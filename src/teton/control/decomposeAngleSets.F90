!***********************************************************************
!                        Last Update:  09/2017, PFN                    *
!                                                                      *
!  decomposeAngleSets -   Decompose all angles into subsets that       *
!                         may be solved concurrently.                  *I
!                                                                      *
!  20250409 BCY - Each angle set created in ConstructPhaseSpaceSets    *
!     will consist of one of more of these angle subsets.              *
!   In the extreme ends, we could have 1 angle set with all of the     *
!   angle subsets (no rank-local parallelism over angle), or           *
!   numAngleSets = maxAngleSets where each angle subset created here   *
!   ends up in its own angle set.                                      *
!                                                                      *
!***********************************************************************


   subroutine decomposeAngleSets

!  Include

   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use BoundaryList_mod
   use Boundary_mod
   use QuadratureList_mod
   use Quadrature_mod


   implicit none

!  Local

   integer         :: set
   integer         :: nReflecting
   integer         :: nReflectGlobal
   integer         :: NumAngles
   integer         :: ndim
   integer         :: bdyID
   integer         :: angle
   integer         :: angle0
   integer         :: angleRef
   integer         :: angleStride
   integer         :: angleSet
   integer         :: bin
   integer         :: binRef
   integer         :: octant
   integer         :: octantRef
   integer         :: anglesPerOctant
   integer         :: setSize
   integer         :: nSets

   integer         :: List(8)
   integer         :: muList(8)
   integer         :: etaList(8)
   integer         :: xsiList(8)
   integer         :: etaXsiList(8)
   integer         :: dirCoupled(Size%ndim)

   real(adqt), parameter :: eps=1.0e-15_adqt
   real(adqt)            :: OmegaDotA
   real(adqt)            :: Area(Size% ndim)

   integer, allocatable  :: octantCoupled(:,:)

   logical(kind=1) :: muCoupled
   logical(kind=1) :: etaCoupled
   logical(kind=1) :: xsiCoupled

   data  muList    /1,2,3,4,5,6,7,8/
   data etaList    /1,4,2,3,5,8,6,7/
   data xsiList    /1,5,2,6,3,7,4,8/
   data etaXsiList /1,4,5,8,2,3,6,7/

!  Constants 

   nReflecting    = getNumberOfReflecting(RadBoundary)
   nReflectGlobal = 0
   nSets          = getNumberOfSets(Quad)
   ndim           = Size% ndim
   angleStride    = 1

   allocate( octantCoupled(8,nReflecting) )

   dirCoupled(1:Size%ndim) = 0
   octantCoupled(:,:)      = 0
   List(:)                 = 0

!  The code supports one high-order angle set and one GTA set

   QuadSetLoop: do set=1,2

     QuadSet   => getQuadrature(Quad, set)
     NumAngles =  QuadSet% NumAngles

     ! Skip angle decomposition for following cases:
     ! - 1d mesh (not supported)
     ! - quadrature set is level symmetric (not supported)
     ! - user specified to use only one angle set, or code is not threaded
     if (ndim == 1 .or. (set /= 2 .and. QuadSet% TypeName == 'levelsym') .or. nSets == 1 ) then

       QuadSet% maxAngleSets = 1

       do angle=1,NumAngles
         QuadSet% angleList(angle) = angle
       enddo

       QuadSet% angleSetSize(1) = NumAngles

     elseif (ndim == 2) then

       xsiCoupled = .FALSE.

!  Loop over reflecting-boundary sets:

       BoundaryLoop2D: do bdyID=1,nReflecting

         Bdy     => getReflecting(RadBoundary, bdyID)
         Area(:) =  Bdy% A_bdy(:,1)

         AngleLoop2D: do angle=1,NumAngles

           OmegaDotA = DOT_PRODUCT( QuadSet% omega(:,angle),Area(:) )

           if (OmegaDotA < -eps) then

             call reflectAxis(ndim, angle, NumAngles, angleRef, &
                              QuadSet% omega, Area)

             bin    = QuadSet% AngleToBin(angle)
             binRef = QuadSet% AngleToBin(angleRef)

             if (bin /= binRef) then
               dirCoupled(2) = 1
             endif

           endif

         enddo AngleLoop2D

       enddo BoundaryLoop2D

       call MPIAllReduce(dirCoupled, "max", MY_COMM_GROUP)

       if ( dirCoupled(2) > 0 ) then
         xsiCoupled = .TRUE.
       endif

       if ( xsiCoupled ) then
         QuadSet% maxAngleSets = QuadSet% nPolarAngles/2

         angle0 = 0
         do angleSet=1,QuadSet% maxAngleSets
           setSize = QuadSet% NangBinList(2*angleSet-1) +  &
                     QuadSet% NangBinList(2*angleSet)

           QuadSet% angleSetSize(angleSet) = setSize

           do angle=1,setSize
             QuadSet% angleList(angle0+angle) = angle0 + angle
           enddo

           angle0 = angle0 + setSize
         enddo

       else
         QuadSet% maxAngleSets = QuadSet% nPolarAngles

         angle0 = 0
         do angleSet=1,QuadSet% maxAngleSets
           setSize = QuadSet% NangBinList(angleSet)
           QuadSet% angleSetSize(angleSet) = setSize

           do angle=1,setSize
             QuadSet% angleList(angle0+angle) = angle0 + angle
           enddo

           angle0 = angle0 + setSize
         enddo

       endif

     elseif (ndim == 3) then

       anglesPerOctant =  NumAngles/8
       muCoupled       = .FALSE.
       etaCoupled      = .FALSE.
       xsiCoupled      = .FALSE.

!  Loop over reflecting-boundary sets.  In XYZ, the angles are
!  numbered consecutively within an octant so we only need to
!  test one angle per octant.

       BoundaryLoop3D: do bdyID=1,nReflecting

         Bdy     => getReflecting(RadBoundary, bdyID)
         Area(:) =  Bdy% A_bdy(:,1)

         octantLoop3D: do octant=1,8

           angle    = octant

           OmegaDotA = DOT_PRODUCT( QuadSet% omega(:,angle),Area(:) )

           if (OmegaDotA < -eps) then

             call reflectAxis(ndim, angle, NumAngles, angleRef, &
                              QuadSet% omega, Area)

             octantRef = mod(angleRef,8)

             if (octantRef == 0) then
               octantRef = 8
             endif

             octantCoupled(octant,bdyID)    = octantRef
             octantCoupled(octantRef,bdyID) = octant

           endif

         enddo octantLoop3D

         if (octantCoupled(1,bdyID) == 2) then
           dirCoupled(1) = 1
         endif

         if (octantCoupled(1,bdyID) == 4) then
           dirCoupled(2) = 1
         endif

         if (octantCoupled(1,bdyID) == 5) then
           dirCoupled(3) = 1 
         endif

       enddo BoundaryLoop3D

       call MPIAllReduce(dirCoupled, "max", MY_COMM_GROUP)

       nReflectGlobal = sum ( dirCoupled(:) )
   
       if ( dirCoupled(1) == 1 ) then
         muCoupled = .TRUE.
       endif

       if ( dirCoupled(2) == 1 ) then 
         etaCoupled = .TRUE.
       endif

       if ( dirCoupled(3) == 1 ) then
         xsiCoupled = .TRUE.
       endif

       if ( nReflectGlobal == 3) then

         QuadSet% maxAngleSets = anglesPerOctant
         angleStride           = 8

         List(:) = muList(:)

       elseif ( nReflectGlobal == 2) then

         QuadSet% maxAngleSets = 2*anglesPerOctant
         angleStride           = 4

         if ( muCoupled .and. etaCoupled ) then
            List(:) = muList(:)
         elseif (muCoupled .and. xsiCoupled ) then
            List(:) = xsiList(:)
         elseif (etaCoupled .and. xsiCoupled ) then
            List(:) = etaXsiList(:)
         endif

       elseif ( nReflectGlobal == 1) then

         QuadSet% maxAngleSets = 4*anglesPerOctant
         angleStride           = 2

         if ( muCoupled ) then
           List(:) = muList(:)
         elseif ( etaCoupled ) then
           List(:) = etaList(:)
         elseif ( xsiCoupled ) then
           List(:) = xsiList(:)
         endif

       elseif ( nReflectGlobal == 0) then

         QuadSet% maxAngleSets = NumAngles
         angleStride           = 1 

         List(:) = muList(:)

       endif

       angle0 = 0
       do angle=1,anglesPerOctant
         do octant=1,8
           QuadSet% angleList(angle0+octant) = angle0 + List(octant)
         enddo
         angle0 = angle0 + 8
       enddo

       QuadSet% angleSetSize(1:QuadSet% maxAngleSets) = angleStride

     endif

     ! Reduce GTA angle sets to one if it exceeded nSets.

     if ( set == 2 .and. QuadSet%maxAngleSets > nSets) then
       QuadSet% maxAngleSets = 1

       do angle=1, QuadSet%NumAngles
         QuadSet% angleList(angle) = angle
       enddo

       QuadSet% angleSetSize(1) = QuadSet%NumAngles
     endif

   enddo QuadSetLoop

   deallocate( octantCoupled )



   return
   end subroutine decomposeAngleSets
