!***********************************************************************
!                        Last Update:  07/2017, TSH                   *
!                                                                      *
!   ConstructPhaseSpaceSets - Builds group-angle sets to expose more   *
!                             parallelism for threading and/or GPU     *
!                             implementations.                         *
!                                                                      *
!***********************************************************************


   subroutine ConstructPhaseSpaceSets(fromRestart) &
                       BIND(C,NAME="teton_constructphasespacesets")

!  Include

   use, intrinsic :: ISO_C_BINDING
   use, intrinsic:: iso_fortran_env, only: stdout=>output_unit
   use kind_mod
   use constant_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use Geometry_mod
   use ZoneData_mod
   use QuadratureList_mod
   use Quadrature_mod
   use BoundaryList_mod
   use Boundary_mod
   use GreyAcceleration_mod
   use SetData_mod
   use AngleSet_mod
   use GroupSet_mod
   use CommSet_mod
   use Options_mod

   implicit none

!  Arguments

   logical(C_BOOL), intent(in)    :: fromRestart

!  Local

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet
   type(GroupSet),   pointer :: GSet
   type(CommSet),    pointer :: CSet

   integer :: nSets
   integer :: nBalancedSets
   integer :: nSetsMax
   integer :: nGTASets
   integer :: nAngleSets
   integer :: nGroupSets
   integer :: nCommSets
   integer :: NumSnSets
   integer :: nZoneSets
   integer :: qset
   integer :: set1
   integer :: setID
   integer :: QuadID
   integer :: groupSetID
   integer :: Groups
   integer :: angleSetID
   integer :: commSetID
   integer :: NumAngles
   integer :: nZones
   integer :: nCorner
   integer :: g
   integer :: g0
   integer :: angle
   integer :: angle0
   integer :: angle1
   integer :: s
   integer :: subset

   integer :: totalWork
   integer :: totalSets
   integer :: workPerSet
   integer :: setsUsed
   integer :: setsLeft
   integer :: extraWork
   integer :: extraGroups
   integer :: n
   integer :: nPerSet
   integer :: nAngles
   integer :: nReflecting

   integer :: sharedID
   integer :: nShared
   integer :: NumBdyElem
   integer :: NEW_COMM_GROUP
   integer :: new_group
   integer :: new_comm 
   integer :: ierror

   integer :: reflID
   integer :: nExit

   integer :: cSet1
   integer :: cSet2
   integer :: z1
   integer :: z2

   logical (kind=1) :: GTASet

   real(adqt) :: workFraction
   real(adqt) :: dot

   integer, dimension (1) :: setMaxWork 

   integer, allocatable   :: AngleSetWork(:)
   integer, allocatable   :: setsPerAngleSet(:)
   integer, allocatable   :: offset(:)
   integer, allocatable   :: workPerQuad(:)
   integer, allocatable   :: setQuadID(:)
   integer, allocatable   :: setGroups(:)
   integer, allocatable   :: setAngles(:)
   integer, allocatable   :: setGroup0(:)
   integer, allocatable   :: setGroupID(:)
   integer, allocatable   :: setAngle0(:)

   logical :: verbose

!  Construct Set Data

   nShared     = getNumberOfShared(RadBoundary)
   nReflecting = getNumberOfReflecting(RadBoundary)
   nSets       = getNumberOfSets(Quad)
   NumSnSets   = getNumSnSets(Quad)
   nGroupSets  = 1
   nZoneSets   = getNumberOfZoneSets(Quad)

   verbose = nSets > 1 .AND. Options%isRankVerbose() > 0

   nSetsMax = 0
   nAngleSets = 0

!  Decompose angle sets (this finds the maximum number of angle sets
!  allowed respecting angular dependencies). Decomposition in angle
!  minimizes run time and gives the best thread scaling so we do this first.
   call decomposeAngleSets

!  Determine maximum number of phase-space sets problem will support.
!  Until we add "zone sets", the maximum number of sets we can use for the
!  sweep is the number of angle sets multiplied by the number of groups
   do qset=1,NumSnSets
     QuadSet    => getQuadrature(Quad, qset)
     nAngleSets =  nAngleSets + QuadSet% maxAngleSets
     nSetsMax   =  nSetsMax   + QuadSet% maxAngleSets*QuadSet% Groups

     if (verbose) then
        print "(A, I0, A, I0, A, I0, A)", "Teton: Quadrature set ", qset, " supports sweeping up to ", QuadSet%maxAngleSets, " sweep directions concurrently and has ", QuadSet%Groups, " energy group bins."
     endif
   enddo

   allocate( setQuadID(nSets) )
   allocate( setGroups(nSets) )
   allocate( setAngles(nSets) )
   allocate( setGroup0(nSets) )
   allocate( setGroupID(nSets) )
   allocate( setAngle0(nSets) )

   ! Reduce nSets if it exceeds what problem will support.
   if (nSets > nSetsMax) then
     if (verbose) then
        print "(A, I0, A, I0, A)", "Teton: This problem lacks enough parallelism to create ", nSets, " phase-space sets.  The maximum available (", nSetsMax, ") will be created."
     endif

     nSets = nSetsMax

     ! TODO - Ask Paul for explanation on what this does.  Resizing some
     ! structures to account for new value of nSets? -- Aaron
     s = 0
     do qset=1,NumSnSets
       QuadSet => getQuadrature(Quad, qset)

       angle0 = 0

       do subSet=1,QuadSet% maxAngleSets
         g0 = 0
         do g=1,QuadSet% Groups
           s             = s + 1
           setQuadID(s)  = qset
           setGroups(s)  = 1
           setGroup0(s)  = g0
           setGroupID(s) = 1
           setAngles(s)  = QuadSet% angleSetSize(subSet)
           setAngle0(s)  = angle0
           g0            = g0 + 1
         enddo
         angle0 = angle0 + QuadSet% angleSetSize(subSet)
       enddo
     enddo
   endif

   ! Create only one phase-space set
   MaxSetTest: if ( nSets == 1) then

     QuadSet       => getQuadrature(Quad, 1)
     setQuadID(1)  =  1
     setGroups(1)  =  QuadSet% Groups
     setAngles(1)  =  QuadSet% numAngles
     setGroup0(1)  =  0
     setGroupID(1) =  1
     setAngle0(1)  =  0

   ! Create multiple phase-space sets
   else
     allocate( AngleSetWork(nAngleSets) )
     allocate( setsPerAngleSet(nAngleSets) )
     allocate( offset(nAngleSets) )
     allocate( workPerQuad(NumSnSets) )

     AngleSetWork(:) = 0
     totalWork       = 0
     s               = 0

     do qset=1,NumSnSets
       QuadSet           => getQuadrature(Quad, qset)
       workPerQuad(qset) = 0

       do subSet=1,QuadSet% maxAngleSets
         s                 = s + 1
         AngleSetWork(s)   = QuadSet% angleSetSize(subSet)*QuadSet% Groups
         workPerQuad(qset) = workPerQuad(qset) + AngleSetWork(s)
         totalWork         = totalWork         + AngleSetWork(s)
       enddo
     enddo

     if (nSets >= nAngleSets) then

       ! Note: Only one Sn Quad set is supported in Teton.  The rest of the
       ! code in ConstructPhaseSpaceSets needs to be refactored to only support
       ! just one Sn Quad set.
       QuadSet           => getSNQuadrature(Quad)

!      If the number of sets desired is greater than the number of angle sets,
!      decompose further in energy and distribute the work as balanced as possible
       nBalancedSets = nAngleSets
!

       ! Keep doubling the sets until we get the closest we can to nSets
       ! without exceeding it.
       nGroupSets=1

       do while ( (nBalancedSets * 2 <= nSets) .AND. (nGroupSets *2 <= QuadSet%maxGroupSets))
         nGroupSets = nGroupSets * 2
         nBalancedSets = nBalancedSets * 2
       enddo

       nSets = nBalancedSets

       workPerSet = totalWork/nSets
       extraWork  = totalWork - workPerSet*nSets

!      Start with one set per angle subset
       setsPerAngleSet(:)  = 1
       setsUsed            = nAngleSets
       setsLeft            = nSets - nAngleSets

       AngleSetWork(:) = AngleSetWork(:) - workPerSet

       SetIteration: do s=1,setsLeft
         setMaxWork             = maxloc( AngleSetWork(:) )
         setID                  = setMaxWork(1)
         setsPerAngleSet(setID) = setsPerAngleSet(setID) + 1
         AngleSetWork(setID)    = AngleSetWork(setID) - workPerSet

         setsUsed = setsUsed + 1

         if (setsUsed == nSets) then
           exit SetIteration
         else
           cycle SetIteration
         endif

       enddo SetIteration

       totalSets = 0
       s         = 0

       do qset=1,NumSnSets
         QuadSet => getQuadrature(Quad, qset)

         angle0  =  0

         do subSet=1,QuadSet% maxAngleSets
           s           = s + 1
           Groups      = QuadSet% Groups/setsPerAngleSet(s)
           extraGroups = QuadSet% Groups - Groups*setsPerAngleSet(s)

           g0 = 0
           do setID=1,setsPerAngleSet(s)
             setQuadID(totalSets+setID)  = qset
             setGroups(totalSets+setID)  = Groups
             setGroup0(totalSets+setID)  = g0
             setGroupID(totalSets+setID) = setID 
             setAngles(totalSets+setID)  = QuadSet% angleSetSize(subSet)
             setAngle0(totalSets+setID)  = angle0
             g0                          = g0 + Groups
           enddo

           do setID=1,extraGroups
             setGroups(totalSets+setID) = setGroups(totalSets+setID) + 1
           enddo

           g0 = 0
           do setID=1,setsPerAngleSet(s)
             setGroup0(totalSets+setID) = g0
             g0                         = g0 + setGroups(totalSets+setID)
           enddo

           totalSets = totalSets + setsPerAngleSet(s)
           angle0    = angle0    + QuadSet% angleSetSize(subSet)

         enddo
       enddo

     elseif (nSets <= nAngleSets) then

!      At present, phase-space sets must contain the same number of angles in each set for the gpu.
!      Adjust nSets so that nSets divides evenly into nAngleSets
       if (nSets /= nAngleSets) then
         do while ( mod(nAngleSets, nSets) /= 0)
           nSets = nSets - 1
         enddo
       endif

!      Combine angle sets and balance the work
       do qset=1,NumSnSets
         workFraction          = real(workPerQuad(qset)/totalWork, adqt)
         setsPerAngleSet(qset) = max(int(workFraction*nSets), 1)
       enddo

       totalSets = 0
       s         = 0

       do qset=1,NumSnSets
         QuadSet => getQuadrature(Quad, qset)

         angle0   = 0
         nPerSet  = QuadSet% maxAngleSets/setsPerAngleSet(qset)
         setsUsed = nPerSet*setsPerAngleSet(qset)
         setsLeft = QuadSet% maxAngleSets - setsUsed
         setID    = 0

         do subSet=1,setsPerAngleSet(qset)
           s = s + 1
           setQuadID(s)  = qset
           setGroups(s)  = QuadSet% Groups
           setGroup0(s)  = 0
           setGroupID(s) = 1
           setAngle0(s)  = angle0

           nAngles = 0
           do n=1,nPerSet
             setID   = setID + 1
             nAngles = nAngles + QuadSet% angleSetSize(setID)
           enddo

           setAngles(s) = nAngles
           angle0       = angle0 + nAngles

         enddo

         offset(1) = 0
         do s=2,QuadSet% maxAngleSets
           offset(s) = offset(s-1) + QuadSet% angleSetSize(s-1)
         enddo

         set1   = setsPerAngleSet(qset) - setsLeft + 1

         if (set1 <= setsPerAngleSet(qset) ) then
           angle0 = setAngle0(set1)
           angle1 = 0

           do s=1,set1-1
             do angle=1,QuadSet% angleSetSize(s)
               QuadSet% angleList(angle1+angle) = offset(s) + angle
             enddo
             angle1 = angle1 + QuadSet% angleSetSize(s)
           enddo

           do s=set1,setsPerAngleSet(qset)
             setID        = setID + 1

             do angle=1,QuadSet% angleSetSize(s)
               QuadSet% angleList(angle1+angle) = offset(s) + angle
             enddo
             angle1 = angle1 + QuadSet% angleSetSize(s)

             do angle=1,QuadSet% angleSetSize(setID)
               QuadSet% angleList(angle1+angle) = offset(setID) + angle
             enddo
             angle1 = angle1 + QuadSet% angleSetSize(setID)

             nAngles      = setAngles(s) + QuadSet% angleSetSize(setID)
             setAngles(s) = nAngles
             setAngle0(s) = angle0
             angle0       = angle0 + nAngles
           enddo
         endif

         totalSets = totalSets + setsPerAngleSet(qset)

       enddo

       nAngleSets = totalSets

     endif

     deallocate( AngleSetWork )
     deallocate( setsPerAngleSet )
     deallocate( offset )
     deallocate( workPerQuad )

!    Error check

     if (totalSets /= nSets) then
       call f90fatal("ConstructPhaseSpaceSets: totalSets /= nSets")
     endif

   endif MaxSetTest

   if (verbose) then
      print "(A,I0,A)", "Teton: Distributing angles and groups across ", nSets, " phase-space sets..."
   endif

   do setID=1,nSets
     if (verbose) then
        write(stdout,100) setID,setQuadID(setID),setAngles(setID),setAngle0(setID)+1,setAngle0(setID)+setAngles(setID),setGroups(setID),setGroup0(setID)+1,setGroup0(setID)+setGroups(setID)
     endif
     ! Remove after support is added for phase-space sets with different numbers of
     ! angles and groups and we have tests exercising this in the suite.
     if ( setID > 1) then
       if ( setGroups(setID) /= setGroups(setID-1) ) then
         call f90fatal("Teton: Unable to evenly distribute energy groups across phase-space sets.  This is currently a requirement. Contact the Teton team for tips on adjusting your energy groups to allow even distribution over the phase-space sets.")
       endif
       if ( setAngles(setID) /= setAngles(setID-1) ) then
         call f90fatal("Teton: Unable to evenly distribute angles across phase-space sets. This is currently a requirement.  Contact the Teton team for tips on adjusting your angle setup to allow even distribution over the phase-space sets.")
       endif
     endif
   enddo

!  Allocate set pointers for the SetData and RadIntensity modules

   if (Size% ndim == 1) then
     nGTASets  = 0
     nCommSets = nSets
   else
     QuadSet   => getGTAQuadrature(Quad)
     nGTASets  =  QuadSet% maxAngleSets
     nCommSets =  nAngleSets
   endif

   call constructSetPointers(Quad, nSets, nAngleSets, nGroupSets,  &
                             nCommSets, nGTASets)

!  Construct the phase-space sets

   GTASet = .FALSE.

   SetLoop: do setID=1,nSets

     Set => getSetData(Quad, setID)

     QuadID     =  setQuadID(setID)
     Groups     =  setGroups(setID)
     NumAngles  =  setAngles(setID)
     g0         =  setGroup0(setID)
     groupSetID =  setGroupID(setID)
     angle0     =  setAngle0(setID)
     angleSetID =  angle0/NumAngles + 1

     if (Size% ndim == 1) then
       commSetID = setID
     else
       commSetID = angleSetID 
     endif

!    In the near future, the number of corners could be different for each set
     nZones     =  Size% nzones
     nCorner    =  Size% ncornr

     QuadSet    => getQuadrature(Quad, QuadID)
     ASet       => getAngleSetData(Quad, angleSetID)
     GSet       => getGroupSetData(Quad, groupSetID)
     CSet       => getCommSetData(Quad, commSetID)

     Quad% groupID(setID) = groupSetID
     Quad% angleID(setID) = angleSetID
     Quad% commID(setID)  = commSetID

!    Construct the set
     call Set%construct(setID, groupSetID, angleSetID, QuadID,      &
                        Groups, NumAngles, g0, angle0, nZones, nCorner,  &
                        QuadSet, GTASet, fromRestart)

!    Construct group sets, but only for the first angle set
     if (angle0 == 0) then
       call construct(GSet, Groups, g0, nZones, nCorner)
     endif

!    Construct angle sets, but only for the first group set
     if (groupSetID == 1) then
       call construct(ASet, NumAngles, angle0, nZones,  &
                      nReflecting, QuadSet)
     endif

!    Construct communication sets - a unique communication group
!    is created for each communication set

     if (groupSetID == 1 .or. Size% ndim == 1) then

!      duplicate the existing communicator
       call MPI_COMM_DUP(MY_COMM_GROUP, NEW_COMM_GROUP, ierror)

!      extract the original group handle
       call MPI_COMM_GROUP(NEW_COMM_GROUP, new_group, ierror)

!      create new communicator
       call MPI_COMM_CREATE(NEW_COMM_GROUP, new_group, new_comm, ierror)

       if (ierror /= MPI_SUCCESS) then
          call f90fatal("MPI COMM Create Failed")
       endif

     endif

     if (Size% ndim > 1) then

       if (groupSetID == 1) then
         cSet1 = setID
         cSet2 = setID + nGroupSets - 1

         call construct(CSet, NumAngles, cSet1, cSet2, Groups, new_comm, ASet)
       endif

     else

       cSet1 = setID
       cSet2 = setID

       call construct(CSet, NumAngles, cSet1, cSet2, Groups, new_comm, ASet)

     endif

     do angle=1,NumAngles
        do g=1,Groups
          Quad% SetIDList(g0+g,angle0+angle) = setID
        enddo
     enddo

   enddo SetLoop


   AngleSetLoop: do setID=1,nAngleSets

     ASet      => getAngleSetData(Quad, setID)
     NumAngles =  ASet% NumAngles

     do reflID=1,nReflecting

       Bdy   => getReflecting(RadBoundary, reflID)
       nExit =  0

       do angle=1,NumAngles
         dot = DOT_PRODUCT( ASet% omega(:,angle),Bdy% A_bdy(:,1) )

         if (dot > zero) then
           nExit = nExit + 1
           ASet% ExitAngleList(nExit,reflID) = angle
         endif
       enddo

       ASet% nExit(reflID) = nExit

     enddo

   enddo AngleSetLoop

!  GTA Set

   if (Size% ndim > 1) then

     GTASet       = .TRUE.
     angle0       =  0
     angleSetID   =  nAngleSets
     groupSetID   =  1

     if (verbose) then
        print "(A)", "Teton: Angle and energy group distribution breakdown (grey acceleration sweep):"
     endif

     do setID=1,nGTASets

       angleSetID =  angleSetID + 1
       commSetID  =  nAngleSets + setID

       Set        => getGTASetData(Quad, setID)
       ASet       => getAngleSetData(Quad, angleSetID)
       CSet       => getCommSetData(Quad, commSetID)
       QuadSet    => getGTAQuadrature(Quad)

       Quad% angleID(nSets+setID) = angleSetID
       Quad% commID(nSets+setID)  = commSetID

       QuadID     =  NumSnSets + 1 
       Groups     =  1 
       NumAngles  =  QuadSet% angleSetSize(setID) 
       g0         =  0
       nZones     =  Size% nZones
       nCorner    =  Size% ncornr

       if (verbose) then
         write(stdout,100) setID,QuadID,NumAngles,angle0+1,angle0+NumAngles,Groups,g0+1,g0+Groups
       endif

!      construct the GTA set
       call Set%construct(setID, groupSetID, angleSetID, QuadID, &
                          Groups, NumAngles, g0, angle0, nZones,      &
                          nCorner, QuadSet, GTASet, fromRestart)

!      construct an angle set for every GTA set
       call construct(ASet, NumAngles, angle0, nZones,  &
                      nReflecting, QuadSet)

!      construct a communication set for every GTA set

!      duplicate the existing communicator
       call MPI_COMM_DUP(MY_COMM_GROUP, NEW_COMM_GROUP, ierror)

!      extract the original group handle
       call MPI_COMM_GROUP(NEW_COMM_GROUP, new_group, ierror)

!      create new communicator
       call MPI_COMM_CREATE(NEW_COMM_GROUP, new_group, new_comm, ierror)

       cSet1 = nSets + setID 
       cSet2 = nSets + setID 

       call construct(CSet, NumAngles, cSet1, cSet2, Groups, new_comm, ASet)

       angle0    = angle0 + NumAngles

     enddo

!    Construct and incident test on shared boundaries

     call initFindExit(nAngleSets, nGTASets)

   endif

 100 format("       Phase-Angle Set ID =",i3,2x," | Quadrature Set ID =",i2,2x, " | # Angles = ",i3," | Angle IDs =",i3," -",i3, " | # Groups =",i3," | Group IDs = ",i3," -",i3)

!  Grey Acceleration Module
!  Moving this constructor here because (in the near future)
!  we will taylor an acceleration method for each set

   allocate (GTA)

   call construct(GTA)

!  Zone Sets - create a start and end for the corner lists

   do setID=1,nZoneSets
     z1 = Geom% zone1(setID)
     z2 = Geom% zone2(setID)

     Geom% corner1(setID) = Geom% cOffSet(z1) + 1
     Geom% corner2(setID) = Geom% cOffSet(z2) + Geom% numCorner(z2)
   enddo

!  Release memory

   deallocate( setQuadID )
   deallocate( setGroups )
   deallocate( setAngles )
   deallocate( setGroup0 )
   deallocate( setGroupID )
   deallocate( setAngle0 )

   return
   end subroutine ConstructPhaseSpaceSets

