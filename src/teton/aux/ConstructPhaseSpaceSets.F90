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
   use cmake_defines_mod,            only: omp_device_num_processors

   use kind_mod
   use constant_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use Quadrature_mod
   use BoundaryList_mod
   use Boundary_mod
   use GreyAcceleration_mod
   use SetData_mod
   use AngleSet_mod
   use GroupSet_mod
   use CommSet_mod
   use ZoneSet_mod
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
   integer :: nZoneSets
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
   integer :: s

   integer :: totalSets
   integer :: angleSetsUsed
   integer :: angleSetsLeft
   integer :: groupsPerSet
   integer :: extraGroups
   integer :: n
   integer :: angleSetsPerSet 
   integer :: nPerSet
   integer :: nAngles
   integer :: nReflecting

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

   integer :: nHypDomMin
   integer :: nHypDomMax

   logical (kind=1) :: GTASet

   real(adqt) :: dot

   integer, allocatable   :: setGroups(:)
   integer, allocatable   :: setAngles(:)
   integer, allocatable   :: setGroup0(:)
   integer, allocatable   :: setGroupID(:)
   integer, allocatable   :: setAngle0(:)

   logical :: verbose

!  Construct Set Data

   nReflecting = getNumberOfReflecting(RadBoundary)
   nSets       = getNumberOfSets(Quad)
   nGroupSets  = 1
   nZoneSets   = getNumberOfZoneSets(Quad)

   verbose = nSets > 1 .AND. Options%isRankVerbose() > 0

!  Decompose angle sets (this finds the maximum number of angle sets
!  allowed respecting angular dependencies). Decomposition in angle
!  minimizes run time and gives the best thread scaling so we do this first.
   call decomposeAngleSets

!  Determine maximum number of phase-space sets problem will support.
!  Until we add "zone sets", the maximum number of sets we can use for the
!  sweep is the number of angle sets multiplied by the number of groups

   QuadSet    => getQuadrature(Quad, 1)
   nAngleSets =  QuadSet% maxAngleSets
   nSetsMax   =  QuadSet% maxAngleSets*QuadSet% Groups
   QuadID     =  1

   if (verbose) then
      print "(A, I0, A, I0, A, I0, A)", "Teton: Quadrature set ", QuadID, " supports sweeping up to ", QuadSet%maxAngleSets, " sweep directions concurrently and has ", QuadSet%Groups, " energy group bins."
   endif

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
   endif

   totalSets = 0

   ! Create only one phase-space set
   MaxSetTest: if ( nSets == 1) then

     setGroups(1)  =  QuadSet% Groups
     setAngles(1)  =  QuadSet% numAngles
     setGroup0(1)  =  0
     setGroupID(1) =  1
     setAngle0(1)  =  0

     totalSets = 1

   else

   ! Create multiple phase-space sets

     if (nSets == nAngleSets) then

!  Assign one angle set to each "set". Each set has the same number of groups (ngr). 

       angle0 = 0

       do s=1,nAngleSets
         setGroups(s)  = QuadSet% Groups
         setGroup0(s)  = 0
         setGroupID(s) = 1
         setAngles(s)  = QuadSet% angleSetSize(s)
         setAngle0(s)  = angle0
         angle0        = angle0 + QuadSet% angleSetSize(s)
       enddo

       totalSets  = nAngleSets
       nGroupSets = 1

     elseif (nSets > nAngleSets) then

!      If the number of sets desired is greater than the number of angle sets,
!      decompose further in energy and distribute the work as balanced as possible
!

       ! Keep doubling the sets until we get the closest we can to nSets
       ! without exceeding it. The assumption here is that each group set
       ! contains the same number of groups. This should be relaxed in the
       ! future.

       nBalancedSets = nAngleSets
       nGroupSets    = 1

       do while ( (nBalancedSets * 2 <= nSets) .AND. (nGroupSets *2 <= QuadSet%maxGroupSets))
         nGroupSets    = nGroupSets * 2
         nBalancedSets = nBalancedSets * 2
       enddo

!      Here nSets = nAngleSets*nGroupSets
       nSets = nBalancedSets

!      The following code block handles the case where the groups sets are
!      unbalanced (i.e. not all group sets contain the same number of groups).

       groupsPerSet = int(QuadSet% Groups/nGroupSets)
       extraGroups  = QuadSet% Groups - nGroupSets*groupsPerSet 

!      If there are extra groups assign one more to the first
!      "extraGroups" group sets 

       totalSets = 0
       angle0    = 0

       do s=1,nAngleSets

         g0 = 0
         do groupSetID=1,nGroupSets
           if (groupSetID <= extraGroups) then
             Groups = groupsPerSet + 1
           else
             Groups = groupsPerSet
           endif

           setGroups(totalSets+groupSetID)  = Groups
           setGroup0(totalSets+groupSetID)  = g0
           setGroupID(totalSets+groupSetID) = groupSetID 
           setAngles(totalSets+groupSetID)  = QuadSet% angleSetSize(s)
           setAngle0(totalSets+groupSetID)  = angle0
           g0                               = g0 + Groups
         enddo

         totalSets = totalSets + nGroupSets
         angle0    = angle0    + QuadSet% angleSetSize(s)

       enddo

     elseif (nSets < nAngleSets) then

!      At present, phase-space sets must contain the same number of angles in each set for the gpu.
!      Adjust nSets so that nSets divides evenly into nAngleSets
       if (nSets /= nAngleSets) then
         do while ( mod(nAngleSets, nSets) /= 0)
           nSets = nSets - 1
         enddo
       endif

!      The following code block allows for the case where the
!      number of angle sets is not an even multiple of the number
!      of sets. 

       angleSetsPerSet = max( int(QuadSet% maxAngleSets/nSets), 1)
       angleSetsUsed   = angleSetsPerSet*nSets
       angleSetsLeft   = QuadSet% maxAngleSets - angleSetsUsed 

!      If the angle sets cannot be evenly divided (angleSetsLeft /= 0) 
!      we assign the extra angle sets to the first few sets until
!      they are gone

       angle0 = 0
       setID  = 0

       do s=1,nSets
         setGroups(s)  = QuadSet% Groups
         setGroup0(s)  = 0
         setGroupID(s) = 1
         setAngle0(s)  = angle0

         if (s <= angleSetsLeft) then
           nPerSet = angleSetsPerSet + 1
         else
           nPerSet = angleSetsPerSet
         endif

         nAngles = 0
         do n=1,nPerSet
           setID   = setID + 1
           nAngles = nAngles + QuadSet% angleSetSize(setID)
         enddo

         setAngles(s) = nAngles
         angle0       = angle0 + nAngles
       enddo

       nAngleSets = nSets
       totalSets  = nSets

     endif

!    Error check

     if (totalSets /= nSets) then
       call f90fatal("ConstructPhaseSpaceSets: totalSets /= nSets")
     endif

   endif MaxSetTest

   if (verbose) then
      print "(A,I0,A)", "Teton: Distributing angles and groups across ", nSets, " phase-space sets..."
   endif

   do setID=1,nSets
     QuadID = 1
     if (verbose) then
        write(stdout,100) setID,QuadID,setAngles(setID),setAngle0(setID)+1,setAngle0(setID)+setAngles(setID),setGroups(setID),setGroup0(setID)+1,setGroup0(setID)+setGroups(setID)
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

!  Allocate pointers for the Set, Angle Set, Group Set, Communication Set
!  and GTA Set modules

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

     QuadID     =  1 
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
                      nReflecting, GTASet, QuadSet)
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

       QuadID     =  2 
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
                      nReflecting, GTASet, QuadSet)

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

   allocate(ZSet)

   call construct(ZSet, nZoneSets)

!   Determine the number of "hyper-domains" to increase parallelism
!   for GTA. We need to be careful for very small zone counts so we
!   estimate a minimum number based on the number of zones.

!   Note that the use of "hyper-domains" will be deprecated once
!   we support sub-meshes per MPI rank.  PFN 09/22/2022

    nHypDomMin = int( 2*sqrt( real(Size%nzones) ) - 1 )
    nHypDomMin = max( nHypDomMin, 1 )
    nHypDomMax = int( omp_device_num_processors/max(nGTASets,1) )
    nHypDomMax = min( nHypDomMax, 20 )

    if (Size% useGPU) then
      Quad% nHyperDomains = min(nHypDomMax, nHypDomMin)
    else
      Quad% nHyperDomains = min(nSets, nHypDomMin)
    endif

!  Release memory

   deallocate( setGroups )
   deallocate( setAngles )
   deallocate( setGroup0 )
   deallocate( setGroupID )
   deallocate( setAngle0 )

   return
   end subroutine ConstructPhaseSpaceSets

