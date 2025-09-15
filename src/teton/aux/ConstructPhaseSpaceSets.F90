#include "macros.h"
!***********************************************************************
!                        Last Update:  04/2025 BCY                     *
!                                                                      *
!   ConstructPhaseSpaceSets - Builds group-angle sets to expose more   *
!                             parallelism for threading and/or GPU     *
!                             implementations.                         *
!                                                                      *
! 20250414 BCY -- I removed the logic for supporting imbalanced phase  *
!                 space sets.  All sets must now have the same number  *
!                 of groups and angles.  It wasn't working in the code *
!                 anyway because it was incompatible with some         *
!                 of the code downstream.  If we want to bring back    *
!                 this logic in the future, see the tag                *
!                         imbalanced-sets                              *
!                 or check out the hash                                *
!                     a263dc521fa92a0057802c6bc71cfa1887dc6e4b         *
!                                                                      *
!***********************************************************************


   subroutine ConstructPhaseSpaceSets(fromRestart) &
                       BIND(C,NAME="teton_constructphasespacesets")

!  Include

   use, intrinsic :: ISO_C_BINDING
   use, intrinsic:: iso_fortran_env, only: stdout=>output_unit
   use Datastore_mod, only : theDatastore

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
   integer :: nSetsMax
   integer :: nSetsMaxUser
   integer :: nGTASets
   integer :: nAngleSets
   integer :: maxAngleSets
   integer :: nGroupSets
   integer :: maxGroupSets
   integer :: nCommSets
   integer :: nZoneSets
   integer :: nHyperDomains
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

   integer :: totalSets
   integer :: nReflecting

   integer :: minGroupSetSize
   integer :: groupsPerSet
   integer :: extraGroups

   integer :: angleSubset, angleSubset0, nAngleSubsets
   integer :: angleSubsetsPerSet
   integer :: anglesInSet, anglesPerSet
   integer :: extraAngles

   integer :: new_comm
   integer :: ierror

   integer :: reflID
   integer :: nExit

   integer :: cSet1
   integer :: cSet2
   integer :: z1
   integer :: z2

   integer :: nHypDomMin
   integer :: nSweepHypDomMax
   integer :: nGTAHypDomMax
   integer :: omp_device_num_processors
   integer :: min_groupset_size


   logical (kind=1) :: GTASet

   real(adqt) :: dot

   integer, allocatable   :: setGroups(:)
   integer, allocatable   :: setAngles(:)
   integer, allocatable   :: setGroup0(:)
   integer, allocatable   :: setGroupID(:)
   integer, allocatable   :: setAngle0(:)

   logical :: verbose

!  Construct Set Data

   nReflecting  = getNumberOfReflecting(RadBoundary)
   nSetsMaxUser = getNumberOfSets(Quad)
   nZoneSets    = getNumberOfZoneSets(Quad)

   verbose = nSetsMaxUser > 1 .AND. Options%isRankVerbose() > 0

!  decomposeAngleSets breaks up the angles as much as possible into subsets
!  Each angle set will have one or more of these subsets
!    i.e., maximum decomposition: 1 angle set = 1 angle subset from decomposeAngleSets
!          no angular decomposition: 1 angle set total with all the angle subsets (i.e., all the angles in one angle set)
   call decomposeAngleSets
!  The above call computed Quad%maxAngleSets.

!  Determine the minimum group set size:
   minGroupSetSize = Options%getMinGroupSetSize()

!!!!!!!
!  Determine maximum number of phase-space sets problem will support.
   QuadSet    => getQuadrature(Quad, 1)
   nSetsMax   =  QuadSet% maxAngleSets*max(int(QuadSet% Groups/minGroupSetSize),1)
   QuadID     =  1

   if (verbose) then
      print "(A, I0, A, I0, A, I0, A)", "Teton: Quadrature set ", QuadID, " supports sweeping up to ", QuadSet%maxAngleSets, " sweep directions concurrently and has ", QuadSet%Groups, " energy group bins."
   endif

   ! Reduce nSets if it exceeds what problem will support.
   if (nSetsMaxUser > nSetsMax) then
     if (verbose) then
        print "(A, I0, A, I0, A)", "Teton: This problem lacks enough parallelism to create ", nSetsMaxUser, " phase-space sets.  The maximum available (up to ", nSetsMax, ") will be created."
     endif

     nSets = nSetsMax
   else
     if (verbose) then
        print "(A, I0, A)", "Teton: Will create up to ", nSetsMaxUser, " phase-space sets (limit requested by user)."
     endif

     nSets = nSetsMaxUser
   endif
!  At the end of this logic, nSets represents a maximum possible number of phase-space sets.  We'll further constrain it.
!!!!!!!

!!!!!!!
!  Now, let's figure out the decomposition in angle.
!  Decomposition in angle minimizes run time and gives the best thread scaling,
!  so we want the highest number of angle sets possible.
    nAngleSubsets = QuadSet%maxAngleSets

    ! default to using nAngleSubsets if option isn't set or if something nonsensical is passed in
    maxAngleSets  = theDatastore%fetchIfExists("options/global/max_num_anglesets", nAngleSubsets)
    if (maxAngleSets < 1) then
      maxAngleSets = nAngleSubsets
    else if (maxAngleSets > nAngleSubsets) then
      maxAngleSets = nAngleSubsets
    endif
    maxAngleSets = min(maxAngleSets, nSets)

    ! Start with most decomposed option (most angle sets/fewest angles per set),
    ! and iterate until we find a decomposition that
    !    1. evenly divides the angles, and
    !    2. is less than the maximum allowed number of angle sets
    nAngleSets = 1
    do angleSubsetsPerSet=1,nAngleSubsets
       if (mod(nAngleSubsets,angleSubsetsPerSet) == 0) then
          ! If angleSubsetsPerSet evenly divides nAngleSubsets
          nAngleSets = nAngleSubsets/angleSubsetsPerSet
          if (nAngleSets <= maxAngleSets) then
            exit
          endif
       endif
    enddo
!  Note: If nAngleSubsets < nSets, then angleSubsetsPerSet = 1 and nAngleSets = nAngleSubsets
!        If nSets == 1, then nAngleSets = 1, angleSubsetsPerSet = nAngleSubsets, and anglesPerSet = nAngles
!!!!!!!

!!!!!!!
!  Next, the decomposition in groups.
! There are two constraints on the number of group sets:
!  1.  nAngleSets*nGroupSets <= nSets
!  2.  groups per set >= minGroupSetSize, which translates to nGroupSets <= nGroups/minGroupSetSize
    maxGroupSets = min( int(nSets/nAngleSets), int(QuadSet% Groups/minGroupSetSize) )
    maxGroupSets = max( maxGroupSets, 1 ) ! This is needed if QuadSet%Groups < minGroupSetSize

! Iterate through the possibilities for group set decomposition, stop at the first
!   integer that evenly divides the groups
    do nGroupSets=maxGroupSets,1,-1
       if (mod(QuadSet%Groups,nGroupSets) == 0) then
         exit
       endif
    enddo
! After this loop nGroupSets will be set to the value we want.
!!!!!!!

! nSets is now the number of phase space sets we will actually use!
   nSets = nGroupSets*nAngleSets
   anglesPerSet = int(QuadSet%NumAngles/nAngleSets)
   groupsPerSet = int(QuadSet%Groups/nGroupSets)

   extraAngles  = QuadSet%NumAngles - nAngleSets*anglesPerSet
   TETON_VERIFY(extraAngles == 0, 'At this time, Teton does not support having a different number of angles in each phase-space set.')
   extraGroups  = QuadSet% Groups - nGroupSets*groupsPerSet
   TETON_VERIFY(extraGroups == 0, 'At this time, Teton does not support having a different number of groups in each phase-space set.')

   allocate( setGroups(nSets) )
   allocate( setAngles(nSets) )
   allocate( setGroup0(nSets) )
   allocate( setGroupID(nSets) )
   allocate( setAngle0(nSets) )

   totalSets = 0

   ! Create multiple phase-space sets

   nSets = nAngleSets*nGroupSets

!  The following code block handles the case where the groups sets are
!  unbalanced (i.e. not all group sets contain the same number of groups).

!  If there are extra groups assign one more to the first
!  "extraGroups" group sets

   totalSets = 0
   angle0    = 0

   do angleSetID=1,nAngleSets

     g0 = 0
     do groupSetID=1,nGroupSets
       setGroups(totalSets+groupSetID)  = groupsPerSet
       setGroup0(totalSets+groupSetID)  = g0
       setGroupID(totalSets+groupSetID) = groupSetID
       setAngles(totalSets+groupSetID)  = anglesPerSet
       setAngle0(totalSets+groupSetID)  = angle0
       g0                               = g0 + groupsPerSet
     enddo

     totalSets = totalSets + nGroupSets
     angle0    = angle0    + anglesPerSet

     anglesInSet = 0
     angleSubset0 = angleSubsetsPerSet*(angleSetID-1)
     do angleSubset=1,angleSubsetsPerSet
       anglesInSet = AnglesInSet + QuadSet% angleSetSize(angleSubset+angleSubset0)
     enddo
     TETON_VERIFY(anglesInSet == anglesPerSet, "Consistency check between ConstructPhaseSpaceSets.F90 and decomposeAngleSets.F90")
   enddo

!  One more sanity check:
   if (totalSets /= nSets) then
     call f90fatal("ConstructPhaseSpaceSets: totalSets /= nSets")
   endif

! Some sanity checks for the simple one-set case:
   if ( nSets == 1) then
     TETON_VERIFY(setGroups(1) == QuadSet% Groups,    'Error: nSets == 1, but nGroupsPerSet != nGroups')
     TETON_VERIFY(setAngles(1) == QuadSet% numAngles, 'Error: nSets == 1, but nAnglesPerSet != nAngles')
     TETON_VERIFY(setGroup0(1) == 0,                  'Error: nSets == 1, but Set%g0 != 0')
     TETON_VERIFY(setGroupID(1) == 1,                 'Error: nSets == 1, but Set%groupID != 1')
     TETON_VERIFY(setAngle0(1) == 0,                  'Error: nSets == 1, but Set%angle0 != 0')
   endif

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
       if ( setAngles(setID) /= setAngles(setID-1) ) then
         call f90fatal("Teton: Unable to evenly distribute angles across phase-space sets. This is currently a requirement.  Contact the Teton team for tips on adjusting your angle setup to allow even distribution over the phase-space sets.")
       endif
       if ( setGroups(setID) /= setGroups(setID-1) ) then
         call f90fatal("Teton: Unable to evenly distribute groups across phase-space sets. This is currently a requirement.  Contact the Teton team for tips on adjusting your angle setup to allow even distribution over the phase-space sets.")
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

!  Determine the number of "hyper-domains" to increase parallelism
!  in the high-order sweeps and "new" GTA. We need to be
!  careful for very small zone counts so we estimate a
!  maximum number based on the number of zones. We also limit
!  the maximum # based on performance observations. This value is set
!  in the Options mod based on cmake/GetGPUInfo.cmake.  PFN 08/28/2024 ACB 04/15/2025

!  Note that the use of "hyper-domains" will be deprecated once
!  we support sub-meshes per MPI rank.  Also, hyper-domains are
!  not used on the CPU.   PFN 02/14/2023

   nSweepHypDomMax = Options%getSweepMaxHyperDomains()
   nSweepHypDomMax = min (nSweepHypDomMax, int( sqrt( real(Size%nzones) )/2 ))
   nSweepHypDomMax = max( nSweepHypDomMax, 1 )

   nGTAHypDomMax = Options%getGTAMaxHyperDomains()
   nGTAHypDomMax = min (nGTAHypDomMax, int( sqrt( real(Size%nzones) )/2 ))
   nGTAHypDomMax = max( nGTAHypDomMax, 1 )

   if (Size% useGPU) then

     omp_device_num_processors = Options%getNumDeviceProcessors()

!    High-order sweep
     nHypDomMin = int( min(omp_device_num_processors,nSetsMaxUser)/max(nSets,1) )
     nHypDomMin = max( nHypDomMin, 1 )

     Quad% nHyperDomains(1) = min(nSweepHypDomMax, nHypDomMin)

!    GTA Sweep
     nHypDomMin = int( min(omp_device_num_processors,nSetsMaxUser)/max(nGTASets,1) )
     nHypDomMin = max( nHypDomMin, 1 )

     Quad% nHyperDomains(2) = min(nGTAHypDomMax, nHypDomMin)

   else
     Quad% nHyperDomains(1) = 1
     Quad% nHyperDomains(2) = 1
   endif

!  Construct the phase-space sets

   GTASet        = .FALSE.
   nHyperDomains = Quad% nHyperDomains(1)

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
     call Set%construct(setID, groupSetID, angleSetID, QuadID,        &
                        Groups, NumAngles, g0, angle0, nZones, nCorner,  &
                        nHyperDomains, QuadSet, GTASet, fromRestart)

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
       call MPI_COMM_DUP(MY_COMM_GROUP, new_comm, ierror)

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

     GTASet        = .TRUE.
     nHyperDomains = Quad% nHyperDomains(2)
     angle0        =  0
     angleSetID    =  nAngleSets
     groupSetID    =  1

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

       QuadID        = 2 
       Groups        = 1 
       NumAngles     = QuadSet% angleSetSize(setID) 
       g0            = 0
       nZones        = Size% nZones
       nCorner       = Size% ncornr

       if (verbose) then
         write(stdout,100) setID,QuadID,NumAngles,angle0+1,angle0+NumAngles,Groups,g0+1,g0+Groups
       endif

!      construct the GTA set
       call Set%construct(setID, groupSetID, angleSetID, QuadID,        &
                          Groups, NumAngles, g0, angle0, nZones, nCorner,  &
                          nHyperDomains, QuadSet, GTASet, fromRestart)

!      construct an angle set for every GTA set
       call construct(ASet, NumAngles, angle0, nZones,  &
                      nReflecting, GTASet, QuadSet)

!      construct a communication set for every GTA set

!      duplicate the existing communicator
       call MPI_COMM_DUP(MY_COMM_GROUP, new_comm, ierror)

       cSet1 = nSets + setID 
       cSet2 = nSets + setID 

       call construct(CSet, NumAngles, cSet1, cSet2, Groups, new_comm, ASet)

       angle0    = angle0 + NumAngles

     enddo

     if (verbose) then
       write(stdout, 300)
       write(stdout, 200) Size% myRankInGroup,Quad% nHyperDomains(1),Quad% nHyperDomains(2)
       write(stdout, 300)
     endif

!    Construct and incident test on shared boundaries

     call initFindExit(nAngleSets, nGTASets)

   endif

 100 format("       Phase-Angle Set ID =",i3,2x," | Quadrature Set ID =",i2,2x, " | # Angles = ",i3," | Angle IDs =",i3," -",i3, " | # Groups =",i3," | Group IDs = ",i3," -",i3)

 200 format( "hyper-domains for rank = ",i4,": high-order = ",i4,", GTA = ",i4)
 300 format(" ")

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


!  Release memory

   deallocate( setGroups )
   deallocate( setAngles )
   deallocate( setGroup0 )
   deallocate( setGroupID )
   deallocate( setAngle0 )

   return
   end subroutine ConstructPhaseSpaceSets

