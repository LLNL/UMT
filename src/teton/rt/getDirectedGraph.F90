!***********************************************************************
!                        Last Update:  04/2024, PFN                    *
!                                                                      *
!   getDirectedGraph - This routine builds an ordered list of corners  *
!                      or zones (depending on sweeping method) for     *
!                      each unique direction.                          *
!                                                                      *
!***********************************************************************
   subroutine getDirectedGraph(aSetID)

   use kind_mod
   use Size_mod
   use constant_mod
   use QuadratureList_mod
   use AngleSet_mod
   use Options_mod

   implicit none

!  Arguments

   integer,        intent(in) :: aSetID

!  Local Variables

   type(AngleSet), pointer    :: ASet
   type(HypPlane), pointer    :: HypPlanePtr

   integer                    :: angle
   integer                    :: offSet
   integer                    :: mCycle
   integer                    :: nHyperDomains
   integer                    :: nAngleSets
   integer                    :: nGroupSets
   integer                    :: setID 
   integer                    :: gSetID
   integer                    :: maxPerPlane
   integer                    :: numAngles
   integer                    :: sweepVersion 

!  Dynamic

   integer,  allocatable      :: badCornerList(:)

!  Constants
   ASet => getAngleSetData(Quad, aSetID)

   nAngleSets    =  getNumberOfAngleSets(Quad)
   nGroupSets    =  getNumberOfGroupSets(Quad)
   numAngles     =  ASet% NumAngles
   sweepVersion  =  Options% getSweepVersion()

!  Determine the sweep order for each angle (i.e. the order in which the 
!  zones are solved: "nextZ") 

   ASet% numCycles(:)    = 0
   ASet% cycleOffSet(:)  = 0
   ASet% nHyperPlanes(:) = 0
   maxPerPlane           = 1

   AngleLoop: do angle=1,numAngles

     if ( .not. ASet% FinishingDirection(angle) ) then

       HypPlanePtr => ASet% HypPlanePtr(angle)

       if ( ASet% GTASet ) then
         nHyperDomains =  getNumberOfHyperDomains(Quad,2)

         call getZoneGraph(aSetID, angle, nHyperDomains)
       else
         nHyperDomains =  getNumberOfHyperDomains(Quad,1)

         if ( sweepVersion == 0 ) then
           call getZoneGraph(aSetID, angle, nHyperDomains)
           maxPerPlane = max(maxPerPlane, HypPlanePtr% maxZones)
         else
           call getCornerGraph(aSetID, angle, nHyperDomains)
           maxPerPlane = max(maxPerPlane, HypPlanePtr% maxCorners)
         endif

       endif
           
     endif

   enddo AngleLoop

   ASet% totalCycles = ASet% numCycles(1) 

   do angle=2,numAngles
     ASet% cycleOffSet(angle) = ASet% cycleOffSet(angle-1) +  &
                                ASet% numCycles(angle-1)
     ASet% totalCycles        = ASet% totalCycles +           &
                                ASet% numCycles(angle)
   enddo

!  Construct cycle List

   ! Cray is having trouble with zero length arrays.  If this array is
   ! zero length it appears to cause code corruption in the cycle list constructor
   ! below.
   ! See https://rzlc.llnl.gov/gitlab/deterministic-transport/TRT/Teton/-/issues/429
   allocate( badCornerList(ASet% totalCycles +1 ) )

   offSet             = 0
   ASet% maxInterface = 1

   do angle=1,numAngles
     HypPlanePtr => ASet% HypPlanePtr(angle)

     do mCycle=1,ASet% numCycles(angle)
       badCornerList(offSet+mCycle) = HypPlanePtr% badCornerList(mCycle)
     enddo

     offSet = offSet + ASet% numCycles(angle)

     if ( .not. ASet% FinishingDirection(angle) ) then
        ASet% maxInterface = max( ASet% maxInterface , HypPlanePtr% interfaceLen )
     endif

   enddo

   call constructCycleList(ASet, badCornerList)

!  Allocate dynamic memory that can change size each cycle

   if ( .not. ASet% GTASet ) then

     offSet = (aSetID - 1)*nGroupSets

     do gSetID=1,nGroupSets
       setID = offSet + gSetID

       call constructDynMemory(setID, maxPerPlane)
     enddo

   endif


   deallocate( badCornerList )


   return
   end subroutine getDirectedGraph 


