!***********************************************************************
!                        Last Update:  02/2012, PFN                    *
!                                                                      *
!   RTORDER - This routine builds an ordered list of corners for each  *
!             unique direction.                                        *
!                                                                      *
!***********************************************************************
   subroutine rtorder(aSetID) 

   use kind_mod
   use Size_mod
   use constant_mod
   use QuadratureList_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,        intent(in) :: aSetID

!  Local Variables

   type(AngleSet), pointer    :: ASet
   type(HypPlane), pointer    :: HypPlanePtr

   integer                    :: angle
   integer                    :: offSet
   integer                    :: mCycle
   integer                    :: nDomains
   integer                    :: nAngleSets
   integer                    :: nGroupSets
   integer                    :: setID 
   integer                    :: gSetID
   integer                    :: maxZonesPerPlane

!  Dynamic

   integer,  allocatable      :: badCornerList(:)

!  Constants
   ASet       => getAngleSetData(Quad, aSetID)
   nAngleSets =  getNumberOfAngleSets(Quad)
   nGroupSets =  getNumberOfGroupSets(Quad)
   nDomains   =  Quad% nHyperDomains

!  Determine the sweep order for each angle (i.e. the order in which the 
!  zones are solved: "nextZ") 

   ASet% numCycles(:)    = 0
   ASet% cycleOffSet(:)  = 0
   ASet% nHyperPlanes(:) = 0

   AngleLoop: do angle=1,ASet% NumAngles

     if ( .not. ASet% FinishingDirection(angle) ) then
       call snnext(aSetID, angle, nDomains) 
     endif

   enddo AngleLoop

   ASet% totalCycles = ASet% numCycles(1) 

   do angle=2,ASet% NumAngles
     ASet% cycleOffSet(angle) = ASet% cycleOffSet(angle-1) +  &
                                ASet% numCycles(angle-1)
     ASet% totalCycles        = ASet% totalCycles +           &
                                ASet% numCycles(angle)
   enddo

!  Construct cycle List

   allocate( badCornerList(ASet% totalCycles) )

   offSet           = 0
   maxZonesPerPlane = 1

   do angle=1,ASet% NumAngles
     HypPlanePtr => ASet% HypPlanePtr(angle)

     do mCycle=1,ASet% numCycles(angle)
       badCornerList(offSet+mCycle) = HypPlanePtr% badCornerList(mCycle)
     enddo

     offSet           = offSet + ASet% numCycles(angle)

     if ( .not. ASet% FinishingDirection(angle) ) then
        maxZonesPerPlane = max(maxZonesPerPlane, HypPlanePtr%maxZones)
     endif

   enddo

   call constructCycleList(ASet, badCornerList)

!  Allocate dynamic memory that can change size each cycle

   if (aSetID <= nAngleSets) then

     offSet = (aSetID - 1)*nGroupSets

     do gSetID=1,nGroupSets
       setID = offSet + gSetID

       call constructDynMemory(setID, maxZonesPerPlane)
     enddo

   endif
   deallocate( badCornerList )

   return
   end subroutine rtorder


