!***********************************************************************
!                        Last Update:  04/2024, PFN                    *
!                                                                      *
!   getDirectedGraph - This routine builds an ordered list of corners  *
!                      or zones (depending on sweeping method) for     *
!                      each unique direction.                          *
!                                                                      *
!***********************************************************************
   subroutine getDirectedGraph

   use kind_mod
   use Size_mod
   use constant_mod
   use QuadratureList_mod
   use AngleSet_mod
   use Options_mod

   implicit none

!  Local Variables

   type(AngleSet), pointer    :: ASet
   type(HypPlane), pointer    :: HypPlanePtr

   integer                    :: aSetID
   integer                    :: angle
   integer                    :: n
   integer                    :: offSet
   integer                    :: mCycle
   integer                    :: nHyperDomains
   integer                    :: nAngleSets
   integer                    :: nGroupSets
   integer                    :: nGTASets
   integer                    :: setID 
   integer                    :: gSetID
   integer                    :: maxPerPlane
   integer                    :: numAngles
   integer                    :: sweepVersion 
   integer                    :: totalAngles

!  Dynamic

   integer,  allocatable      :: badCornerList(:)
   integer,  allocatable      :: angleList(:,:)
   integer,  allocatable      :: elementsPerPlane(:)

!  Constants

   nAngleSets    =  getNumberOfAngleSets(Quad)
   nGroupSets    =  getNumberOfGroupSets(Quad)
   nGTASets      =  getNumberOfGTASets(Quad)
   sweepVersion  =  Options% getSweepVersion()

!  Create a list of angleSet/angle number pairs

   totalAngles = 0

   do aSetID=1,nAngleSets+nGTASets
     ASet => getAngleSetData(Quad, aSetID)
     totalAngles = totalAngles + ASet% NumAngles
   enddo

   allocate( angleList(2,totalAngles) )
   allocate( elementsPerPlane(totalAngles) )

   elementsPerPlane(:) = 0

   n = 0
   do aSetID=1,nAngleSets+nGTASets
     ASet => getAngleSetData(Quad, aSetID)

     ASet% numCycles(:)    = 0
     ASet% cycleOffSet(:)  = 0
     ASet% nHyperPlanes(:) = 0

     do angle=1,ASet% NumAngles
       n              = n + 1
       angleList(1,n) = aSetID
       angleList(2,n) = angle
     enddo
   enddo

!  Determine the sweep order for each angle (i.e. the order in which the 
!  zones or corners are solved: "next") 
   !$omp parallel do default(none) schedule(static) &
   !$omp& shared(totalAngles,sweepVersion,angleList,elementsPerPlane,Quad) &
   !$omp& private(ASet,HypPlanePtr,aSetID,angle,nHyperDomains)
   AngleListLoop: do n=1,totalAngles

     aSetID = angleList(1,n)
     angle  = angleList(2,n)

     ASet   => getAngleSetData(Quad, aSetID)

     if ( .not. ASet% FinishingDirection(angle) ) then

       HypPlanePtr => ASet% HypPlanePtr(angle)

       if ( ASet% GTASet ) then
         nHyperDomains =  getNumberOfHyperDomains(Quad,2)

         call getZoneGraph(aSetID, angle, nHyperDomains)
       else
         nHyperDomains =  getNumberOfHyperDomains(Quad,1)

         if ( sweepVersion == 1 ) then
           call getZoneGraph(aSetID, angle, nHyperDomains)
           elementsPerPlane(n) = HypPlanePtr% maxZones
         else if ( sweepVersion == 2 ) then
           call getCornerGraph(aSetID, angle, nHyperDomains)
           elementsPerPlane(n) = HypPlanePtr% maxCorners
         endif

       endif
           
     endif

   enddo AngleListLoop

   !$omp end parallel do

   maxPerPlane = maxval( elementsPerPlane(1:totalAngles) )
   maxPerPlane = max( maxPerPlane, 1 )

   AngleSetLoop: do aSetID=1,nAngleSets+nGTASets
     ASet      => getAngleSetData(Quad, aSetID)
     numAngles =  ASet% NumAngles

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
     allocate( badCornerList(ASet% totalCycles+1 ) )

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

   enddo AngleSetLoop

!  Find the maximum number of hyper-elements (high-order and GTA)

   do aSetID=1,nAngleSets
     Quad% nHyperElements(1) = max( Quad% nHyperElements(1), Quad% AngSetPtr(aSetID)% maxInterface )
   enddo

   do aSetID=nAngleSets+1,nAngleSets+nGTASets
     Quad% nHyperElements(2) = max( Quad% nHyperElements(2), Quad% AngSetPtr(aSetID)% maxInterface )
   enddo



   deallocate( angleList )
   deallocate( elementsPerPlane )


   return
   end subroutine getDirectedGraph 


