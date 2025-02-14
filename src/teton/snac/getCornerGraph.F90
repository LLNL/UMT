!***********************************************************************
!                        Last Update:  04/2019, PFN                    *
!                                                                      *
!   getCornerGraph - This routine builds the sweep ordering array      *
!                    (a directed graph) of CORNERS for a single        *
!                    direction.                                        *
!                                                                      *
!***********************************************************************
   subroutine getCornerGraph(aSetID, angle, nHyperDomains) 

   use kind_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: aSetID
   integer,    intent(in)    :: angle 
   integer,    intent(in)    :: nHyperDomains

!  Local Variables

   type(AngleSet), pointer   :: ASet

   integer    :: i
   integer    :: c
   integer    :: cExit
   integer    :: newCorners
   integer    :: addedCorners
   integer    :: lastCorner
   integer    :: nextCorner
   integer    :: meshCycles
   integer    :: ncornr
   integer    :: nzones

   integer    :: cID
   integer    :: ndone
   integer    :: nHyperPlanes
   integer    :: nZoneSets

   real(adqt) :: omega(Size%ndim)

!  Dynamic

   integer,          allocatable :: need(:)
   integer,          allocatable :: cornerList(:)
   integer,          allocatable :: cornersInPlane(:)
   integer,          allocatable :: cycleList(:)
   integer,          allocatable :: CToHypPlane(:)
   integer,          allocatable :: nDSC(:)
   integer,          allocatable :: DSC(:,:)

   logical (kind=1), allocatable :: done(:)
   logical (kind=1), allocatable :: onCycleList(:)

!  Constants

   ASet      => getAngleSetData(Quad, aSetID) 

   nZoneSets = getNumberOfZoneSets(Quad)
   nzones    =  Size% nzones
   ncornr    =  Size% ncornr
   omega(:)  =  ASet% Omega(:,angle)

!  Allocate arrays

   allocate( need(ncornr) )
   allocate( cornerList(ncornr) )
   allocate( cornersInPlane(ncornr) )
   allocate( cycleList(ncornr) )
   allocate( CToHypPlane(ncornr) )
   allocate( nDSC(ncornr) )
   allocate( DSC(2*Size%maxcf,ncornr) )
   allocate( done(ncornr+1) )
   allocate( onCycleList(ncornr) )

   done(:)        = .FALSE.
   onCycleList(:) = .FALSE.
   done(ncornr+1) = .TRUE. 
   meshCycles     = 0

!  Build NEED array by computing Outward_Normal dot Omega(m)

   call getCornerDependency(omega, need, nDSC, DSC)

!  Create a list of corners to start the sweep

   call getNewCorners(newCorners, meshCycles, need,   &
                      cornerList, cycleList, nDSC, DSC) 

!  Create the "next" array. 

   ndone        = 0
   nextCorner   = 0
   lastCorner   = 0
   nHyperPlanes = 0

   OuterIteration: do

!  Advance to a new hyper-plane

     nHyperPlanes                 = nHyperPlanes + 1
     cornersInPlane(nHyperPlanes) = newCorners
     nextCorner                   = lastCorner + newCorners 
     addedCorners                 = 0

!    Loop over all corners in the current list

     CornerLoop: do cID=1,newCorners

       c       =  cornerList(lastCorner+cID)
       ndone   =  ndone + 1
       done(c) = .TRUE.

!  Loop over the down-stream corners for the corner just added
!  to the next list, decrementing the need array for these
!  neighboring corners 

       do i=1,nDSC(c)
         cExit = DSC(i,c)

         if ( .not. done(cExit) ) then
           need(cExit) = need(cExit) - 1

           if ( need(cExit) == 0 ) then
             nextCorner             = nextCorner   + 1
             addedCorners           = addedCorners + 1
             cornerList(nextCorner) = cExit
           elseif ( need(cExit) < 0 ) then
             write(6,100) Size% myRankInGroup,angle,cExit,c
             call f90fatal("need < 0 in getCornerGraph")
           endif

         endif
       enddo

       ASet% nextC(ndone,angle) = c
       CToHypPlane(c)           = nHyperPlanes

     enddo CornerLoop

     lastCorner = lastCorner + newCorners 

     if (lastCorner == ncornr) then

       exit OuterIteration

     else

       if (addedCorners > 0) then

         newCorners = addedCorners

       elseif (addedCorners == 0) then

!        Break a cycle to add a corner to the list

         call cycleBreakerCorner(ndone, meshCycles, nextCorner,  &
                                 addedCorners, need, cornerList, &
                                 cycleList, nDSC, DSC, onCycleList) 

         newCorners = addedCorners

       endif

       cycle OuterIteration

     endif

   enddo OuterIteration

!  End of Outer Loop, save the number of hyperplanes

   if (meshCycles > Size% ncornr) then
     call f90fatal("MeshCycles exceeds the number of corners in getCornerGraph") 
   endif

   ASet% numCycles(angle) = meshCycles

   call constructHyperPlane( ASet, angle, nHyperPlanes, meshCycles, &
                             nHyperDomains, nZoneSets,              &
                             cornersInPlane(1:nHyperPlanes),        &
                             CToHypPlane, cycleList(1:meshCycles) )

!  Set the number of hyperplanes in the set module for this angle

   ASet% nHyperPlanes(angle) = nHyperPlanes

!  Final error check

   if (ndone /= ncornr) then
     call f90fatal("Wrong number of corners in getCornerGraph")
   endif

 100  format("On Process ",i7," angle = ",i4," cExit = ",i9,  & 
             " has already been done and is down stream of ", &
             " c = ",i9," in getCornerGraph")


!  Release memory

   deallocate( need )
   deallocate( cornerList )
   deallocate( cornersInPlane )
   deallocate( cycleList )
   deallocate( CToHypPlane )
   deallocate( nDSC )
   deallocate( DSC )
   deallocate( done )
   deallocate( onCycleList )

 
   return
   end subroutine getCornerGraph 

