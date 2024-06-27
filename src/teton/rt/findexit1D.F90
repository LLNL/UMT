!***********************************************************************
!                        Version 1:  12/98, PFN                        *
!                                                                      *
!   FINDEXIT - Find the angles that are exiting for each boundary      *
!              element on a shared surface.                            *
!                                                                      *
!***********************************************************************

   subroutine findexit1D(cSetID)

   use kind_mod
   use constant_mod
   use Size_mod
   use Communicator_mod
   use QuadratureList_mod
   use Quadrature_mod
   use BoundaryList_mod
   use Boundary_mod
   use AngleSet_mod
   use CommSet_mod

   implicit none

!  Arguments

   integer, intent(in)             :: cSetID

!  Local

   type(CommSet),         pointer  :: CSet
   type(AngleSet),        pointer  :: ASet
   type(Boundary),        pointer  :: BdyT
   type(Communicator),    pointer  :: CommT
   type(CommunicateFlux), pointer  :: CommX

   integer    :: sharedID
   integer    :: nSend,nRecv
   integer    :: zone
   integer    :: ib, cornerID

   integer    :: angle
   integer    :: angleSend1, angleSend2
   integer    :: angleRecv1, angleRecv2

   integer    :: nShared
   integer    :: nGroupSets
   integer    :: Groups
   integer    :: NumAngles 


   DecompTest: if (Size% ncomm > 0) then

     nShared    =  getNumberOfShared(RadBoundary)
     nGroupSets =  getNumberOfGroupSets(Quad)

     CSet       => getCommSetData(Quad, cSetID)
     ASet       => getAngleSetFromSetID(Quad, cSetID) 

     NumAngles  =  CSet% NumAngles
     Groups     =  CSet% Groups

     cornerID   =  0 
     angleSend1 =  0
     angleSend2 =  0
     angleRecv1 =  0
     angleRecv2 =  0

!  Now create the lists of incident and exiting boundary elements 
!  for each angle in the quadrature set

!  If the zoneID = 1:
!      Incident angles  mu>0   NumAngles/2 + 1 <= angle <= NumAngles
!      Exiting angles   mu<0   1 <= angle <= NumAngles/2
!  If the zoneID = nzones:
!      Incident angles  mu<0   1 <= angle <= NumAngles/2
!      Exiting angles   mu>0   NumAngles/2 + 1 <= angle <= NumAngles


     CommunicatorLoop: do sharedID=1,nShared

       BdyT => getShared(RadBoundary, sharedID)
       zone =  BdyT% BdyToZone(1)
       ib   =  getFirstBdyElement(BdyT)

       if (zone == 1) then
         cornerID   = 1
         angleSend1 = 1
         angleSend2 = NumAngles/2
         angleRecv1 = NumAngles/2 + 1
         angleRecv2 = NumAngles
       elseif (zone == Size% nzones) then
         cornerID   = 2*Size% nzones 
         angleRecv1 = 1
         angleRecv2 = NumAngles/2
         angleSend1 = NumAngles/2 + 1
         angleSend2 = NumAngles
       endif

!  Allocate send and receive buffers

       SendLoop: do angle=angleSend1,angleSend2 

         CommT => getMessage(CSet, sharedID, angle)

         if ( ASet% FinishingDirection(Angle) ) then

           nSend = 0 
           nRecv = 0

           call constructBuffer(CommT, nSend, nRecv, Groups, nGroupSets)

         else

           nSend = 1 
           nRecv = 0

           call constructBuffer(CommT, nSend, nRecv, Groups, nGroupSets)

           CommT% ListSend(1,1) = ib
           CommT% ListSend(2,1) = cornerID 

         endif

       enddo SendLoop

       ReceiveLoop: do angle=angleRecv1,angleRecv2 

         CommT => getMessage(CSet, sharedID, angle)

         if ( ASet% FinishingDirection(Angle) ) then

           nSend = 0
           nRecv = 0 

           call constructBuffer(CommT, nSend, nRecv, Groups, nGroupSets)

         else

           nSend = 0 
           nRecv = 1 

           call constructBuffer(CommT, nSend, nRecv, Groups, nGroupSets)

           CommT% ListRecv(1) = ib

         endif

       enddo ReceiveLoop

       CommX => getSharedFlux(CSet, sharedID)
       call constructFluxBuffer(CommX, NumAngles)

     enddo CommunicatorLoop

   endif DecompTest

   return
   end subroutine findexit1D

