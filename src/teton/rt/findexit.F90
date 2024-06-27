!***********************************************************************
!                        Version 1:  12/98, PFN                        *
!                                                                      *
!   FINDEXIT - Find the angles that are exiting for each boundary      *
!              element on a shared surface.                            *
!                                                                      *
!***********************************************************************

   subroutine findexit(aSetID)

   use kind_mod
   use constant_mod
   use flags_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use Geometry_mod
   use Communicator_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use Boundary_mod
   use CommSet_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,               intent(in) :: aSetID

!  Local

   type(CommSet),         pointer    :: CSet
   type(AngleSet),        pointer    :: ASet
   type(IncidentTest),    pointer    :: TSet
   type(Communicator),    pointer    :: CommT
   type(CommunicateFlux), pointer    :: CommX
   type(Boundary),        pointer    :: BdyT

   integer    :: sharedID
   integer    :: nShared
   integer    :: neighborRank
   integer    :: nSend
   integer    :: nRecv
   integer    :: izero
   integer    :: NumBdyElem
   integer    :: bin
   integer    :: b
   integer    :: b0
   integer    :: BCType

   integer    :: angle 
   integer    :: angle1 
   integer    :: angle2
   integer    :: angleRecv1 
   integer    :: angleRecv2

   integer    :: myRankInGroup 
   integer    :: Groups
   integer    :: NumAngles
   integer    :: NumBin0
   integer    :: offset

   integer    :: nGroupSets
   integer    :: nAngleSets

   integer    :: COMM_GROUP

   integer    :: n
   integer    :: nxBdy
   integer    :: nBoundary

   real(adqt) :: dot

!  Dynamic

   integer,  allocatable :: bdyList(:,:)
   integer,  allocatable :: ListSend(:,:)
   integer,  allocatable :: ListRecv(:) 

!  Constants

   parameter (izero=0)

!  There is a one-to-one correspondence between angle sets
!  and communication sets 

   ASet       => getAngleSetData(Quad, aSetID)
   CSet       => getCommSetData(Quad, aSetID)

   NumAngles  =  ASet% NumAngles
   nShared    =  getNumberOfShared(RadBoundary)
   nAngleSets =  getNumberOfAngleSets(Quad)


   if (aSetID <= nAngleSets) then
     nGroupSets = getNumberOfGroupSets(Quad)
   else
     nGroupSets = 1
   endif

   DecompTest: if (nShared > 0) then

     myRankInGroup =  Size% myRankInGroup 
     COMM_GROUP    =  CSet% COMM_GROUP
     NumBin0       =  ASet% NumBin0

!    Wait for all nodes to arrive

     call MPIBarrier(COMM_GROUP)

!  Loop over all shared boundary elements and decide  
!  which angles are exiting and which are incident on 
!  shared surfaces (we only need to check the unique 
!  angle set).  Each process on a shared surface computes
!  dot products for half of the angles and then
!  the results are exchanged.

     PostReceiveLoop: do sharedID=1,nShared

!  Start receives

       TSet => getIncidentTest(ASet, sharedID)

       TSet% IncTestR(:) = 0
       TSet% IncTest(:)  = 0

       call MPIStart(TSet% request(2)) 

     enddo PostReceiveLoop

!  Start sends

     PostSendLoop: do sharedID=1,nShared

       TSet         => getIncidentTest(ASet, sharedID)
       BdyT         => getShared(RadBoundary, sharedID)
       neighborRank =  getNeighborID(BdyT)
       NumBdyElem   =  getNumberOfBdyElements(BdyT)

       if (myRankInGroup < neighborRank) then
         angle1     = 1
         angle2     = NumAngles/2
       else
         angle1     = NumAngles/2 + 1
         angle2     = NumAngles
       endif

!  Loop over boundary elements

       AngleLoop1: do angle=angle1,angle2

         offset = NumBdyElem*(angle - angle1)

         BoundaryElements1: do b=1,NumBdyElem

           dot = DOT_PRODUCT( ASet%omega(:,angle),BdyT%A_bdy(:,b) ) 

           if (dot < zero) then
             TSet% IncTest(offset+b) = -1
           elseif (dot > zero) then
             TSet% IncTest(offset+b) =  1
           endif

         enddo BoundaryElements1 

       enddo AngleLoop1


       call MPIStart(TSet% request(1)) 

     enddo PostSendLoop 

!  Check that all sends are complete

     do sharedID=1,nShared
       TSet => getIncidentTest(ASet, sharedID)
       call MPIWait(TSet% request(1))
     enddo

!  Now create the lists of incident and exiting boundary elements 
!  for each angle in the quadrature set

     CSet% NetFlux(:,:) = zero

!  Check that all receives are complete

     do sharedID=1,nShared
       TSet => getIncidentTest(ASet, sharedID)
       call MPIWait(TSet% request(2))
     enddo

     ConstructBufferLoop: do sharedID=1,nShared

       TSet         => getIncidentTest(ASet, sharedID)
       BdyT         => getShared(RadBoundary, sharedID)
       NumBdyElem   =  getNumberOfBdyElements(BdyT) 
       neighborRank =  getNeighborID(BdyT)
       b0           =  getFirstBdyElement(BdyT) - 1 

       if (myRankInGroup < neighborRank) then
         angle1     = 1
         angleRecv1 = NumAngles/2 + 1
         angleRecv2 = NumAngles
       else
         angleRecv1 = 1
         angleRecv2 = NumAngles/2
         angle1     = NumAngles/2 + 1
       endif

!  Dot products received from neighbor have the opposite sign

       TSet% IncTestR(:) = -TSet% IncTestR(:)

       allocate( ListSend(2,NumBdyElem) )
       allocate( ListRecv(NumBdyElem) )

       AngleLoop2: do angle=1,NumAngles

         bin   = CSet% AngleToBin(angle)
         nSend = 0
         nRecv = 0

         if (angle >= angleRecv1 .and. angle <= angleRecv2) then

           offset = NumBdyElem*(angle - angleRecv1)

           BoundaryElements2: do b=1,NumBdyElem

             if (TSet% IncTestR(offset+b) < izero) then
               nRecv             = nRecv + 1
               ListRecv(nRecv)   = b0 + b
             elseif (TSet% IncTestR(offset+b) > izero) then
               nSend             = nSend + 1
               ListSend(1,nSend) = b0 + b
               ListSend(2,nSend) = BdyT% BdyToC(b)
             endif

           enddo BoundaryElements2

         else

           offset = NumBdyElem*(angle - angle1)

           BoundaryElements3: do b=1,NumBdyElem

             if (TSet% IncTest(offset+b) < izero) then
               nRecv             = nRecv + 1
               ListRecv(nRecv)   = b0 + b
             elseif (TSet% IncTest(offset+b) > izero) then
               nSend             = nSend + 1
               ListSend(1,nSend) = b0 + b
               ListSend(2,nSend) = BdyT% BdyToC(b)
             endif

           enddo BoundaryElements3

         endif

!  Allocate send/receive buffers for this message

         CommT  => getMessage(CSet, sharedID, angle)
         Groups =  CSet% Groups

         call constructBuffer(CommT, nSend, nRecv, Groups, nGroupSets)

         if (nSend > 0) then
           CommT%ListSend(:,1:nSend) = ListSend(:,1:nSend)
         endif

         if (nRecv > 0) then
           CommT%ListRecv(1:nRecv) = ListRecv(1:nRecv)
         endif

         CSet% NetFlux(sharedID,bin) = CSet% NetFlux(sharedID,bin) + &
                                       real( nRecv*Groups*nGroupSets, adqt )
       enddo AngleLoop2


       CommX => getSharedFlux(CSet, sharedID)

       call constructFluxBuffer(CommX, NumBin0)

       deallocate( ListSend )
       deallocate( ListRecv )

     enddo ConstructBufferLoop 

   endif DecompTest

!  Loop over all boundaries and create a list of exiting boundary elements

   nBoundary =  getNumberOfBoundaries(RadBoundary)
   nxBdy     =  0

   allocate( bdyList(2,Size% nbelem) )

   AngleLoop0: do angle=1,NumAngles

     nxBdy = 0

!    All boundaries except shared

     do n=1,nBoundary
       BdyT         => getBoundary(RadBoundary, n)
       BCType       =  getBCType(BdyT)

       if (BCType /= bcType_shared) then

         NumBdyElem   =  getNumberOfBdyElements(BdyT)
         b0           =  getFirstBdyElement(BdyT) - 1

         do b=1,NumBdyElem
           dot = DOT_PRODUCT( ASet%omega(:,angle),BdyT%A_bdy(:,b) )

           if (dot > zero) then
             nxBdy            = nxBdy + 1
             bdyList(1,nxBdy) = b0 + b
             bdyList(2,nxBdy) = BdyT% BdyToC(b)
           endif
         enddo

       endif

     enddo

!    Shared boundaries

     do sharedID=1,nShared
       CommT => getMessage(CSet, sharedID, angle)

       do n=1,CommT% nSend
         nxBdy            = nxBdy + 1
         bdyList(1,nxBdy) = CommT% ListSend(1,n)
         bdyList(2,nxBdy) = CommT% ListSend(2,n)
       enddo

     enddo

     call constructBdyExitList(ASet, angle, nxBdy, bdyList)

   enddo AngleLoop0

   deallocate( bdyList )

   return
   end subroutine findexit

