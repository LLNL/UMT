#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Version 1:  11/98, PFN                        *
!                                                                      *
!   INITCOMM - Initializes persistent communication objects that are   *
!              used to exchange boundary fluxes.                       *
!                                                                      *
!   Input:   psibsend - dedicated send buffer                          *
!            psibrecv - dedicated receive buffer                       *
!                                                                      *
!   Output:  IREQUEST - communication "handles"                        *
!                                                                      *
!***********************************************************************

   subroutine initcomm(cSetID)

   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use QuadratureList_mod
   use Communicator_mod
   use BoundaryList_mod
   use Boundary_mod
   use CommSet_mod
   use Options_mod
   use, intrinsic :: iso_c_binding, only : c_bool

   implicit none

!  Arguments

   integer,               intent(in) :: cSetID

!  Local

   type(CommSet),         pointer    :: CSet
   type(Communicator),    pointer    :: CommT
   type(CommunicateFlux), pointer    :: CommX
   type(Boundary),        pointer    :: BdyT

   integer :: sharedID
   integer :: nShared
   integer :: neighbor
   integer :: tag
   integer :: tagBase
   integer :: offset
   integer :: nrecv
   integer :: nsend

   integer :: COMM_GROUP
   integer :: myRank

   integer :: Angle
   integer :: NumAngles
   integer :: numSets
   integer :: Groups

!  Wait for all nodes to arrive

   CSet       => getCommSetData(Quad, cSetID)
   COMM_GROUP =  getCommunicationGroup(CSet)
   myRank     =  Size% myRankInGroup
   numSets    =  CSet% set2 - CSet% set1 + 1

   nShared    =  getNumberOfShared(RadBoundary)
   tagBase    =  10000*cSetID

   call MPIBarrier(COMM_GROUP)

!  Each message is attached to an angle 

   Groups    =  CSet% Groups 
   NumAngles =  CSet% NumAngles

   AngleLoop: do Angle=1,NumAngles


!  Loop over the number of shared surfaces

     CommunicatorLoop: do sharedID=1,nShared

       BdyT     => getShared(RadBoundary, sharedID)
       neighbor =  getNeighborID(BdyT)
       CommT    => getMessage(CSet, sharedID, Angle)

!  Initialize Persistant Communication Objects (receives are "even"
!  requests and sends are "odd" requests)
       nsend  = CommT% nSend*Groups*numSets 
       nrecv  = CommT% nRecv*Groups*numSets
       offset = (myRank + neighbor)*NumAngles
       tag    = tagBase + offset + Angle

       if (nrecv > 0) then
         call MPIRecvInit(CommT%psibrecv, nrecv, neighbor, tag,   &
                          COMM_GROUP, CommT% irequest(2)) 
       endif

       if (nsend > 0) then
         call MPISendInit(CommT%psibsend, nsend, neighbor, tag,  & 
                          COMM_GROUP, CommT% irequest(1)) 
       endif

     enddo CommunicatorLoop

   enddo AngleLoop

!  Initialize Persistant Communication Objects for shared boundary fluxes

   tagBase = 200*cSetID

   do sharedID=1,nShared
     BdyT     => getShared(RadBoundary, sharedID)
     neighbor =  getNeighborID(BdyT)
     CommX    => getSharedFlux(CSet, sharedID)

     nsend    =  CommX% nAngleBins
     nrecv    =  CommX% nAngleBins
     offset   =  (myRank + neighbor)*CommX% nAngleBins
     tag      =  tagBase + offset 

     call MPIRecvInit1(CommX% IncFlux, nrecv, neighbor, tag,   &
                       COMM_GROUP, CommX% irequest(2))

     call MPISendInit1(CommX% ExitFlux, nsend, neighbor, tag,  &
                       COMM_GROUP, CommX% irequest(1))
   enddo


   return
   end subroutine initcomm

