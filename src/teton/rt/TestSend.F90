!***********************************************************************
!                        Version 1:  11/98, PFN                        *
!                                                                      *
!   TestSend - Tests for the completion of Sends.                      *
!                                                                      *
!***********************************************************************

   subroutine TestSend(cSetID, sendIndex) 

   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use Communicator_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use CommSet_mod

   implicit none

!  Arguments

   integer, intent(in)         :: cSetID
   integer, intent(in)         :: sendIndex

!  Local

   type(CommSet),      pointer :: CSet
   type(Communicator), pointer :: CommT

   integer                     :: Angle
   integer                     :: sharedID
   integer                     :: nShared

!  Constants

   CSet    => getCommSetData(Quad, cSetID)
   nShared =  getNumberOfShared(RadBoundary)

!  Test for the completion of sends  

   do sharedID=1,nShared
     Angle =  CSet% RecvOrder(sendIndex,sharedID)
     CommT  => getMessage(CSet, sharedID, Angle)

     if (CommT% nSend > 0) then
       call MPIWait(CommT% irequest(1))
     endif
   enddo


   return
   end subroutine TestSend 

