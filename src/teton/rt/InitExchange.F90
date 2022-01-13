!***********************************************************************
!                        Version 1:  11/98, PFN                        *
!                                                                      *
!   INITEXCHANGE - Posts receives for boundary fluxes.                 *
!                                                                      *
!   Input:                                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!***********************************************************************

   subroutine InitExchange(cSetID)


   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use Communicator_mod
   use CommSet_mod

   implicit none

!  Arguments

   integer,    intent(in)      :: cSetID

!  Local

   type(CommSet),      pointer :: CSet
   type(Communicator), pointer :: CommT

   integer                     :: angle
   integer                     :: bin
   integer                     :: NumAngles
   integer                     :: sharedID
   integer                     :: nShared
   integer                     :: COMM_GROUP

!  Constants

   nShared    =  getNumberOfShared(RadBoundary)

   if (nShared > 0) then

   CSet       => getCommSetData(Quad, cSetID)
   COMM_GROUP =  getCommunicationGroup(CSet)
   NumAngles  =  CSet% NumAngles

!  Wait for all nodes to arrive

   call MPIBarrier(COMM_GROUP)

!  First, start all of the receives (receives use even numbered handles)

   do angle=1,NumAngles
     bin = CSet% AngleToBin(angle)

     if ( .not. CSet% Converged(bin) ) then

       do sharedID=1,nShared
         CommT => getMessage(CSet, sharedID, angle)
         if (CommT% nRecv > 0) then
           call MPIStart(CommT% irequest(2))
         endif
       enddo

     endif

   enddo


   call MPIBarrier(COMM_GROUP)

   endif


   return
   end subroutine InitExchange

