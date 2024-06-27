#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Version 1:  11/98, PFN                        *
!                                                                      *
!   RecvFlux - Receives boundary fluxes on shared surfaces.            *
!                                                                      *
!***********************************************************************

   subroutine RecvFlux(SnSweep, cSetID, Angle, PsiB) 

   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use AngleSet_mod
   use Communicator_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use CommSet_mod
   use SetData_mod
   use Options_mod
   use, intrinsic :: iso_c_binding, only : c_bool

   implicit none

!  Arguments

   logical (kind=1),     intent(in)    :: SnSweep
   integer,              intent(in)    :: cSetID
   integer,              intent(in)    :: Angle
   real(adqt), optional, intent(inout) :: PsiB(Size%nSurfElem,Size%nangGTA)

!  Local

   type(SetData),        pointer       :: Set
   type(CommSet),        pointer       :: CSet
   type(AngleSet),       pointer       :: ASet
   type(Communicator),   pointer       :: CommT

   integer                             :: setID
   integer                             :: n
   integer                             :: i
   integer                             :: b
   integer                             :: angle0
   integer                             :: sharedID
   integer                             :: nShared

!  Constants

   nShared =  getNumberOfShared(RadBoundary)
   CSet    => getCommSetData(Quad, cSetID)

!  Process data received from neighbors 

   ReceiveLoop: do sharedID=1,nShared

     CommT => getMessage(CSet, sharedID, Angle)

     if ( SnSweep) then

       if (CommT% nRecv > 0) then

!      Check for completion of the receive

         call MPIWait(CommT% irequest(2))

!  Loop over boundary elements that are incident for this communicator

         n = 0
         do setID=CSet% set1,CSet%set2
           Set => getSetData(Quad, setID)

           do i=1,CommT% nRecv
             b                    = CommT% ListRecv(i)
             Set% PsiB(:,b,Angle) = CommT% psibrecv(:,n+i)
           enddo

           n = n + CommT% nRecv
         enddo

       endif

     else

       if (CommT% nRecv > 0) then

         ASet   => CSet% angleSetPtr
         angle0 =  ASet% angle0

!      Check for completion of the receive

         call MPIWait(CommT% irequest(2))

!  Loop over boundary elements that are incident for this communicator

         do i=1,CommT% nRecv
           b                    = CommT% ListRecv(i)
           PsiB(b,angle0+Angle) = CommT% psibrecv(1,i)
         enddo

       endif

     endif

!  End loop over shared surfaces

   enddo ReceiveLoop


   return
   end subroutine RecvFlux 

