#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Version 1:  11/98, PFN                        *
!                                                                      *
!   SendFlux - Sends boundary fluxes on shared surfaces.               *
!                                                                      *
!***********************************************************************

   subroutine SendFlux(SnSweep, cSetID, sendIndex, PsiB) 

   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use Communicator_mod
   use AngleSet_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use CommSet_mod
   use SetData_mod
   use Options_mod
   use, intrinsic :: iso_c_binding, only : c_bool

   implicit none

!  Arguments

   logical (kind=1),     intent(in) :: SnSweep
   integer,              intent(in) :: cSetID
   integer,              intent(in) :: sendIndex
   real(adqt), optional, intent(in) :: PsiB(Size%nSurfElem,Size%nangGTA)

!  Local

   type(SetData),        pointer    :: Set
   type(CommSet),        pointer    :: CSet
   type(AngleSet),       pointer    :: ASet
   type(Communicator),   pointer    :: CommT

   integer                          :: setID
   integer                          :: n
   integer                          :: angle0
   integer                          :: Angle
   integer                          :: i
   integer                          :: b
   integer                          :: sharedID
   integer                          :: nShared

!  Constants

   nShared =  getNumberOfShared(RadBoundary)
   CSet    => getCommSetData(Quad, cSetID)


!  Send boundary fluxes needed by neighboring domains for this sweep 

   SendLoop: do sharedID=1,nShared

     Angle =  CSet% RecvOrder(sendIndex,sharedID)
     CommT => getMessage(CSet, sharedID, Angle)

     if ( SnSweep) then

       if ( CommT% nSend > 0 ) then

!  Loop over boundary elements that are exiting for this communicator

         n = 0
         do setID=CSet% set1,CSet% set2
           Set => getSetData(Quad, setID)

           do i=1,CommT% nSend
             b                      = CommT% ListSend(1,i)
             CommT% psibsend(:,n+i) = Set% PsiB(:,b,Angle)
           enddo

           n = n + CommT% nSend 
         enddo

!  Start send for this communicator (odd numbered handle)

         call MPIStart(CommT% irequest(1))

       endif

     else

       if ( CommT% nSend > 0 ) then

         ASet   => CSet% angleSetPtr
         angle0 =  ASet% angle0

         do i=1,CommT% nSend
           b                    = CommT% ListSend(1,i)
           CommT% psibsend(1,i) = PsiB(b,angle0+Angle)
         enddo

!  Start send for this communicator (odd numbered handle)

         call MPIStart(CommT% irequest(1))

       endif

     endif

   enddo SendLoop


   return
   end subroutine SendFlux 

