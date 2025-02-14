#include "macros.h"
!***********************************************************************
!                         Version 0: 04/06 PFN                         *
!                                                                      *
!    setCommunicationGroup - set the communication group for this      *
!                            MPI process.                              *
!                                                                      *
!    Updated 07/2017 by maginot1                                       *
!      - Added ISO_C_BINDING to improve Teton library interoperability *
!                                                                      *
!***********************************************************************

   subroutine setCommunicationGroup(comm) &
        BIND(C,NAME="teton_setcommunicationgroup")

   use mpif90_mod 
   use mpi_param_mod
   use ISO_C_BINDING
   
   implicit none

!  Passed Variables
   integer(C_INT), intent(in)    :: comm

#if defined (TETON_ENABLE_OPENMP)
!  Local Variables
   integer                       :: provided ! level of thread support
   integer                       :: errorcode,  ierr     ! error code from MPI
#endif

   MY_COMM_GROUP = comm

! Note - attempts to use TETON_VERIFY macro resulted in crashes.  Calling
! MPI_Abort directly works.
#if defined (TETON_ENABLE_OPENMP)
   call MPI_Query_thread(provided, ierr)
   if (ierr /= MPI_SUCCESS) then
      print *, "Teton: MPI thread support query failed."
      errorcode = 1
      call MPI_Abort(MPI_COMM_WORLD, errorcode, ierr)
   endif
   if (provided /= MPI_THREAD_MULTIPLE) then
      print *, "Teton: MPI was not initialized with thread support level MPI_THREAD_MULTIPLE.  A thread-safe MPI is required when OpenMP is enabled in Teton."
      errorcode = 1
      call MPI_Abort(MPI_COMM_WORLD, errorcode, ierr)
   endif
#endif



   return
   end subroutine setCommunicationGroup 
