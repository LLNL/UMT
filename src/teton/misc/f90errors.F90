!-----------------------------------------------------------------------
! f90fatal2
!   This routine emits a message that a fatal error has occurred and
!   shuts down the code.
!
! routine   routine at which error was encountered (provided by cpp)
! line      line number at which error was encountered (provided by cpp)
! message   error message printed
!-----------------------------------------------------------------------
subroutine f90fatal2(routine,line,message)

   use io_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod

!  variable declarations
   implicit none

!  passed variables
   character(*), intent(in) :: routine, message
   integer,      intent(in) :: line

   print '(a,a,a,i6,a,a,a,i6)', &
         "Teton fatal error: '", message, "' on rank: ", &
         Size%myRankInGroup, " in file: ", routine, " at line: ", line

      call MPIAbort(MY_COMM_GROUP)

   return

end subroutine f90fatal2

!-----------------------------------------------------------------------
! f90fatal - older version.
!   Try to use the newer f90fatal2 version.  This older version does
!   not support providing the source file and line where the error
!   occurred.
!-----------------------------------------------------------------------
subroutine f90fatal(message)

   use io_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod

!  variable declarations
   implicit none

!  passed variables
   character(*) :: message

!-----------------------------------------------------------------------

!  issue the fatal error message and exit
   write(nout,500) Size%myRankInGroup, message

   call MPIAbort(MY_COMM_GROUP)

!-----------------------------------------------------------------------
!  format statements
!-----------------------------------------------------------------------
500 format(/1x,"The following fatal error has occurred on process: ", &
           i6,/3x,a)
!-----------------------------------------------------------------------

   return
end subroutine f90fatal

!-----------------------------------------------------------------------
! f90verify This routine performs verification checking.
!
! bool      scalar logical (boolean)
! routine   routine in which verification is written (provided by cpp)
! line      line number at which verificion is written (provided by cpp)
! message   error message printed upon failed verification
!-----------------------------------------------------------------------
subroutine f90verify(bool,routine,line,message)

   use io_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod

!  variable declarations
   implicit none

!  passed variables
   logical,      intent(in) :: bool
   character(*), intent(in) :: routine, message
   integer,      intent(in) :: line

   if ( .not. bool ) then
      print '(a,a,a,i6,a,a,a,i6)', &
         "Teton verification failed: '", message, "' on rank: ", &
         Size%myRankInGroup, " in file: ", routine, " at line: ", line

      call MPIAbort(MY_COMM_GROUP)
   endif

   return

end subroutine f90verify



!-----------------------------------------------------------------------
! f90assert This routine performs assertion checking.
!
! bool      scalar logical (boolean) assertion
! routine   routine in which assertion is written (provided by cpp)
! line      line number at which assertion is written (provided by cpp)
! message   error message printed upon failed assertion
!-----------------------------------------------------------------------
subroutine f90assert(bool,routine,line,message)
   use io_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod

!  variable declarations
   implicit none

!  passed variables
   logical,      intent(in) :: bool
   character(*), intent(in) :: routine, message
   integer,      intent(in) :: line

   if ( .not. bool ) then
      print '(a,a,a,i6,a,a,a,i6)', &
         "Teton assertion failed: '", message, "' on rank: ", &
         Size%myRankInGroup, " in file: ", routine, " at line: ", line

      call MPIAbort(MY_COMM_GROUP)
   endif

   return

end subroutine f90assert
