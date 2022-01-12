!*******************************************************************************
!  Array Bounds Checking - A set of subroutines that help detect array         *
!     operation errors, such as out-of-bounds accesses.                        *
!*******************************************************************************

module ArrayChecks_mod
   use, intrinsic :: iso_c_binding, only : c_double, c_int
   implicit none

   interface is_legal_access
      module procedure is_legal_access_double_array_1d, &
                       is_legal_access_double_array_2d, &
                       is_legal_access_double_array_3d, &
                       is_legal_access_int_array_1d, &
                       is_legal_access_int_array_2d, &
                       is_legal_access_int_array_3d, &
                       is_legal_access_AngleSet_array_1d, &
                       is_legal_access_SetData_array_1d
                     
   end interface

contains

   !--------------------------------------------------------------
   ! Checks that an array access is within the bounds of the array
   !--------------------------------------------------------------
   function is_legal_access_double_array_1d(ptr, i) result(isLegal)

      real(kind=c_double), intent(in), pointer, contiguous, dimension(:) :: ptr
      integer :: i
      logical :: isLegal

      if (i > SIZE(ptr, 1)) then
         isLegal = .FALSE.
      else
         isLegal = .TRUE.
      endif
         
      return
   end function

   function is_legal_access_double_array_2d(ptr, i, j) result(isLegal)

      real(kind=c_double), intent(in), pointer, contiguous, dimension(:,:) :: ptr
      integer :: i, j
      logical :: isLegal

      if (i > SIZE(ptr, 1) .OR. j > SIZE(ptr, 2) ) then
         isLegal = .FALSE.
      else
         isLegal = .TRUE.
      endif
         
      return
   end function

   function is_legal_access_double_array_3d(ptr, i, j, k) result(isLegal)

      real(kind=c_double), intent(in), pointer, contiguous, dimension(:,:,:) :: ptr
      integer :: i, j, k
      logical :: isLegal

      if (i > SIZE(ptr, 1) .OR. j > SIZE(ptr, 2) .OR. k > SIZE(ptr, 3) ) then
         isLegal = .FALSE.
      else
         isLegal = .TRUE.
      endif
         
      return
   end function

   function is_legal_access_int_array_1d(ptr, i) result(isLegal)

      integer(kind=c_int), intent(in), pointer, contiguous, dimension(:) :: ptr
      integer :: i
      logical :: isLegal

      if (i > SIZE(ptr, 1)) then
         isLegal = .FALSE.
      else
         isLegal = .TRUE.
      endif
         
      return
   end function

   function is_legal_access_int_array_2d(ptr, i, j) result(isLegal)

      integer(kind=c_int), intent(in), pointer, contiguous, dimension(:,:) :: ptr
      integer :: i, j
      logical :: isLegal

      if (i > SIZE(ptr, 1) .OR. j > SIZE(ptr, 2) ) then
         isLegal = .FALSE.
      else
         isLegal = .TRUE.
      endif
         
      return
   end function

   function is_legal_access_int_array_3d(ptr, i, j, k) result(isLegal)

      integer(kind=c_int), intent(in), pointer, contiguous, dimension(:,:,:) :: ptr
      integer :: i, j, k
      logical :: isLegal

      if (i > SIZE(ptr, 1) .OR. j > SIZE(ptr, 2) .OR. k > SIZE(ptr, 3) ) then
         isLegal = .FALSE.
      else
         isLegal = .TRUE.
      endif
         
      return
   end function

   function is_legal_access_AngleSet_array_1d(ptr, i) result(isLegal)
      use AngleSet_mod

      type(AngleSet), intent(in), pointer, dimension(:) :: ptr
      integer :: i
      logical :: isLegal

      if (i > SIZE(ptr, 1) ) then
         isLegal = .FALSE.
      else
         isLegal = .TRUE.
      endif
         
      return
   end function

   function is_legal_access_SetData_array_1d(ptr, i) result(isLegal)
      use SetData_mod

      type(SetData), intent(in), pointer, dimension(:) :: ptr
      integer :: i
      logical :: isLegal

      if (i > SIZE(ptr, 1) ) then
         isLegal = .FALSE.
      else
         isLegal = .TRUE.
      endif
         
      return
   end function

end module ArrayChecks_mod
