function factorial(n)

!  This function computes the factorial of an integer:

!    factorial(n) == n*(n-1)*(n-2)* ... *2*1

!  Variable declarations
implicit none

!  Arguments
integer :: factorial, n

!  Local variables
integer :: i

!  Compute the factorial

      factorial = 1
      do i=1,n
        factorial = factorial*i
      enddo

return

end function factorial

