function kronecker(m,n)

!***********************************************************************
!                        Version 1:  07/98, MRZ                        *
!                                                                      *
!  This function computes the Dirac delta function:                    *
!                                                                      *
!                     { 1  m==n                                        *
!   kronecker(m,n) =  {                                                *
!                     { 0  m /= n                                      *
!                                                                      *
!***********************************************************************

!  Variable declarations
implicit none

!  Arguments
integer  ::  kronecker, m, n

!  Compute the Kronecker delta function

      if (m == n) then
        kronecker = 1
      else
        kronecker = 0
      endif


return
end function kronecker

