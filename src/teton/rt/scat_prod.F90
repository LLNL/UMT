function scat_prod(x,y) result(sInnerProduct)

!=======================================================================
!                       Version 1.0: 03/99, MRZ
!-----------------------------------------------------------------------
! scat_prod
!   This function computes the scattering inner product, defined as:
!
!                  I
!     <x,y>   :=  SUM  x(i)*y(i)*sigmaS(i)*vol(i)
!          s      i=1
!
! The result is protected against underflow exceptions.
!
!   Routine is used for 2D and 3D; called by RTACCELMD.
!
! nconr        number of corners
!
! x(i)         arbitrary vector x
! y(i)         arbitrary vector y
! sigmaS(i)    within-group scattering cross section, volume i
! vol(i)       volume of spatial element i
!-----------------------------------------------------------------------
! v1.0: Original implementation
!=======================================================================

   use kind_mod
   use constant_mod
   use Size_mod
   use GreyAcceleration_mod
   use mpi_param_mod
   use mpif90_mod

!  variable declarations
   implicit none

!  passed variables
   real(adqt), intent(in) :: x(Size%ncornr)
   real(adqt), intent(in) :: y(Size%ncornr)

   real(adqt)             :: sInnerProduct

!  local variables
#ifdef UNDERFLOW
   real(adqt) :: tinySquareRoot, prodSquareRoot
   integer    :: i
#endif

!-----------------------------------------------------------------------
#ifdef UNDERFLOW

!  Compute the scattering inner product, protected against underflow
!  errors

!  initialize
   tinySquareRoot = sqrt(adqtTiny)

!  compute the scattering inner product, protecting against underflow

   sInnerProduct = zero
   do i = 1, Size%ncornr 
      prodSquareRoot = sign(one,x(i))*sqrt(abs(x(i))) * &
                       sign(one,y(i))*sqrt(abs(y(i))) * &
                       sqrt(GTA%GreySigScatVol(i))

      if (abs(prodSquareRoot) > tinySquareRoot) then
         sInnerProduct = sInnerProduct + prodSquareRoot*prodSquareRoot
      endif
   enddo

   call MPIAllReduce(sInnerProduct, "sum", MY_COMM_GROUP)

#else

!  Compute the scattering inner product, assuming that the compiler
!  sets any underflow to zero (the proper behavior in a sum such as
!  this)

   sInnerProduct = sum(x(:)*y(:)*GTA%GreySigScatVol(:))

   call MPIAllReduce(sInnerProduct, "sum", MY_COMM_GROUP)

#endif
!-----------------------------------------------------------------------

   return
end function scat_prod
