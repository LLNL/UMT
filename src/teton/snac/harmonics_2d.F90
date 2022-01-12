subroutine harmonics_2d(YLK,mu,eta,num_harm,num_degree,num_angle)

!  This subroutine assigns the values of the surface harmonics functions
!  for degrees l=0,...num_degree and orders k=-l...+l evaluated at each of
!  the direction cosines in the quadrature set for two-dimensional
!  geometries.  The surface harmonics are an orthonormal set when
!  integrated over the unit sphere (quadrature set normalized to 2*pi).

!  We do not want to store the surface harmonics as a triply-dimensioned
!  array, so we use the following storage scheme:

!     k
!    Y (OMEGA[m]) =: YLK( l*(l+1)/2 + k + 1, m )
!     l

!  The direction of particle travel is given:

!    OMEGA[m] == (mu[m],eta[m])

!  where,

!    mu[m]  == cos(theta[m])
!    eta[m] == sin(theta[m])*cos(w[m])

!  mu, eta may be measured relative to different axes for
!  different geometries.

!  The surface harmonics are expressed as a function of mu[m] and w[m]
!  (see Refs. 1 & 2)

!   k
!  Y (OMEGA[m]) = Y[l,k](mu[m],w[m]) = a[l,k]*P[l,k](mu[m])*cos(k*w[m])
!   l

!  for degree l=0...num_degree, and order k=0...l, where

!    a[l,k] = SQRT{ [(2l+1)*(l-k)!] / [pi*(1+delta(k,0))*(l+k)!]}
!    P[l,k](mu[m]) = associated Legendre polynomial of degree l and
!                    order k evaluated at mu[m]

!  References:

!  (1)  E.E. Lewis and W.F. Miller, Jr., Computational Methods of
!       Neutron Transport, American Nuclear Society, Inc., Illinois (1993).

!  (2)  R.D. ODell and R.E. Alcouffe, Transport Calculations for
!       Nuclear Analyses:  Theory and Guidelines for Effective Use of
!       Transport Codes, LANL (1987).

use kind_mod
use constant_mod

!  Variable declarations
implicit none

!  External subprograms
integer, external :: factorial, kronecker

!  Arguments
integer,    intent(in)    :: num_harm,num_degree,num_angle

real(adqt), intent(in)    :: mu(num_angle),eta(num_angle)

real(adqt), intent(inout) :: ylk(num_harm,num_angle)

!  Local variables
real(adqt), allocatable :: plk(:,:)
real(adqt)              :: a,w
integer                 :: l,k,l_k,m,alloc_err


!***********************************************************************

!  Declare storage for the harmonics calculation

      allocate ( plk(num_angle,num_harm) )

!  Assign the associate Legendre polynomials

      call assoc_legendre(plk,mu,num_degree,num_harm,num_angle)

!  Loop over all degrees

      DegreeLoop: do l=0,num_degree

!     Loop over all orders

        OrderLoop: do k=0,l

!       the following indices map to the indicated functions
!         l_k --> P[l,k] and Y[l,k]

          l_k = l*(l+1)/2 + k + 1

          a = sqrt( (real((2*l+1)*factorial(l-k),adqt)) / &
                    (pi*real((1+kronecker(k,0))*factorial(l+k),adqt)) )

!       Load the surface harmonic for every direction

          AngleLoop: do m=1,num_angle

!         compute the angle "w" 

            w = acos( eta(m)/sqrt(one - mu(m)*mu(m)) )

!         store the value of the surface harmonic

            ylk(l_k,m) = a*plk(m,l_k)*cos(k*w)

          enddo AngleLoop

        enddo OrderLoop

      enddo DegreeLoop

!  Free memory

      deallocate(plk, stat=alloc_err)


return
end subroutine harmonics_2d

 
