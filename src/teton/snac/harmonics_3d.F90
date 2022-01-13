subroutine harmonics_3d(YLK,mu,eta,xi,num_harm,num_degree,num_angle)

!  This subroutine assigns the values of the surface harmonics functions
!  for degrees l=0,...num_degree and orders k=-l...+l evaluated at each of
!  the direction cosines in the quadrature set for three-dimensional
!  geometries.  The surface harmonics are an orthonormal set when
!  integrated over the unit sphere (quadrature set normalized to 4*pi).

!  We do not want to store the surface harmonics as a triply-dimensioned
!  array, so we use the following storage scheme:

!     k
!    Y (OMEGA[m]) =: YLK( l*(l+1) + k + 1, m )
!     l

!  The direction of particle travel is given:

!    OMEGA[m] == (mu[m],eta[m],xi[m])

!  where,

!    mu[m]  == cos(theta[m])
!    eta[m] == sin(theta[m])*cos(w[m])
!    xi[m]  == sin(theta[m])*sin(w[m])

!  mu, eta, xi may be measured relative to different axes for
!  different geometries.

!  The surface harmonics are expressed as a function of mu[m] and w[m]
!  (see Refs. 1 & 2)

!   k
!  Y (OMEGA[m]) = Y[l,k](mu[m],w[m]) = a[l,k]*P[l,|k|](mu[m])*tau[k](w[m])
!   l

!  for degree l=0...num_degree, and order k=-l...+l, where

!    a[l,k] = SQRT{ [(2l+1)*(l-|k|)!] / [2*pi*(1+delta(k,0))*(l+|k|)!]}
!    P[l,k](mu[m]) = associated Legendre polynomial of degree l and
!                    order k evaluated at mu[m]
!    tau[k](w[m])  = cos(k*w[m])     for k >= 0
!                    sin(|k|*w[m])   for k < 0

!  References:

!  (1)  R.D. ODell and R.E. Alcouffe, Transport Calculations for
!       Nuclear Analyses:  Theory and Guidelines for Effective Use of
!       Transport Codes, LANL (1987).

!  (2)  P.N. Brown, A Linear Algebraic Development of Diffusion
!       Synthetic Acceleration for Three-Dimensional Transport
!       Equations, SIAM J. Numer. Anal., 32, 179 (1995).

use kind_mod
use constant_mod

!  Variable declarations
implicit none

!  External subprograms
integer, external :: factorial, kronecker

!  Arguments
integer,    intent(in)    :: num_harm,num_degree,num_angle

real(adqt), intent(in)    :: mu(num_angle),eta(num_angle),xi(num_angle)

real(adqt), intent(inout) :: ylk(num_harm,num_angle)

!  Local variables
real(adqt), allocatable :: plk(:,:)
real(adqt)              :: a,w,tau
integer                 :: num_plk_mom,l,k,l_k_ylk,l_k_plk,m,alloc_err


!***********************************************************************

!  Declare storage for the harmonics calculation

      num_plk_mom = ( (num_degree+1)*(num_degree+2) )/2

      allocate ( plk(num_angle,num_plk_mom) )

!  Assign the associate Legendre polynomials

      call assoc_legendre(plk,mu,num_degree,num_plk_mom,num_angle)

!  Loop over all degrees

      DegreeLoop: do l=0,num_degree

!     Loop over all orders

        OrderLoop: do k=-l,l

!       the following indices map to the indicated functions
!         l_k_plk --> P[l,|k|]
!         l_k_ylk --> Y[l,k]

          l_k_ylk = l*(l+1) + k + 1
          l_k_plk = l*(l+1)/2 + abs(k) + 1

          a = sqrt( (real((2*l+1)*factorial(l-abs(k)),adqt)) / &
                    (pi*real(2*(1+kronecker(k,0))*factorial(l+abs(k)),adqt)) )

!       Load the surface harmonic for every direction

          AngleLoop: do m=1,num_angle

!         compute the angle "w" and "tau"

            w = atan2( xi(m),eta(m) )

            if (k >= 0) then
              tau = cos(k*w)
            else
              tau = sin(abs(k)*w)
            endif

!         store the value of the surface harmonic

            ylk(l_k_ylk,m) = a*plk(m,l_k_plk)*tau

          enddo AngleLoop

        enddo OrderLoop

      enddo DegreeLoop

!  Free memory

      deallocate(plk, stat=alloc_err)


return
end subroutine harmonics_3d

 
