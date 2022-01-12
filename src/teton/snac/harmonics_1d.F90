subroutine harmonics_1d(YLK,mu,num_harm,num_degree,num_angle)

!  This subroutine assigns the values of the surface harmonics functions
!  for degrees l=0,...num_degree and orders k=-l...+l evaluated at each of
!  the direction cosines in the quadrature set for one-dimensional
!  geometries.  The surface harmonics are an orthonormal set when
!  integrated over the unit sphere (quadrature set normalized to 2).

!  We do not want to store the surface harmonics as a triply-dimensioned
!  array, so we use the following storage scheme:

!     0 
!    Y (OMEGA[m]) =: YLK( l+1, m )
!     l

!  The direction of particle travel is given:

!    OMEGA[m] == (mu[m])

!  where,

!    mu[m]  == cos(theta[m])

!  mu may be measured relative to different axes for
!  different geometries.

!  The surface harmonics are expressed as a function of mu[m]
!  (see Ref. 1)

!   0 
!  Y (OMEGA[m]) = Y[l](mu[m]) = a[l]*P[l](mu[m])
!   l

!  for degree l=0...num_degree, where

!    a[l]        = SQRT{ (2l+1)/2 }
!    P[l](mu[m]) = Legendre polynomial of degree l evaluated at mu[m] 

!  The first two Legendre polynomials are given: (Refs. 1 & 2)

!    P[0](mu[m]) = 1
!    P[1](mu[m]) = mu[m]

!  The remaining Legendre polynomials are computed using recurrence
!  relations:  (Refs. 1 & 2)

!    P[l](mu[m]) = (1/l)*{ (2l-1)*mu[m]*P[l-1](mu[m]) -
!                          (l-1)*P[l-2](mu[m]) }

!  References:

!  (1)  E.E. Lewis and W.F. Miller, Jr., Computational Methods of
!       Neutron Transport, American Nuclear Society, Inc., Illinois (1993).

!  (2)  U.W. Hochstrasser, Orthogonal Polynomials, in Handbook of
!       Mathematical Functions with Formulas, Graphs and Mathematical
!       Tables, M. Abramowitz & I.A. Stegun (eds.), Dover
!       Publications, Inc., New York (1972).

use kind_mod
use constant_mod

!  Variable declarations
implicit none

!  Arguments
integer,    intent(in)    :: num_harm,num_degree,num_angle

real(adqt), intent(in)    :: mu(num_angle)

real(adqt), intent(inout) :: ylk(num_harm,num_angle)

!  Local variables
integer                   :: l, ll, l1, l2 

!***********************************************************************

!  Loop over all degrees

      DegreeLoop: do l=0,num_degree

!     the following indices map to the indicated functions
!       ll --> P[l]
!       l1 --> P[l-1]
!       l2 --> P[l-2]

        ll = l + 1

!     construct the Legendre polynomial of degree 1 for all angles

        if (l == 0) then

          ylk(ll,:) = one

        elseif (l == 1) then

          ylk(ll,:) = mu(:)

        else

          l1 = l
          l2 = l - 1

          ylk(ll,:) = (one/real(l,adqt))*                     &
                         ( real(2*l-1,adqt)*mu(:)*ylk(l1,:) - &
                           real(l-1,adqt)*ylk(l2,:) )

        endif

      enddo DegreeLoop

!  Normalize to obtain the surface harmonics in each direction

      do l=0,num_degree
        ll = l + 1
        ylk(ll,:) = sqrt( real(2*l+1,adqt)/two )*ylk(ll,:)
      enddo


return
end subroutine harmonics_1d

 
