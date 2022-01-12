subroutine assoc_legendre(PLK,mu,num_degree,num_harm,num_angle)

!=======================================================================
!                       Version 1.0: 07/98, MRZ
!-----------------------------------------------------------------------
! This subroutine assigns the values of the associated Legendre 
! polynomials for degrees l=0,...num_degree, order k=0..l, evaluated at
! each of the x-direction cosines, mu(m).

! We do not want to store the associated Legendre polynomials as a
! triply-dimensioned array, so we use the following storage scheme:

!     k
!    P (mu[m]) == P[l,k](mu[m]) =: plk( m, l*(l+1)/2 + k + 1 )
!     l

! The first two associated Legendre polynomials are given: (Ref. 1)

!      P[0,0](mu[m]) == 1
!      P[1,0](mu[m]) == mu[m]

! The remaining associated Legendre polynomials are computed using
! recurrence relations.  For l>0, k>0 we use:  (Ref. 2)

!      P[l,k](mu[m]) == (1/sqrt(1-mu[m]^2))*{
!                                   (l-k+1)*mu[m]*P[l-1,k](mu[m]) -
!                                   (l+k-1)*P[l-1,k-1](mu[m]) }

! For l>1, k=0 we use:  (Ref. 2)

!      P[l,k](mu[m]) == (1/(l-l))*{ (2l-1)*mu[m]*P[l-1,k](mu[m]) -
!                                   (l+k-1)*P[l-2,k](mu[m]) }

!  References:

!  (1)  E.E. Lewis and W.F. Miller, Jr., Computational Methods of
!       Neutron Transport, American Nuclear Society, Inc., Illinois (1993).

!  (2)  I.A. Stegun, Legendre Functions, in Handbook of
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

real(adqt), intent(inout) :: plk(num_angle,num_harm)

!  Local variables
integer :: l, k, l_k, l_k1, l1_k1, l1_k, l2_k


!***********************************************************************

!  Loop over all degrees

      DegreeLoop: do l=0,num_degree

!     Loop over all orders

        OrderLoop: do k=0,l

!       the following indices map to the indicated polynomials,
!       each is computed only if needed 
!         l_k   --> P[l,k]
!         l_k1  --> P[l,k-1]
!         l1_k1 --> P[l-1,k-1]
!         l1_k  --> P[l-1,k]
!         l2_k  --> P[l-2,k]

!       compute the associated Legendre polynomial of degree l,
!       order k, evaluated at mu[m].

          l_k = (( l )*(( l )+1))/2 + ( k ) + 1 

          if (l==0 .and. k==0) then
!           P[0,0] = 1
            plk(:,l_k) = one

          elseif (l==1 .and. k==0) then
!           P[1,0] = mu
            plk(:,l_k) = mu(:)

          elseif (l>0 .and. k>0) then
!           P[l,k] defined using the first recurrence relation above

            l_k1  = (( l )*(( l )+1))/2 + (k-1) + 1
            l1_k1 = ((l-1)*((l-1)+1))/2 + (k-1) + 1

            plk(:,l_k) = (one/sqrt(one-mu(:)*mu(:))) *   &
                           ( real(l-k+1,adqt)*mu(:)*plk(:,l_k1) - &
                             real(l+k-1,adqt)*plk(:,l1_k1) )

          elseif (l>1 .and. k==0) then
!           P[l,k] defined using the second recurrence relation above

            l1_k  = ((l-1)*((l-1)+1))/2 + ( k ) + 1
            l2_k  = ((l-2)*((l-2)+1))/2 + ( k ) + 1

            plk(:,l_k) = (one/real(l-k,adqt)) *    &
                           ( real(2*l-1,adqt)*mu(:)*plk(:,l1_k) - &
                             real(l+k-1,adqt)*plk(:,l2_k) )
          endif

        enddo OrderLoop

      enddo DegreeLoop



return
end subroutine assoc_legendre 

 
