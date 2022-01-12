!***********************************************************************
   subroutine snpnmset(self, ndim, isctp1)

   use kind_mod
   use constant_mod
   use Quadrature_mod

   implicit none

!  Arguments

   type(Quadrature)          :: self

   integer,    intent(in)    :: ndim, isctp1

!  Local

   integer    :: m, d, momentCount 

   real(adqt) :: pifactor

!***********************************************************************
!
!     version 1 1/94 by marvin adams
!
!
!     This routine, called by sngsrcxy, sets the pnm array.  Pnm(n,m)
!     is the nth spherical harmonic evaluated at the mth quadrature
!     direction, multiplied by a few constants.
!
!     If iscat = isctp1-1 is the order of anisotropic scattering, then
!     we need the (l,k)th spherical harmonic for l = 0,...,iscat and
!     k = 0,...,l.
!
!     We do not want to mess with doubly-dimensioned spherical-harmonic
!     functions, so we map from (l,k) --> n as follows:
!
!                         n = 1 + l*(l+1)/2 + k  .
!
!     Given that mapping, our pnm(n,m) is given by:
!
!                             2l+1 (l-k)! 2-delta(k,0)
!                Pnm(n,m)  =  ---- ------ ------------ Ylk(m) ,
!                              pi  (l+k)!      2
!
!     where
!
!          Ylk(m)     = Real part of (l,k)th spherical harmonic,
!                       evaluated at mth direction,
!
!          delta(k,0) = 0 if k.ne.0;  = 1 if k=0,
!
!          (j)!       = j factorial.
!
!***********************************************************************
 
!  Check to make sure that the moment count and scattering order agree

   momentCount = 0
   pifactor    = zero

   if (ndim == 1) then
     momentCount = isctp1
     pifactor    = half
   elseif (ndim == 2) then
     momentCount = isctp1*(isctp1+1)/2
     pifactor    = one/(two*pi)
   elseif (ndim == 3) then
     momentCount = isctp1*isctp1
     pifactor    = one/(four*pi)
   else
     call f90fatal("Incorrect spatial dimension in SNPNMSET")
   endif
 
   if (self% NumMoments /= momentCount) then
     call f90fatal("Moment count wrong in SNPNMSET")
   endif
 
!  Loop over quadrature angles:

   do m=1,self% NumAngles
 
!    l=k=0   =>  Ylk = 1.0, n=1
 
     self% pnm(1,m) = pifactor 
 
!    l = 1, k=0  =>  Y=u, n=2;  l=k=1 => Y=eta, n=3.
 
     if (isctp1 > 1) then

       do d=1,ndim
         self% pnm(d+1,m) = three*pifactor*self% omega(d,m)
       enddo

     endif

     if (isctp1 > 2) then
       call f90fatal("Not ready for scat ord > 1 in SNPNMSET")
     endif

   enddo

 
   return
   end subroutine snpnmset


