!***********************************************************************
!                                                                      *
!                       Last Update: 09/2017 by PFN                    *
!                                                                      *
!   reflectArbitrary - This routine, given an incident angle and       *
!                      outward normal of a boundary surface, finds     *
!                      the reflected angles in the set of angles       *
!                      provided.                                       *
!                                                                      *
!                                                                      *
!     The reflection algorithm works with angular currents instead     *
!     of angular fluxes.  (We believe currents are the important       *
!     quantities at boundaries.)  For a given incoming angle, the      *
!     algorithm picks one outgoing angle for its partner, and          *
!     requires that the net current caused by the pair vanish:         *
!                                                                      *
!      w(mm)*|cosine(mm)|*psi(mm) = w(mref)*|cosine(mref)|*psi(mref)   *
!                                                                      *
!     This routine returns mref and the ratio:                         *
!                                                                      *
!        cosrat  =  (w(mref)*|cosine(mref)|) / (w(mm)*|cosine(mm)|)    *
!                                                                      *
!     The main 'flaw' in our algorithm is that an isotropic angular    *
!     flux will not be isotropically reflected.  (The net current will *
!     always be zero, however.)                                        *
!                                                                      *
!     Input:                                                           *
!       Minc    index of the angle whose partner we must find,         *
!       Area    area vector,                                           *
!       omega   array of quadrature direction cosines,                 *
!                                                                      *
!     Output:                                                          *
!       mref    index of the partner angle,                            *
!       cosrat  ratio of cosines;  see above                           *
!                                                                      *
!***********************************************************************

   subroutine reflectArbitrary(ndim,AngleInc, NumAngles, nReflect, &
                               ReflAngle, omega, weight, Area, cosRatio)

   use kind_mod
   use constant_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: ndim
   integer,    intent(in)    :: AngleInc
   integer,    intent(in)    :: NumAngles
   integer,    intent(inout) :: nReflect
   integer,    intent(inout) :: ReflAngle(NumAngles) 

   real(adqt), intent(in)    :: omega(ndim,NumAngles) 
   real(adqt), intent(in)    :: weight(NumAngles)
   real(adqt), intent(in)    :: Area(ndim)
   real(adqt), intent(inout) :: cosratio(NumAngles) 

!  Local

   integer    :: ia,d,m1,m2

   real(adqt) :: dot,dotWW,dot_Minc,  &
                 AreaMag,AreaMagInv,wf_sum,eps,wf,wf_m

   real(adqt) :: omega_Inc(ndim)
   real(adqt) :: omega_refl(ndim,NumAngles)

!  Constants

   real(adqt), parameter :: tol=1.0e-10_adqt


!  Most of the time, we will have reflection in one of the directions
!  of the coordinate system.  The component of OMEGA that changes sign
!  corresponds to the component of the normal to the side (or face)
!  that is non-zero ( e.g. if Area(1).ne.0 then we have reflection in 
!  mu=omega(1)  ).  First we check for the component that changes
!  sign ( omega(d) ), and then we check to see that the other
!  components are unchanged.

   omega_Inc(:) = omega(:,AngleInc)

   AreaMag = zero

   do d=1,ndim
     AreaMag = AreaMag + Area(d)*Area(d)
   enddo

!------------------------------------------------------------------------------
! THEORY: REFLECTIVE B.C. on ARBITRARILY ORIENTED SURFACES
!
! If the surface is not orthogonal then we have to play some tricks to get
! reflected intensities.
!
! The exact reflected direction can be computed as:
!
!                omega_xit .dot. A_fep
!      A-prime = --------------------- A_fep
!                   A_fep .dot A_fep
!
!      Omega_refl = omega_xit - 2*A-prime
!
! Since we know the intensity along this "exact direction" (which is probably
! not in our set of discrete ordinate directions) we must "spread" this
! intensity (psi_w) among nearby d.o. directions.  This must be done in a
! fashion that preserves the net current flow across the boundary.  For a
! reflective surface the net flow should be zero.  We also want all of our
! reflected angular fluxes to be positive.
!
! Given an exiting direction we would compute Omega_refl and then "spread" this
! angular flux among "nearby" discrete ordinate directions.  We choose this
! set to be:
!
!     M = {k :: omega(k) .dot. A_fep < 0}
!               \----------------------/
!                 incoming directions
!
! To preserve the net current flow across A_fep we require (let j be the
! exiting direction of interest):
!
!    ( A .dot. omega_j ) * psi_j * wt_j
!             = - SUM ( A .dot. omega(k) ) * psi(k) * wt(k)
!               [k in {M}]
!
! This leads us to the following form for the j-th component of Psi(k),
! k is in {M}.
!
!           / wt_j * A .dot. omega_j \           delta(k,j)
! psi_k = - | ---------------------- | * psi_j * ----------
!           \ wt_k * A .dot. omega_k /              D(j)
!
! For now we let the weight delta(k,j) be defined as (1/sin(theta_k)).
! Where theta_k is the angle between omega_k and the reflection of omega_j
!
!     Theta_k = 1/2 * acos (omega_k .dot. omega_refl_j)
!
! We also requrie that D(j) = SUM [ delta(k,j) ]
!                           [k in {M}]
!
! For balance to hold     SUM [ delta(k,j) / D(j) ] == 1
!                       [k in {M}]
!
! Note that our angle of interest Omega_m is a member of {M}
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! IMPLEMENTATION:
!
! The theory above assumes that we have an exiting direction and intensity,
! psi(iw,m).  However, in SNAC, we need to build an incoming intensity in a
! specified direction.  To do this we will establish a set of exiting
! directions Omega_j in the set {M-tilde} that can contribute in building
! the reflected direction psi_m (k=m is a member of set {M}).  We define this
! set of exiting directions to be:
!
! {M-tilde}
!      = { j:: omega(j) .dot. A_fep > 0 }
!              \----------------------/
!                 exiting directions
!
! We define omega_refl(j) to be the exact reflected direction (not necessarily
! a discrete ordinate direction) of omega(j).  Omega(k=m) is the direction for
! which we are trying to build psi(m).
!
! It should be noted that the set {M}_j from the Theory section depends on
! omega_refl(j).
!
! To build psi(k=m) we will loop over the directions omega(j) in {M-tilde}.  For
! each of these directions we determine the contribution from psi(omega(j)) to
! psi(omega(k=m)).  Finally, we keep a running sum to build psi(k=m) from each of
! these contributing directions.
!------------------------------------------------------------------------------

! Determine the components for the set {M-tilde}.
! These are the exiting directions (omega_j) that can contribute to psi_m.
! Find the exact reflected direction omega_refl
!------------------------------------------------------------------------------

   dot_Minc   = DOT_PRODUCT( omega_inc(:),Area(:) )
   AreaMagInv = one/AreaMag
   nreflect   = 0

   do ia=1,NumAngles
     dot = DOT_PRODUCT( omega(:,ia),Area(:) )

     if (dot > zero) then
       nreflect               = nreflect + 1
       ReflAngle(nreflect)    = ia
       omega_refl(:,nreflect) = omega(:,ia) - two*dot*AreaMagInv*Area(:)
     endif
   enddo

! Loop over each of the directions j in the set {M_tilde}.  For each j:
!    1. Build the set {M}_j, which is the set of directions that omega_j can
!        contribute to.  Omega(k=m) is a member of {M}.
!    2. Build the "sine sum" as we build {M}.
!    3. With the "sine sum" we can now define the contribution coefficient for
!       psi(M_tilde_set(j)) toward psi_m.  This is the variable "cosratio".
!------------------------------------------------------------------------------

   ReflectAngles: do ia=1,nreflect

     wf_sum = zero
     m2     = ReflAngle(ia)

! If the dot product of the direction of interest (Minc) and an
! exact reflected direction is close to one (i.e. the angle between
! the two vectors is close to zero) the weight will be very large. 

     dotWW = min( DOT_PRODUCT( omega_inc(:),omega_refl(:,ia) ), one )
     eps   = one - dotWW + tol
     wf_m  = one/(eps*eps*eps*eps)

! Build {M} and compute the sum of the weights (wf_sum).

     do m1=1,NumAngles
       dot = DOT_PRODUCT( omega(:,m1),Area(:) )
       if (dot < zero) then

! Is the angle in the same hemisphere as omega_refl(j)?
! Omega_m1 .dot. Omega_refl(j) > 0?

         dotWW = min( DOT_PRODUCT( omega(:,m1),omega_refl(:,ia) ), one )
         eps   = one - dotWW + tol
         wf    = one/(eps*eps*eps*eps)
         wf_sum = wf_sum + wf

       endif
     enddo

! Find the fraction of psi_j that contributes to the building of psi(k=m).

     dot          = DOT_PRODUCT( omega(:,m2),Area(:) )
     cosratio(ia) = -(weight(m2)  *dot*wf_m)/  &
                     (weight(AngleInc)*dot_Minc*wf_sum)

   enddo ReflectAngles



   return
   end subroutine reflectArbitrary 

