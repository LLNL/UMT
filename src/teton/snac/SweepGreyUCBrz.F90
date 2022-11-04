!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   SweepGreyUCBrz  - This routine calculates angular fluxes for a     *
!                     single direction and single energy groups for    *
!                     for an upstream corner-balance (UCB) spatial     *
!                     in rz-geometry. It is only used for Grey         *
!                     Transport Acceleration (GTA).                    *
!                                                                      *
!***********************************************************************

   subroutine SweepGreyUCBrz(setID, Angle, PsiB)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use SetData_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: setID
   integer,    intent(in)    :: Angle
   real(adqt), intent(inout) :: PsiB(Size%nbelem,Size%nangGTA)

!  Local

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet

   integer    :: c
   integer    :: c0
   integer    :: nCorner
   integer    :: ib
   integer    :: ii
   integer    :: angle0
   integer    :: zone
   integer    :: nzones
   integer    :: ndoneZ
   integer    :: hyperPlane
   integer    :: nHyperPlanes

   real(adqt) :: fac
   real(adqt) :: quadwt
   real(adqt) :: quadTauW1
   real(adqt) :: quadTauW2
   real(adqt) :: omega(2)

   logical(kind=1) :: StartingDirection

!  Constants

   Set         => getSetData(Quad,setID)
   ASet        => getAngleSetFromSetID(Quad, setID)

   nHyperPlanes      = getNumberOfHyperPlanes(ASet,Angle)
   angle0            = ASet% angle0
   omega(:)          = ASet% omega(:,Angle)
   quadwt            = ASet% weight(Angle)
   fac               = ASet% angDerivFac(Angle)
   quadTauW1         = ASet% quadTauW1(Angle)
   quadTauW2         = ASet% quadTauW2(Angle)
   StartingDirection = ASet% StartingDirection(Angle)

!  Initialize the temporary Psi to account for mesh cycles

   if (Size% useNewGTASolver ) then
     do c=1,Set% nCorner
       Set% tPsi(c) = zero
       Set% pInc(c) = zero
     enddo
   else
     do c=1,Set% nCorner
       Set% tPsi(c) = zero
     enddo
   endif

   do ib=1,Set%nbelem
     Set% tPsi(Set%nCorner+ib) = PsiB(ib,angle0+Angle)
   enddo

!  Loop through all of the corners using the NEXT list

   ndoneZ = 0 

   HyperPlaneLoop: do hyperPlane=1,nHyperPlanes

     nzones = getZonesInPlane(ASet,Angle,hyperPlane)

     if (Size% useNewGTASolver ) then

!$omp  parallel do default(none) schedule(static) &
!$omp& private(zone, nCorner, c0)  &
!$omp& shared(ASet, Geom, Psib, angle0, nzones, ndoneZ, Angle, setID, omega)  &
!$omp& shared(fac, quadTauW1, quadTauW2, StartingDirection)
       ZoneLoop1: do ii=1,nzones

         zone    = iabs( ASet% nextZ(ndoneZ+ii,Angle) )
         nCorner = Geom% numCorner(zone)
         c0      = Geom% cOffSet(zone)

         call SweepGreyUCBrzKernelNew(setID, PsiB, nCorner, c0, angle0,  &
                                      Angle, omega, fac, quadTauW1,      &
                                      quadTauW2, StartingDirection)

       enddo ZoneLoop1
!$omp end parallel do

     else

!$omp  parallel do default(none) schedule(static) &
!$omp& private(zone, nCorner, c0)  &
!$omp& shared(ASet, Geom, Psib, angle0, nzones, ndoneZ, Angle, setID, omega, fac, quadwt)
       ZoneLoop2: do ii=1,nzones

         zone    = iabs( ASet% nextZ(ndoneZ+ii,Angle) )
         nCorner = Geom% numCorner(zone)
         c0      = Geom% cOffSet(zone)

         call SweepGreyUCBrzKernel(setID, PsiB, nCorner, c0, angle0,  &
                                   Angle, omega, fac, quadwt)

     enddo ZoneLoop2
!$omp end parallel do

     endif

     ndoneZ = ndoneZ + nzones

   enddo HyperPlaneLoop

   if (Size% useNewGTASolver ) then
     GTA% PhiInc(:) = GTA% PhiInc(:) + quadwt*Set% pInc(:)
   endif

!  Set the "half-angle" angular intensity (PsiM) for the next angle

   if ( StartingDirection ) then
     do c=1,Set%nCorner
       Set% tPsiM(c) = Set% tPsi(c)
     enddo
   else
     do c=1,Set%nCorner
       Set% tPsiM(c) = quadTauW1*Set% tPsi(c) - quadTauW2*Set% tPsiM(c)
     enddo
   endif


   return
   end subroutine SweepGreyUCBrz



   subroutine SweepGreyUCBrzKernelNew(setID, PsiB, nCorner, c0, angle0,  &
                                      Angle, omega, fac, quadTauW1,      &
                                      quadTauW2, StartingDirection)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use SetData_mod
   use AngleSet_mod

   implicit none

!  Arguments

   real(adqt), intent(inout)   :: PsiB(Size%nbelem,Size%nangGTA)

   integer,    intent(in)      :: setID
   integer,    intent(in)      :: nCorner
   integer,    intent(in)      :: c0
   integer,    intent(in)      :: angle0
   integer,    intent(in)      :: Angle

   real(adqt), intent(in)      :: omega(2)
   real(adqt), intent(in)      :: fac
   real(adqt), intent(in)      :: quadTauW1
   real(adqt), intent(in)      :: quadTauW2

   logical(kind=1), intent(in) :: StartingDirection

!  Local

   type(SetData),    pointer   :: Set
   type(AngleSet),   pointer   :: ASet

   integer    :: b
   integer    :: i
   integer    :: ib
   integer    :: cface
   integer    :: cfp
   integer    :: cez
   integer    :: c
   integer    :: nxBdy 

   integer    :: nxez(Size% maxCorner)
   integer    :: ez_exit(2,Size% maxCorner)
   integer    :: bdy_exit(2,2*Size% maxCorner)

   real(adqt) :: sigA
   real(adqt) :: sigA2

   real(adqt) :: sez
   real(adqt) :: gnum
   real(adqt) :: gtau

   real(adqt) :: aez
   real(adqt) :: afp
   real(adqt) :: R
   real(adqt) :: R_afp

   real(adqt) :: denom(Size% maxCorner)
   real(adqt) :: coefpsi(2,Size% maxCorner)
   real(adqt) :: Sigt(Size% maxCorner)
   real(adqt) :: Q(Size% maxCorner)
   real(adqt) :: src(Size% maxCorner)

   real(adqt), parameter :: fouralpha=1.82d0

!  Constants

   Set  => getSetData(Quad,setID)
   ASet => getAngleSetFromSetID(Quad, setID)

   nxBdy   = 0
   nxez(:) = 0

!  Contributions from volume terms (if a starting direction add angular derivative)

   do c=1,nCorner
     Q(c)            = GTA%GreySigtInv(c0+c)*GTA%TsaSource(c0+c)
     src(c)          = Geom% Volume(c0+c)*GTA%TsaSource(c0+c) + &
                       fac*Geom% Area(c0+c)*Set% tPsiM(c0+c)
     Sigt(c)         = GTA%GreySigTotal(c0+c)
     denom(c)        = Sigt(c)*Geom% Volume(c0+c) + fac*Geom% Area(c0+c)
     Set% pInc(c0+c) = fac*Geom% Area(c0+c)*Set% tInc(c0+c)
   enddo

   CornerLoop: do c=1,nCorner

     CornerFaceLoop: do cface=1,2

       afp = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c0+c) )
       aez = DOT_PRODUCT( omega(:),Geom% A_ez(:,cface,c0+c) )

       cfp = Geom% cFP(cface,c0+c)

       if ( afp < zero ) then

         R_afp           = Geom% RadiusFP(cface,c0+c)*afp
         denom(c)        = denom(c)        - R_afp
         src(c)          = src(c)          - R_afp*Set% tPsi(cfp)
         Set% pInc(c0+c) = Set% pInc(c0+c) - R_afp*Set% tPsi(cfp)

       else

         if (cfp > Set% nCorner) then
           nxBdy             = nxBdy + 1
           bdy_exit(1,nxBdy) = c
           bdy_exit(2,nxBdy) = cfp - Set% nCorner
         endif

       endif

       if ( aez > zero ) then

         R    = Geom% RadiusEZ(cface,c0+c)
         cez  = Geom% cEZ(cface,c0+c)

         nxez(c)            = nxez(c) + 1
         ez_exit(nxez(c),c) = cez
         coefpsi(nxez(c),c) = R*aez
         denom(cez)         = denom(cez) + R*aez

         if ( afp < zero ) then

           sigA      = Sigt(c)*Geom% Area(c0+c)
           sigA2     = sigA*sigA

           gnum      = aez*aez*( fouralpha*sigA2 +  &
                       aez*(four*sigA + three*aez) )

           gtau      = gnum/                                           &
                     ( gnum + four*sigA2*sigA2 + aez*sigA*(six*sigA2 + &
                       two*aez*(two*sigA + aez)) )

           sez       = R*( gtau*sigA*( Set% tPsi(cfp) - Q(c) ) +  &
                       half*aez*(one - gtau)*( Q(c) - Q(cez) ) )

           src(c)    = src(c)    + sez
           src(cez)  = src(cez)  - sez

           Set% pInc(c0+c)   = Set% pInc(c0+c)   + R*gtau*sigA*Set% tPsi(cfp)
           Set% pInc(c0+cez) = Set% pInc(c0+cez) - R*gtau*sigA*Set% tPsi(cfp)

         else

           sez      = half*R*aez*( Q(c) - Q(cez) )
           src(c)   = src(c)    + sez
           src(cez) = src(cez)  - sez

         endif

       endif

     enddo CornerFaceLoop

   enddo CornerLoop


   do i=1,nCorner

     c = ASet% nextC(c0+i,angle)

!    Corner angular flux
     Set% tPsi(c0+c) = src(c)/denom(c)
     Set% pInc(c0+c) = Set% pInc(c0+c)/denom(c)

!    Calculate the contribution of this flux to the sources of
!    downstream corners in this zone. The downstream corner index is
!    "ez_exit."

     do cface=1,nxez(c)
       cez               = ez_exit(cface,c)
       src(cez)          = src(cez)  + coefpsi(cface,c)*Set% tPsi(c0+c)
       Set% pInc(c0+cez) = Set% pInc(c0+cez) + coefpsi(cface,c)*Set% pInc(c0+c)
     enddo

   enddo

!  Set exiting boundary fluxes

   do ib=1,nxBdy
     c                    = bdy_exit(1,ib)
     b                    = bdy_exit(2,ib)
     PsiB(b,angle0+Angle) = Set% tPsi(c0+c)
   enddo

   if ( StartingDirection ) then
     do c=1,nCorner
       Set% tInc(c0+c) = Set% pInc(c0+c)
     enddo
   else
     do c=1,nCorner
       Set% tInc(c0+c) = quadTauW1*Set% pInc(c0+c) - quadTauW2*Set% tInc(c0+c)
     enddo
   endif


   return
   end subroutine SweepGreyUCBrzKernelNEW 


   subroutine SweepGreyUCBrzKernel(setID, PsiB, nCorner, c0, angle0, &
                                   Angle, omega, fac, quadwt)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use SetData_mod
   use AngleSet_mod

   implicit none

!  Arguments

   real(adqt), intent(inout)   :: PsiB(Size%nbelem,Size%nangGTA)

   integer,    intent(in)      :: setID
   integer,    intent(in)      :: nCorner
   integer,    intent(in)      :: c0
   integer,    intent(in)      :: angle0
   integer,    intent(in)      :: Angle

   real(adqt), intent(in)      :: omega(2)
   real(adqt), intent(in)      :: fac
   real(adqt), intent(in)      :: quadwt

!  Local

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet

   integer    :: b
   integer    :: i
   integer    :: ib
   integer    :: cface
   integer    :: cez
   integer    :: cfp
   integer    :: ID
   integer    :: c
   integer    :: nxBdy

   integer    :: nxez(Size% maxCorner)
   integer    :: ez_exit(2,Size% maxCorner)
   integer    :: bdy_exit(2,2*Size% maxCorner)

   real(adqt) :: sigA
   real(adqt) :: sigA2
   real(adqt) :: sez
   real(adqt) :: gnum
   real(adqt) :: gtau
   real(adqt) :: aez
   real(adqt) :: afp
   real(adqt) :: R
   real(adqt) :: R_afp

   real(adqt) :: denom(Size% maxCorner)
   real(adqt) :: coefpsi(2,Size% maxCorner)
   real(adqt) :: Sigt(Size% maxCorner)
   real(adqt) :: Q(Size% maxCorner)
   real(adqt) :: src(Size% maxCorner)

   real(adqt), parameter :: fouralpha=1.82d0

!  Constants

   Set    => getSetData(Quad,setID)
   ASet   => getAngleSetFromSetID(Quad, setID) 
   ID     =  GTA% ID

   nxBdy   = 0
   nxez(:) = 0

!  Contributions from volume terms (if a starting direction add angular
!  derivative)

   do c=1,nCorner
     Q(c)     = GTA%GreySigtInv2(c0+c,ID)*GTA%TsaSource(c0+c)
     src(c)   = Geom% Volume(c0+c)*GTA%TsaSource(c0+c) + & 
                Geom% Area(c0+c)*fac*Set% tPsiM(c0+c)
     Sigt(c)  = GTA%GreySigTotal2(c0+c,ID)
     denom(c) = Sigt(c)*Geom% Volume(c0+c) + Geom% Area(c0+c)*fac
   enddo

   CornerLoop: do c=1,nCorner

     CornerFaceLoop: do cface=1,2

       afp = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c0+c) )
       aez = DOT_PRODUCT( omega(:),Geom% A_ez(:,cface,c0+c) )

       cfp = Geom% cFP(cface,c0+c)

       if ( afp < zero ) then

         R_afp    = Geom% RadiusFP(cface,c0+c)*afp
         denom(c) = denom(c) - R_afp
         src(c)   = src(c)   - R_afp*Set% tPsi(cfp)

       else

         if (cfp > Set% nCorner) then
           nxBdy             = nxBdy + 1
           bdy_exit(1,nxBdy) = c
           bdy_exit(2,nxBdy) = cfp - Set% nCorner
         endif

       endif

       if ( aez > zero ) then

         R    = Geom% RadiusEZ(cface,c0+c)
         cez  = Geom% cEZ(cface,c0+c)

         nxez(c)            = nxez(c) + 1
         ez_exit(nxez(c),c) = cez
         coefpsi(nxez(c),c) = R*aez
         denom(cez)         = denom(cez) + R*aez

         if ( afp < zero ) then

           sigA      = Sigt(c)*Geom% Area(c0+c)
           sigA2     = sigA*sigA

           gnum      = aez*aez*( fouralpha*sigA2 +  &
                       aez*(four*sigA + three*aez) )

           gtau      = gnum/                                           &
                     ( gnum + four*sigA2*sigA2 + aez*sigA*(six*sigA2 + &
                       two*aez*(two*sigA + aez)) )

           sez       = R*( gtau*sigA*( Set% tPsi(cfp) - Q(c) ) +  &
                       half*aez*(one - gtau)*( Q(c) - Q(cez) ) )

           src(c)    = src(c)    + sez
           src(cez)  = src(cez)  - sez

         else

           sez      = half*R*aez*( Q(c) - Q(cez) )
           src(c)   = src(c)    + sez
           src(cez) = src(cez)  - sez

         endif

       endif

     enddo CornerFaceLoop

   enddo CornerLoop


   do i=1,nCorner

     c = ASet% nextC(c0+i,angle)

!    Corner angular flux
     Set% tPsi(c0+c) = src(c)/denom(c)
     Set% tPhi(c0+c) = Set% tPhi(c0+c) + quadwt*Set% tPsi(c0+c)

!    Calculate the contribution of this flux to the sources of
!    downstream corners in this zone. The downstream corner index is
!    "ez_exit."

     do cface=1,nxez(c)
       cez      = ez_exit(cface,c)
       src(cez) = src(cez) + coefpsi(cface,c)*Set% tPsi(c0+c)
     enddo

   enddo

!  Set exiting boundary fluxes

   do ib=1,nxBdy
     c                    = bdy_exit(1,ib)
     b                    = bdy_exit(2,ib)
     PsiB(b,angle0+Angle) = Set% tPsi(c0+c)
   enddo


   return
   end subroutine SweepGreyUCBrzKernel


