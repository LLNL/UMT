!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   SweepGreyUCBxyz  - This routine calculates angular fluxes for a    *
!                      single direction and single energy groups for   *
!                      for an upstream corner-balance (UCB) spatial    *
!                      in xyz-geometry. It is only used for Grey       *
!                      Transport Acceleration (GTA).                   *
!                                                                      *
!***********************************************************************

   subroutine SweepGreyUCBxyz(setID, Angle, PsiB)

   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use Geometry_mod
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

   real(adqt) :: quadwt
   real(adqt) :: omega(3)

!  Constants

   Set          => getSetData(Quad,setID)
   ASet         => getAngleSetFromSetID(Quad, setID)
   angle0       =  Set% angle0
   nHyperPlanes =  getNumberOfHyperPlanes(ASet,Angle)

   omega(:)     =  ASet% omega(:,Angle)
   quadwt       =  ASet% weight(Angle)

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
!$omp& shared(ASet, Geom, Psib, angle0, nzones, ndoneZ, Angle, setID, omega)
     ZoneLoop1: do ii=1,nzones

       zone    = iabs( ASet% nextZ(ndoneZ+ii,Angle) ) 
       nCorner = Geom% numCorner(zone)
       c0      = Geom% cOffSet(zone)

       call SweepGreyUCBxyzKernelNew(setID, PsiB, nCorner, c0, angle0, &
                                     Angle, omega)

     enddo ZoneLoop1
!$omp end parallel do

     else

!$omp  parallel do default(none) schedule(static) &
!$omp& private(zone, nCorner, c0)  &
!$omp& shared(ASet, Geom, Psib, angle0, nzones, ndoneZ, Angle, setID, omega, quadwt)
     ZoneLoop2: do ii=1,nzones

       zone    = iabs( ASet% nextZ(ndoneZ+ii,Angle) )
       nCorner = Geom% numCorner(zone)
       c0      = Geom% cOffSet(zone)

       call SweepGreyUCBxyzKernel(setID, PsiB, nCorner, c0, angle0, &
                                  Angle, omega, quadwt)

     enddo ZoneLoop2
!$omp end parallel do

     endif

     ndoneZ = ndoneZ + nzones

   enddo HyperPlaneLoop

   if (Size% useNewGTASolver ) then
     GTA% PhiInc(:) = GTA% PhiInc(:) + quadwt*Set% pInc(:)
   endif



   return
   end subroutine SweepGreyUCBxyz



   subroutine SweepGreyUCBxyzKernelNew(setID, PsiB, nCorner, c0, angle0, &
                                       Angle, omega) 

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

   real(adqt), intent(inout) :: PsiB(Size%nbelem,Size%nangGTA)

   integer,    intent(in)    :: setID
   integer,    intent(in)    :: nCorner
   integer,    intent(in)    :: c0
   integer,    intent(in)    :: angle0
   integer,    intent(in)    :: Angle

   real(adqt), intent(in)    :: omega(3)

!  Local Variables

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet

   integer    :: b
   integer    :: i
   integer    :: ib
   integer    :: cface
   integer    :: ifp
   integer    :: k
   integer    :: nxBdy
   integer    :: c 
   integer    :: cez 
   integer    :: cfp
   integer    :: nCFaces

   integer    :: nxez(Size% maxCorner)
   integer    :: ez_exit(Size%maxcf,Size% maxCorner)
   integer    :: bdy_exit(2,Size%maxcf*Size% maxCorner)

   real(adqt) :: aez
   real(adqt) :: area_opp
   real(adqt) :: sigv
   real(adqt) :: sigv2
   real(adqt) :: gnum
   real(adqt) :: gtau
   real(adqt) :: sez
   real(adqt) :: psi_opp

   real(adqt) :: denom(Size% maxCorner)
   real(adqt) :: afp(Size%maxcf)
   real(adqt) :: coefpsi(Size%maxcf,Size% maxCorner)
   real(adqt) :: psifp(Size%maxcf)
   real(adqt) :: src(Size% maxCorner) 
   real(adqt) :: Q(Size% maxCorner)

   real(adqt), parameter :: fouralpha=1.82_adqt

!  Constants

   Set     => getSetData(Quad,setID)
   ASet    => getAngleSetFromSetID(Quad, setID)

   nxBdy   =  0

!  Contributions from volume terms

   do c=1,nCorner
     Q(c)   = GTA%GreySigtInv(c0+c)*GTA%TsaSource(c0+c) 
     src(c) = Geom% Volume(c0+c)*GTA%TsaSource(c0+c)
   enddo

   nxez(:) = 0

   CornerLoop: do c=1,nCorner

     sigv     = Geom% Volume(c0+c)*GTA%GreySigTotal(c0+c)
     denom(c) = sigv 
     nCFaces  = Geom% nCFacesArray(c0+c)

!    Contributions from external corner faces (FP faces)

     do cface=1,nCFaces

       afp(cface) = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c0+c) )

       cfp = Geom% cFP(cface,c0+c)

       if ( afp(cface) > zero ) then
         denom(c) = denom(c) + afp(cface)

         if (cfp > Set% nCorner) then
           nxBdy             = nxBdy + 1
           bdy_exit(1,nxBdy) = c
           bdy_exit(2,nxBdy) = cfp - Set% nCorner
         endif

       elseif ( afp(cface) < zero ) then

         psifp(cface) = Set% tPsi(cfp)

         src(c)          = src(c)          - afp(cface)*psifp(cface)
         Set% pInc(c0+c) = Set% pInc(c0+c) - afp(cface)*psifp(cface)
       endif
     enddo

!    Contributions from interior corner faces (EZ faces)

     do cface=1,ncfaces

       aez = DOT_PRODUCT( omega(:),Geom% A_ez(:,cface,c0+c) )
       cez = Geom% cEZ(cface,c0+c)

       if (cez > c) then
         if (aez > zero ) then
           nxez(c)                = nxez(c)   + 1
           ez_exit(nxez(c),c)     = cez
           coefpsi(nxez(c),c)     = aez
         elseif (aez < zero) then
           nxez(cez)              = nxez(cez) + 1
           ez_exit(nxez(cez),cez) = c
           coefpsi(nxez(cez),cez) = -aez
         endif
       endif

       if (aez > zero ) then

         psi_opp  = zero
         area_opp = zero
         denom(c) = denom(c) + aez

         ifp = mod(cface,nCFaces) + 1
                                                                                                   
         if ( afp(ifp) < zero ) then
           area_opp = -afp(ifp)
           psi_opp  = -afp(ifp)*psifp(ifp)
         endif
                                                                                                   
         do k=2,nCFaces-2
           ifp = mod(ifp,nCFaces) + 1
           if ( afp(ifp) < zero ) then
             area_opp = area_opp - afp(ifp)
             psi_opp  = psi_opp  - afp(ifp)*psifp(ifp)
           endif
         enddo
        
         TestOppositeFace: if (area_opp > zero) then

           psi_opp   = psi_opp/area_opp
           sigv2     = sigv*sigv

           gnum      = aez*aez*( fouralpha*sigv2 +    &
                       aez*(four*sigv + three*aez) )

           gtau      = gnum/    &
                     ( gnum + four*sigv2*sigv2 + aez*sigv*(six*sigv2 + &
                       two*aez*(two*sigv + aez)) )

           sez       = gtau*sigv*( psi_opp - Q(c) ) +  &
                       half*aez*(one - gtau)*( Q(c) - Q(cez) )
           src(c)    = src(c)    + sez
           src(cez)  = src(cez)  - sez

           Set% pInc(c0+c)   = Set% pInc(c0+c)   + gtau*sigv*psi_opp
           Set% pInc(c0+cez) = Set% pInc(c0+cez) - gtau*sigv*psi_opp
          
         else

           sez       = half*aez*( Q(c) - Q(cez) )
           src(c)    = src(c)    + sez
           src(cez)  = src(cez)  - sez
          
         endif TestOppositeFace

       endif

     enddo

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
       src(cez)          = src(cez)          + coefpsi(cface,c)*Set% tPsi(c0+c)
       Set% pInc(c0+cez) = Set% pInc(c0+cez) + coefpsi(cface,c)*Set% pInc(c0+c)
     enddo

   enddo
       
!  Set exiting boundary fluxes

   do ib=1,nxBdy
     c                    = bdy_exit(1,ib)
     b                    = bdy_exit(2,ib)
     Psib(b,angle0+Angle) = Set% tPsi(c0+c)
   enddo


   return
   end subroutine SweepGreyUCBxyzKernelNew 

   subroutine SweepGreyUCBxyzKernel(setID, PsiB, nCorner, c0, angle0,  & 
                                    Angle, omega, quadwt)

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

   real(adqt), intent(inout) :: PsiB(Size%nbelem,Size%nangGTA)

   integer,    intent(in)    :: setID
   integer,    intent(in)    :: nCorner
   integer,    intent(in)    :: c0
   integer,    intent(in)    :: angle0
   integer,    intent(in)    :: Angle

   real(adqt), intent(in)    :: omega(3)
   real(adqt), intent(in)    :: quadwt

!  Local Variables

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet

   integer    :: b
   integer    :: i
   integer    :: ib
   integer    :: cface
   integer    :: ifp
   integer    :: k
   integer    :: nxBdy

   integer    :: c
   integer    :: cez
   integer    :: cfp
   integer    :: ID
   integer    :: nCFaces

   integer    :: nxez(Size% maxCorner)
   integer    :: ez_exit(Size%maxcf,Size% maxCorner)
   integer    :: bdy_exit(2,Size%maxcf*Size% maxCorner)

   real(adqt) :: fouralpha,aez,area_opp,sigv,sigv2,  &
                 gnum,gtau,sez,psi_opp

   real(adqt) :: coef(Size% maxCorner)
   real(adqt) :: afp(Size%maxcf)
   real(adqt) :: coefpsi(Size%maxcf,Size% maxCorner)
   real(adqt) :: psifp(Size%maxcf)
   real(adqt) :: src(Size% maxCorner)
   real(adqt) :: Q(Size% maxCorner)
   real(adqt) :: SigtVol(Size% maxCorner)

!  Constants

   parameter (fouralpha=1.82d0)


   Set     => getSetData(Quad,setID)
   ASet    => getAngleSetFromSetID(Quad, setID) 
   ID      =  GTA% ID

   psi_opp = zero
   nxBdy   = 0
   nxez(:) = 0

!  Contributions from volume terms

   do c=1,nCorner
     Q(c)       = GTA%GreySigtInv2(c0+c,ID)*GTA%TsaSource(c0+c)
     src(c)     = Geom% Volume(c0+c)*GTA%TsaSource(c0+c)
     SigtVol(c) = Geom% Volume(c0+c)*GTA%GreySigTotal2(c0+c,ID)
   enddo

   CornerLoop: do c=1,nCorner

     sigv    = SigtVol(c)
     coef(c) = SigtVol(c)
     nCFaces = Geom% nCFacesArray(c0+c)

!    Contributions from external corner faces (FP faces)

     do cface=1,nCFaces

       afp(cface) = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c0+c) )

       cfp = Geom% cFP(cface,c0+c)

       if ( afp(cface) > zero ) then
         coef(c) = coef(c) + afp(cface)

         if (cfp > Set% nCorner) then
           nxBdy             = nxBdy + 1
           bdy_exit(1,nxBdy) = c
           bdy_exit(2,nxBdy) = cfp - Set% nCorner
         endif

       elseif ( afp(cface) < zero ) then

         psifp(cface) = Set% tPsi(cfp)

         src(c) = src(c) - afp(cface)*psifp(cface)
       endif
     enddo

!    Contributions from interior corner faces (EZ faces)

     do cface=1,ncfaces

       aez = DOT_PRODUCT( omega(:),Geom% A_ez(:,cface,c0+c) )
       cez = Geom% cEZ(cface,c0+c)

       if (cez > c) then
         if (aez > zero ) then
           nxez(c)                = nxez(c)   + 1
           ez_exit(nxez(c),c)     = cez
           coefpsi(nxez(c),c)     = aez
         elseif (aez < zero) then
           nxez(cez)              = nxez(cez) + 1
           ez_exit(nxez(cez),cez) = c
           coefpsi(nxez(cez),cez) = -aez
         endif
       endif

       if (aez > zero ) then

         area_opp = zero
         coef(c)  = coef(c) + aez

         if (nCFaces == 3) then

           ifp = mod(cface,nCFaces) + 1

           if ( afp(ifp) < zero ) then
             area_opp = -afp(ifp)
             psi_opp  =  psifp(ifp)
           endif

         else

           ifp      = cface
           area_opp = zero
           psi_opp  = zero

           do k=1,nCFaces-2
             ifp = mod(ifp,nCFaces) + 1
             if ( afp(ifp) < zero ) then
               area_opp = area_opp - afp(ifp)
               psi_opp  = psi_opp  - afp(ifp)*psifp(ifp)
             endif
           enddo

           psi_opp = psi_opp/area_opp

         endif

         TestOppositeFace: if (area_opp > zero) then

           sigv2     = sigv*sigv

           gnum      = aez*aez*( fouralpha*sigv2 +    &
                       aez*(four*sigv + three*aez) )

           gtau      = gnum/    &
                     ( gnum + four*sigv2*sigv2 + aez*sigv*(six*sigv2 + &
                       two*aez*(two*sigv + aez)) )

           sez       = gtau*sigv*( psi_opp - Q(c) ) +  &
                       half*aez*(one - gtau)*( Q(c) - Q(cez) )
           src(c)    = src(c)    + sez
           src(cez)  = src(cez)  - sez

         else

           sez       = half*aez*( Q(c) - Q(cez) )
           src(c)    = src(c)    + sez
           src(cez)  = src(cez)  - sez

         endif TestOppositeFace

       endif

     enddo

   enddo CornerLoop


   do i=1,nCorner

     c = ASet% nextC(c0+i,angle) 

!    Corner angular flux
     Set% tPsi(c0+c) = src(c)/coef(c)
     Set% tPhi(c0+c) = Set% tPhi(c0+c) + quadwt*Set% tPsi(c0+c)

!    Calculate the contribution of this flux to the sources of
!    downstream corners in this zone. The downstream corner index is
!    "ez_exit."

     do cface=1,nxez(c)
       cez       = ez_exit(cface,c)
       src(cez)  = src(cez) + coefpsi(cface,c)*Set% tPsi(c0+c)
     enddo

   enddo

!  Set exiting boundary fluxes

   do ib=1,nxBdy
     c                    = bdy_exit(1,ib)
     b                    = bdy_exit(2,ib)
     PsiB(b,angle0+Angle) = Set% tPsi(c0+c)
   enddo


   return
   end subroutine SweepGreyUCBxyzKernel



