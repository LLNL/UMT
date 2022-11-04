!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   SweepUCBxyz  - This routine calculates angular fluxes for a        *
!                 single direction and multiple energy groups for      *
!                 for an upstream corner-balance (UCB) spatial         *
!                 in xyz-geometry.                                     *
!                                                                      *
!***********************************************************************

   subroutine SweepUCBxyz(Set, setID, Groups, Angle, savePsi) 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use GroupSet_mod

   implicit none

!  Arguments

   type(SetData),    intent(inout) :: Set
   integer,          intent(in)    :: setID
   integer,          intent(in)    :: Groups 
   integer,          intent(in)    :: Angle 

   logical (kind=1), intent(in)    :: savePsi

!  Local Variables

   type(AngleSet),   pointer       :: ASet
   type(GroupSet),   pointer       :: GSet

   integer    :: b, i, ib, cfp, cface, ifp, g, k
   integer    :: nzones, zone, zone0, c, c0, cez, ii
   integer    :: nCorner, nCFaces
   integer    :: nxBdy
   integer    :: ndoneZ
   integer    :: hyperPlane
   integer    :: nHyperplanes

   integer    :: nxez(Size% maxCorner)
   integer    :: ez_exit(Size%maxcf,Size% maxCorner)
   integer    :: bdy_exit(2,Size%maxcf*Size% maxCorner)

   real(adqt) :: fouralpha 
   real(adqt) :: aez 
   real(adqt) :: aez2 
   real(adqt) :: area_opp 
   real(adqt) :: area_inv

   real(adqt) :: source 
   real(adqt) :: sig
   real(adqt) :: vol
   real(adqt) :: sigv 
   real(adqt) :: sigv2 
   real(adqt) :: gnum 
   real(adqt) :: gden 
   real(adqt) :: sez
   real(adqt) :: quadwt
   real(adqt) :: tau

   real(adqt) :: sigInv
   real(adqt) :: sumArea(Size% maxCorner)
   real(adqt) :: omega(3)
   real(adqt) :: psi_opp(Groups)
   real(adqt) :: SigtVol(Groups,Size% maxCorner)
   real(adqt) :: src(Groups,Size% maxCorner)
   real(adqt) :: Q(Groups,Size% maxCorner)
   real(adqt) :: afp(Size%maxcf)
   real(adqt) :: coefpsi(Size%maxcf,Size% maxCorner)
   real(adqt) :: psifp(Groups,Size%maxcf)

!  Constants

   parameter (fouralpha=1.82d0)

   ASet         => getAngleSetFromSetID(Quad, setID)
   GSet         => getGroupSetFromSetID(Quad, setID)

   nHyperPlanes =  getNumberOfHyperPlanes(ASet,Angle)
   omega(:)     =  ASet% omega(:,Angle)
   quadwt       =  ASet% weight(Angle)
   tau          =  Size% tau

   psi_opp(:) = zero
   psifp(:,:) = zero

!  Initialize boundary values in Psi1

   do ib=1,Set%nbelem
     Set% Psi1(:,Set%nCorner+ib) = Set% PsiB(:,ib,Angle)
   enddo

!  Track the location of the start of the hyperplane zone in nextZ list.

   ndoneZ = 0

   HyperPlaneLoop: do hyperPlane=1,nHyperPlanes

     nzones = getZonesInPlane(ASet,Angle,hyperPlane)

     ZoneLoop: do ii=1,nzones
 
       zone0   = ASet% nextZ(ndoneZ+ii,Angle)

       zone    = iabs( zone0 )
       nCorner = Geom% numCorner(zone)
       c0      = Geom% cOffSet(zone)

       nxBdy   = 0

!      Contributions from volume terms

       do c=1,nCorner
         do g=1,Groups
           source       = GSet% STotal(g,c0+c) + tau*Set% Psi(g,c0+c,Angle)
           Q(g,c)       = source 
           src(g,c)     = Geom% Volume(c0+c)*source
           SigtVol(g,c) = GSet% Sigt(g,zone)*Geom% Volume(c0+c)
         enddo
       enddo

       nxez(:) = 0

       CornerLoop: do c=1,nCorner

!      Calculate Area_CornerFace dot Omega to determine the 
!      contributions from incident fluxes across external 
!      corner faces (FP faces)

         sumArea(c) = zero
         nCFaces    = Geom% nCFacesArray(c0+c)
 
         do cface=1,nCFaces

           afp(cface) = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c0+c) )

           cfp = Geom% cFP(cface,c0+c)

           if ( afp(cface) > zero ) then

             sumArea(c)    = sumArea(c) + afp(cface)

             if (cfp > Set% nCorner) then
               nxBdy             = nxBdy + 1
               bdy_exit(1,nxBdy) = c
               bdy_exit(2,nxBdy) = cfp - Set% nCorner
             endif

           elseif ( afp(cface) < zero ) then

             psifp(:,cface) = Set% Psi1(:,cfp)
             src(:,c)       = src(:,c) - afp(cface)*Set% Psi1(:,cfp)

           endif
         enddo

!        Contributions from interior corner faces (EZ faces)

         do cface=1,nCFaces

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

             sumArea(c) = sumArea(c) + aez
             area_opp   = zero

             if (nCFaces == 3) then

               ifp = mod(cface,nCFaces) + 1

               if ( afp(ifp) < zero ) then
                 psi_opp(:) =  psifp(:,ifp)
                 area_opp   = -afp(ifp)
               endif

             else

               ifp        = cface
               area_opp   = zero
               psi_opp(:) = zero

               do k=1,nCFaces-2
                 ifp = mod(ifp,nCFaces) + 1
                 if ( afp(ifp) < zero ) then
                   area_opp   = area_opp   - afp(ifp)
                   psi_opp(:) = psi_opp(:) - afp(ifp)*psifp(:,ifp)
                 endif
               enddo

               if (area_opp > zero) then
                 area_inv   = one/area_opp
                 psi_opp(:) = psi_opp(:)*area_inv
               endif

             endif

             TestOppositeFace: if ( area_opp > zero ) then

               aez2 = aez*aez
               vol  = Geom% Volume(c0+c)

               do g=1,Groups
  
                 sig        = GSet% Sigt(g,zone)
                 sigv       = sig*vol 
                 sigv2      = sigv*sigv
 
                 gnum       = aez2*( fouralpha*sigv2 +              &
                              aez*(four*sigv + three*aez) )

                 gden       = vol*( four*sigv*sigv2 + aez*(six*sigv2 + &
                              two*aez*(two*sigv + aez)) )

                 sez        = ( vol*gnum*( sig*psi_opp(g) - Q(g,c) ) +   &
                                half*aez*gden*( Q(g,c) - Q(g,cez) ) )/ &
                              ( gnum + gden*sig)

                 src(g,c)   = src(g,c)   + sez
                 src(g,cez) = src(g,cez) - sez

               enddo

             else

               do g=1,Groups
                 sigInv     = one/GSet% Sigt(g,zone)
                 sez        = half*aez*sigInv*( Q(g,c) - Q(g,cez) )
                 src(g,c)   = src(g,c)   + sez
                 src(g,cez) = src(g,cez) - sez
               enddo

             endif TestOppositeFace

           endif

         enddo

       enddo CornerLoop


       if ( zone0 > 0 ) then
          
         do i=1,nCorner
           c = ASet% nextC(c0+i,angle)

!          Corner angular flux
           Set% Psi1(:,c0+c) = src(:,c)/(sumArea(c) + SigtVol(:,c))

!          Scalar Flux
           Set% Phi(:,c0+c)  = Set% Phi(:,c0+c) + quadwt*Set% Psi1(:,c0+c)

!          Calculate the contribution of this flux to the sources of
!          downstream corners in this zone. The downstream corner index is
!          "ez_exit."

           do cface=1,nxez(c)
             cez        = ez_exit(cface,c)
             src(:,cez) = src(:,cez) + coefpsi(cface,c)*Set% Psi1(:,c0+c)
           enddo

         enddo
        
       else
        
!      Direct Solve (non-lower triangular, use old values of Psi1 for now)
         do c=1,nCorner
           do cface=1,nxez(c)
             cez        = ez_exit(cface,c)
             src(:,cez) = src(:,cez) + coefpsi(cface,c)*Set% Psi1(:,c0+c)
           enddo
         enddo

         do c=1,nCorner
           Set% Psi1(:,c0+c) = src(:,c)/(sumArea(c) + SigtVol(:,c))
           Set% Phi(:,c0+c)  = Set% Phi(:,c0+c) + quadwt*Set% Psi1(:,c0+c)
         enddo
          
       endif

!      Set exiting boundary fluxes

       do ib=1,nxBdy
         c                    = bdy_exit(1,ib)
         b                    = bdy_exit(2,ib)
         Set% PsiB(:,b,Angle) = Set% Psi1(:,c0+c)
       enddo

     enddo ZoneLoop

     ndoneZ = ndoneZ + nzones

   enddo HyperPlaneLoop

   if ( savePsi ) then
     do c=1,Set% nCorner
       Set% Psi(:,c,Angle) = Set% Psi1(:,c)
     enddo
   endif


   return
   end subroutine SweepUCBxyz 

