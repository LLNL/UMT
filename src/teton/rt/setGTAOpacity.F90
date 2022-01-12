!***********************************************************************
!                        Version 1:  04/2008, PFN                      *
!                                                                      *
!   setGTAOpacity - calculates grey opacities for grey-transport       *
!                   acceleration (GTA) in 2D/3D or grey-diffusion      *
!                   acceleration (GDA) in 1D.                          *
!                                                                      *
!***********************************************************************
   subroutine setGTAOpacityNEW(zone) 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use Material_mod
   use GreyAcceleration_mod

   implicit none

!  Arguments

   integer, intent(in)     :: zone 

!  Local Variables

   integer    :: g
   integer    :: c
   integer    :: c0
   integer    :: cID
   integer    :: nCorner

   real(adqt) :: tau
   real(adqt) :: greysigt
   real(adqt) :: greysiga
   real(adqt) :: greysigs
   real(adqt) :: SigtInv
   real(adqt) :: SigtInvAve
   real(adqt) :: Sigt2InvAve
   real(adqt) :: SigaAve
   real(adqt) :: SigsAve 

!  Constants

   tau     =  Size%tau
   nCorner =  Geom% numCorner(zone)
   c0      =  Geom% cOffSet(zone)
 
   do c=1,nCorner

     cID = c0 + c

!    Initialize spectral sums
 
     SigtInvAve  = zero
     Sigt2InvAve = zero
     SigaAve     = zero
     SigsAve     = zero

     do g=1,Size%ngr
       SigtInv     = one/(Mat% Siga(g,zone) + Mat% Sigs(g,zone) + tau)
       SigtInvAve  = SigtInvAve   + GTA% Chi(g,cID)*SigtInv
       Sigt2InvAve = Sigt2InvAve  + GTA% Chi(g,cID)*SigtInv*SigtInv
       SigaAve     = SigaAve      + GTA% Chi(g,cID)*Mat% Siga(g,zone)*SigtInv
       SigsAve     = SigsAve      + GTA% Chi(g,cID)*Mat% Sigs(g,zone)*SigtInv

       GTA% Chi(g,cID) = GTA% Chi(g,cID)*SigtInv
     enddo
 
     GTA% Chi(:,cID) = GTA% Chi(:,cID)/SigtInvAve 
 
!    In 2D and 3D, we require slightly different opacities for
!    acceleration.  Compute grey opacities and keep GSIGSC positive
!    or zero

     greysigt = SigtInvAve/Sigt2InvAve
     greysiga = tau + (one - Mat% Eta(cID))*SigaAve/SigtInvAve

     greysigs = greysigt - greysiga
 
!    GSIGSC < 0 is a degenerate case and we reset GSIGTC.

     if (greysigs <= zero) then
       GTA%GreySigScat(cID)  = zero
       GTA%GreySigTotal(cID) = greysiga
     else
       GTA%GreySigScat(cID)  = greysigs
       GTA%GreySigTotal(cID) = greysigt
     endif

     GTA%GreySigScatVol(cID) = GTA%GreySigScat(cID)*Geom% Volume(cID)
     GTA%GreySigtInv(cID)    = one/GTA%GreySigTotal(cID)

   enddo

 
   return
   end subroutine setGTAOpacityNEW 


   subroutine setGTAOpacity(zone)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use Material_mod
   use GreyAcceleration_mod

   implicit none

!  Arguments

   integer, intent(in)     :: zone 

!  Local Variables

   integer    :: g
   integer    :: c
   integer    :: c0
   integer    :: cID
   integer    :: nCorner

   real(adqt) :: tau
   real(adqt) :: scatratio
   real(adqt) :: greysigt
   real(adqt) :: greysiga
   real(adqt) :: greysigs
   real(adqt) :: SigtInv
   real(adqt) :: SigtInvAve
   real(adqt) :: Sigt2InvAve
   real(adqt) :: SigaAve
   real(adqt) :: SigsAve

!  Constants

   tau     =  Size%tau
   nCorner =  Geom% numCorner(zone)
   c0      =  Geom% cOffSet(zone)

   do c=1,nCorner

     cID = c0 + c

!    Initialize spectral sums

     SigtInvAve  = zero
     Sigt2InvAve = zero
     SigaAve     = zero
     SigsAve     = zero

     do g=1,Size%ngr
       SigtInv     = one/(Mat% Siga(g,zone) + Mat% Sigs(g,zone) + tau)
       SigtInvAve  = SigtInvAve   + GTA% Chi(g,cID)*SigtInv
       Sigt2InvAve = Sigt2InvAve  + GTA% Chi(g,cID)*SigtInv*SigtInv
       SigaAve     = SigaAve      + GTA% Chi(g,cID)*Mat% Siga(g,zone)*SigtInv
       SigsAve     = SigsAve      + GTA% Chi(g,cID)*Mat% Sigs(g,zone)*SigtInv

       GTA% Chi(g,cID) = GTA% Chi(g,cID)*SigtInv
     enddo

     GTA% Chi(:,cID) = GTA% Chi(:,cID)/SigtInvAve

!  Compute grey opacities

     if (Size%ndim == 1) then

!    In 1D we need a diffusion coefficient

       GTA%GreyDiffCoef(cID) = Sigt2InvAve/(three*SigtInvAve)
       GTA%GreySigEff(cID)   = (one - Mat% Eta(cID)*SigaAve - SigsAve)/SigtInvAve

     elseif (Size%ndim >= 2) then

!    In 2D and 3D, we require slightly different opacities for
!    acceleration.  Compute grey opacities and keep GSIGSC positive
!    or zero

       greysigt = SigtInvAve/Sigt2InvAve
       greysiga = tau + (one - Mat% Eta(cID))*SigaAve/SigtInvAve

       greysigs = greysigt - greysiga

!    GSIGSC < 0 is a degenerate case and we reset GSIGTC.
!    Compute the grey opacities for the stretched-TSA step.
!   "Epsilon" is chosen so that the scattering term vanishes
!    in the stretched problem.

       if (greysigs <= zero) then
         GTA%GreySigScat(cID)     = zero
         GTA%GreySigTotal2(cID,1) = greysiga
       else
         GTA%GreySigScat(cID)     = greysigs
         GTA%GreySigTotal2(cID,1) = greysigt
       endif

       scatratio                = GTA%GreySigScat(cID)/GTA%GreySigTotal2(cID,1)
       GTA%eps(cID)             = one/sqrt(one - scatratio)
       GTA%GreySigTotal2(cID,2) = GTA%GreySigTotal2(cID,1)/GTA%eps(cID)
       GTA%GreySigScatVol(cID)  = GTA%GreySigScat(cID)*Geom% Volume(cID)

       GTA%GreySigtInv2(cID,1)  = one/GTA%GreySigTotal2(cID,1)
       GTA%GreySigtInv2(cID,2)  = one/GTA%GreySigTotal2(cID,2)

     endif

   enddo


   return
   end subroutine setGTAOpacity
