!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   SweepUCBrz  - This routine calculates angular fluxes for a         *
!                 single direction and multiple energy groups for      *
!                 for an upstream corner-balance (UCB) spatial         *
!                 in rz-geometry.                                      *
!                                                                      *
!***********************************************************************

   subroutine SweepUCBrz(Set, setID, Groups, Angle, savePsi)

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

!  Local

   type(AngleSet),   pointer       :: ASet
   type(GroupSet),   pointer       :: GSet
   type(BdyExit),    pointer       :: BdyExitPtr

   integer    :: b
   integer    :: i
   integer    :: ib
   integer    :: g
   integer    :: cface
   integer    :: cez
   integer    :: cfp
   integer    :: c0

   integer    :: nzones
   integer    :: zone
   integer    :: zone0
   integer    :: nCorner
   integer    :: c
   integer    :: ii
   integer    :: ndoneZ
   integer    :: hyperPlane
   integer    :: nHyperplanes

   integer    :: nxez(Size% maxCorner)
   integer    :: ez_exit(2,Size% maxCorner)

   real(adqt) :: fac
   real(adqt) :: quadwt
   real(adqt) :: quadTauW1
   real(adqt) :: quadTauW2
   real(adqt) :: sigA
   real(adqt) :: sigA2
   real(adqt) :: source
   real(adqt) :: tau

   real(adqt) :: area
   real(adqt) :: sig
   real(adqt) :: sez
   real(adqt) :: gnum
   real(adqt) :: gden

   real(adqt) :: aez
   real(adqt) :: afp
   real(adqt) :: R
   real(adqt) :: R_afp

   real(adqt) :: sumArea(Size% maxCorner)
   real(adqt) :: omega(2)
   real(adqt) :: coefpsi(2,Size% maxCorner)
   real(adqt) :: Q(Groups,Size% maxCorner)
   real(adqt) :: src(Groups,Size% maxCorner)

   real(adqt), parameter :: fouralpha=1.82d0

   logical (kind=1) :: StartingDirection
   logical (kind=1) :: setFinishingDirection


!  Constants

   ASet                  => getAngleSetFromSetID(Quad, setID)
   GSet                  => getGroupSetFromSetID(Quad, setID)

   nHyperPlanes          =  getNumberOfHyperPlanes(ASet,Angle)
   omega(:)              =  ASet% omega(:,Angle)
   quadwt                =  ASet% weight(Angle)
   fac                   =  ASet% angDerivFac(Angle)
   quadTauW1             =  ASet% quadTauW1(Angle)
   quadTauW2             =  ASet% quadTauW2(Angle)
   tau                   =  Size% tau
   StartingDirection     =  ASet% StartingDirection(Angle)
   setFinishingDirection =  ASet% FinishingDirection(Angle+1)

!  Initialize boundary values in Psi1

   do ib=1,Set%nbelem
     Set% Psi1(:,Set%nCorner+ib) = Set% PsiB(:,ib,Angle)
   enddo

!  Loop through all of the corners using the NEXT list

   ndoneZ = 0
   cfp = -1

   HyperPlaneLoop: do hyperPlane=1,nHyperPlanes

     nzones = getZonesInPlane(ASet,Angle,hyperPlane)
 
     ZoneLoop: do ii=1,nzones
 
       zone0   = ASet% nextZ(ndoneZ+ii,Angle)
       zone    = iabs( zone0 ) 
       nCorner = Geom% numCorner(zone)
       c0      = Geom% cOffSet(zone)

       nxez(:) = 0

!      Contributions from volume terms (if a starting direction add angular derivative) 

       do c=1,nCorner 
         do g=1,Groups
           source   = GSet% STotal(g,c0+c) + tau*Set% Psi(g,c0+c,Angle)
           Q(g,c)   = source
           src(g,c) = Geom% Volume(c0+c)*source
         enddo
         sumArea(c) = fac*Geom% Area(c0+c)
       enddo

       CornerLoop: do c=1,nCorner

         CornerFaceLoop: do cface=1,2

           afp = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c0+c) )
           aez = DOT_PRODUCT( omega(:),Geom% A_ez(:,cface,c0+c) )

           if ( afp < zero ) then
             cfp        = Geom% cFP(cface,c0+c)
             R_afp      = Geom% RadiusFP(cface,c0+c)*afp
             sumArea(c) = sumArea(c) - R_afp

             src(:,c)   = src(:,c)   - R_afp*Set% Psi1(:,cfp)
           endif

           if ( aez > zero ) then

             R    = Geom% RadiusEZ(cface,c0+c)
             cez  = Geom% cEZ(cface,c0+c)
             area = Geom% Area(c0+c)

             nxez(c)            = nxez(c) + 1
             ez_exit(nxez(c),c) = cez
             coefpsi(nxez(c),c) = R*aez
             sumArea(cez)       = sumArea(cez) + R*aez
            
             if ( afp < zero ) then

               do g=1,Groups
                 sig        = GSet% Sigt(g,zone)
                 sigA       = sig*area
                 sigA2      = sigA*sigA

                 gnum       = aez*aez*( fouralpha*sigA2 +  &
                              aez*(four*sigA + three*aez) )

                 gden       = area*(four*sigA*sigA2 + aez*(six*sigA2 + &
                              two*aez*(two*sigA + aez)))

                 sez        = R*( area*gnum*( sig*Set% Psi1(g,cfp) - Q(g,c) ) + &
                              half*aez*gden*( Q(g,c) - Q(g,cez) ) )/  &
                              ( gnum + gden*sig )

                 src(g,c)   = src(g,c)   + sez
                 src(g,cez) = src(g,cez) - sez
               enddo

             else

               do g=1,Groups
                 sez        = half*R*aez*( Q(g,c) - Q(g,cez) )/GSet% Sigt(g,zone)
                 src(g,c)   = src(g,c)    + sez
                 src(g,cez) = src(g,cez)  - sez
               enddo

             endif

           endif

         enddo CornerFaceLoop

       enddo CornerLoop


       do i=1,nCorner

         c = ASet% nextC(c0+i,angle) 

!        Corner angular flux
         Set% Psi1(:,c0+c) = (src(:,c) + Geom% Area(c0+c)*fac*    &
                              Set% PsiM(:,c0+c))/( sumArea(c) +  &
                              GSet% Sigt(:,zone)*Geom% Volume(c0+c) )

!        Calculate the contribution of this flux to the sources of
!        downstream corners in this zone. The downstream corner index is
!        "ez_exit."

         do cface=1,nxez(c)
           cez        = ez_exit(cface,c)
           src(:,cez) = src(:,cez) + coefpsi(cface,c)*Set% Psi1(:,c0+c)
         enddo

       enddo

!      Set the "half-angle" angular intensity (PsiM) for the next angle

       if ( StartingDirection ) then
         do c=1,nCorner
           Set% PsiM(:,c0+c) = Set% Psi1(:,c0+c) 
         enddo
       else
         do c=1,nCorner
           Set% PsiM(:,c0+c) = quadTauW1*Set%Psi1(:,c0+c) -  &
                               quadTauW2*Set% PsiM(:,c0+c)

           Set% Phi(:,c0+c)  = Set% Phi(:,c0+c) + quadwt*Set% Psi1(:,c0+c)
         enddo
       endif

     enddo ZoneLoop

     ndoneZ = ndoneZ + nzones

   enddo HyperPlaneLoop

!  Set exiting boundary fluxes

   BdyExitPtr => ASet% BdyExitPtr(Angle)

   do ib=1,BdyExitPtr% nxBdy
     b = BdyExitPtr% bdyList(1,ib)
     c = BdyExitPtr% bdyList(2,ib)

     Set% PsiB(:,b,Angle) = Set% Psi1(:,c)
   enddo

   if ( setFinishingDirection ) then
     do ib=1,BdyExitPtr% nxBdy
       b = BdyExitPtr% bdyList(1,ib)
       c = BdyExitPtr% bdyList(2,ib)

       Set% PsiB(:,b,Angle+1) = Set% PsiM(:,c)
     enddo
   endif

   if ( savePsi ) then
     do c=1,Set% nCorner
       Set% Psi(:,c,Angle) = Set% Psi1(:,c)
     enddo

     if ( setFinishingDirection ) then
       do c=1,Set% nCorner
         Set% Psi(:,c,Angle+1) = Set% PsiM(:,c)
       enddo
     endif
   endif


   return
   end subroutine SweepUCBrz  

