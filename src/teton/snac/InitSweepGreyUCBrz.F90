!***********************************************************************
!                        Last Update:  09/2018, PFN                    *
!                                                                      *
!   InitGreySweepUCBrz  - This routine calculates the transfer         *
!                         matrices used for the direct solve for the   *
!                         scalar corrections.                          *
!                                                                      *
!***********************************************************************

   subroutine InitGreySweepUCBrz(zone)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,        intent(in) :: zone

!  Local

   type(AngleSet), pointer    :: ASet

   integer    :: nSets
   integer    :: nGTASets
   integer    :: setID
   integer    :: i
   integer    :: cface
   integer    :: cez
   integer    :: c0
   integer    :: nCorner
   integer    :: c
   integer    :: c1
   integer    :: Angle
   integer    :: numAngles

   integer    :: nxez(Size% maxCorner)
   integer    :: ez_exit(2,Size% maxCorner)

   real(adqt) :: fac
   real(adqt) :: quadwt
   real(adqt) :: sigA
   real(adqt) :: sigA2
   real(adqt) :: gnum
   real(adqt) :: gtau
   real(adqt) :: quadTauW1
   real(adqt) :: quadTauW2
   real(adqt) :: afp
   real(adqt) :: R_afp
   real(adqt) :: aez
   real(adqt) :: R
   real(adqt) :: dInv
   real(adqt) :: B0 
   real(adqt) :: B1
   real(adqt) :: B2
   real(adqt) :: coef

   real(adqt) :: omega(2)
   real(adqt) :: denom(Size% maxCorner)
   real(adqt) :: coefpsi(2,Size% maxCorner)
   real(adqt) :: Sigt(Size% maxCorner)

   real(adqt) :: Tvv(Size% maxCorner,Size% maxCorner)
   real(adqt) :: Pvv(Size% maxCorner,Size% maxCorner)

   real(adqt), parameter :: fouralpha=1.82_adqt

   logical (kind=1) :: StartingDirection
   logical (kind=1) :: FinishingDirection

!  Constants

   nSets    = getNumberOfSets(Quad)
   nGTASets = getNumberOfGTASets(Quad)
   nCorner  = Geom% numCorner(zone)
   c0       = Geom% cOffSet(zone)

   Tvv(:,:) = zero

   do c=1,nCorner
     GTA% TT(:,c0+c) = zero
     Sigt(c)         = GTA%GreySigTotal(c0+c)
   enddo

   GTASetLoop: do setID=nSets+1,nSets+nGTASets

     ASet      => getAngleSetData(Quad, setID)
     numAngles =  ASet% numAngles

     AngleLoop: do Angle=1,numAngles

       StartingDirection  = ASet% StartingDirection(Angle)
       FinishingDirection = ASet% FinishingDirection(Angle)

       omega(:)           = ASet% omega(:,Angle)
       quadwt             = ASet% weight(Angle)

       if ( .not. FinishingDirection ) then

         fac       = ASet% angDerivFac(Angle)
         quadTauW1 = ASet% quadTauW1(Angle)
         quadTauW2 = ASet% quadTauW2(Angle)

         nxez(:)  = 0
         Pvv(:,:) = zero

!        Contributions from volume terms (if a starting direction add angular derivative) 

         do c=1,nCorner 
           Pvv(c,c) = Geom% Volume(c0+c)
           denom(c) = Sigt(c)*Geom% Volume(c0+c) + fac*Geom% Area(c0+c)

           do c1=1,nCorner
             Pvv(c1,c) = Pvv(c1,c) + fac*Geom% Area(c0+c)*Tvv(c1,c)
           enddo
         enddo

         CornerLoop: do c=1,nCorner

           CornerFaceLoop: do cface=1,2

             afp = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c0+c) )
             aez = DOT_PRODUCT( omega(:),Geom% A_ez(:,cface,c0+c) )

             if ( afp < zero ) then
               R_afp    = Geom% RadiusFP(cface,c0+c)*afp
               denom(c) = denom(c) - R_afp
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

                 B0           = half*aez*(one - gtau)*R
                 B1           = (B0 - R*gtau*sigA)/Sigt(c)
                 B2           =  B0/Sigt(cez)

!                Pvv(column,row)
                 Pvv(c,c)     = Pvv(c,c)     + B1
                 Pvv(cez,c)   = Pvv(cez,c)   - B2
                 Pvv(c,cez)   = Pvv(c,cez)   - B1
                 Pvv(cez,cez) = Pvv(cez,cez) + B2

               else

                 B1           = half*R*aez/Sigt(c)
                 B2           = half*R*aez/Sigt(cez)

!                Pvv(column,row)
                 Pvv(c,c)     = Pvv(c,c)     + B1
                 Pvv(cez,c)   = Pvv(cez,c)   - B2
                 Pvv(c,cez)   = Pvv(c,cez)   - B1
                 Pvv(cez,cez) = Pvv(cez,cez) + B2

               endif

             endif

           enddo CornerFaceLoop

         enddo CornerLoop

         do i=1,nCorner

           c    = ASet% nextC(c0+i,angle)
           dInv = one/denom(c)

           do c1=1,nCorner
             Pvv(c1,c) = dInv*Pvv(c1,c)
           enddo

!          Calculate the contribution of this flux to the sources of
!          downstream corners in this zone. The downstream corner index is
!          "ez_exit."

           do cface=1,nxez(c)
             cez       = ez_exit(cface,c)
             coef      = coefpsi(cface,c)

             do c1=1,nCorner
               Pvv(c1,cez) = Pvv(c1,cez) + coef*Pvv(c1,c)
             enddo

           enddo

         enddo

!        Increment matrix elements for this angle 

         if ( StartingDirection ) then
           Tvv(:,:) = Pvv(:,:)
         else
           do c=1,nCorner
             do c1=1,nCorner
               GTA% TT(c1,c0+c) = GTA% TT(c1,c0+c)    + quadwt*Pvv(c1,c)
               Tvv(c1,c)        = quadTauW1*Pvv(c1,c) - quadTauW2*Tvv(c1,c)
             enddo
           enddo
         endif

       endif

     enddo AngleLoop

   enddo GTASetLoop



   return
   end subroutine InitGreySweepUCBrz  

