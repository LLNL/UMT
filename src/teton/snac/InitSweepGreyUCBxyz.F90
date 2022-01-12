!***********************************************************************
!                        Last Update:  09/2018, PFN                    *
!                                                                      *
!   InitGreySweepUCBxyz - This routine calculates the transfer         *
!                         matrices used for the direct solve for the   *
!                         scalar corrections.                          *
!                                                                      *
!***********************************************************************

   subroutine InitGreySweepUCBxyz(zone)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: zone 

!  Local Variables

   type(AngleSet),   pointer :: ASet

   integer, dimension (1)    :: cfirst

   integer    :: nSets
   integer    :: nGTASets
   integer    :: setID
   integer    :: angle
   integer    :: numAngles
   integer    :: i
   integer    :: k
   integer    :: cface
   integer    :: ifp

   integer    :: c
   integer    :: c0 
   integer    :: c1
   integer    :: cez 
   integer    :: nCorner 
   integer    :: nCFaces

   integer    :: nxez(Size% maxCorner)
   integer    :: need(Size% maxCorner)
   integer    :: ez_exit(Size%maxcf,Size% maxCorner)

   real(adqt) :: aez
   real(adqt) :: area_opp
   real(adqt) :: sigv
   real(adqt) :: sigv2
   real(adqt) :: gnum
   real(adqt) :: gtau
   real(adqt) :: B0
   real(adqt) :: B1
   real(adqt) :: B2
   real(adqt) :: coef
   real(adqt) :: dInv
   real(adqt) :: quadwt

   real(adqt) :: omega(3)
   real(adqt) :: denom(Size% maxCorner)
   real(adqt) :: afp(Size%maxcf)
   real(adqt) :: coefpsi(Size%maxcf,Size% maxCorner)
   real(adqt) :: Sigt(Size% maxCorner)
   real(adqt) :: Pvv(Size% maxCorner,Size% maxCorner)

   real(adqt), parameter :: fouralpha=1.82_adqt

!  Constants

   nSets    = getNumberOfSets(Quad)
   nGTASets = getNumberOfGTASets(Quad)
   nCorner  = Geom% numCorner(zone)
   c0       = Geom% cOffSet(zone)

   do c=1,nCorner
     GTA% TT(:,c0+c) = zero
   enddo

   GTASetLoop: do setID=nSets+1,nSets+nGTASets

     ASet      => getAngleSetData(Quad, setID)
     numAngles =  ASet% numAngles

     AngleLoop: do Angle=1,numAngles

       omega(:) = ASet% omega(:,Angle)
       quadwt   = ASet% weight(Angle)

       nxez(:)  = 0
       need(:)  = 0
       Pvv(:,:) = zero

       do c=1,nCorner
         Pvv(c,c) = Geom% Volume(c0+c)
         Sigt(c)  = GTA%GreySigTotal(c0+c)
       enddo

!      Loop over corners 

       CornerLoop: do c=1,nCorner

         sigv     = Geom% Volume(c0+c)*Sigt(c) 
         denom(c) = sigv 
         nCFaces  = Geom% nCFacesArray(c0+c)

!        Contributions from external corner faces (FP faces)

         do cface=1,nCFaces

           afp(cface) = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c0+c) )

           if ( afp(cface) > zero ) then
             denom(c) = denom(c) + afp(cface)
           endif
         enddo

!        Contributions from interior corner faces (EZ faces)

         do cface=1,ncfaces

           aez = DOT_PRODUCT( omega(:),Geom% A_ez(:,cface,c0+c) )
           cez = Geom% cEZ(cface,c0+c)

           if (cez > c) then
             if (aez > zero ) then
               need(cez)              = need(cez) + 1
               nxez(c)                = nxez(c)   + 1
               ez_exit(nxez(c),c)     = cez
               coefpsi(nxez(c),c)     = aez
             elseif (aez < zero) then
               need(c)                = need(c)   + 1
               nxez(cez)              = nxez(cez) + 1
               ez_exit(nxez(cez),cez) = c
               coefpsi(nxez(cez),cez) = -aez
             endif
           endif

           if (aez > zero ) then

             denom(c) = denom(c) + aez
             area_opp = zero

             if (nCFaces == 3) then

               ifp = mod(cface,nCFaces) + 1

               if ( afp(ifp) < zero ) then
                 area_opp = -afp(ifp)
               endif

             else

               ifp      = cface
               area_opp = zero

               do k=1,nCFaces-2
                 ifp = mod(ifp,nCFaces) + 1
                 if ( afp(ifp) < zero ) then
                   area_opp = area_opp - afp(ifp)
                 endif
               enddo

             endif

             TestOppositeFace: if (area_opp > zero) then

               sigv2     = sigv*sigv

               gnum      = aez*aez*( fouralpha*sigv2 +    &
                           aez*(four*sigv + three*aez) )

               gtau      = gnum/    &
                         ( gnum + four*sigv2*sigv2 + aez*sigv*(six*sigv2 + &
                           two*aez*(two*sigv + aez)) )

               B0        = half*aez*(one - gtau)
               B1        = (B0 - gtau*sigv)/Sigt(c) 
               B2        =  B0/Sigt(cez) 

!              Pvv(column,row)
               Pvv(c,c)     = Pvv(c,c)     + B1
               Pvv(cez,c)   = Pvv(cez,c)   - B2
               Pvv(c,cez)   = Pvv(c,cez)   - B1
               Pvv(cez,cez) = Pvv(cez,cez) + B2

             else

               B1           = half*aez/Sigt(c)
               B2           = half*aez/Sigt(cez)

!              Pvv(column,row)
               Pvv(c,c)     = Pvv(c,c)     + B1
               Pvv(cez,c)   = Pvv(cez,c)   - B2
               Pvv(c,cez)   = Pvv(c,cez)   - B1
               Pvv(cez,cez) = Pvv(cez,cez) + B2

             endif TestOppositeFace

           endif

         enddo

       enddo CornerLoop

       do i=1,nCorner

         cfirst = minloc( need(1:nCorner) )
         c      = cfirst(1)

         dInv = one/denom(c)

         do c1=1,nCorner
           Pvv(c1,c) = dInv*Pvv(c1,c)
         enddo

!        Calculate the contribution of this flux to the sources of
!        downstream corners in this zone. The downstream corner index is
!        "ez_exit."

         do cface=1,nxez(c)
           cez       = ez_exit(cface,c)
           coef      = coefpsi(cface,c)
           need(cez) = need(cez) - 1

           do c1=1,nCorner
             Pvv(c1,cez) = Pvv(c1,cez) + coef*Pvv(c1,c)
           enddo
         enddo

         need(c) = 99

       enddo

       do c1=1,nCorner
         do c=1,nCorner
           GTA% TT(c,c0+c1) = GTA% TT(c,c0+c1) + quadwt*Pvv(c,c1)
         enddo
       enddo

     enddo AngleLoop

   enddo GTASetLoop


   return
   end subroutine InitGreySweepUCBxyz 


