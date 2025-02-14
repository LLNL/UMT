#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  09/2018, PFN                    *
!                                                                      *
!   InitGreySweepUCBrz  - This routine calculates the transfer         *
!                         matrices used for the direct solve for the   *
!                         scalar corrections.                          *
!                                                                      *
!***********************************************************************

   subroutine InitGreySweepUCBrz_GPU

   use, intrinsic :: iso_c_binding, only : c_int
   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use AngleSet_mod
   use OMPWrappers_mod
   use CodeChecks_mod

   implicit none

!  Local

   type(AngleSet),    pointer :: ASet

   integer    :: nAngleSets
   integer    :: nGTASets
   integer    :: nZoneSets
   integer    :: zSetID
   integer    :: setID
   integer    :: zone
   integer    :: i
   integer    :: cface
   integer    :: cez
   integer    :: c0
   integer    :: nCorner
   integer    :: c
   integer    :: c1
   integer    :: angle
   integer    :: aSetID
   integer    :: angGTA
   integer    :: n
   integer    :: numAngles

   real(adqt) :: sigA
   real(adqt) :: sigA2
   real(adqt) :: gnum
   real(adqt) :: gtau
   real(adqt) :: afp
   real(adqt) :: R_afp
   real(adqt) :: R_afp2
   real(adqt) :: aez
   real(adqt) :: R
   real(adqt) :: R2
   real(adqt) :: dInv
   real(adqt) :: B0 
   real(adqt) :: B1
   real(adqt) :: B2
   real(adqt) :: quadwt
   real(adqt) :: Sigt
   real(adqt) :: SigtEZ

   real(adqt), parameter :: fouralpha=1.82_adqt

!  Dynamic

   integer,         allocatable :: angleList(:,:)
   real(adqt),      allocatable :: omega(:,:)
   real(adqt),      allocatable :: fac(:)
   real(adqt),      allocatable :: quadTauW1(:)
   real(adqt),      allocatable :: quadTauW2(:)
   logical(kind=1), allocatable :: Starting(:)

!  Constants
   nZoneSets  =  getNumberOfZoneSets(Quad)
   nAngleSets =  getNumberOfAngleSets(Quad)
   nGTASets   =  getNumberOfGTASets(Quad)

!  The quadrature weight is a constant

   ASet   => getAngleSetData(Quad, nAngleSets+1)
   quadwt =  ASet% weight(2)

   numAngles = 0
   do setID=1,nGTASets
     ASet => getAngleSetData(Quad, nAngleSets+setID)

     do angle=1,ASet% numAngles
       if ( .not. ASet% FinishingDirection(angle) ) then
         numAngles = numAngles + 1 
       endif
     enddo
   enddo

   allocate( angleList(2,numAngles) )
   allocate( omega(2,numAngles) )
   allocate( fac(numAngles) )
   allocate( quadTauW1(numAngles) )
   allocate( quadTauW2(numAngles) )
   allocate( Starting(numAngles) )

   angle = 0

   do setID=1,nGTASets
     ASet => getAngleSetData(Quad, nAngleSets+setID)

     do n=1,ASet% numAngles
       if ( .not. ASet% FinishingDirection(n) ) then
         angle              = angle + 1
         angleList(1,angle) = nAngleSets + setID
         angleList(2,angle) = n 
         omega(1:2,angle)   = ASet% omega(1:2,n)
         fac(angle)         = ASet% angDerivFac(n)
         quadTauW1(angle)   = ASet% quadTauW1(n)
         quadTauW2(angle)   = ASet% quadTauW2(n)
         Starting(angle)    = ASet% StartingDirection(n)
       endif
     enddo
   enddo


     ! Verify we won't get out-of-bounds access in angle loop below.
     TETON_CHECK_BOUNDS1(Geom%corner1, nZoneSets)
     TETON_CHECK_BOUNDS1(Geom%corner2, nZoneSets)

     TOMP_MAP(target enter data map(to: numAngles, angleList, omega, quadwt, fac, quadTauW1, quadTauW2, Starting))

     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nZoneSets, numAngles, omega, Geom, GTA) &)
     TOMPC(private(angle))

     ZoneSetLoop1: do zSetID=1,nZoneSets

       do angle=1,numAngles

         !$omp parallel do collapse(2) default(none)  &
         !$omp& shared(Geom, GTA, angle, omega, zSetID)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           do cface=1,2
             GTA% AfpNorm(cface,c,angle) = DOT_PRODUCT( omega(:,angle),Geom% A_fp(:,cface,c) )
             GTA% AezNorm(cface,c,angle) = DOT_PRODUCT( omega(:,angle),Geom% A_ez(:,cface,c) )
           enddo
         enddo

         !$omp end parallel do

       enddo

     enddo ZoneSetLoop1

     TOMP(end target teams distribute)


     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nZoneSets, numAngles, GTA, Geom, fac)&)
     TOMPC(private(R_afp, R_afp2, R, R2, angle))

     ZoneSetLoop2: do zSetID=1,nZoneSets

       do angle=1,numAngles

         !$omp  parallel do default(none)  &
         !$omp& shared(angle, Geom, GTA, fac, zSetID) private(R_afp,R_afp2,R,R2)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           R_afp  = Geom% RadiusFP(1,c)
           R_afp2 = Geom% RadiusFP(2,c)
           R      = Geom% RadiusEZ(1,c)
           R2     = Geom% RadiusEZ(2,c)

           GTA% ANormSum(c,angle) = fac(angle)*Geom% Area(c) - half*(  &
              R_afp *(GTA% AfpNorm(1,c,angle) - abs(GTA% AfpNorm(1,c,angle))) +  &
              R_afp2*(GTA% AfpNorm(2,c,angle) - abs(GTA% AfpNorm(2,c,angle))) +  &
              R     *(GTA% AezNorm(1,c,angle) - abs(GTA% AezNorm(1,c,angle))) +  &
              R2    *(GTA% AezNorm(2,c,angle) - abs(GTA% AezNorm(2,c,angle))) )
         enddo

         !$omp end parallel do

       enddo

     enddo ZoneSetLoop2

     TOMP(end target teams distribute)


     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
     TOMPC(shared(nZoneSets, Geom, GTA, Quad, numAngles, quadwt, angleList, quadTauW1, quadTauW2, fac, Starting)  &)
     TOMPC(private(aSetID,angGTA,c0,cez,nCorner,R,aez,afp,sigA,sigA2,gnum,gtau,B0,B1,B2,dInv,Sigt,SigtEZ))

     ZoneSetLoop: do zSetID=1,nZoneSets

     !$omp  parallel do default(none)  &
     !$omp& shared(Geom, GTA, Quad, numAngles, quadwt, zSetID)      &
     !$omp& shared(angleList, quadTauW1, quadTauW2, fac, Starting)  &
     !$omp& private(aSetID,angGTA,c0,cez,nCorner,R,aez,afp)         &
     !$omp& private(sigA,sigA2,gnum,gtau,B0,B1,B2,dInv,Sigt,SigtEZ, c1, c, cface)

       ZoneLoop: do zone=Geom% zone1(zSetID),Geom% zone2(zSetID)

         nCorner = Geom% numCorner(zone)
         c0      = Geom% cOffSet(zone)

         do c1=1,nCorner
           do c=1,nCorner
             GTA% TT(c,c0+c1)  = zero
             GTA% Tvv(c,c0+c1) = zero
           enddo
         enddo

         AngleLoop: do angle=1,numAngles

           aSetID = angleList(1,angle)
           angGTA = angleList(2,angle)

           do c=1,nCorner
             do c1=1,nCorner
               GTA% Pvv(c1,c0+c) = zero
             enddo
             GTA% Pvv(c,c0+c) = Geom% Volume(c0+c)
           enddo

!          Contributions from angular derivative 

           do c=1,nCorner 
             do c1=1,nCorner
               GTA% Pvv(c1,c0+c) = GTA% Pvv(c1,c0+c) + fac(angle)*  &
                                   Geom% Area(c0+c)*GTA% Tvv(c1,c0+c)
             enddo
           enddo

           CornerLoop: do c=1,nCorner

             Sigt = GTA%GreySigTotal(c0+c)
             sigA = Sigt*Geom% Area(c0+c)

             CornerFaceLoop: do cface=1,2

               afp = GTA% AfpNorm(cface,c0+c,angle) 
               aez = GTA% AezNorm(cface,c0+c,angle) 

               if ( aez > zero ) then

                 R      = Geom% RadiusEZ(cface,c0+c)
                 cez    = Geom% cEZ(cface,c0+c)
                 SigtEZ = GTA%GreySigTotal(c0+cez)

                 if ( afp < zero ) then

                   sigA2    = sigA*sigA

                   gnum     = aez*aez*( fouralpha*sigA2 +  &
                              aez*(four*sigA + three*aez) )

                   gtau      = gnum/                                           &
                             ( gnum + four*sigA2*sigA2 + aez*sigA*(six*sigA2 + &
                               two*aez*(two*sigA + aez)) )

                   B0        = half*aez*(one - gtau)*R
                   B1        = (B0 - R*gtau*sigA)/Sigt
                   B2        =  B0/SigtEZ

                   ! Pvv(column,row)
                   GTA% Pvv(c,c0+c)     = GTA% Pvv(c,c0+c)     + B1
                   GTA% Pvv(cez,c0+c)   = GTA% Pvv(cez,c0+c)   - B2
                   GTA% Pvv(c,c0+cez)   = GTA% Pvv(c,c0+cez)   - B1
                   GTA% Pvv(cez,c0+cez) = GTA% Pvv(cez,c0+cez) + B2

                 else

                   B1        = half*R*aez/Sigt
                   B2        = half*R*aez/SigtEZ

                   ! Pvv(column,row)
                   GTA% Pvv(c,c0+c)     = GTA% Pvv(c,c0+c)     + B1
                   GTA% Pvv(cez,c0+c)   = GTA% Pvv(cez,c0+c)   - B2
                   GTA% Pvv(c,c0+cez)   = GTA% Pvv(c,c0+cez)   - B1
                   GTA% Pvv(cez,c0+cez) = GTA% Pvv(cez,c0+cez) + B2

                 endif

               endif

             enddo CornerFaceLoop

           enddo CornerLoop

           do i=1,nCorner

             c    = Quad% AngSetPtr(aSetID)% nextC(c0+i,angGTA)

             dInv = one/(GTA% ANormSum(c0+c,angle) +   &
                         Geom% Volume(c0+c)*GTA%GreySigTotal(c0+c))

             do c1=1,nCorner
               GTA% Pvv(c1,c0+c) = dInv*GTA% Pvv(c1,c0+c)
             enddo

!            Calculate the contribution of this flux to the sources of
!            downstream corners in this zone.

             do cface=1,2
               aez = GTA% AezNorm(cface,c0+c,angle)

               if (aez > zero) then
                 R   = Geom% RadiusEZ(cface,c0+c)
                 cez = Geom% cEZ(cface,c0+c)

                 do c1=1,nCorner
                   GTA% Pvv(c1,c0+cez) = GTA% Pvv(c1,c0+cez) + R*aez*GTA% Pvv(c1,c0+c)
                 enddo
               endif

             enddo

           enddo

!          Increment matrix elements for this angle 

           if ( Starting(angle) ) then
             do c=1,nCorner
               do c1=1,nCorner
                 GTA% Tvv(c1,c0+c) = GTA% Pvv(c1,c0+c)
               enddo
             enddo
           else
             do c=1,nCorner
               do c1=1,nCorner
                 GTA% TT(c1,c0+c) = GTA% TT(c1,c0+c)  +  &
                                    quadwt*GTA% Pvv(c1,c0+c)
                 GTA% Tvv(c1,c0+c)= quadTauW1(angle)*GTA% Pvv(c1,c0+c) -  &
                                    quadTauW2(angle)*GTA% Tvv(c1,c0+c)
               enddo
             enddo
           endif

         enddo AngleLoop

       enddo ZoneLoop

       !$omp end parallel do

     enddo ZoneSetLoop

     TOMP(end target teams distribute)

     TOMP_MAP(target exit data map(release: numAngles, angleList, omega, quadwt, fac, quadTauW1, quadTauW2, Starting))


   deallocate( angleList )
   deallocate( omega )
   deallocate( fac )
   deallocate( quadTauW1 )
   deallocate( quadTauW2 )
   deallocate( Starting )


   return
   end subroutine InitGreySweepUCBrz_GPU

