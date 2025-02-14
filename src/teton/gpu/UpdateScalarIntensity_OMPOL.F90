#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  09/2018, PFN                    *
!                                                                      *
!   UpdateScalarIntensity:                                             *
!                                                                      *
!   Invert the transport operator to solve for an updated value of     *
!   the scalar intensity.                                              *
!                                                                      *
!                                                                      *
!***********************************************************************
   subroutine ScalarIntensityDecompose_GPU

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use GreyAcceleration_mod
   use QuadratureList_mod
   use OMPWrappers_mod

   implicit none

!  Arguments

!  Local

   integer    :: zSetID
   integer    :: nZoneSets
   integer    :: zone
   integer    :: c
   integer    :: cc 
   integer    :: c0
   integer    :: nCorner
   integer    :: i, j, k

   real(adqt) :: wtiso
   real(adqt) :: diagInv
   real(adqt) :: t
   real(adqt) :: v

!  Constants

   nZoneSets = getNumberOfZoneSets(Quad)
   wtiso     = Size% wtiso


   TOMP_MAP(target enter data map(to: wtiso))

   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
   TOMPC(shared(nZoneSets, Geom, GTA, wtiso)&)
   TOMPC(private(cc,c0,nCorner,diagInv,t,v))

   ZoneSetLoop: do zSetID=1,nZoneSets

     !$omp  parallel do default(none)  &
     !$omp& shared(Geom, GTA, wtiso, zSetID)  &
     !$omp& private(cc,c0,nCorner,diagInv,t,v)

     ZoneLoop: do zone=Geom% zone1(zSetID),Geom% zone2(zSetID)

!      Update Scalar Intensity (solve A*Phi = y)

       nCorner = Geom% numCorner(zone)
       c0      = Geom% cOffSet(zone)

       do c=1,nCorner
         GTA% Sscat(c0+c) = zero 
         do cc=1,nCorner
           GTA% Sscat(c0+c)  =  GTA% Sscat(c0+c) + GTA% TT(cc,c0+c)*  &
                               wtiso*GTA%GreySource(c0+cc)
           GTA% TT(cc,c0+c) = -wtiso*GTA%GreySigScat(c0+cc)* &
                               GTA% TT(cc,c0+c)
         enddo
         GTA% TT(c,c0+c)    =  one + GTA% TT(c,c0+c) 
       enddo

!      Decompose:  A = LU 

       do i=1,nCorner

         t = zero

         do k=1,i-1
           t = t + GTA% TT(k,c0+i)*GTA% TT(i,c0+k)
         enddo

         GTA% TT(i,c0+i) = GTA% TT(i,c0+i) - t
         diagInv         = one/GTA% TT(i,c0+i)

         do j=i+1,nCorner

           t = zero
           v = zero

           do k=1,i-1
             t = t + GTA% TT(k,c0+i)*GTA% TT(j,c0+k)
             v = v + GTA% TT(k,c0+j)*GTA% TT(i,c0+k)
           enddo

           GTA% TT(j,c0+i) = GTA% TT(j,c0+i) - t
           GTA% TT(i,c0+j) = diagInv*(GTA% TT(i,c0+j) - v)

         enddo

       enddo

     enddo ZoneLoop

!$omp end parallel do

   enddo ZoneSetLoop

   TOMP(end target teams distribute)

   TOMP_MAP(target exit data map(release: wtiso))

   return
   end subroutine ScalarIntensityDecompose_GPU


   subroutine ScalarIntensitySolve_GPU(P)

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use GreyAcceleration_mod
   use Options_mod
   use QuadratureList_mod

   implicit none

!  Arguments

   real(adqt), intent(inout) :: P(Size% ncornr)

!  Local

   integer    :: zSetID
   integer    :: nZoneSets
   integer    :: zone
   integer    :: c
   integer    :: c0
   integer    :: nCorner
   integer    :: i, k

   real(adqt) :: t

!  Constants

   nZoneSets = getNumberOfZoneSets(Quad)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, P, Geom, GTA)&)
   TOMPC(private(c0,nCorner,t))

   ZoneSetLoop: do zSetID=1,nZoneSets

     !$omp  parallel do default(none)  &
     !$omp& shared(P, Geom, GTA, zSetID)  &
     !$omp& private(c0,nCorner,t)

     ZoneLoop: do zone=Geom% zone1(zSetID),Geom% zone2(zSetID)

!      Update Scalar Intensity (solve A*Phi = y)

       nCorner = Geom% numCorner(zone)
       c0      = Geom% cOffSet(zone)

       do c=1,nCorner
         P(c0+c) = GTA% PhiInc(c0+c)
       enddo

!      Solve Ly = S

       do k=2,nCorner
         t = zero
         do i=1,k-1
           t = t - GTA% TT(i,c0+k)*P(c0+i)
         enddo
         P(c0+k) = P(c0+k) + t
       enddo

!      Solve Ux = y

       P(c0+nCorner) = P(c0+nCorner)/GTA% TT(nCorner,c0+nCorner)

       do k=nCorner-1,1,-1
         t = zero

         do i=k+1,nCorner
           t = t + P(c0+i)*GTA% TT(i,c0+k)
         enddo

         P(c0+k) = (P(c0+k) - t)/GTA% TT(k,c0+k)
       enddo

     enddo ZoneLoop

!$omp end parallel do

   enddo ZoneSetLoop

   TOMP(end target teams distribute)


   return
   end subroutine ScalarIntensitySolve_GPU

