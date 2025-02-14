#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                         Last Update: 10/2016 PFN                     *
!                                                                      *
!    getPhiTotal  -  Called from within Teton to retrieve the total    *
!                    scalar radiation intensity (PhiTotal) by summing  *
!                    over Sets.                                        *
!                                                                      *
!***********************************************************************

   subroutine getPhiTotal(sendIndex)

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use Size_mod
   use constant_mod
   use Geometry_mod
   use RadIntensity_mod
   use SetData_mod
   use AngleSet_mod
   use QuadratureList_mod

   implicit none

!  Arguments

   integer, optional, intent(in) :: sendIndex

!  Local

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet

   integer    :: Groups
   integer    :: g
   integer    :: c
   integer    :: g0
   integer    :: setID
   integer    :: nSets
   integer    :: zSetID
   integer    :: nZoneSets
   integer    :: ngr
   integer    :: Angle
   real(adqt) :: quadwt

!  Constants
   nSets     = getNumberOfSets(Quad)
   nZoneSets = getNumberOfZoneSets(Quad)
   ngr       = Size% ngr

   if (Size%useGPU) then

     TOMP_MAP(target enter data map(to: nSets, ngr, sendIndex))

     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
     TOMPC(shared(Geom, Rad, ngr, nZoneSets, sendIndex, nSets, Quad)&)
     TOMPC(private(setID, Set, ASet, g0, Groups, Angle, quadwt))

     ZoneSetLoop1: do zSetID=1,nZoneSets

       if (sendIndex == 1) then

         !$omp  parallel do collapse(2) default(none)  &
         !$omp& shared(Geom, Rad, ngr, zSetID)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           do g=1,ngr
             Rad% PhiTotal(g,c) = zero 
           enddo
         enddo

         !$omp end parallel do

       endif

!      Sum over phase-space sets

       do setID=1,nSets

         Set    => Quad% SetDataPtr(setID) 
         ASet   => Quad% AngSetPtr(Set% angleSetID)

         Angle  =  Set% AngleOrder(sendIndex)
         Groups =  Set% Groups
         g0     =  Set% g0
         quadwt =  ASet% weight(Angle)

         if (.not. ASet% StartingDirection(Angle) ) then

         !$omp  parallel do collapse(2) default(none)  &
         !$omp& shared(Geom, Rad, Set, g0, Groups, quadwt, zSetID)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           do g=1,Groups
             Rad% PhiTotal(g0+g,c) = Rad% PhiTotal(g0+g,c) + quadwt*Set% Psi1(g,c)
           enddo
         enddo

         !$omp end parallel do

         endif

       enddo

     enddo ZoneSetLoop1

     TOMP(end target teams distribute)

     TOMP_MAP(target exit data map(release: nSets, ngr, sendIndex))

   else

     Rad% PhiTotal(:,:) = zero

     !$omp parallel do default(none) schedule(static) &
     !$omp& shared(Quad, Geom, Rad, nZoneSets, nSets) &
     !$omp private(Set, zSetID, g0, Groups)

     ZoneSetLoop2: do zSetID=1,nZoneSets

       do setID=1,nSets

         Set    => getSetData(Quad, setID)
         Groups =  Set% Groups
         g0     =  Set% g0

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           do g=1,Groups
             Rad% PhiTotal(g0+g,c) = Rad% PhiTotal(g0+g,c) + Set% Phi(g,c)
           enddo
         enddo

       enddo

     enddo ZoneSetLoop2

     !$omp end parallel do

   endif


   return
   end subroutine getPhiTotal 
