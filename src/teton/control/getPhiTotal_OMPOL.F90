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
   use, intrinsic :: iso_c_binding, only : c_int
   use Options_mod
   use kind_mod
   use Size_mod
   use constant_mod
   use Geometry_mod
   use SetData_mod
   use AngleSet_mod
   use QuadratureList_mod

   implicit none

!  Arguments

   integer, optional, intent(in) :: sendIndex

!  Local

   type(SetData),    pointer :: Set

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
   type(AngleSet),   pointer :: ASet
   real(adqt) :: quadwt
   integer(kind=c_int) :: nOmpMaxTeamThreads

!  Constants

   nOmpMaxTeamThreads = Options%getNumOmpMaxTeamThreads()
   nSets     = getNumberOfSets(Quad)
   nZoneSets = getNumberOfZoneSets(Quad)
   ngr       = Size% ngr

   if (Size%useGPU) then

     TOMP(target data map(to: nSets, ngr, sendIndex))

     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(nOmpMaxTeamThreads) &)
     TOMPC(private(Set, ASet, zSetID, setID, g0, Groups, Angle, quadwt))

     ZoneSetLoop1: do zSetID=1,nZoneSets

       if (sendIndex == 1) then

         !$omp  parallel do collapse(2) default(none)  &
         !$omp& shared(Geom, ngr, zSetID) private(c,g)
         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           do g=1,ngr
             Geom% PhiTotal(g,c) = zero 
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
         !$omp& shared(Geom, Set, g0, Groups, quadwt, zSetID) private(c,g)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           do g=1,Groups
             Geom% PhiTotal(g0+g,c) = Geom% PhiTotal(g0+g,c) + quadwt*Set% Psi1(g,c)
           enddo
         enddo

         !$omp end parallel do

         endif

       enddo

     enddo ZoneSetLoop1


     TOMP(end target teams distribute)
     TOMP(end target data)

   else

     Geom% PhiTotal(:,:) = zero

     !$omp parallel do private(Set, zSetID, setID, c, g, g0, Groups) &
     !$omp& shared(Quad, Geom) schedule(static)

     ZoneSetLoop2: do zSetID=1,nZoneSets

       do setID=1,nSets

         Set    => getSetData(Quad, setID)
         Groups =  Set% Groups
         g0     =  Set% g0

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           do g=1,Groups
             Geom% PhiTotal(g0+g,c) = Geom% PhiTotal(g0+g,c) + Set% Phi(g,c)
           enddo
         enddo

       enddo

     enddo ZoneSetLoop2

     !$omp end parallel do

   endif


   return
   end subroutine getPhiTotal 
