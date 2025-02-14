#include "macros.h"
#include "omp_wrappers.h"

!***********************************************************************
!                        Version 1:  04/2008, PFN                      *
!                                                                      *
!   setGTAOpacity - calculates grey opacities for grey-transport       *
!                   acceleration (GTA) in 2D/3D.                       *
!                                                                      *
!***********************************************************************
   subroutine setGTAOpacityNEW_GPU 

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use Material_mod
   use GreyAcceleration_mod
   use QuadratureList_mod
   use ZoneSet_mod

   implicit none

!  Local Variables

   integer    :: zone
   integer    :: g
   integer    :: c
   integer    :: ngr
   integer    :: zSetID
   integer    :: nZoneSets

   real(adqt) :: tau
   real(adqt) :: greysigs
   real(adqt) :: SigtInv
   real(adqt) :: ChiSigt
   real(adqt) :: scatRatio

   real(adqt), parameter :: minRatio = 1.0e-10_adqt

!  Constants

   nZoneSets = getNumberOfZoneSets(Quad)
   tau       = Size% tau
   ngr       = Size% ngr

!  If the CUDA solver is used we need to map Eta/Chi

   if ( Size%useCUDASolver ) then
     TOMP_UPDATE(target update to(GTA% Chi))
     TOMP_UPDATE(target update to(Mat% Eta))
   endif

   TOMP_MAP(target enter data map(to: tau, ngr))


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, ZSet, Geom, Mat, tau, ngr))

   do zSetID=1,nZoneSets

     !$omp parallel do collapse(2) default(none) schedule(dynamic)  &
     !$omp& shared(ZSet, Geom, Mat, tau, ngr, zSetID)

     do zone=Geom% zone1(zSetID),Geom% zone2(zSetID)
       do g=1,ngr
         ZSet% B(g,zone) = one/(Mat% Siga(g,zone) + Mat% Sigs(g,zone) + tau)
       enddo
     enddo

     !$omp end parallel do

   enddo

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, ZSet, Geom))

   do zSetID=1,nZoneSets

     !$omp parallel do default(none) schedule(dynamic)  &
     !$omp& shared(zSetID, ZSet, Geom) 

     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       ZSet% sumT(c)    = zero 
       ZSet% delta(c)   = zero 
       ZSet% netRate(c) = zero 
     enddo

     !$omp end parallel do

   enddo

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, ZSet, Geom, Mat, GTA, ngr) &)
   TOMPC(private(zone, SigtInv, ChiSigt) )

   do zSetID=1,nZoneSets

     !$omp parallel do default(none) schedule(dynamic)  &
     !$omp& shared(ZSet, Geom, Mat, GTA, ngr, zSetID)  &
     !$omp& private(zone, SigtInv, ChiSigt) 

     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       do g=1,ngr
         zone             = Geom% CToZone(c)
         SigtInv          = ZSet% B(g,zone)
         ChiSigt          = GTA% Chi(g,c)*SigtInv

         ZSet% sumT(c)    = ZSet% sumT(c)    + ChiSigt
         ZSet% delta(c)   = ZSet% delta(c)   + ChiSigt*SigtInv
         ZSet% netRate(c) = ZSet% netRate(c) + ChiSigt*Mat% Siga(g,zone)
         GTA% Chi(g,c)    = ChiSigt 
       enddo
     enddo

     !$omp end parallel do

   enddo

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, ZSet, GTA, Geom, ngr))

   do zSetID=1,nZoneSets

     !$omp parallel do default(none) schedule(dynamic)  &
     !$omp& shared(zSetID, ZSet, GTA, Geom, ngr)

     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       if ( ZSet% sumT(c) > zero ) then
         do g=1,ngr
           GTA% Chi(g,c) = GTA% Chi(g,c)/ZSet% sumT(c)
         enddo
       endif
     enddo

     !$omp end parallel do

   enddo

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
   TOMPC(shared(nZoneSets, ZSet, Geom, Mat, tau))

   do zSetID=1,nZoneSets

     !$omp  parallel do default(none) schedule(dynamic)  &
     !$omp& shared(zSetID, ZSet, Geom, Mat, tau)

     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       if ( ZSet% sumT(c) > zero ) then
         ZSet% dTCompton(c) = ZSet% sumT(c)/ZSet% delta(c)
         ZSet% comptonSe(c) = tau + (one - Mat% Eta(c))*ZSet% netRate(c)/ZSet% sumT(c)
       else
         ZSet% dTCompton(c) = tau
         ZSet% comptonSe(c) = tau
       endif
     enddo

     !$omp end parallel do

   enddo

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
   TOMPC(shared(nZoneSets, ZSet, GTA, Geom)&)
   TOMPC(private(greysigs, scatRatio))

   do zSetID=1,nZoneSets

     !$omp  parallel do default(none) schedule(dynamic)  &
     !$omp& shared(zSetID, ZSet, GTA, Geom) private(greysigs, scatRatio)

     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       greysigs  = ZSet% dTCompton(c) - ZSet% comptonSe(c)
       scatRatio = greysigs/ZSet% dTCompton(c)

       if (scatRatio <= minRatio) then
         GTA%GreySigScat(c)  = zero
         GTA%GreySigTotal(c) = ZSet% comptonSe(c)
       else
         GTA%GreySigScat(c)  = greysigs
         GTA%GreySigTotal(c) = ZSet% dTCompton(c)
       endif
     enddo

     !$omp end parallel do

   enddo

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
   TOMPC(shared(nZoneSets, GTA, Geom))

   do zSetID=1,nZoneSets

!$omp parallel do default(none) schedule(dynamic)  &
!$omp& shared(zSetID, GTA, Geom)

     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       GTA%GreySigScatVol(c) = GTA%GreySigScat(c)*Geom% Volume(c)
       GTA%GreySigtInv(c)    = one/GTA%GreySigTotal(c)
     enddo

     !$omp end parallel do

   enddo

   TOMP(end target teams distribute)

   TOMP_MAP(target exit data map(release: tau, ngr))

   TOMP_UPDATE(target update from(GTA% GreySigScatVol))
   TOMP_UPDATE(target update from(GTA% Chi))

 
   return
   end subroutine setGTAOpacityNEW_GPU 

