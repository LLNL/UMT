#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  03/2013, PFN                    *
!                                                                      *
!   getCollisionRate - Computes the total collision rate. This         *
!                      quantity is used to compute the residual        *
!                      source for grey-transport acceleration (GTA).   *
!                                                                      *
!***********************************************************************
   subroutine getCollisionRate_GPU(residualFlag) 

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use Material_mod
   use GreyAcceleration_mod
   use QuadratureList_mod
   use ZoneSet_mod

   implicit none

!  Arguments

   integer,  intent(in) :: residualFlag

!  Local

   integer    :: c
   integer    :: g
   integer    :: ngr
   integer    :: zone 
   integer    :: zSetID
   integer    :: nZoneSets

!  Constants

   nZoneSets = getNumberOfZoneSets(Quad)
   ngr       = Size% ngr

!  Calculate the total energy absorption rate density 


TOMP_MAP(target enter data map(to: ngr))

#ifdef TETON_ENABLE_OPENACC
   !$acc parallel loop gang num_gangs(nZoneSets) &
   !$acc& vector_length(omp_device_team_thread_limit) &
   !$acc& private(zone)
#else
   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
   TOMPC(shared(nZoneSets, ZSet, Geom, Rad, Mat, ngr) &)
   TOMPC(private(zone))
#endif

   do zSetID=1,nZoneSets

#ifdef TETON_ENABLE_OPENACC
     !$acc loop vector collapse(2) &
     !$acc& private(zone)
#else
     !$omp  parallel do default(none) collapse(2) schedule(dynamic) &
     !$omp& shared(nZoneSets, ZSet, zSetID, Geom, Rad, Mat, ngr)   &
     !$omp& private(zone) 
#endif
     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       do g=1,ngr
         zone          =  Geom% CToZone(c)
         ZSet% ex(g,c) = (Mat%Eta(c)*Mat%siga(g,zone) + Mat%sigs(g,zone))* &
                          Rad% PhiTotal(g,c)
       enddo
     enddo
#ifndef TETON_ENABLE_OPENACC
     !$omp end parallel do
#endif

   enddo

#ifdef TETON_ENABLE_OPENACC
   !$acc end parallel loop
#else
   TOMP(end target teams distribute)
#endif


TOMP_MAP(target exit data map(release: ngr))

   if (residualFlag == 0) then

#ifdef TETON_ENABLE_OPENACC
     !$acc parallel loop gang num_gangs(nZoneSets) &
     !$acc& vector_length(omp_device_team_thread_limit)
#else
     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
     TOMPC(shared(nZoneSets, ZSet, Geom, GTA))
#endif

     do zSetID=1,nZoneSets

#ifdef TETON_ENABLE_OPENACC
       !$acc loop vector
#else
       !$omp  parallel do default(none) schedule(dynamic) &
       !$omp& shared(nZoneSets, ZSet, zSetID, Geom, GTA)
#endif
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         GTA% GreySource(c) = sum( ZSet% ex(:,c) ) 
       enddo
#ifndef TETON_ENABLE_OPENACC
       !$omp end parallel do
#endif

     enddo

#ifdef TETON_ENABLE_OPENACC
     !$acc end parallel loop
#else
     TOMP(end target teams distribute)
#endif

   else

#ifdef TETON_ENABLE_OPENACC
     !$acc parallel loop gang num_gangs(nZoneSets) &
     !$acc& vector_length(omp_device_team_thread_limit)
#else
     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
     TOMPC(shared(nZoneSets, ZSet, Geom, GTA))
#endif

     do zSetID=1,nZoneSets

#ifdef TETON_ENABLE_OPENACC
       !$acc loop vector
#else
       !$omp  parallel do default(none) schedule(dynamic) &
       !$omp& shared(nZoneSets, ZSet, zSetID, Geom, GTA)
#endif
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         GTA% GreySource(c) = sum( ZSet% ex(:,c) ) - GTA% GreySource(c) 
       enddo
#ifndef TETON_ENABLE_OPENACC
       !$omp end parallel do
#endif

     enddo

#ifdef TETON_ENABLE_OPENACC
     !$acc end parallel loop
#else
     TOMP(end target teams distribute)
#endif

   endif

 
   return
   end subroutine getCollisionRate_GPU 

