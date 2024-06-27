#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                         Last Update: 10/2016 PFN                     *
!                                                                      *
!   initPhiTotal  -  Called from within Teton to initialize the total  *
!                    scalar radiation intensity (PhiTotal) by          *
!                    integrating Psi over Sets.                        *
!                                                                      *
!***********************************************************************

   subroutine initPhiTotal

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use Size_mod
   use constant_mod
   use Geometry_mod
   use SetData_mod
   use AngleSet_mod
   use RadIntensity_mod
   use QuadratureList_mod
   use OMPWrappers_mod

   implicit none

!  Local

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet

   integer    :: angle
   integer    :: numAngles
   integer    :: c
   integer    :: g
   integer    :: g0
   integer    :: Groups
   integer    :: setID
   integer    :: nSets
   integer    :: zSetID
   integer    :: nZoneSets
   integer    :: ngr

   real(adqt) :: quadwt
   real(adqt) :: volRatio

!  Constants

   nSets     = getNumberOfSets(Quad)
   nZoneSets = getNumberOfZoneSets(Quad)
   ngr       = Size% ngr

   if ( Size% useGPU ) then

     TOMP_MAP(target enter data map(to: nSets, ngr))

#ifdef TETON_ENABLE_OPENACC
     !$acc parallel loop gang num_gangs(nZoneSets) &
     !$acc& vector_length(omp_device_team_thread_limit) &
     !$acc& private(Set, ASet, g0, Groups, NumAngles, quadwt, volRatio, Angle, setID)
#else
     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nZoneSets, Geom, Rad, ngr, Quad, nSets )&)
     TOMPC(private(Set, ASet, g0, Groups, NumAngles, quadwt, volRatio, Angle, setID))
#endif

     ZoneSetLoop1: do zSetID=1,nZoneSets

#ifdef TETON_ENABLE_OPENACC
       !$acc loop vector collapse(2)
#else
       !$omp  parallel do collapse(2) default(none)  &
       !$omp& shared(Geom, Rad, ngr, zSetID)
#endif
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         do g=1,ngr
           Rad% PhiTotal(g,c) = zero 
         enddo
       enddo
#ifndef TETON_ENABLE_OPENACC
       !$omp end parallel do
#endif

!      Sum over phase-space sets

       do setID=1,nSets

         Set       => Quad% SetDataPtr(setID) 
         ASet      => Quad% AngSetPtr(Set% angleSetID)

         NumAngles =  Set% NumAngles
         Groups    =  Set% Groups
         g0        =  Set% g0

         AngleLoop: do Angle=1,NumAngles
           quadwt =  ASet% weight(Angle)

#ifdef TETON_ENABLE_OPENACC
           !$acc loop vector collapse(2) &
           !$acc& private(volRatio)
#else
           !$omp  parallel do collapse(2) default(none)  &
           !$omp& shared(Geom, Rad, Set, g0, Groups, Angle, quadwt, zSetID) &
           !$omp& private(volRatio)
#endif
           do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
             do g=1,Groups
               volRatio              = Geom% VolumeOld(c)/Geom% Volume(c)
               Set% Psi(g,c,Angle)   = volRatio*Set% Psi(g,c,Angle)

               Rad% PhiTotal(g0+g,c) = Rad% PhiTotal(g0+g,c) + &
                                       quadwt*Set% Psi(g,c,Angle)
             enddo
           enddo

#ifndef TETON_ENABLE_OPENACC
          !$omp end parallel do
#endif

         enddo AngleLoop

       enddo

     enddo ZoneSetLoop1

#ifdef TETON_ENABLE_OPENACC
     !$acc end parallel loop
#else
     TOMP(end target teams distribute)
#endif

     TOMP_MAP(target exit data map(release: nSets, ngr))

   else

     Rad% PhiTotal(:,:) = zero

     !$omp parallel do default(none) schedule(static) &
     !$omp& private(Set, ASet, zSetID, setID)  &
     !$omp& private(c, g, g0, Groups, Angle, NumAngles, quadwt, volRatio) &
     !$omp& shared(Quad, Geom, Rad, nZoneSets, nSets)

     ZoneSetLoop2: do zSetID=1,nZoneSets

       do setID=1,nSets

         Set       => getSetData(Quad, setID)
         ASet      => getAngleSetFromSetID(Quad, setID)

         NumAngles =  Set% NumAngles
         Groups    =  Set% Groups
         g0        =  Set% g0

         AngleLoop2: do Angle=1,NumAngles
           quadwt =  ASet% weight(Angle)

           do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
             volRatio = Geom% VolumeOld(c)/Geom% Volume(c)
             do g=1,Groups
               Set% Psi(g,c,Angle)   = volRatio*Set% Psi(g,c,Angle)

               Rad% PhiTotal(g0+g,c) = Rad% PhiTotal(g0+g,c) + &
                                       quadwt*Set% Psi(g,c,Angle)
             enddo
           enddo
         enddo AngleLoop2

       enddo

     enddo ZoneSetLoop2

     !$omp end parallel do

   endif


   return
   end subroutine initPhiTotal 

