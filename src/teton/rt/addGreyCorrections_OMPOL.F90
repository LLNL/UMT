#include "macros.h"
#include "omp_wrappers.h"

!***********************************************************************
!                       Last Update:  10/2016, PFN                     *
!                                                                      *
!   addGreyCorrections - Update scalar and angle-dependent intensities *
!                        withgrey corrections.                         *
!                                                                      *
!***********************************************************************
   subroutine addGreyCorrections_GPU 

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use GreyAcceleration_mod
   use Communicator_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use Boundary_mod
   use SetData_mod
   use AngleSet_mod
   use CommSet_mod
   use ZoneSet_mod

   implicit none

!  Local

   type(SetData),          pointer  :: Set
   type(AngleSet),         pointer  :: ASet
   type(CommSet),          pointer  :: CSet
   type(Communicator),     pointer  :: CommT
   type(Boundary),         pointer  :: BdyT
   type(HypPlane),         pointer  :: HypPlanePtr

   integer    :: zone
   integer    :: nzones
   integer    :: b
   integer    :: b0
   integer    :: c
   integer    :: c0
   integer    :: nCorner
   integer    :: angle
   integer    :: exit_angle
   integer    :: NumAngles
   integer    :: g
   integer    :: g0
   integer    :: ngr
   integer    :: Groups
   integer    :: i
   integer    :: setID
   integer    :: nSets
   integer    :: zSetID
   integer    :: nZoneSets
   integer    :: sharedID
   integer    :: nShared
   integer    :: reflID
   integer    :: nReflecting
   integer    :: nBdyElem

   real(adqt) :: wtiso

!  Constants

   nzones      = Size% nzones 
   wtiso       = Size% wtiso
   ngr         = Size% ngr
   nShared     = getNumberOfShared(RadBoundary)
   nReflecting = getNumberOfReflecting(RadBoundary)
   nSets       = getNumberOfSets(Quad)
   nZoneSets   = getNumberOfZoneSets(Quad)

!  Compute the group-dependent corrections

   TOMP_UPDATE(target update to(GTA% GreyCorrection))

   TOMP_MAP(target enter data map(to: ngr, wtiso))

#ifdef TETON_ENABLE_OPENACC
   !$acc parallel loop gang num_gangs(nZoneSets) vector_length(omp_device_team_thread_limit)
#else
   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, GTA, Geom, Rad, ngr))
#endif
   do zSetID=1,nZoneSets

#ifdef TETON_ENABLE_OPENACC
     !$acc loop vector collapse(2)
#else
     !$omp parallel do collapse(2) default(none) schedule(dynamic)  &
     !$omp& shared(zSetID, GTA, Geom, Rad, ngr)
#endif
     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       do g=1,ngr
         Rad% PhiTotal(g,c) = Rad% PhiTotal(g,c) + GTA%GreyCorrection(c)*GTA% Chi(g,c)
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


#ifdef TETON_ENABLE_OPENACC
   !$acc parallel loop gang num_gangs(nZoneSets) vector_length(omp_device_team_thread_limit)
#else
   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, ZSet, Geom, Rad))
#endif
   do zSetID=1,nZoneSets

#ifdef TETON_ENABLE_OPENACC
     !$acc loop vector
#else
     !$omp parallel do default(none) schedule(dynamic)  &
     !$omp& shared(zSetID, ZSet, Geom, Rad)
#endif
     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       ZSet% sumT(c) = sum( Rad% PhiTotal(:,c) )
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


#ifdef TETON_ENABLE_OPENACC
   !$acc  parallel loop gang num_gangs(nZoneSets) vector_length(omp_device_team_thread_limit) &
   !$acc& private(c0, nCorner)
#else
   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, Geom, Rad, ZSet)&)
   TOMPC(private(c0, nCorner))
#endif
   do zSetID=1,nZoneSets

#ifdef TETON_ENABLE_OPENACC
     !$acc loop vector &
     !$acc& private(c0, nCorner)
#else
     !$omp parallel do default(none) schedule(dynamic)  &
     !$omp& shared(zSetID, Geom, Rad, ZSet) private(c0, nCorner) 
#endif
     do zone=Geom% zone1(zSetID),Geom% zone2(zSetID)
       nCorner              = Geom% numCorner(zone)
       c0                   = Geom% cOffSet(zone)
       Rad% radEnergy(zone) = zero

       do c=1,nCorner
         Rad% radEnergy(zone) = Rad% radEnergy(zone) + &
                                Geom% Volume(c0+c)*ZSet% sumT(c0+c)
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

#ifdef TETON_ENABLE_OPENACC
   !$acc  parallel loop gang num_gangs(nSets) vector_length(omp_device_team_thread_limit) &
   !$acc& private(Set, ASet, HypPlanePtr, Groups, NumAngles, c, g0)
#else
   TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nSets, Quad, GTA, wtiso)&)
   TOMPC(private(Set, ASet, HypPlanePtr, Groups, NumAngles, c, g0))
#endif
   do setID=1,nSets
     Set        => Quad% SetDataPtr(setID)
     ASet       => Quad% AngSetPtr(Set% angleSetID)
     Groups     =  Set% Groups
     g0         =  Set% g0
     NumAngles  =  Set% NumAngles

     do angle=1,NumAngles
       HypPlanePtr => ASet% HypPlanePtr(angle)

#ifdef TETON_ENABLE_OPENACC
       !$acc  loop vector collapse(2) &
       !$acc& private(c)
#else
       !$omp  parallel do collapse(2) default(none) &
       !$omp& shared(Set, HypPlanePtr, GTA, angle, g0, Groups, wtiso) &
       !$omp& private(c)
#endif
       do i=1,HypPlanePtr% interfaceLen
         do g=1,Groups
           c = HypPlanePtr% interfaceList(i)
           Set% PsiInt(g,i,angle) = Set% PsiInt(g,i,angle) + wtiso*  &
                                    GTA%GreyCorrection(c)*GTA% Chi(g0+g,c) 
         enddo
       enddo
#ifndef TETON_ENABLE_OPENACC
       !$omp end parallel do
#endif

     enddo

   enddo

#ifdef TETON_ENABLE_OPENACC
   !$acc end parallel loop
#else
   TOMP(end target teams distribute)
#endif
 

TOMP_MAP(target exit data map(release: ngr, wtiso))

TOMP_UPDATE(target update from(Rad% radEnergy))

!  Update Set dependent boundary fluxes

   !$omp  parallel do default(none) schedule(static) &
   !$omp& private(NumAngles, g0, Groups, b, c, reflID, nBdyElem, exit_angle) &
   !$omp& private(b0, Set, ASet, CSet, CommT, BdyT) &
   !$omp& shared(nSets, nShared, nReflecting, Quad, GTA, RadBoundary, wtiso)

   SetLoop: do setID=1,nSets

     Set       => getSetData(Quad, setID)
     ASet      => getAngleSetFromSetID(Quad, setID)
     CSet      => getCommSetFromSetID(Quad, setID)

     NumAngles =  Set% NumAngles
     g0        =  Set% g0
     Groups    =  Set% Groups 

!    In spatially decomposed runs, also correct the boundary flux

     SharedLoop: do sharedID=1,nShared

       do angle=1,NumAngles
         CommT => getMessage(CSet, sharedID, angle)
         do i=1,CommT% nSend
           b =  CommT% ListSend(1,i)
           c =  CommT% ListSend(2,i)

           do g=1,Groups
             Set% PsiB(g,b,angle) = Set% PsiB(g,b,angle) + wtiso*  & 
                                    GTA%GreyCorrection(c)*GTA% Chi(g0+g,c)
           enddo
         enddo
       enddo

     enddo SharedLoop

!    Opposing Reflecting Boundaries

     ReflLoop: do reflID=1,nReflecting

       BdyT      => getReflecting(RadBoundary, reflID)
       nBdyElem  =  getNumberOfBdyElements(BdyT)
       b0        =  getFirstBdyElement(BdyT) - 1

       do i=1,ASet% nExit(reflID)
         exit_angle = ASet% ExitAngleList(i,reflID)
         do b=1,nBdyElem
           c = BdyT% BdyToC(b)

           do g=1,Groups
             Set% PsiB(g,b0+b,exit_angle) = Set% PsiB(g,b0+b,exit_angle) + wtiso*  &
                                            GTA%GreyCorrection(c)*GTA% Chi(g0+g,c) 
           enddo
         enddo
       enddo

     enddo ReflLoop

   enddo SetLoop
   !$omp end parallel do

 
   return
   end subroutine addGreyCorrections_GPU 

