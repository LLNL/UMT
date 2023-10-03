#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  03/2023, PFN                    *
!                                                                      *
!  CornerSweepUCBxyz_GPU:                                              *
!                                                                      * 
!                    This routine calculates angular fluxes for a      *
!                    single direction and multiple energy groups for   *
!                    for an upstream corner-balance (UCB) spatial      *
!                    in xyz-geometry.                                  *
!                                                                      *
!                    The work is offloaded to a GPU in which each      *
!                    computational "set" is a GPU block. The threads   *
!                    assigned to the block compute one group in one    *
!                    corner in a hyperplane (by definition, all of     *
!                    the corners in a hyperplane are independent).     *
!                                                                      *
!***********************************************************************

   subroutine CornerSweepUCBxyz_GPU(nSets, sendIndex, savePsi)

   use, intrinsic :: iso_c_binding, only : c_int
   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use GroupSet_mod
   use ArrayChecks_mod

   implicit none

!  Arguments

   integer,          intent(in)    :: nSets
   integer,          intent(in)    :: sendIndex
   logical (kind=1), intent(in)    :: savePsi

!  Local

   type(SetData),    pointer       :: Set
   type(AngleSet),   pointer       :: ASet
   type(GroupSet),   pointer       :: GSet
   type(HypPlane),   pointer       :: HypPlanePtr
   type(BdyExit),    pointer       :: BdyExitPtr

   integer            :: setID
   integer            :: zSetID
   integer            :: Angle
   integer            :: g
   integer            :: Groups

   integer            :: mCycle
   integer            :: offSet
   integer            :: nAngleSets
   integer            :: nZoneSets

   integer            :: nzones
   integer            :: ii
   integer            :: ndoneZ
   integer            :: hyperPlane
   integer            :: nHyperplanes

   real(adqt)         :: tau

!  Local

   integer    :: b 
   integer    :: i 
   integer    :: cfp 
   integer    :: cface 
   integer    :: ifp
   integer    :: c
   integer    :: c0 
   integer    :: cez
   integer    :: nCFaces

   integer    :: zone
   integer    :: zone0
   integer    :: nCorner

   real(adqt), parameter :: fouralpha=1.82d0

   real(adqt) :: aez
   real(adqt) :: area_opp

   real(adqt) :: source
   real(adqt) :: sig
   real(adqt) :: vol
   real(adqt) :: sigv
   real(adqt) :: sigv2
   real(adqt) :: gnum
   real(adqt) :: gden
   real(adqt) :: sez
   real(adqt) :: psi_opp
   real(adqt) :: afp

!  Dynamic

   integer, allocatable :: angleList(:)

!  Constants

   tau        = Size% tau
   nAngleSets = getNumberOfAngleSets(Quad)
   nZoneSets  = getNumberOfZoneSets(Quad)

   allocate( angleList(nAngleSets) )

   do setID=1,nSets
     Set => Quad% SetDataPtr(setID)
     angleList(Set% angleSetID) = Set% AngleOrder(sendIndex)
   enddo

!  Here the maximum block size is the product of the maximum
!  number of zones in a hyperplane and the number of groups;
!  The maximum value is the over all teams  

!  Note: num_blocks = nSets and the number of threads per 
!  team (a.k.a. "block") <= block_threads

!  Initialize

   ! Verify we won't get out-of-bounds accesses below.
   TETON_CHECK_BOUNDS1(Quad%AngSetPtr, nAngleSets)
   TETON_CHECK_BOUNDS1(Geom%corner1, nZoneSets)
   TETON_CHECK_BOUNDS1(Geom%corner2, nZoneSets)

#ifdef TETON_ENABLE_OPENACC
   !$acc data copyin(tau, sendIndex, angleList)
#else
   TOMP(target enter data map(to: tau, sendIndex, angleList))
#endif

#ifdef TETON_ENABLE_OPENACC
   !$acc parallel loop gang num_gangs(nZoneSets) &
   !$acc& vector_length(omp_device_team_thread_limit) &
   !$acc& private(ASet, setID, Angle)
#else
   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(private(ASet, setID, Angle) &)
   TOMPC(shared(nZoneSets, angleList, Quad, Geom, nAngleSets) )
#endif

   ZoneSetLoop: do zSetID=1,nZoneSets

!    Loop over angle sets

     do setID=1,nAngleSets

       ASet  => Quad% AngSetPtr(setID)
       angle =  angleList(setID)

! NOTE: This loop doesn't support a collapse(2), as its not a canonical form
! loop ( the inner loop bounds can not be predetermined ); it's significantly
! faster to split into two loops as below

#ifdef TETON_ENABLE_OPENACC
!$acc loop vector collapse(2)
#else
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Geom, ASet, Angle, zSetID)
#endif
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         do cface=1,3
           ASet% AfpNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_fp(:,cface,c) )
           ASet% AezNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_ez(:,cface,c) )
         enddo
       enddo
#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

#ifdef TETON_ENABLE_OPENACC
!$acc loop vector
#else
!$omp  parallel do default(none)  &
!$omp& shared(Geom, ASet, Angle, zSetID)
#endif
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         do cface=4,Geom% nCFacesArray(c)
           ASet% AfpNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_fp(:,cface,c) )
           ASet% AezNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_ez(:,cface,c) )
         enddo
       enddo
#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

     enddo

   enddo ZoneSetLoop

#ifdef TETON_ENABLE_OPENACC
!$acc end parallel loop
#else
TOMP(end target teams distribute)
#endif


#ifdef TETON_ENABLE_OPENACC
!$acc parallel loop gang num_gangs(nZoneSets) &
!$acc& vector_length(omp_device_team_thread_limit) &
!$acc& private(ASet)
#else
TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
TOMPC(private(ASet) &)
TOMPC(shared(nZoneSets, nAngleSets, Quad, Geom))
#endif

   ZoneSetLoop2: do zSetID=1,nZoneSets

!    Loop over angle sets

     do setID=1,nAngleSets

       ASet  => Quad% AngSetPtr(setID)

#ifdef TETON_ENABLE_OPENACC
!$acc loop vector
#else
!$omp  parallel do default(none)  &
!$omp& shared(Geom, ASet, zSetID)
#endif

       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         ASet% ANormSum(c) = zero
         do cface=1,Geom% nCFacesArray(c)
           ASet% ANormSum(c) = ASet% ANormSum(c) + half*   &
                              (ASet% AfpNorm(cface,c) + abs( ASet% AfpNorm(cface,c) ) +  & 
                               ASet% AezNorm(cface,c) + abs( ASet% AezNorm(cface,c) ) )
         enddo
       enddo

#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

     enddo

   enddo ZoneSetLoop2

#ifdef TETON_ENABLE_OPENACC
!$acc end parallel loop
#else
TOMP(end target teams distribute)
#endif


#ifdef TETON_ENABLE_OPENACC
!$acc parallel loop gang num_gangs(nSets) &
!$acc& vector_length(omp_device_team_thread_limit) &
!$acc& private(Set, ASet, Angle, Groups, offSet)
#else
TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
TOMPC(shared(sendIndex, Quad, nSets) &)
TOMPC(private(Set, ASet, Angle, Groups, offSet))
#endif


   SetLoop0: do setID=1,nSets

     Set          => Quad% SetDataPtr(setID)
     ASet         => Quad% AngSetPtr(Set% angleSetID)

     Groups       =  Set% Groups
     Angle        =  Set% AngleOrder(sendIndex)
     offSet       =  ASet% cycleOffSet(angle)

!  Initialize boundary values in Psi1 and interior values on the cycle
!  list

#ifdef TETON_ENABLE_OPENACC
!$acc loop vector collapse(2) &
!$acc& private(c)
#else
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Angle, Set, ASet, offSet, Groups) private(c)
#endif
     do mCycle=1,ASet% numCycles(Angle)
       do g=1,Groups
         c              = ASet% cycleList(offSet+mCycle)
         Set% Psi1(g,c) = Set% cyclePsi(g,offSet+mCycle)
       enddo
     enddo
#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif


#ifdef TETON_ENABLE_OPENACC
!$acc loop vector collapse(2)
#else
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, Groups, Angle)
#endif
     do c=1,Set%nbelem
       do g=1,Groups
         Set% Psi1(g,Set%nCorner+c) = Set% PsiB(g,c,Angle)
       enddo
     enddo
#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

   enddo SetLoop0

#ifdef TETON_ENABLE_OPENACC
!$acc end parallel loop
#else
TOMP(end target teams distribute)
#endif


! TODO:
! IBM XLF segfaults if 'hyperPlane' is not scoped to private below.
! This should not be necessary, as this is a loop control variables which the runtime should automatically scope to
! private.
!
! Relevant portions of OpenMP spec:
! `The loop iteration variable in any associated loop of a for, parallel for,
! taskloop, or distribute construct is private.`
!
! `A loop iteration variable for a sequential loop in a parallel or task
! generating construct is private in the innermost such construct that encloses
! the loop.`
! 
! Look into reporting this bug to IBM, using UMT as a reproducer.

#ifdef TETON_ENABLE_OPENACC
!$acc parallel loop gang num_gangs(nSets) &
!$acc& vector_length(omp_device_team_thread_limit) &
!$acc& private(Set, ASet, GSet, HypPlanePtr, Angle, Groups) &
!$acc& private(nHyperPlanes, ndoneZ, nzones, hyperPlane)
#else
TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) &)
TOMPC(private(Set, ASet, GSet, HypPlanePtr, Angle, Groups) &)
TOMPC(private(nHyperPlanes, ndoneZ, nzones, hyperPlane)) 
#endif

   SetLoop: do setID=1,nSets

     Set          => Quad% SetDataPtr(setID)
     ASet         => Quad% AngSetPtr(Set% angleSetID) 
     GSet         => Quad% GrpSetPtr(Set% groupSetID) 

     Groups       =  Set% Groups
     Angle        =  Set% AngleOrder(sendIndex)
     nHyperPlanes =  ASet% nHyperPlanes(Angle)
     ndoneZ       =  0
     HypPlanePtr  => ASet% HypPlanePtr(Angle)

     HyperPlaneLoop: do hyperPlane=1,nHyperPlanes

       nzones = HypPlanePtr% zonesInPlane(hyperPlane)

! NOTE: Two loop sections that follow do not support a collapse(3), 
! as its not in canonical form (all loop bounds are not predetermined); 
! it's significantly faster to split into two loops as below

#ifdef TETON_ENABLE_OPENACC
!$acc  loop vector collapse(3) &
!$acc& private(c0,zone,zone0,nCorner,source,nCFaces,afp,cfp) 
#else
!$omp  parallel do collapse(3) default(none) &
!$omp& shared(Set, Geom, ASet, GSet, Angle, nzones, Groups) &
!$omp& shared(ndoneZ, tau) &
!$omp& private(c0,zone,zone0,source,nCFaces,afp,cfp)
#endif

       ZoneLoop0: do ii=1,nzones
         CornerLoop0: do c=1,8
           GroupLoop0: do g=1,Groups

!          Loop through the zones using the NEXTZ list

             zone0   = ASet% nextZ(ndoneZ+ii,Angle)
             zone    = iabs( zone0 )
             c0      = Geom% cOffSet(zone)

!            Contributions from volume terms 

             source         = GSet% STotal(g,c0+c) + tau*Set% Psi(g,c0+c,Angle)
             Set% Q(g,c,ii) = source
             Set% S(g,c,ii) = Geom% Volume(c0+c)*source

             nCFaces = Geom% nCFacesArray(c0+c)

             do cface=1,nCFaces

               afp = ASet% AfpNorm(cface,c0+c)
               cfp = Geom% cFP(cface,c0+c)

               if ( afp < zero ) then
                 Set% S(g,c,ii) = Set% S(g,c,ii) - afp*Set% Psi1(g,cfp)
               endif
             enddo

           enddo GroupLoop0
         enddo CornerLoop0
       enddo ZoneLoop0                                                                                                         

#ifndef TETON_ENABLE_OPENACC                                                                                 
!$omp end parallel do                                                                                                  
#endif

#ifdef TETON_ENABLE_OPENACC
!$acc  loop vector collapse(2) &
!$acc& private(c0,zone,zone0,nCorner,source,nCFaces,afp,cfp) 
#else
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, Geom, ASet, GSet, Angle, nzones, Groups) &
!$omp& shared(ndoneZ, tau) &
!$omp& private(c0,zone,zone0,source,nCFaces,afp,cfp)
#endif

       ZoneLoop00: do ii=1,nzones
         GroupLoop00: do g=1,Groups

!        Loop through the zones using the NEXTZ list

           zone0   = ASet% nextZ(ndoneZ+ii,Angle)
           zone    = iabs( zone0 )

           CornerLoop00: do c=9,Geom% numCorner(zone)

             c0 = Geom% cOffSet(zone)

!            Contributions from volume terms 

             source         = GSet% STotal(g,c0+c) + tau*Set% Psi(g,c0+c,Angle)
             Set% Q(g,c,ii) = source
             Set% S(g,c,ii) = Geom% Volume(c0+c)*source

             nCFaces = Geom% nCFacesArray(c0+c)

             do cface=1,nCFaces

               afp = ASet% AfpNorm(cface,c0+c)
               cfp = Geom% cFP(cface,c0+c)

               if ( afp < zero ) then
                 Set% S(g,c,ii) = Set% S(g,c,ii) - afp*Set% Psi1(g,cfp)
               endif
             enddo

           enddo CornerLoop00
         enddo GroupLoop00
       enddo ZoneLoop00                                                                                                   

#ifndef TETON_ENABLE_OPENACC                                                                                 
!$omp end parallel do                                                                                                  
#endif


#ifdef TETON_ENABLE_OPENACC
!$acc  loop vector collapse(3) &
!$acc& private(c0,cfp,ifp,cez,zone,zone0,nCorner,nCFaces) &
!$acc& private(aez,area_opp,sig,vol) &
!$acc& private(sigv,sigv2,sez,gnum,gden,psi_opp,afp)
#else
!$omp  parallel do collapse(3) default(none) &
!$omp& shared(Set, Geom, ASet, GSet, Angle, nzones, Groups, ndoneZ) &
!$omp& private(c0,cfp,ifp,cez,zone,zone0,nCFaces) &
!$omp& private(aez,area_opp,sig,vol) &
!$omp& private(sigv,sigv2,sez,gnum,gden,psi_opp,afp)
#endif

       ZoneLoop1: do ii=1,nzones
         CornerLoop1: do c=1,8
           GroupLoop1: do g=1,Groups

!            Loop through the zones using the NEXTZ list

             zone0   = ASet% nextZ(ndoneZ+ii,Angle)
             zone    = iabs( zone0 )
             c0      = Geom% cOffSet(zone)

             sig     = GSet% Sigt(g,zone)

!            Calculate Area_CornerFace dot Omega to determine the
!            contributions from incident fluxes across external
!            corner faces (FP faces)

             nCFaces = Geom% nCFacesArray(c0+c)

!            Contributions from interior corner faces (EZ faces)

             do cface=1,nCFaces

               aez = ASet% AezNorm(cface,c0+c) 
               cez = Geom% cEZ(cface,c0+c)

               if (aez > zero ) then

                 area_opp = zero
                 psi_opp  = zero

                 ifp = mod(cface,nCFaces) + 1
                 afp = ASet% AfpNorm(ifp,c0+c)

                 if ( afp < zero ) then
                   cfp      =  Geom% cFP(ifp,c0+c)
                   area_opp = -afp
                   psi_opp  = -afp*Set% Psi1(g,cfp)
                 endif

                 do i=2,nCFaces-2
                   ifp = mod(ifp,nCFaces) + 1
                   afp = ASet% AfpNorm(ifp,c0+c)
                   if ( afp < zero ) then
                     cfp      = Geom% cFP(ifp,c0+c)
                     area_opp = area_opp - afp
                     psi_opp  = psi_opp  - afp*Set% Psi1(g,cfp)
                   endif
                 enddo

                 TestOppositeFace: if ( area_opp > zero ) then

                   psi_opp  = psi_opp/area_opp
 
                   vol      = Geom% Volume(c0+c)
 
                   sigv     = sig*vol
                   sigv2    = sigv*sigv

                   gnum     = aez*aez*( fouralpha*sigv2 +              &
                              aez*(four*sigv + three*aez) )

                   gden     = vol*( four*sigv*sigv2 + aez*(six*sigv2 + &
                              two*aez*(two*sigv + aez)) )

                   sez      = ( vol*gnum*( sig*psi_opp - Set% Q(g,c,ii) ) +   &
                                half*aez*gden*( Set% Q(g,c,ii) - Set% Q(g,cez,ii) ) )/ &
                              ( gnum + gden*sig)

                   Set% S(g,c,ii)   = Set% S(g,c,ii)   + sez
                   Set% S(g,cez,ii) = Set% S(g,cez,ii) - sez

                 else

                   sez              = half*aez*( Set% Q(g,c,ii) - Set% Q(g,cez,ii) )/sig
                   Set% S(g,c,ii)   = Set% S(g,c,ii)   + sez
                   Set% S(g,cez,ii) = Set% S(g,cez,ii) - sez

                 endif TestOppositeFace

               endif

             enddo

           enddo GroupLoop1
         enddo CornerLoop1
       enddo ZoneLoop1

#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

#ifdef TETON_ENABLE_OPENACC
!$acc  loop vector collapse(2) &
!$acc& private(c0,cfp,ifp,cez,zone,zone0,nCorner,nCFaces) &
!$acc& private(aez,area_opp,sig,vol) &
!$acc& private(sigv,sigv2,sez,gnum,gden,psi_opp,afp)
#else
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, Geom, ASet, GSet, Angle, nzones, Groups, ndoneZ) &
!$omp& private(c0,cfp,ifp,cez,zone,zone0,nCFaces) &
!$omp& private(aez,area_opp,sig,vol) &
!$omp& private(sigv,sigv2,sez,gnum,gden,psi_opp,afp)
#endif

       ZoneLoop11: do ii=1,nzones
         GroupLoop11: do g=1,Groups

!          Loop through the zones using the NEXTZ list

           zone0   = ASet% nextZ(ndoneZ+ii,Angle)
           zone    = iabs( zone0 )

           CornerLoop11: do c=9,Geom% numCorner(zone)

             c0      = Geom% cOffSet(zone)
             sig     = GSet% Sigt(g,zone)

!            Calculate Area_CornerFace dot Omega to determine the
!            contributions from incident fluxes across external
!            corner faces (FP faces)

             nCFaces = Geom% nCFacesArray(c0+c)

!            Contributions from interior corner faces (EZ faces)

             do cface=1,nCFaces

               aez = ASet% AezNorm(cface,c0+c)
               cez = Geom% cEZ(cface,c0+c)

               if (aez > zero ) then

                 area_opp = zero
                 psi_opp  = zero

                 ifp = mod(cface,nCFaces) + 1
                 afp = ASet% AfpNorm(ifp,c0+c)

                 if ( afp < zero ) then
                   cfp      =  Geom% cFP(ifp,c0+c)
                   area_opp = -afp
                   psi_opp  = -afp*Set% Psi1(g,cfp)
                 endif

                 do i=2,nCFaces-2
                   ifp = mod(ifp,nCFaces) + 1
                   afp = ASet% AfpNorm(ifp,c0+c)
                   if ( afp < zero ) then
                     cfp      = Geom% cFP(ifp,c0+c)
                     area_opp = area_opp - afp
                     psi_opp  = psi_opp  - afp*Set% Psi1(g,cfp)
                   endif
                 enddo

                 TestOppositeFace1: if ( area_opp > zero ) then

                   psi_opp  = psi_opp/area_opp

                   vol      = Geom% Volume(c0+c)

                   sigv     = sig*vol
                   sigv2    = sigv*sigv

                   gnum     = aez*aez*( fouralpha*sigv2 +              &
                              aez*(four*sigv + three*aez) )

                   gden     = vol*( four*sigv*sigv2 + aez*(six*sigv2 + &
                              two*aez*(two*sigv + aez)) )

                   sez      = ( vol*gnum*( sig*psi_opp - Set% Q(g,c,ii) ) +   &
                                half*aez*gden*( Set% Q(g,c,ii) - Set% Q(g,cez,ii) ) )/ &
                              ( gnum + gden*sig)

                   Set% S(g,c,ii)   = Set% S(g,c,ii)   + sez
                   Set% S(g,cez,ii) = Set% S(g,cez,ii) - sez

                 else

                   sez              = half*aez*( Set% Q(g,c,ii) - Set% Q(g,cez,ii) )/sig
                   Set% S(g,c,ii)   = Set% S(g,c,ii)   + sez
                   Set% S(g,cez,ii) = Set% S(g,cez,ii) - sez

                 endif TestOppositeFace1

               endif

             enddo

           enddo CornerLoop11
         enddo GroupLoop11
       enddo ZoneLoop11

#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif


#ifdef TETON_ENABLE_OPENACC
!$acc  loop vector collapse(2) &
!$acc& private(c,c0,cez,zone,zone0,nCorner,nCFaces) &
!$acc& private(aez,sig)
#else
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, GSet, Geom, ASet, Angle, nzones, Groups, ndoneZ) &
!$omp& private(c,c0,cez,zone,zone0,nCorner,nCFaces) &
!$omp& private(aez,sig)
#endif

       ZoneLoop: do ii=1,nzones
         GroupLoop: do g=1,Groups

!          Loop through the zones using the NEXTZ list

           zone0   = ASet% nextZ(ndoneZ+ii,Angle)
           zone    = iabs( zone0 )
           nCorner = Geom% numCorner(zone)
           c0      = Geom% cOffSet(zone)

           sig     = GSet% Sigt(g,zone)


           if ( zone0 > 0 ) then

             do i=1,nCorner

               c = ASet% nextC(c0+i,angle) 

!              Corner angular flux
               Set% Psi1(g,c0+c) = Set% S(g,c,ii)/(ASet% ANormSum(c0+c) + sig*Geom% Volume(c0+c))

!              Calculate the contribution of this flux to the sources of
!              downstream corners in this zone. The downstream corner index is
!              "ez_exit."

               nCFaces = Geom% nCFacesArray(c0+c)

               do cface=1,nCFaces
                 aez = ASet% AezNorm(cface,c0+c)

                 if (aez > zero) then
                   cez              = Geom% cEZ(cface,c0+c) 
                   Set% S(g,cez,ii) = Set% S(g,cez,ii) + aez*Set% Psi1(g,c0+c)
                 endif
               enddo

             enddo

           else

!          Direct Solve (non-lower triangular, use old values of Psi1)
             do c=1,nCorner
               do cface=1,Geom% nCFacesArray(c0+c)
                 aez = ASet% AezNorm(cface,c0+c)

                 if (aez > zero) then
                   cez              = Geom% cEZ(cface,c0+c)
                   Set% S(g,cez,ii) = Set% S(g,cez,ii) + aez*Set% Psi1(g,c0+c)
                 endif
               enddo
             enddo

             do c=1,nCorner
               Set% Psi1(g,c0+c) = Set% S(g,c,ii)/(ASet% ANormSum(c0+c) + sig*Geom% Volume(c0+c))
             enddo

           endif

         enddo GroupLoop
       enddo ZoneLoop

#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

       ndoneZ = ndoneZ + nzones

     enddo HyperPlaneLoop

   enddo SetLoop

#ifdef TETON_ENABLE_OPENACC
!$acc end parallel loop
#else
TOMP(end target teams distribute)
#endif


!  Update Boundary data

#ifdef TETON_ENABLE_OPENACC
!$acc parallel loop gang num_gangs(nSets) &
!$acc& vector_length(omp_device_team_thread_limit) &
!$acc& private(Set, ASet, BdyExitPtr, offSet, Angle, Groups, b, c)
#else
TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
TOMPC(shared(nSets, Quad, sendIndex)&)
TOMPC(private(Set, ASet, BdyExitPtr, offSet, Angle, Groups, b, c))
#endif

     SetLoop3: do setID=1,nSets

       Set        => Quad% SetDataPtr(setID)
       ASet       => Quad% AngSetPtr(Set% angleSetID)
       Groups     =  Set% Groups
       Angle      =  Set% AngleOrder(sendIndex)
       offSet     =  ASet% cycleOffSet(angle)
       BdyExitPtr => ASet% BdyExitPtr(Angle)

#ifdef TETON_ENABLE_OPENACC
!$acc  loop vector collapse(2) &
!$acc& private(b,c)
#else
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, BdyExitPtr, Groups, Angle) private(b,c)
#endif

       do i=1,BdyExitPtr% nxBdy
         do g=1,Groups
           b = BdyExitPtr% bdyList(1,i)
           c = BdyExitPtr% bdyList(2,i)

           Set% PsiB(g,b,Angle) = Set% Psi1(g,c)
         enddo
       enddo

#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

!      Update Psi in the cycle list

#ifdef TETON_ENABLE_OPENACC
!$acc  loop vector collapse(2) &
!$acc& private(c)
#else
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Angle, Set, ASet, offSet, Groups) private(c)
#endif
       do mCycle=1,ASet% numCycles(angle)
         do g=1,Groups
           c                              = ASet% cycleList(offSet+mCycle)
           Set% cyclePsi(g,offSet+mCycle) = Set% Psi1(g,c)
         enddo
       enddo
#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

     enddo SetLoop3

#ifdef TETON_ENABLE_OPENACC
!$acc end parallel loop
#else
TOMP(end target teams distribute)
#endif

!  We only store Psi if this is the last transport sweep in the time step

   if ( savePsi ) then

#ifdef TETON_ENABLE_OPENACC
!$acc parallel loop gang num_gangs(nSets) &
!$acc& vector_length(omp_device_team_thread_limit) &
!$acc& private(Set, Angle, Groups)
#else
TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none)&)
TOMPC(shared(nSets, Quad, sendIndex)&)
TOMPC(private(Set, Angle, Groups))
#endif

     SetLoop2: do setID=1,nSets

       Set    => Quad% SetDataPtr(setID)
       Groups =  Set% Groups
       Angle  =  Set% AngleOrder(sendIndex)

#ifdef TETON_ENABLE_OPENACC
!$acc  loop vector collapse(2)
#else
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, Angle, Groups)
#endif

       CornerLoop2: do c=1,Set% nCorner
         GroupLoop2: do g=1,Groups

           Set% Psi(g,c,Angle) = Set% Psi1(g,c)

         enddo GroupLoop2
       enddo CornerLoop2

#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

     enddo SetLoop2

#ifdef TETON_ENABLE_OPENACC
!$acc end parallel loop
#else
TOMP(end target teams distribute)
#endif

   endif

#ifdef TETON_ENABLE_OPENACC
!$acc end data
#else
TOMP(target exit data map(release: tau, sendIndex, angleList))
#endif


   deallocate( angleList )

   return
   end subroutine CornerSweepUCBxyz_GPU

