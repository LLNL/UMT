#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   CornerSweepUCBrz_GPU - This routine calculates angular fluxes for a      *
!                    single direction and multiple energy groups for   *
!                    for an upstream corner-balance (UCB) spatial      *
!                    in rz-geometry. The mesh is traversed corner
!                    by corner.                                        *
!                                                                      *
!                    The work is offloaded to a GPU in which each      *
!                    computational "set" is a GPU block. The threads   *
!                    assigned to the block compute one group in one    *
!                    zone in a hyperplane (by definition, all of the   *
!                    zones in a hyperplane are independent).           *
!                                                                      *
!***********************************************************************

   subroutine CornerSweepUCBrz_GPU(nSets, sendIndex, savePsi)

   use, intrinsic :: iso_c_binding, only : c_int
   use cmake_defines_mod, only : omp_device_team_thread_limit
   use Options_mod
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use GroupSet_mod
   use CodeChecks_mod

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

   integer    :: setID
   integer    :: zSetID
   integer    :: Angle
   integer    :: g
   integer    :: Groups

   integer    :: mCycle
   integer    :: offset
   integer    :: nAngleSets
   integer    :: nZoneSets
   integer    :: nHyperDomains

   integer    :: nzones
   integer    :: ii
   integer    :: ndone
   integer    :: hyperPlane
   integer    :: domID
   integer    :: hplane1
   integer    :: hplane2

   real(adqt) :: tau

!  Local

   integer    :: b
   integer    :: i
   integer    :: cface
   integer    :: cez
   integer    :: cfp
   integer    :: ifp
   integer    :: c
   integer    :: c0
   integer    :: c1
   integer    :: c2
   integer    :: zone
   integer    :: nCorner

   real(adqt), parameter :: fouralpha=1.82d0

   real(adqt) :: fac
   real(adqt) :: sigA
   real(adqt) :: sigA2
   real(adqt) :: source
   real(adqt) :: SS
   real(adqt) :: QQ
   real(adqt) :: Qez
   real(adqt) :: mult

   real(adqt) :: area
   real(adqt) :: sig
   real(adqt) :: sez
   real(adqt) :: gnum
   real(adqt) :: gden
   real(adqt) :: quadTauW1
   real(adqt) :: quadTauW2

   real(adqt) :: afp
   real(adqt) :: aez
   real(adqt) :: R
   real(adqt) :: R2
   real(adqt) :: R_afp
   real(adqt) :: R_afp2

!  Dynamic

   integer, allocatable :: angleList(:)

!  Constants
   tau           = Size% tau
   nAngleSets    = getNumberOfAngleSets(Quad)
   nZoneSets     = getNumberOfZoneSets(Quad)
   nHyperDomains = getNumberOfHyperDomains(Quad,1)

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

   ! Verify we won't get out-of-bounds accesses below.
   TETON_CHECK_BOUNDS1(Quad%AngSetPtr, nAngleSets)
   TETON_CHECK_BOUNDS1(Geom%corner1, nZoneSets)
   TETON_CHECK_BOUNDS1(Geom%corner2, nZoneSets)

   TOMP(target enter data map(to: tau, sendIndex, angleList))


#ifdef TETON_ENABLE_OPENACC
   !$acc parallel loop gang num_gangs(nZoneSets) vector_length(omp_device_team_thread_limit) &
   !$acc& private(ASet, angle)
#else
   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, nAngleSets,Geom, angleList, Quad)&)
   TOMPC(private(ASet, angle))
#endif

   ZoneSetLoop: do zSetID=1,nZoneSets

!    Loop over angle sets

     do setID=1,nAngleSets

       ASet  => Quad% AngSetPtr(setID)
       angle =  angleList(setID)

#ifdef TETON_ENABLE_OPENACC
       !$acc loop vector collapse(2)
#else
       !$omp  parallel do collapse(2) default(none)  &
       !$omp& shared(Geom, ASet, Angle, zSetID)
#endif
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         do cface=1,2
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
   !$acc parallel loop gang num_gangs(nZoneSets) vector_length(omp_device_team_thread_limit) &
   !$acc& private(ASet, angle, fac, R_afp, R_afp2, R, R2)
#else
   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
   TOMPC(shared(nZoneSets, nAngleSets, Geom, angleList, Quad)&)
   TOMPC(private(ASet, angle, fac, R_afp, R_afp2, R, R2))
#endif

   ZoneSetLoop2: do zSetID=1,nZoneSets

!    Loop over angle sets

     do setID=1,nAngleSets

       ASet  => Quad% AngSetPtr(setID)
       angle =  angleList(setID)
       fac   =  ASet% angDerivFac(Angle)

#ifdef TETON_ENABLE_OPENACC
       !$acc parallel loop gang num_gangs(nZoneSets) vector_length(omp_device_team_thread_limit) &
       !$acc& private(R_afp, R_afp2, R, R2)
#else
       !$omp  parallel do default(none)  &
       !$omp& shared(Geom, ASet, angle, fac, zSetID) private(R_afp, R_afp2, R, R2)
#endif

       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         R_afp  = Geom% RadiusFP(1,c)
         R_afp2 = Geom% RadiusFP(2,c)
         R      = Geom% RadiusEZ(1,c)
         R2     = Geom% RadiusEZ(2,c)

         ASet% ANormSum(c) = fac*Geom% Area(c) - half*(  &
            R_afp *(ASet% AfpNorm(1,c) - abs(ASet% AfpNorm(1,c))) +    &
            R_afp2*(ASet% AfpNorm(2,c) - abs(ASet% AfpNorm(2,c))) +    &
            R     *(ASet% AezNorm(1,c) - abs(ASet% AezNorm(1,c))) +    &
            R2    *(ASet% AezNorm(2,c) - abs(ASet% AezNorm(2,c))) ) 
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
   !$acc& private(Set, ASet, HypPlanePtr, Angle, Groups, offSet, c)
#else
   TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(sendIndex, Quad, nSets) &)
   TOMPC(private(Set, ASet, HypPlanePtr, Angle, Groups, offSet, c))
#endif

   SetLoop0: do setID=1,nSets

     Set          => Quad% SetDataPtr(setID)
     ASet         => Quad% AngSetPtr(Set% angleSetID)

     Groups       =  Set% Groups
     Angle        =  Set% AngleOrder(sendIndex)
     offSet       =  ASet% cycleOffSet(angle)
     HypPlanePtr  => ASet% HypPlanePtr(angle)

!  Initialize boundary values in Psi1 and interior values on the cycle
!  list

#ifdef TETON_ENABLE_OPENACC
     !$acc  loop vector collapse(2) &
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
     do b=1,Set%nbelem
       do g=1,Groups
         Set% Psi1(g,Set%nCorner+b) = Set% PsiB(g,b,Angle)
       enddo
     enddo
#ifndef TETON_ENABLE_OPENACC
     !$omp end parallel do
#endif

!    Initialize values at hyper-domain interfaces

#ifdef TETON_ENABLE_OPENACC
     !$acc loop vector collapse(2) &
     !$acc& private(c)
#else
     !$omp  parallel do collapse(2) default(none) &
     !$omp& shared(Set, HypPlanePtr, Groups, angle) private(c)
#endif
     do b=1,HypPlanePtr% interfaceLen
       do g=1,Groups
         c              = HypPlanePtr% interfaceList(b)
         Set% Psi1(g,c) = Set% PsiInt(g,b,angle)
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
! IBM XLF segfaults if 'mCycle', 'b', and 'g' are not scoped to private below.  This should not
! be necessary, as these are loop control variables which the runtime should automatically scope to private.
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
   !$acc  parallel loop gang collapse(2) num_gangs(nSets*nHyperDomains) &
   !$acc& vector_length(omp_device_team_thread_limit) &
   !$acc& private(Set, ASet, GSet, HypPlanePtr, Swp, Angle, Groups, hplane1, hplane2, ndone) &
   !$acc& private(hyperPlane, nzones, fac, c, c0, c1, c2, cfp, cez, zone, nCorner) &
   !$acc& private(sigA, sigA2, source, area, sig, sez, SS, QQ, Qez) &
   !$acc& private(gnum, gden, aez, afp, R, R_afp, ifp, mult)
#else
   TOMP(target teams distribute collapse(2) num_teams(nSets*nHyperDomains) default(none) &)
   TOMPC(thread_limit(omp_device_team_thread_limit) &)
   TOMPC(shared(nSets, nHyperDomains, Quad, Geom, sendIndex, tau)&)
   TOMPC(private(Set, ASet, GSet, HypPlanePtr, Angle, Groups, hplane1, hplane2, ndone) &)
   TOMPC(private(b, g, hyperPlane, nzones, fac, c, c0, c1, c2, cfp, cez, zone, nCorner)&)
   TOMPC(private(sigA, sigA2, source, area, sig, sez, SS, QQ, Qez)&)
   TOMPC(private(gnum, gden, aez, afp, R, R_afp, ifp, mult))
#endif

   SetLoop: do setID=1,nSets
     DomainLoop: do domID=1,nHyperDomains

     Set          => Quad% SetDataPtr(setID)
     ASet         => Quad% AngSetPtr(Set% angleSetID)
     GSet         => Quad% GrpSetPtr(Set% groupSetID)

     Groups       =  Set% Groups
     Angle        =  Set% AngleOrder(sendIndex)
     HypPlanePtr  => ASet% HypPlanePtr(Angle)
     hplane1      =  HypPlanePtr% hplane1(domID)
     hplane2      =  HypPlanePtr% hplane2(domID)
     ndone        =  HypPlanePtr% ndone(domID)

!    Angle Constants

     fac          = ASet% angDerivFac(Angle)

     HyperPlaneLoop: do hyperPlane=hplane1,hplane2

       nCorner = HypPlanePtr% cornersInPlane(hyperPlane)

#ifdef TETON_ENABLE_OPENACC
       !$acc  loop vector collapse(2) &
       !$acc& private(c, c0, c1, c2, cez, cfp, cface, ifp, zone, sigA, sigA2, source) &
       !$acc& private(SS, QQ, Qez, area, sig, sez, gnum, gden, aez, afp, R, R_afp)
#else
       !$omp  parallel do collapse(2) default(none) &
       !$omp& shared(Set, Geom, ASet, GSet, Angle, nCorner, Groups, ndone, tau, fac) &
       !$omp& private(c, c0, c1, c2, cez, cfp, cface, ifp, zone, sigA, sigA2, source, mult) &
       !$omp& private(SS, QQ, Qez, area, sig, sez, gnum, gden, aez, afp, R, R_afp)
#endif

       CornerLoop: do ii=1,nCorner
         GroupLoop: do g=1,Groups

!  Loop through all of the corners using the NEXT list

           c       = ASet% nextC(ndone+ii,Angle)
           zone    = Geom% CToZone(c)
           c0      = Geom% cOffSet(zone)

           sig     = GSet% Sigt(g,zone)

!  Contributions from volume terms (if a starting direction add angular
!  derivative)

           source  = GSet% STotal(g,c) + tau*Set% Psi(g,c,Angle)
           SS      = Geom% Volume(c)*source + fac*Geom% Area(c)*Set% PsiM(g,c) 

           do cface=1,2
             afp   = ASet% AfpNorm(cface,c)
             cfp   = Geom% cFP(cface,c)
             R_afp = Geom% RadiusFP(cface,c)*afp

             if ( afp < zero ) then
               SS = SS - R_afp*Set% Psi1(g,cfp)
             endif
           enddo

           CornerFaceLoop: do cface=1,2

             aez = ASet% AezNorm(cface,c) 

             if (aez > zero ) then

               c1   = c
               c2   = c0 + Geom% cEZ(cface,c)
               QQ   = source
               Qez  = GSet% STotal(g,c2) + tau*Set% Psi(g,c2,Angle)
               mult = one
               ifp  = cface

             elseif (aez < zero ) then

               c2   = c
               c1   = c0 + Geom% cEZ(cface,c)
               QQ   = GSet% STotal(g,c1) + tau*Set% Psi(g,c1,Angle)
               Qez  = source
               mult = -one
               ifp  = mod(cface,2) + 1
               R    = Geom% RadiusEZ(cface,c)

!              Contributions from upsteam fluxes in the same zone
               SS   = SS - R*aez*Set% Psi1(g,c1)
               aez  = -aez

             endif

             afp = ASet% AfpNorm(ifp,c1)
             R   = Geom% RadiusEZ(ifp,c1)

             if ( afp < zero ) then

               cfp   = Geom% cFP(ifp,c1)
               area  = Geom% Area(c1)
               sigA  = sig*area
               sigA2 = sigA*sigA

               gnum  = aez*aez*( fouralpha*sigA2 +  &
                       aez*(four*sigA + three*aez) )

               gden  = area*(four*sigA*sigA2 + aez*(six*sigA2 + &
                       two*aez*(two*sigA + aez)))

               sez   = R*( area*gnum*( sig*Set% Psi1(g,cfp) - QQ ) +  &
                           half*aez*gden*( QQ - Qez ) )/  &
                         ( gnum + gden*sig )

               SS    = SS + mult*sez

             else

               sez   = half*R*aez*( QQ - Qez )/sig
               SS    = SS + mult*sez

             endif

           enddo CornerFaceLoop

!          Corner angular flux
           Set% Psi1(g,c) = SS/(ASet% ANormSum(c) + sig*Geom% Volume(c))

         enddo GroupLoop
       enddo CornerLoop
#ifndef TETON_ENABLE_OPENACC
       !$omp end parallel do 
#endif

       ndone = ndone + nCorner 

     enddo HyperPlaneLoop

     enddo DomainLoop
   enddo SetLoop

#ifdef TETON_ENABLE_OPENACC
   !$acc end parallel loop
#else
   TOMP(end target teams distribute)
#endif


#ifdef TETON_ENABLE_OPENACC
     !$acc  parallel loop gang num_gangs(nSets) vector_length(omp_device_team_thread_limit) &
     !$acc& private(Set, ASet, Angle, Groups, quadTauW1, quadTauW2)
#else
     TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nSets, Quad, sendIndex)&)
     TOMPC(private(Set, ASet, Angle, Groups, quadTauW1, quadTauW2))
#endif

     SetLoop2: do setID=1,nSets

       Set    => Quad% SetDataPtr(setID)
       ASet   => Quad% AngSetPtr(Set% angleSetID)

       Groups =  Set% Groups
       Angle  =  Set% AngleOrder(sendIndex)

!      Set the "half-angle" angular intensity (PsiM) for the next angle

       if ( ASet% StartingDirection(Angle) ) then

#ifdef TETON_ENABLE_OPENACC
         !$acc loop vector collapse(2)
#else
         !$omp  parallel do collapse(2) default(none) &
         !$omp& shared(Set, Groups)
#endif
         do c=1,Set% nCorner
           do g=1,Groups
             Set% PsiM(g,c) = Set% Psi1(g,c)
           enddo 
         enddo 
#ifndef TETON_ENABLE_OPENACC
         !$omp end parallel do
#endif

       else

         quadTauW1 = ASet% quadTauW1(Angle)
         quadTauW2 = ASet% quadTauW2(Angle)

#ifdef TETON_ENABLE_OPENACC
         !$acc loop vector collapse(2)
#else
         !$omp  parallel do collapse(2) default(none) &
         !$omp& shared(Set, Groups, quadTauW1, quadTauW2)
#endif

         do c=1,Set% nCorner
           do g=1,Groups
             Set% PsiM(g,c) = quadTauW1*Set% Psi1(g,c) -  &
                              quadTauW2*Set% PsiM(g,c)
           enddo 
         enddo
#ifndef TETON_ENABLE_OPENACC
         !$omp end parallel do
#endif

       endif

     enddo SetLoop2
#ifdef TETON_ENABLE_OPENACC
     !$acc end parallel loop
#else
     TOMP(end target teams distribute)
#endif


#ifdef TETON_ENABLE_OPENACC
     !$acc  parallel loop gang num_gangs(nSets) vector_length(omp_device_team_thread_limit) &
     !$acc& private(Set, ASet, BdyExitPtr, HypPlanePtr, offSet, Angle, Groups, b, c)
#else
     TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none)&)
     TOMPC(shared(nSets, Quad, sendIndex)&)
     TOMPC(private(Set, ASet, BdyExitPtr, HypPlanePtr, offSet, Angle, Groups, b, c))
#endif

     SetLoop3: do setID=1,nSets

       Set         => Quad% SetDataPtr(setID)
       ASet        => Quad% AngSetPtr(Set% angleSetID)
       Groups      =  Set% Groups
       Angle       =  Set% AngleOrder(sendIndex)
       offSet      =  ASet% cycleOffSet(angle)
       BdyExitPtr  => ASet% BdyExitPtr(Angle)
       HypPlanePtr => ASet% HypPlanePtr(angle)

#ifdef TETON_ENABLE_OPENACC
       !$acc  loop vector collapse(2) &
       !$acc& private(b, c)
#else
       !$omp  parallel do collapse(2) default(none) &
       !$omp& shared(Set, BdyExitPtr, Angle, Groups) private(b, c)
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

!      Update Interface Elements

#ifdef TETON_ENABLE_OPENACC
       !$acc  loop vector collapse(2) &
       !$acc& private(c)
#else
       !$omp  parallel do collapse(2) default(none) &
       !$omp& shared(Set, HypPlanePtr, Groups, angle) private(c)
#endif

       do i=1,HypPlanePtr% interfaceLen
         do g=1,Groups
           c = HypPlanePtr% interfaceList(i)
           Set% PsiInt(g,i,angle) = Set% Psi1(g,c)
         enddo
       enddo

#ifndef TETON_ENABLE_OPENACC
       !$omp end parallel do
#endif

!    Update Psi in the cycle list

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


       if ( ASet% FinishingDirection(Angle+1) ) then

#ifdef TETON_ENABLE_OPENACC
         !$acc  loop vector collapse(2) &
         !$acc& private(b, c)
#else
         !$omp  parallel do collapse(2) default(none) &
         !$omp& shared(Set, BdyExitPtr, Angle, Groups) private(b, c)
#endif
         do i=1,BdyExitPtr% nxBdy
           do g=1,Groups
             b = BdyExitPtr% bdyList(1,i)
             c = BdyExitPtr% bdyList(2,i)

             Set% PsiB(g,b,Angle+1) = Set% PsiM(g,c)
           enddo
         enddo
#ifndef TETON_ENABLE_OPENACC
!$omp end parallel do
#endif

       endif 

     enddo SetLoop3

#ifdef TETON_ENABLE_OPENACC
     !$acc end parallel loop
#else
     TOMP(end target teams distribute)
#endif


!  We only store Psi if this is the last transport sweep in the time step

   if ( savePsi ) then

#ifdef TETON_ENABLE_OPENACC
     !$acc  parallel loop gang num_gangs(nSets) vector_length(omp_device_team_thread_limit) &
     !$acc& private(Set, ASet, Angle, Groups)
#else
     TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nSets, sendIndex, Quad)&)
     TOMPC(private(Set, ASet, Angle, Groups))
#endif
     SetLoop4: do setID=1,nSets

       Set    => Quad% SetDataPtr(setID)
       ASet   => Quad% AngSetPtr(Set% angleSetID)

       Groups =  Set% Groups
       Angle  =  Set% AngleOrder(sendIndex)

#ifdef TETON_ENABLE_OPENACC
       !$acc loop vector collapse(2)
#else
       !$omp  parallel do collapse(2) default(none) &
       !$omp& shared(Set, ASet, Angle, Groups)
#endif
       CornerLoop4: do c=1,Set% nCorner
         GroupLoop4: do g=1,Groups

           Set% Psi(g,c,Angle) = Set% Psi1(g,c)

           if ( ASet% FinishingDirection(Angle+1) ) then
             Set% Psi(g,c,Angle+1) = Set% PsiM(g,c)
           endif

         enddo GroupLoop4
       enddo CornerLoop4
#ifndef TETON_ENABLE_OPENACC
       !$omp end parallel do
#endif

     enddo SetLoop4

#ifdef TETON_ENABLE_OPENACC
     !$acc end parallel loop
#else
     TOMP(end target teams distribute)
#endif

   endif


   TOMP(target exit data map(always,release: tau, sendIndex, angleList))

   deallocate( angleList )


   return
   end subroutine CornerSweepUCBrz_GPU

