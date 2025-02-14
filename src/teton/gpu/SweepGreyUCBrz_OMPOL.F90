#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   SweepGreyUCBrz  - This routine calculates angular fluxes for a     *
!                     single direction and single energy groups for    *
!                     for an upstream corner-balance (UCB) spatial     *
!                     in rz-geometry. It is only used for Grey         *
!                     Transport Acceleration (GTA).                    *
!                                                                      *
!***********************************************************************

   subroutine SweepGreyUCBrzNEW_GPU(sendIndex, PsiB)

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use SetData_mod
   use AngleSet_mod
   use CodeChecks_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: sendIndex 
   real(adqt), intent(inout) :: PsiB(Size%nSurfElem,Size%nangGTA)

!  Local

   type(SetData),    pointer :: Set
   type(AngleSet),   pointer :: ASet
   type(HypPlane),   pointer :: HypPlanePtr
   type(BdyExit),    pointer :: BdyExitPtr

   integer    :: domID
   integer    :: setID
   integer    :: zSetID
   integer    :: nSets
   integer    :: nGTASets
   integer    :: nZoneSets
   integer    :: nHyperDomains

   integer    :: c
   integer    :: cc
   integer    :: c0
   integer    :: nCorner
   integer    :: ii
   integer    :: angle
   integer    :: angle0
   integer    :: angle2
   integer    :: zone
   integer    :: nzones
   integer    :: ndoneZ
   integer    :: hplane1
   integer    :: hplane2
   integer    :: hyperPlane
   integer    :: b
   integer    :: i
   integer    :: cfp
   integer    :: cez
   integer    :: cface

   real(adqt) :: fac
   real(adqt) :: quadwt
   real(adqt) :: quadTauW1
   real(adqt) :: quadTauW2
   real(adqt) :: sigA
   real(adqt) :: sigA2
   real(adqt) :: sez
   real(adqt) :: gnum
   real(adqt) :: gtau
   real(adqt) :: aez
   real(adqt) :: afp
   real(adqt) :: R
   real(adqt) :: R_afp
   real(adqt) :: dInv

   real(adqt), parameter :: fouralpha=1.82d0

   logical(kind=1) :: FinishingDirection

!  Dynamic

   integer, allocatable :: angleList(:,:)

!  Constants

   nSets         = getNumberOfSets(Quad)
   nGTASets      = getNumberOfGTASets(Quad)
   nZoneSets     = getNumberOfZoneSets(Quad)
   nHyperDomains = getNumberOfHyperDomains(Quad,2)

!  Quiet the compiler 'variable may not be initialized' warning.
   FinishingDirection = .FALSE.

   allocate( angleList(3,nGTASets) )

   do setID=1,nGTASets
     Set  => Quad% SetDataPtr(nSets+setID)
     ASet => Quad% AngSetPtr(Set% angleSetID)

     angle              = Set% AngleOrder(sendIndex)
     angleList(1,setID) = Set% AngleOrder(sendIndex)
     angleList(2,setID) = Set% angle0

     ! Compute the global non-finishing angle index:
     if (nGTASets == 1) then
        ! This else-logic is fine if there is exactly one finishing direction per set.
        ! But it fails if there's only one gta set.  In that case, we have
        !    angles 4 and 8 as finishing directions.

        ! Assume that the finishing directions are indexed ASet%numAngles/2 and
        !    ASet%numAngles
        if (angle > ASet%numAngles/2) then
           angleList(3,setID) = angle-1
        else
           angleList(3,setID) = angle
        endif
     else
        angleList(3,setID) = Set% angle0 + angle + 1 - setID
     endif

     FinishingDirection = ASet% FinishingDirection(angle)
   enddo

   if ( FinishingDirection ) then
     deallocate( angleList )

     return
   endif


   ! Verify we won't get out-of-bounds accesses below.
   TETON_CHECK_BOUNDS1(Quad%SetDataPtr, nSets+nGTASets)
   TETON_CHECK_BOUNDS1(Geom%corner1, nZoneSets)
   TETON_CHECK_BOUNDS1(Geom%corner2, nZoneSets)

   TOMP_MAP(target enter data map(to: angleList, nSets, nGTASets, nHyperDomains, PsiB))

   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets,nGTASets,Geom,nSets,Quad )&)
   TOMPC(private(Set))

   ZoneSetLoop0: do zSetID=1,nZoneSets

     do setID=1,nGTASets
       Set => Quad% SetDataPtr(nSets+setID)

       !$omp  parallel do default(none)  &
       !$omp& shared(Set, Geom, zSetID)

       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         Set% tPsi(c) = zero
       enddo

       !$omp end parallel do

     enddo

   enddo ZoneSetLoop0

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nGTASets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nGTASets, nSets,PsiB, Quad, angleList)&)
   TOMPC(private(Set, ASet, angle0, angle, c))

   GTASetLoop0: do setID=1,nGTASets

     Set         => Quad% SetDataPtr(nSets+setID)
     ASet        => Quad% AngSetPtr(Set% angleSetID)
     angle       =  angleList(1,setID)
     angle0      =  angleList(2,setID)

!    Initialize Boundary Values

     !$omp  parallel do default(none) &
     !$omp& shared(Set, PsiB, angle0, angle)

     do b=1,Set%nbelem
       Set% tPsi(Set%nCorner+b) = PsiB(b,angle0+angle)
     enddo

     !$omp end parallel do

   enddo GTASetLoop0

   TOMP(end target teams distribute)

   if ( nHyperDomains > 1 ) then

     TOMP(target teams distribute collapse(2) num_teams(nZoneSets*nGTASets) &)
     TOMPC(thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nZoneSets, nGTASets, Quad, nSets, angleList, PsiB)&)
     TOMPC(private(Set, ASet, HypPlanePtr, angle0, angle, cc))

     ZoneSetLoop1: do zSetID=1,nZoneSets
       do setID=1,nGTASets
         Set         => Quad% SetDataPtr(nSets+setID)
         ASet        => Quad% AngSetPtr(Set% angleSetID)
         angle       =  angleList(1,setID)
         angle0      =  angleList(2,setID)
         HypPlanePtr => ASet% HypPlanePtr(angle)

!        Initialize values at hyper-domain interfaces

         !$omp  parallel do default(none) &
         !$omp& shared(Set, HypPlanePtr, PsiB, angle0, angle, zSetID) private(cc)

         do i=HypPlanePtr% c1(zSetID),HypPlanePtr% c2(zSetID)
           cc = HypPlanePtr% interfaceList(i)
           Set% tPsi(cc) = PsiB(Set%nbelem+i,angle0+angle)
         enddo

     !$omp end parallel do

       enddo
     enddo ZoneSetLoop1

     TOMP(end target teams distribute)

   endif

   TOMP(target teams distribute collapse(2) num_teams(nHyperDomains*nGTASets) &)
   TOMPC(thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nGTASets, nHyperDomains, angleList, Geom, GTA, Quad, nSets)&)
   TOMPC(private(Set, ASet, HypPlanePtr, setID, domID, angle, angle0, angle2, hyperPlane) &)
   TOMPC(private(hplane1, hplane2, ndoneZ, nzones, fac, zone, nCorner, c0, cez, cfp) &)
   TOMPC(private(aez, afp, R, R_afp, sigA, sigA2, gnum, gtau, sez, dInv))

   GTASetLoop: do setID=1,nGTASets
     DomainLoop: do domID=1,nHyperDomains

       Set          => Quad% SetDataPtr(nSets+setID)
       ASet         => Quad% AngSetPtr(Set% angleSetID)
       angle        =  angleList(1,setID)
       angle0       =  angleList(2,setID)
       angle2       =  angleList(3,setID) 

       HypPlanePtr  => ASet% HypPlanePtr(angle)
       hplane1      =  HypPlanePtr% hplane1(domID)
       hplane2      =  HypPlanePtr% hplane2(domID)
       fac          =  ASet% angDerivFac(angle)
       ndoneZ       =  HypPlanePtr% ndone(domID)

!      Loop through all of the zones using the NEXT list

       HyperPlaneLoop: do hyperPlane=hplane1,hplane2

         nzones = HypPlanePtr% zonesInPlane(hyperPlane)

         !$omp  parallel do default(none) schedule(static)   &
         !$omp& shared(Geom, GTA, ASet, Set, Quad)           &
         !$omp& shared(ndoneZ, nzones, angle, angle2, fac)   &
         !$omp& private(zone, nCorner, c, c0, cez, cfp, aez, afp) &
         !$omp& private(R, R_afp, sigA, sigA2, gnum, gtau, sez, dInv)

         ZoneLoop: do ii=1,nzones

           zone    = iabs( ASet% nextZ(ndoneZ+ii,angle) )
           nCorner = Geom% numCorner(zone)
           c0      = Geom% cOffSet(zone)

!          Contributions from volume terms and angular derivative

           do c=1,nCorner
             Set% src(c0+c)  = GTA%TsaSource(c0+c) + &
                               fac*Geom% Area(c0+c)*Set% tPsiM(c0+c)
             Set% pInc(c0+c) = fac*Geom% Area(c0+c)*Set% tInc(c0+c)
           enddo

           CornerLoop: do c=1,nCorner

             sigA = GTA%GreySigTotal(c0+c)*Geom% Area(c0+c)

             CornerFaceLoop: do cface=1,2

               afp = GTA% AfpNorm(cface,c0+c,angle2)
               aez = GTA% AezNorm(cface,c0+c,angle2)

               if ( afp < zero ) then
                 cfp             = Geom% cFP(cface,c0+c)
                 R_afp           = Geom% RadiusFP(cface,c0+c)*afp
                 Set% src(c0+c)  = Set% src(c0+c)  - R_afp*Set% tPsi(cfp)
                 Set% pInc(c0+c) = Set% pInc(c0+c) - R_afp*Set% tPsi(cfp)
               else
                 cfp = -1
               endif

               if ( aez > zero ) then

                 R   = Geom% RadiusEZ(cface,c0+c)
                 cez = Geom% cEZ(cface,c0+c)

                 if ( afp < zero ) then

                   sigA2 = sigA*sigA

                   gnum  = aez*aez*( fouralpha*sigA2 +  &
                           aez*(four*sigA + three*aez) )

                   gtau  = gnum/                                           &
                         ( gnum + four*sigA2*sigA2 + aez*sigA*(six*sigA2 + &
                           two*aez*(two*sigA + aez)) )

                   sez   = R*( gtau*sigA*( Set% tPsi(cfp) - GTA% Q(c0+c) ) +  &
                           half*aez*(one - gtau)*( GTA% Q(c0+c) - GTA% Q(c0+cez) ) )

                   Set% src(c0+c)    = Set% src(c0+c)    + sez
                   Set% src(c0+cez)  = Set% src(c0+cez)  - sez

                   Set% pInc(c0+c)   = Set% pInc(c0+c)   + R*gtau*sigA*Set% tPsi(cfp)
                   Set% pInc(c0+cez) = Set% pInc(c0+cez) - R*gtau*sigA*Set% tPsi(cfp)

                 else

                   sez               = half*R*aez*( GTA% Q(c0+c) - GTA% Q(c0+cez) )
                   Set% src(c0+c)    = Set% src(c0+c)   + sez
                   Set% src(c0+cez)  = Set% src(c0+cez) - sez

                 endif

               endif

             enddo CornerFaceLoop

           enddo CornerLoop

!          Solve the corners in the correct order

           do i=1,nCorner

             c    = ASet% nextC(c0+i,angle) 

             dInv = one/(GTA% ANormSum(c0+c,angle2) +   &
                    Geom% Volume(c0+c)*GTA%GreySigTotal(c0+c))

!            Corner angular flux
             Set% tPsi(c0+c) = dInv*Set% src(c0+c)
             Set% pInc(c0+c) = dInv*Set% pInc(c0+c)

!            Calculate the contribution of this flux to the sources of
!            downstream corners in this zone. 

             do cface=1,2
               aez = GTA% AezNorm(cface,c0+c,angle2)

               if (aez > zero) then
                 R                 = Geom% RadiusEZ(cface,c0+c)
                 cez               = Geom% cEZ(cface,c0+c)

                 Set% src(c0+cez)  = Set% src(c0+cez)  + R*aez*Set% tPsi(c0+c)
                 Set% pInc(c0+cez) = Set% pInc(c0+cez) + R*aez*Set% pInc(c0+c)
               endif
             enddo

           enddo

         enddo ZoneLoop

         !$omp end parallel do

         ndoneZ = ndoneZ + nzones

       enddo HyperPlaneLoop

     enddo DomainLoop
   enddo GTASetLoop

   TOMP(end target teams distribute)

!  Update exiting boundary fluxes 

   TOMP(target teams distribute num_teams(nGTASets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nGTASets, PsiB, Quad, angleList, nSets)&)
   TOMPC(private(Set, ASet, BdyExitPtr, HypPlanePtr, angle0, angle, b, c))

   do setID=1,nGTASets

     Set         => Quad% SetDataPtr(nSets+setID)
     ASet        => Quad% AngSetPtr(Set% angleSetID)
     angle       =  angleList(1,setID)
     angle0      =  angleList(2,setID)
     BdyExitPtr  => ASet% BdyExitPtr(angle)
     HypPlanePtr => ASet% HypPlanePtr(angle)

     !$omp  parallel do default(none) &
     !$omp& shared(Set, BdyExitPtr, angle0, angle, PsiB) &
     !$omp& private(b,c)

     do i=1,BdyExitPtr% nxBdy
       b = BdyExitPtr% bdyList(1,i)
       c = BdyExitPtr% bdyList(2,i)

       PsiB(b,angle0+angle) = Set% tPsi(c)
     enddo

     !$omp end parallel do

   enddo

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
   TOMPC(shared(nZoneSets, nGTASets, Geom, GTA, Quad, angleList, nSets)&)
   TOMPC(private(Set, ASet, angle, angle0, quadwt, quadTauW1, quadTauW2))

   ZoneSetLoop2: do zSetID=1,nZoneSets
     do setID=1,nGTASets

       Set         => Quad% SetDataPtr(nSets+setID)
       ASet        => Quad% AngSetPtr(Set% angleSetID)

       angle       =  angleList(1,setID)
       angle0      =  angleList(2,setID)
       quadwt      =  ASet% weight(angle)
       quadTauW1   =  ASet% quadTauW1(angle)
       quadTauW2   =  ASet% quadTauW2(angle)

       if ( ASet% StartingDirection(angle) ) then

         !$omp  parallel do default(none) &
         !$omp& shared(Geom,Set,zSetID)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           Set% tPsiM(c)  = Set% tPsi(c)
           Set% tInc(c)   = Set% pInc(c)
         enddo

         !$omp end parallel do

       else

         !$omp  parallel do default(none)  &
         !$omp& shared(Geom, Set, GTA, quadwt, quadTauW1, quadTauW2, zSetID)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           Set% tPsiM(c)  = quadTauW1*Set% tPsi(c) - quadTauW2*Set% tPsiM(c)
           Set% tInc(c)   = quadTauW1*Set% pInc(c) - quadTauW2*Set% tInc(c)
           GTA% PhiInc(c) = GTA% PhiInc(c)         + quadwt*Set% pInc(c)
         enddo

         !$omp end parallel do

       endif

     enddo
   enddo ZoneSetLoop2

   TOMP(end target teams distribute)

   if ( nHyperDomains > 1 ) then

#ifdef TETON_ENABLE_OPENACC
     !$acc parallel loop gang collapse(2) num_gangs(nZoneSets*nGTASets)  &
     !$acc& vector_length(omp_device_team_thread_limit) &
     !$acc& private(Set, ASet, HypPlanePtr, angle, angle0)
#else
     TOMP(target teams distribute collapse(2) num_teams(nZoneSets*nGTASets) &)
     TOMPC(thread_limit(omp_device_team_thread_limit) default(none)&)
     TOMPC(shared(nZoneSets, nGTASets, Quad, angleList, nSets, PsiB)&)
     TOMPC(private(Set, ASet, HypPlanePtr, angle, angle0, cc))
#endif

     ZoneSetLoop3: do zSetID=1,nZoneSets
       do setID=1,nGTASets

         Set         => Quad% SetDataPtr(nSets+setID)
         ASet        => Quad% AngSetPtr(Set% angleSetID)
         angle       =  angleList(1,setID)
         angle0      =  angleList(2,setID)
         HypPlanePtr => ASet% HypPlanePtr(angle)

!        Update Interface values

#ifdef TETON_ENABLE_OPENACC
         !$acc  loop vector &
         !$acc& private(cc)
#else
         !$omp  parallel do default(none) &
         !$omp& shared(Set, HypPlanePtr, angle0, angle, zSetID, PsiB) private(cc)
#endif
         do i=HypPlanePtr% c1(zSetID),HypPlanePtr% c2(zSetID)
           cc = HypPlanePtr% interfaceList(i)
           PsiB(Set%nbelem+i,angle0+angle) = Set% tPsi(cc)
         enddo
#ifndef TETON_ENABLE_OPENACC
         !$omp end parallel do
#endif

       enddo
     enddo ZoneSetLoop3

#ifdef TETON_ENABLE_OPENACC
     !$acc end parallel loop
#else
     TOMP(end target teams distribute)
#endif

   endif

   TOMP_MAP(target exit data map(release: angleList, nSets, nGTASets, nHyperDomains))
   TOMP_MAP(target exit data map(from: PsiB))

   deallocate( angleList )


   return
   end subroutine SweepGreyUCBrzNEW_GPU
