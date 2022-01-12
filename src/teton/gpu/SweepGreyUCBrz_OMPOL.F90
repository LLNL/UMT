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

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use SetData_mod
   use AngleSet_mod
   use ArrayChecks_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: sendIndex 
   real(adqt), intent(inout) :: PsiB(Size%nbelem,Size%nangGTA)

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
   integer    :: iter

   integer    :: c
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
   nHyperDomains = Quad% nHyperDomains

!  Quiet the compiler 'variable may not be initialized' warning.
   FinishingDirection = .FALSE.

   allocate( angleList(3,nGTASets) )

   do setID=1,nGTASets
     Set  => Quad% SetDataPtr(nSets+setID)
     ASet => Quad% AngSetPtr(Set% angleSetID)

     Angle              = Set% AngleOrder(sendIndex)
     angleList(1,setID) = Set% AngleOrder(sendIndex)
     angleList(2,setID) = Set% angle0
     angleList(3,setID) = Set% angle0 + Angle + 1 - setID

     FinishingDirection = ASet% FinishingDirection(Angle)
   enddo

   if ( FinishingDirection ) then
     return
   endif

   if (Size%useGPU) then

     ! Verify we won't get out-of-bounds accesses below.
     TETON_CHECK_BOUNDS1(Quad%SetDataPtr, nSets+nGTASets)
     TETON_CHECK_BOUNDS1(Geom%corner1, nZoneSets)
     TETON_CHECK_BOUNDS1(Geom%corner2, nZoneSets)

     TOMP(target data map(to: angleList, nSets, nGTASets, nHyperDomains) map(tofrom: PsiB))

     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(1024) private(Set, setID, zSetID))

     ZoneSetLoop0: do zSetID=1,nZoneSets

       do setID=1,nGTASets
         Set => Quad% SetDataPtr(nSets+setID)

!$omp  parallel do default(none)  &
!$omp& shared(Set, Geom, zSetID) private(c)
         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           Set% tPsi(c) = zero
         enddo
!$omp end parallel do

       enddo

     enddo ZoneSetLoop0

     TOMP(end target teams distribute)

     SweepIteration: do iter=1,GTA% nGreySweepIters

       TOMP(target teams distribute num_teams(nGTASets) thread_limit(1024) private(Set, setID, angle0, Angle))

       GTASetLoop0: do setID=1,nGTASets

         Set    => Quad% SetDataPtr(nSets+setID)
         Angle  =  angleList(1,setID)
         angle0 =  angleList(2,setID)

         !$omp  parallel do default(none) &
         !$omp& shared(Set, PsiB, angle0, Angle) private(b)
         do b=1,Set%nbelem
           Set% tPsi(Set%nCorner+b) = PsiB(b,angle0+Angle)
         enddo
         !$omp end parallel do

       enddo GTASetLoop0

       TOMP(end target teams distribute)


       TOMP(target teams distribute collapse(2) num_teams(nHyperDomains*nGTASets) thread_limit(1024) private(Set, ASet) &)
       TOMPC(private(HypPlanePtr, setID, domID, angle, angle0, angle2, hyperPlane, hplane1, hplane2, ndoneZ, nzones, fac))

       GTASetLoop: do setID=1,nGTASets
         DomainLoop: do domID=1,nHyperDomains

           Set          => Quad% SetDataPtr(nSets+setID)
           ASet         => Quad% AngSetPtr(Set% angleSetID)
           Angle        =  angleList(1,setID)
           angle0       =  angleList(2,setID)
           angle2       =  angleList(3,setID)

           HypPlanePtr  => ASet% HypPlanePtr(Angle)
           hplane1      =  HypPlanePtr% hplane1(domID)
           hplane2      =  HypPlanePtr% hplane2(domID)
           fac          =  ASet% angDerivFac(Angle)
           ndoneZ       =  HypPlanePtr% ndone(domID)

!          Loop through all of the zones using the NEXT list

           HyperPlaneLoop: do hyperPlane=hplane1,hplane2

             nzones = HypPlanePtr% zonesInPlane(hyperPlane)

             !$omp  parallel do default(none)                    &
             !$omp& shared(Geom, GTA, ASet, Set, Quad)           &
             !$omp& shared(ndoneZ, nzones, Angle, angle2, fac)   &
             !$omp& private(ii, zone, nCorner, c, c0, cface)     &
             !$omp& private(i, cez, cfp, aez, afp, R, R_afp)     &
             !$omp& private(sigA, sigA2, gnum, gtau, sez, dInv)  &
             !$omp& schedule(static) 

             ZoneLoop: do ii=1,nzones

               zone    = ASet% nextZ(ndoneZ+ii,Angle)
               nCorner = Geom% numCorner(zone)
               c0      = Geom% cOffSet(zone)

!              Contributions from volume terms and angular derivative

               do c=1,nCorner
                 Set% src(c0+c)  = GTA%TsaSource(c0+c) + &
                                   fac*Geom% Area(c0+c)*Set% tPsiM(c0+c)
                 Set% pInc(c0+c) = fac*Geom% Area(c0+c)*Set% tInc(c0+c)
               enddo

               CornerLoop: do c=1,nCorner

                 sigA = GTA%GreySigTotal(c0+c)*Geom% Area(c0+c)

                 CornerFaceLoop: do cface=1,2

                   afp = Quad% AngSetPtr(angle2)% AfpNorm(cface,c0+c)
                   aez = Quad% AngSetPtr(angle2)% AezNorm(cface,c0+c)

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

!              Solve the corners in the correct order

               do i=1,nCorner

                 c    = ASet% nextC(c0+i,angle) 

                 dInv = one/(Quad% AngSetPtr(angle2)% ANormSum(c0+c) +   &
                        Geom% Volume(c0+c)*GTA%GreySigTotal(c0+c))

!                Corner angular flux
                 Set% tPsi(c0+c) = dInv*Set% src(c0+c)
                 Set% pInc(c0+c) = dInv*Set% pInc(c0+c)

!                Calculate the contribution of this flux to the sources of
!                downstream corners in this zone. 

                 do cface=1,2
                   aez = Quad% AngSetPtr(Angle2)% AezNorm(cface,c0+c)

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

     enddo SweepIteration

!    Update exiting boundary fluxes 

     TOMP(target teams distribute num_teams(nGTASets) thread_limit(1024) private(Set, ASet, BdyExitPtr, setID, angle0, Angle))

     do setID=1,nGTASets

       Set        => Quad% SetDataPtr(nSets+setID)
       ASet       => Quad% AngSetPtr(Set% angleSetID)
       Angle      =  angleList(1,setID)
       angle0     =  angleList(2,setID)
       BdyExitPtr => ASet% BdyExitPtr(Angle)

       !$omp  parallel do default(none) &
       !$omp& shared(Set, BdyExitPtr, angle0, Angle, PsiB) private(i,b,c)

       do i=1,BdyExitPtr% nxBdy
         b = BdyExitPtr% bdyList(1,i)
         c = BdyExitPtr% bdyList(2,i)

         PsiB(b,angle0+Angle) = Set% tPsi(c)
       enddo

       !$omp end parallel do

     enddo
     TOMP(end target teams distribute)


     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(1024) &)
     TOMPC(private(Set, ASet, setID, zSetID, Angle, quadwt, quadTauW1, quadTauW2))

     ZoneSetLoop1: do zSetID=1,nZoneSets

       do setID=1,nGTASets

         Set       => Quad% SetDataPtr(nSets+setID)
         ASet      => Quad% AngSetPtr(Set% angleSetID)

         Angle     =  angleList(1,setID)
         quadwt    =  ASet% weight(Angle)
         quadTauW1 =  ASet% quadTauW1(Angle)
         quadTauW2 =  ASet% quadTauW2(Angle)

         if ( ASet% StartingDirection(Angle) ) then

           !$omp  parallel do default(none) shared(Geom,Set,zSetID) private(c)
           do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
             Set% tPsiM(c)  = Set% tPsi(c)
             Set% tInc(c)   = Set% pInc(c)
           enddo
           !$omp end parallel do

         else

           !$omp  parallel do default(none)  &
           !$omp& shared(Geom, Set, GTA, quadwt, quadTauW1, quadTauW2, zSetID) private(c)
           do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
             Set% tPsiM(c)  = quadTauW1*Set% tPsi(c) - quadTauW2*Set% tPsiM(c)
             Set% tInc(c)   = quadTauW1*Set% pInc(c) - quadTauW2*Set% tInc(c)
             GTA% PhiInc(c) = GTA% PhiInc(c)         + quadwt*Set% pInc(c)
           enddo
           !$omp end parallel do

         endif

       enddo

     enddo ZoneSetLoop1

     TOMP(end target teams distribute)
     TOMP(end target data)
   endif ! if Size%useGPU

   deallocate( angleList )


   return
   end subroutine SweepGreyUCBrzNEW_GPU
