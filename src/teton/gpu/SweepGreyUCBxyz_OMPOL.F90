#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   SweepGreyUCBxyz  - This routine calculates angular fluxes for a    *
!                      single direction and single energy groups for   *
!                      for an upstream corner-balance (UCB) spatial    *
!                      in xyz-geometry. It is only used for Grey       *
!                      Transport Acceleration (GTA).                   *
!                                                                      *
!***********************************************************************

   subroutine SweepGreyUCBxyzNEW_GPU(sendIndex, PsiB)

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use Geometry_mod
   use SetData_mod
   use AngleSet_mod
   use GreyAcceleration_mod
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

   integer    :: setID
   integer    :: zSetID
   integer    :: nSets
   integer    :: nGTASets
   integer    :: nZoneSets
   integer    :: nHyperDomains
   integer    :: iter

   integer    :: angle0
   integer    :: angle
   integer    :: angGTA

   integer    :: c
   integer    :: c0
   integer    :: nCorner
   integer    :: ii
   integer    :: zone
   integer    :: nzones
   integer    :: ndoneZ
   integer    :: hyperPlane
   integer    :: domID

   integer    :: hplane1
   integer    :: hplane2

   integer    :: b
   integer    :: i
   integer    :: cface
   integer    :: ifp
   integer    :: cez
   integer    :: cfp
   integer    :: nCFaces

   real(adqt) :: aez
   real(adqt) :: area_opp
   real(adqt) :: sigv
   real(adqt) :: sigv2
   real(adqt) :: gnum
   real(adqt) :: gtau
   real(adqt) :: sez
   real(adqt) :: psi_opp

   real(adqt) :: denom
   real(adqt) :: afp

   real(adqt), parameter :: fouralpha=1.82_adqt

   real(adqt) :: quadwt

!  Dynamic

   integer, allocatable :: angleList(:,:)

!  Constants

   nSets         = getNumberOfSets(Quad)
   nGTASets      = getNumberOfGTASets(Quad)
   nZoneSets     = getNumberOfZoneSets(Quad)
   nHyperDomains = Quad% nHyperDomains

   allocate( angleList(2,nGTASets) )

   do setID=1,nGTASets
     Set => Quad% SetDataPtr(nSets+setID)
     angleList(1,setID) = Set% AngleOrder(sendIndex)
     angleList(2,setID) = Set% angle0
   enddo

   ! Verify we won't get out-of-bounds accesses below.
   TETON_CHECK_BOUNDS1(Quad%SetDataPtr, nSets+nGTASets)
   TETON_CHECK_BOUNDS1(Geom%corner1, nZoneSets)
   TETON_CHECK_BOUNDS1(Geom%corner2, nZoneSets)

   TOMP(target enter data map(to: angleList, nSets, nGTASets, nHyperDomains, PsiB))
   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, nGTASets, Geom, Quad, nSets)&)
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


   TOMP(target teams distribute num_teams(nGTASets) thread_limit(omp_device_team_thread_limit) default(none)&)
   TOMPC(shared(nGTASets, PsiB, Quad, nSets, angleList)&)
   TOMPC(private(Set, angle0, angle))

   GTASetLoop0: do setID=1,nGTASets

     Set    => Quad% SetDataPtr(nSets+setID)
     angle  =  angleList(1,setID)
     angle0 =  angleList(2,setID)

!    Initialize Boundary Values 

     !$omp  parallel do default(none) &
     !$omp& shared(Set, PsiB, angle0, angle)
     do b=1,Set%nbelem
       Set% tPsi(Set%nCorner+b) = PsiB(b,angle0+angle)
     enddo
     !$omp end parallel do

   enddo GTASetLoop0

   TOMP(end target teams distribute)


   SweepIteration: do iter=1,GTA% nGreySweepIters

     TOMP(target teams distribute collapse(2) num_teams(nHyperDomains*nGTASets) &)
     TOMPC(thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nGTASets, nHyperDomains, Geom, GTA, Quad, nSets, angleList)&)
     TOMPC(private(Set, ASet, HypPlanePtr, angle0, angle, angGTA, hyperPlane, hplane1, hplane2, ndoneZ, nzones, zone, nCorner) &)
     TOMPC(private(c0, ifp, cez, cfp, nCFaces, aez, area_opp, sigv, sigv2, gnum, gtau, sez, psi_opp, denom, afp))

   GTASetLoop: do setID=1,nGTASets
     DomainLoop: do domID=1,nHyperDomains

     Set          => Quad% SetDataPtr(nSets+setID)
     ASet         => Quad% AngSetPtr(Set% angleSetID)
     angle        =  angleList(1,setID)
     angle0       =  angleList(2,setID)
     angGTA       =  angle0 + angle

     HypPlanePtr  => ASet% HypPlanePtr(angle)
     hplane1      =  HypPlanePtr% hplane1(domID)
     hplane2      =  HypPlanePtr% hplane2(domID)
     ndoneZ       =  HypPlanePtr% ndone(domID) 

!  Loop through all of the zones using the NEXT list

     HyperPlaneLoop: do hyperPlane=hplane1,hplane2

       nzones = HypPlanePtr% zonesInPlane(hyperPlane) 

!$omp  parallel do default(none) schedule(static) &
!$omp& shared(Geom, GTA, ASet, Set, ndoneZ, nzones, angle, angGTA) &
!$omp& private(zone, nCorner, c0, ifp, cez, cfp, nCFaces, aez, area_opp, sigv, sigv2, gnum, gtau, sez, psi_opp) &
!$omp& private(denom, afp)

       ZoneLoop: do ii=1,nzones

         zone    = iabs( ASet% nextZ(ndoneZ+ii,angle) )
         nCorner = Geom% numCorner(zone)
         c0      = Geom% cOffSet(zone)

!        Contributions from volume terms

         do c=1,nCorner
           Set% src(c0+c)  = GTA%TsaSource(c0+c)
           Set% pInc(c0+c) = zero
         enddo

         CornerLoop: do c=1,nCorner

           sigv     = Geom% Volume(c0+c)*GTA%GreySigTotal(c0+c)
           nCFaces  = Geom% nCFacesArray(c0+c)

!          Contributions from external corner faces (FP faces)

           do cface=1,nCFaces

             afp = GTA% AfpNorm(cface,c0+c,angGTA) 
             cfp = Geom% cFP(cface,c0+c)

             if ( afp < zero ) then
               Set% src(c0+c)  = Set% src(c0+c)  - afp*Set% tPsi(cfp)
               Set% pInc(c0+c) = Set% pInc(c0+c) - afp*Set% tPsi(cfp)
             endif
           enddo

!          Contributions from interior corner faces (EZ faces)

           do cface=1,ncfaces

             aez = GTA% AezNorm(cface,c0+c,angGTA) 
             cez = Geom% cEZ(cface,c0+c)

             if (aez > zero ) then

               psi_opp  = zero
               area_opp = zero

               ifp = mod(cface,nCFaces) + 1
               afp = GTA% AfpNorm(ifp,c0+c,angGTA)

               if ( afp < zero ) then
                 cfp      =  Geom% cFP(ifp,c0+c)
                 area_opp = -afp
                 psi_opp  = -afp*Set% tPsi(cfp)
               endif

               do i=2,nCFaces-2
                 ifp = mod(ifp,nCFaces) + 1
                 afp = GTA% AfpNorm(ifp,c0+c,angGTA)
                 if ( afp < zero ) then
                   cfp      = Geom% cFP(ifp,c0+c)
                   area_opp = area_opp - afp
                   psi_opp  = psi_opp  - afp*Set% tPsi(cfp)
                 endif
               enddo

               TestOppositeFace: if (area_opp > zero) then

                 psi_opp   = psi_opp/area_opp
                 sigv2     = sigv*sigv

                 gnum      = aez*aez*( fouralpha*sigv2 +    &
                             aez*(four*sigv + three*aez) )

                 gtau      = gnum/    &
                           ( gnum + four*sigv2*sigv2 + aez*sigv*(six*sigv2 + &
                             two*aez*(two*sigv + aez)) )

                 sez       = gtau*sigv*( psi_opp - GTA% Q(c0+c) ) +  &
                             half*aez*(one - gtau)*( GTA% Q(c0+c) - GTA% Q(c0+cez) )

                 Set% src(c0+c)    = Set% src(c0+c)    + sez
                 Set% src(c0+cez)  = Set% src(c0+cez)  - sez

                 Set% pInc(c0+c)   = Set% pInc(c0+c)   + gtau*sigv*psi_opp
                 Set% pInc(c0+cez) = Set% pInc(c0+cez) - gtau*sigv*psi_opp

               else

                 sez               = half*aez*( GTA% Q(c0+c) - GTA% Q(c0+cez) )
                 Set% src(c0+c)    = Set% src(c0+c)   + sez
                 Set% src(c0+cez)  = Set% src(c0+cez) - sez

               endif TestOppositeFace

             endif

           enddo

         enddo CornerLoop

         do i=1,nCorner
           c = ASet% nextC(c0+i,angle) 

!          Corner angular flux
           denom           = GTA% ANormSum(c0+c,angGTA) +  &
                             Geom% Volume(c0+c)*GTA%GreySigTotal(c0+c) 
           Set% tPsi(c0+c) = Set% src(c0+c)/denom
           Set% pInc(c0+c) = Set% pInc(c0+c)/denom

!          Calculate the contribution of this flux to the sources of
!          downstream corners in this zone. The downstream corner index is
!          "ez_exit."

           nCFaces = Geom% nCFacesArray(c0+c)

           do cface=1,nCFaces
             aez = GTA% AezNorm(cface,c0+c,angGTA)

             if (aez > zero) then
               cez               = Geom% cEZ(cface,c0+c)
               Set% src(c0+cez)  = Set% src(c0+cez)  + aez*Set% tPsi(c0+c)
               Set% pInc(c0+cez) = Set% pInc(c0+cez) + aez*Set% pInc(c0+c)
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


!  Update exiting boundary fluxes 

   TOMP(target teams distribute num_teams(nGTASets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nGTASets, nSets, PsiB, angleList, Quad)&)
   TOMPC(private(Set, ASet, BdyExitPtr, angle0, angle, b, c))

   do setID=1,nGTASets

     Set        => Quad% SetDataPtr(nSets+setID)
     ASet       => Quad% AngSetPtr(Set% angleSetID)
     angle      =  angleList(1,setID) 
     angle0     =  angleList(2,setID) 
     BdyExitPtr => ASet% BdyExitPtr(angle)

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


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, nGTASets, nSets, GTA, angleList, Quad, Geom)&)
   TOMPC(private(Set, ASet, zSetID, angle, quadwt)) 

     ZoneSetLoop3: do zSetID=1,nZoneSets

       do setID=1,nGTASets

         Set    => Quad% SetDataPtr(nSets+setID)
         ASet   => Quad% AngSetPtr(Set% angleSetID)

         angle  =  angleList(1,setID) 
         quadwt =  ASet% weight(angle)

!$omp  parallel do default(none)  &
!$omp& shared(Geom, Set, GTA, quadwt, zSetID)
         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           GTA% PhiInc(c) = GTA% PhiInc(c) + quadwt*Set% pInc(c) 
         enddo
!$omp end parallel do

       enddo

     enddo ZoneSetLoop3

   TOMP(end target teams distribute)
   TOMP(target exit data map(release: angleList, nSets, nGTASets, nHyperDomains))
   TOMP(target exit data map(from: PsiB))

   deallocate( angleList )


   return
   end subroutine SweepGreyUCBxyzNEW_GPU

