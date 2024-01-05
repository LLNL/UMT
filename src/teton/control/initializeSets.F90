#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   initializeSet   - Constructs a total opacity for transport         *
!                     comprised of absorption, time-absorption and     *
!                     scattering.                                      *
!                                                                      *
!***********************************************************************
 
   subroutine initializeSets 


   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use TimeStepControls_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use GroupSet_mod
   use CommSet_mod
   use ZoneSet_mod
   use Geometry_mod
   use GreyAcceleration_mod
   use Material_mod
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   use ComptonControl_mod
   use flags_mod
#endif
   use RadIntensity_mod
   use MemoryAllocator_mod
   use OMPWrappers_mod
   use Options_mod
   use, intrinsic :: iso_c_binding, only : C_SIZE_T
   use system_info_mod

#if defined(TETON_ENABLE_CUDA)
   use cuda_utils_mod
#endif

#if defined(TETON_ENABLE_CALIPER)
   use caliper_mod
#endif

   implicit none

!  Local

   type(CommSet),  pointer  :: CSet

   integer                  :: setID
   integer                  :: aSetID
   integer                  :: gSetID
   integer                  :: cSetID
   integer                  :: commID
   integer                  :: nSets
   integer                  :: nAngleSets
   integer                  :: nGroupSets
   integer                  :: nGTASets
   integer                  :: nCommSets
   integer                  :: angle
   logical(kind=1)          :: useBoltzmannCompton

   real(adqt)               :: dtrad

!  Constants

   dtrad      = getRadTimeStep(DtControls)
   nSets      = getNumberOfSets(Quad)
   nAngleSets = getNumberOfAngleSets(Quad)
   nGroupSets = getNumberOfGroupSets(Quad)
   nGTASets   = getNumberOfGTASets(Quad)
   nCommSets  = getNumberOfCommSets(Quad)
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   useBoltzmannCompton = getUseBoltzmann(Compton)
#endif

   Size% tau  = one/(speed_light*dtrad)

!  Update the Geometry

   call getGeometry

!  Set Total opacity

   GroupSetLoop: do gSetID=1,nGroupSets
     call setTotalOpacity(gSetID)
   enddo GroupSetLoop

   !$omp parallel do default(none) schedule(static) &
   !$omp& shared(Size,nGTASets,nAngleSets)
   AngleSetLoop: do aSetID=1,nAngleSets+nGTASets
     ! Find reflected angles on all symmetry boundaries
     call findReflectedAngles(aSetID)

     ! Establish angle order for transport sweeps (rtorder) and create
     ! a list of exiting boundary elements by angle (findexit)
     if (Size% ndim >= 2) then
       call rtorder(aSetID)
       call findexit(aSetID)
     endif
   enddo AngleSetLoop
   !$omp end parallel do

!  Create an exiting list for shared boundary elements in 1D

   if (Size% ndim == 1) then
     do cSetID=1,nCommSets
       call findexit1D(cSetID)
     enddo
   endif

!  If we are using the GPU, we need to map some data before the set loop
   if ( Size% useGPU ) then

     ! Use UMPIRE pinned memory allocation size as an estimator for amount of device memory needed.
     ! TODO - The NLsolver is not currently using UMPIRE, this won't be taken
     ! into account on this estimate until that is done.
#if defined(TETON_ENABLE_UMPIRE)
     if ( Allocator%umpire_host_allocator_id >= 0 .AND. Options%isRankVerbose() > 0 ) then
        call printGPUMemRequired(Size%myRankInGroup)
     endif
#endif

!    Map Quadrature List
     TOMP(target enter data map(to:Quad))

     UMPIRE_DEVICE_POOL_ALLOC(Quad%SetDataPtr)
     TOMP(target enter data map(always,to:Quad%SetDataPtr))

     UMPIRE_DEVICE_POOL_ALLOC(Quad%GrpSetPtr)
     TOMP(target enter data map(always,to:Quad%GrpSetPtr))

     UMPIRE_DEVICE_POOL_ALLOC(Quad%AngSetPtr)
     TOMP(target enter data map(always,to:Quad%AngSetPtr))


!    Map Group Sets
     do gSetID=1,nGroupSets
       UMPIRE_DEVICE_POOL_ALLOC(Quad% GrpSetPtr(gSetID)% STotal)
       TOMP(target enter data map(always,to:Quad% GrpSetPtr(gSetID)% STotal))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% GrpSetPtr(gSetID)% Sigt)
       TOMP(target enter data map(always,to:Quad% GrpSetPtr(gSetID)% Sigt))

     enddo

!    Map ZoneSets

     TOMP(target enter data map(to:ZSet))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% AL)
     TOMP(target enter data map(always,to:ZSet% AL))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% AU)
     TOMP(target enter data map(always,to:ZSet% AU))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% nCornerSet)
     TOMP(target enter data map(always,to:ZSet% nCornerSet))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% nCornerBatch)
     TOMP(target enter data map(always,to:ZSet% nCornerBatch))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% offset)
     TOMP(target enter data map(always,to:ZSet% offset))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% cornerList)
     TOMP(target enter data map(always,to:ZSet% cornerList))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% cornerMap)
     TOMP(target enter data map(always,to:ZSet% cornerMap))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% zoneList)
     TOMP(target enter data map(always,to:ZSet% zoneList))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% cornerConverged)
     TOMP(target enter data map(always,to:ZSet% cornerConverged))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% Te)
     TOMP(target enter data map(always,to:ZSet% Te))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% TeOld)
     TOMP(target enter data map(always,to:ZSet% TeOld))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% delta)
     TOMP(target enter data map(always,to:ZSet% delta))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% sumT)
     TOMP(target enter data map(always,to:ZSet% sumT))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% netRate)
     TOMP(target enter data map(always,to:ZSet% netRate))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% dTCompton)
     TOMP(target enter data map(always,to:ZSet% dTCompton))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% B)
     TOMP(target enter data map(always,to:ZSet% B))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% dBdT)
     TOMP(target enter data map(always,to:ZSet% dBdT))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% Snu0)
     TOMP(target enter data map(always,to:ZSet% Snu0))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% dSnu0dT)
     TOMP(target enter data map(always,to:ZSet% dSnu0dT))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% AD)
     TOMP(target enter data map(always,to:ZSet% AD))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% z)
     TOMP(target enter data map(always,to:ZSet% z))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% fk2)
     TOMP(target enter data map(always,to:ZSet% fk2))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% nI)
     TOMP(target enter data map(always,to:ZSet% nI))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% nS)
     TOMP(target enter data map(always,to:ZSet% nS))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% ex)
     TOMP(target enter data map(always,to:ZSet% ex))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% expPH)
     TOMP(target enter data map(always,to:ZSet% expPH))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% comptonDeltaEr)
     TOMP(target enter data map(always,to:ZSet% comptonDeltaEr))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% dComptonDT)
     TOMP(target enter data map(always,to:ZSet% dComptonDT))

     UMPIRE_DEVICE_POOL_ALLOC(ZSet% comptonSe)
     TOMP(target enter data map(always,to:ZSet% comptonSe))


!    Map Angle Sets
     do aSetID=1,nAngleSets+nGTASets

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% nextZ)

       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% nextZ))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% nextC)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% nextC))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% StartingDirection)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% StartingDirection))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% FinishingDirection)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% FinishingDirection))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% Omega)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% Omega))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% Weight)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% Weight))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% numCycles)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% numCycles))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% cycleOffSet)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% cycleOffSet))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% cycleList)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% cycleList))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% nHyperPlanes)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% nHyperPlanes))


       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% HypPlanePtr)

       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% HypPlanePtr))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% BdyExitPtr)
       TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% BdyExitPtr))


       if ( aSetID <= nAngleSets ) then
         UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% AfpNorm)
         TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% AfpNorm))

         UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% AezNorm)
         TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% AezNorm))

         UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% ANormSum)
         TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% ANormSum))

       endif

       do angle=1,Quad% AngSetPtr(aSetID)% numAngles

         ! Unable to map this to UMPIRE device pool, causes a segfault.
         TOMP(target enter data map(to: Quad% AngSetPtr(aSetID)% BdyExitPtr(angle)% bdyList))

         if ( .not. Quad% AngSetPtr(aSetID)% FinishingDirection(angle) ) then
           ! Unable to map these to UMPIRE device pool, causes a segfault or wrong answers.
           TOMP(target enter data map(to:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% zonesInPlane))
           TOMP(target enter data map(to:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane1))
           TOMP(target enter data map(to:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane2))
           TOMP(target enter data map(to:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% ndone))
         endif

       enddo

       if (Size% ndim == 2) then
         UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% angDerivFac)
         TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% angDerivFac))

         UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% quadTauW1)
         TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% quadTauW1))

         UMPIRE_DEVICE_POOL_ALLOC(Quad% AngSetPtr(aSetID)% quadTauW2)
         TOMP(target enter data map(always,to:Quad% AngSetPtr(aSetID)% quadTauW2))

       endif

     enddo

!    Geometry

     TOMP(target enter data map(to:Geom))
     UMPIRE_DEVICE_POOL_ALLOC(Geom% Volume)
     TOMP(target enter data map(always,to:Geom% Volume))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% VolumeOld)
     TOMP(target enter data map(always,to:Geom% VolumeOld))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% VolumeZone)
     TOMP(target enter data map(always,to:Geom% VolumeZone))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% cOffSet)
     TOMP(target enter data map(always,to:Geom% cOffSet))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% numCorner)
     TOMP(target enter data map(always,to:Geom% numCorner))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% CToZone)
     TOMP(target enter data map(always,to:Geom% CToZone))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% corner1)
     TOMP(target enter data map(always,to:Geom% corner1))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% corner2)
     TOMP(target enter data map(always,to:Geom% corner2))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% zone1)
     TOMP(target enter data map(always,to:Geom% zone1))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% zone2)
     TOMP(target enter data map(always,to:Geom% zone2))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% cEZ)
     TOMP(target enter data map(always,to:Geom% cEZ))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% cFP)
     TOMP(target enter data map(always,to:Geom% cFP))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% A_ez)
     TOMP(target enter data map(always,to:Geom% A_ez))

     UMPIRE_DEVICE_POOL_ALLOC(Geom% A_fp)
     TOMP(target enter data map(always,to:Geom% A_fp))


     if (Size% ndim == 2) then
       UMPIRE_DEVICE_POOL_ALLOC(Geom% Area)
       TOMP(target enter data map(always,to:Geom% Area))

       UMPIRE_DEVICE_POOL_ALLOC(Geom% RadiusEZ)
       TOMP(target enter data map(always,to:Geom% RadiusEZ))

       UMPIRE_DEVICE_POOL_ALLOC(Geom% RadiusFP)
       TOMP(target enter data map(always,to:Geom% RadiusFP))

     elseif (Size% ndim == 3) then
       UMPIRE_DEVICE_POOL_ALLOC(Geom% nCFacesArray)
       TOMP(target enter data map(always,to:Geom% nCFacesArray))

     endif

!    Radiation Intensity

     TOMP(target enter data map(to:Rad))
     UMPIRE_DEVICE_POOL_ALLOC(Rad% PhiTotal)
     TOMP(target enter data map(always,to:Rad% PhiTotal))

     UMPIRE_DEVICE_POOL_ALLOC(Rad% radEnergy)
     TOMP(target enter data map(always,to:Rad% radEnergy))


#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     TOMP(target enter data map(to:Compton))

     if (getComptonFlag(Compton) /= comptonType_None) then
       UMPIRE_DEVICE_POOL_ALLOC(Compton% gamMean)
       TOMP(target enter data map(always,to:Compton% gamMean))

       UMPIRE_DEVICE_POOL_ALLOC(Compton% gamSqdDGam)
       TOMP(target enter data map(always,to:Compton% gamSqdDGam))

       UMPIRE_DEVICE_POOL_ALLOC(Compton% gamCubedDGam)
       TOMP(target enter data map(always,to:Compton% gamCubedDGam))

       UMPIRE_DEVICE_POOL_ALLOC(Compton% gamD)
       TOMP(target enter data map(always,to:Compton% gamD))

     endif
#endif

!    GTA

     TOMP(target enter data map(to:GTA))
     UMPIRE_DEVICE_POOL_ALLOC(GTA% GreySource)
     TOMP(target enter data map(always,to:GTA% GreySource))

     UMPIRE_DEVICE_POOL_ALLOC(GTA% GreyCorrection)
     TOMP(target enter data map(always,to:GTA% GreyCorrection))

     UMPIRE_DEVICE_POOL_ALLOC(GTA% Chi)
     TOMP(target enter data map(always,to:GTA% Chi))


     if (Size%useNewGTASolver) then
        UMPIRE_DEVICE_POOL_ALLOC(GTA% TT)
        TOMP(target enter data map(always,to:GTA% TT))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% Pvv)
        TOMP(target enter data map(always,to:GTA% Pvv))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% GreySigTotal)
        TOMP(target enter data map(always,to:GTA% GreySigTotal))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% GreySigScat)
        TOMP(target enter data map(always,to:GTA% GreySigScat))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% GreySigScatVol)
        TOMP(target enter data map(always,to:GTA% GreySigScatVol))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% GreySigtInv)
        TOMP(target enter data map(always,to:GTA% GreySigtInv))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% PhiInc)
        TOMP(target enter data map(always,to:GTA% PhiInc))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% Q)
        TOMP(target enter data map(always,to:GTA% Q))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% TsaSource)
        TOMP(target enter data map(always,to:GTA% TsaSource))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% AfpNorm)
        TOMP(target enter data map(always,to:GTA% AfpNorm))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% AezNorm)
        TOMP(target enter data map(always,to:GTA% AezNorm))

        UMPIRE_DEVICE_POOL_ALLOC(GTA% ANormSum)
        TOMP(target enter data map(always,to:GTA% ANormSum))


        if (Size% ndim == 2) then
          UMPIRE_DEVICE_POOL_ALLOC(GTA% Tvv)
          TOMP(target enter data map(always,to:GTA% Tvv))

        endif
     endif
   endif ! Size%useGPU

!  Initialize communication handles for persistent communicators

!  QUESTION - We're passing in angle set IDs, but inside the initcomm the
!  parameter is 'cSetID'.  Should this be a loop over comm sets or angle sets?? -black27
   do aSetID=1,nAngleSets+nGTASets
     call initcomm(aSetID)
   enddo

!  Begin Initialize Phase

!$omp parallel do schedule(static) default(none) &
!$omp& shared(nSets)
   do setID=1,nSets
!    Initialize Boundary Flux
     call setBoundarySources(setID)
   enddo
!$omp end parallel do

   !Allocate GPU Memory
   if (Size%useGPU) then
     call initializeGPUMemory
   endif

!  Initialize the radiation field (Psi, PsiB, PhiTotal)
   call initPhiTotal
   call initializeRadiationField

!  Map PsiB back to the CPU
   if (Size%useGPU) then
     do cSetID=1,nCommSets
       CSet => getCommSetData(Quad, cSetID)
       do setID=CSet% set1,CSet% set2
         TOMP(target update from(Quad% SetDataPtr(setID)% PsiB))
       enddo
     enddo
   endif

!  Establish angle order for transport sweeps
!$omp  parallel do default(none) schedule(static) &
!$omp& shared(nCommSets, Size)
   do cSetID=1,nCommSets
     call setNetFlux(cSetID)

     if (Size% ndim >= 2) then
       call SweepScheduler(cSetID)
     endif
   enddo
!$omp end parallel do

   if (Size%useGPU) then
     do setID=1,nSets
       UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% AngleOrder)
       TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% AngleOrder))

     enddo
   endif

!  Initialize zonal material properties, 
!  Contains a threaded loop over zones

   call initializeZones

!  Initialize GTA set

   if (Size% ndim >= 2) then
     GTASetLoop: do cSetID=nCommSets+1,nCommSets+nGTASets
       call SweepScheduler(cSetID)
     enddo GTASetLoop
   endif

!    Material

   if ( Size% useGPU ) then
     TOMP(target enter data map(to:Mat))
     UMPIRE_DEVICE_POOL_ALLOC(Mat% Tec)
     TOMP(target enter data map(always,to:Mat% Tec))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% Tecn)
     TOMP(target enter data map(always,to:Mat% Tecn))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% denec)
     TOMP(target enter data map(always,to:Mat% denec))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% cve)
     TOMP(target enter data map(always,to:Mat% cve))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% rho)
     TOMP(target enter data map(always,to:Mat% rho))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% nez)
     TOMP(target enter data map(always,to:Mat% nez))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% stimComptonMult)
     TOMP(target enter data map(always,to:Mat% stimComptonMult))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% Siga)
     TOMP(target enter data map(always,to:Mat% Siga))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% Sigs)
     TOMP(target enter data map(always,to:Mat% Sigs))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% Eta)
     TOMP(target enter data map(always,to:Mat% Eta))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% EmissionRate)
     TOMP(target enter data map(always,to:Mat% EmissionRate))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% SMatEff)
     TOMP(target enter data map(always,to:Mat% SMatEff))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% PowerEmitted)
     TOMP(target enter data map(always,to:Mat% PowerEmitted))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% PowerCompton)
     TOMP(target enter data map(always,to:Mat% PowerCompton))

     UMPIRE_DEVICE_POOL_ALLOC(Mat% nonLinearIterations)
     TOMP(target enter data map(always,to:Mat% nonLinearIterations))

   endif

!  Map GTA set variables

   if ( Size% useGPU ) then

     do setID=nSets+1,nSets+nGTASets
       UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% AngleOrder)
       TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% AngleOrder))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% tPsi)
       TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% tPsi))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% pInc)
       TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% pInc))

       UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% src)
       TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% src))


       if (Size% ndim == 2) then
         UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% tPsiM)
         TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% tPsiM))

         UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% tInc)
         TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% tInc))

       endif
     enddo

   endif

!  Initialize memory for non-linear solver

#if defined(TETON_ENABLE_CUDA)
#  if !defined(TETON_ENABLE_MINIAPP_BUILD)
   if (useBoltzmannCompton .AND. Size% useCUDASolver .AND. Size% ngr >= 16) then
     call fallocateGpuMemory(Size%ngr, Size%nBCITabG2Gs, Size%nBCITabTaus,  &
                             Size% zoneBatchSize, Size%maxCorner)
   endif
#  endif
#endif

#if defined(TETON_ENABLE_UMPIRE)
   if ( Allocator%umpire_host_allocator_id >= 0 .AND. Options%isRankVerbose() > 0 ) then
      call printGPUMemInfo(Size%myRankInGroup)
   endif
#endif

   return
   end subroutine initializeSets 
