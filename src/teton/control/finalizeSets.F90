#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   finalizeSets    - Constructs a total opacity for transport         *
!                     comprised of absorption, time-absorption and     *
!                     scattering.                                      *
!                                                                      *
!***********************************************************************
 
   subroutine finalizeSets 

   use kind_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use SetData_mod
   use CommSet_mod
   use AngleSet_mod
   use GroupSet_mod
   use ZoneSet_mod
   use GreyAcceleration_mod
   use Material_mod
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   use ComptonControl_mod
   use flags_mod
#endif
   use RadIntensity_mod
   use OMPWrappers_mod
   use MemoryAllocator_mod
   use Options_mod
   use Size_mod
   use, intrinsic :: iso_c_binding
   use system_info_mod

#if defined(TETON_ENABLE_CUDA)
   use cuda_utils_mod
#endif

   implicit none

!  Local

   type(SetData),  pointer  :: Set
   type(CommSet),  pointer  :: CSet
   type(AngleSet), pointer  :: ASet

   integer                  :: setID
   integer                  :: aSetID
   integer                  :: gSetID
   integer                  :: cSetID
   integer                  :: nSets
   integer                  :: nCommSets
   integer                  :: nAngleSets
   integer                  :: nGroupSets
   integer                  :: nGTASets
   integer                  :: nHyperDomains
   integer                  :: angle
   integer                  :: sweepVersion
   logical(kind=1)          :: useBoltzmannCompton
   logical(kind=1)          :: startCycle

!  Constants

   nSets         = getNumberOfSets(Quad)
   nCommSets     = getNumberOfCommSets(Quad)
   nAngleSets    = getNumberOfAngleSets(Quad)
   nGroupSets    = getNumberOfGroupSets(Quad)
   nGTASets      = getNumberOfGTASets(Quad)
   nHyperDomains = getNumberOfHyperDomains(Quad, 1)
   sweepVersion  = Options% getSweepVersion()

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   useBoltzmannCompton = getUseBoltzmann(Compton)
#else
   useBoltzmannCompton = .FALSE.
#endif

!  Release GPU Memory
   if (Size%useGPU) then

!    Update PhiTotal and edits on the CPU

     TOMP_UPDATE(target update from(Rad% PhiTotal))

     if ( useBoltzmannCompton .and. Size%useCUDASolver .and. Size%ngr >= 16) then
!  In this case these edits are already on the CPU
     else
       TOMP_UPDATE(target update from(Mat% denec))
       TOMP_UPDATE(target update from(Mat% nonLinearIterations))
       TOMP_UPDATE(target update from(Mat% PowerEmitted))
       TOMP_UPDATE(target update from(Mat% PowerCompton))
     endif

     ! Unmap zone sets

     UMPIRE_DEVICE_POOL_FREE(ZSet% nCornerSet)
     TOMP_MAP(target exit data map(always,release:ZSet% nCornerSet))

     UMPIRE_DEVICE_POOL_FREE(ZSet% nCornerBatch)
     TOMP_MAP(target exit data map(always,release:ZSet% nCornerBatch))

     UMPIRE_DEVICE_POOL_FREE(ZSet% offset)
     TOMP_MAP(target exit data map(always,release:ZSet% offset))

     UMPIRE_DEVICE_POOL_FREE(ZSet% cornerList)
     TOMP_MAP(target exit data map(always,release:ZSet% cornerList))

     UMPIRE_DEVICE_POOL_FREE(ZSet% cornerMap)
     TOMP_MAP(target exit data map(always,release:ZSet% cornerMap))

     UMPIRE_DEVICE_POOL_FREE(ZSet% zoneList)
     TOMP_MAP(target exit data map(always,release:ZSet% zoneList))

     UMPIRE_DEVICE_POOL_FREE(ZSet% cornerConverged)
     TOMP_MAP(target exit data map(always,release:ZSet% cornerConverged))

     UMPIRE_DEVICE_POOL_FREE(ZSet% Te)
     TOMP_MAP(target exit data map(always,release:ZSet% Te))

     UMPIRE_DEVICE_POOL_FREE(ZSet% TeOld)
     TOMP_MAP(target exit data map(always,release:ZSet% TeOld))

     UMPIRE_DEVICE_POOL_FREE(ZSet% delta)
     TOMP_MAP(target exit data map(always,release:ZSet% delta))

     UMPIRE_DEVICE_POOL_FREE(ZSet% sumT)
     TOMP_MAP(target exit data map(always,release:ZSet% sumT))

     UMPIRE_DEVICE_POOL_FREE(ZSet% netRate)
     TOMP_MAP(target exit data map(always,release:ZSet% netRate))

     UMPIRE_DEVICE_POOL_FREE(ZSet% dTCompton)
     TOMP_MAP(target exit data map(always,release:ZSet% dTCompton))

     UMPIRE_DEVICE_POOL_FREE(ZSet% B)
     TOMP_MAP(target exit data map(always,release:ZSet% B))

     UMPIRE_DEVICE_POOL_FREE(ZSet% dBdT)
     TOMP_MAP(target exit data map(always,release:ZSet% dBdT))

     UMPIRE_DEVICE_POOL_FREE(ZSet% Snu0)
     TOMP_MAP(target exit data map(always,release:ZSet% Snu0))

     UMPIRE_DEVICE_POOL_FREE(ZSet% dSnu0dT)
     TOMP_MAP(target exit data map(always,release:ZSet% dSnu0dT))

     UMPIRE_DEVICE_POOL_FREE(ZSet% AD)
     TOMP_MAP(target exit data map(always,release:ZSet% AD))

     UMPIRE_DEVICE_POOL_FREE(ZSet% z)
     TOMP_MAP(target exit data map(always,release:ZSet% z))

     UMPIRE_DEVICE_POOL_FREE(ZSet% fk2)
     TOMP_MAP(target exit data map(always,release:ZSet% fk2))

     UMPIRE_DEVICE_POOL_FREE(ZSet% nI)
     TOMP_MAP(target exit data map(always,release:ZSet% nI))

     UMPIRE_DEVICE_POOL_FREE(ZSet% nS)
     TOMP_MAP(target exit data map(always,release:ZSet% nS))

     UMPIRE_DEVICE_POOL_FREE(ZSet% ex)
     TOMP_MAP(target exit data map(always,release:ZSet% ex))

     UMPIRE_DEVICE_POOL_FREE(ZSet% expPH)
     TOMP_MAP(target exit data map(always,release:ZSet% expPH))

     UMPIRE_DEVICE_POOL_FREE(ZSet% comptonDeltaEr)
     TOMP_MAP(target exit data map(always,release:ZSet% comptonDeltaEr))

     UMPIRE_DEVICE_POOL_FREE(ZSet% dComptonDT)
     TOMP_MAP(target exit data map(always,release:ZSet% dComptonDT))

     UMPIRE_DEVICE_POOL_FREE(ZSet% comptonSe)
     TOMP_MAP(target exit data map(always,release:ZSet% comptonSe))

     UMPIRE_DEVICE_POOL_FREE(ZSet% AU)
     TOMP_MAP(target exit data map(always,release:ZSet% AU))

     UMPIRE_DEVICE_POOL_FREE(ZSet% AL)
     TOMP_MAP(target exit data map(always,release:ZSet% AL))

     TOMP_MAP(target exit data map(release: ZSet))

     ! Unmap group sets
     do gSetID=1,nGroupSets
       UMPIRE_DEVICE_POOL_FREE(Quad% GrpSetPtr(gSetID)% STotal)
       TOMP_MAP(target exit data map(always,release:Quad% GrpSetPtr(gSetID)% STotal))

       UMPIRE_DEVICE_POOL_FREE(Quad% GrpSetPtr(gSetID)% Sigt)
       TOMP_MAP(target exit data map(always,release:Quad% GrpSetPtr(gSetID)% Sigt))

     enddo

     do aSetID=1,nAngleSets+nGTASets

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% nextZ)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% nextZ))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% nextC)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% nextC))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% StartingDirection)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% StartingDirection))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% FinishingDirection)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% FinishingDirection))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% Omega)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% Omega))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% Weight)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% Weight))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% numCycles)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% numCycles))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% cycleOffSet)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% cycleOffSet))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% cycleList)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% cycleList))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% nHyperPlanes)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% nHyperPlanes))


       ! This loop unmaps internal components of HypPlanePtr and BdyExitPtr.
       ! Delay unmapping these until this loop is done.
       do angle=1,Quad%AngSetPtr(aSetID)% numAngles
         ! Unable to map this to UMPIRE device pool, causes segfault.
         TOMP_MAP(target exit data map(release:Quad% AngSetPtr(aSetID)% BdyExitPtr(angle)%bdyList))

         if ( .not. Quad%AngSetPtr(aSetID)% FinishingDirection(angle) ) then
           ! Unable to map these to UMPIRE device pool, causes segfault or wrong answers.

           if (aSetID > nAngleSets) then
             TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% zonesInPlane))
           else
             if ( sweepVersion == 0 ) then
               TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% zonesInPlane))
             else
               TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% cornersInPlane))
             endif
           endif

           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane1))
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane2))
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% ndone))
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% interfaceList))
         endif
       enddo

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% HypPlanePtr)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% HypPlanePtr))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% BdyExitPtr)
       TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% BdyExitPtr))


       if ( aSetID <= nAngleSets ) then
         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% AfpNorm)
         TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% AfpNorm))

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% AezNorm)
         TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% AezNorm))

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% ANormSum)
         TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% ANormSum))

       endif


       if (Size% ndim == 2) then

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% angDerivFac)
         TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% angDerivFac))

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% quadTauW1)
         TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% quadTauW1))

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% quadTauW2)
         TOMP_MAP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% quadTauW2))

       endif

     enddo

!    Geometry

     UMPIRE_DEVICE_POOL_FREE(Geom% Volume)
     TOMP_MAP(target exit data map(always,release:Geom% Volume))

     UMPIRE_DEVICE_POOL_FREE(Geom% VolumeOld)
     TOMP_MAP(target exit data map(always,release:Geom% VolumeOld))

     UMPIRE_DEVICE_POOL_FREE(Geom% VolumeZone)
     TOMP_MAP(target exit data map(always,release:Geom% VolumeZone))

     UMPIRE_DEVICE_POOL_FREE(Geom% cOffSet)
     TOMP_MAP(target exit data map(always,release:Geom% cOffSet))

     UMPIRE_DEVICE_POOL_FREE(Geom% numCorner)
     TOMP_MAP(target exit data map(always,release:Geom% numCorner))

     UMPIRE_DEVICE_POOL_FREE(Geom% CToZone)
     TOMP_MAP(target exit data map(always,release:Geom% CToZone))

     UMPIRE_DEVICE_POOL_FREE(Geom% corner1)
     TOMP_MAP(target exit data map(always,release:Geom% corner1))

     UMPIRE_DEVICE_POOL_FREE(Geom% corner2)
     TOMP_MAP(target exit data map(always,release:Geom% corner2))

     UMPIRE_DEVICE_POOL_FREE(Geom% zone1)
     TOMP_MAP(target exit data map(always,release:Geom% zone1))

     UMPIRE_DEVICE_POOL_FREE(Geom% zone2)
     TOMP_MAP(target exit data map(always,release:Geom% zone2))

     UMPIRE_DEVICE_POOL_FREE(Geom% cEZ)
     TOMP_MAP(target exit data map(always,release:Geom% cEZ))

     UMPIRE_DEVICE_POOL_FREE(Geom% cFP)
     TOMP_MAP(target exit data map(always,release:Geom% cFP))

     UMPIRE_DEVICE_POOL_FREE(Geom% A_ez)
     TOMP_MAP(target exit data map(always,release:Geom% A_ez))

     UMPIRE_DEVICE_POOL_FREE(Geom% A_fp)
     TOMP_MAP(target exit data map(always,release:Geom% A_fp))


     if (Size% ndim == 2) then
       UMPIRE_DEVICE_POOL_FREE(Geom% Area)
       TOMP_MAP(target exit data map(always,release:Geom% Area))

       UMPIRE_DEVICE_POOL_FREE(Geom% RadiusEZ)
       TOMP_MAP(target exit data map(always,release:Geom% RadiusEZ))

       UMPIRE_DEVICE_POOL_FREE(Geom% RadiusFP)
       TOMP_MAP(target exit data map(always,release:Geom% RadiusFP))

     elseif (Size% ndim == 3) then
       UMPIRE_DEVICE_POOL_FREE(Geom% nCFacesArray)
       TOMP_MAP(target exit data map(always,release:Geom% nCFacesArray))

     endif

     TOMP_MAP(target exit data map(release:Geom))

!    Radiation Intensity

     UMPIRE_DEVICE_POOL_FREE(Rad% PhiTotal)
     TOMP_MAP(target exit data map(always,release:Rad% PhiTotal))

     UMPIRE_DEVICE_POOL_FREE(Rad% radEnergy)
     TOMP_MAP(target exit data map(always,release:Rad% radEnergy))

     TOMP_MAP(target exit data map(release:Rad))

!    GTA

     if (Size%useNewGTASolver) then
       UMPIRE_DEVICE_POOL_FREE(GTA% TT)
       TOMP_MAP(target exit data map(always,release:GTA% TT))

       UMPIRE_DEVICE_POOL_FREE(GTA% Pvv)
       TOMP_MAP(target exit data map(always,release:GTA% Pvv))

       UMPIRE_DEVICE_POOL_FREE(GTA% GreySigTotal)
       TOMP_MAP(target exit data map(always,release:GTA% GreySigTotal))

       UMPIRE_DEVICE_POOL_FREE(GTA% GreySigScat)
       TOMP_MAP(target exit data map(always,release:GTA% GreySigScat))

       UMPIRE_DEVICE_POOL_FREE(GTA% GreySigScatVol)
       TOMP_MAP(target exit data map(always,release:GTA% GreySigScatVol))

       UMPIRE_DEVICE_POOL_FREE(GTA% GreySigtInv)
       TOMP_MAP(target exit data map(always,release:GTA% GreySigtInv))

       UMPIRE_DEVICE_POOL_FREE(GTA% PhiInc)
       TOMP_MAP(target exit data map(always,release:GTA% PhiInc))

       UMPIRE_DEVICE_POOL_FREE(GTA% Sscat)
       TOMP_MAP(target exit data map(always,release:GTA% Sscat))

       UMPIRE_DEVICE_POOL_FREE(GTA% Q)
       TOMP_MAP(target exit data map(always,release:GTA% Q))

       UMPIRE_DEVICE_POOL_FREE(GTA% TsaSource)
       TOMP_MAP(target exit data map(always,release:GTA% TsaSource))

       UMPIRE_DEVICE_POOL_FREE(GTA% AfpNorm)
       TOMP_MAP(target exit data map(always,release:GTA% AfpNorm))

       UMPIRE_DEVICE_POOL_FREE(GTA% AezNorm)
       TOMP_MAP(target exit data map(always,release:GTA% AezNorm))

       UMPIRE_DEVICE_POOL_FREE(GTA% ANormSum)
       TOMP_MAP(target exit data map(always,release:GTA% ANormSum))


       if (Size% ndim == 2) then
         UMPIRE_DEVICE_POOL_FREE(GTA% Tvv)
         TOMP_MAP(target exit data map(always,release:GTA% Tvv))

       endif
     endif

     UMPIRE_DEVICE_POOL_FREE(GTA% GreySource)
     TOMP_MAP(target exit data map(always,release:GTA% GreySource))

     UMPIRE_DEVICE_POOL_FREE(GTA% GreyCorrection)
     TOMP_MAP(target exit data map(always,release:GTA% GreyCorrection))

     UMPIRE_DEVICE_POOL_FREE(GTA% Chi)
     TOMP_MAP(target exit data map(always,release:GTA% Chi))

     TOMP_MAP(target exit data map(release:GTA))

     do setID=nSets+1,nSets+nGTASets
       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% AngleOrder)
       TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% AngleOrder))

       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% tPsi)
       TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% tPsi))

       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% pInc)
       TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% pInc))

       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% src)
       TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% src))


       if (Size% ndim == 2) then
         UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% tPsiM)
         TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% tPsiM))

         UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% tInc)
         TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% tInc))

       endif
     enddo

!    Material

     UMPIRE_DEVICE_POOL_FREE(Mat% Tec)
     TOMP_MAP(target exit data map(always,release:Mat% Tec))

     UMPIRE_DEVICE_POOL_FREE(Mat% Tecn)
     TOMP_MAP(target exit data map(always,release:Mat% Tecn))

     UMPIRE_DEVICE_POOL_FREE(Mat% denec)
     TOMP_MAP(target exit data map(always,release:Mat% denec))

     UMPIRE_DEVICE_POOL_FREE(Mat% cve)
     TOMP_MAP(target exit data map(always,release:Mat% cve))

     UMPIRE_DEVICE_POOL_FREE(Mat% rho)
     TOMP_MAP(target exit data map(always,release:Mat% rho))

     UMPIRE_DEVICE_POOL_FREE(Mat% nez)
     TOMP_MAP(target exit data map(always,release:Mat% nez))

     UMPIRE_DEVICE_POOL_FREE(Mat% stimComptonMult)
     TOMP_MAP(target exit data map(always,release:Mat% stimComptonMult))

     UMPIRE_DEVICE_POOL_FREE(Mat% Siga)
     TOMP_MAP(target exit data map(always,release:Mat% Siga))

     UMPIRE_DEVICE_POOL_FREE(Mat% Sigs)
     TOMP_MAP(target exit data map(always,release:Mat% Sigs))

     UMPIRE_DEVICE_POOL_FREE(Mat% Eta)
     TOMP_MAP(target exit data map(always,release:Mat% Eta))

     UMPIRE_DEVICE_POOL_FREE(Mat% EmissionRate)
     TOMP_MAP(target exit data map(always,release:Mat% EmissionRate))

     UMPIRE_DEVICE_POOL_FREE(Mat% SMatEff)
     TOMP_MAP(target exit data map(always,release:Mat% SMatEff))

     UMPIRE_DEVICE_POOL_FREE(Mat% PowerEmitted)
     TOMP_MAP(target exit data map(always,release:Mat% PowerEmitted))

     UMPIRE_DEVICE_POOL_FREE(Mat% PowerCompton)
     TOMP_MAP(target exit data map(always,release:Mat% PowerCompton))

     UMPIRE_DEVICE_POOL_FREE(Mat% nonLinearIterations)
     TOMP_MAP(target exit data map(always,release:Mat% nonLinearIterations))

     TOMP_MAP(target exit data map(release:Mat))

! IF THESE ARE UNALLOCATED WILL CRASH?
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     if (getComptonFlag(Compton) /= comptonType_None) then
       UMPIRE_DEVICE_POOL_FREE(Compton% gamMean)
       TOMP_MAP(target exit data map(always,release:Compton% gamMean))

       UMPIRE_DEVICE_POOL_FREE(Compton% gamSqdDGam)
       TOMP_MAP(target exit data map(always,release:Compton% gamSqdDGam))

       UMPIRE_DEVICE_POOL_FREE(Compton% gamCubedDGam)
       TOMP_MAP(target exit data map(always,release:Compton% gamCubedDGam))

       UMPIRE_DEVICE_POOL_FREE(Compton% gamD)
       TOMP_MAP(target exit data map(always,release:Compton% gamD))

     endif

     TOMP_MAP(target exit data map(release:Compton))
#endif

   endif !endif useGPU

!  Deallocation for Communication Sets 
   CommSetLoop: do cSetID=1,nCommSets+nGTASets
     CSet  => getCommSetData(Quad, cSetID)
     call destructComm(CSet)
   enddo CommSetLoop

   SetLoop: do setID=1,nSets

     Set => getSetData(Quad, setID)

!    Update Psi on the host and Release GPU Memory

     if ( Size% useGPU ) then
       call finalizeGPUMemory(setID)
       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% AngleOrder)
       TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% AngleOrder))
     endif

!  Release Dynamic Memory allocated at the beginning of the time step

     if (Size% ndim > 1) then
       call Set%destructDynMemory(nHyperDomains)
     endif

   enddo SetLoop

!  Update boundary edits

!$omp parallel do default(none) schedule(static) &
!$omp& shared(nSets)
   do setID=1,nSets
     call BoundaryEdit(setID)
   enddo
!$omp end parallel do

!  Update end-of-cycle material properties used by the host code

   startCycle = .FALSE.
   call advanceMaterialProperties(startCycle)

!  Release set pointers

   if ( Size% useGPU ) then

     UMPIRE_DEVICE_POOL_FREE(Quad%AngSetPtr)

     TOMP_MAP(target exit data map(always,release:Quad%AngSetPtr))

     UMPIRE_DEVICE_POOL_FREE(Quad%GrpSetPtr)
     TOMP_MAP(target exit data map(always,release:Quad%GrpSetPtr))

     UMPIRE_DEVICE_POOL_FREE(Quad%SetDataPtr)
     TOMP_MAP(target exit data map(always,release:Quad%SetDataPtr))

     TOMP_MAP(target exit data map(release:Quad))

   endif

!  Update radiation energy density

   call setEnergyDensity

!  Release Communication buffers and hyperplanes

   if (Size% ndim >= 2) then

     AngleSetLoop: do setID=1,nAngleSets+nGTASets
       ASet => getAngleSetData(Quad, setID)

       if (setID > nAngleSets) then
         sweepVersion = 0
       endif

       call destructHyperPlane(ASet, sweepVersion)
       call destructBdyExitList(ASet)
       call destructCycleList(ASet)
     enddo AngleSetLoop

   endif

#if defined(TETON_ENABLE_CUDA)
#  if !defined(TETON_ENABLE_MINIAPP_BUILD)
   if (useBoltzmannCompton .AND. Size% useCUDASolver .and. Size%ngr >= 16) then
     call freeGpuMemory ()
   endif
#  endif
#endif

   return
   end subroutine finalizeSets 

