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
   integer                  :: commID
   integer                  :: nSets
   integer                  :: nCommSets
   integer                  :: nAngleSets
   integer                  :: nGroupSets
   integer                  :: nGTASets
   integer                  :: angle
   logical(kind=1)          :: useBoltzmannCompton
   logical(kind=1)          :: startCycle

!  Constants

   nSets      = getNumberOfSets(Quad)
   nCommSets  = getNumberOfCommSets(Quad)
   nAngleSets = getNumberOfAngleSets(Quad)
   nGroupSets = getNumberOfGroupSets(Quad)
   nGTASets   = getNumberOfGTASets(Quad)
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   useBoltzmannCompton = getUseBoltzmann(Compton)
#else
   useBoltzmannCompton = .FALSE.
#endif

!  Release GPU Memory
   if (Size%useGPU) then

!    Update PhiTotal and edits on the CPU

     TOMP(target update from(Rad% PhiTotal))

     if ( useBoltzmannCompton .and. Size%useCUDASolver .and. Size%ngr >= 16) then
!  In this case these edits are already on the CPU
     else
       TOMP(target update from(Mat% denec))
       TOMP(target update from(Mat% nonLinearIterations))
       TOMP(target update from(Mat% PowerEmitted))
       TOMP(target update from(Mat% PowerCompton))
     endif

     ! Unmap zone sets

     UMPIRE_DEVICE_POOL_FREE(ZSet% nCornerSet)
     TOMP(target exit data map(always,release:ZSet% nCornerSet))

     UMPIRE_DEVICE_POOL_FREE(ZSet% nCornerBatch)
     TOMP(target exit data map(always,release:ZSet% nCornerBatch))

     UMPIRE_DEVICE_POOL_FREE(ZSet% offset)
     TOMP(target exit data map(always,release:ZSet% offset))

     UMPIRE_DEVICE_POOL_FREE(ZSet% cornerList)
     TOMP(target exit data map(always,release:ZSet% cornerList))

     UMPIRE_DEVICE_POOL_FREE(ZSet% cornerMap)
     TOMP(target exit data map(always,release:ZSet% cornerMap))

     UMPIRE_DEVICE_POOL_FREE(ZSet% zoneList)
     TOMP(target exit data map(always,release:ZSet% zoneList))

     UMPIRE_DEVICE_POOL_FREE(ZSet% cornerConverged)
     TOMP(target exit data map(always,release:ZSet% cornerConverged))

     UMPIRE_DEVICE_POOL_FREE(ZSet% Te)
     TOMP(target exit data map(always,release:ZSet% Te))

     UMPIRE_DEVICE_POOL_FREE(ZSet% TeOld)
     TOMP(target exit data map(always,release:ZSet% TeOld))

     UMPIRE_DEVICE_POOL_FREE(ZSet% delta)
     TOMP(target exit data map(always,release:ZSet% delta))

     UMPIRE_DEVICE_POOL_FREE(ZSet% sumT)
     TOMP(target exit data map(always,release:ZSet% sumT))

     UMPIRE_DEVICE_POOL_FREE(ZSet% netRate)
     TOMP(target exit data map(always,release:ZSet% netRate))

     UMPIRE_DEVICE_POOL_FREE(ZSet% dTCompton)
     TOMP(target exit data map(always,release:ZSet% dTCompton))

     UMPIRE_DEVICE_POOL_FREE(ZSet% B)
     TOMP(target exit data map(always,release:ZSet% B))

     UMPIRE_DEVICE_POOL_FREE(ZSet% dBdT)
     TOMP(target exit data map(always,release:ZSet% dBdT))

     UMPIRE_DEVICE_POOL_FREE(ZSet% Snu0)
     TOMP(target exit data map(always,release:ZSet% Snu0))

     UMPIRE_DEVICE_POOL_FREE(ZSet% dSnu0dT)
     TOMP(target exit data map(always,release:ZSet% dSnu0dT))

     UMPIRE_DEVICE_POOL_FREE(ZSet% AD)
     TOMP(target exit data map(always,release:ZSet% AD))

     UMPIRE_DEVICE_POOL_FREE(ZSet% z)
     TOMP(target exit data map(always,release:ZSet% z))

     UMPIRE_DEVICE_POOL_FREE(ZSet% fk2)
     TOMP(target exit data map(always,release:ZSet% fk2))

     UMPIRE_DEVICE_POOL_FREE(ZSet% nI)
     TOMP(target exit data map(always,release:ZSet% nI))

     UMPIRE_DEVICE_POOL_FREE(ZSet% nS)
     TOMP(target exit data map(always,release:ZSet% nS))

     UMPIRE_DEVICE_POOL_FREE(ZSet% ex)
     TOMP(target exit data map(always,release:ZSet% ex))

     UMPIRE_DEVICE_POOL_FREE(ZSet% expPH)
     TOMP(target exit data map(always,release:ZSet% expPH))

     UMPIRE_DEVICE_POOL_FREE(ZSet% comptonDeltaEr)
     TOMP(target exit data map(always,release:ZSet% comptonDeltaEr))

     UMPIRE_DEVICE_POOL_FREE(ZSet% dComptonDT)
     TOMP(target exit data map(always,release:ZSet% dComptonDT))

     UMPIRE_DEVICE_POOL_FREE(ZSet% comptonSe)
     TOMP(target exit data map(always,release:ZSet% comptonSe))

     UMPIRE_DEVICE_POOL_FREE(ZSet% AU)
     TOMP(target exit data map(always,release:ZSet% AU))

     UMPIRE_DEVICE_POOL_FREE(ZSet% AL)
     TOMP(target exit data map(always,release:ZSet% AL))

     TOMP(target exit data map(release: ZSet))

     ! Unmap group sets
     do gSetID=1,nGroupSets
       UMPIRE_DEVICE_POOL_FREE(Quad% GrpSetPtr(gSetID)% STotal)
       TOMP(target exit data map(always,release:Quad% GrpSetPtr(gSetID)% STotal))

       UMPIRE_DEVICE_POOL_FREE(Quad% GrpSetPtr(gSetID)% Sigt)
       TOMP(target exit data map(always,release:Quad% GrpSetPtr(gSetID)% Sigt))

     enddo

     do aSetID=1,nAngleSets+nGTASets

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% nextZ)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% nextZ))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% nextC)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% nextC))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% StartingDirection)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% StartingDirection))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% FinishingDirection)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% FinishingDirection))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% Omega)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% Omega))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% Weight)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% Weight))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% numCycles)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% numCycles))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% cycleOffSet)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% cycleOffSet))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% cycleList)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% cycleList))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% nHyperPlanes)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% nHyperPlanes))


       ! This loop unmaps internal components of HypPlanePtr and BdyExitPtr.
       ! Delay unmapping these until this loop is done.
       do angle=1,Quad%AngSetPtr(aSetID)% numAngles
         ! Unable to map this to UMPIRE device pool, causes segfault.
         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% BdyExitPtr(angle)%bdyList))

         if ( .not. Quad%AngSetPtr(aSetID)% FinishingDirection(angle) ) then
           ! Unable to map these to UMPIRE device pool, causes segfault or wrong answers.
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% zonesInPlane))
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane1))
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane2))
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% ndone))
         endif
       enddo

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% HypPlanePtr)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% HypPlanePtr))

       UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% BdyExitPtr)
       TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% BdyExitPtr))


       if ( aSetID <= nAngleSets ) then
         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% AfpNorm)
         TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% AfpNorm))

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% AezNorm)
         TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% AezNorm))

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% ANormSum)
         TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% ANormSum))

       endif


       if (Size% ndim == 2) then

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% angDerivFac)
         TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% angDerivFac))

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% quadTauW1)
         TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% quadTauW1))

         UMPIRE_DEVICE_POOL_FREE(Quad% AngSetPtr(aSetID)% quadTauW2)
         TOMP(target exit data map(always,release:Quad% AngSetPtr(aSetID)% quadTauW2))

       endif

     enddo

!    Geometry

     UMPIRE_DEVICE_POOL_FREE(Geom% Volume)
     TOMP(target exit data map(always,release:Geom% Volume))

     UMPIRE_DEVICE_POOL_FREE(Geom% VolumeOld)
     TOMP(target exit data map(always,release:Geom% VolumeOld))

     UMPIRE_DEVICE_POOL_FREE(Geom% VolumeZone)
     TOMP(target exit data map(always,release:Geom% VolumeZone))

     UMPIRE_DEVICE_POOL_FREE(Geom% cOffSet)
     TOMP(target exit data map(always,release:Geom% cOffSet))

     UMPIRE_DEVICE_POOL_FREE(Geom% numCorner)
     TOMP(target exit data map(always,release:Geom% numCorner))

     UMPIRE_DEVICE_POOL_FREE(Geom% CToZone)
     TOMP(target exit data map(always,release:Geom% CToZone))

     UMPIRE_DEVICE_POOL_FREE(Geom% corner1)
     TOMP(target exit data map(always,release:Geom% corner1))

     UMPIRE_DEVICE_POOL_FREE(Geom% corner2)
     TOMP(target exit data map(always,release:Geom% corner2))

     UMPIRE_DEVICE_POOL_FREE(Geom% zone1)
     TOMP(target exit data map(always,release:Geom% zone1))

     UMPIRE_DEVICE_POOL_FREE(Geom% zone2)
     TOMP(target exit data map(always,release:Geom% zone2))

     UMPIRE_DEVICE_POOL_FREE(Geom% cEZ)
     TOMP(target exit data map(always,release:Geom% cEZ))

     UMPIRE_DEVICE_POOL_FREE(Geom% cFP)
     TOMP(target exit data map(always,release:Geom% cFP))

     UMPIRE_DEVICE_POOL_FREE(Geom% A_ez)
     TOMP(target exit data map(always,release:Geom% A_ez))

     UMPIRE_DEVICE_POOL_FREE(Geom% A_fp)
     TOMP(target exit data map(always,release:Geom% A_fp))


     if (Size% ndim == 2) then
       UMPIRE_DEVICE_POOL_FREE(Geom% Area)
       TOMP(target exit data map(always,release:Geom% Area))

       UMPIRE_DEVICE_POOL_FREE(Geom% RadiusEZ)
       TOMP(target exit data map(always,release:Geom% RadiusEZ))

       UMPIRE_DEVICE_POOL_FREE(Geom% RadiusFP)
       TOMP(target exit data map(always,release:Geom% RadiusFP))

     elseif (Size% ndim == 3) then
       UMPIRE_DEVICE_POOL_FREE(Geom% nCFacesArray)
       TOMP(target exit data map(always,release:Geom% nCFacesArray))

     endif

     TOMP(target exit data map(release:Geom))

!    Radiation Intensity

     UMPIRE_DEVICE_POOL_FREE(Rad% PhiTotal)
     TOMP(target exit data map(always,release:Rad% PhiTotal))

     UMPIRE_DEVICE_POOL_FREE(Rad% radEnergy)
     TOMP(target exit data map(always,release:Rad% radEnergy))

     TOMP(target exit data map(release:Rad))

!    GTA

     if (Size%useNewGTASolver) then
       UMPIRE_DEVICE_POOL_FREE(GTA% TT)
       TOMP(target exit data map(always,release:GTA% TT))

       UMPIRE_DEVICE_POOL_FREE(GTA% Pvv)
       TOMP(target exit data map(always,release:GTA% Pvv))

       UMPIRE_DEVICE_POOL_FREE(GTA% GreySigTotal)
       TOMP(target exit data map(always,release:GTA% GreySigTotal))

       UMPIRE_DEVICE_POOL_FREE(GTA% GreySigScat)
       TOMP(target exit data map(always,release:GTA% GreySigScat))

       UMPIRE_DEVICE_POOL_FREE(GTA% GreySigScatVol)
       TOMP(target exit data map(always,release:GTA% GreySigScatVol))

       UMPIRE_DEVICE_POOL_FREE(GTA% GreySigtInv)
       TOMP(target exit data map(always,release:GTA% GreySigtInv))

       UMPIRE_DEVICE_POOL_FREE(GTA% PhiInc)
       TOMP(target exit data map(always,release:GTA% PhiInc))

       UMPIRE_DEVICE_POOL_FREE(GTA% Q)
       TOMP(target exit data map(always,release:GTA% Q))

       UMPIRE_DEVICE_POOL_FREE(GTA% TsaSource)
       TOMP(target exit data map(always,release:GTA% TsaSource))

       UMPIRE_DEVICE_POOL_FREE(GTA% AfpNorm)
       TOMP(target exit data map(always,release:GTA% AfpNorm))

       UMPIRE_DEVICE_POOL_FREE(GTA% AezNorm)
       TOMP(target exit data map(always,release:GTA% AezNorm))

       UMPIRE_DEVICE_POOL_FREE(GTA% ANormSum)
       TOMP(target exit data map(always,release:GTA% ANormSum))


       if (Size% ndim == 2) then
         UMPIRE_DEVICE_POOL_FREE(GTA% Tvv)
         TOMP(target exit data map(always,release:GTA% Tvv))

       endif
     endif

     UMPIRE_DEVICE_POOL_FREE(GTA% GreySource)
     TOMP(target exit data map(always,release:GTA% GreySource))

     UMPIRE_DEVICE_POOL_FREE(GTA% GreyCorrection)
     TOMP(target exit data map(always,release:GTA% GreyCorrection))

     UMPIRE_DEVICE_POOL_FREE(GTA% Chi)
     TOMP(target exit data map(always,release:GTA% Chi))

     TOMP(target exit data map(release:GTA))

     do setID=nSets+1,nSets+nGTASets
       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% AngleOrder)
       TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% AngleOrder))

       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% tPsi)
       TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% tPsi))

       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% pInc)
       TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% pInc))

       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% src)
       TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% src))


       if (Size% ndim == 2) then
         UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% tPsiM)
         TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% tPsiM))

         UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% tInc)
         TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% tInc))

       endif
     enddo

!    Material

     UMPIRE_DEVICE_POOL_FREE(Mat% Tec)
     TOMP(target exit data map(always,release:Mat% Tec))

     UMPIRE_DEVICE_POOL_FREE(Mat% Tecn)
     TOMP(target exit data map(always,release:Mat% Tecn))

     UMPIRE_DEVICE_POOL_FREE(Mat% denec)
     TOMP(target exit data map(always,release:Mat% denec))

     UMPIRE_DEVICE_POOL_FREE(Mat% cve)
     TOMP(target exit data map(always,release:Mat% cve))

     UMPIRE_DEVICE_POOL_FREE(Mat% rho)
     TOMP(target exit data map(always,release:Mat% rho))

     UMPIRE_DEVICE_POOL_FREE(Mat% nez)
     TOMP(target exit data map(always,release:Mat% nez))

     UMPIRE_DEVICE_POOL_FREE(Mat% stimComptonMult)
     TOMP(target exit data map(always,release:Mat% stimComptonMult))

     UMPIRE_DEVICE_POOL_FREE(Mat% Siga)
     TOMP(target exit data map(always,release:Mat% Siga))

     UMPIRE_DEVICE_POOL_FREE(Mat% Sigs)
     TOMP(target exit data map(always,release:Mat% Sigs))

     UMPIRE_DEVICE_POOL_FREE(Mat% Eta)
     TOMP(target exit data map(always,release:Mat% Eta))

     UMPIRE_DEVICE_POOL_FREE(Mat% EmissionRate)
     TOMP(target exit data map(always,release:Mat% EmissionRate))

     UMPIRE_DEVICE_POOL_FREE(Mat% SMatEff)
     TOMP(target exit data map(always,release:Mat% SMatEff))

     UMPIRE_DEVICE_POOL_FREE(Mat% PowerEmitted)
     TOMP(target exit data map(always,release:Mat% PowerEmitted))

     UMPIRE_DEVICE_POOL_FREE(Mat% PowerCompton)
     TOMP(target exit data map(always,release:Mat% PowerCompton))

     UMPIRE_DEVICE_POOL_FREE(Mat% nonLinearIterations)
     TOMP(target exit data map(always,release:Mat% nonLinearIterations))

     TOMP(target exit data map(release:Mat))

! IF THESE ARE UNALLOCATED WILL CRASH?
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     if (getComptonFlag(Compton) /= comptonType_None) then
       UMPIRE_DEVICE_POOL_FREE(Compton% gamMean)
       TOMP(target exit data map(always,release:Compton% gamMean))

       UMPIRE_DEVICE_POOL_FREE(Compton% gamSqdDGam)
       TOMP(target exit data map(always,release:Compton% gamSqdDGam))

       UMPIRE_DEVICE_POOL_FREE(Compton% gamCubedDGam)
       TOMP(target exit data map(always,release:Compton% gamCubedDGam))

       UMPIRE_DEVICE_POOL_FREE(Compton% gamD)
       TOMP(target exit data map(always,release:Compton% gamD))

     endif

     TOMP(target exit data map(release:Compton))
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
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
       call finalizeGPUMemory(setID)
       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% AngleOrder)
       TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% AngleOrder))

#endif
     endif

!  Release Dynamic Memory allocated at the beginning of the time step

     if (Size% ndim > 1) then
       call Set%destructDynMemory()
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

     TOMP(target exit data map(always,release:Quad%AngSetPtr))

     UMPIRE_DEVICE_POOL_FREE(Quad%GrpSetPtr)
     TOMP(target exit data map(always,release:Quad%GrpSetPtr))

     UMPIRE_DEVICE_POOL_FREE(Quad%SetDataPtr)
     TOMP(target exit data map(always,release:Quad%SetDataPtr))

     TOMP(target exit data map(release:Quad))

   endif

!  Update radiation energy density

   call setEnergyDensity

!  Release Communication buffers and hyperplanes

   if (Size% ndim >= 2) then

     AngleSetLoop: do setID=1,nAngleSets+nGTASets
       ASet => getAngleSetData(Quad, setID)

       call destructHyperPlane(ASet)
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

