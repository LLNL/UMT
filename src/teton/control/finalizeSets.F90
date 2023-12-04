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

     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% nCornerSet)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% nCornerBatch)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% offset)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% cornerList)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% cornerMap)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% zoneList)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% cornerConverged)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% Te)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% TeOld)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% delta)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% sumT)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% netRate)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% dTCompton)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% B)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% dBdT)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% Snu0)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% dSnu0dT)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% AD)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% z)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% fk2)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% nI)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% nS)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% ex)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% expPH)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% comptonDeltaEr)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% dComptonDT)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% comptonSe)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% AU)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(ZSet% AL)
     TOMP(target exit data map(release: ZSet))

     ! Unmap group sets
     do gSetID=1,nGroupSets
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% GrpSetPtr(gSetID)% STotal)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% GrpSetPtr(gSetID)% Sigt)
     enddo

     do aSetID=1,nAngleSets+nGTASets

       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% nextZ)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% nextC)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% StartingDirection)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% FinishingDirection)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% Omega)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% Weight)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% numCycles)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% cycleOffSet)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% cycleList)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% nHyperPlanes)

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

       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% HypPlanePtr)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% BdyExitPtr)

       if ( aSetID <= nAngleSets ) then
         TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% AfpNorm)
         TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% AezNorm)
         TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% ANormSum)
       endif


       if (Size% ndim == 2) then

         TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% angDerivFac)
         TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% quadTauW1)
         TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% AngSetPtr(aSetID)% quadTauW2)

       endif

     enddo

!    Geometry

     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% Volume)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% VolumeOld)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% VolumeZone)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% cOffSet)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% numCorner)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% CToZone)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% corner1)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% corner2)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% zone1)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% zone2)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% cEZ)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% cFP)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% A_ez)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% A_fp)

     if (Size% ndim == 2) then
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% Area)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% RadiusEZ)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% RadiusFP)
     elseif (Size% ndim == 3) then
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Geom% nCFacesArray)
     endif

     TOMP(target exit data map(release:Geom))

!    Radiation Intensity

     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Rad% PhiTotal)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Rad% radEnergy)
     TOMP(target exit data map(release:Rad))

!    GTA

     if (Size%useNewGTASolver) then
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% TT)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% Pvv)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% GreySigTotal)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% GreySigScat)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% GreySigScatVol)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% GreySigtInv)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% PhiInc)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% Q)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% TsaSource)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% AfpNorm)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% AezNorm)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% ANormSum)

       if (Size% ndim == 2) then
         TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% Tvv)
       endif
     endif

     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% GreySource)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% GreyCorrection)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(GTA% Chi)
     TOMP(target exit data map(release:GTA))

     do setID=nSets+1,nSets+nGTASets
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% AngleOrder)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% tPsi)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% pInc)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% src)

       if (Size% ndim == 2) then
         TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% tPsiM)
         TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% tInc)
       endif
     enddo

!    Material

     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% Tec)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% Tecn)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% denec)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% cve)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% rho)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% nez)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% stimComptonMult)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% Siga)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% Sigs)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% Eta)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% EmissionRate)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% SMatEff)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% PowerEmitted)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% PowerCompton)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Mat% nonLinearIterations)
     TOMP(target exit data map(release:Mat))

! IF THESE ARE UNALLOCATED WILL CRASH?
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     if (getComptonFlag(Compton) /= comptonType_None) then
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Compton% gamMean)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Compton% gamSqdDGam)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Compton% gamCubedDGam)
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Compton% gamD)
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
       TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% AngleOrder)
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

     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad%AngSetPtr)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad%GrpSetPtr)
     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad%SetDataPtr)
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

