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
     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad%SetDataPtr)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad%GrpSetPtr)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad%AngSetPtr)

!    Map Group Sets
     do gSetID=1,nGroupSets
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% GrpSetPtr(gSetID)% STotal)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% GrpSetPtr(gSetID)% Sigt)
     enddo

!    Map ZoneSets

     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% AL)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% AU)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% nCornerSet)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% nCornerBatch)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% offset)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% cornerList)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% cornerMap)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% zoneList)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% cornerConverged)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% Te)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% TeOld)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% delta)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% sumT)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% netRate)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% dTCompton)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% B)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% dBdT)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% Snu0)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% dSnu0dT)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% AD)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% z)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% fk2)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% nI)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% nS)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% ex)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% expPH)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% comptonDeltaEr)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% dComptonDT)
     TOMP_TARGET_ENTER_DATA_MAP_TO(ZSet% comptonSe)


     if (Options%getMPIUseDeviceAddresses()) then
  !    Map Comm Sets
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr)

       CommSetLoop1: do cSetID=1,nCommSets+nGTASets

!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% NangBinList)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% AngleToBin)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% AngleOrder)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% AngleOrder0)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% RecvOrder0)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% RecvOrder)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% request)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% IncFlux)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% IncFluxOld)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% NetFlux)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% relError)
!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% Converged)

         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% CommPtr)

         ! Loop over communicators
         do commID=1,Size%ncomm
           do angle=1,Quad%CommSetPtr(cSetID)%NumAngles
             if (Quad%CommSetPtr(cSetID)% CommPtr(commID, angle)%nSend > 0) then
                 TOMP_TARGET_ENTER_DATA_MAP_ALLOC(Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% psibsend)
!              TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% ListSend)
             endif

             if (Quad%CommSetPtr(cSetID)% CommPtr(commID, angle)%nRecv > 0) then
!             TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% ListRecv)
               TOMP_TARGET_ENTER_DATA_MAP_ALLOC(Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% psibrecv)
             endif

!           TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% irequest)
           end do
         end do

!       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% CommFluxPtr)
       ! Loops over communicator fluxes
!         do commID=1,Size%ncomm
!           TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% CommFluxPtr(commID)% irequest)
!           TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% CommFluxPtr(commID)% IncFlux)
!           TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% CommSetPtr(cSetID)% CommFluxPtr(commID)% ExitFlux)
!         enddo
       enddo CommSetLoop1

     endif ! Options%getMPIUseDeviceAddresses()

!    Map Angle Sets
     do aSetID=1,nAngleSets+nGTASets

       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% nextZ)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% nextC)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% StartingDirection)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% FinishingDirection)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% Omega)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% Weight)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% numCycles)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% cycleOffSet)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% cycleList)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% nHyperPlanes)

       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% HypPlanePtr)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% BdyExitPtr)

       if ( aSetID <= nAngleSets ) then
         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% AfpNorm)
         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% AezNorm)
         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% ANormSum)
       endif

       do angle=1,Quad% AngSetPtr(aSetID)% numAngles

         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% BdyExitPtr(angle)% bdyList)

         if ( .not. Quad% AngSetPtr(aSetID)% FinishingDirection(angle) ) then
           TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% zonesInPlane)
           TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane1)
           TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane2)
           TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% ndone)
         endif

       enddo

       if (Size% ndim == 2) then
         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% angDerivFac)
         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% quadTauW1)
         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% AngSetPtr(aSetID)% quadTauW2)
       endif

     enddo

!    Geometry

     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% Volume)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% VolumeOld)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% VolumeZone)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% cOffSet)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% numCorner)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% CToZone)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% corner1)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% corner2)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% zone1)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% zone2)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% cEZ)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% cFP)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% A_ez)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% A_fp)

     if (Size% ndim == 2) then
       TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% Area)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% RadiusEZ)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% RadiusFP)
     elseif (Size% ndim == 3) then
       TOMP_TARGET_ENTER_DATA_MAP_TO(Geom% nCFacesArray)
     endif

!    Radiation Intensity

     TOMP_TARGET_ENTER_DATA_MAP_TO(Rad)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Rad% PhiTotal)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Rad% radEnergy)

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Compton)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Compton% gamMean)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Compton% gamSqdDGam)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Compton% gamCubedDGam)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Compton% gamD)
#endif

!    GTA

     TOMP_TARGET_ENTER_DATA_MAP_TO(GTA)
     TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% GreySource)
     TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% GreyCorrection)
     TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% Chi)

     if (Size%useNewGTASolver) then
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% TT)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% Pvv)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% GreySigTotal)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% GreySigScat)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% GreySigtInv)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% PhiInc)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% Q)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% TsaSource)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% AfpNorm)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% AezNorm)
        TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% ANormSum)

        if (Size% ndim == 2) then
          TOMP_TARGET_ENTER_DATA_MAP_TO(GTA% Tvv)
        endif
     endif
   endif ! Size%useGPU

!  Initialize communication handles for persistent communicators

!  Moved from findexit routine.   This needs to be done after the comm data is
!  mapped to the GPU, as these are set up with the device addresses of these buffers.
!  -- Aaron
!  QUESTION - We're passing in angle set IDs, but inside the initcomm the
!  parameter is 'cSetID'.  Should this be a loop over comm sets or angle sets??
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
     do cSetID=1,nCommSets
       CSet => getCommSetData(Quad, cSetID)
       do setID=CSet% set1,CSet% set2
         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% AngleOrder)
       enddo
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
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% Tec)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% Tecn)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% denec)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% cve)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% rho)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% nez)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% stimComptonMult)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% Siga)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% Sigs)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% Eta)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% EmissionRate)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% SMatEff)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% PowerEmitted)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% PowerCompton)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Mat% nonLinearIterations)
   endif

!  Map GTA set variables

   if ( Size% useGPU ) then

     do setID=nSets+1,nSets+nGTASets
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% AngleOrder)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% tPsi)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% pInc)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% src)

       if (Size% ndim == 2) then
         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% tPsiM)
         TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% tInc)
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
