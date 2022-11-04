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

     TOMP(target exit data map(release: ZSet% nCornerSet))
     TOMP(target exit data map(release: ZSet% nCornerBatch))
     TOMP(target exit data map(release: ZSet% offset))
     TOMP(target exit data map(release: ZSet% cornerList))
     TOMP(target exit data map(release: ZSet% cornerMap))
     TOMP(target exit data map(release: ZSet% zoneList))
     TOMP(target exit data map(release: ZSet% cornerConverged))
     TOMP(target exit data map(release: ZSet% Te))
     TOMP(target exit data map(release: ZSet% TeOld))
     TOMP(target exit data map(release: ZSet% delta))
     TOMP(target exit data map(release: ZSet% sumT))
     TOMP(target exit data map(release: ZSet% netRate))
     TOMP(target exit data map(release: ZSet% dTCompton))
     TOMP(target exit data map(release: ZSet% B))
     TOMP(target exit data map(release: ZSet% dBdT))
     TOMP(target exit data map(release: ZSet% Snu0))
     TOMP(target exit data map(release: ZSet% dSnu0dT))
     TOMP(target exit data map(release: ZSet% AD))
     TOMP(target exit data map(release: ZSet% z))
     TOMP(target exit data map(release: ZSet% fk2))
     TOMP(target exit data map(release: ZSet% nI))
     TOMP(target exit data map(release: ZSet% nS))
     TOMP(target exit data map(release: ZSet% ex))
     TOMP(target exit data map(release: ZSet% expPH))
     TOMP(target exit data map(release: ZSet% comptonDeltaEr))
     TOMP(target exit data map(release: ZSet% dComptonDT))
     TOMP(target exit data map(release: ZSet% comptonSe))
     TOMP(target exit data map(release: ZSet% AU))
     TOMP(target exit data map(release: ZSet% AL))
     TOMP(target exit data map(release: ZSet))

     ! Unmap group sets
     do gSetID=1,nGroupSets
       TOMP(target exit data map(release:Quad% GrpSetPtr(gSetID)% STotal))
       TOMP(target exit data map(release:Quad% GrpSetPtr(gSetID)% Sigt))
     enddo

! Leave as-of-yet unneeded maps in here as comments.
! I anticipate some of these will be needed as more code is pushed to the GPU for the MPI.
! -- Aaron
     if (Options%getMPIUseDeviceAddresses()) then
       ! Unmap comm sets
       do cSetID=1,nCommSets+nGTASets

         ! Loop over communicator fluxes
!         do commID=1,Size%ncomm
!           TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommFluxPtr(commID)% irequest))
!           TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommFluxPtr(commID)% IncFlux))
!           TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommFluxPtr(commID)% ExitFlux))
!         enddo

!        TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommFluxPtr))

         ! Loop over communicators
         do commID=1,Size%ncomm
           do angle=1,Quad%CommSetPtr(cSetID)% NumAngles
             if (Quad%CommSetPtr(cSetID)% CommPtr(commID, angle)%nSend > 0) then
!              TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% ListSend))
               TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% psibsend))
             endif
             if (Quad%CommSetPtr(cSetID)% CommPtr(commID, angle)%nRecv > 0) then
!               TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% ListRecv))
               TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% psibrecv))
             endif
!             TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommPtr(commID, angle)% irequest))
           end do
         end do
         TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% CommPtr))

!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% NangBinList))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% AngleToBin))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% AngleOrder))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% AngleOrder0))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% RecvOrder0))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% RecvOrder))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% request))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% IncFlux))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% IncFluxOld))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% NetFlux))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% relError))
!       TOMP(target exit data map(release:Quad% CommSetPtr(cSetID)% Converged))

         TOMP(target exit data map(release:Quad% CommSetPtr))
       enddo
     endif ! Options%getMPIUseDeviceAddresses()

     do aSetID=1,nAngleSets+nGTASets

       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% nextZ))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% nextC))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% Omega))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% Weight))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% numCycles))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% cycleOffSet))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% cycleList))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% nHyperPlanes))

       do angle=1,Quad%AngSetPtr(aSetID)% numAngles
         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% BdyExitPtr(angle)%bdyList))

         if ( .not. Quad%AngSetPtr(aSetID)% FinishingDirection(angle) ) then
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% zonesInPlane))
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane1))
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% hplane2))
           TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr(angle)% ndone))
         endif
       enddo

       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% StartingDirection))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% FinishingDirection))

       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% HypPlanePtr))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% BdyExitPtr))

       if ( aSetID <= nAngleSets ) then
         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% AfpNorm))
         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% AezNorm))
         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% ANormSum))
       endif


       if (Size% ndim == 2) then

         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% angDerivFac))
         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% quadTauW1))
         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% quadTauW2))

       endif

     enddo

!    Geometry

     TOMP(target exit data map(release:Geom% Volume))
     TOMP(target exit data map(release:Geom% VolumeOld))
     TOMP(target exit data map(release:Geom% VolumeZone))
     TOMP(target exit data map(release:Geom% cOffSet))
     TOMP(target exit data map(release:Geom% numCorner))
     TOMP(target exit data map(release:Geom% CToZone))
     TOMP(target exit data map(release:Geom% corner1))
     TOMP(target exit data map(release:Geom% corner2))
     TOMP(target exit data map(release:Geom% zone1))
     TOMP(target exit data map(release:Geom% zone2))
     TOMP(target exit data map(release:Geom% cEZ))
     TOMP(target exit data map(release:Geom% cFP))
     TOMP(target exit data map(release:Geom% A_ez))
     TOMP(target exit data map(release:Geom% A_fp))

     if (Size% ndim == 2) then
       TOMP(target exit data map(release:Geom% Area))
       TOMP(target exit data map(release:Geom% RadiusEZ))
       TOMP(target exit data map(release:Geom% RadiusFP))
     elseif (Size% ndim == 3) then
       TOMP(target exit data map(release:Geom% nCFacesArray))
     endif

     TOMP(target exit data map(release:Geom))

!    Radiation Intensity

     TOMP(target exit data map(release:Rad% PhiTotal))
     TOMP(target exit data map(release:Rad% radEnergy))
     TOMP(target exit data map(release:Rad))

!    GTA

     if (Size%useNewGTASolver) then
       TOMP(target exit data map(release:GTA% TT))
       TOMP(target exit data map(release:GTA% Pvv))
       TOMP(target exit data map(release:GTA% GreySigTotal))
       TOMP(target exit data map(release:GTA% GreySigScat))
       TOMP(target exit data map(release:GTA% GreySigtInv))
       TOMP(target exit data map(release:GTA% PhiInc))
       TOMP(target exit data map(release:GTA% Q))
       TOMP(target exit data map(release:GTA% TsaSource))
       TOMP(target exit data map(release:GTA% AfpNorm))
       TOMP(target exit data map(release:GTA% AezNorm))
       TOMP(target exit data map(release:GTA% ANormSum))

       if (Size% ndim == 2) then
         TOMP(target exit data map(release:GTA% Tvv))
       endif
     endif

     TOMP(target exit data map(release:GTA% GreySource))
     TOMP(target exit data map(release:GTA% GreyCorrection))
     TOMP(target exit data map(release:GTA% Chi))
     TOMP(target exit data map(release:GTA))

     do setID=nSets+1,nSets+nGTASets
       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% AngleOrder))
       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% tPsi))
       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% pInc))
       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% src))

       if (Size% ndim == 2) then
         TOMP(target exit data map(release:Quad% SetDataPtr(setID)% tPsiM))
         TOMP(target exit data map(release:Quad% SetDataPtr(setID)% tInc))
       endif
     enddo

!    Material

     TOMP(target exit data map(release:Mat% Tec))
     TOMP(target exit data map(release:Mat% Tecn))
     TOMP(target exit data map(release:Mat% denec))
     TOMP(target exit data map(release:Mat% cve))
     TOMP(target exit data map(release:Mat% rho))
     TOMP(target exit data map(release:Mat% nez))
     TOMP(target exit data map(release:Mat% stimComptonMult))
     TOMP(target exit data map(release:Mat% Siga))
     TOMP(target exit data map(release:Mat% Sigs))
     TOMP(target exit data map(release:Mat% Eta))
     TOMP(target exit data map(release:Mat% EmissionRate))
     TOMP(target exit data map(release:Mat% SMatEff))
     TOMP(target exit data map(release:Mat% PowerEmitted))
     TOMP(target exit data map(release:Mat% PowerCompton))
     TOMP(target exit data map(release:Mat% nonLinearIterations))
     TOMP(target exit data map(release:Mat))

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     TOMP(target exit data map(release:Compton% gamMean))
     TOMP(target exit data map(release:Compton% gamSqdDGam))
     TOMP(target exit data map(release:Compton% gamCubedDGam))
     TOMP(target exit data map(release:Compton% gamD))
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

     TOMP(target exit data map(release:Quad%AngSetPtr))
     TOMP(target exit data map(release:Quad%GrpSetPtr))
     TOMP(target exit data map(release:Quad%SetDataPtr))
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

