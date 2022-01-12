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
   use GreyAcceleration_mod
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   use ComptonControl_mod
#endif
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

!$omp parallel do private(aSetID, angle) shared(nGTASets, nAngleSets, Quad, Size) schedule(static) default(none)

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

       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% AfpNorm))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% AezNorm))
       TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% ANormSum))


       if (Size% ndim == 2) then

         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% angDerivFac))
         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% quadTauW1))
         TOMP(target exit data map(release:Quad% AngSetPtr(aSetID)% quadTauW2))

       endif

     enddo
!$omp end parallel do

     TOMP(target exit data map(release:Geom% Volume))
     TOMP(target exit data map(release:Geom% VolumeOld))
     TOMP(target exit data map(release:Geom% cOffSet))
     TOMP(target exit data map(release:Geom% numCorner))
     TOMP(target exit data map(release:Geom% corner1))
     TOMP(target exit data map(release:Geom% corner2))
     TOMP(target exit data map(release:Geom% zone1))
     TOMP(target exit data map(release:Geom% zone2))
     TOMP(target exit data map(release:Geom% cEZ))
     TOMP(target exit data map(release:Geom% cFP))
     TOMP(target exit data map(release:Geom% A_ez))
     TOMP(target exit data map(release:Geom% A_fp))
     TOMP(target exit data map(release:Geom% PhiTotal))

     if (Size% ndim == 2) then

       TOMP(target exit data map(release:Geom% Area))
       TOMP(target exit data map(release:Geom% RadiusEZ))
       TOMP(target exit data map(release:Geom% RadiusFP))

     elseif (Size% ndim == 3) then

       TOMP(target exit data map(release:Geom% nCFacesArray))

     endif

     TOMP(target exit data map(release:Geom))

!    GTA
     if (Size%useNewGTASolver) then
       TOMP(target exit data map(release:GTA% TT))
       TOMP(target exit data map(release:GTA% Pvv))
       TOMP(target exit data map(release:GTA% GreySigTotal))
       TOMP(target exit data map(release:GTA% GreySigScat))
       TOMP(target exit data map(release:GTA% GreySigtInv))
       TOMP(target exit data map(release:GTA% GreySource))
       TOMP(target exit data map(release:GTA% PhiInc))
       TOMP(target exit data map(release:GTA% Q))
       TOMP(target exit data map(release:GTA% TsaSource))

       if (Size% ndim == 2) then
         TOMP(target exit data map(release:GTA% Tvv))
       endif

       TOMP(target exit data map(release:GTA))
     endif

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

   endif !endif useGPU

!  Deallocation for Communication Sets 
   CommSetLoop: do cSetID=1,nCommSets+nGTASets
     CSet  => getCommSetData(Quad, cSetID)
     call destructComm(CSet)
   enddo CommSetLoop

!$omp parallel do private(setID, Set) schedule(static)

   SetLoop: do setID=1,nSets

     Set => getSetData(Quad, setID)

!    Update Psi on the host and Release GPU Memory

     if ( Size% useGPU ) then
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
       call finalizeGPUMemory(setID)
#endif
     endif

     if (Size% ndim > 1) then
       call Set%destructCyclePsi()
     endif

   enddo SetLoop

!$omp end parallel do

!  Update boundary edits

!$omp parallel do private(setID) schedule(static)
   do setID=1,nSets
     call BoundaryEdit(setID)
   enddo
!$omp end parallel do

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
   if (useBoltzmannCompton .AND. Size% useCUDASolver) then
     call freeGpuMemory ()
   endif
#  endif
#endif

   return
   end subroutine finalizeSets 

