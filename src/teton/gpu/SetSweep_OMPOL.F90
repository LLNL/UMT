#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Version 1:  09/2017, PFN                      *
!                                                                      *
!   SetSweep_GPU  - This routine controls the computational set        *
!                   sweeps and communication when using a GPU          *
!                   accelerator for the sweeps.                        *
!                                                                      *
!***********************************************************************

   subroutine SetSweep_GPU(savePsi)


   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use Geometry_mod
   use SetData_mod
   use CommSet_mod
   use AngleSet_mod
   use GroupSet_mod
   use mpi_param_mod
   use mpif90_mod
   use iter_control_list_mod
   use iter_control_mod
   use OMPWrappers_mod
   use Options_mod
   use MemoryAllocator_mod, only : Allocator

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   use ComptonControl_mod
#endif

#if defined(TETON_ENABLE_CALIPER)
   use caliper_mod
#endif

   implicit none

!  Arguments

   logical (kind=1), intent(in)  :: savePsi

!  Local

   type(SetData),     pointer    :: Set  => null()
   type(CommSet),     pointer    :: CSet => null()
   type(AngleSet),    pointer    :: ASet => null()
   type(GroupSet),    pointer    :: GSet => null()
   type(IterControl), pointer    :: incidentFluxControl => NULL()

   integer                       :: nSets
   integer                       :: setID
   integer                       :: cSetID
   integer                       :: Angle
   integer                       :: ndim 
   integer                       :: sendIndex
   integer                       :: fluxIter
   integer                       :: NumAnglesDyn
   integer                       :: nGroupSets
   integer                       :: nCommSets
   integer                       :: nConv
   integer                       :: nNotConv
   integer                       :: maxIters
   integer                       :: sweepVersion

   real(adqt)                    :: time1
   real(adqt)                    :: time2
   real(adqt)                    :: dtime

   logical (kind=1)              :: SnSweep 

   logical (kind=1), allocatable :: FluxConverged(:)
   logical(kind=1)               :: useBoltzmannCompton

   logical                       :: commDataOnCPU

!  Constants

   incidentFluxControl => getIterationControl(IterControls,"incidentFlux")
   maxIters            =  getMaxNumberOfIterations(incidentFluxControl)
   nSets               =  getNumberOfSets(Quad)
   nGroupSets          =  getNumberOfGroupSets(Quad)
   nCommSets           =  getNumberOfCommSets(Quad)
   ndim                =  Size% ndim
   SnSweep             = .TRUE.
   sweepVersion        =  Options%getSweepVersion()

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   useBoltzmannCompton = getUseBoltzmann(Compton)
#endif

   allocate( FluxConverged(nCommSets) )

!  If the CUDA solver is used the source needs to be mapped to the GPU

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
#if !defined(TETON_OPENMP_HAS_UNIFIED_MEMORY)
   if ( useBoltzmannCompton .and. Size% useCUDASolver .and. Size% ngr >= 16) then
   START_RANGE("Teton_OpenMP_data_movement")
     do setID=1,nGroupSets
       GSet   => getGroupSetData(Quad, setID)
       TOMP_UPDATE(target update to (GSet% STotal))
     enddo
   END_RANGE("Teton_OpenMP_data_movement")
   endif
#endif
#endif

!  At this point, all sets must have the same number of angles

   Set          => getSetData(Quad, 1)
   NumAnglesDyn =  Set% NumAnglesDyn

!  If this is the first flux iteration, initialize the communication
!  order and incident flux on shared boundaries
   do cSetID=1,nCommSets

     CSet => getCommSetData(Quad, cSetID)

     START_RANGE("Teton_Comm_Init_Order")
     call restoreCommOrder(CSet)
     END_RANGE("Teton_Comm_Init_Order")
     START_RANGE("Teton_SetIncidentFlux")
     call setIncidentFlux(cSetID)
     END_RANGE("Teton_SetIncidentFlux")
   enddo

!  Begin Flux Iteration 

   START_RANGE("Teton_Flux_Iter_Loop")
   fluxIter         = 0
   FluxConverged(:) = .FALSE.

   FluxIteration: do

     fluxIter = fluxIter + 1

!    Post receives for all data

     START_RANGE("Teton_Comm_Post_Receives")
     do cSetID=1,nCommSets
       call InitExchange(cSetID)
     enddo
     END_RANGE("Teton_Comm_Post_Receives")

!    Loop over angles, solving for each in turn:

     AngleLoop: do sendIndex=1,NumAnglesDyn
       START_RANGE("Teton_Comm_Send_Recv_Fluxes")
!$omp parallel do default(none) schedule(dynamic) &
!$omp& shared(nCommSets,Quad,SnSweep, sendIndex) &
!$omp& private(CSet,Angle)
       do cSetID=1,nCommSets

         CSet  => getCommSetData(Quad, cSetID)
         Angle =  CSet% AngleOrder(sendIndex)

!        Send the boundary information needed by my neighbors
         call SendFlux(SnSweep, cSetID, sendIndex)

!        Test for completion of the sends needed by my neighbors
         call TestSend(cSetID, sendIndex)

!        Receive the boundary information needed to compute this angle
         call RecvFlux(SnSweep, cSetID, Angle)

       enddo
!$omp end parallel do
       END_RANGE("Teton_Comm_Send_Recv_Fluxes")

       do setID=1,nSets

         Set   => getSetData(Quad, setID)
         Angle =  Set% AngleOrder(sendIndex)

!        Update incident fluxes on reflecting boundaries

         START_RANGE("Teton_Update_Fluxes_On_Refl_Boundaries")
         call snreflect(SnSweep, setID, Angle)
         END_RANGE("Teton_Update_Fluxes_On_Refl_Boundaries")

!  Map the latest boundary values

#if !defined(TETON_OPENMP_HAS_UNIFIED_MEMORY)
         START_RANGE("Teton_OpenMP_data_movement")
         TOMP_UPDATE(target update to( Set%PsiB(:,:,Angle) ) )
         END_RANGE("Teton_OpenMP_data_movement")
#endif
       enddo

!      Sweep the mesh, calculating PSI for each corner; the 
!      boundary flux array PSIB is also updated here. 

       ASet  => getAngleSetData(Quad, 1)

       AngleType: if ( .not. ASet% FinishingDirection(Angle) ) then

         time1 = MPIWtime()
         START_RANGE("Teton_Sweep_Kernel")

         if (sweepVersion == 1) then

           if (ndim == 3) then
             call SweepUCBxyz_GPU(nSets, sendIndex, savePsi)
           elseif (ndim == 2) then
             call SweepUCBrz_GPU(nSets, sendIndex, savePsi)
           endif

         elseif (sweepVersion == 2) then

           if (ndim == 3) then
             call CornerSweepUCBxyz_GPU(nSets, sendIndex, savePsi)
           elseif (ndim == 2) then
             call CornerSweepUCBrz_GPU(nSets, sendIndex, savePsi)
           endif

         else
           TETON_FATAL("Invalid value set for Sweep kernel version to use.")
         endif
         END_RANGE("Teton_Sweep_Kernel")

!    Update the total scalar intensity on the GPU

         START_RANGE("Teton_Total_Scalar_Intensity")
         call getPhiTotal(sendIndex)
         END_RANGE("Teton_Total_Scalar_Intensity")

         time2 = MPIWtime()
         dtime = (time2 - time1)/sixty
         Size%GPUSweepTimeCycle = Size%GPUSweepTimeCycle + dtime

       endif AngleType

       do setID=1,nSets

         Set   => getSetData(Quad, setID)
         Angle =  Set% AngleOrder(sendIndex)

#if !defined(TETON_OPENMP_HAS_UNIFIED_MEMORY)
         START_RANGE("Teton_OpenMP_Updates")
         TOMP_UPDATE(target update from( Set%PsiB(:,:,Angle) ))
         END_RANGE("Teton_OpenMP_Updates")
#endif
       enddo

     enddo AngleLoop

!    Test convergence of incident fluxes

     START_RANGE("Teton_Test_For_Conv")
!$omp parallel do default(none) schedule(static) &
!$omp& shared(nCommSets, FluxConverged)
     do cSetID=1,nCommSets
       call setIncidentFlux(cSetID)
       call testFluxConv(cSetID, FluxConverged(cSetID))

!  Do not reduce the the number of sweep angles as we would for the dynamic
!  iteration; all angles are repeated.

!         CSet => getCommSetData(Quad, cSetID)
!         call setCommOrder(CSet)

     enddo
!$omp end parallel do

!    If this is the end of the radiation step and we are saving Psi do
!    not perform additional sweeps

     END_RANGE("Teton_Test_For_Conv")
     if (savePsi) then
       exit FluxIteration
     endif
     START_RANGE("Teton_Test_For_Conv")

     nConv = 0
     do cSetID=1,nCommSets
       if ( FluxConverged(cSetID) ) then
         nConv = nConv + 1
       endif
     enddo

     nNotConv = nCommSets - nConv

!    Make sure all processes are in sync
     START_RANGE("Teton_Comm_All_Reduce_On_Conv")
     call MPIAllReduce(nNotConv, "max", MY_COMM_GROUP)
     END_RANGE("Teton_Comm_All_Reduce_On_Conv")

     END_RANGE("Teton_Test_For_Conv")

     if ( nNotConv == 0 .or. fluxIter >= maxIters ) then
       exit FluxIteration
     else
       cycle FluxIteration
     endif

   enddo FluxIteration
   END_RANGE("Teton_Flux_Iter_Loop")


   deallocate( FluxConverged )

   return
   end subroutine SetSweep_GPU 


