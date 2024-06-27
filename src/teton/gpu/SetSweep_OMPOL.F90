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
   if ( useBoltzmannCompton .and. Size% useCUDASolver .and. Size% ngr >= 16) then
     do setID=1,nGroupSets
       GSet   => getGroupSetData(Quad, setID)
       TOMP_UPDATE(target update to (GSet% STotal))
     enddo
   endif
#endif

!  At this point, all sets must have the same number of angles

   Set          => getSetData(Quad, 1)
   NumAnglesDyn =  Set% NumAnglesDyn

!  If this is the first flux iteration, initialize the communication
!  order and incident flux on shared boundaries

   START_RANGE("Teton_Init_Comm_Order")
   do cSetID=1,nCommSets

     CSet => getCommSetData(Quad, cSetID)

     call restoreCommOrder(CSet)
     call setIncidentFlux(cSetID)
   enddo
   END_RANGE("Teton_Init_Comm_Order")

!  Begin Flux Iteration 

   fluxIter         = 0
   FluxConverged(:) = .FALSE.

   FluxIteration: do

     fluxIter = fluxIter + 1

!    Post receives for all data

     START_RANGE("Teton_Comm_Boundary_Fluxes")
     do cSetID=1,nCommSets
       call InitExchange(cSetID)
     enddo
     END_RANGE("Teton_Comm_Boundary_Fluxes")

!    Loop over angles, solving for each in turn:

     AngleLoop: do sendIndex=1,NumAnglesDyn

       START_RANGE("Teton_Comm_Boundary_Fluxes")

!!$omp parallel do default(none) schedule(dynamic) &
!!$omp& shared(nCommSets,Quad,SnSweep, sendIndex) &
!!$omp& private(CSet,Angle)
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
!!$omp end parallel do
       END_RANGE("Teton_Comm_Boundary_Fluxes")

       do setID=1,nSets

         Set   => getSetData(Quad, setID)
         Angle =  Set% AngleOrder(sendIndex)

!        Update incident fluxes on reflecting boundaries

         START_RANGE("Teton_Update_Reflecting_Fluxes")
         call snreflect(SnSweep, setID, Angle)
         END_RANGE("Teton_Update_Reflecting_Fluxes")

!  Map the latest boundary values

         START_RANGE("Teton_OpenMP_Updates")
         TOMP_UPDATE(target update to( Set%PsiB(:,:,Angle) ) )
         END_RANGE("Teton_OpenMP_Updates")

       enddo

!      Sweep the mesh, calculating PSI for each corner; the 
!      boundary flux array PSIB is also updated here. 

       ASet  => getAngleSetData(Quad, 1)

       AngleType: if ( .not. ASet% FinishingDirection(Angle) ) then

         time1 = MPIWtime()
         START_RANGE("Teton_Calc_Rad_Energy_Density")

         if (sweepVersion == 0) then

           if (ndim == 3) then
             call SweepUCBxyz_GPU(nSets, sendIndex, savePsi)
           elseif (ndim == 2) then
             call SweepUCBrz_GPU(nSets, sendIndex, savePsi)
           endif

         elseif (sweepVersion == 1) then

           if (ndim == 3) then
             call CornerSweepUCBxyz_GPU(nSets, sendIndex, savePsi)
           elseif (ndim == 2) then
             call CornerSweepUCBrz_GPU(nSets, sendIndex, savePsi)
           endif

         else
           TETON_FATAL("Invalid value set for Sweep kernel version to use.")
         endif
         END_RANGE("Teton_Calc_Rad_Energy_Density")

!    Update the total scalar intensity on the GPU

         START_RANGE("Teton_Update_Total_Scalar_Intensity")
         call getPhiTotal(sendIndex)
         END_RANGE("Teton_Update_Total_Scalar_Intensity")

         time2 = MPIWtime()
         dtime = (time2 - time1)/sixty
         Size%GPUSweepTimeCycle = Size%GPUSweepTimeCycle + dtime

       endif AngleType

       do setID=1,nSets

         Set   => getSetData(Quad, setID)
         Angle =  Set% AngleOrder(sendIndex)

         START_RANGE("Teton_OpenMP_Updates")
         TOMP_UPDATE(target update from( Set%PsiB(:,:,Angle) ))
         END_RANGE("Teton_OpenMP_Updates")

       enddo

     enddo AngleLoop

!    Test convergence of incident fluxes

     START_RANGE("Teton_Sweep_Test_Conv")
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
     END_RANGE("Teton_Sweep_Test_Conv")

!    If this is the end of the radiation step and we are saving Psi do
!    not perform additional sweeps

     if (savePsi) then
       exit FluxIteration
     endif

     nConv = 0
     do cSetID=1,nCommSets
       if ( FluxConverged(cSetID) ) then
         nConv = nConv + 1
       endif
     enddo

     nNotConv = nCommSets - nConv

!    Make sure all processes are in sync
     call MPIAllReduce(nNotConv, "max", MY_COMM_GROUP)

     if ( nNotConv == 0 .or. fluxIter >= maxIters ) then
       exit FluxIteration
     else
       cycle FluxIteration
     endif

   enddo FluxIteration


   deallocate( FluxConverged )


   return
   end subroutine SetSweep_GPU 


