#include "macros.h"
!***********************************************************************
!                       Last Update:  03/2012, PFN                     *
!                                                                      *
!   RTMAINSN - Control program for SN radiation transport in 1D,       *
!              2D and 3D geometries.                                   *
!                                                                      *
!***********************************************************************

   subroutine rtmainsn

   use, intrinsic :: iso_c_binding, only : c_int
   use Options_mod
   use kind_mod
   use iter_control_list_mod
   use iter_control_mod
   use Size_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use GreyAcceleration_mod
   use constant_mod
   use mpi_param_mod
   use mpif90_mod
   use default_iter_controls_mod, only : outer_slow_conv_threshold

#if defined(TETON_ENABLE_CALIPER)
   use caliper_mod
#endif
#if defined(TETON_ENABLE_OPENMP)
   use omp_lib
#endif

   implicit none

!  Local

   integer    :: tempIter, nTotalSweeps, izero
   integer    :: ndim
   integer    :: nThreadsInitial
   integer    :: maxIterCheck
   integer    :: sweepVersion

   real(adqt) :: maxEnergyDensityError, maxTempError 
   real(adqt) :: time1, time2, dtime, epsilonCheck

   logical(kind=1) :: savePsi
   logical(kind=1) :: SnSweep

   type(IterControl) , pointer :: temperatureControl  => NULL()
   type(IterControl) , pointer :: intensityControl    => NULL()
   type(IterControl) , pointer :: greyControl         => NULL()
   type(IterControl) , pointer :: incidentFluxControl => NULL()
   type(IterControl) , pointer :: nonlinearControl => NULL()

   integer(kind=c_int) :: numOMPThreads

   character(len=512)   :: descriptor

!  Constants:

   parameter (izero=0)

!  Set some scalars used for dimensioning

   ndim                   = Size%ndim

   Size%MatCoupTimeCycle  = zero
   Size%SweepTimeCycle    = zero
   Size%GPUSweepTimeCycle = zero
   Size%GTATimeCycle      = zero

   sweepVersion           = Options%getSweepVersion()

   savePsi = .FALSE.

!  Iteration Controls

   temperatureControl  => getIterationControl(IterControls, "temperature")
   intensityControl    => getIterationControl(IterControls, "intensity")
   greyControl         => getIterationControl(IterControls, "grey")
   incidentFluxControl => getIterationControl(IterControls, "incidentFlux")
   nonlinearControl    => getIterationControl(IterControls, "nonLinear")

!  Check that tolerances have been set identically on all ranks.
   epsilonCheck = getEpsilonPoint(temperatureControl)
   call MPIAllReduce(epsilonCheck, "max", MY_COMM_GROUP)
   TETON_VERIFY(epsilonCheck == getEpsilonPoint(temperatureControl), "Teton: Outer temperature tolerance not identical on all ranks.")
   epsilonCheck = getEpsilonPoint(intensityControl)
   call MPIAllReduce(epsilonCheck, "max", MY_COMM_GROUP)
   TETON_VERIFY(epsilonCheck == getEpsilonPoint(intensityControl), "Teton: Outer energy density tolerance not identical on all ranks.")
   epsilonCheck = getEpsilonPoint(greyControl)
   call MPIAllReduce(epsilonCheck, "max", MY_COMM_GROUP)
   TETON_VERIFY(epsilonCheck == getEpsilonPoint(greyControl), "Teton: Grey tolerance not identical on all ranks.")
   epsilonCheck = getEpsilonPoint(incidentFluxControl)
   call MPIAllReduce(epsilonCheck, "max", MY_COMM_GROUP)
   TETON_VERIFY(epsilonCheck == getEpsilonPoint(incidentFluxControl), "Teton: Incident flux tolerance not identical on all ranks.")
   epsilonCheck = getEpsilonPoint(nonlinearControl)
   call MPIAllReduce(epsilonCheck, "max", MY_COMM_GROUP)
   TETON_VERIFY(epsilonCheck == getEpsilonPoint(nonlinearControl), "Teton: Inner nonlinear tolerance not identical on all ranks.")

   maxIterCheck = getMaxNumberOfIterations(temperatureControl)
   call MPIAllReduce(maxIterCheck, "max", MY_COMM_GROUP)
   TETON_VERIFY(maxIterCheck == getMaxNumberOfIterations(temperatureControl), "Teton: Outer max iters not identical on all ranks.")
   maxIterCheck = getMaxNumberOfIterations(greyControl)
   call MPIAllReduce(maxIterCheck, "max", MY_COMM_GROUP)
   TETON_VERIFY(maxIterCheck == getMaxNumberOfIterations(greyControl), "Teton: Grey max iters not identical on all ranks.")
   maxIterCheck = getMaxNumberOfIterations(incidentFluxControl)
   call MPIAllReduce(maxIterCheck, "max", MY_COMM_GROUP)
   TETON_VERIFY(maxIterCheck == getMaxNumberOfIterations(incidentFluxControl), "Teton: Incident flux max iters not identical on all ranks.")
   maxIterCheck = getMaxNumberOfIterations(nonlinearControl)
   call MPIAllReduce(maxIterCheck, "max", MY_COMM_GROUP)
   TETON_VERIFY(maxIterCheck == getMaxNumberOfIterations(nonlinearControl), "Teton: Inner nonlienar max iters not identical on all ranks.")

!  Initialize counters for this time step
   call setNumberOfIterations(temperatureControl,izero)
   call setNumberOfIterations(intensityControl,izero)
   call setNumberOfIterations(greyControl,izero)

   call setGlobalError(temperatureControl,0.1d0)
   call setGlobalError(intensityControl,0.1d0)
   call setSweepCounters(Quad)

   GTA%epsGrey = 1.0e-5_adqt

!  Initialize threading.
!  Get the initial max threads value, as we will restore that after
!  finishing.
#if defined(TETON_ENABLE_OPENMP)
   nThreadsInitial = omp_get_max_threads()
   numOMPThreads = Options%getNumOmpMaxThreads()
   call omp_set_num_threads(numOMPThreads)
   TETON_VERIFY(omp_get_max_threads() == numOMPThreads, "Teton: Unable to set max threads in OpenMP runtime.")
#endif
!***********************************************************************
!     Nonlinear Solver                                                 *
!***********************************************************************

!  Calculate new electron temperature and energy change

   write(descriptor,'(A63)') "In rtmainsn, before any NL solves or temperature iterations."
   call PrintEnergies(trim(descriptor))

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   time1 = MPIWtime()
   START_RANGE("Teton_NonLinearSolver")
   call NonLinearSolver
   END_RANGE("Teton_NonLinearSolver")

   time2 = MPIWtime()
   dtime = (time2 - time1)/sixty
   Size%MatCoupTimeCycle = Size%MatCoupTimeCycle + dtime
#else
   Size%MatCoupTimeCycle = 0.0
#endif

!***********************************************************************
!     BEGIN IMPLICIT ELECTRON/RADIATION COUPLING ITERATION             *
!***********************************************************************
   tempIter = 0

   TemperatureIteration: do

      tempIter = tempIter + 1

!***********************************************************************
!     BEGIN PHOTON INTENSITY LINEAR ITERATION (only one step for now)  *
!***********************************************************************

     call LinearSolver(savePsi)

     write(descriptor,'(A45,I5,A35)') "In rtmainsn, TemperatureIteration number ", tempIter, ", after call to linear solver"
     call PrintEnergies(trim(descriptor))

!***********************************************************************
!     END PHOTON INTENSITY LINEAR ITERATION                            *
!***********************************************************************

!  Calculate new electron temperature and energy change, including
!  the effect of Compton scattering (if present)

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     time1 = MPIWtime()

     START_RANGE("Teton_NonLinearSolver")
     call NonLinearSolver
     END_RANGE("Teton_NonLinearSolver")

     time2 = MPIWtime()
     dtime = (time2 - time1)/sixty
     Size%MatCoupTimeCycle = Size%MatCoupTimeCycle + dtime
#else
     Size%MatCoupTimeCycle = 0.0
#endif

!***********************************************************************
!     CHECK CONVERGENCE OF RADIATION ENERGY DENSITY AND TEMPERATURE    *
!***********************************************************************
     call ConvergenceTest(maxEnergyDensityError, maxTempError)

     if (tempIter >= getMaxNumberOfIterations(temperatureControl)) then

       exit TemperatureIteration

     else if (tempIter > 1 .and.  &
         maxTempError <  getEpsilonPoint(temperatureControl) .and.  &
         maxEnergyDensityError <  getEpsilonPoint(intensityControl)) then

       if (associated(GTA) .eqv. .true. .and. GTA% forceExtraOuter .and. .not. getConvergenceState(greyControl)) then

          if (Options%isRankVerbose() > 1) then
             print *, "In outer iteration ", tempIter, ", GTA is not converged so Teton will do another outer iteration."
          endif

          GTA%epsGrey = max(maxEnergyDensityError,maxTempError)/20.d0
          cycle TemperatureIteration

       endif

       exit TemperatureIteration

     else

       GTA%epsGrey = max(maxEnergyDensityError,maxTempError)/20.d0

       cycle TemperatureIteration

     endif

   enddo TemperatureIteration

!  We did not save the angle-dependent intensity so generate it here

   SnSweep = .TRUE.
   savePsi = .TRUE.

   call PrintEnergies("In rtmainsn, after exiting TemperatureIteration")

   START_RANGE("Teton_Sweep")
   if (ndim > 1) then
     call ControlSweep(SnSweep, savePsi)
   else
     call ControlSweep1D(savePsi)
   endif
   END_RANGE("Teton_Sweep")

   call PrintEnergies("In rtmainsn, after final ControlSweep to save psi")

!***********************************************************************
!     END IMPLICIT ELECTRON/RADIATION COUPLING ITERATION (OUTER)       *
!***********************************************************************

!***********************************************************************
!     EDITS                                                            *
!***********************************************************************

!  Hook to call optional pre-intensity iteration functions for
!  user-defined sources for test problems (i.e., MMS tests)
#if defined(TETON_ENABLE_INTENSITY_SOLVE_PREPOST_HOOKS)
   call postIntensityIterAction
#endif

!  Update Iteration Counts and sweep timer. Here we subtract the
!  GPU portion of the total sweep time so we can report separate
!  CPU and GPU times

   nTotalSweeps         = getTotalSweeps(Quad)
   Size% SweepTimeCycle = Size% SweepTimeCycle - Size% GPUSweepTimeCycle

   call setNumberOfIterations(temperatureControl,tempIter)
   call setNumberOfIterations(intensityControl,nTotalSweeps)
 
!  Restore openmp max threads to what it was before teton rtmain called
!  so teton doesn't impact thread behavior of any other libraries that 
!  run after teton.
#if defined(TETON_ENABLE_OPENMP)
   call omp_set_num_threads(nThreadsInitial)
   TETON_VERIFY(omp_get_max_threads() == nThreadsInitial, "Teton: Unable to restore max threads in OpenMP runtime to pre-Teton value.")
#endif

   return
   end subroutine rtmainsn
