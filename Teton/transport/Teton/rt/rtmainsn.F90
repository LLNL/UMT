!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   RTMAINSN - Control program for SN radiation transport in 1D,       *
!              2D and 3D geometries.                                   *
!                                                                      *
!                                                                      *
!   Units:   E/e/T/m/L/A/V/t -                                         *
!        energy/photon energy/temperature/mass/length/area/volume/time *
!***********************************************************************

   subroutine rtmainsn(dtrad, PSIR, PHI, psib, angleLoopTime)

   use, intrinsic :: iso_c_binding
   use kind_mod
   use iter_control_list_mod
   use iter_control_mod
   use Size_mod
   use Material_mod
   use QuadratureList_mod
   use Quadrature_mod
   use constant_mod
   use radconstant_mod
   use cudafor
   use GPUhelper_mod

!!#include "cudaProfiler.h"

   implicit none
   

!  Arguments

   real(adqt), intent(in)    :: dtrad

   real(adqt), intent(inout) :: psib(Size%ngr,Size%nbelem,Size%nangSN), &
                                psir(Size%ngr,Size%ncornr,Size%nangSN), &
                                Phi(Size%ngr,Size%ncornr), angleLoopTime

!  Local

   integer    :: NumSnSets

   integer    :: noutrt, ninrt, intensityIter, izero

   integer    :: set, NumQuadSets, NumBin

   real(adqt) :: maxEnergyDensityError, maxTempError 

   integer :: mm1, buffer

!  Dynamic Arrays


#ifdef PROFILING_ON
   integer profiler(2) / 0, 0 /
   save profiler
#endif

!  Constants:

   parameter (izero=0)

#ifdef PROFILING_ON
   call TAU_PROFILE_TIMER(profiler, 'rtmainsn')
   call TAU_PROFILE_START(profiler)
#endif

   ! start cuda profiler
   call cudaProfilerStart()


   NumSnSets = getNumSnSets(Quad)

   Size%CommTimeCycle = zero

!  Iteration Controls

   temperatureControl  => getIterationControl(IterControls, "temperature")
   intensityControl    => getIterationControl(IterControls, "intensity")
   scatteringControl   => getIterationControl(IterControls, "scattering")
   incidentFluxControl => getIterationControl(IterControls, "incidentFlux")

!  Initialize counters for this time step

   call setNumberOfIterations(temperatureControl,izero)
   call setNumberOfIterations(intensityControl,izero)
   call setNumberOfIterations(scatteringControl,izero)
   call setControls(incidentFluxControl,maxNumberOfIterations=2)
   call setGlobalError(temperatureControl,0.1d0)


!***********************************************************************
!     SWEEP ORDER                                                      *
!***********************************************************************
 
!  Find reflected angles on all reflecting boundaries 

   call timer_beg('reflect')
   call nvtxStartRange("reflect",6)
   call findReflectedAngles
   call nvtxEndRange
   call timer_end('reflect')

!  Calculate ordering for grid sweeps

   call timer_beg('rtorder')
   call nvtxStartRange("rtorder",5)
   call rtorder 
   call nvtxEndRange
   call timer_end('rtorder')

!***********************************************************************
!     INITIALIZE COMMUNICATION                                         *
!***********************************************************************

!  Create an incident and exiting list for shared boundary elements
   call timer_beg('findexit')
   call nvtxStartRange("findexit",6)
   call findexit
   call nvtxEndRange
   call timer_end('findexit')

!  Establish angle order for transport sweeps
   ! once angle order is established, could move in psi(first_angle)

   call timer_beg('scheduler')
   call nvtxStartRange("scheduler",7)
   call SweepScheduler
   call nvtxEndRange
   call timer_end('scheduler')

!***********************************************************************
!     SAVE ZONE AVERAGE TEMPERATURES FOR TIME STEP CALCULATION         *
!*********************************************************************** 

   call timer_beg('advanceRT')
   call nvtxStartRange("advanceRT",5)
   ! calls snmoments to consume psir, produce phi
   ! scales psi (and phi too).
   call advanceRT(dtrad, PSIR, PHI, psib)
   call nvtxEndRange
   call timer_end('advanceRT')



!***********************************************************************
!     SAVE PREVIOUS CYCLE INFORMATION AND BEGIN MATERIAL COUPLING      *
!***********************************************************************
 
!  Save various quantities from previous time step and calculate
!  the time-dependent source


   call timer_beg('rtstrtsn')
   ! in: psir, phi
   ! out: psib from set boundary
   ! removed psir and psib touching routines to advanceRT.
   call rtstrtsn( Phi )
   call timer_end('rtstrtsn')

!  Energy Change due to Compton scattering

   call timer_beg('compton')
   ! in: Phi (not even used in UMT)
   ! out: nothing
   call rtcompton(Phi) 
   call timer_end('compton')

!***********************************************************************
!     EXCHANGE BOUNDARY FLUXES                                         *
!***********************************************************************


!  Initialize Absorption Rate

   call timer_beg('absorbrate')
   call nvtxStartRange("absorbrate")
   ! in: Phi
   ! out: absorbrate (1 value per corner)
   call getAbsorptionRate(Phi) 
   call nvtxEndRange
   call timer_end('absorbrate')

   call timer_beg('material')
   call nvtxStartRange("material")
   call UpdateMaterialCoupling(dtrad)
   call nvtxEndRange
   call timer_end('material')

!***********************************************************************
!     BEGIN IMPLICIT ELECTRON/RADIATION COUPLING ITERATION (OUTER)     *
!***********************************************************************

   ! ! debugging optimized code:                                                

   ! print *, "psib before starting sweeps: ", psib(1,1,1), psib(1,1,Size%nangSN)

   ! print *, "phi before starting sweeps: ", phi(1,1), phi(1,Size%ncornr)

   ! print *, "STime before starting sweeps: ", Geom%ZDataSoA%STime(1,1,1), Geom%ZDataSoA%STime(1,1,Size%nangSN)
 

   ! print *, "d_psi(1)%owner = ", d_psi(1)%owner, "d_psi(2)%owner = ", d_psi(2)%owner

   ! print *, "d_STime(1)%owner = ", d_STime(1)%owner, "d_STime(2)%owner = ", d_STime(2)%owner

   noutrt = 0
   ninrt  = 0
 
   TemperatureIteration: do
 
     noutrt = noutrt + 1

!***********************************************************************
!     BEGIN PHOTON INTENSITY ITERATION (INNER)                         *
!***********************************************************************

! it looks like only 1 inner per outer for SuOlson, i.e. this loop is not a loop.

     intensityIter = 0
 
     IntensityIteration: do
 
       intensityIter = intensityIter + 1
 
!***********************************************************************
!     BEGIN LOOP OVER BATCHES (meaning batches of quadsets here)       *
!***********************************************************************
 
       GroupSetLoop: do set=1,NumSnSets

         QuadSet => getQuadrature(Quad, set)

!  Sweep all angles in all groups in this "batch"

 
         ! CPU code should wait until psib is on the host before exchanging.
         istat=cudaEventSynchronize( psib_OnHost( current%batch ) )

         call timer_beg('exch')
         call nvtxStartRange("exch all bins")
         call InitExchange
         ! exchange over all angle bins. Maybe could be done a bin at at time, overlapped with advanceRT.
         call exchange(PSIB, izero, izero) 
         call nvtxEndRange
         call timer_end('exch')

         call timer_beg('rswpmd')
         call nvtxStartRange("rswpmd",2)
         call rswpmd(PSIB, PSIR, PHI, angleLoopTime, intensityIter, noutrt)
         call nvtxEndRange
         call timer_end('rswpmd')

       enddo GroupSetLoop
 
!***********************************************************************
!     END FREQUENCY GROUP LOOP                                         *
!***********************************************************************

!  Update Absorption Rate

       Mat%AbsorptionRateOld(:) = Mat%AbsorptionRate(:)

       call timer_beg('absorbrate')
       call nvtxStartRange("getAbsorptionRate")
       ! in: Phi
       ! out: absorbrate (1 value per corner)
       call getAbsorptionRate(Phi)
       call nvtxEndRange
       call timer_end('absorbrate')

!***********************************************************************
!     CHECK CONVERGENCE OF SCALAR INTENSITIES                          *
!***********************************************************************

       call timer_beg('rtconi')
       call nvtxStartRange("rtconi")
       call rtconi(maxEnergyDensityError, Phi)
       call nvtxEndRange
       call timer_end('rtconi')
 
       if (maxEnergyDensityError < getEpsilonPoint(intensityControl) .or. &
           intensityIter >= getMaxNumberOfIterations(intensityControl)) then
         exit IntensityIteration
       else
         cycle IntensityIteration
       endif
 
     enddo IntensityIteration

     ninrt = ninrt + intensityIter
 
!***********************************************************************
!     END PHOTON INTENSITY ITERATION (INNER)                           *
!***********************************************************************
 
!  Calculate new electron temperature and energy change
 
     call timer_beg('material')
     call nvtxStartRange("material")
     call UpdateMaterialCoupling(dtrad)
     call nvtxEndRange
     call timer_end('material')

!  Check convergence of electron temperature
 
     call timer_beg('rtconv')
     call nvtxStartRange("rtconv")
     call rtconv(maxTempError) 
     call nvtxEndRange
     call timer_end('rtconv')

     if ((maxTempError <  getEpsilonPoint(temperatureControl) .and.  &
          maxEnergyDensityError <  getEpsilonPoint(intensityControl))  .or.   &
          noutrt >= getMaxNumberOfIterations(temperatureControl)) then

       exit TemperatureIteration

     else

       if (maxTempError <  getEpsilonPoint(temperatureControl)) then
         call setControls(incidentFluxControl,maxNumberOfIterations=4)
       else
         call setControls(incidentFluxControl,maxNumberOfIterations=3)
       endif

       cycle TemperatureIteration

     endif
 
   enddo TemperatureIteration

!  Update Iteration Counts

   call setNumberOfIterations(temperatureControl,noutrt)
   call setNumberOfIterations(intensityControl,ninrt)
 
!***********************************************************************
!     END IMPLICIT ELECTRON/RADIATION COUPLING ITERATION (OUTER)       *
!***********************************************************************
 

!***********************************************************************
!     BOUNDARY EDITS                                                   *
!***********************************************************************

   call timer_beg('bdyedt')
   call nvtxStartRange("bdyedt")
   ! in: psib
   call bdyedt(psib)
   call nvtxEndRange
   call timer_end('bdyedt')

!***********************************************************************
!     RELEASE MEMORY                                                   *
!***********************************************************************
 
!  Photon Intensities on the problem boundary

   !! These are saved and reused across calls
   !deallocate( psib )
   !istat = cudaHostUnregister(C_LOC(psir(1,1,1)))
   !istat = cudaHostUnregister(C_LOC(phi(1,1)))
 
!  Structures for communicating boundary fluxes and sweeps

   NumQuadSets = getNumQuadSets(Quad)
   do set=1,NumQuadSets
     QuadSet => getQuadrature(Quad, set)
     call destructComm(QuadSet)
     call destructExitList(QuadSet)
   enddo


   call cudaProfilerStop()

#ifdef PROFILING_ON
   call TAU_PROFILE_STOP(profiler)
#endif

   return
   end subroutine rtmainsn



