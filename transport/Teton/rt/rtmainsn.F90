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

   subroutine rtmainsn(dtrad, PSIR, PHI, angleLoopTime)

   use kind_mod
   use iter_control_list_mod
   use iter_control_mod
   use Size_mod
   use Material_mod
   use QuadratureList_mod
   use Quadrature_mod
   use constant_mod
   use radconstant_mod

   implicit none

!  Arguments

   real(adqt), intent(in)    :: dtrad

   real(adqt), intent(inout) :: psir(Size%ngr,Size%ncornr,Size%nangSN), &
                                Phi(Size%ngr,Size%ncornr), angleLoopTime

!  Local

   integer    :: NumSnSets

   integer    :: noutrt, ninrt, intensityIter, izero
   integer    :: nbelem, ngr, nangSN
   integer    :: set, NumQuadSets, NumBin

   real(adqt) :: maxEnergyDensityError, maxTempError 

!  Dynamic Arrays
 
!  Photon Intensities on the problem boundary

   real(adqt), allocatable :: psib(:,:,:)



!  Constants:

   parameter (izero=0)



!  Set some scalars used for dimensioning

   nbelem   = Size%nbelem
   ngr      = Size%ngr
   nangSN   = Size%nangSN

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
!                                                                      *
!     ALLOCATE MEMORY                                                  *
!                                                                      *
!***********************************************************************
 
!  Photon Intensities on the problem boundary

   allocate( psib(ngr,nbelem,nangSN) )

!***********************************************************************
!     SWEEP ORDER                                                      *
!***********************************************************************
 
!  Find reflected angles on all reflecting boundaries

   call timer_beg('reflect')
   call findReflectedAngles
   call timer_end('reflect')

!  Calculate ordering for grid sweeps

   call timer_beg('rtorder')
   call rtorder 
   call timer_end('rtorder')

!***********************************************************************
!     INITIALIZE COMMUNICATION                                         *
!***********************************************************************

!  Create an incident and exiting list for shared boundary elements

   call findexit

!***********************************************************************
!     SAVE PREVIOUS CYCLE INFORMATION AND BEGIN MATERIAL COUPLING      *
!***********************************************************************
 
!  Save various quantities from previous time step and calculate
!  the time-dependent source

   call timer_beg('rtstrtsn')
   call rtstrtsn(psir, Phi, PSIB)
   call timer_end('rtstrtsn')

!  Energy Change due to Compton scattering

   call timer_beg('compton')
   call rtcompton(Phi) 
   call timer_end('compton')

!***********************************************************************
!     EXCHANGE BOUNDARY FLUXES                                         *
!***********************************************************************

!  Establish angle order for transport sweeps

   call timer_beg('scheduler')
   call SweepScheduler
   call timer_end('scheduler')

!  Initialize Absorption Rate

   call timer_beg('absorbrate')
   call getAbsorptionRate(Phi) 
   call timer_end('absorbrate')

   call timer_beg('material')
   call UpdateMaterialCoupling(dtrad)
   call timer_end('material')

!***********************************************************************
!     BEGIN IMPLICIT ELECTRON/RADIATION COUPLING ITERATION (OUTER)     *
!***********************************************************************
 
   noutrt = 0
   ninrt  = 0
 
   TemperatureIteration: do
 
     noutrt = noutrt + 1

!***********************************************************************
!     BEGIN PHOTON INTENSITY ITERATION (INNER)                         *
!***********************************************************************

     intensityIter = 0
 
     IntensityIteration: do
 
       intensityIter = intensityIter + 1
 
!***********************************************************************
!     BEGIN LOOP OVER BATCHES                                          *
!***********************************************************************
 
       GroupSetLoop: do set=1,NumSnSets

         QuadSet => getQuadrature(Quad, set)

!  Sweep all angles in all groups in this "batch"
 
         call timer_beg('exch')
         call InitExchange
         call exchange(PSIB, izero, izero)
         call timer_end('exch')

         call timer_beg('rswpmd')
         call rswpmd(PSIB, PSIR, PHI, angleLoopTime)
         call timer_end('rswpmd')

       enddo GroupSetLoop
 
!***********************************************************************
!     END FREQUENCY GROUP LOOP                                         *
!***********************************************************************

!  Update Absorption Rate

       Mat%AbsorptionRateOld(:) = Mat%AbsorptionRate(:)

       call timer_beg('absorbrate')
       call getAbsorptionRate(Phi)
       call timer_end('absorbrate')

!***********************************************************************
!     CHECK CONVERGENCE OF SCALAR INTENSITIES                          *
!***********************************************************************

       call timer_beg('rtconi')
       call rtconi(maxEnergyDensityError, Phi)
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
     call UpdateMaterialCoupling(dtrad)
     call timer_end('material')

!  Check convergence of electron temperature
 
     call timer_beg('rtconv')
     call rtconv(maxTempError) 
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
   call bdyedt(psib)
   call timer_end('bdyedt')

!***********************************************************************
!     RELEASE MEMORY                                                   *
!***********************************************************************
 
!  Photon Intensities on thr problem boundary

   deallocate( psib )
 
!  Structures for communicating boundary fluxes and sweeps

   NumQuadSets = getNumQuadSets(Quad)
   do set=1,NumQuadSets
     QuadSet => getQuadrature(Quad, set)
     call destructComm(QuadSet)
     call destructExitList(QuadSet)
   enddo




   return
   end subroutine rtmainsn



