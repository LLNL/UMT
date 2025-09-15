#include "macros.h"
!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   DTNEW - Calculates a new time step based on the maximum changes in *
!           the electron and radiation temperatures.                   *
!                                                                      *
!***********************************************************************
   subroutine dtnew(maxOSComptonChangeCorner,  &
                    maxOSComptonChange) BIND(C,NAME="teton_dtnew")


   USE ISO_C_BINDING
   use flags_mod
   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use constant_mod
   use iter_control_list_mod
   use iter_control_mod
   use TimeStepControls_mod
   use Size_mod
   use Geometry_mod
   use Material_mod

   use Datastore_mod, only: theDatastore

   implicit none

!  Arguments

   integer(C_INT), intent(in) :: maxOSComptonChangeCorner
   real(C_DOUBLE), intent(in) :: maxOSComptonChange

!  Local

   integer, dimension (1) :: indexDtRec

   integer                :: zone
   integer                :: zoneMaxChangeTe
   integer                :: zoneMaxChangeEr
   integer                :: zoneOSCompton
   integer                :: zoneConvControl
   integer                :: numTempIterations 
   integer                :: ZoneConst(2) 
   integer                :: Rank
   integer                :: myRankInGroup 
   integer                :: nsend 

   ! time step votes are stored accoridng to these indices
   integer, parameter     :: indexEr=1
   integer, parameter     :: indexTe=2
   integer, parameter     :: indexCompton=3
   integer, parameter     :: indexSlowConv=4
   integer, parameter     :: indexNoConv=5

   ! this array is like a dictionary to return the correct flag
   ! it must match the same order as the parameters above
   integer, dimension(5) :: dtConstraint = &
        (/ dtControl_radTemp , &
           dtControl_elecTemp, &
           dtControl_compton, &
           dtControl_slowConv, &
           dtControl_noConv /)

   real(adqt)       :: delta_te, delta_Er, deltaTeMax, deltaErMax
   real(adqt)       :: Er_offset
   real(adqt)       :: facTe, facEr
   real(adqt)       :: tez, tezn, tr, trn, Er, Ern

   ! User-settable controls with default values:
   real(adqt)       :: maxChangeTe, maxChangeEr      ! maximum allowed change in Te and Er ("delte" and "deltr" in inputs)
   real(adqt)       :: Te_cutoff                     ! T_e values below this temperautre do not participate in the dt voting
   real(adqt)       :: Er_offset_frac                ! The denominator offset for the radiation energy dt vote is Er_offset_frac*(maximum zonal radiation energy)
   real(adqt)       :: dtReductionLimit              ! Votes for T_r and T_e cannot be smaller than this fraction of the old time step
   real(adqt)       :: dtIncreaseLimit_t             ! Maximum growth for T_r and T_e based votes
   real(adqt)       :: dtIncreaseLimit_i             ! Maximum growth for iteration count based votes
   integer          :: outerSlowConvergenceThreshold ! Above this value, dtSlowConv = dtold
   integer          :: outerNoConvergenceThreshold   ! Above this value, dtNoConv = dtold*dtNoConvReductionFactor
   real(adqt)       :: dtNoConvReductionFactor

   character(len=20) :: temp_str ! for printouts in VERIFY

   real(adqt)       :: dtMin, dtMax, dtRec, dtRad, my_dtRad 

   real(adqt)       :: dtRecList(5)

   type(IterControl) , pointer :: temperatureControl => NULL()
   
   ! votes for dt flags
   integer :: constraint = dtControl_invalid  ! default  (overwritten or this routine is broken)

   myRankInGroup = Size% myRankInGroup 
   !  Iteration Control
   !  Time step is only controlled by temperature (outer) iteration count
   temperatureControl => getIterationControl(IterControls,"temperature")

   ! Defaults from before the refactor:
   maxChangeTe                   = getMaxChangeTe(DtControls) ! currently defined by iteration/delte
   maxChangeEr                   = getMaxChangeEr(DtControls) ! currently defined by iteration/deltr
   ! ^ You can change these in the middle of the run by changing iteration/dtcontrols/delte or iteration/dtcontrols/deltr
   Te_cutoff                     = 0.008d0
   Er_offset_frac                = 0.01
   dtReductionLimit              = half
   dtIncreaseLimit_t             = one+ninth ! i.e., 1.11111...
   dtIncreaseLimit_i             = two
   outerSlowConvergenceThreshold = 16_C_INT
   outerNoConvergenceThreshold   = MIN(20_C_INT, getMaxNumberOfIterations(temperatureControl))
   dtNoConvReductionFactor       = 0.8_adqt

   ! Optional inputs from datastore:
   Te_cutoff                     = theDatastore%fetchIfExists("options/iteration/dtcontrols/delte_cutoff", Te_cutoff)
   Er_offset_frac                = theDatastore%fetchIfExists("options/iteration/dtcontrols/deltr_offset_frac", Er_offset_frac)
   maxChangeTe                   = theDatastore%fetchIfExists("options/iteration/dtcontrols/delte", maxChangeTe)
   maxChangeEr                   = theDatastore%fetchIfExists("options/iteration/dtcontrols/deltr", maxChangeEr)
   dtReductionLimit              = theDatastore%fetchIfExists("options/iteration/dtcontrols/dt_reduction_limit", dtReductionLimit)
   dtIncreaseLimit_t             = theDatastore%fetchIfExists("options/iteration/dtcontrols/dt_increase_limit_t", dtIncreaseLimit_t)
   dtIncreaseLimit_i             = theDatastore%fetchIfExists("options/iteration/dtcontrols/dt_increase_limit_i", dtIncreaseLimit_i)
   outerSlowConvergenceThreshold = theDatastore%fetchIfExists("options/iteration/dtcontrols/outer_slow_convergence_threshold", outerSlowConvergenceThreshold)
   outerNoConvergenceThreshold   = theDatastore%fetchIfExists("options/iteration/dtcontrols/outer_no_convergence_threshold", outerNoConvergenceThreshold)
   dtNoConvReductionFactor       = theDatastore%fetchIfExists("options/iteration/dtcontrols/dt_noconv_reduction_factor", dtNoConvReductionFactor)

   WRITE(temp_str, '(F8.4)') Te_cutoff
   TETON_VERIFY(Te_cutoff > zero, "Electron temperature cutoff for time step vote must be positive. Bad options/iteration/dtcontrols/delte_cutoff value:"//trim(adjustl(temp_str)))
   WRITE(temp_str, '(F8.4)') Er_offset_frac
   TETON_VERIFY(Er_offset_frac > zero .and. Er_offset_frac < 1.000001_adqt, "options/iteration/dtcontrols/deltr_offset_frac must be > 0 and <= 1, current value:"//trim(adjustl(temp_str)))
   WRITE(temp_str, '(F8.4)') maxChangeTe
   TETON_VERIFY(maxChangeTe > zero, "options/iteration/dtcontrols/delte must be positive, current value:"//trim(adjustl(temp_str)))
   WRITE(temp_str, '(F8.4)') maxChangeEr
   TETON_VERIFY(maxChangeEr > zero, "options/iteration/dtcontrols/deltr must be positive, current value:"//trim(adjustl(temp_str)))
   WRITE(temp_str, '(F8.4)') dtReductionLimit
   TETON_VERIFY(dtReductionLimit > zero .and. dtReductionLimit < 1.000001_adqt, "options/iteration/dtcontrols/dt_reduction_limit must be > 0 and <= 1, current value:"//trim(adjustl(temp_str)))
   WRITE(temp_str, '(F8.4)') dtIncreaseLimit_t
   TETON_VERIFY(dtIncreaseLimit_t > 0.999999_adqt, "options/iteration/dtcontrols/dt_increase_limit_t must be >= 1, current value:"//trim(adjustl(temp_str)))
   WRITE(temp_str, '(F8.4)') dtIncreaseLimit_i
   TETON_VERIFY(dtIncreaseLimit_i > 0.999999_adqt, "options/iteration/dtcontrols/dt_increase_limit_i must be >= 1, current value:"//trim(adjustl(temp_str)))
   WRITE(temp_str, '(F8.4)') dtNoConvReductionFactor
   TETON_VERIFY(dtNoConvReductionFactor > zero .and. dtNoConvReductionFactor < 1.000001_adqt, "options/iteration/dtcontrols/dt_noconv_reduction_factor must be > 0 and <= 1, current value:"//trim(adjustl(temp_str)))

!  offset for time step control for E_r is 1% of maximum zonal radiation energy
   tr = Size%tfloor
   Er_offset = tr*tr*tr*tr*minval(Geom%VolumeZone)
   do zone=1,Size%nzones
     tr        = Mat%Trz(zone)
     Er        = tr*tr*tr*tr*Geom%VolumeZone(zone)
     Er_offset = max(Er_offset, Er)
   enddo
   Er_offset = Er_offset_frac*Er_offset

   call MPIAllReduce(Er_offset, "max", MY_COMM_GROUP)

!  We test on the radiation energy density and electron temperature

   zoneMaxChangeTe  = 1 
   zoneMaxChangeEr  = 1
   deltaTeMax       = zero
   deltaErMax       = zero

   ZoneLoop: do zone=1,Size%nzones

     if (.not. Mat% isVoid(zone)) then

       tr   = Mat%trz(zone)
       trn  = Mat%trzn(zone)
       Er   = tr*tr*tr*tr*Geom%VolumeZone(zone)
       Ern  = trn*trn*trn*trn*Geom%VolumeZone(zone)

       delta_Er = abs(Er - Ern)/(Ern + Er_offset)

       if (delta_Er > deltaErMax) then
         zoneMaxChangeEr = zone
         deltaErMax = delta_Er
       endif

       tez  = Mat%tez(zone)
       tezn = Mat%tezn(zone)

       if (tez > Te_cutoff .and. tezn > Te_cutoff) then
          delta_te  = abs(tez - tezn)/tezn

          if (delta_te > deltaTeMax) then
            zoneMaxChangeTe = zone
            deltaTeMax = delta_te
          endif
       endif

     endif

   enddo ZoneLoop

!  What is the controlling time step

   constraint = dtControl_none

!  Temperature-based time step vote can decrease by only a factor 2 per cycle

   dtRad  = getRadTimeStep(DtControls)

   facTe  = min(one/dtReductionLimit,deltaTeMax/maxChangeTe)
   facEr = min(one/dtReductionLimit,deltaErMax/maxChangeEr)

!  Temperature-based time step vote can increase by only a factor dtIncreaseLimit_t per cycle

   dtRecList(indexEr)  = dtRad/max(one/dtIncreaseLimit_t, facEr)
   dtRecList(indexTe)  = dtRad/max(one/dtIncreaseLimit_t, facTe )

!  Operator-split Compton

   dtRecList(indexCompton) = one/(maxOSComptonChange + half)*dtRad
   zoneOSCompton           = Geom% CToZone(maxOSComptonChangeCorner)

!  If the iteration count is approaching the maximum allowed,
!  do not increase the time step further. 
   
   numTempIterations           = getNumberOfIterations(temperatureControl)
   zoneConvControl             = getZoneOfMax(temperatureControl)

   if (numTempIterations >= outerSlowConvergenceThreshold) then
     dtRecList(indexSlowConv) = dtRad
   else
     dtRecList(indexSlowConv) = dtIncreaseLimit_i*dtRad
   endif

   if (numTempIterations >= outerNoConvergenceThreshold) then
     dtRecList(indexNoConv) = dtNoConvReductionFactor*dtRad
   else
     dtRecList(indexNoConv) = dtIncreaseLimit_i*dtRad
   endif

   dtMax = getMaxTimeStep(DtControls)
   dtMin = getMinTimeStep(DtControls)

   ! which of the 5 possible voting indices limits us
   ! minloc returns a vector of length 1 (since dtRecList is a 1-D array)
   indexDtRec = minloc( dtRecList(:) )
   dtRec      = dtRecList( indexDtRec(1) )
   constraint = dtConstraint ( indexDtRec(1) )


   dtRad  = max(dtRec,dtMin)
   dtRad  = min(dtRad,dtMax)

!  Choose the minimum time step over all domains

   my_dtRad = dtRad

   call MPIAllReduce(DTRAD, "min", MY_COMM_GROUP)

!  Select the process controlling the time step first (in the event of a tie),
!  and use it to provide the constraint and zone

   Rank = -1
   if (my_dtRad == dtRad) then
     Rank = myRankInGroup
   endif

   call MPIAllReduce(Rank, "max", MY_COMM_GROUP)

!  For the controlling process, broadcast the constraint and zone

   if (myRankInGroup == Rank) then

     if (constraint == dtControl_elecTemp) then
       ZoneConst(1) = zoneMaxChangeTe
     elseif (constraint == dtControl_radTemp) then
       ZoneConst(1) = zoneMaxChangeEr
     elseif (constraint == dtControl_slowConv .or. &
             constraint == dtControl_noConv) then
       ZoneConst(1) = zoneConvControl
     elseif (constraint == dtControl_compton) then
       ZoneConst(1) = zoneOSCompton
     endif

     ZoneConst(2) = constraint

   endif
                                                                                       
   nsend = 2
   call MPIBcast(ZoneConst, nsend, Rank, MY_COMM_GROUP)

   ! The rank and zone are already synced in rt/ConvergenceTest.F90 in the case
   !   of slow/no convergence controlling dt.  If reason is slow/no conv.,
   !   rank has the highest processor ID stored, which generally does not
   !   correspond to the process of the zone in ZoneConst(1)/zoneConvControl.
   if (ZoneConst(2) == dtControl_slowConv .or. &
       ZoneConst(2) == dtControl_noConv) then
      Rank = getProcessOfMax(temperatureControl)
   endif

!  Update controls

   call setDtControls(DtControls,                        &
                      ControlProcess=Rank,               &
                      ControlZone=ZoneConst(1),          &
                      ZoneMaxChangeEr=zoneMaxChangeEr,   &
                      ZoneMaxChangeTe=zoneMaxChangeTe,   &
                      RecTimeStep=dtRad,                 &
                      MaxFracChangeEr=deltaErMax,        &
                      MaxFracChangeTe=deltaTeMax,        &
                      dtConstraint=ZoneConst(2) )


   return
   end subroutine dtnew


