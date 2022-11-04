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
   use default_iter_controls_mod, only : outer_slow_conv_threshold

   implicit none

!  Arguments

   integer(C_INT), intent(in) :: maxOSComptonChangeCorner
   real(C_DOUBLE), intent(in) :: maxOSComptonChange

!  Local

   integer, dimension (1) :: indexDtRec

   integer                :: zone
   integer                :: zoneMaxChangeTe
   integer                :: zoneMaxChangeTr4
   integer                :: zoneOSCompton
   integer                :: zoneConvControl
   integer                :: numTempIterations 
   integer                :: maxTempIterations
   integer                :: ZoneConst(2) 
   integer                :: Rank
   integer                :: myRankInGroup 
   integer                :: nsend 

   ! time step votes are stored accoridng to these indices
   integer, parameter     :: indexTr4=1
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

   real(adqt)       :: delta_te, delta_tr4, deltaTeMax, deltaTr4Max
   real(adqt)       :: facTe, facTr4
   real(adqt)       :: tmin, tez, tezn, tr, trn, tr4, tr4n, TrMax, Tr4Max
   real(adqt)       :: maxChangeTe, maxChangeTr4
   real(adqt)       :: threshold

   real(adqt)       :: dtMin, dtMax, dtRec, dtRad, my_dtRad 

   real(adqt)       :: dtRecList(5)

   type(IterControl) , pointer :: temperatureControl => NULL()
   
   ! votes for dt flags
   integer :: constraint = dtControl_invalid  ! default  (overwritten or this routine is broken)

!  Constants
! TODO: is this tmin an undocumented floor?
   myRankInGroup = Size% myRankInGroup 
   tmin        = 0.008d0

   !  Iteration Control
   !  Time step is only controlled by temperature iteration
   temperatureControl => getIterationControl(IterControls,"temperature")

!  We test on the radiation energy density and electron temperature

   zoneMaxChangeTe  = 1 
   zoneMaxChangeTr4 = 1 
   deltaTeMax       = zero
   deltaTr4Max      = zero
   TrMax            = zero
   Tr4Max           = zero

   do zone=1,Size%nzones
     Tr = Mat%Trz(zone)
     if (Tr > tmin) then
       Tr4    = Tr*Tr*Tr*Tr*Geom% VolumeZone(zone)
       Tr4Max = max(Tr4, Tr4Max)
     endif
   enddo

   threshold = cutoff*Tr4Max

   call MPIAllReduce(THRESHOLD, "max", MY_COMM_GROUP)

   ZoneLoop: do zone=1,Size%nzones

     tr   = Mat%trz(zone)
     trn  = Mat%trzn(zone)
     tr4  = tr*tr*tr*tr
     tr4n = trn*trn*trn*trn 

     tez  = Mat%tez(zone)
     tezn = Mat%tezn(zone)

     if (Mat% isVoid(zone)) then

       delta_tr4 = zero
       delta_te  = zero

     else

       delta_tr4 = abs(tr4 - tr4n)/tr4n

       if (tr4*Geom% VolumeZone(zone) > threshold) then
         delta_tr4 = abs(tr4 - tr4n)/tr4n

         if (delta_tr4 > deltaTr4Max) then
           zoneMaxChangeTr4 = zone 
           deltaTr4Max = delta_tr4
         endif
       endif

       if (tez > tmin .and. tezn > tmin) then
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

!  Time step can decrease by only a factor 2 per cycle

   maxChangeTe   = getMaxChangeTe(DtControls)
   maxChangeTr4  = getMaxChangeTr4(DtControls)
   dtRad         = getRadTimeStep(DtControls)

   facTe  = min(two,deltaTeMax/maxChangeTe)
   facTr4 = min(two,deltaTr4Max/maxChangeTr4)

!  Time step can increase by only a factor 1/TempFraction per cycle

   dtRecList(indexTr4) = dtRad/max(TempFraction,facTr4)
   dtRecList(indexTe)  = dtRad/max(TempFraction,facTe)

!  Operator-split Compton

   dtRecList(indexCompton) = one/(maxOSComptonChange + half)*dtRad
   zoneOSCompton           = Geom% CToZone(maxOSComptonChangeCorner)

!  If the iteration count is approaching the maximum allowed,
!  do not increase the time step further. 

!  20+ iterations is considered slow even if you're allowing for many more outers
   maxTempIterations = MIN(outer_slow_conv_threshold,getMaxNumberOfIterations(temperatureControl))
   numTempIterations = getNumberOfIterations(temperatureControl)
   zoneConvControl   = getZoneOfMax(temperatureControl)

   if (numTempIterations >= IterFraction*maxTempIterations) then
     dtRecList(indexSlowConv) = dtRad
   else
     dtRecList(indexSlowConv) = two*dtRad
   endif

!  If the iteration did not converge cut the timestep by IterFraction
!    We may want to consider cutting it in half instead!

   if (numTempIterations >= maxTempIterations) then
     dtRecList(indexNoConv) = IterFraction*dtRad
   else
     dtRecList(indexNoConv) = two*dtRad
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
       ZoneConst(1) = zoneMaxChangeTr4
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
                      ZoneMaxChangeTr4=zoneMaxChangeTr4, &
                      ZoneMaxChangeTe=zoneMaxChangeTe,   &
                      RecTimeStep=dtRad,                 &
                      MaxFracChangeTr4=deltaTr4Max,      &
                      MaxFracChangeTe=deltaTeMax,        &
                      Tr4Threshold=threshold,            &
                      dtConstraint=ZoneConst(2) )


   return
   end subroutine dtnew


