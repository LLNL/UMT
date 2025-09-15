!***********************************************************************
!                     Version 0: 01/2005 PFN                           *
!                                                                      *
!    GETRUNSTATS -  Called from host to get information for            *
!                   convergence and time step control.                 *
!                                                                      *
!***********************************************************************

   subroutine getRunStats(MatCoupTimeTotal, SweepTimeTotal,      &
                          GPUSweepTimeTotal, GTATimeTotal,       &
                          RadtrTimeTotal, InitTimeTotal,         &
                          FinalTimeTotal, timeNonRad, timeOther) &
                          BIND(C,NAME="teton_getrunstats")
   
   use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                             stdout=>output_unit, &
                                             stderr=>error_unit
   use, intrinsic :: iso_c_binding, only : c_long, c_double, c_int

   use flags_mod
   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use constant_mod
   use Material_mod
   use Geometry_mod
   use Quadrature_mod, only : Quadrature
   use QuadratureList_mod
   use iter_control_list_mod
   use iter_control_mod
   use TimeStepControls_mod
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   use ComptonControl_mod
#endif
   use Options_mod
   use Datastore_mod, only : theDatastore
   use OMPUtilities_mod
   use system_info_mod
#if defined(TETON_ENABLE_UMPIRE)
   use umpire_mod, only : get_process_memory_usage, get_process_memory_usage_hwm
#endif

   implicit none

!  Arguments

   real(C_DOUBLE) :: MatCoupTimeTotal
   real(C_DOUBLE) :: SweepTimeTotal
   real(C_DOUBLE) :: GPUSweepTimeTotal
   real(C_DOUBLE) :: GTATimeTotal
   real(C_DOUBLE) :: RadtrTimeTotal
   real(C_DOUBLE) :: InitTimeTotal
   real(C_DOUBLE) :: FinalTimeTotal
   real(C_DOUBLE) :: timeNonRad
   real(C_DOUBLE) :: timeOther

!  Local

   real(C_DOUBLE) :: timeNonRadCycle
   real(C_DOUBLE) :: timeNonRadTotal

   real(adqt)     :: ConvControlError
   real(adqt)     :: errorTemp
   real(adqt)     :: errorPsi
   real(adqt)     :: dtrad
   real(adqt)     :: currentDtRad
   real(adqt)     :: DtControlChange
   real(adqt)     :: ConvState(5+Size% ndim)
   real(adqt)     :: DtState(8+Size% ndim)
   real(adqt)     :: zoneCenter(Size% ndim)
   real(adqt)     :: throughputCycle

   integer        :: ncycle
   integer        :: ConvControlProcess
   integer        :: ConvControlZone
   integer        :: DtControlProcess
   integer        :: DtControlZone
   integer        :: myRankInGroup 
   integer        :: ConvControlReason
   integer        :: DtConstraint
   integer        :: indexDt
   integer        :: indexCaveat
   integer        :: i
   integer        :: nsend
   integer        :: nsendDt
   integer        :: ndim

   integer (kind=c_long) :: globalNumUnknowns

   ! Variables affecting threading
   integer        :: numOmpCPUThreads
   integer        :: nZoneSets
   integer        :: nSets
   integer        :: nGTASets
   integer        :: nSweepHyperDomains
   integer        :: nGreySweepHyperDomains
   integer        :: sweepVersion
   integer(c_int) :: nGPUProcessors

   character(len=26), parameter :: Tformat1 = "(1X,A16,1X,F14.6,3X,F14.6)" 
   character(len=48), parameter :: Tformat2 = "(1X,A16,1X,F14.6,3X,F14.6,5X,F5.1,A1)" 
   character(len=14), parameter :: Sformat = "(A21,1pe18.11)" 
   character(len=13), parameter :: format1D = "(A20,1pe12.4)"
   character(len=21), parameter :: format2D = "(A20,1pe12.4,1pe12.4)"
   character(len=29), parameter :: format3D = "(A20,1pe12.4,1pe12.4,1pe12.4)"

   character(len=19), dimension(5) :: DtControlReasonC =  &
                                     (/" Rad Energy Density",  &
                                       " Electron Temp     ",  &
                                       " Compton Scattering",  &
                                       " Slow Convergence  ",  &
                                       " No Convergence    " /)

   character(len=21), dimension(2) :: DtCaveat = &
                                     (/" At Minimum time step", &
                                       " At Maximum time step" /)

   character(len=57), dimension(5) :: messageRoot = &
                     (/"Radiation Energy Density change is                       ", &
                       "Electron Temperature change is                           ", &
                       "Operator-split Compton change is                         ", &
                       "Solver converged, approaching maximum iterations allowed ", &
                       "Solver did not converge within maximum iterations allowed" /)
        
   character(len=11)  :: ConvControlReasonC
   character(len=13)  :: changeStr
   character(len=8)   :: zoneStr
   character(len=7)   :: procStr
   character(len=7)   :: cycleStr
   character(len=36)  :: coordStr
   character(len=160) :: messageStr

   type(IterControl) , pointer :: temperatureControl => NULL() 
   type(IterControl) , pointer :: intensityControl   => NULL()
   
!  Threading information
   sweepVersion              = Options%getSweepVersion()
   numOmpCPUThreads          = Options%getNumOmpMaxThreads()
   nZoneSets                 = getNumberOfZoneSets(Quad)
   nSets                     = getNumberOfSets(Quad)
   nGTASets                  = getNumberOfGTASets(Quad)
   nSweepHyperDomains        = getNumberOfHyperDomains(Quad,1) 
   nGreySweepHyperDomains    = getNumberOfHyperDomains(Quad,2) 

!  Iteration Controls

   temperatureControl => getIterationControl(IterControls,"temperature")
   intensityControl   => getIterationControl(IterControls,"intensity")

   myRankInGroup = Size% myRankInGroup 
   ndim          = Size% ndim
   nsend         = 5 + ndim
   nsendDt       = 8 + ndim
   ConvState(:)  = zero
   DtState(:)    = zero

!  Iteration Statistics
   errorTemp = getGlobalError(temperatureControl)
   errorPsi  = getGlobalError(intensityControl)

   if (errorPsi >= errorTemp) then
     ConvControlProcess = getProcessOfMax(intensityControl)
     ConvControlReason  = convControl_fluxIter
     ConvControlZone    = getZoneOfMax(intensityControl)
     ConvControlError   = errorPsi
     ConvControlReasonC = " Intensity "

     if (myRankInGroup == ConvControlProcess) then
       ConvState(1)     = Mat%trz(ConvControlZone)
       ConvState(2)     = Mat%tez(ConvControlZone)
       ConvState(3)     = Mat%rho(ConvControlZone)
       ConvState(4)     = Mat%cve(ConvControlZone)
       ConvState(5)     = Mat%SMatEff(ConvControlZone)

       zoneCenter(:)    = getZoneCenter(Geom, ConvControlZone)

       do i=1,ndim
         ConvState(5+i) = zoneCenter(i)
       enddo

     endif
   else
     ConvControlProcess = getProcessOfMax(temperatureControl)
     ConvControlReason  = convControl_tempIter
     ConvControlZone    = getZoneOfMax(temperatureControl)
     ConvControlError   = errorTemp
     ConvControlReasonC = "Temperature"

     if (myRankInGroup == ConvControlProcess) then
       ConvState(1)     = Mat%trz(ConvControlZone)
       ConvState(2)     = Mat%tez(ConvControlZone)
       ConvState(3)     = Mat%rho(ConvControlZone)
       ConvState(4)     = Mat%cve(ConvControlZone)
       ConvState(5)     = Mat%SMatEff(ConvControlZone)

       zoneCenter(:)    = getZoneCenter(Geom, ConvControlZone)

       do i=1,ndim
         ConvState(5+i) = zoneCenter(i)
       enddo

     endif
   endif

   flush(stdout)

   call MPIBcast(ConvState, nsend, ConvControlProcess, MY_COMM_GROUP)

!  Time Step Statistics 

   DtControlProcess = getControlProcess(DtControls)
   DtControlZone    = getControlZone(DtControls)
   DtConstraint     = getDtConstraint(DtControls)
   dtrad            = getRecTimeStep(DtControls)
   currentDtRad     = getRadTimeStep(DtControls)

   indexCaveat = 0 

   if (dtrad == getMinTimeStep(DtControls)) then
     indexCaveat = 1 
   elseif (dtrad == getMaxTimeStep(DtControls)) then
     indexCaveat = 2 
   endif

   if (myRankInGroup == DtControlProcess) then

     if (DtConstraint == dtControl_radTemp) then
       DtControlChange  = getMaxFracChangeEr(DtControls)
     elseif (DtConstraint == dtControl_elecTemp) then
       DtControlChange  = getMaxFracChangeTe(DtControls)
     else
       DtControlChange  = zero
     endif

     zoneCenter(:) = getZoneCenter(Geom, DtControlZone)

     DtState(1) = DtControlChange
     DtState(2) = Mat%trz(DtControlZone)
     DtState(3) = Mat%trzn(DtControlZone)
     DtState(4) = Mat%tez(DtControlZone)
     DtState(5) = Mat%tezn(DtControlZone)
     DtState(6) = Mat%rho(DtControlZone)
     DtState(7) = Mat%cve(DtControlZone)
     DtState(8) = Mat%SMatEff(DtControlZone)

     do i=1,ndim
       DtState(8+i) = zoneCenter(i)
     enddo

   endif

   call MPIBcast(DtState, nsendDt, DtControlProcess, MY_COMM_GROUP)

   DtControlChange = DtState(1)
   indexDt         = DtConstraint - 20

!  Timings

   MatCoupTimeTotal  = Size% MatCoupTimeTotal
   SweepTimeTotal    = Size% SweepTimeTotal
   GPUSweepTimeTotal = Size% GPUSweepTimeTotal
   GTATimeTotal      = Size% GTATimeTotal
   RadtrTimeTotal    = max(Size% RadtrTimeTotal,adqtSmall)
   InitTimeTotal     = Size% InitTimeTotal
   FinalTimeTotal    = Size% FinalTimeTotal

   timeNonRadCycle   = timeNonRad/sixty
   timeNonRadTotal   = timeOther/sixty

   if ( Options%isRankVerbose() > 0 ) then

     ncycle = getRadCycle(DtControls) 

     print *,"************     Configuration Info    *************"
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     if (Size%useNewGTASolver) then
       print *, " GTA solver version: 2"
     else
       print *, " GTA solver version: 1"
     endif

     if (getUseBoltzmann(Compton)) then
       print *, " Compton scattering kernel: Boltzmann"
     endif

     if (getUseFokkerPlanck(Compton)) then
       print *, " Compton scattering kernel: Fokker Planck"
     endif
#endif

     if (sweepVersion == 1) then
       print *, " Sweep scheduled over zones"
     elseif (sweepVersion == 2) then
       print *, " Sweep scheduled over corners"
     endif

     if (Size%useGPU) then
       print *," Device : GPU"
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
       print *, " # GPU processor units = ", Options%getNumDeviceProcessors()
       print *, " # GPU thread teams utilized by zone sets = ", nZoneSets
       print *, " # GPU thread teams utilized by sweep  = ", nSets*nSweepHyperDomains
       print *, " # GPU thread teams utilized by grey sweep = ", nGTASets*nGreySweepHyperDomains
       ! Note: Umpire does not yet support querying the device memory usage for ROCM, only CUDA.
       ! I created an issue requesting this on the Umpire github site and provided code for ROCM
       ! support. (This went in April 29 2025 in Umpire)
       ! https://github.com/LLNL/Umpire/issues/959
       ! In general I like seeing both the used and free memory so am using our own
       ! internal implementation.  -- A. Black
       call printGPUMemInfo(Size%myRankInGroup)
#endif
     else
       print *," Device : CPU"
     endif

#if defined(TETON_ENABLE_UMPIRE)
     print *, " Process memory usage (CPU) = ", get_process_memory_usage() / 1024 / 1024, " MB"
     print *, " Process memory usage high water mark (CPU) = ", get_process_memory_usage_hwm() / 1024 / 1024, " MB"
#endif

#if defined(TETON_ENABLE_OPENMP)
     print '(A30,1X,I3)', "  # CPU threads per mpi rank = ", numOmpCPUThreads
! Indent the next print a couple spaces to match this block of text
     write(*, fmt="(a)", advance="no") "  "
     call print_thread_bindings()
#endif
     flush(stdout)

     print *," "
     print *,"******************     Run Time (minutes) ******************"
     print *,"                          Cycle      Accumulated  % of RADTR"
     print Tformat2, "RADTR          =", Size% RadtrTimeCycle,    RadtrTimeTotal,    100.0, "%"
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     print Tformat2, "Mat. Coupling  =", Size% MatCoupTimeCycle,  MatCoupTimeTotal,  MatCoupTimeTotal/RadtrTimeTotal*100.0, "%"
#endif
     print Tformat2, "Sweep(CPU)     =", Size% SweepTimeCycle,    SweepTimeTotal,    SweepTimeTotal/RadtrTimeTotal*100.0, "%"
     print Tformat2, "Sweep(GPU)     =", Size% GPUSweepTimeCycle, GPUSweepTimeTotal, GPUSweepTimeTotal/RadtrTimeTotal*100.0, "%"
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     print Tformat2, "Grey Tr. Accel =", Size% GTATimeCycle,      GTATimeTotal,      GTATimeTotal/RadtrTimeTotal*100.00, "%"
#endif
     print Tformat2, "Initialization =", Size% InitTimeCycle,     InitTimeTotal,     InitTimeTotal/RadtrTimeTotal*100.0, "%"
     print Tformat2, "Finalization   =", Size% FinalTimeCycle,    FinalTimeTotal,    FinalTimeTotal/RadtrTimeTotal*100.0, "%"
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     print Tformat1, "Non-Rad        =",       timeNonRadCycle,   timeNonRadTotal
#endif

     print *," "
! Print the 'throughput' for use in performance measurements.
! Throughput is defined as # unknowns solved per second in most performance
! reports.  In Teton's case for the cycle throughput I'm using:
! the # elements in PSI * timestep / RADTR walltime.
! -- black27
     if ( Options%isRankVerbose() > 1 ) then
       if (theDatastore%root%has_path("metrics/global/sweep/number_of_unknowns")) then
          globalNumUnknowns = theDatastore%root%fetch_path_as_int64("metrics/global/sweep/number_of_unknowns")
          ! These can be very large values.  Size up everything to doubles to avoid overflows.
          throughputCycle = dble(globalNumUnknowns) * currentDtRad / (Size%RadtrTimeCycle * 60.0 + adqtEpsilon)
          Size%throughputTotal =  Size%throughputTotal + throughputCycle
          print '(1X,A52,1X,ES14.6)', "Cycle throughput (# elements in PSI * dt / radtr ) =", throughputCycle
          print '(1X,A27,1X,ES14.6)', "Average cycle throughput = ", Size%throughputTotal / dble(max(1,ncycle))
       else
          print *, "Cycle throughput: Unavailable, requires teton to be initialized with the C++ Teton::initialize() call."
       endif
     endif

     write(zoneStr, "(i7)") ConvControlZone
     write(procStr, "(i7)") ConvControlProcess

     print *,"*****************   Convergence    *****************"
     print *,"    Controlled by = ", ConvControlReasonC
     print *,"    ProcessID     = ", trim(procStr) 
     print *,"    Zone          = ", trim(zoneStr) 
     print Sformat,"     Rel Error     = ", ConvControlError
     print Sformat,"     Tr            = ", ConvState(1)
     print Sformat,"     Te            = ", ConvState(2)
     print Sformat,"     Rho           = ", ConvState(3)
     print Sformat,"     Cv            = ", ConvState(4)
     print Sformat,"     Source Rate   = ", ConvState(5)
     if (ndim == 1) then
       print format1D,"     Coordinates   =", ConvState(6)
     elseif (ndim == 2) then
       print format2D,"     Coordinates   =", ConvState(6),ConvState(7)
     elseif (ndim == 3) then
       print format3D,"     Coordinates   =", ConvState(6),ConvState(7),ConvState(8)
     endif
     print *," "

     write(cycleStr, "(i7)") ncycle + 1

     print *,"*****************  Time Step Vote  *****************"
     print *,"    For Cycle     = ", trim(cycleStr)
     print *,"    Controlled by = ", DtControlReasonC(indexDt)

     if (indexCaveat > 0) then
       print *,"    Caveat        = ", DtCaveat(indexCaveat)
     endif

     write(zoneStr, "(i7)") DtControlZone
     write(procStr, "(i7)") DtControlProcess

     print *,"    ProcessID     = ", trim(procStr) 
     print *,"    Control Zone  = ", trim(zoneStr) 
     print Sformat,"     Recommend Dt  = ", dtrad
     print Sformat,"     Max Change    = ", DtState(1)
     print Sformat,"     Tr            = ", DtState(2)
     print Sformat,"     Tr Old        = ", DtState(3)
     print Sformat,"     Te            = ", DtState(4)
     print Sformat,"     Te Old        = ", DtState(5)
     print Sformat,"     Rho           = ", DtState(6)
     print Sformat,"     Cv            = ", DtState(7)
     print Sformat,"     Source Rate   = ", DtState(8)

     if (ndim == 1) then
       print format1D,"     Coordinates   =", DtState(9)
     elseif (ndim == 2) then
       print format2D,"     Coordinates   =", DtState(9),DtState(10)
     elseif (ndim == 3) then
       print format3D,"     Coordinates   =", DtState(9),DtState(10),DtState(11)
     endif

     print *," "
     flush(stdout)

   endif

!  Construct a message for the host code

   messageStr = messageRoot(indexDt)

   write(changeStr,"(1pe13.6)") DtState(1) 
   write(zoneStr, "(i8)")       DtControlZone
   write(procStr, "(i6)")       DtControlProcess

   if (ndim == 1) then
     write(coordStr,"(1pe12.4)") DtState(9)
   elseif (ndim == 2) then
     write(coordStr,"(1pe12.4,1pe12.4)") DtState(9),DtState(10)
   elseif (ndim == 3) then
     write(coordStr,"(1pe12.4,1pe12.4,1pe12.4)") DtState(9),DtState(10),DtState(11)
   endif

   if (indexDt <= 3) then
     messageStr = trim(messageStr)//changeStr//" in zone"//trim(zoneStr)  &
                  //" on process"//trim(procStr)//" @ ("//trim(coordStr)//" )"
   else
     messageStr = trim(messageStr)//" in zone"//trim(zoneStr)  &
                  //" on process"//trim(procStr)//" @ ("//trim(coordStr)//" )"
   endif

!   print *, messageStr

   DtControls% dtMessage = messageStr


   return
   end subroutine getRunStats


