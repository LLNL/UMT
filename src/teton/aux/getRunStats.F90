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
   
   use ISO_C_BINDING
   use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                             stdout=>output_unit, &
                                             stderr=>error_unit
   use flags_mod
   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use constant_mod
   use Material_mod
   use Geometry_mod
   use QuadratureList_mod
   use iter_control_list_mod
   use iter_control_mod
   use TimeStepControls_mod
   use Options_mod

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
   real(adqt)     :: DtControlChange
   real(adqt)     :: ConvState(5+Size% ndim)
   real(adqt)     :: DtState(8+Size% ndim)
   real(adqt)     :: zoneCenter(Size% ndim)

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

   ! Variables affecting threading
   integer        :: numOmpCPUThreads
   integer        :: nZoneSets
   integer        :: nSets
   integer        :: nHyperDomains

   character(len=26), parameter :: Tformat = "(1X,A16,1X,F14.8,5X,F14.8)" 
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
   numOmpCPUThreads = Options%getNumOmpMaxThreads()
   nZoneSets        = getNumberOfZoneSets(Quad)
   nSets            = getNumberOfSets(Quad)
   nHyperDomains    = getNumberOfHyperDomains(Quad,1) 

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

   indexCaveat = 0 

   if (dtrad == getMinTimeStep(DtControls)) then
     indexCaveat = 1 
   elseif (dtrad == getMaxTimeStep(DtControls)) then
     indexCaveat = 2 
   endif

   if (myRankInGroup == DtControlProcess) then

     if (DtConstraint == dtControl_radTemp) then
       DtControlChange  = getMaxFracChangeTr4(DtControls)
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
   RadtrTimeTotal    = Size% RadtrTimeTotal
   InitTimeTotal     = Size% InitTimeTotal
   FinalTimeTotal    = Size% FinalTimeTotal

   timeNonRadCycle   = timeNonRad/sixty
   timeNonRadTotal   = timeOther/sixty

   if ( Options%isRankVerbose() > 0 ) then

     ncycle = getRadCycle(DtControls) 

#if defined(TETON_ENABLE_OPENMP)
     print *,"*****************     Threading     ****************"
     print '(A,i5)', " # threads per rank, cpu        = ", numOmpCPUThreads
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
     if (Size%useGPU) then
! Number of thread teams used for kernels iterating over zone sets.
       print '(A,i5)', " # thread teams over zone sets  = ", nZoneSets
! Number of thread teams used for kernels iterating over phase-angle
! sets.  May comment this line out later, as Paul intends to migrate all
! kernels to be over zone sets.
       print '(A,i5)', " # thread teams over sweep sets = ", nSets*nHyperDomains
     endif
     print *," "
#endif
#endif
     print *,"*****************     Run Time     *****************"
     print *,"                    Cycle (min)     Accumulated (min)"
     print Tformat, "RADTR          =", Size% RadtrTimeCycle,    RadtrTimeTotal 
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     print Tformat, "Mat. Coupling  =", Size% MatCoupTimeCycle,  MatCoupTimeTotal 
#endif
     print Tformat, "Sweep(CPU)     =", Size% SweepTimeCycle,    SweepTimeTotal
     print Tformat, "Sweep(GPU)     =", Size% GPUSweepTimeCycle, GPUSweepTimeTotal
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     print Tformat, "Grey Tr. Accel =", Size% GTATimeCycle,      GTATimeTotal
#endif
     print Tformat, "Initialization =", Size% InitTimeCycle,     InitTimeTotal
     print Tformat, "Finalization   =", Size% FinalTimeCycle,    FinalTimeTotal
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
     print Tformat, "Non-Rad        =",       timeNonRadCycle,   timeNonRadTotal
#endif
     print *," "

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


