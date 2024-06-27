!***********************************************************************
!                                                                      *
!    PUBLISHEDITS -  Called from host to get current values for        *
!                various edits so that they can be printed             *
!                to the log file.  These values were computed in       *
!                BoundaryEdit at the end of the TRT solve (i.e., at    *
!                the end of the teton_radtr() call)                    *
!                                                                      *
!***********************************************************************

   subroutine publishEdits(dtrad)  BIND(C,NAME="teton_publishedits")

   use, intrinsic :: iso_c_binding, only: c_int, c_bool, c_ptr, c_null_ptr, c_double
   use, intrinsic :: iso_fortran_env, only : stdin=>input_unit,   &
                                             stdout=>output_unit, &
                                             stderr=>error_unit

   use kind_mod
   use Size_mod
   use iter_control_list_mod
   use iter_control_mod
   use TimeStepControls_mod
   use Editor_mod
   use Datastore_mod
   use Options_mod, only : Options

   implicit none


!  Arguments

   real(C_DOUBLE), intent(inout) :: dtrad

!  local arguments

   integer(C_INT) :: noutrt
   integer(C_INT) :: ninrt
   integer(C_INT) :: ngdart
   integer(C_INT) :: nNLIters
   integer(C_INT) :: maxNLIters
   integer(C_INT), dimension (1) :: TrMaxZone
   integer(C_INT), dimension (1) :: TeMaxZone
   integer(C_INT) :: TrMaxProcess
   integer(C_INT) :: TeMaxProcess

   real(C_DOUBLE) :: dtused
   real(C_DOUBLE) :: TrMax
   real(C_DOUBLE) :: TeMax
   real(C_DOUBLE) :: deltaEMat
   real(C_DOUBLE) :: EnergyRadiation
   real(C_DOUBLE) :: PowerIncident
   real(C_DOUBLE) :: PowerEscape
   real(C_DOUBLE) :: PowerAbsorbed
   real(C_DOUBLE) :: PowerEmitted
   real(C_DOUBLE) :: PowerExtSources
   real(C_DOUBLE) :: PowerCompton
   real(C_DOUBLE) :: EnergyCheck

   type(IterControl), pointer    :: temperatureControl => NULL()
   type(IterControl), pointer    :: intensityControl   => NULL()
   type(IterControl), pointer    :: incidentFluxControl=> NULL()
   type(IterControl), pointer    :: greyControl        => NULL()
   type(IterControl), pointer    :: nonLinearControl   => NULL()

!  Local

   integer(C_INT) :: ncycle
   integer(C_INT) :: outerMaxIts, greyMaxIts, incidentFluxMaxIts, innerNLMaxIts
   real(C_DOUBLE) :: timerad
   real(C_DOUBLE) :: outerTempRelTol, outerEDRelTol, greyRelTol, incidentFluxRelTol, innerNLRelTol

   character(len=38), parameter :: Cformat = "(1X,A10,i8,A12,F18.10,A10,1pe18.10)"
   character(len=25), parameter :: Iformat = "(1X,A13,i6,A15,i6,A15,i6)"
   character(len=18), parameter :: Jformat = "(1X,A20,i6,A21,i6)"
   character(len=30), parameter :: Tformat = "(1X,A7,1X,F18.10,A9,i7,A12,i5)"
   character(len=55), parameter :: Eformat = "(1X,A31,1X,1pe18.10,1X,A13,1X,1pe18.10,1X,A15,1pe18.10)"
   character(len=17), parameter :: Dformat = "(1X,A43,1pe18.10)"

 ! Iteration Controls

   temperatureControl  => getIterationControl(IterControls,"temperature")
   intensityControl    => getIterationControl(IterControls,"intensity")
   incidentFluxControl => getIterationControl(IterControls, "incidentFlux")
   greyControl         => getIterationControl(IterControls,"grey")
   nonLinearControl    => getIterationControl(IterControls,"nonLinear")

   noutrt     = getNumberOfIterations(temperatureControl)
   ninrt      = getNumberOfIterations(intensityControl)
   ngdart     = getNumberOfIterations(greyControl)
   nNLIters   = getNumberOfIterations(nonLinearControl)
   maxNLIters = getGlobalMaxIterationsTaken(nonLinearControl)

!  Tolerance control settings
   outerMaxIts = getMaxNumberOfIterations(temperatureControl)
   greyMaxIts = getMaxNumberOfIterations(greyControl)
   incidentFluxMaxIts = getMaxNumberOfIterations(incidentFluxControl)
   innerNLMaxIts = getMaxNumberOfIterations(nonLinearControl)
   outerTempRelTol = getEpsilonPoint(temperatureControl)
   outerEDRelTol = getEpsilonPoint(intensityControl)
   greyRelTol = getEpsilonPoint(greyControl)
   incidentFluxRelTol = getEpsilonPoint(incidentFluxControl)
   innerNLRelTol = getEpsilonPoint(nonLinearControl)


!  Time Step Controls

   ncycle     = getRadCycle(DtControls)
   timerad    = getRadTime(DtControls)
   dtrad      = getRecTimeStep(DtControls)
   dtused     = getRadTimeStep(DtControls)

!  Energy and Temperature Edits

   TrMaxZone       = getTrMaxZone(RadEdit)
   TeMaxZone       = getTeMaxZone(RadEdit)
   TrMaxProcess    = getTrMaxProcess(RadEdit)
   TeMaxProcess    = getTeMaxProcess(RadEdit)
   TrMax           = getTrMax(RadEdit)
   TeMax           = getTeMax(RadEdit)
   deltaEMat       = getDeltaEMat(RadEdit)
   EnergyRadiation = getEnergyRadiation(RadEdit)
   PowerIncident   = getPowerIncident(RadEdit)
   PowerEscape     = getPowerEscape(RadEdit)
   PowerAbsorbed   = getPowerAbsorbed(RadEdit)
   PowerEmitted    = getPowerEmitted(RadEdit)
   PowerExtSources = getPowerExtSources(RadEdit)
   PowerCompton    = getPowerCompton(RadEdit)
   EnergyCheck     = getEnergyCheck(RadEdit)

   ! store this in the datastore
   call theDatastore%root%set_path("rtedits/noutrt", noutrt)
   call theDatastore%root%set_path("rtedits/ninrt", ninrt)
   call theDatastore%root%set_path("rtedits/ngdart", ngdart)
   call theDatastore%root%set_path("rtedits/nNLIters", nNLIters)
   call theDatastore%root%set_path("rtedits/maxNLIters", maxNLIters)
   call theDatastore%root%set_path("rtedits/TrMaxZone", TrMaxZone(1))
   call theDatastore%root%set_path("rtedits/TeMaxZone", TeMaxZone(1))
   call theDatastore%root%set_path("rtedits/TrMaxProcess", TrMaxProcess)
   call theDatastore%root%set_path("rtedits/TeMaxProcess", TeMaxProcess)
   call theDatastore%root%set_path("rtedits/TeMax", TeMax)
   call theDatastore%root%set_path("rtedits/TrMax", TrMax)
   call theDatastore%root%set_path("rtedits/deltaEMat", deltaEMat)
   call theDatastore%root%set_path("rtedits/EnergyRadiation", EnergyRadiation)
   call theDatastore%root%set_path("rtedits/PowerIncident", PowerIncident)
   call theDatastore%root%set_path("rtedits/PowerEscape", PowerEscape)
   call theDatastore%root%set_path("rtedits/PowerAbsorbed", PowerAbsorbed)
   call theDatastore%root%set_path("rtedits/PowerEmitted", PowerEmitted)
   call theDatastore%root%set_path("rtedits/PowerExtSources", PowerExtSources)
   call theDatastore%root%set_path("rtedits/PowerCompton", PowerCompton)
   call theDatastore%root%set_path("rtedits/EnergyCheck", EnergyCheck)

   if ( Options%isRankVerbose() > 0 ) then
      print *, ""
      print *, ">>>>>>>>>>>>>>>     End of Radiation Step Report    <<<<<<<<<<<<<<<"
      print Cformat, "TIME STEP ", ncycle,"  timerad = ", timerad,"  dtrad = ", dtused
      print *, ""
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
      print Iformat,"TempIters = ", noutrt, "  FluxIters = ", ninrt, "  GTASweeps = ",ngdart
      print Jformat,"AveNonLinearIters = ", nNLIters, "  MaxNonLinearIters = ",maxNLIters
#else
      print *, "FluxIters = ", ninrt
#endif
      if( Options%isRankVerbose() > 1 ) then
         print *, "  *** max outer iterations = ", outerMaxIts
         print *, "  *** max outer temperature rel tol = ", outerTempRelTol
         print *, "  *** max outer energy density rel tol = ", outerEDRelTol
         print *, "  *** max incident flux iterations = ", incidentFluxMaxIts
         print *, "  *** max incident flux rel tol = ", incidentFluxRelTol
         print *, "  *** max grey iterations = ", greyMaxIts
         print *, "  *** max grey rel tol = ", greyRelTol
         print *, "  *** max inner nonlinear iterations = ", innerNLMaxIts
         print *, "  *** max inner nonlinear rel tol = ", innerNLRelTol
      endif
      print Tformat, "TrMax =", TrMax, " in Zone ", &
                                TrMaxZone," on Process ",TrMaxProcess
      print Tformat, "TeMax =", TeMax, " in Zone ", &
                                TeMaxZone," on Process ",TeMaxProcess
      print Eformat, "Energy deposited in material = ", deltaEMat, "ERad total = ", EnergyRadiation, "Energy check = ",EnergyCheck
      print Dformat, "Recommended time step for next rad cycle = ", dtrad
      print *, ""
      flush(stdout)
   endif

   return
   end subroutine publishEdits

