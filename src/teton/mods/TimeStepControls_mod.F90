! Time Step Control Module:  Contains controls and edits for the time step 
                                                                                 
module TimeStepControls_mod 

  use kind_mod
  use constant_mod
  use flags_mod

  private

! public interfaces

  public construct 
  public destruct 
  public setDtControls
  public getRadCycle
  public getRadTime 
  public getRadTimeStep 
  public getRecTimeStep
  public getMaxChangeTe 
  public getMaxChangeEr
  public getMinTimeStep 
  public getMaxTimeStep
  public getMaxFracChangeEr
  public getMaxFracChangeTe
  public getZoneMaxChangeEr
  public getZoneMaxChangeTe
  public getControlProcess 
  public getControlZone
  public getDtConstraint 

  type, public :: TimeStepControls 

     integer            :: RadCycle
     integer            :: ZoneMaxChangeEr   ! Zone with max change in Tr**4
     integer            :: ZoneMaxChangeTe   ! Zone with max change in Te
     integer            :: ControlProcess    ! Process controlling time step
     integer            :: ControlZone       ! Zone controlling time step
                                                                                 
     real(adqt)         :: RadTimeStep       ! Current radiation time step
     real(adqt)         :: RecTimeStep       ! Recommended time step for next cycle
     real(adqt)         :: MaxFracChangeEr   ! Max fractional change in Tr**4
     real(adqt)         :: MaxFracChangeTe   ! Max fractional change in Te

     real(adqt)         :: RadTime           ! Current radiation time
     real(adqt)         :: MaxChangeTe       ! Max allowed change in Te per cycle 
     real(adqt)         :: MaxChangeEr       ! Max allowed change in Tr**4 per cycle
     real(adqt)         :: MinTimeStep       ! Minimum allowed time step
     real(adqt)         :: MaxTimeStep       ! Maximum allowed time step

     integer            :: DtConstraint      ! flag for what is controlling the time step 
     character(len=180) :: dtMessage

  end type TimeStepControls 

  type(TimeStepControls), pointer, public :: DtControls => null()

  interface construct
    module procedure TimeStepControls_ctor
  end interface

  interface destruct
    module procedure TimeStepControls_dtor
  end interface

  interface setDtControls
    module procedure TimeStepControls_set
  end interface

  interface getRadCycle
    module procedure TimeStepControls_get_RadCycle
  end interface

  interface getRadTime
    module procedure TimeStepControls_get_RadTime
  end interface

  interface getRadTimeStep
    module procedure TimeStepControls_get_RadTimeStep
  end interface

  interface getRecTimeStep
    module procedure TimeStepControls_get_RecTimeStep
  end interface

  interface getMaxChangeTe
    module procedure TimeStepControls_get_MaxChangeTe
  end interface

  interface getMaxChangeEr
    module procedure TimeStepControls_get_MaxChangeEr
  end interface

  interface getMinTimeStep
    module procedure TimeStepControls_get_MinTimeStep
  end interface

  interface getMaxTimeStep
    module procedure TimeStepControls_get_MaxTimeStep
  end interface

  interface getMaxFracChangeEr
    module procedure TimeStepControls_get_MaxFracChangeEr
  end interface

  interface getMaxFracChangeTe
    module procedure TimeStepControls_get_MaxFracChangeTe
  end interface

  interface getZoneMaxChangeEr
    module procedure TimeStepControls_get_ZoneMaxChangeEr
  end interface

  interface getZoneMaxChangeTe
    module procedure TimeStepControls_get_ZoneMaxChangeTe
  end interface

  interface getControlProcess
    module procedure TimeStepControls_get_ControlProcess
  end interface

  interface getControlZone
    module procedure TimeStepControls_get_ControlZone
  end interface

  interface getDtConstraint
    module procedure TimeStepControls_get_DtConstraint
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================
                                                                                   
  subroutine TimeStepControls_ctor(self,RadTimeStep,MaxChangeTe,  &
                                   MaxChangeEr,      &
                                   MinTimeStep,MaxTimeStep)

    implicit none

!   Passed variables

    type(TimeStepControls),    intent(inout) :: self
    real(adqt), optional, intent(in)         :: RadTimeStep 
    real(adqt), optional, intent(in)         :: MaxChangeTe
    real(adqt), optional, intent(in)         :: MaxChangeEr
    real(adqt), optional, intent(in)         :: MinTimeStep
    real(adqt), optional, intent(in)         :: MaxTimeStep 

!   Construct the Timestep control object

    if (present(RadTimeStep)) then
      self % RadTimeStep = RadTimeStep
      self % RecTimeStep = RadTimeStep
    else
      self % RadTimeStep = 1.0e-3_adqt
      self % RecTimeStep = 1.0e-3_adqt 
    endif

    if (present(MaxChangeTe)) then
      self % MaxChangeTe = MaxChangeTe 
    else
      self % MaxChangeTe = 4.0e-1_adqt
    endif

    if (present(MaxChangeEr)) then
      self % MaxChangeEr = MaxChangeEr
    else
      self % MaxChangeEr = 4.0e-1_adqt
    endif

    if (present(MinTimeStep)) then
      self % MinTimeStep = MinTimeStep
    else
      self % MinTimeStep = 1.0e-4_adqt
    endif

    if (present(MaxTimeStep)) then
      self % MaxTimeStep = MaxTimeStep
    else
      self % MaxTimeStep = 1.0e-1_adqt
    endif

    self % ZoneMaxChangeEr  = -1
    self % ZoneMaxChangeTe  = -1
    self % ControlProcess   = -1

    self % MaxFracChangeEr  = zero
    self % MaxFracChangeTe  = zero
    self % RadTime          = zero

    self % DtConstraint     = dtControl_none


    return

  end subroutine TimeStepControls_ctor
                                                      
!=======================================================================
! destruct interface
!=======================================================================
                                                                                    
  subroutine TimeStepControls_dtor(self)

    implicit none

!   Passed variables

    type(TimeStepControls),    intent(inout) :: self
                                                                
    self % ZoneMaxChangeEr  = -1
    self % ZoneMaxChangeTe  = -1
    self % ControlProcess   = -1
                                                                                         
    self % RadTimeStep      = zero
    self % RecTimeStep      = zero
    self % MaxFracChangeEr  = zero
    self % MaxFracChangeTe  = zero
                                                                                         
    self % RadTime          = zero
    self % MaxChangeTe      = zero
    self % MaxChangeEr      = zero
    self % MinTimeStep      = zero
    self % MaxTimeStep      = zero
                                                                                         
    self % DtConstraint     = dtControl_invalid

    WRITE(6,*) "Teton has destroyed the time step control"


    return

  end subroutine TimeStepControls_dtor

!=======================================================================
! set interface
!=======================================================================
                                                                                         
  subroutine TimeStepControls_set(self,             &
                                  ControlProcess,   &
                                  ControlZone,      &
                                  ZoneMaxChangeEr,  &
                                  ZoneMaxChangeTe,  &
                                  RadCycle,         &
                                  RadTimeStep,      &
                                  RecTimeStep,      &
                                  MaxFracChangeEr,  &
                                  MaxFracChangeTe,  &
                                  RadTime,          &
                                  DtConstraint      )
                                                                                         
    implicit none
                                                                                         
!   Passed variables
                                                                                         
    type(TimeStepControls),    intent(inout) :: self
    integer,    optional, intent(in)         :: ControlProcess
    integer,    optional, intent(in)         :: ControlZone
    integer,    optional, intent(in)         :: ZoneMaxChangeEr
    integer,    optional, intent(in)         :: ZoneMaxChangeTe
    integer,    optional, intent(in)         :: RadCycle

    real(adqt), optional, intent(in)         :: RadTimeStep
    real(adqt), optional, intent(in)         :: RecTimeStep 
    real(adqt), optional, intent(in)         :: MaxFracChangeEr
    real(adqt), optional, intent(in)         :: MaxFracChangeTe
    real(adqt), optional, intent(in)         :: RadTime

    integer, optional, intent(in)            :: DtConstraint
                                                                                         
!   Update the Timestep control object

    if (present(ControlProcess)) then
      self % ControlProcess = ControlProcess 
    endif

    if (present(ControlZone)) then
      self % ControlZone = ControlZone
    endif

    if (present(ZoneMaxChangeEr)) then
      self % ZoneMaxChangeEr = ZoneMaxChangeEr
    endif

    if (present(ZoneMaxChangeTe)) then
      self % ZoneMaxChangeTe = ZoneMaxChangeTe
    endif

    if (present(RadCycle)) then
      self% RadCycle = RadCycle
    endif

    if (present(RadTimeStep)) then
      self % RadTimeStep = RadTimeStep
    endif

    if (present(RecTimeStep)) then
      self % RecTimeStep = RecTimeStep
    endif

    if (present(MaxFracChangeEr)) then
      self % MaxFracChangeEr = MaxFracChangeEr
    endif

    if (present(MaxFracChangeTe)) then
      self % MaxFracChangeTe = MaxFracChangeTe
    endif

    if (present(RadTime)) then
      self % RadTime = RadTime 
    endif

    if (present(DtConstraint)) then
      self % DtConstraint = DtConstraint 
    else
      self % DtConstraint = dtControl_none
    endif

    return
 
  end subroutine TimeStepControls_set

!=======================================================================
! getRadTime interface
!=======================================================================
  function TimeStepControls_get_RadTime(self) result(RadTime)

!    Return the current radiation time (RadTime)  

!    variable declarations
     implicit none

!    passed variables
     type(TimeStepControls), intent(in) :: self
     real(adqt)                         :: RadTime 

     RadTime = self % RadTime 

     return
  end function TimeStepControls_get_RadTime

!=======================================================================
! getRadCycle interface
!=======================================================================
  function TimeStepControls_get_RadCycle(self) result(RadCycle)

!    Return the current radiation time (RadTime)  

!    variable declarations
     implicit none

!    passed variables
     type(TimeStepControls), intent(in) :: self
     integer                            :: RadCycle

     RadCycle = self% RadCycle

     return
  end function TimeStepControls_get_RadCycle

!=======================================================================
! getRadTimeStep interface
!=======================================================================
  function TimeStepControls_get_RadTimeStep(self) result(RadTimeStep)

!    Return the current radiation time step (RadTimeStep)

!    variable declarations
     implicit none

!    passed variables
     type(TimeStepControls), intent(in) :: self
     real(adqt)                         :: RadTimeStep

     RadTimeStep = self % RadTimeStep

     return
  end function TimeStepControls_get_RadTimeStep

!=======================================================================
! getRecTimeStep interface
!=======================================================================
  function TimeStepControls_get_RecTimeStep(self) result(RecTimeStep)
                                                                                       
!    Return the recommended radiation time step (RecTimeStep)
                                                                                       
!    variable declarations
     implicit none
                                                                                       
!    passed variables
     type(TimeStepControls), intent(in) :: self
     real(adqt)                         :: RecTimeStep
       
     RecTimeStep = self % RecTimeStep
       
     return
  end function TimeStepControls_get_RecTimeStep

!=======================================================================
! getMaxChangeTe interface
!=======================================================================
  function TimeStepControls_get_MaxChangeTe(self) result(MaxChangeTe)

!    Return the Maximum allowed change in Te per cycle (MaxChangeTe)

!    variable declarations
     implicit none

!    passed variables
     type(TimeStepControls), intent(in) :: self
     real(adqt)                         :: MaxChangeTe 

     MaxChangeTe = self % MaxChangeTe 

     return
  end function TimeStepControls_get_MaxChangeTe

!=======================================================================
! getMaxChangeEr interface
!=======================================================================
  function TimeStepControls_get_MaxChangeEr(self) result(MaxChangeEr)

!    Return the Maximum allowed change in Er per cycle (MaxChangeEr)

!    variable declarations
     implicit none

!    passed variables
     type(TimeStepControls), intent(in) :: self
     real(adqt)                         :: MaxChangeEr

     MaxChangeEr = self % MaxChangeEr

     return
  end function TimeStepControls_get_MaxChangeEr

!=======================================================================
! getMinTimeStep interface
!=======================================================================
  function TimeStepControls_get_MinTimeStep(self) result(MinTimeStep)

!    Return the minimum allowed time step (MinTimeStep)

!    variable declarations
     implicit none

!    passed variables
     type(TimeStepControls), intent(in) :: self
     real(adqt)                         :: MinTimeStep 

     MinTimeStep = self % MinTimeStep 

     return
  end function TimeStepControls_get_MinTimeStep

!=======================================================================
! getMaxTimeStep interface
!=======================================================================
  function TimeStepControls_get_MaxTimeStep(self) result(MaxTimeStep)

!    Return the maximum allowed time step (MaxTimeStep)

!    variable declarations
     implicit none

!    passed variables
     type(TimeStepControls), intent(in) :: self
     real(adqt)                         :: MaxTimeStep

     MaxTimeStep = self % MaxTimeStep

     return
  end function TimeStepControls_get_MaxTimeStep

!=======================================================================
! getMaxFracChangeEr interface
!=======================================================================
  function TimeStepControls_get_MaxFracChangeEr(self) result(MaxFracChangeEr)
                                                                                       
!    Return the maximum observed fractional change in Er (MaxFracChangeEr)
                                                                                       
!    variable declarations
     implicit none
                                                                                       
!    passed variables
     type(TimeStepControls), intent(in) :: self
     real(adqt)                         :: MaxFracChangeEr
                                                                                       
     MaxFracChangeEr = self % MaxFracChangeEr
                                                                                       
     return
  end function TimeStepControls_get_MaxFracChangeEr

!=======================================================================
! getMaxFracChangeTe interface
!=======================================================================
  function TimeStepControls_get_MaxFracChangeTe(self) result(MaxFracChangeTe)
                                                                                       
!    Return the maximum observed fractional change in Te (MaxFracChangeTe)
                                                                                       
!    variable declarations
     implicit none
                                                                                       
!    passed variables
     type(TimeStepControls), intent(in) :: self
     real(adqt)                         :: MaxFracChangeTe
                                                                                       
     MaxFracChangeTe = self % MaxFracChangeTe
                                                                                       
     return
  end function TimeStepControls_get_MaxFracChangeTe

!=======================================================================
! getZoneMaxChangeEr interface
!=======================================================================
  function TimeStepControls_get_ZoneMaxChangeEr(self) result(ZoneMaxChangeEr)
                                                                                       
!    Return the zone with the maximum observed fractional 
!    change in Er (ZoneMaxChangeEr)
                                                                                       
!    variable declarations
     implicit none
                                                                                       
!    passed variables
     type(TimeStepControls), intent(in) :: self
     integer                            :: ZoneMaxChangeEr
                                                                                       
     ZoneMaxChangeEr = self % ZoneMaxChangeEr
                                                                                       
     return
  end function TimeStepControls_get_ZoneMaxChangeEr

!=======================================================================
! getZoneMaxChangeTe interface
!=======================================================================
  function TimeStepControls_get_ZoneMaxChangeTe(self) result(ZoneMaxChangeTe)
                                                                                       
!    Return the zone with the maximum observed fractional
!    change in Te (ZoneMaxChangeTe)
                                                                                       
!    variable declarations
     implicit none
                                                                                       
!    passed variables
     type(TimeStepControls), intent(in) :: self
     integer                            :: ZoneMaxChangeTe
                                                                                       
     ZoneMaxChangeTe = self % ZoneMaxChangeTe
                                                                                       
     return
  end function TimeStepControls_get_ZoneMaxChangeTe

!=======================================================================
! getControlProcess interface
!=======================================================================
  function TimeStepControls_get_ControlProcess(self) result(ControlProcess)
                                                                                       
!    Return the process controlling the time step (ControlProcess)

!    variable declarations
     implicit none
                                                                                       
!    passed variables
     type(TimeStepControls), intent(in) :: self
     integer                            :: ControlProcess 
                                                                                       
     ControlProcess = self % ControlProcess 
                                                                                       
     return
  end function TimeStepControls_get_ControlProcess

!=======================================================================
! getControlZone interface
!=======================================================================
  function TimeStepControls_get_ControlZone(self) result(ControlZone)

!    Return the zone controlling the time step (ControlZone)

!    variable declarations
     implicit none

!    passed variables
     type(TimeStepControls), intent(in) :: self
     integer                            :: ControlZone

     ControlZone = self% ControlZone

     return
  end function TimeStepControls_get_ControlZone

!=======================================================================
! getDtConstraint interface
!=======================================================================
  function TimeStepControls_get_DtConstraint(self) result(DtConstraint)
                                                                                       
!    Return the process controlling the time step (DtConstraint)
                                                                                       
!    variable declarations
     implicit none
                                                                                       
!    passed variables
     type(TimeStepControls), intent(in) :: self
     integer                            :: DtConstraint 
                                                                                       
     DtConstraint = self % DtConstraint 
                                                                                       
     return
  end function TimeStepControls_get_DtConstraint


end module TimeStepControls_mod

