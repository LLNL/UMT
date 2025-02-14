! Editor Module:  Contains various problem edits 
                                                                                 
module Editor_mod 

  use kind_mod
  use constant_mod

  private

! public interfaces

  public construct 
  public destruct
  public setEdits
  public getTrMaxZone 
  public getTeMaxZone
  public getTrMaxProcess 
  public getTeMaxProcess
  public getTrMax
  public getTeMax
  public getNumberOfSpectrumAngleBins
  public getDeltaEMat
  public getEnergyRadiation
  public getPowerIncident 
  public getPowerEscape
  public getPowerAbsorbed 
  public getPowerEmitted
  public getPowerExtSources 
  public getPowerCompton
  public getEnergyCheck

  type, public :: Editor 

     integer                 :: numSpectrumAngleBins

     integer, dimension (1)  :: TrMaxZone           ! Zone with max value of Tr
     integer, dimension (1)  :: TeMaxZone           ! Zone with max value of Te

     integer                 :: TrMaxProcess        ! Process with max Tr
     integer                 :: TeMaxProcess        ! Process with max Te

     real(adqt)              :: TrMax               ! Max value of Tr
     real(adqt)              :: TeMax               ! Max value of Te

!    Instantaneous (per cycle) global edits
     real(adqt)              :: deltaEMat           ! Energy deposited in material (this is what gets returned to the host code)
     real(adqt)              :: EnergyRadiation     ! Energy in radiation field
     real(adqt)              :: PowerIncident       ! Power incident on external boundaries
     real(adqt)              :: PowerEscape         ! Power escaping from external boundaries
     real(adqt)              :: PowerAbsorbed       ! Power absorbed
     real(adqt)              :: PowerEmitted        ! Power emitted
     real(adqt)              :: PowerExtSources     ! Power from external sources
     real(adqt)              :: PowerCompton        ! Net Power absorbed by the material from Compton exchange 
     real(adqt)              :: EnergyRadBOC        ! Beginning of cycle radiation energy
     real(adqt)              :: EnergyCheck         ! Energy check or error

!    Group and Angle/Group dependent Edits on external boundaries 
     real(adqt), contiguous, pointer :: spectrumAngleBinBoundaries(:) => null() ! polar sector bin boundaries 
     real(adqt), contiguous, pointer :: RadPowerEscape(:)             => null() ! Power escaping by group 
     real(adqt), contiguous, pointer :: RadPowerIncident(:)           => null() ! Power incident by group 
     real(adqt), contiguous, pointer :: PolarSectorPowerEscape(:)     => null() ! Power escaping by polar sector and group

  end type Editor 

  type(Editor), pointer, public :: RadEdit => null()

  interface construct
    module procedure Editor_ctor
  end interface

  interface destruct
    module procedure Editor_dtor
  end interface

  interface setEdits
    module procedure Editor_set
  end interface

  interface getTrMaxZone
    module procedure Editor_get_TrMaxZone
  end interface

  interface getTeMaxZone
    module procedure Editor_get_TeMaxZone
  end interface

  interface getTrMaxProcess
    module procedure Editor_get_TrMaxProcess
  end interface
                                                                                       
  interface getTeMaxProcess
    module procedure Editor_get_TeMaxProcess
  end interface
                                                                                       
  interface getTrMax
    module procedure Editor_get_TrMax
  end interface

  interface getTeMax
    module procedure Editor_get_TeMax
  end interface

  interface getNumberOfSpectrumAngleBins
    module procedure Editor_get_NumSpectrumAngleBins
  end interface

  interface getDeltaEMat
    module procedure Editor_get_DeltaEMat
  end interface

  interface getEnergyRadiation
    module procedure Editor_get_EnergyRadiation
  end interface

  interface getPowerIncident
    module procedure Editor_get_PowerIncident
  end interface

  interface getPowerEscape
    module procedure Editor_get_PowerEscape
  end interface

  interface getPowerAbsorbed
    module procedure Editor_get_PowerAbsorbed
  end interface

  interface getPowerEmitted
    module procedure Editor_get_PowerEmitted
  end interface

  interface getPowerExtSources
    module procedure Editor_get_PowerExtSources
  end interface

  interface getPowerCompton
    module procedure Editor_get_PowerCompton
  end interface

  interface getEnergyCheck
    module procedure Editor_get_EnergyCheck
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================
                                                                                   
  subroutine Editor_ctor(self,                       &
                         ngr,                        &
                         numSpectrumAngleBins,       &
                         spectrumAngleBinBoundaries, &
                         RadPowerEscape,             &
                         RadPowerIncident,           &
                         PolarSectorPowerEscape)

    implicit none

!   Passed variables

    type(Editor),    intent(inout) :: self
    integer,         intent(in)    :: ngr
    integer,         intent(in)    :: numSpectrumAngleBins
    real(adqt), target, intent(in) :: spectrumAngleBinBoundaries(numSpectrumAngleBins+1)
    real(adqt), target, intent(in) :: RadPowerEscape(ngr)
    real(adqt), target, intent(in) :: RadPowerIncident(ngr)
    real(adqt), target, intent(in) :: PolarSectorPowerEscape(numSpectrumAngleBins*ngr)

!   Local

    integer    :: bin

!   Initialize scalers

    self% numSpectrumAngleBins =  numSpectrumAngleBins 
    self% TrMaxZone            =  0
    self% TeMaxZone            =  0
    self% TrMaxProcess         = -1
    self% TeMaxProcess         = -1

    self% TrMax                = zero
    self% TeMax                = zero

    self% deltaEMat            = zero
    self% EnergyRadiation      = zero
    self% PowerIncident        = zero
    self% PowerEscape          = zero
    self% PowerAbsorbed        = zero
    self% PowerEmitted         = zero
    self% PowerExtSources      = zero
    self% PowerCompton         = zero
    self% EnergyRadBOC         = zero
    self% EnergyCheck          = zero

    allocate( self% spectrumAngleBinBoundaries(numSpectrumAngleBins+1) )

    self% RadPowerEscape             => RadPowerEscape
    self% RadPowerIncident           => RadPowerIncident
    self% PolarSectorPowerEscape     => PolarSectorPowerEscape

    do bin=1,numSpectrumAngleBins+1
      self% spectrumAngleBinBoundaries(bin) = spectrumAngleBinBoundaries(bin)
    enddo


    return

  end subroutine Editor_ctor

!=======================================================================
! destruct interface
!=======================================================================
                                                                                          
  subroutine Editor_dtor(self)
                                                                                          
    implicit none
                                                                                          
!   Passed variables
                                                                                          
    type(Editor),    intent(inout) :: self
                                                                                          
!   Zero scalers
                                                                                          
    self% numSpectrumAngleBins =  0 
    self% TrMaxZone            =  0
    self% TeMaxZone            =  0
    self% TrMaxProcess         = -1
    self% TeMaxProcess         = -1

    self% TrMax                = zero
    self% TeMax                = zero

    self% deltaEMat            = zero
    self% EnergyRadiation      = zero
    self% PowerIncident        = zero
    self% PowerEscape          = zero
    self% PowerAbsorbed        = zero
    self% PowerEmitted         = zero
    self% PowerExtSources      = zero
    self% PowerCompton         = zero
    self% EnergyRadBOC         = zero
    self% EnergyCheck          = zero

    deallocate( self% spectrumAngleBinBoundaries )

    return
          
  end subroutine Editor_dtor

!=======================================================================
! set interface
!=======================================================================
                                                                                         
  subroutine Editor_set(self,             &
                        TrMaxZone,        &
                        TeMaxZone,        &
                        TrMaxProcess,     &
                        TeMaxProcess,     &
                        TrMax,            &
                        TeMax,            &
                        deltaEMat,        &
                        EnergyRadiation,  &
                        PowerIncident,    &
                        PowerEscape,      &
                        PowerAbsorbed,    &
                        PowerEmitted,     &
                        PowerExtSources,  &
                        PowerCompton,     &
                        EnergyCheck)

    implicit none
                                                                                         
!   Passed variables
                                                                                         
    type(Editor),         intent(inout)             :: self
    integer,    optional, intent(in), dimension (1) :: TrMaxZone
    integer,    optional, intent(in), dimension (1) :: TeMaxZone
    integer,    optional, intent(in)                :: TrMaxProcess
    integer,    optional, intent(in)                :: TeMaxProcess
    real(adqt), optional, intent(in)                :: TrMax 
    real(adqt), optional, intent(in)                :: TeMax 
    real(adqt), optional, intent(in)                :: deltaEMat
    real(adqt), optional, intent(in)                :: EnergyRadiation
    real(adqt), optional, intent(in)                :: PowerIncident
    real(adqt), optional, intent(in)                :: PowerEscape
    real(adqt), optional, intent(in)                :: PowerAbsorbed
    real(adqt), optional, intent(in)                :: PowerEmitted
    real(adqt), optional, intent(in)                :: PowerExtSources
    real(adqt), optional, intent(in)                :: PowerCompton
    real(adqt), optional, intent(in)                :: EnergyCheck

!   Update the Timestep control object

    if (present(TrMaxZone)) then
      self% TrMaxZone = TrMaxZone 
    endif

    if (present(TeMaxZone)) then
      self% TeMaxZone = TeMaxZone
    endif

    if (present(TrMaxProcess)) then
      self% TrMaxProcess = TrMaxProcess
    endif

    if (present(TeMaxProcess)) then
      self% TeMaxProcess = TeMaxProcess
    endif

    if (present(TrMax)) then
      self% TrMax = TrMax
    endif

    if (present(TeMax)) then
      self% TeMax = TeMax
    endif

    if (present(deltaEMat)) then
      self% deltaEMat = deltaEMat
    endif

    if (present(EnergyRadiation)) then
      self% EnergyRadiation = EnergyRadiation
    endif

    if (present(PowerIncident)) then
      self% PowerIncident = PowerIncident 
    endif

    if (present(PowerEscape)) then
      self% PowerEscape = PowerEscape 
    endif

    if (present(PowerAbsorbed)) then
      self% PowerAbsorbed = PowerAbsorbed    
    endif

    if (present(PowerEmitted)) then
      self% PowerEmitted = PowerEmitted
    endif

    if (present(PowerExtSources)) then
      self% PowerExtSources = PowerExtSources 
    endif

    if (present(PowerCompton)) then
      self% PowerCompton = PowerCompton
    endif

    if (present(EnergyCheck)) then
      self% EnergyCheck = EnergyCheck 
    endif


    return
 
  end subroutine Editor_set

!=======================================================================
! getTrMaxZone interface
!=======================================================================
  function Editor_get_TrMaxZone(self) result(TrMaxZone)

!    Return the zone with maximum radiation temperature (TrMaxZone)  

!    variable declarations
     implicit none

!    passed variables
     type(Editor), intent(in) :: self
     integer, dimension (1)   :: TrMaxZone 

     TrMaxZone = self% TrMaxZone 

     return
  end function Editor_get_TrMaxZone

!=======================================================================
! getTeMaxZone interface
!=======================================================================
  function Editor_get_TeMaxZone(self) result(TeMaxZone)

!    Return the zone with maximum electron temperature (TeMaxZone)

!    variable declarations
     implicit none

!    passed variables
     type(Editor), intent(in) :: self
     integer, dimension (1)   :: TeMaxZone

     TeMaxZone = self% TeMaxZone

     return
  end function Editor_get_TeMaxZone

!=======================================================================
! getTrMaxProcess interface
!=======================================================================
  function Editor_get_TrMaxProcess(self) result(TrMaxProcess)
                                                                                       
!    Return the node with maximum radiation temperature (TrMaxProcess)
                                                                                       
!    variable declarations
     implicit none
                                                                                       
!    passed variables
     type(Editor), intent(in) :: self
     integer                  :: TrMaxProcess
                                                                                       
     TrMaxProcess = self% TrMaxProcess
                                                                                       
     return
  end function Editor_get_TrMaxProcess

!=======================================================================
! getTeMaxProcess interface
!=======================================================================
  function Editor_get_TeMaxProcess(self) result(TeMaxProcess)
                                                                                       
!    Return the node with maximum electron temperature (TeMaxProcess)
                                                                                       
!    variable declarations
     implicit none
                                                                                       
!    passed variables
     type(Editor), intent(in) :: self
     integer                  :: TeMaxProcess
                                                                                       
     TeMaxProcess = self% TeMaxProcess
                                                                                       
     return
  end function Editor_get_TeMaxProcess

!=======================================================================
! getTrMax interface
!=======================================================================
  function Editor_get_TrMax(self) result(TrMax)

!    Return the maximum radiation temperature (TrMax)

!    variable declarations
     implicit none

!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: TrMax

     TrMax = self% TrMax

     return
  end function Editor_get_TrMax

!=======================================================================
! getTeMax interface
!=======================================================================
  function Editor_get_TeMax(self) result(TeMax)

!    Return the maximum electron temperature (TeMax)

!    variable declarations
     implicit none
 
!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: TeMax
 
     TeMax = self% TeMax
 
     return
  end function Editor_get_TeMax

!=======================================================================
! getNumberOfSpectrumAngleBins interface
!=======================================================================
  function Editor_get_NumSpectrumAngleBins(self) result(numSpectrumAngleBins)

!    Return the number of spectrum angle bins (numSpectrumAngleBins)

!    variable declarations
     implicit none

!    passed variables
     type(Editor), intent(in) :: self
     integer                  :: numSpectrumAngleBins 

     numSpectrumAngleBins = self% numSpectrumAngleBins 

     return
  end function Editor_get_NumSpectrumAngleBins

!=======================================================================
! getDeltaEMat interface
!=======================================================================

  function Editor_get_DeltaEMat(self) result(deltaEMat)

!    Return the material energy deposited (deltaEMat).
!       This is the quantity that the host code sees.

!    variable declarations
     implicit none

!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: deltaEMat

     deltaEMat = self% deltaEMat

     return
  end function Editor_get_DeltaEMat

!=======================================================================
! getEnergyRadiation interface
!=======================================================================
                                                                                         
  function Editor_get_EnergyRadiation(self) result(EnergyRadiation)
                                                                                         
                                                                                         
!    Return the radiation energy (EnergyRadiation)
                                                                                         
!    variable declarations
     implicit none
                                                                                         
!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: EnergyRadiation

     EnergyRadiation = self% EnergyRadiation

     return
  end function Editor_get_EnergyRadiation

!=======================================================================
! getPowerIncident interface
!=======================================================================
                                                                                          
  function Editor_get_PowerIncident(self) result(PowerIncident)
                                                                                          
                                                                                          
!    Return the incident power on the boundary (PowerIncident)
                                                                                          
!    variable declarations
     implicit none
                                                                                          
!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: PowerIncident
                                                                                          
     PowerIncident = self% PowerIncident
                                                                                          
     return
  end function Editor_get_PowerIncident

!=======================================================================
! getPowerEscape interface
!=======================================================================
                                                                                          
  function Editor_get_PowerEscape(self) result(PowerEscape)
                                                                                          
                                                                                          
!    Return the escaping power from the boundary (PowerEscape)
                                                                                          
!    variable declarations
     implicit none
                                                                                          
!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: PowerEscape
                                                                                          
     PowerEscape = self% PowerEscape
                                                                                          
     return
  end function Editor_get_PowerEscape

!=======================================================================
! getPowerAbsorbed interface
!=======================================================================

  function Editor_get_PowerAbsorbed(self) result(PowerAbsorbed)


!    Return the power absorbed by the material (PowerAbsorbed)

!    variable declarations
     implicit none

!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: PowerAbsorbed

     PowerAbsorbed = self% PowerAbsorbed

     return
  end function Editor_get_PowerAbsorbed

!=======================================================================
! getPowerEmitted interface
!=======================================================================

  function Editor_get_PowerEmitted(self) result(PowerEmitted)


!    Return the power emitted by the material (PowerEmitted)

!    variable declarations
     implicit none

!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: PowerEmitted

     PowerEmitted = self% PowerEmitted

     return
  end function Editor_get_PowerEmitted

!=======================================================================
! getPowerExtSources interface
!=======================================================================
                                                                                          
  function Editor_get_PowerExtSources(self) result(PowerExtSources)
                                                                                          
                                                                                          
!    Return the power from external sources (PowerExtSources)
                                                                                          
!    variable declarations
     implicit none
                                                                                          
!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: PowerExtSources
                                                                                          
     PowerExtSources = self% PowerExtSources
                                                                                          
     return
  end function Editor_get_PowerExtSources

!=======================================================================
! getPowerCompton interface
!=======================================================================

  function Editor_get_PowerCompton(self) result(PowerCompton)


!    Return the power transfered from radiation to electrons through Compton 
!    scattering (PowewrCompton)

!    variable declarations
     implicit none

!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: PowerCompton

     PowerCompton = self% PowerCompton

     return
  end function Editor_get_PowerCompton

!=======================================================================
! getEnergyCheck interface
!=======================================================================
                                                                                          
  function Editor_get_EnergyCheck(self) result(EnergyCheck)
                                                                                          
                                                                                          
!    Return the energy check or error (EnergyCheck)
                                                                                          
!    variable declarations
     implicit none
                                                                                          
!    passed variables
     type(Editor), intent(in) :: self
     real(adqt)               :: EnergyCheck
                                                                                          
     EnergyCheck = self% EnergyCheck
                                                                                          
     return
  end function Editor_get_EnergyCheck


end module Editor_mod

