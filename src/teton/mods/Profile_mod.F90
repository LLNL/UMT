! Profile Module:  Contains data structures that describe source profiles 

module Profile_mod 

  use kind_mod
  use Size_mod

  private

! public interfaces

  public construct
  public reset
  public destruct
                                                                                 
  type, public :: Profile 

     integer                  :: NumTimes          ! number of profile times 
     integer                  :: NumValues         ! number of profile values 

     real(adqt)               :: Multiplier        ! profile multiplier

     logical(kind=1)          :: profileOn         ! profile on/off flag
     logical(kind=1)          :: BlackBody         ! frequency-dependence flag
     logical(kind=1)          :: Isotropic         ! angle-dependence flag

     real(adqt), allocatable  :: Times(:)          ! profile times
     real(adqt), allocatable  :: Values(:)         ! profile values
     real(adqt), allocatable  :: InterpValues(:)   ! interpolated values

  end type Profile 


  interface construct
    module procedure Profile_ctor
  end interface

  interface reset
    module procedure Profile_reset
  end interface

  interface destruct
    module procedure Profile_dtor
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================
                                                                                   
  subroutine Profile_ctor(self,            &
                          NumTimes,        &
                          NumValues,       &
                          Multiplier,      &
                          BlackBody,       &
                          Isotropic,       &
                          Times,           &
                          Values)

    implicit none

!   Passed variables

    type(Profile), intent(inout)           :: self

    integer, intent(in)                    :: NumTimes         
    integer, intent(in)                    :: NumValues       

    real(adqt), optional, intent(in)       :: Multiplier     

    logical(kind=1), intent(in)            :: BlackBody
    logical(kind=1),  optional, intent(in) :: Isotropic

    real(adqt), intent(in)                 :: Times(NumTimes)         
    real(adqt), intent(in)                 :: Values(NumValues)       

!   Set Properties

    self% NumTimes        = NumTimes
    self% NumValues       = NumValues

    if (present(Multiplier)) then
      self% Multiplier = Multiplier 
    else
      self% Multiplier = 1.0_adqt
    endif

    self% BlackBody = BlackBody 

    if (present(Isotropic)) then
      self% Isotropic = Isotropic 
    else
      self% Isotropic = .TRUE. 
    endif

    self% profileOn  = .FALSE. 

    allocate( self% Times(self % NumTimes) )
    allocate( self% Values(self % NumValues) )
    allocate( self% InterpValues(Size % ngr) )

    self% Times(:)  = Times(:)
    self% Values(:) = Values(:) 


    return

  end subroutine Profile_ctor

!=======================================================================
! construct interface
!=======================================================================

  subroutine Profile_reset(self,            &
                           NumTimes,        &
                           NumValues,       &
                           Multiplier,      &
                           Times,           &
                           Values)

    implicit none

!   Passed variables

    type(Profile), intent(inout)           :: self

    integer, intent(in)                    :: NumTimes
    integer, intent(in)                    :: NumValues

    real(adqt), optional, intent(in)       :: Multiplier

    real(adqt), intent(in)                 :: Times(NumTimes)
    real(adqt), intent(in)                 :: Values(NumValues)

!   Set Properties

    if (NumTimes /= self%NumTimes) then
      self% NumTimes = NumTimes
      deallocate( self% Times                  )
      allocate  ( self% Times(self % NumTimes) )
    endif
    if (NumValues /= self%NumValues) then
      self% NumValues = NumValues
      deallocate( self% Values                   )
      allocate  ( self% Values(self % NumValues) )
    endif

    if (present(Multiplier)) then
      self% Multiplier = Multiplier
    else
      self% Multiplier = 1.0_adqt
    endif

    self% profileOn  = .FALSE.

    self% Times(:)  = Times(:)
    self% Values(:) = Values(:)

    return

  end subroutine Profile_reset
                                                      
!=======================================================================
! destruct interface
!=======================================================================
                                                                                    
  subroutine Profile_dtor(self)


    implicit none

!   Passed variables
                                                                                     
    type(Profile),  intent(inout) :: self

    deallocate( self% Times )
    deallocate( self% Values )
    deallocate( self% InterpValues )

!   Local


    return

  end subroutine Profile_dtor

end module Profile_mod

