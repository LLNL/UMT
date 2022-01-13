! BoundaryList Module:  Contains attributes for problem boundaries 

module BoundaryList_mod

  use kind_mod
  use Boundary_mod
  use Profile_mod

  private

! public interfaces

  public construct 
  public setBoundary 
  public getNumberOfReflecting
  public getNumberOfBoundaries 
  public getNumberOfVacuum
  public getNumberOfShared 
  public getNumberOfSource
  public getReflecting 
  public getVacuum 
  public getSource
  public getShared 
  public getBoundary 
  public getProfile 
  public setProfile
  public resetSourceProfile
  public destruct

  public innerBdyReflect 
  public outerBdyReflect
  public innerBdyShared  
  public outerBdyShared
  public getInnerBdyID
  public getOuterBdyID
  public getInnerSharedID 
  public getOuterSharedID

  type, public :: BoundaryList

     integer                 :: NumBoundary        ! number of total boundaries
     integer                 :: NumReflecting      ! number of reflecting boundaries 
     integer                 :: NumVacuum          ! number of vacuum boundaries
     integer                 :: NumSource          ! number of source boundaries
     integer                 :: NumShared          ! number of shared boundaries

     integer                 :: PtrVac             ! points to vacuum boundary
     integer                 :: PtrSrc             ! points to source boundary
     integer                 :: PtrShared          ! points to shared boundary

     integer                 :: ProfileCounter     ! counter for adding source profiles

!    Special for 1D

     integer                 :: innerBdyID
     integer                 :: outerBdyID
     integer                 :: innerSharedID
     integer                 :: outerSharedID

     logical(kind=1)         :: inner1DBdyReflect
     logical(kind=1)         :: outer1DBdyReflect
     logical(kind=1)         :: inner1DBdyShared
     logical(kind=1)         :: outer1DBdyShared

     type(Boundary), pointer :: iBoundary(:) => null() ! boundary flux pointers
     type(Profile),  pointer :: iProfile(:)  => null() ! profile pointers

  end type BoundaryList

  type(BoundaryList), pointer, public :: RadBoundary => null()

  interface construct
    module procedure BoundaryList_ctor
  end interface

  interface setBoundary 
    module procedure BoundaryList_set
  end interface

  interface getReflecting
    module procedure BoundaryList_getReflBdy
  end interface

  interface getVacuum
    module procedure BoundaryList_getVacBdy
  end interface

  interface getSource
    module procedure BoundaryList_getSrcBdy
  end interface

  interface getShared
    module procedure BoundaryList_getSharBdy
  end interface

  interface getBoundary
    module procedure BoundaryList_getBdy
  end interface

  interface getProfile
    module procedure BoundaryList_getProf
  end interface

  interface setProfile
    module procedure BoundaryList_setProf
  end interface

  interface resetSourceProfile
    module procedure BoundaryList_resetProf
  end interface

  interface getNumberOfBoundaries
    module procedure BoundaryList_getNum
  end interface

  interface getNumberOfReflecting
    module procedure BoundaryList_getNumRefl
  end interface

  interface getNumberOfVacuum
    module procedure BoundaryList_getNumVac
  end interface

  interface getNumberOfShared
    module procedure BoundaryList_getNumShared
  end interface

  interface getNumberOfSource
    module procedure BoundaryList_getNumSrc
  end interface

  interface destruct
    module procedure BoundaryList_dtor
  end interface

  interface innerBdyReflect
    module procedure BoundaryList_innerBdyReflect
  end interface

  interface outerBdyReflect
    module procedure BoundaryList_outerBdyReflect
  end interface

  interface innerBdyShared
    module procedure BoundaryList_innerBdyShared
  end interface

  interface outerBdyShared
    module procedure BoundaryList_outerBdyShared
  end interface

  interface getInnerBdyID
    module procedure BoundaryList_getInnerBdyID
  end interface

  interface getOuterBdyID
    module procedure BoundaryList_getOuterBdyID
  end interface

  interface getInnerSharedID
    module procedure BoundaryList_getInnerSharedID
  end interface

  interface getOuterSharedID
    module procedure BoundaryList_getOuterSharedID
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine BoundaryList_ctor(self, NumReflecting, &
                                     NumVacuum,     &
                                     NumSource,     &
                                     NumShared)

    implicit none

!   Passed variables

    type(BoundaryList), intent(inout) :: self
    integer, intent(in)               :: NumReflecting
    integer, intent(in)               :: NumVacuum
    integer, intent(in)               :: NumSource
    integer, intent(in)               :: NumShared 

!   Initialize type counters
                                                                                                    
    self% NumReflecting  = NumReflecting 
    self% NumVacuum      = NumVacuum
    self% NumSource      = NumSource
    self% NumShared      = NumShared

    self% NumBoundary    = NumReflecting + NumVacuum +  &
                           NumSource     + NumShared 

    self% PtrVac         = NumReflecting
    self% PtrSrc         = NumReflecting + NumVacuum
    self% PtrShared      = NumReflecting + NumVacuum + NumSource

    self% ProfileCounter = 0

    allocate( self% iBoundary(self% NumBoundary) )
    allocate( self% iProfile(self% NumSource) )

!   Special for 1D

    self% innerBdyID    = -1
    self% outerBdyID    = -1
    self% innerSharedID = -1
    self% outerSharedID = -1

    self% inner1DBdyReflect = .FALSE.
    self% outer1DBdyReflect = .FALSE.
    self% inner1DBdyShared  = .FALSE.
    self% outer1DBdyShared  = .FALSE.

    return

  end subroutine BoundaryList_ctor

!=======================================================================
! set interface
!
!  PGM: updated 7/27/2017- propagated integer flags through for BcType
!=======================================================================

  subroutine BoundaryList_set(self,        &
                              BoundaryID,  &
                              BCType,      &
                              NumBdyElem,  &
                              BdyElem1,    &
                              NeighborID)

    implicit none

!   Passed variables

    type(BoundaryList), intent(inout)      :: self
    integer, intent(in)                    :: BoundaryID
    integer, intent(in)                    :: BCType
    integer, intent(in)                    :: NumBdyElem 
    integer, intent(in)                    :: BdyElem1 
    integer, intent(in)                    :: NeighborID

    

    call construct(self % iBoundary(BoundaryID), &
                          BCType,                &
                          NumBdyElem,            &
                          BdyElem1,              &
                          NeighborID)

    return

  end subroutine BoundaryList_set

!=======================================================================
! set Profile interface
!=======================================================================

  subroutine BoundaryList_setProf(self,            &
                                  NumTimes,        &
                                  NumValues,       &
                                  Multiplier,      &
                                  BlackBody,       &
                                  Isotropic,       &
                                  Times,           &
                                  Values,          &
                                  ProfileID)

    implicit none

!   Passed variables

    type(BoundaryList), intent(inout)      :: self
    integer, intent(in)                    :: NumTimes
    integer, intent(in)                    :: NumValues

    real(adqt), optional, intent(in)       :: Multiplier

    logical(kind=1), intent(in)            :: BlackBody
    logical(kind=1), optional, intent(in)  :: Isotropic

    real(adqt), intent(in)                 :: Times(NumTimes)
    real(adqt), intent(in)                 :: Values(NumValues)
    integer, intent(inout)                 :: ProfileID


    self% ProfileCounter = self% ProfileCounter + 1

    call construct(self% iProfile(self% ProfileCounter), &
                         NumTimes,                       &
                         NumValues,                      &
                         Multiplier,                     &
                         BlackBody,                      &
                         Isotropic,                      &
                         Times,                          &
                         Values)

    ProfileID = self%ProfileCounter

    return

  end subroutine BoundaryList_setProf

!=======================================================================
! reset Profile interface
!=======================================================================

  subroutine BoundaryList_resetProf(self,          &
                                  ProfileID,       &
                                  NumTimes,        &
                                  NumValues,       &
                                  Multiplier,      &
                                  Times,           &
                                  Values)

    implicit none

!   Passed variables

    type(BoundaryList), intent(inout)      :: self

    integer, intent(in)                    :: ProfileID
    integer, intent(in)                    :: NumTimes
    integer, intent(in)                    :: NumValues

    real(adqt), optional, intent(in)       :: Multiplier

    real(adqt), intent(in)                 :: Times(NumTimes)
    real(adqt), intent(in)                 :: Values(NumValues)


    call reset(self% iProfile(ProfileID),            &
                     NumTimes,                       &
                     NumValues,                      &
                     Multiplier,                     &
                     Times,                          &
                     Values)

    return

  end subroutine BoundaryList_resetProf

!=======================================================================
! getBoundary interface
!=======================================================================
  function BoundaryList_getBdy(self,BdyID) result(iBoundary)

!    Return a pointer to a boundary definition
                                                                                                  
!    variable declarations
     implicit none
     
!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer,            intent(in)   :: BdyID
     type(Boundary),     pointer      :: iBoundary
     
     
     iBoundary => self % iBoundary(BdyID)
     
     return
     
  end function BoundaryList_getBdy

!=======================================================================
! getBoundary interface
!=======================================================================
  function BoundaryList_getProf(self,BdyID) result(iProfile)

!    Return a pointer to a profile definition

!    variable declarations
     implicit none


!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer,            intent(in)   :: BdyID
     type(Profile),      pointer      :: iProfile


     iProfile => self% iProfile(BdyID)

     return

  end function BoundaryList_getProf

!=======================================================================
! getReflectingBoundary interface
!=======================================================================
  function BoundaryList_getReflBdy(self,BdyID) result(iBoundary)
                                                                                                    
!    Return a pointer to a boundary definition 

!    variable declarations
     implicit none
                                                                                                    
!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer,            intent(in)   :: BdyID
     type(Boundary),     pointer      :: iBoundary
                                                                                                    
 
     iBoundary => self % iBoundary(BdyID)
 
     return
 
  end function BoundaryList_getReflBdy

!=======================================================================
! getVacuumBoundary interface
!=======================================================================
  function BoundaryList_getVacBdy(self,BdyID) result(iBoundary)

!    Return a pointer to a boundary definition

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer,            intent(in)   :: BdyID
     type(Boundary),     pointer      :: iBoundary


     iBoundary => self % iBoundary(self%PtrVac + BdyID)

     return

  end function BoundaryList_getVacBdy

!=======================================================================
! getSourceBoundary interface
!=======================================================================
  function BoundaryList_getSrcBdy(self,BdyID) result(iBoundary)

!    Return a pointer to a boundary definition

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer,            intent(in)   :: BdyID
     type(Boundary),     pointer      :: iBoundary


     iBoundary => self % iBoundary(self%PtrSrc + BdyID)

     return

  end function BoundaryList_getSrcBdy

!=======================================================================
! getSharedBoundary interface
!=======================================================================
  function BoundaryList_getSharBdy(self,BdyID) result(iBoundary)

!    Return a pointer to a boundary definition

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer,            intent(in)   :: BdyID
     type(Boundary),     pointer      :: iBoundary


     iBoundary => self % iBoundary(self%PtrShared + BdyID)

     return

  end function BoundaryList_getSharBdy

!=======================================================================
! getNumberOfReflecting interface
!=======================================================================
  function BoundaryList_getNum(self) result(NumBoundary)

!    Return the number of reflecting boundaries
                                                                                                 
!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer                          :: NumBoundary
                                                                                                 
     NumBoundary = self % NumBoundary
                                                                                                 
     return
                                                                                                 
  end function BoundaryList_getNum

!=======================================================================
! getNumberOfReflecting interface
!=======================================================================
  function BoundaryList_getNumRefl(self) result(NumReflecting)
                                                                                                    
!    Return the number of reflecting boundaries 
                                                                                                    
!    variable declarations
     implicit none
                                                                                                    
!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer                          :: NumReflecting
                                                                                                    
     NumReflecting = self % NumReflecting
                                                                                                    
     return
                                                                                                    
  end function BoundaryList_getNumRefl

!=======================================================================
! getNumberOfVacuum interface
!=======================================================================
  function BoundaryList_getNumVac(self) result(NumVacuum)
                                                                                                    
!    Return the number of vacuum boundaries
                                                                                                    
!    variable declarations
     implicit none
                                                                                                    
!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer                          :: NumVacuum
                                                                                                    
     NumVacuum = self % NumVacuum 
                                                                                                    
     return
                                                                                                    
  end function BoundaryList_getNumVac

!=======================================================================
! getNumberOfShared interface
!=======================================================================
  function BoundaryList_getNumShared(self) result(NumShared)
                                                                                                    
!    Return the number of shared boundaries
                                                                                                    
!    variable declarations
     implicit none
                                                                                                    
!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer                          :: NumShared
                                                                                                    
     NumShared = self % NumShared
                                                                                                    
     return
                                                                                                    
  end function BoundaryList_getNumShared

!=======================================================================
! getNumberOfSource interface
!=======================================================================
  function BoundaryList_getNumSrc(self) result(NumSource)
                                                                                                    
!    Return the number of source boundaries
                                                                                                    
!    variable declarations
     implicit none
                                                                                                    
!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer                          :: NumSource
                                                                                                    
     NumSource = self % NumSource
                                                                                                    
     return
                                                                                                    
  end function BoundaryList_getNumSrc

!=======================================================================
! innerBdyReflect interface
!=======================================================================
  function BoundaryList_innerBdyReflect(self) result(inner1DBdyReflect)

!    Return the inner boundary element number for a 1D mesh 

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     logical(kind=1)                  :: inner1DBdyReflect 

     inner1DBdyReflect = self% inner1DBdyReflect 

     return

  end function BoundaryList_innerBdyReflect

!=======================================================================
! outerBdyReflect interface
!=======================================================================
  function BoundaryList_outerBdyReflect(self) result(outer1DBdyReflect)

!    Return the inner boundary element number for a 1D mesh 

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     logical(kind=1)                  :: outer1DBdyReflect

     outer1DBdyReflect = self% outer1DBdyReflect

     return

  end function BoundaryList_outerBdyReflect

!=======================================================================
! innerBdyShared interface
!=======================================================================
  function BoundaryList_innerBdyShared(self) result(inner1DBdyShared)

!    Return the inner boundary element number for a 1D mesh 

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     logical(kind=1)                  :: inner1DBdyShared

     inner1DBdyShared = self% inner1DBdyShared

     return

  end function BoundaryList_innerBdyShared

!=======================================================================
! outerBdyShared interface
!=======================================================================
  function BoundaryList_outerBdyShared(self) result(outer1DBdyShared)

!    Return the inner boundary element number for a 1D mesh 

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     logical(kind=1)                  :: outer1DBdyShared

     outer1DBdyShared = self% outer1DBdyShared

     return

  end function BoundaryList_outerBdyShared

!=======================================================================
! getInnerBdyID interface
!=======================================================================
  function BoundaryList_getInnerBdyID(self) result(innerBdyID)

!    Return the inner boundary element number for a 1D mesh 

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer                          :: innerBdyID 

     innerBdyID = self% innerBdyID 

     return

  end function BoundaryList_getInnerBdyID

!=======================================================================
! getOuterBdyID interface
!=======================================================================
  function BoundaryList_getOuterBdyID(self) result(outerBdyID)

!    Return the outer boundary element number for a 1D mesh 

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer                          :: outerBdyID

     outerBdyID = self% outerBdyID

     return

  end function BoundaryList_getOuterBdyID

!=======================================================================
! getInnerSharedID interface
!=======================================================================
  function BoundaryList_getInnerSharedID(self) result(innerSharedID)

!    Return the inner boundary element number for a 1D mesh 

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer                          :: innerSharedID

     innerSharedID = self% innerSharedID

     return

  end function BoundaryList_getInnerSharedID

!=======================================================================
! getOuterSharedID interface
!=======================================================================
  function BoundaryList_getOuterSharedID(self) result(outerSharedID)

!    Return the outer boundary element number for a 1D mesh 

!    variable declarations
     implicit none

!    passed variables
     type(BoundaryList), intent(in)   :: self
     integer                          :: outerSharedID

     outerSharedID = self% outerSharedID

     return

  end function BoundaryList_getOuterSharedID

!=======================================================================
! destruct interface
!=======================================================================

  subroutine BoundaryList_dtor(self)

    implicit none

!   Passed variables

    type(BoundaryList),  intent(inout) :: self

!   Local

    integer :: n

    do n=1,self% NumBoundary
      call destruct( self% iBoundary(n) ) 
    enddo

    do n=1,self% NumSource
      call destruct( self% iProfile(n) )
    enddo

    deallocate( self% iBoundary )
    deallocate( self% iProfile  )

    return

  end subroutine BoundaryList_dtor


end module BoundaryList_mod

