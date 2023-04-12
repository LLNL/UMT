! Boundary Module:  Contains data structures for boundary fluxes 

module Boundary_mod 

  use kind_mod

  private

! public interfaces

  public construct
  public destruct 
  public getNumberOfBdyElements
  public getFirstBdyElement 
  public getNeighborID
  public getBCType
  public constructReflectedAngle 
  public destructReflectedAngle
  public setReflectedAngle
  public getReflectedAngle

  type, public :: ReflectedAngle
     integer, allocatable  :: ReflAngle(:)        ! list of angle IDs
  end type ReflectedAngle

  type, public :: Boundary 

     integer                  :: NumBdyElem          ! number of boundary elements 
     integer                  :: BdyElemCtr          ! used for counting new boundary elements
     integer                  :: BdyElem1            ! index of first boundary element
     integer                  :: NeighborID          ! shared process ID
     integer                  :: BCType              ! boundary type (defined in flags_mod)

     integer,    allocatable  :: BdyToC(:)           ! BdyToC(NumBdyElem)
     integer,    allocatable  :: BdyToZone(:)        ! BdyToZone(NumBdyElem)

     real(adqt), allocatable  :: A_bdy(:,:)          ! A_bdy(ndim,NumBdyElem)
     real(adqt), allocatable  :: Radius(:)           ! Radius(NumBdyElem)

     real(adqt), allocatable  :: nodePosition(:,:)   ! Used to check alignment of
     real(adqt), allocatable  :: nodePosition2(:,:)  ! shared boundary elements

     type(ReflectedAngle), pointer :: iRef(:) => null() ! Pointers to reflected angle data

  end type Boundary 

  type(Boundary), pointer, public :: Bdy => null()

  interface construct
    module procedure Boundary_ctor
  end interface

  interface getNumberOfBdyElements
    module procedure Boundary_getNumBdyElem
  end interface

  interface getFirstBdyElement
    module procedure Boundary_getFirstBE
  end interface

  interface getNeighborID
    module procedure Boundary_getNeigh
  end interface

  interface getBCType
    module procedure Boundary_getBCType
  end interface

  interface constructReflectedAngle
    module procedure Boundary_ctorRef
  end interface

  interface destructReflectedAngle
    module procedure Boundary_dtorRef
  end interface
                                                                                                   
  interface setReflectedAngle
    module procedure Boundary_setRef
  end interface
                                                                                                   
  interface getReflectedAngle
    module procedure Boundary_getRef
  end interface

  interface destruct
    module procedure Boundary_dtor
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================
                                                                                   
  subroutine Boundary_ctor(self,        &
                           BCType,      &
                           NumBdyElem,  &
                           BdyElem1,    &
                           NeighborID)

    use flags_mod
    use Size_mod
    use QuadratureList_mod

    implicit none

!   Passed variables

    type(Boundary), intent(inout)  :: self

    integer, intent(in)           :: BCType
    integer, intent(in)           :: NumBdyElem         
    integer, intent(in)           :: BdyElem1
    integer, intent(in)           :: NeighborID

!   Set Properties

    self% NumBdyElem = NumBdyElem 
    self% BdyElemCtr = 0
    self% BdyElem1   = BdyElem1
    self% NeighborID = NeighborID
    self% BCType     = BCType

    allocate( self % BdyToC(self% NumBdyElem) )
    allocate( self % A_bdy(Size%ndim,self% NumBdyElem) )

    if (Size%ndim <= 2) then
      allocate( self % Radius(self% NumBdyElem) )
    endif

    allocate( self% BdyToZone(self% NumBdyElem) )

    if (self% BCType == bcType_refl) then
      allocate( self % iRef(2) )
    endif


    return

  end subroutine Boundary_ctor

!=======================================================================
! construct Reflected Angle interface
!=======================================================================
                                                                                                   
  subroutine Boundary_ctorRef(self)
                                                                                                   
    use QuadratureList_mod
    use Quadrature_mod

    implicit none
                                                                                                   
!   Passed variables
    type(Boundary), intent(inout)  :: self

!   Local
    integer                        :: set
    type(ReflectedAngle), pointer  :: iRef
                                                                                                   
!   Allocate space

    do set=1,2
      QuadSet => getQuadrature(Quad,set)
      iRef    => self% iRef(set)

      allocate( iRef% ReflAngle(QuadSet% NumAngles) )

      iRef% ReflAngle(:) = -1
    enddo

    return
                                                                                                   
  end subroutine Boundary_ctorRef

!=======================================================================
! destruct Reflected Angle interface
!=======================================================================
                                                                                                   
  subroutine Boundary_dtorRef(self)

    use QuadratureList_mod

    implicit none
                                                                                                   
!   Passed variables
    type(Boundary), intent(inout)  :: self

!   Local
    integer                        :: set
    type(ReflectedAngle), pointer  :: iRef

!   Release space

    do set=1,2
      iRef    => self% iRef(set)
                                                                                                   
      deallocate( iRef% ReflAngle )
    enddo
                                                                                                   
    return
                                                                                                   
  end subroutine Boundary_dtorRef

!=======================================================================
! get Reflected Angle interface
!=======================================================================
                                                                                                   
  function Boundary_getRef(self, QuadID, IncAngle) result(ReflAngle)
                                                                                                   
    use Quadrature_mod

    implicit none
                                                                                                   
!   Passed variables
                                                                                                   
    type(Boundary),       intent(inout)  :: self
    integer,              intent(in)     :: QuadID
    integer,              intent(in)     :: IncAngle
    integer                              :: ReflAngle

!   Local
    type(ReflectedAngle), pointer        :: iRef

                                                                                                   
    iRef => self% iRef(QuadID)

    ReflAngle = iRef% ReflAngle(IncAngle)
                                                                                                   
    return
                                                                                                   
  end function Boundary_getRef

!=======================================================================
! set Reflected Angle interface
!=======================================================================
                                                                                                   
  subroutine Boundary_setRef(self, QuadID, IncAngle, ReflAngle)
                                                                                                   
    use Quadrature_mod

    implicit none
                                                                                                   
!   Passed variables
                                                                                                   
    type(Boundary),       intent(inout)  :: self
    integer,              intent(in)     :: QuadID
    integer,              intent(in)     :: IncAngle
    integer,              intent(in)     :: ReflAngle

!   Local
    type(ReflectedAngle), pointer        :: iRef
                                                                                                   
    iRef => self% iRef(QuadID)                                                                                                 
    iRef% ReflAngle(IncAngle) = ReflAngle
                                                                                                   
    return
                                                                                                   
  end subroutine Boundary_setRef

!=======================================================================
! getNumberOfBdyElements interface
!=======================================================================
  function Boundary_getNumBdyElem(self) result(NumBdyElem)

!    Returns the number of boundary elements 

!    variable declarations
     implicit none

!    passed variables
     type(Boundary), intent(in) :: self
     integer                    :: NumBdyElem

     NumBdyElem = self% NumBdyElem

     return

  end function Boundary_getNumBdyElem 

!=======================================================================
! getFirstBdyElement interface
!=======================================================================
  function Boundary_getFirstBE(self) result(BdyElem1)
                                                                                                    
!    Returns the first boundary element
                                                                                                    
!    variable declarations
     implicit none
                                                                                                    
!    passed variables
     type(Boundary), intent(in) :: self
     integer                    :: BdyElem1 
                                                                                                    
     BdyElem1 = self% BdyElem1 
                                                                                                    
     return
                                                                                                    
  end function Boundary_getFirstBE

!=======================================================================
! getNeighborID interface
!=======================================================================
  function Boundary_getNeigh(self) result(NeighborID)
                                                                                                    
!    Returns the NeighborID
                                                                                                    
!    variable declarations
     implicit none
                                                                                                    
!    passed variables
     type(Boundary), intent(in) :: self
     integer                    :: NeighborID
                                                                                                    
     NeighborID = self% NeighborID
                                                                                                    
     return
                                                                                                    
  end function Boundary_getNeigh

!=======================================================================
! getBCType interface
!=======================================================================
  function Boundary_getBCType(self) result(BCType)

!    Returns the type 
                                                                                                 
!    variable declarations
     implicit none
                                                                                                 
!    passed variables
     type(Boundary), intent(in) :: self
     integer                    :: BCType 
                                                                                                 
     BCType = self% BCType 
                                                                                                 
     return
                                                                                                 
  end function Boundary_getBCType

!=======================================================================
! destruct interface
!=======================================================================
                                                                                    
  subroutine Boundary_dtor(self)

    use Size_mod
    use flags_mod

    implicit none

!   Passed variables
                                                                                     
    type(Boundary),  intent(inout) :: self

!   Local

    deallocate( self% BdyToC )
    deallocate( self% BdyToZone )
    deallocate( self% A_bdy  )

    if (Size% ndim <= 2) then
      deallocate( self% Radius )
    endif

    if (self% BCType == bcType_refl) then
      call destructReflectedAngle(self)

      deallocate( self% iRef )
    endif

!  NOTE: nodePosition and nodePosition2 were deallocated immediatel;y after
!  use in aux/checkSharedBoundary

    return

  end subroutine Boundary_dtor

end module Boundary_mod

