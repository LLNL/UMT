#include "macros.h"
! ZoneData Module:  Contains data structures for zone information 

module ZoneData_mod 

  use kind_mod
  use constant_mod

  private

! public interfaces

  public constructZone
  public destructZone

  type, public :: ZoneData 

     integer                  :: nCorner             ! number of corners
     integer                  :: nSides              ! number of sides
     integer                  :: c0                  ! global ID of first corner
     integer                  :: side0               ! global ID of first side
     integer                  :: nFaces              ! number of zone faces
     real(adqt)               :: EnergyDensityOld    ! old energy density
     logical(kind=1)          :: BoundaryZone        ! is zone on a boundary

!    For 1D-only
     real(adqt)               :: zoneWidth           ! delta radius
     real(adqt)               :: Rave                ! average radius 
     real(adqt)               :: Rmin                ! inner radius
     real(adqt)               :: Rmax                ! outer radius

  end type ZoneData 

  type(ZoneData), pointer, public :: Z => null()

  interface constructZone
    module procedure ZoneData_ctor,  &
                     ZoneData1D_ctor
  end interface

  interface destructZone
    module procedure ZoneData_dtor
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================
                                                                                   
  subroutine ZoneData_ctor(self,         &
                           nCorner,      &
                           corner0,      &
                           nFaces,       &
                           nSides,       &
                           side0,        &
                           cFP)

    use Size_mod

    implicit none

!   Passed variables

    type(ZoneData), intent(inout)    :: self

    integer, intent(in)              :: nCorner
    integer, intent(in)              :: corner0
    integer, intent(in)              :: nFaces
    integer, intent(in)              :: nSides 
    integer, intent(in)              :: side0 
    integer, intent(in)              :: cFP(Size%maxcf,Size%maxCorner) 

!   Local

    integer            :: cID 
    integer            :: i

!   Set Properties

    self% nCorner      = nCorner 
    self% c0           = corner0
    self% nFaces       = nFaces
    self% nSides       = nSides
    self% side0        = side0
    self% BoundaryZone = .FALSE.


!   Tag boundary zones

    do cID=1,self% nCorner
      do i=1,Size%maxcf
        if (cFP(i,cID) > Size% ncornr) then
          self% BoundaryZone = .TRUE.
        endif
      enddo
    enddo

    return

  end subroutine ZoneData_ctor

!=======================================================================
! construct 1D interface
!=======================================================================

  subroutine ZoneData1D_ctor(self,     &
                             corner0,  &  
                             BoundaryZone)

    use Size_mod
    use constant_mod

    implicit none

!   Passed variables

    type(ZoneData), intent(inout)    :: self

    integer, intent(in)              :: corner0
    logical(kind=1), intent(in)      :: BoundaryZone

!   Set Properties

    self% nCorner      = 2
    self% c0           = corner0
    self% nFaces       = 2
    self% BoundaryZone = BoundaryZone 

    self% zoneWidth    = zero
    self% Rave         = zero
    self% Rmin         = zero
    self% Rmax         = zero


    return

  end subroutine ZoneData1D_ctor

!=======================================================================
! destruct interface
!=======================================================================
                                                                                    
  subroutine ZoneData_dtor(self)

    use Size_mod

    implicit none

!   Passed variables

    type(ZoneData), intent(inout)    :: self

    SUPPRESS_UNUSED_VAR_WARNING(self)

    return

  end subroutine ZoneData_dtor


end module ZoneData_mod

