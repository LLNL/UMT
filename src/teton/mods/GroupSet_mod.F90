! Group Set Module:  Contains data structures for a group set 

module GroupSet_mod

  use kind_mod
  use constant_mod
  use, intrinsic :: iso_c_binding, only : c_double
  implicit none

  private

! public interfaces

  public construct
  public destruct

  type, public :: GroupSet

     integer                         :: Groups       ! number of energy groups 
     integer                         :: g0           ! group offset

     real(c_double), pointer, contiguous :: Sigt(:,:) => null()   ! total opacity
     real(c_double), pointer, contiguous :: STotal(:,:) => null() ! fixed + scat source

!    Misc                                                                                                                                                                                                                                 
     character(len=13) :: label ! A string descriptor for this set.

  end type GroupSet 

  interface construct
    module procedure GroupSet_ctor
  end interface

  interface destruct
    module procedure GroupSet_dtor
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine GroupSet_ctor(self,        &
                           Groups,      &
                           g0,          &
                           nZones,      &
                           nCorner)

    use Size_mod
    use constant_mod
    use MemoryAllocator_mod

    implicit none

!   Passed variables

    type(GroupSet), intent(inout)    :: self

    integer,        intent(in)       :: Groups
    integer,        intent(in)       :: g0
    integer,        intent(in)       :: nZones
    integer,        intent(in)       :: nCorner

!   Set Properties

    write(self%label,'(I0.3)') g0
    self% label = "group_set_"//self%label    

    self% Groups = Groups
    self% g0     = g0

!   Allocate Memory
    call Allocator%allocate(Size%usePinnedMemory,self%label, "Sigt",   self% Sigt,   Groups,nZones)
    call Allocator%allocate(Size%usePinnedMemory,self%label, "STotal", self% STotal, Groups,nCorner)

    self% Sigt(:,:)   = zero
    self% STotal(:,:) = zero

    return

  end subroutine GroupSet_ctor

!=======================================================================
! destruct interface
!=======================================================================

  subroutine GroupSet_dtor(self)

    use MemoryAllocator_mod
    use Size_mod

    implicit none

!   Passed variables

    type(GroupSet), intent(inout)    :: self

    call Allocator%deallocate(Size%usePinnedMemory,self%label, "Sigt",   self% Sigt)
    call Allocator%deallocate(Size%usePinnedMemory,self%label, "STotal", self% STotal)

    return

  end subroutine GroupSet_dtor

end module GroupSet_mod
