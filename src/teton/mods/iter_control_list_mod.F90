#include "macros.h"
!=======================================================================
!                       Version 1.0: 03/99, MRZ
!-----------------------------------------------------------------------
! Iteration Control List
!   This class contains the list of iteration controls in the problem
!
! nIterCon     number of iteration controls
! maxIterCon   maximum number of iteration controls
! names        iteration control names
! iControls    iteration controls
!-----------------------------------------------------------------------
! v1.0: Original implementation
!=======================================================================
module iter_control_list_mod
      use kind_mod
      use iter_control_mod
      use default_iter_controls_mod

private

! public interfaces
  public construct 
  public resetNumberOfIterations 
  public destruct
  public getNumberOfIterationControls 
  public getIterationControl

  type, public :: IterControlList
    private
    integer                        :: nIterCon
    integer                        :: maxIterCon
    character(12),     allocatable :: names(:)
    type(IterControl), pointer     :: iControls(:) => null()
  end type IterControlList

  type(IterControlList), pointer, public :: IterControls => null()

  interface construct
    module procedure iter_control_list_ctor
  end interface

  interface resetNumberOfIterations
    module procedure iter_control_list_reset_nIter
  end interface

  interface destruct
    module procedure iter_control_list_dtor
  end interface

  interface getNumberOfIterationControls
    module procedure iter_control_list_get_nIterCon
  end interface

  interface getIterationControl
    module procedure iter_control_list_get_iCon
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine iter_control_list_ctor(self)

!    Construct the iteration control list object

!    variable declarations
     implicit none

!    passed variables
     type(IterControlList), intent(inout) :: self

!    local variables
     integer :: nIterCon

!    initialize the maximum number of iteration controls
     self % maxIterCon = 10

!    allocate memory for the iteration control list
     allocate (self % names(self%maxIterCon))
     allocate (self % iControls(self%maxIterCon))

!    construct each of the default iteration control objects
     nIterCon = 0

     nIterCon = nIterCon + 1
     self%names(nIterCon) = "temperature"
     call construct(self%iControls(nIterCon), &
       epsilonPoint=outer_temp_reltol, &
       maxNumberOfIterations=outer_max_it)

     nIterCon = nIterCon + 1
     self%names(nIterCon) = "intensity"
     call construct(self%iControls(nIterCon), &
       epsilonPoint=outer_intensity_reltol, &
       maxNumberOfIterations=1)

     nIterCon = nIterCon + 1
     self%names(nIterCon) = "grey"
     call construct(self%iControls(nIterCon), &
       epsilonPoint=grey_reltol, &
       maxNumberOfIterations=grey_max_sweeps)

     nIterCon = nIterCon + 1
     self%names(nIterCon) = "incidentFlux"
     call construct(self%iControls(nIterCon), &
       epsilonPoint=incident_flux_reltol, &
       maxNumberOfIterations=incident_flux_max_it)

     nIterCon = nIterCon + 1
     self%names(nIterCon) = "nonLinear"
     call construct(self%iControls(nIterCon), &
       epsilonPoint=inner_nl_reltol, &
       maxNumberOfIterations=inner_nl_max_it)

     self % nIterCon = nIterCon

!    assertions
     TETON_ASSERT(self%nIterCon > 0,"Invalid iter control list")
     TETON_ASSERT(self%nIterCon<=self%maxIterCon,"Invalid iter control list")
     TETON_ASSERT(associated(self%iControls),"Invalid iter control list")

     return
  end subroutine iter_control_list_ctor

!=======================================================================
! resetNumberOfIterations interface
!=======================================================================

  subroutine iter_control_list_reset_nIter(self)

!    Reset the number of iterations in each iteration control object
!    in the list

!    variable declarations
     implicit none

!    passed variables
     type(IterControlList), intent(inout) :: self

!    local variables
     integer :: iIterCon

!    reset the number of iterations
     ListLoop: do iIterCon = 1, self % nIterCon
        call resetNumberOfIterations(self % iControls(iIterCon))
     enddo ListLoop

     return
  end subroutine iter_control_list_reset_nIter

!=======================================================================
! destruct interface
!=======================================================================

  subroutine iter_control_list_dtor(self)

!    Destruct the iteration control list object

!    variable declarations
     implicit none

!    passed variables
     type(IterControlList), intent(inout) :: self

!    local variables
     integer :: iIterCon, alloc_stat

!    destruct the list from the bottom up
     DestructLoop: do iIterCon = 1, self%nIterCon

        call destruct(self%iControls(iIterCon))
        self % names(iIterCon) = ""

     enddo DestructLoop

     deallocate(self % names, stat=alloc_stat)
     deallocate(self % iControls, stat=alloc_stat)

     self % nIterCon = 0

!    assertions
     TETON_ASSERT(.not.allocated(self%names),"Invalid destruct")
     TETON_ASSERT(.not.associated(self%iControls),"Invalid destruct")
     TETON_ASSERT(self%nIterCon==0,"Invalid destruct")

     return
  end subroutine iter_control_list_dtor

!=======================================================================
! external data access routines
!=======================================================================

  function iter_control_list_get_nIterCon(self) result(nIterCon)

!    Return the number of iteration controls in the list
!      nIterCon   number of iteration controls in the list

!    variable declarations
     implicit none

!    passed variables
     type(IterControlList), intent(in) :: self
     integer                           :: nIterCon

     nIterCon = self % nIterCon

!    assertions
     TETON_ASSERT(nIterCon==self%nIterCon,"Improper data access")
     TETON_ASSERT(nIterCon>=0,"Improper data access")

     return
  end function iter_control_list_get_nIterCon

!-----------------------------------------------------------------------
  function iter_control_list_get_iCon(self,iteration) result(iControl)

!    Return a pointer to an iteration control
!      iteration   iteration control name
!      iControl    pointer to the iteration control

!    variable declarations
     implicit none

!    passed variables
     type(IterControlList), intent(in) :: self
     character(*),          intent(in) :: iteration
     type(IterControl),     pointer    :: iControl

!    local variables
     integer :: iIterCon

!    assertions
     TETON_ASSERT(allocated(self%names),"Invalid iter control list")
     TETON_ASSERT(associated(self%iControls),"Invalid iter control list")
     TETON_ASSERT(any(iteration==self%names(:)),"Invalid iteration name")

     iControl => self% iControls(self% maxIterCon)

!    loop through the list to the requested iteration control
     ListLoop: do iIterCon = 1, self%nIterCon
        if (self%names(iIterCon) == iteration) then
           iControl => self % iControls(iIterCon)
           exit ListLoop
        endif
     enddo ListLoop

!    assertions
     TETON_ASSERT(associated(iControl),"Invalid iter control")

     return
  end function iter_control_list_get_iCon

end module iter_control_list_mod
