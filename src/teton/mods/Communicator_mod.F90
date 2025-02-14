! Communication Module:  Contains data structures that tell processors
!                        how to communicate with each other 
                                                                                 
module Communicator_mod 

  use kind_mod
  use constant_mod

  private

! public interfaces

  public constructBuffer
  public destructBuffer
  public constructFluxBuffer
  public destructFluxBuffer

  type, public :: CommunicateFlux

     integer                         :: nAngleBins
     integer                         :: irequest(2)

     real(adqt), pointer, contiguous :: IncFlux(:)
     real(adqt), pointer, contiguous :: ExitFlux(:)

  end type CommunicateFlux
                                                                                 
  type, public :: Communicator 

     integer                         :: nSend                    ! number of boundary elements to send
     integer                         :: nRecv                    ! number of boundary elements to receive

     integer,    pointer, contiguous :: ListRecv(:)   => null()  ! ListRecv(nRecv)
     integer,    pointer, contiguous :: ListSend(:,:) => null()  ! ListSend(2,nSend)
     integer                         :: irequest(2)              ! irequest(2)

     real(adqt), pointer, contiguous :: psibsend(:,:) => null()  ! psibsend(Groups,nSend) - send buffer
     real(adqt), pointer, contiguous :: psibrecv(:,:) => null()  ! psibrecv(Groups,nRecv) - receive buffer

  end type Communicator 

  type(Communicator), pointer, public :: Comm => null()

  interface constructBuffer
    module procedure Communicator_ctor_buffer
  end interface

  interface destructBuffer
    module procedure Communicator_dtor
  end interface

  interface constructFluxBuffer
    module procedure Communicator_ctor_Fluxbuffer
  end interface

  interface destructFluxBuffer
    module procedure Communicator_dtor_Fluxbuffer
  end interface

contains

!=======================================================================
! constructBuffer interface
!=======================================================================
                                                                                    
  subroutine Communicator_ctor_buffer(self, nSend, nRecv, Groups, &
                                      nGroupSets) 

    use MemoryAllocator_mod, only : Allocator
    implicit none

!   Passed variables
    type(Communicator),    intent(inout) :: self
    integer,               intent(in)    :: nSend
    integer,               intent(in)    :: nRecv
    integer,               intent(in)    :: Groups
    integer,               intent(in)    :: nGroupSets

    self% nSend = nSend
    self% nRecv = nRecv

    if (self% nSend > 0) then
      allocate( self% ListSend(2,nSend) )
      call Allocator%allocate(Allocator%use_for_comm_data,"Communicator","psibsend", self% psibsend, Groups, nSend*nGroupSets)
    endif

    if (self% nRecv > 0) then
      allocate( self% ListRecv(nRecv) )
      call Allocator%allocate(Allocator%use_for_comm_data,"Communicator","psibrecv", self% psibrecv, Groups, nRecv*nGroupSets)
    endif

    return

  end subroutine Communicator_ctor_buffer

!=======================================================================
! destruct interface
!=======================================================================
                                                                                    
  subroutine Communicator_dtor(self)

    use mpi_param_mod
    use mpif90_mod
    use MemoryAllocator_mod, only : Allocator

    implicit none

!   Passed variables
                                                                                     
    type(Communicator), intent(inout)    :: self

!   Local

    integer  :: ierr

!   Free Communication Requests

    if (self% nSend > 0) then
      call MPI_Request_Free(self%irequest(1), ierr)

      deallocate( self% ListSend  )
      call Allocator%deallocate(Allocator%use_for_comm_data,"Communicator","psibsend",  self% psibsend)

    endif
                                                                                    
    if (self% nRecv > 0) then
      call MPI_Request_Free(self%irequest(2), ierr)

      deallocate( self% ListRecv  )
      call Allocator%deallocate(Allocator%use_for_comm_data,"Communicator","psibrecv",  self% psibrecv)

    endif

    return

  end subroutine Communicator_dtor
!=======================================================================
! constructBuffer interface
!=======================================================================

  subroutine Communicator_ctor_Fluxbuffer(self, nAngleBins)

    implicit none

!   Passed variables

    type(CommunicateFlux), intent(inout) :: self
    integer,               intent(in)    :: nAngleBins


    self% nAngleBins = nAngleBins

    allocate( self% IncFlux(self% nAngleBins) )
    allocate( self% ExitFlux(self% nAngleBins) )

    return

  end subroutine Communicator_ctor_Fluxbuffer

!=======================================================================
! destruct interface
!=======================================================================

  subroutine Communicator_dtor_Fluxbuffer(self)

    use mpi_param_mod
    use mpif90_mod

    implicit none

!   Passed variables

    type(CommunicateFlux), intent(inout)    :: self

!   Local

    integer  :: ierr

!   Free Communication Requests

    call MPI_Request_Free(self%irequest(1), ierr)

    deallocate( self% IncFlux  )

    call MPI_Request_Free(self%irequest(2), ierr)

    deallocate( self% ExitFlux  )


    return

  end subroutine Communicator_dtor_Fluxbuffer
end module Communicator_mod

