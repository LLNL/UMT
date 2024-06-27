! Communication Set Module:  Contains data structures for communication 

module CommSet_mod

  use kind_mod
  use constant_mod
  use Communicator_mod
  use AngleSet_mod
  USE ISO_C_BINDING

  private

! public interfaces

  public construct
  public destruct
  public getCommunicationGroup
  public getMessage
  public destructComm
  public restoreCommOrder
  public setCommOrder
  public getSharedFlux


  type, public :: CommSet

     integer                              :: COMM_GROUP
     integer                              :: set1
     integer                              :: set2
     integer                              :: Groups
     integer                              :: NumAngles 
     integer                              :: NumAnglesDyn
     integer                              :: NumBin
     integer                              :: NumBin0
     integer                              :: fluxSweeps

!    For monitoring convergence of sweep iterations
     real(adqt)                           :: totalIncFlux

     integer,         pointer, contiguous :: NangBinList(:)
     integer,         pointer, contiguous :: AngleToBin(:)
     integer,         pointer, contiguous :: AngleOrder(:)
     integer,         pointer, contiguous :: AngleOrder0(:)
     integer,         pointer, contiguous :: RecvOrder0(:,:)
     integer,         pointer, contiguous :: RecvOrder(:,:)
     integer,         pointer, contiguous :: request(:)

     real(adqt),      pointer, contiguous :: IncFlux(:)
     real(adqt),      pointer, contiguous :: IncFluxOld(:)
     real(adqt),      pointer, contiguous :: NetFlux(:,:)
     real(adqt),      pointer, contiguous :: relError(:)

     logical(kind=1), pointer, contiguous :: Converged(:)

!    Communication Pointers

     type(AngleSet),           pointer    :: AngleSetPtr    => null()
     type(Communicator),       pointer    :: CommPtr(:,:)   => null()
     type(CommunicateFlux),    pointer    :: CommFluxPtr(:) => null()

  end type CommSet 

  interface construct
    module procedure CommSet_ctor
  end interface

  interface destruct
    module procedure CommSet_dtor
  end interface

  interface getCommunicationGroup
    module procedure CommSet_getCommunicationGroup
  end interface

  interface getMessage
    module procedure CommSet_getMessage
  end interface

  interface destructComm
    module procedure CommSet_dtorComm
  end interface

  interface restoreCommOrder
    module procedure CommSet_resCommOrd
  end interface

  interface setCommOrder
    module procedure CommSet_setCommOrd
  end interface

  interface getSharedFlux
    module procedure CommSet_getSharedFlux
  end interface


contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine CommSet_ctor(self,        &
                          NumAngles,   &
                          set1,        &
                          set2,        &
                          Groups,      &
                          new_group,   &
                          AngleSetPtr)

    use Size_mod
    use constant_mod

    implicit none

!   Passed variables

    type(CommSet),  intent(inout)      :: self

    integer,        intent(in)         :: NumAngles
    integer,        intent(in)         :: set1
    integer,        intent(in)         :: set2
    integer,        intent(in)         :: Groups
    integer,        intent(in)         :: new_group !This is really a comm, not a comm group

    type(AngleSet), target, intent(in) :: AngleSetPtr

!   Local


!   Set Properties

    self% AngleSetPtr  => AngleSetPtr

    self% COMM_GROUP   =  new_group
    self% set1         =  set1
    self% set2         =  set2
    self% Groups       =  Groups
    self% NumAngles    =  NumAngles
    self% NumAnglesDyn =  NumAngles
    self% NumBin       =  AngleSetPtr% NumBin
    self% NumBin0      =  AngleSetPtr% NumBin0
    self% fluxSweeps   =  0
    self% totalIncFlux =  zero

!   Allocate Memory

    allocate( self% NangBinList(self% NumBin) )
    allocate( self% AngleToBin(self% NumAngles) )
    allocate( self% AngleOrder(self% NumAngles) )
    allocate( self% AngleOrder0(self% NumAngles) )

!   For testing boundary flux convergence

    allocate( self% IncFlux(self% NumBin) )
    allocate( self% IncFluxOld(self% NumBin) )
    allocate( self% relError(self% NumBin) )
    allocate( self% Converged(self% NumBin) )

!   Initialize

    self% NangBinList(:) = AngleSetPtr% NangBinList(:)
    self% AngleToBin(:)  = AngleSetPtr% AngleToBin(:)

    self% AngleOrder(:)  = 0
    self% AngleOrder0(:) = 0

    self% IncFlux(:)     = zero
    self% IncFluxOld(:)  = zero
    self% relError(:)    = zero
    self% Converged(:)   = .FALSE.

!   For MPI communiaction

    if (Size% ncomm > 0) then

      allocate( self% CommPtr(Size% ncomm,self% NumAngles) )
      allocate( self% CommFluxPtr(Size% ncomm) )

      allocate( self% RecvOrder(self% NumAngles,Size% ncomm) )
      allocate( self% RecvOrder0(self% NumAngles,Size% ncomm) )
      allocate( self% request(2*Size% ncomm) )

      allocate( self% NetFlux(Size% ncomm,self% NumBin) )

      self% RecvOrder(:,:)  = 0
      self% RecvOrder0(:,:) = 0
      self% request(:)      = 0

      self% NetFlux(:,:)    = zero
    endif



    return

  end subroutine CommSet_ctor

!=======================================================================
! destruct interface
!=======================================================================

  subroutine CommSet_dtor(self)

    use Size_mod

    implicit none

!   Passed variables

    type(CommSet), intent(inout)    :: self

   integer                   :: ierr

    if (Size% ncomm > 0) then

      deallocate( self% CommPtr )
      deallocate( self% CommFluxPtr )

      deallocate( self% RecvOrder )
      deallocate( self% RecvOrder0 )
! ERROR: Actually calling request free here causes MPI errors, but only durring shutdown
!      call MPI_Request_free(self% request(1), ierr)
!      call MPI_Request_free(self% request(2), ierr)
      deallocate( self% request )

      deallocate( self% NetFlux )

    endif

    deallocate( self% IncFlux )
    deallocate( self% IncFluxOld )
    deallocate( self% relError )
    deallocate( self% Converged )

    deallocate( self% NangBinList )
    deallocate( self% AngleToBin )
    deallocate( self% AngleOrder )
    deallocate( self% AngleOrder0 )

    call MPI_COMM_FREE(self% COMM_GROUP, ierr)

    return

  end subroutine CommSet_dtor

!-----------------------------------------------------------------------
  function CommSet_getCommunicationGroup(self) result(COMM_GROUP)

!    Returns the MPI communication group for this set 

!    variable declarations
     implicit none

!    passed variables
     type(CommSet), intent(in) :: self
     integer                   :: COMM_GROUP


     COMM_GROUP = self% COMM_GROUP

     return

  end function CommSet_getCommunicationGroup

!-----------------------------------------------------------------------
  function CommSet_getMessage(self,sharedID,angle) result(CommPtr)

!    Return a pointer to a communicator for boundary and angle 
!      angle      angle number
!      CommPtr    pointer to a communicator

!    variable declarations
     implicit none

!    passed variables
     type(CommSet),       intent(in) :: self
     integer,             intent(in) :: sharedID
     integer,             intent(in) :: angle
     type(Communicator),  pointer    :: CommPtr


     CommPtr => self% CommPtr(sharedID,angle)

     return

  end function CommSet_getMessage

!-----------------------------------------------------------------------
  function CommSet_getSharedFlux(self,sharedID) result(CommFluxPtr)

!    Return a pointer to a communicator for boundary and angle 
!      angle      angle number
!      CommPtr    pointer to a communicator

!    variable declarations
     implicit none

!    passed variables
     type(CommSet),         intent(in) :: self
     integer,               intent(in) :: sharedID
     type(CommunicateFlux), pointer    :: CommFluxPtr


     CommFluxPtr => self% CommFluxPtr(sharedID)

     return

  end function CommSet_getSharedFlux

!=======================================================================
! destruct communicator interface
!=======================================================================
  subroutine CommSet_dtorComm(self)

    use Size_mod

    implicit none

!   Passed variables
    type(CommSet),        intent(inout) :: self

!   Local
    type(CommunicateFlux), pointer      :: CommFlux
    integer                             :: Angle
    integer                             :: sharedID
    integer                             :: nShared

!   Free Communication Requests

    nShared = Size% ncomm

    do Angle=1,self% NumAngles
      do sharedID=1,nShared
        Comm => getMessage(self, sharedID, Angle)
        call destructBuffer(Comm)
      enddo
    enddo

    do sharedID=1,nShared
      CommFlux => getSharedFlux(self, sharedID)
      call destructFluxBuffer(CommFlux)
    enddo

    return

  end subroutine CommSet_dtorComm

!=======================================================================
! restore Communication order
!=======================================================================
  subroutine CommSet_resCommOrd(self)

    use Size_mod

    implicit none

!   Passed variables
    type(CommSet),  intent(inout) :: self


    self% NumBin         = self% NumBin0
    self% NumAnglesDyn   = self% NumAngles
    self% AngleOrder(:)  = self% AngleOrder0(:)
    self% Converged(:)   = .FALSE.

    if (Size% ncomm > 0) then
      self% RecvOrder(:,:) = self% RecvOrder0(:,:)
    endif


    return

  end subroutine CommSet_resCommOrd

!=======================================================================
! set Communication order
!=======================================================================
  subroutine CommSet_setCommOrd(self)

    use Size_mod

    implicit none

!   Passed variables
    type(CommSet),  intent(inout) :: self

!   Local
    integer                       :: angle 
    integer                       :: bin 
    integer                       :: i
    integer                       :: ishared 
    integer                       :: nsend
    integer                       :: nrecv

!   Set the send and receive orders based on a convergence test

    self% NumAnglesDyn = 0
    nsend              = 0

    do i=1,self% NumAngles
      angle = self% AngleOrder0(i)
      bin   = self% AngleToBin(angle)
      if ( .not. self% Converged(bin) ) then
        self% NumAnglesDyn                   = self% NumAnglesDyn + 1
        self% AngleOrder(self% NumAnglesDyn) = angle
        nsend                                = nsend + 1
      endif
    enddo

    do ishared=1,Size% ncomm
      nrecv = 0
      do i=1,self% NumAngles
        angle = self% RecvOrder0(i,ishared)
        bin   = self% AngleToBin(angle)
        if ( .not. self% Converged(bin) ) then
          nrecv                          = nrecv + 1
          self% RecvOrder(nrecv,ishared) = angle
        endif
      enddo
      if (nrecv /= nsend) then
        call f90fatal("Mismatch of send and receive counts in setCommOrder")
      endif
    enddo

!   Check the no-op case

    if (self% NumAnglesDyn == self% NumAngles) then
      do i=1,self% NumAngles
        if (self% AngleOrder(i) /= self% AngleOrder0(i)) then
          call f90fatal("No-op case failed for AngleOrder")
        endif
      enddo

      do ishared=1,Size% ncomm
        do i=1,self% NumAngles
          if (self% RecvOrder(i,ishared) /= self% RecvOrder0(i,ishared)) then
            call f90fatal("No-op case failed for RecvOrder")
          endif
        enddo
      enddo

    endif

    return

  end subroutine CommSet_setCommOrd


end module CommSet_mod
