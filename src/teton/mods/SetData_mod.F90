#include "macros.h"
! SetData Module:  Contains data structures for a group-angle set 

module SetData_mod 

  use kind_mod
  use constant_mod
  use Quadrature_mod
  use, intrinsic :: iso_c_binding, only : c_double

  private

  type, public :: SetData 

     integer                  :: SetID              ! group-angle set ID
     integer                  :: groupSetID
     integer                  :: angleSetID
     integer                  :: QuadID
     integer                  :: Groups             ! number of energy groups 
     integer                  :: NumAngles          ! number of angles
     integer                  :: NumAnglesDyn 
     integer                  :: nPolarAngles
     integer                  :: Order
     integer                  :: nZones
     integer                  :: nCorner
     integer                  :: maxCorner
     integer                  :: nbelem
     integer                  :: g0                 ! offset to global group
     integer                  :: angle0

     integer                  :: NumBin0            !
     integer                  :: NumBin             ! number of angle bins for

!    Quadrature related
     integer,         pointer, contiguous :: AngleOrder(:) => null()
     integer,         pointer, contiguous :: PolarAngle(:) => null()
     real(adqt),      pointer, contiguous :: PolarAngleList(:) => null()

!    For transport Sweeps

     real(C_DOUBLE),  pointer, contiguous :: Psi(:,:,:) => null()
     real(adqt),      pointer, contiguous :: Psi1(:,:) => null()
     real(adqt),      pointer, contiguous :: PsiM(:,:) => null()
     real(C_DOUBLE),  pointer, contiguous :: PsiB(:,:,:) => null()   ! boundary intensities
     real(adqt),      pointer, contiguous :: Phi(:,:) => null()
     real(adqt),      pointer, contiguous :: cyclePsi(:,:) => null()
     real(adqt),      pointer, contiguous :: Q(:,:,:) => null()
     real(adqt),      pointer, contiguous :: S(:,:,:) => null()

!    For GTA Sweeps

     real(adqt),      pointer, contiguous :: tPsi(:) => null()
     real(adqt),      pointer, contiguous :: tPsiM(:) => null()
     real(adqt),      pointer, contiguous :: tPhi(:) => null()
     real(adqt),      pointer, contiguous :: pInc(:) => null()
     real(adqt),      pointer, contiguous :: tInc(:) => null()
     real(adqt),      pointer, contiguous :: src(:)  => null()

!    For Spectral Tallies

     real(adqt),      pointer, contiguous :: RadPowerEscape(:) => null()
     real(adqt),      pointer, contiguous :: RadPowerIncident(:) => null()
     real(adqt),      pointer, contiguous :: PolarSectorPowerEscape(:,:) => null()

!    Misc
     character(len=19) :: label ! A string descriptor for this set.

     contains
       procedure, public :: construct => SetData_ctor
       procedure, public :: destruct => SetData_dtor
       procedure, public :: destructCyclePsi => SetData_dtorCyclePsi

  end type SetData 

contains

!=======================================================================
! construct interface
!=======================================================================
                                                                                   
  subroutine SetData_ctor(self,        &
                           SetID,       &
                           groupSetID,  &
                           angleSetID,  &
                           QuadID,      &
                           Groups,      &
                           NumAngles,   &
                           g0,          &
                           angle0,      &
                           nZones,      &
                           nCorner,     &
                           QuadPtr,     &
                           GTASet,      &
                           fromRestart )

    use Size_mod
    use constant_mod
    use MemoryAllocator_mod
    use, intrinsic :: iso_c_binding, only : c_size_t

    implicit none

!   Passed variables

    class(SetData), intent(inout)        :: self

    integer, intent(in)                  :: SetID
    integer, intent(in)                  :: groupSetID
    integer, intent(in)                  :: angleSetID
    integer, intent(in)                  :: QuadID
    integer, intent(in)                  :: Groups       
    integer, intent(in)                  :: NumAngles
    integer, intent(in)                  :: g0
    integer, intent(in)                  :: angle0
    integer, intent(in)                  :: nZones
    integer, intent(in)                  :: nCorner

    type(Quadrature), target, intent(in) :: QuadPtr

    logical (kind=1), intent(in)         :: GTASet
    logical (kind=1), intent(in)         :: fromRestart

!   Local

    integer :: n
    integer :: ndim
    integer :: angle
    integer :: nLevels
    integer :: maxZones
    logical :: usePinnedMemory

!   Set Properties

    ndim               =  Size% ndim

!   IMPORTANT: "maxZones" is intended to estimate the maximum number of zones
!   in any sweep hyperplane. It is sized to be conservative (i.e. it is much
!   bigger than it needs to be). It is sized here to avoid gpu memory errors.
!   It should really be sized in the hyperplane constructor. (PFN)

    ! This used to be nzones/25, but the max hyperplane is very problem
    ! dependent and that was insufficient for some problem configurations,
    ! resulting in array overruns on the gpu.
    ! Should make another attempt to move this to the hyperplane constructor,
    ! per Paul's comment above.
    maxZones           =  MAX(100, (Size% nzones / 3) )

    self% SetID        =  SetID
    self% groupSetID   =  groupSetID
    self% angleSetID   =  angleSetID
    self% QuadID       =  QuadID
    self% Groups       =  Groups 
    self% NumAngles    =  NumAngles
    self% NumAnglesDyn =  NumAngles
    self% nPolarAngles =  QuadPtr% nPolarAngles
    self% Order        =  QuadPtr% Order
    self% nZones       =  nZones
    self% nCorner      =  nCorner 
    self% maxCorner    =  Size% maxCorner
    self% nbelem       =  Size% nbelem
    self% g0           =  g0
    self% angle0       =  angle0

    usePinnedMemory = Size%useGPU

    ! Make sure loadTetonVariables.F90 stay consistent with this logic that sets
    ! the label.  It also creates this label during a restart load.
    write(self%label,'(I0.3)') setID
    self%label = "phase_space_set_"//self%label

    if ( .not. GTASet ) then

      if ( .not. fromRestart ) then
        call Allocator%allocate(usePinnedMemory,self%label,"Psi", self% Psi, Groups,nCorner,NumAngles)
        self% Psi(:,:,:) = zero
      endif

      call Allocator%allocate(usePinnedMemory,self%label,"PsiB", self% PsiB, Groups,Size% nbelem,NumAngles)
      self% PsiB(:,:,:) = zero

      call Allocator%allocate(usePinnedMemory,self%label,"Psi1", self% Psi1, Groups,nCorner+Size% nbelem)
      self% Psi1(:,:)   = zero

      call Allocator%allocate(usePinnedMemory,self%label,"Q",    self% Q,    Groups,Size% maxCorner,maxZones)
      self% Q(:,:,:)    = zero

      call Allocator%allocate(usePinnedMemory,self%label,"S",    self% S,    Groups,Size% maxCorner,maxZones)
      self% S(:,:,:)    = zero

      if (ndim == 2) then
        call Allocator%allocate(usePinnedMemory,self%label,"PsiM",  self% PsiM,  Groups,nCorner)
        self% PsiM(:,:) = zero
      endif
      ! Using either CPU or CUDA Sweep
      if (Size%useGPU .NEQV. .TRUE.) then
        if (Size%useCUDASweep .EQV. .TRUE.) then
          ! CUDA Sweep solver needs it page-locked for optimal transfer performance.
          call Allocator%allocate(.TRUE.,self%label,"Phi",  self% Phi, Groups,nCorner )
        else
          ! CPU Sweep does not need this page-locked.
          allocate( self% Phi(Groups,nCorner) )
        endif
      endif

    else

      call Allocator%allocate(usePinnedMemory,self%label,"tPsi", self% tPsi, nCorner+Size% nbelem)
      call Allocator%allocate(usePinnedMemory,self%label,"pInc", self% pInc, nCorner)
      call Allocator%allocate(usePinnedMemory,self%label,"src" , self% src,  nCorner)

      self% tPsi(:) = zero
      self% pInc(:) = zero
      self% src(:)  = zero

      if (ndim == 2) then
        call Allocator%allocate(usePinnedMemory,self%label,"tInc",  self% tInc, nCorner)
        call Allocator%allocate(usePinnedMemory,self%label,"tPsiM", self% tPsiM, nCorner)

        self% tInc(:)  = zero
        self% tPsiM(:) = zero
      endif

      if ( .not. Size% useNewGTASolver ) then
        allocate( self% tPhi(nCorner) )
        self% tPhi(:) = zero
      endif

    endif

    allocate( self% AngleOrder(self% NumAngles) )
    allocate( self% PolarAngle(self% NumAngles) )
    allocate( self% PolarAngleList(self% nPolarAngles) )

    self% PolarAngleList(:) = QuadPtr% PolarAngleList(:)

!   For Spectral Tallies - 

    allocate( self% RadPowerEscape(self% Groups) )
    allocate( self% RadPowerIncident(self% Groups) )
    allocate( self% PolarSectorPowerEscape(self% Groups, self% nPolarAngles) )

!   Quadrature Related

    nLevels = 0

    do n=1,NumAngles
      angle               = QuadPtr% angleList(angle0+n)
      self% PolarAngle(n) = QuadPtr% PolarAngle(angle)
      self% AngleOrder(n) = 0

      if ( QuadPtr% StartingDirection(angle) ) then
        nLevels = nLevels + 1
      endif
    enddo

    if (ndim == 2) then
      self% NumBin  = nLevels
      self% NumBin0 = nLevels
    else
      self% NumBin  = NumAngles
      self% NumBin0 = NumAngles
    endif

    return

  end subroutine SetData_ctor

!=======================================================================
! destruct interface
!=======================================================================
  subroutine SetData_dtor(self, GTASet)

    use Size_mod
    use MemoryAllocator_mod

    implicit none

!   Passed variables
                                                                                     
    class(SetData),   intent(inout) :: self
    logical(kind=1), intent(in)     :: GTASet

    logical :: usePinnedMemory
    usePinnedMemory = Size%useGPU

    if ( .not. GTASet ) then
      call Allocator%deallocate(usePinnedMemory,self%label,"Psi",  self% Psi  )
      call Allocator%deallocate(usePinnedMemory,self%label,"PsiB", self% PsiB )
      call Allocator%deallocate(usePinnedMemory,self%label,"Psi1", self% Psi1 )
      call Allocator%deallocate(usePinnedMemory,self%label,"Q",    self% Q)
      call Allocator%deallocate(usePinnedMemory,self%label,"S",    self% S)

      if ( Size% ndim == 2) then
        call Allocator%deallocate(usePinnedMemory,self%label,"PsiM", self% PsiM )
      endif

      ! Using either CPU or CUDA Sweep
      if (Size%useGPU .NEQV. .TRUE.) then
        if (Size%useCUDASweep .EQV. .TRUE.) then
          ! CUDA Sweep solver needs it page-locked for optimal transfer performance.
          call Allocator%deallocate(usePinnedMemory,self%label,"Phi",  self% Phi  )
        else
          ! CPU Sweep does not need this page-locked.
          deallocate( self% Phi )
        endif
      endif

    else

      call Allocator%deallocate(usePinnedMemory,self%label,"tPsi", self% tPsi )
      call Allocator%deallocate(usePinnedMemory,self%label,"pInc", self% pInc )
      call Allocator%deallocate(usePinnedMemory,self%label,"src" , self% src )

      if (Size% ndim == 2) then
        call Allocator%deallocate(usePinnedMemory,self%label,"tInc", self% tInc )
        call Allocator%deallocate(usePinnedMemory,self%label,"tPsiM", self% tPsiM )
      endif


      if ( .not. Size% useNewGTASolver ) then
        deallocate( self% tPhi )
      endif

    endif

    deallocate( self% AngleOrder )
    deallocate( self% PolarAngle ) 
    deallocate( self% PolarAngleList )

    deallocate( self% RadPowerEscape )
    deallocate( self% RadPowerIncident )
    deallocate( self% PolarSectorPowerEscape )

    return

  end subroutine SetData_dtor

!=======================================================================
! destruct cycle Psi    
!=======================================================================
  subroutine SetData_dtorCyclePsi(self)

    use MemoryAllocator_mod
    use Size_mod

    implicit none

!   Passed variables
    class(SetData),    intent(inout) :: self

    logical :: usePinnedMemory
    usePinnedMemory = Size%useGPU

!   Release Memory 

    call Allocator%deallocate(usePinnedMemory,self%label,"cyclePsi", self% cyclePsi)

    return

  end subroutine SetData_dtorCyclePsi

end module SetData_mod
