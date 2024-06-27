! Size Module:  mesh-dependent parameters used by Sn
#include "macros.h"                                                                                
module Size_mod 

  use flags_mod
  use kind_mod
  use constant_mod
  use Options_mod

  private

! public interfaces

  public construct
  public getGeometryFactor
  public getGPUStatus
  public getCUDASweepStatus
                                                                                 
  type, public :: MeshSize
     integer          :: myRankInGroup   ! my rank in my communication group 
     integer          :: myThreadID      ! my thread ID
     integer          :: nThreads        ! number of threads per MPI process
     integer          :: nzones          ! number of zones
     integer          :: ncornr          ! number of corners
     integer          :: nSides
     integer          :: nbelem          ! number of boundary elements
     integer          :: nSurfElem       ! boundary + interface elements
     integer          :: maxcf           ! maximum number of zone faces a corner touches
     integer          :: maxCorner       ! maximum number of corners in a zone
     integer          :: maxFaces        ! maximum number of zone-faces
     integer          :: maxSides        ! maximum number of sides
     integer          :: ncomm           ! number of processors to communicate with

     integer          :: ndim            ! number of spatial dimensions
     integer          :: ngr             ! number of energy groups
     integer          :: nangGTA         ! number of angles used for GTA sweeps

     integer          :: nBCITabTaus     ! # of Boltz Compton intgrl table tau vals
     integer          :: nBCITabG2Gs     ! # of Boltz Compton intgrl table group-to-group vals

     integer          :: functionRNLTE

     real(adqt)       :: tfloor              ! temperature floor
     real(adqt)       :: tr4floor            ! tfloor**4
     real(adqt)       :: wtiso               ! isotropic normalization factor
     real(adqt)       :: radForceMultiplier  ! radiation force multiplier
     real(adqt)       :: betaNLTE
     real(adqt)       :: gammaNLTE
     real(adqt)       :: tau                 ! reciprocal of timestep*SpeedOfLight
     real(adqt)       :: geometryFactor      ! multiplier need for a 3D volume
     real(adqt)       :: MatCoupTimeCycle
     real(adqt)       :: SweepTimeCycle
     real(adqt)       :: GPUSweepTimeCycle
     real(adqt)       :: GTATimeCycle
     real(adqt)       :: RadtrTimeCycle
     real(adqt)       :: InitTimeCycle
     real(adqt)       :: FinalTimeCycle
     real(adqt)       :: MatCoupTimeTotal
     real(adqt)       :: SweepTimeTotal
     real(adqt)       :: GPUSweepTimeTotal
     real(adqt)       :: GTATimeTotal
     real(adqt)       :: RadtrTimeTotal
     real(adqt)       :: InitTimeTotal
     real(adqt)       :: FinalTimeTotal

     logical (kind=1) :: DopplerShiftOn        ! Doppler shift control
     logical (kind=1) :: useNewNonLinearSolver ! Non Linear solver control 
     logical (kind=1) :: useNewGTASolver
     logical (kind=1) :: usePWLD               ! use PWLD spatial differencing
     logical (kind=1) :: useSurfaceMassLumping ! surface and mass lumping for PWLD 
     logical (kind=1) :: useGPU                ! offload computations to the GPU
     logical          :: usePinnedMemory       ! use page-locked cpu memory for
                                               ! data that will be transferred to GPU.
     logical (kind=1) :: useCUDASweep          ! use CUDA sweep on the GPU
     logical (kind=1) :: useCUDASolver         ! use CUDA solver on the GPU

     integer          :: zoneBatchSize
     integer          :: nConcurrentBatches
     integer          :: igeom                 ! geometry flag

  end type MeshSize 

  type(MeshSize), pointer, public :: Size

  interface construct
    module procedure Size_ctor
  end interface

  interface getGeometryFactor
    module procedure Size_getGeometryFactor
  end interface

  interface getGPUStatus
    module procedure Size_getGPUStatus
  end interface

  interface getCUDASweepStatus
    module procedure Size_getCUDASweepStatus
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine Size_ctor(self, myRankInGroup, nzones, ncornr, nSides,   &
                             nbelem, maxcf, maxCorner, ncomm,         &
                             ndim, ngr,                               &
                             functionRNLTE, tfloor,                   &
                             radForceMultiplier, betaNLTE, gammaNLTE, &
                             DopplerShiftOn,                          &
                             useNewNonLinearSolver, useNewGTASolver,  &
                             usePWLD, useSurfaceMassLumping,          &
                             useGPU, useCUDASweep, useCUDASolver,     &
                             zoneBatchSize, nConcurrentBatches,       &
                             igeom) 

    implicit none
                                                                                         
!   Passed variables

    type(MeshSize),    intent(inout)       :: self

    integer, optional, intent(in)          :: myRankInGroup 
    integer, optional, intent(in)          :: nzones
    integer, optional, intent(in)          :: ncornr
    integer, optional, intent(in)          :: nSides
    integer, optional, intent(in)          :: nbelem
    integer, optional, intent(in)          :: maxcf
    integer, optional, intent(in)          :: maxCorner
    integer, optional, intent(in)          :: ncomm
    integer, optional, intent(in)          :: ndim
    integer, optional, intent(in)          :: ngr
    integer, optional, intent(in)          :: functionRNLTE

    real(adqt), optional, intent(in)       :: tfloor
    real(adqt), optional, intent(in)       :: radForceMultiplier
    real(adqt), optional, intent(in)       :: betaNLTE
    real(adqt), optional, intent(in)       :: gammaNLTE

    logical (kind=1), optional, intent(in) :: DopplerShiftOn
    logical (kind=1), optional, intent(in) :: useNewNonLinearSolver
    logical (kind=1), optional, intent(in) :: useNewGTASolver
    logical (kind=1), optional, intent(in) :: usePWLD
    logical (kind=1), optional, intent(in) :: useSurfaceMassLumping
    logical (kind=1), optional, intent(in) :: useGPU
    logical (kind=1), optional, intent(in) :: useCUDASweep
    logical (kind=1), optional, intent(in) :: useCUDASolver

    integer, optional, intent(in)          :: zoneBatchSize
    integer, optional, intent(in)          :: nConcurrentBatches

    integer, optional, intent(in) :: igeom

!   This is a 'hack'.
!   Teton doesn't have an overall 'initialization' subroutine that we can
!   add a one-time initialization of classes to.  Adding it here in the 'Size'
!   module for now, as we've been parking things at least related to runtime
!   options here. --Aaron
    call Options%initialize()
    call Options%check()
   
!   Problem Size Parameters
    self% myRankInGroup      = myRankInGroup 
    self% myThreadID         = 0
    self% nThreads           = 1
    self% nzones             = nzones
    self% ncornr             = ncornr
    self% nSides             = nSides
    self% nbelem             = nbelem
    self% nSurfElem          = nbelem
    self% maxcf              = maxcf
    self% maxCorner          = maxCorner 
! In 1D, 2D maxFaces = maxCorner; In 3D it's an overestimate
    self% maxFaces           = maxCorner 
    self% maxSides           = 1
    self% ncomm              = ncomm
    self% ndim               = ndim
    self% ngr                = ngr
    self% nangGTA            = 8
    self% functionRNLTE      = functionRNLTE

    self% nBCITabTaus        = 21
    self% nBCITabG2Gs        = ngr*(ngr-1)/2

    self% tfloor             = tfloor
    self% tr4floor           = tfloor*tfloor*tfloor*tfloor
    self% radForceMultiplier = radForceMultiplier
    self% betaNLTE           = betaNLTE
    self% gammaNLTE          = gammaNLTE
    self% tau                = zero

    self% MatCoupTimeCycle   = zero
    self% SweepTimeCycle     = zero
    self% GPUSweepTimeCycle  = zero
    self% GTATimeCycle       = zero
    self% RadtrTimeCycle     = zero
    self% InitTimeCycle      = zero
    self% FinalTimeCycle     = zero
    self% MatCoupTimeTotal   = zero
    self% SweepTimeTotal     = zero
    self% GPUSweepTimeTotal  = zero
    self% GTATimeTotal       = zero
    self% RadtrTimeTotal     = zero
    self% InitTimeTotal      = zero
    self% FinalTimeTotal     = zero

!  Check consistency of dopper shift flag

    if (self% ngr == 1) then
      self% DopplerShiftOn = .FALSE.
    else
      self% DopplerShiftOn = DopplerShiftOn
    endif

    self% useNewNonLinearSolver = useNewNonLinearSolver
    self% usePWLD               = usePWLD
    self% useSurfaceMassLumping = useSurfaceMassLumping

    if (ndim == 1) then
    ! For 1D, teton uses the GDA Solver, not GTA Solver.  Set this to False.
       self% useNewGTASolver = .FALSE.
    ! GPU kernels do not support 1d meshes.
       self% useGPU = .FALSE.
    else
       self% useNewGTASolver       = useNewGTASolver
       self% useGPU = useGPU
    endif

! Pinned memory improves CPU<->GPU memory transfer performance, so default to
! preferring to use Umpire for pinned memory allocations on CPU.
#if defined(TETON_ENABLE_UMPIRE)
    self% usePinnedMemory       = self%useGPU
#else
    self% usePinnedMemory       = .FALSE.
#endif

    if (self% useGPU .AND. .NOT. self% usePinnedMemory .AND. myRankInGroup == 0) then
       print *, "TETON WARNING: Detected that GPU kernels are enabled, but UMPIRE support is not enabled.  This may impact CPU<->GPU memory transfer performance."
    endif

    self% useCUDASolver         = useCUDASolver
    self% zoneBatchSize         = zoneBatchSize
    self% nConcurrentBatches    = nConcurrentBatches
    self% igeom                 = igeom

    ! Internal re-setting of useCUDASweep to avoid changing
    ! host code interfaces. useGPU must be FALSE in order to
    ! activate CUDA Sweep.
    if ( useCUDASweep ) then
      if ( .not. self%useGPU ) then
        if ( myRankInGroup == 0) then
           print *, "ACTIVATING CUDA SWEEP"
        endif
        self% useCUDASweep          = useCUDASweep
      else
        call f90fatal( "EXITING: CUDA SWEEP NOT ACTIVATED, please set useGPU=false" )
      endif
    else
      self% useCUDASweep          = .FALSE.
    endif

!  Angle-dependent variables are (by convention) per unit
!  cosine in 1D, per 2pi steradians in 2D and per 4pi
!  steradians in 3D.

    select case (self% ndim)
      case (1)
        self% wtiso = half

        if (self% igeom == geometry_slab) then
          self% geometryFactor = one
        elseif (self% igeom == geometry_sphere) then
          self% geometryFactor = four*pi
        elseif (self% igeom == geometry_cylinder) then
          self% geometryFactor = two*pi 
        endif
      case (2)
        self% wtiso          = one/(two*pi)
        self% geometryFactor = two*pi
      case (3)
        self% wtiso          = one/(four*pi)
        self% geometryFactor = one
    end select

    return
   
  end subroutine Size_ctor

!=======================================================================
! getGeometryFactor interface
!=======================================================================

  function Size_getGeometryFactor(self) result(geometryFactor)

!  Returns a geometry factor for the problem geometry

!    variable declarations
     implicit none

!    passed variables
     type(MeshSize), intent(in)   :: self
     real(adqt)                   :: geometryFactor 

     geometryFactor = self% geometryFactor 

     return

  end function Size_getGeometryFactor

!=======================================================================
! getGPUStatus interface
!=======================================================================

  function Size_getGPUStatus(self) result(useGPU)

!  Returns the GPU control 

!    variable declarations
     implicit none

!    passed variables
     type(MeshSize), intent(in)   :: self
     logical(kind=1)              :: useGPU

     useGPU = self% useGPU

     return

  end function Size_getGPUStatus
                                                                                 
!=======================================================================
! getCUDASweepStatus interface
!=======================================================================

  function Size_getCUDASweepStatus(self) result(useCUDASweep)

!  Returns the CUDA sweep control 

!    variable declarations
     implicit none

!    passed variables
     type(MeshSize), intent(in)   :: self
     logical(kind=1)              :: useCUDASweep

     useCUDASweep = self% useCUDASweep

     return

  end function Size_getCUDASweepStatus
                                                                                 
                                                      
end module Size_mod

