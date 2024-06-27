! GreyAcceleration Module:  module for iterative acceleration 
                                                                                 
module GreyAcceleration_mod
                                                                                 
  use kind_mod

  private

! public interfaces
                                                                                             
  public construct, destruct
                                                                                 
  type, public :: GreyAcceleration 

     integer                 :: ID
     integer                 :: nGreySISweeps
     real(adqt)              :: epsGrey

!    Both GTA solvers
     real(adqt), pointer, contiguous :: GreySource(:)     => null()  ! GreySource(ncornr) 
     real(adqt), pointer, contiguous :: GreyCorrection(:) => null()  ! GreyCorrection(ncornr) 
     real(adqt), pointer, contiguous :: GreySigScat(:)    => null()  ! GreySigScat(ncornr) 
     real(adqt), pointer, contiguous :: GreySigScatVol(:) => null()  ! GreySigScatVol(ncornr)
     real(adqt), pointer, contiguous :: Chi(:,:)          => null()  ! Chi(ngr,ncornr) 
     real(adqt), pointer, contiguous :: TsaSource(:)      => null()  ! TsaSource(ncornr)

!    New GTA Solver only
     real(adqt), pointer, contiguous :: GreySigTotal(:)   => null()
     real(adqt), pointer, contiguous :: GreySigtInv(:)    => null()
     real(adqt), pointer, contiguous :: PhiInc(:)         => null()
     real(adqt), pointer, contiguous :: Sscat(:)          => null()
     real(adqt), pointer, contiguous :: Q(:)              => null()
     real(adqt), pointer, contiguous :: TT(:,:)           => null()
     real(adqt), pointer, contiguous :: Pvv(:,:)          => null()
     real(adqt), pointer, contiguous :: Tvv(:,:)          => null()

     real(adqt), pointer, contiguous :: AfpNorm(:,:,:)    => null()
     real(adqt), pointer, contiguous :: AezNorm(:,:,:)    => null()
     real(adqt), pointer, contiguous :: ANormSum(:,:)     => null()

!    Old GTA Solver only
     real(adqt), pointer, contiguous :: GreySigTotal2(:,:)   => null()
     real(adqt), pointer, contiguous :: GreySigtInv2(:,:)    => null()
     real(adqt), pointer, contiguous :: TsaPsib(:,:)         => null()
     real(adqt), pointer, contiguous :: OldGreyCorrection(:) => null()
     real(adqt), pointer, contiguous :: eps(:)               

!  The following are used in 1D for Grey Diffusion Acceleration
     real(adqt), pointer, contiguous :: GreyDiffCoef(:)   => null() ! GreyDiffCoef(ncornr) 
     real(adqt), pointer, contiguous :: GreySigEff(:)     => null() ! GreySigEff(ncornr)
     real(adqt), pointer, contiguous :: MatrixL(:,:)      => null() ! MatrixL(2,ncornr)
     real(adqt), pointer, contiguous :: MatrixU(:,:)      => null() ! MatrixU(3,ncornr)
     real(adqt), pointer, contiguous :: MatrixA(:,:)      => null() ! MatrixA(4,ncornr)

!    Misc
     character(len=3) :: label ! A string descriptor for this set.

!    For getting new GTA working:

     logical(kind=1)         :: enforceHardGTAIterMax
     logical(kind=1)         :: forceExtraOuter

  end type GreyAcceleration 

  type(GreyAcceleration), pointer, public :: GTA => null()

  interface construct
    module procedure GreyAcceleration_ctor
  end interface
 
  interface destruct
    module procedure GreyAcceleration_dtor
  end interface
 
contains
 
!=======================================================================
! construct interface
!=======================================================================

  subroutine GreyAcceleration_ctor(self)
 
    use Size_mod
    use constant_mod
    use MemoryAllocator_mod
 
    implicit none
 
!   Passed variables
    type(GreyAcceleration),  intent(inout) :: self

!   Local

    integer :: nAngGTA

    self%label = "gta"

!  Initialize a minimum convergence tolerance for the grey iteration
                                                                                             
    self% ID              = 1
!   Number of Grey source iteration sweeps made before starting BiCG
    self% nGreySISweeps   = 5 
    self% epsGrey         = 0.1_adqt

    if (Size% ndim == 2) then
      nAngGTA = 6
    elseif (Size% ndim == 3) then
      nAngGTA = 8
    endif

    call Allocator%allocate(Size%usePinnedMemory,self%label,"GreySource", self% GreySource, Size%ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"GreyCorrection", self% GreyCorrection, Size%ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"Chi", self% Chi, Size%ngr, Size%ncornr)

    self% GreyCorrection(:) = zero

!   Defaults for most solvers
    self% forceExtraOuter       = .false.
    self% enforceHardGTAIterMax = .false.

    if (Size% ndim == 1) then

!     Grey Diffusion Acceleration in 1D
      allocate( self% GreyDiffCoef(Size%ncornr) )
      allocate( self% GreySigEff(Size%ncornr) )
      allocate( self% MatrixL(2,Size%ncornr) )
      allocate( self% MatrixU(3,Size%ncornr) )
      allocate( self% MatrixA(4,Size%ncornr) )

      self% MatrixL(:,:) = zero
      self% MatrixU(:,:) = zero
      self% MatrixA(:,:) = zero

    else

!     Both GTA Solvers

      call Allocator%allocate(Size%usePinnedMemory,self%label,"GreySigScat",     self% GreySigScat,    Size%ncornr)
      call Allocator%allocate(Size%usePinnedMemory,self%label,"TsaSource",       self% TsaSource,      Size%ncornr)
      call Allocator%allocate(Size%usePinnedMemory,self%label,"GreySigScatVol",  self% GreySigScatVol, Size%ncornr )

      if (Size% useNewGTASolver) then

!       New GTA Solver only

        call Allocator%allocate(Size%usePinnedMemory,self%label,"GreySigTotal", self% GreySigTotal, Size%ncornr)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"GreySigtInv",  self% GreySigtInv,  Size%ncornr)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"PhiInc",       self% PhiInc,       Size%ncornr)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"Sscat",        self% Sscat,        Size%ncornr)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"Q",            self% Q,            Size%ncornr)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"TT",           self% TT,           Size%maxCorner,Size%ncornr)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"Pvv",          self% Pvv,          Size%maxCorner,Size%ncornr)

        call Allocator%allocate(Size%usePinnedMemory,self%label,"AfpNorm",  self% AfpNorm,  Size%maxcf, Size%ncornr, nAngGTA)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"AezNorm",  self% AezNorm,  Size%maxcf, Size%ncornr, nAngGTA)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"ANormSum", self% ANormSum, Size%ncornr, nAngGTA)

        if (Size% ndim == 2) then
          call Allocator%allocate(Size%usePinnedMemory,self%label,"Tvv", self% Tvv, Size%maxCorner,Size%ncornr)
        endif

        self% enforceHardGTAIterMax = .true. ! This is needed so that BiCGSTAB doesn't accumulate roundoff errors that lead to a divergent iteration scheme

      else

!       Old GTA Solver only
        allocate( self% GreySigTotal2(Size%ncornr,2) )
        allocate( self% GreySigtInv2(Size%ncornr,2) )
        allocate( self% TsaPsib(Size%nbelem,Size%nangGTA) )
        allocate( self% OldGreyCorrection(Size%ncornr) )
        allocate( self% eps(Size%ncornr) )

      endif

    endif

    return
 
  end subroutine GreyAcceleration_ctor

!=======================================================================
! destruct interface
!=======================================================================
                                                            
  subroutine GreyAcceleration_dtor(self)
                                      
    use Size_mod
    use MemoryAllocator_mod
                
    implicit none
                 
!   Passed variables
                    
    type(GreyAcceleration),  intent(inout) :: self
                                            
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"GreySource", self% GreySource)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"GreyCorrection", self% GreyCorrection )
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"Chi", self% Chi )

    if (Size% ndim == 1) then

!     Grey Diffusion Acceleration in 1D
      deallocate( self% GreyDiffCoef )
      deallocate( self% GreySigEff )
      deallocate( self% MatrixL )
      deallocate( self% MatrixU )
      deallocate( self% MatrixA )

    else

!     Both GTA Solvers

      call Allocator%deallocate(Size%usePinnedMemory,self%label,"GreySigScat", self% GreySigScat)
      call Allocator%deallocate(Size%usePinnedMemory,self%label,"GreySigScatVol", self% GreySigScatVol)
      call Allocator%deallocate(Size%usePinnedMemory,self%label,"TsaSource",   self% TsaSource)

      if (Size% useNewGTASolver) then

!       New GTA Solver only

        call Allocator%deallocate(Size%usePinnedMemory,self%label,"GreySigTotal", self% GreySigTotal)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"GreySigtInv",  self% GreySigtInv)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"PhiInc",       self% PhiInc)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"Sscat",        self% Sscat)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"Q",            self% Q)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"TT",           self% TT)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"Pvv",          self% Pvv)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"AfpNorm",      self% AfpNorm)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"AezNorm",      self% AezNorm)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"ANormSum",     self% ANormSum)

        if (Size% ndim == 2) then
          call Allocator%deallocate(Size%usePinnedMemory,self%label,"Tvv",        self% Tvv)
        endif

      else

!       Old GTA Solver only
        deallocate( self% GreySigTotal2 )
        deallocate( self% GreySigtInv2 )
        deallocate( self% TsaPsib )
        deallocate( self% OldGreyCorrection )
        deallocate( self% eps )

      endif

    endif


    return
          
  end subroutine GreyAcceleration_dtor

end module GreyAcceleration_mod
