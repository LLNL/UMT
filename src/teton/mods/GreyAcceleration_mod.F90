! GreyAcceleration Module:  module for iterative acceleration 
                                                                                 
module GreyAcceleration_mod
                                                                                 
  use kind_mod

  private

! public interfaces
                                                                                             
  public construct, destruct
                                                                                 
  type, public :: GreyAcceleration 

     integer                 :: ID
     integer                 :: nGreySweepIters
     real(adqt)              :: epsGrey

!    Both GTA solvers
     real(adqt), pointer, contiguous :: GreySource(:)     => null()  ! GreySource(ncornr) 
     real(adqt), pointer, contiguous :: GreyCorrection(:) => null()  ! GreyCorrection(ncornr) 
     real(adqt), pointer, contiguous :: GreySigScat(:)    => null()  ! GreySigScat(ncornr) 
     real(adqt), pointer, contiguous :: GreySigScatVol(:) => null()  ! GreySigScatVol(ncornr)
     real(adqt), pointer, contiguous :: Chi(:,:)          => null()  ! Chi(ngr,ncornr) 
     real(adqt), pointer, contiguous :: TsaSource(:)      => null()  ! TsaSource(ncornr)

     real(adqt), pointer, contiguous :: CGDirectionB(:,:) => null()  ! CGDirectionB(nbelem,nangGTA)
     real(adqt), pointer, contiguous :: CGResidualB(:,:)  => null()  ! CGResidualB(nbelem,nangGTA)
     real(adqt), pointer, contiguous :: CGActionB(:,:)    => null()  ! CGActionB(nbelem,nangGTA)
     real(adqt), pointer, contiguous :: CGActionSB(:,:)   => null()  ! CGActionSB(nbelem,nangGTA)

!    New GTA Solver only
     real(adqt), pointer, contiguous :: GreySigTotal(:)   => null()
     real(adqt), pointer, contiguous :: GreySigtInv(:)    => null()
     real(adqt), pointer, contiguous :: PhiInc(:)         => null()
     real(adqt), pointer, contiguous :: Q(:)              => null()
     real(adqt), pointer, contiguous :: TT(:,:)           => null()
     real(adqt), pointer, contiguous :: Pvv(:,:)          => null()
     real(adqt), pointer, contiguous :: Tvv(:,:)          => null()

!    Old GTA Solver only
     real(adqt), pointer, contiguous :: GreySigTotal2(:,:)   => null()
     real(adqt), pointer, contiguous :: GreySigtInv2(:,:)    => null()
     real(adqt), pointer, contiguous :: TsaCorrection(:)     => null()
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

    logical :: usePinnedMemory
    usePinnedMemory = Size%useGPU

    self%label = "gta"

!  Initialize a minimum convergence tolerance for the grey iteration
                                                                                             
    self% ID              = 1
    self% nGreySweepIters = 2
    self% epsGrey         = 0.1_adqt

    call Allocator%allocate(usePinnedMemory,self%label,"GreySource", self% GreySource, Size%ncornr)

    allocate( self% GreyCorrection(Size%ncornr) )
    allocate( self% Chi(Size%ngr,Size%ncornr) )

    self% GreyCorrection(:) = zero

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

      call Allocator%allocate(usePinnedMemory,self%label,"GreySigScat", self% GreySigScat, Size%ncornr)
      call Allocator%allocate(usePinnedMemory,self%label,"TsaSource",   self% TsaSource,   Size%ncornr)

      allocate( self% GreySigScatVol(Size%ncornr) )
      allocate( self% CGDirectionB(Size%nbelem,Size%nangGTA) )
      allocate( self% CGResidualB(Size%nbelem,Size%nangGTA) )
      allocate( self% CGActionB(Size%nbelem,Size%nangGTA) )
      allocate( self% CGActionSB(Size%nbelem,Size%nangGTA) )

      if (Size% useNewGTASolver) then

!       New GTA Solver only

        call Allocator%allocate(usePinnedMemory,self%label,"GreySigTotal", self% GreySigTotal, Size%ncornr)
        call Allocator%allocate(usePinnedMemory,self%label,"GreySigtInv",  self% GreySigtInv,  Size%ncornr)
        call Allocator%allocate(usePinnedMemory,self%label,"PhiInc",       self% PhiInc,       Size%ncornr)
        call Allocator%allocate(usePinnedMemory,self%label,"Q",            self% Q,            Size%ncornr)
        call Allocator%allocate(usePinnedMemory,self%label,"TT",           self% TT,           Size%maxCorner,Size%ncornr)
        call Allocator%allocate(usePinnedMemory,self%label,"Pvv",          self% Pvv,          Size%maxCorner,Size%ncornr)

        if (Size% ndim == 2) then
          call Allocator%allocate(usePinnedMemory,self%label,"Tvv", self% Tvv, Size%maxCorner,Size%ncornr)
        endif

      else

!       Old GTA Solver only
        allocate( self% GreySigTotal2(Size%ncornr,2) )
        allocate( self% GreySigtInv2(Size%ncornr,2) )
        allocate( self% TsaCorrection(Size%ncornr) )
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
                                            
    logical :: usePinnedMemory
    usePinnedMemory = Size%useGPU

    call Allocator%deallocate(usePinnedMemory,self%label,"GreySource", self% GreySource)

    deallocate( self% GreyCorrection )
    deallocate( self% Chi )

    if (Size% ndim == 1) then

!     Grey Diffusion Acceleration in 1D
      deallocate( self% GreyDiffCoef )
      deallocate( self% GreySigEff )
      deallocate( self% MatrixL )
      deallocate( self% MatrixU )
      deallocate( self% MatrixA )

    else

!     Both GTA Solvers

      call Allocator%deallocate(usePinnedMemory,self%label,"GreySigScat", self% GreySigScat)
      call Allocator%deallocate(usePinnedMemory,self%label,"TsaSource",   self% TsaSource)

      deallocate( self% GreySigScatVol )
      deallocate( self% CGDirectionB )
      deallocate( self% CGResidualB )
      deallocate( self% CGActionB )
      deallocate( self% CGActionSB )

      if (Size% useNewGTASolver) then

!       New GTA Solver only

        call Allocator%deallocate(usePinnedMemory,self%label,"GreySigTotal", self% GreySigTotal)
        call Allocator%deallocate(usePinnedMemory,self%label,"GreySigtInv",  self% GreySigtInv)
        call Allocator%deallocate(usePinnedMemory,self%label,"PhiInc",       self% PhiInc)
        call Allocator%deallocate(usePinnedMemory,self%label,"Q",            self% Q)
        call Allocator%deallocate(usePinnedMemory,self%label,"TT",           self% TT)
        call Allocator%deallocate(usePinnedMemory,self%label,"Pvv",          self% Pvv)

        if (Size% ndim == 2) then
          call Allocator%deallocate(usePinnedMemory,self%label,"Tvv",        self% Tvv)
        endif

      else

!       Old GTA Solver only
        deallocate( self% GreySigTotal2 )
        deallocate( self% GreySigtInv2 )
        deallocate( self% TsaCorrection )
        deallocate( self% TsaPsib )
        deallocate( self% OldGreyCorrection )
        deallocate( self% eps )

      endif

    endif


    return
          
  end subroutine GreyAcceleration_dtor

end module GreyAcceleration_mod
