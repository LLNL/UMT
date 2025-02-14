! Material Module:  material properties 
                                                                                 
module Material_mod
  use kind_mod, only : adqt
  implicit none

  private

  type, public :: Material

     real(adqt),      pointer, contiguous :: cve(:)              => null()  ! "effective" specific heat 
     real(adqt),      pointer, contiguous :: rho(:)              => null()  ! Material Density 
     real(adqt),      pointer, contiguous :: nez(:)              => null()  ! Electron Density 
     real(adqt),      pointer, contiguous :: SMatEff(:)          => null()  ! Source to the material 
     real(adqt),      pointer, contiguous :: denez(:)            => null()  ! denez
     real(adqt),      pointer, contiguous :: Trz(:)              => null()  ! trz
     real(adqt),      pointer, contiguous :: Tez(:)              => null()  ! tez
     real(adqt),      pointer, contiguous :: Trzn(:)             => null()  ! trzn
     real(adqt),      pointer, contiguous :: Tezn(:)             => null()  ! tezn
     real(adqt),      pointer, contiguous :: Tezold(:)           => null()  ! tezold
     real(adqt),      pointer, contiguous :: EnergyDensityOld(:) => null()
     real(adqt),      pointer, contiguous :: stimComptonMult(:)  => null()  ! stimulated Compton mult
     logical(kind=1), pointer, contiguous :: isVoid(:)           => null()  ! void zone
     logical(kind=1), pointer, contiguous :: nonLTE(:)           => null()  ! non-local thermodynamic equilibrium

!    Thread-Safe Tallies
     integer,    pointer, contiguous :: nonLinearIterations(:)   => null()
     real(adqt), pointer, contiguous :: PowerEmitted(:)          => null()
     real(adqt), pointer, contiguous :: PowerCompton(:)          => null()

     real(adqt), pointer, contiguous :: denec(:)                 => null()  ! denec(nCorner)
     real(adqt), pointer, contiguous :: Tec(:)                   => null()  ! tec(nCorner)
     real(adqt), pointer, contiguous :: Tecn(:)                  => null()  ! tecn(nCorner)

     real(adqt), pointer, contiguous :: sigA(:,:)                => null()  ! absorption coefficient 
     real(adqt), pointer, contiguous :: sigS(:,:)                => null()  ! scattering coefficient 
     real(adqt), pointer, contiguous :: Eta(:)                   => null()  ! Eta(nCorner)

!    For nonLTE
     real(adqt), pointer, contiguous :: emis(:,:)                => null()  ! emissivity
     real(adqt), pointer, contiguous :: demisdT(:,:)             => null()

!    Non-linear Solver
     real(adqt), pointer, contiguous :: EmissionRate(:,:)        => null()

!    Misc
     character(len=8) :: label

     contains

       procedure, public :: construct => Material_ctor
       procedure, public :: destruct  => Material_dtor

  end type Material

  type(Material), pointer, public :: Mat => null()

contains
 
!=======================================================================
! construct interface
!=======================================================================

  subroutine Material_ctor(self, nonLTE, fromRestart)
    use Size_mod, only : Size
    use constant_mod, only : zero
    use, intrinsic :: iso_c_binding, only : c_size_t
    use MemoryAllocator_mod, only : Allocator
    use Datastore_mod, only : theDatastore
 
!   Passed variables
 
    class(Material),   intent(inout) :: self
    logical(kind=1),   intent(in)    :: nonLTE
    logical(kind=1),   intent(in)    :: fromRestart

!   Some of these data structures are not mapped to the GPU.  Do not use pinned memory.
!   We still use the MemoryAllocator here, as it adds entries in the conduit
!   datastore for each allocation.  Adding these to the conduit datastore makes
!   them easily accessible from C++, or dumping to files for troubleshooting.

    self%label = "material"

    call Allocator%allocate(Size%usePinnedMemory, self%label, "cve", self%cve, Size%nzones)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "rho", self%rho, Size%nzones)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "nez", self%nez, Size%nzones)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "SMatEff", self%SMatEff, Size%nzones)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "stimComptonMult", self%stimComptonMult, Size%nzones)

    call Allocator%allocate(Size%usePinnedMemory, self%label, "nonLinearIterations", self%nonLinearIterations, Size%ncornr)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "PowerEmitted", self%PowerEmitted, Size%ncornr)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "PowerCompton", self%PowerCompton, Size%ncornr)

    call Allocator%allocate(Size%usePinnedMemory, self%label, "denec", self%denec, Size%ncornr)
    if (.not. fromRestart) then
       call Allocator%allocate(Size%usePinnedMemory, self%label, "Tec", self%Tec, Size%ncornr)
    endif
    call Allocator%allocate(Size%usePinnedMemory, self%label, "Tecn", self%Tecn, Size%ncornr)

    call Allocator%allocate(Size%usePinnedMemory, self%label, "sigA", self%sigA, Size%ngr, Size%nzones)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "sigS", self%sigS, Size%ngr, Size%nzones)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "EmissionRate", self%EmissionRate, Size%ngr, Size%ncornr)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "Eta", self%Eta, Size%ncornr)

    call Allocator%allocate(.FALSE., self%label, "denez", self%denez, Size%nzones)
    call Allocator%allocate(.FALSE., self%label, "Trz", self%Trz, Size%nzones)
    call Allocator%allocate(.FALSE., self%label, "Tez", self%Tez, Size%nzones)
    call Allocator%allocate(.FALSE., self%label, "Trzn", self%Trzn, Size%nzones)
    call Allocator%allocate(.FALSE., self%label, "Tezn", self%Tezn, Size%nzones)
    call Allocator%allocate(.FALSE., self%label, "Tezold", self%Tezold, Size%nzones)
    call Allocator%allocate(.FALSE., self%label, "EnergyDensityOld", self%EnergyDensityOld, Size%nzones)

    allocate( self% isVoid(Size% nzones) )
    self% isVoid(:)              = .FALSE.
    allocate( self% nonLTE(Size% nzones) )
    self% nonLTE(:)              = .FALSE.

    if ( nonLTE ) then
      self% nonLTE(:) = .TRUE.
      call Allocator%allocate(.FALSE., self%label, "emis", self%emis, Size%ngr, Size%nzones)
      call Allocator%allocate(.FALSE., self%label, "demisdT", self%demisdT, Size%ngr, Size%nzones)
    endif

    return
 
  end subroutine Material_ctor

!=======================================================================
! destruct interface
!=======================================================================
                                                            
  subroutine Material_dtor(self, nonLTE)

    use Size_mod, only : Size
    use MemoryAllocator_mod, only : Allocator
    use Datastore_mod, only : theDatastore
                 
!   Passed variables
                    
    class(Material),   intent(inout) :: self
    logical(kind=1),   intent(in)    :: nonLTE

    call Allocator%deallocate(Size%usePinnedMemory, self%label, "cve", self%cve)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "rho", self%rho)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "nez", self%nez)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "SMatEff", self%SMatEff)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "stimComptonMult", self%stimComptonMult)

    call Allocator%deallocate(Size%usePinnedMemory, self%label, "nonLinearIterations", self%nonLinearIterations)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "PowerEmitted", self%PowerEmitted)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "PowerCompton", self%PowerCompton)

    call Allocator%deallocate(Size%usePinnedMemory, self%label, "denec", self%denec)

    call Allocator%deallocate(Size%usePinnedMemory, self%label, "Tec", self%Tec)

    call Allocator%deallocate(Size%usePinnedMemory, self%label, "Tecn", self%Tecn)

    call Allocator%deallocate(Size%usePinnedMemory, self%label, "sigA", self%sigA)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "sigS", self%sigS)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "EmissionRate", self%EmissionRate)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "Eta", self%Eta)

    call Allocator%deallocate(.FALSE., self%label, "denez", self%denez)
    call Allocator%deallocate(.FALSE., self%label, "Trz", self%Trz)
    call Allocator%deallocate(.FALSE., self%label, "Tez", self%Tez)
    call Allocator%deallocate(.FALSE., self%label, "Trzn", self%Trzn)
    call Allocator%deallocate(.FALSE., self%label, "Tezn", self%Tezn)
    call Allocator%deallocate(.FALSE., self%label, "Tezold", self%Tezold)
    call Allocator%deallocate(.FALSE., self%label, "EnergyDensityOld", self%EnergyDensityOld)

    deallocate( self% isVoid )
    deallocate( self% nonLTE )

    if ( nonLTE ) then
      call Allocator%deallocate(.FALSE., self%label, "emis", self%emis)
      call Allocator%deallocate(.FALSE., self%label, "demisdT", self%demisdT)
    endif

    return
          
  end subroutine Material_dtor

end module Material_mod
