! Material Module:  material properties 
                                                                                 
module Material_mod
  use kind_mod, only : adqt
  implicit none

  private

  type, public :: Material

     real(adqt),      pointer, contiguous :: cve(:)             => null()  ! "effective" specific heat 
     real(adqt),      pointer, contiguous :: rho(:)             => null()  ! Material Density 
     real(adqt),      pointer, contiguous :: nez(:)             => null()  ! Electron Density 
     real(adqt),      pointer, contiguous :: SMatEff(:)         => null()  ! Source to the material 
     real(adqt),      pointer, contiguous :: denez(:)           => null()  ! denez
     real(adqt),      pointer, contiguous :: Trz(:)             => null()  ! trz
     real(adqt),      pointer, contiguous :: Tez(:)             => null()  ! tez
     real(adqt),      pointer, contiguous :: Trzn(:)            => null()  ! trzn
     real(adqt),      pointer, contiguous :: Tezn(:)            => null()  ! tezn
     real(adqt),      pointer, contiguous :: Tezold(:)          => null()  ! tezold
     real(adqt),      pointer, contiguous :: stimComptonMult(:) => null()  ! stimulated Compton mult
     logical(kind=1), pointer, contiguous :: isVoid(:)          => null()  ! void zone
     logical(kind=1), pointer, contiguous :: nonLTE(:)          => null()  ! non-local thermodynamic equilibrium

!    Thread-Safe Tallies
     integer,    pointer, contiguous :: nonLinearIterations(:)  => null()
     real(adqt), pointer, contiguous :: PowerEmitted(:)         => null()
     real(adqt), pointer, contiguous :: PowerCompton(:)         => null()

     real(adqt), pointer, contiguous :: denec(:)                => null()  ! denec(nCorner)
     real(adqt), pointer, contiguous :: Tec(:)                  => null()  ! tec(nCorner)
     real(adqt), pointer, contiguous :: Tecn(:)                 => null()  ! tecn(nCorner)

     real(adqt), pointer, contiguous :: sigA(:,:)               => null()  ! absorption coefficient 
     real(adqt), pointer, contiguous :: sigS(:,:)               => null()  ! scattering coefficient 
     real(adqt), pointer, contiguous :: Eta(:)                  => null()  ! Eta(nCorner)

!    For nonLTE
     real(adqt), pointer, contiguous :: emis(:,:)               => null()  ! emissivity
     real(adqt), pointer, contiguous :: demisdT(:,:)            => null()

!    Misc
     character(len=8) :: label

     contains

       procedure, public :: construct => Material_ctor
       procedure, public :: destruct => Material_dtor

  end type Material

  type(Material), pointer, public :: Mat                        => null()

contains
 
!=======================================================================
! construct interface
!=======================================================================

  subroutine Material_ctor(self, nonLTE)
    use Size_mod, only : Size
    use constant_mod, only : zero
    use MemoryAllocator_mod
 
!   Passed variables
 
    class(Material),    intent(inout) :: self
    logical(kind=1),   intent(in)    :: nonLTE

    logical :: usePinnedMemory
    usePinnedMemory = .FALSE.

    self%label = "material"

    call Allocator%allocate(usePinnedMemory, self%label, "cve", self%cve, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "rho", self%rho, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "nez", self%nez, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "SMatEff", self%SMatEff, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "denez", self%denez, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "Trz", self%Trz, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "Tez", self%Tez, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "Trzn", self%Trzn, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "Tezn", self%Tezn, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "Tezold", self%Tezold, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "stimComptonMult", self%stimComptonMult, Size%nzones)

    allocate( self% isVoid(Size% nzones) )
    self% isVoid(:)              = .FALSE.
    allocate( self% nonLTE(Size% nzones) )
    self% nonLTE(:)              = .FALSE.

    call Allocator%allocate(usePinnedMemory, self%label, "nonLinearIterations", self%nonLinearIterations, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "PowerEmitted", self%PowerEmitted, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "PowerCompton", self%PowerCompton, Size%nzones)

    call Allocator%allocate(usePinnedMemory, self%label, "denec", self%denec, Size%ncornr)
    call Allocator%allocate(usePinnedMemory, self%label, "Tec", self%Tec, Size%ncornr)
    call Allocator%allocate(usePinnedMemory, self%label, "Tecn", self%Tecn, Size%ncornr)

    call Allocator%allocate(usePinnedMemory, self%label, "sigA", self%sigA, Size%ngr, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "sigS", self%sigS, Size%ngr, Size%nzones)
    call Allocator%allocate(usePinnedMemory, self%label, "Eta", self%Eta, Size%ncornr)

    if ( nonLTE ) then
      self% nonLTE(:) = .TRUE.
      call Allocator%allocate(usePinnedMemory, self%label, "emis", self%emis, Size%ngr, Size%nzones)
      call Allocator%allocate(usePinnedMemory, self%label, "demisdT", self%demisdT, Size%ngr, Size%nzones)
    endif

    return
 
  end subroutine Material_ctor

!=======================================================================
! destruct interface
!=======================================================================
                                                            
  subroutine Material_dtor(self, nonLTE)
    use MemoryAllocator_mod
                 
!   Passed variables
                    
    class(Material),    intent(inout) :: self
    logical(kind=1),   intent(in)    :: nonLTE

    logical :: usePinnedMemory
    usePinnedMemory = .FALSE.

    call Allocator%deallocate(usePinnedMemory, self%label, "cve", self%cve)
    call Allocator%deallocate(usePinnedMemory, self%label, "rho", self%rho)
    call Allocator%deallocate(usePinnedMemory, self%label, "nez", self%nez)
    call Allocator%deallocate(usePinnedMemory, self%label, "SMatEff", self%SMatEff)
    call Allocator%deallocate(usePinnedMemory, self%label, "denez", self%denez)
    call Allocator%deallocate(usePinnedMemory, self%label, "Trz", self%Trz)
    call Allocator%deallocate(usePinnedMemory, self%label, "Tez", self%Tez)
    call Allocator%deallocate(usePinnedMemory, self%label, "Trzn", self%Trzn)
    call Allocator%deallocate(usePinnedMemory, self%label, "Tezn", self%Tezn)
    call Allocator%deallocate(usePinnedMemory, self%label, "Tezold", self%Tezold)
    call Allocator%deallocate(usePinnedMemory, self%label, "stimComptonMult", self%stimComptonMult)

    deallocate( self% isVoid )
    deallocate( self% nonLTE )

    call Allocator%deallocate(usePinnedMemory, self%label, "nonLinearIterations", self%nonLinearIterations)
    call Allocator%deallocate(usePinnedMemory, self%label, "PowerEmitted", self%PowerEmitted)
    call Allocator%deallocate(usePinnedMemory, self%label, "PowerCompton", self%PowerCompton)

    call Allocator%deallocate(usePinnedMemory, self%label, "denec", self%denec)
    call Allocator%deallocate(usePinnedMemory, self%label, "Tec", self%Tec)
    call Allocator%deallocate(usePinnedMemory, self%label, "Tecn", self%Tecn)

    call Allocator%deallocate(usePinnedMemory, self%label, "sigA", self%sigA)
    call Allocator%deallocate(usePinnedMemory, self%label, "sigS", self%sigS)
    call Allocator%deallocate(usePinnedMemory, self%label, "Eta", self%Eta)

    if ( nonLTE ) then
      call Allocator%deallocate(usePinnedMemory, self%label, "emis", self%emis)
      call Allocator%deallocate(usePinnedMemory, self%label, "demisdT", self%demisdT)
    endif

    return
          
  end subroutine Material_dtor

end module Material_mod
