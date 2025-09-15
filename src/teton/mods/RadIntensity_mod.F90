! Radiation Intensity Module:  angle-dependent and scalar intensity
                                                                                 
module RadIntensity_mod
                                                                                 
  use kind_mod

  private

! public interfaces
                                                                                             
  public construct
  public destruct

  type, public :: RadIntensity

     real(adqt), pointer, contiguous  :: radEnergy(:)             => null()
     real(adqt), pointer, contiguous  :: RadiationForce(:,:)      => null()
     real(adqt), pointer, contiguous  :: RadiationFlux(:,:,:)     => null()
     real(adqt), pointer, contiguous  :: EddingtonTensorDiag(:,:) => null()
     real(adqt), pointer, contiguous  :: PhiTotal(:,:)            => null()
     real(adqt), pointer, contiguous  :: RadEnergyDensity(:,:)    => null()

!    Misc
     character(len=12) :: label ! A string descriptor for this module.

  end type RadIntensity 

  type(RadIntensity), pointer, public :: Rad => null() 

  interface construct
    module procedure RadIntensity_ctor
  end interface

 
  interface destruct
    module procedure RadIntensity_dtor
  end interface
 
contains
 
!=======================================================================
! construct interface
!=======================================================================

  subroutine RadIntensity_ctor(self)
 
    use Size_mod
    use constant_mod
    use MemoryAllocator_mod
    use Datastore_mod, only : theDatastore
    use conduit_obj, only : node
    use, intrinsic :: iso_c_binding, only : c_size_t

    implicit none
 
!   Passed variables
 
    type(RadIntensity), intent(inout) :: self

!   Local

    integer(kind=c_size_t) :: num_elements_size_t
    type(node) :: blueprint_node
    character(len=60) :: field_path

    self%label = "radintensity"


!   The following are used on the GPU and require pinned memory

    call Allocator%allocate(Size%usePinnedMemory,self%label,"PhiTotal", self% PhiTotal, Size% ngr, Size% ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"radEnergy", self% radEnergy, Size%nZones)

!   The following are only used on the CPU

    call Allocator%allocate(.FALSE.,self%label, "RadiationForce", self% RadiationForce, Size%ndim, Size%ncornr)
    call Allocator%allocate(.FALSE.,self%label, "RadiationFlux", self% RadiationFlux, Size%ndim, Size%ngr, Size%nzones)
    call Allocator%allocate(.FALSE.,self%label, "EddingtonTensorDiag", self% EddingtonTensorDiag, Size%ndim, Size%nzones)
    call Allocator%allocate(.FALSE.,self%label, "RadEnergyDensity", self% RadEnergyDensity, Size%nzones, Size%ngr)

    ! Describe the field in the mesh blueprint.  This is a field over the main
    ! mesh (zone mesh), over zones (elements) as opposed to vertices (nodes).
    ! Multi-value fields are not support for visualization, however.
    num_elements_size_t = Size% nzones * Size% ngr

    blueprint_node = theDatastore%get_blueprint_node()
    field_path = "fields/radiation_energy_density/"
    call blueprint_node%set_path_external_float64_ptr(trim(field_path)//"values", self% RadEnergyDensity,num_elements_size_t )
    call blueprint_node%set_path(trim(field_path)//"association", "element")
    call blueprint_node%set_path(trim(field_path)//"topology", "main")

    num_elements_size_t = Size% nzones

    field_path = "fields/rad_energy/"
    call blueprint_node%set_path_external_float64_ptr(trim(field_path)//"values", self% radEnergy,num_elements_size_t )
    call blueprint_node%set_path(trim(field_path)//"association", "element")
    call blueprint_node%set_path(trim(field_path)//"topology", "main")

!   Initialize

    self% radEnergy(:)              = zero
    self% PhiTotal(:,:)             = zero
    self% RadiationForce(:,:)       = zero
    self% RadiationFlux(:,:,:)      = zero
    self% EddingtonTensorDiag(:,:)  = zero
    self% RadEnergyDensity(:,:)     = zero

    return
 
  end subroutine RadIntensity_ctor


!=======================================================================
! destruct interface
!=======================================================================
                                                            
  subroutine RadIntensity_dtor(self)

    use Size_mod
    use MemoryAllocator_mod
    use Datastore_mod, only : theDatastore
                                      
    implicit none
                 
!   Passed variables
                    
    type(RadIntensity),    intent(inout) :: self

    call Allocator%deallocate(Size%usePinnedMemory,self%label, "PhiTotal", self% PhiTotal)
    call Allocator%deallocate(Size%usePinnedMemory,self%label, "radEnergy", self% radEnergy)
    call Allocator%deallocate(.FALSE.,             self%label, "RadiationForce", self% RadiationForce)
    call Allocator%deallocate(.FALSE.,             self%label, "RadiationFlux", self% RadiationFlux)
    call Allocator%deallocate(.FALSE.,             self%label, "EddingtonTensorDiag", self% EddingtonTensorDiag)

    if (.NOT. theDatastore%partitioning()) then
      call theDatastore%root%remove_path("blueprint/fields/radiation_energy_density")
    endif
    call Allocator%deallocate(.FALSE.,self%label,"RadEnergyDensity", self%RadEnergyDensity)


    return
          
  end subroutine RadIntensity_dtor

end module RadIntensity_mod

