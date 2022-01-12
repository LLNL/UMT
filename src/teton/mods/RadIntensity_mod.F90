! Radiation Intensity Module:  angle-dependent and scalar intensity
                                                                                 
module RadIntensity_mod
                                                                                 
  use kind_mod

  private

! public interfaces
                                                                                             
  public construct
  public destruct

  ! TODO: These are stored per set.  But that really isn't needed.
  type, public :: RadIntensity

     real(adqt), pointer, contiguous  :: radEnergy(:)
     real(adqt), pointer, contiguous  :: RadiationForce(:,:)
     real(adqt), pointer, contiguous  :: RadiationFlux(:,:,:)
     real(adqt), pointer, contiguous  :: EddingtonTensorDiag(:,:)

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

  subroutine RadIntensity_ctor(self, Groups)
 
    use Size_mod
    use constant_mod

    implicit none
 
!   Passed variables
 
    type(RadIntensity), intent(inout) :: self
    integer,            intent(in)    :: Groups


    allocate( self% radEnergy(Size%nzones) )
    allocate( self% RadiationForce(Size%ndim,Size%ncornr) )
    allocate( self% RadiationFlux(Size%ndim,Groups,Size%nzones) )
    allocate( self% EddingtonTensorDiag(Size%ndim,Size%nzones) )

    self% radEnergy(:)              = zero
    self% RadiationForce(:,:)       = zero
    self% RadiationFlux(:,:,:)      = zero
    self% EddingtonTensorDiag(:,:)  = zero

    return
 
  end subroutine RadIntensity_ctor



!=======================================================================
! destruct interface
!=======================================================================
                                                            
  subroutine RadIntensity_dtor(self)
                                      
    implicit none
                 
!   Passed variables
                    
    type(RadIntensity),    intent(inout) :: self


    deallocate( self% radEnergy )
    deallocate( self% RadiationForce )
    deallocate( self% RadiationFlux )
    deallocate( self% EddingtonTensorDiag )


    return
          
  end subroutine RadIntensity_dtor

end module RadIntensity_mod

