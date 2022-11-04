!***********************************************************************
!                        Version 1:  10/2016, PFN                      *
!                                                                      *
!   getRadiationEnergyDensity                                          *
!                                                                      *
!***********************************************************************
 
   subroutine getRadiationEnergyDensity(RadEnergyDensity) &
        BIND(C,NAME="teton_getradiationenergydensity")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use constant_mod
   use RadIntensity_mod

   implicit none 

!  Arguments

   real(C_DOUBLE), intent(inout)  :: RadEnergyDensity(Size% nzones,Size%ngr) 

!  Local

!***********************************************************************
!  Update Radiation Energy Density                                     * 
!***********************************************************************
 
   RadEnergyDensity(:,:) = Rad% RadEnergyDensity(:,:) 


   return
   end subroutine getRadiationEnergyDensity 
