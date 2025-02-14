!***********************************************************************
!                        Version 1:  08/2024, BCY                      *
!                                                                      *
!   getRadiationEnergyDensityPointer                                   *
!                                                                      *
!***********************************************************************
 
   subroutine getRadiationEnergyDensityPointer(RadEnergyDensityPtr) &
        BIND(C,NAME="teton_getradiationenergydensityptr")

   USE ISO_C_BINDING
   use kind_mod
   use RadIntensity_mod

   implicit none 

!  Arguments

   type(C_PTR), intent(out)  :: RadEnergyDensityPtr

!  Local

!***********************************************************************
!  Update Radiation Energy Density                                     * 
!***********************************************************************
 
   RadEnergyDensityPtr = C_LOC(Rad% RadEnergyDensity)

   return
   end subroutine getRadiationEnergyDensityPointer
