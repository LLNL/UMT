!***********************************************************************
!                        Version 1:  03/02, PFN                        *
!                                                                      *
!   setEnergyDensity - Calculate the radiation energy density for      *
!                      each group and each zone.  These quantities     *
!                      are remapped and used to construct new corner   *
!                      intensties after remap.                         *
!                                                                      *
!   Output:      RadEnergyDensity  -  zonal rad energy density (E/V)   *         
!                                                                      *
!***********************************************************************
 
   subroutine setEnergyDensity 

   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod

   implicit none 

!  Local

   integer    :: c 
   integer    :: c0
   integer    :: g
   integer    :: nCorner
   integer    :: zone 

   real(adqt) :: factor 
   real(adqt) :: VolFrac

!***********************************************************************
!  Update Radiation Energy Density                                     * 
!***********************************************************************
 
   factor = one/speed_light

   Rad% RadEnergyDensity(:,:) = zero

   ZoneLoop: do zone=1,Size%nzones

     nCorner = Geom% numCorner(zone)
     c0      = Geom% cOffSet(zone)

     do c=1,nCorner

       VolFrac = factor*Geom% Volume(c0+c)/Geom% VolumeZone(zone)
       
       do g=1,Size% ngr
         Rad% RadEnergyDensity(zone,g) = Rad% RadEnergyDensity(zone,g) +  &
                                         VolFrac*Rad% PhiTotal(g,c0+c)
       enddo

     enddo

   enddo ZoneLoop


   return
   end subroutine setEnergyDensity 
