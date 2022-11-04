!***********************************************************************
!                        Last Update:  03/2012, PFN                    *
!                                                                      *
!    SetNodeVelocity -  Called from host to set node velocities        *
!                       in the Geometry module.                        *
!                                                                      *
!***********************************************************************

   subroutine setNodeVelocity(zoneID, velocity) &
        BIND(C,NAME="teton_setnodevelocity")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use Geometry_mod
   use radconstant_mod

   implicit none

!  Arguments

   integer(C_INT), intent(in) :: zoneID 
   real(C_DOUBLE), intent(in) :: velocity(Size% ndim,Size% maxCorner)

!  Local

   integer  :: c, c0, nCorner


   nCorner = Geom% numCorner(zoneID)
   c0      = Geom% cOffSet(zoneID) 

   do c=1,nCorner
     Geom% VoC(:,c0+c) = velocity(:,c)/speed_light
   enddo


   return
   end subroutine setNodeVelocity 
