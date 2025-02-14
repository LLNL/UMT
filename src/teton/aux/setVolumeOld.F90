!***********************************************************************
!                        Version 1:  01/25, BCY                        *
!                                                                      *
!   setVolumeOld   - This routine sets the VolumeOld variable equal    *
!                    to the current Volume variable.                   *
!                                                                      *
!***********************************************************************

   subroutine setVolumeOld() BIND(C,NAME="teton_setvolumeold")

!  Include

   use Geometry_mod

   implicit none

!  Update the "old" corner volumes so that dV = Volume-VolumeOld = 0

   Geom% VolumeOld(:) = Geom% Volume(:)

   return
   end subroutine setVolumeOld 


