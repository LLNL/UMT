!***********************************************************************
!                        Version 1:  02/06, PFN                        *
!                                                                      *
!   getVolume   - This routine gets mesh volumes and areas.            *
!                                                                      *
!***********************************************************************

   subroutine getVolume() BIND(C,NAME="teton_getvolume")

!  Include

   use kind_mod
   use Size_mod
   use Geometry_mod

   implicit none


!  Save the "old" corner volumes

   Geom% VolumeOld(:) = Geom% Volume(:)

!  Call the appropriate function based on the spatial dimensionality

   if (Size% ndim == 3) then

     call volumeUCBxyz 

   else if (Size% ndim == 2) then

     if (Size% usePWLD) then
       call volumePWLDrz
     else
       call volumeUCBrz
     endif

   else if (Size% ndim == 1) then

     call volumeSCB1D 

   endif 



   return
   end subroutine getVolume 


