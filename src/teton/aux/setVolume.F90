!***********************************************************************
!                        Version 1:  01/25, BCY                        *
!                                                                      *
!   setVolume   - This routine updates mesh volumes and areas.         *
!                 Unlike the old deprecated "getVolume", it does NOT   *
!                 update Geom%VolumeOld so it is safe to call multiple *
!                 times.                                               *
!                                                                      *
!***********************************************************************

   subroutine setVolume() BIND(C,NAME="teton_setvolume")

!  Include

   use Size_mod
   use Geometry_mod

   implicit none

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
   end subroutine setVolume 


