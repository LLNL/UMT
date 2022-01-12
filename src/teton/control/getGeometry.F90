!***********************************************************************
!                        Version 1:  02/06, PFN                        *
!                                                                      *
!   getGeometry - This routine gets mesh volumes and areas.            *
!                                                                      *
!***********************************************************************

   subroutine getGeometry

!  Include

   use kind_mod
   use Size_mod

   implicit none

!  Call the appropriate function based on the spatial dimensionality

   if (Size% ndim == 3) then

     call geometryUCBxyz 

   else if (Size% ndim == 2) then

     if ( Size% usePWLD ) then
       call geometryPWLDrz
     else
       call geometryUCBrz
     endif

   else if (Size% ndim == 1) then

     call geometrySCB1D 

   endif 



   return
   end subroutine getGeometry 


