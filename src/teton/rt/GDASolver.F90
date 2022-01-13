!***********************************************************************
!                       Last Update:  03/2012, PFN                     *
!                                                                      *
!   GDASolver   - Performs a grey-diffusion acceleration for 1D        *
!                 problems.                                            *
!                                                                      *
!***********************************************************************

   subroutine GDASolver 

   use kind_mod
   use Size_mod

   implicit none

!  Local


!  Constants


!***********************************************************************
!  Grey Diffusion Synthetic Acceleration                               *
!***********************************************************************


   call rtgdac
   call pentfbsb


   return
   end subroutine GDASolver 
