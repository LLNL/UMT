!***********************************************************************
!                         Version 0: 04/06 PFN                         *
!                                                                      *
!    SetMaterialSource -  Called from host to update the               *
!                         energy deposition rate to the material.      *
!                                                                      *
!***********************************************************************

   subroutine setMaterialSource(SMatEff) &
        BIND(C,NAME="teton_setmaterialsource")

   use ISO_C_BINDING
   use kind_mod
   use Size_mod 
   use Material_mod

   implicit none

!  Arguments

   real(C_DOUBLE), intent(in) :: SMatEff(Size%nzones)

!  Local

   integer  :: zone

   do zone=1,Size%nzones
     Mat% SMatEff(zone) = SMatEff(zone)
   enddo


   return
   end subroutine setMaterialSource 
