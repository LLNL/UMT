!***********************************************************************
!                         Last Update: 02/2018 PFN                     *
!                                                                      *
!    SetOpacity -  Called from host to update opacities                *
!                  after they have been changed.                       *
!                                                                      *
!***********************************************************************

   subroutine setOpacity(zone, siga, sigs, useTableSigmaS) &
        BIND(C,NAME="teton_setopacity")

   use ISO_C_BINDING
   use kind_mod
   use Size_mod 
   use Material_mod

   use radconstant_mod

   implicit none

!  Arguments

   integer(C_INT),  intent(in) :: zone
   real(C_DOUBLE),  intent(in) :: siga(Size% ngr) 
   real(C_DOUBLE),  intent(in) :: sigs(Size% ngr)
   logical(C_BOOL), intent(in) :: useTableSigmaS 
   integer                     :: grp

!  Local

!  TODO: Make sigmaCeiling from radconstant_mod tunable by the host codes? See TETON-131.
   do grp=1,Size% ngr
     Mat% siga(grp,zone) = Mat% siga(grp,zone) + min(siga(grp),sigmaCeiling)

     if (useTableSigmaS) then
       Mat% sigs(grp,zone) = Mat% sigs(grp,zone) + min(sigs(grp),sigmaCeiling)
     end if
   enddo

   return
   end subroutine setOpacity
