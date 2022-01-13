!***********************************************************************
!                         Last Update: 02/2012 PFN                     *
!                                                                      *
!    ResetProfile -  Called to alter an existing source profile        *
!                    in the profile list.                              *
!                                                                      *
!***********************************************************************
   subroutine resetProfile(TetonProfileID,               &
                         NumTimes, NumValues,            &
                         Multiplier, Times, Values)      &
                         BIND(C,NAME="teton_resetprofile_internal")

   use ISO_C_BINDING
   use kind_mod
   use BoundaryList_mod

   implicit none


!  Arguments

   integer(C_INT), intent(in)       :: TetonProfileID
   integer(C_INT), intent(in)       :: NumTimes
   integer(C_INT), intent(in)       :: NumValues
   real(C_DOUBLE), intent(in)       :: Multiplier
   real(C_DOUBLE), intent(in)       :: Times(NumTimes)
   real(C_DOUBLE), intent(in)       :: Values(NumValues)

!  Add this profile to the list 

   call resetSourceProfile(RadBoundary,       &
                           TetonProfileID,    &
                           NumTimes,          &
                           NumValues,         &
                           Multiplier,        &
                           Times,             &
                           Values)


   return
   end subroutine resetProfile 



