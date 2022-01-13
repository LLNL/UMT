!***********************************************************************
!                         Last Update: 02/2012 PFN                     *
!                                                                      *
!    AddProfile   -  Called from host to add a source profile          *
!                    to the profile list.                              *
!                                                                      *
!    Output:  TetonProfileID - integer index correpsonding to profile  *
!                                                                      *
!***********************************************************************
   subroutine addProfile(NumTimes, NumValues,            &
                         Multiplier,                     &
                         BlackBody, Isotropic,           &
                         Times, Values, TetonProfileID)  &
                         BIND(C,NAME="teton_addprofile_internal")

   use ISO_C_BINDING
   use kind_mod
   use BoundaryList_mod

   implicit none


!  Arguments

   integer(C_INT), intent(in)       :: NumTimes
   integer(C_INT), intent(in)       :: NumValues
              
   real(C_DOUBLE), intent(in)       :: Multiplier

   logical(C_BOOL), intent(in)      :: BlackBody
   logical(C_BOOL), intent(in)      :: Isotropic
              
   real(C_DOUBLE), intent(in)       :: Times(NumTimes)
   real(C_DOUBLE), intent(in)       :: Values(NumValues)

   integer(C_INT), intent(inout)    :: TetonProfileID

!  Add this profile to the list 

   call setProfile(RadBoundary,       &
                   NumTimes,          &
                   NumValues,         &
                   Multiplier,        &
                   BlackBody,         &
                   Isotropic,         &
                   Times,             &
                   Values,            &
                   TetonProfileID)


   return
   end subroutine addProfile 



