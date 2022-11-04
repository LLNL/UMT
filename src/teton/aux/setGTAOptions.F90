!***********************************************************************
!                         Last Update: 03/2022 BCY                     *
!                                                                      *
!    setGTAOptions  -  Called from host to update two GTA options      *
!                                                                      *
!    Input:   enforceHardGTAIterMax:                                   *
!                 Force the GTA iteration to exit,                     *
!                 even if it hasn't converged to the minimum tolerance *
!                 True by default for new GTA.  False for old GTA.     *
!             forceExtraOuter:                                         *
!                 Force another outer iteration if                     *
!                 GTA isn't converged.  False by default.              *
!                                                                      *
!***********************************************************************
   subroutine setGTAOptions(enforceHardGTAIterMax, forceExtraOuter) &
     BIND(C,NAME="teton_setgtaoptions")

   use GreyAcceleration_mod
   use Size_mod
   use ISO_C_BINDING 
   
   implicit none


!  Arguments

   logical(C_BOOL), intent(in) :: enforceHardGTAIterMax
   logical(C_BOOL), intent(in) :: forceExtraOuter

!  Update controls

   if (associated(GTA)) then
      GTA% enforceHardGTAIterMax = enforceHardGTAIterMax
      GTA% forceExtraOuter       = forceExtraOuter
   else if (Size%myRankInGroup == 0) then
      print *, "You're trying to set GTA options before GTA is constructed! Teton will ignore this API call!"
   endif


   return
   end subroutine setGTAOptions 



