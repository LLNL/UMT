!***********************************************************************
!                        Version 1:  02/94, PFN                        *
!                                                                      *
!   PROFINT - Performs an interpolation of all profiles at the current *
!             time.                                                    *
!                                                                      *
!***********************************************************************

   subroutine profint

   use kind_mod
   use constant_mod
   use radconstant_mod
   use TimeStepControls_mod
   use BoundaryList_mod
   use Profile_mod
   use QuadratureList_mod
   use Size_mod
!  From PhysicsUtils library
   use emissionutils

   implicit none

!  Local

   type(Profile),  pointer  :: ProfID
                                                                                         
   integer    :: prof,ig,n,nl,nh,nlow,nhigh,ngr,NumTimes,nSource

   real(adqt) :: ac,timel,timeh,dtime,timerad,Mult,kb
   real(adqt) :: Te

   real(adqt) :: gnu(Size%ngr+1)
   real(adqt) :: planck(Size%ngr+1)

   logical(kind=1)  :: BlackBody

!  Constants

   ac      = rad_constant*speed_light
   kb      = one
   ngr     = Size% ngr

   timerad = getRadTime(DtControls)

   gnu(:) = getEnergyGroups(Quad,ngr)

!  We need to check each profile to see if it is
!  on (set to zero if off).

   nSource = getNumberOfSource(RadBoundary) 

   ProfileLoop: do prof=1,nSource

     ProfID => getProfile(RadBoundary, prof)

     BlackBody         = ProfID% BlackBody
     NumTimes          = ProfID% NumTimes
     Mult              = ProfID% Multiplier
     ProfID% profileOn = .FALSE. 
     timel             = zero
     nl                = 1

!    Interpolate in time

     ProfID% InterpValues(1:ngr) = zero

     TestTime: if ( (timerad >= ProfID% Times(1)) .and.  &
                    (timerad <= ProfID% Times(NumTimes)) ) then

       ProfID% profileOn = .TRUE.

!      Find time factor

       n = 0

       TimeIteration: do 
         n = n + 1

         if (timerad >= ProfID% Times(n) .and. n < NumTimes) then
           timel = ProfID% Times(n)
           nl    = n
           cycle TimeIteration
         else
           exit TimeIteration 
         endif

       enddo TimeIteration

       nh    = nl + 1
       timeh = ProfID% Times(nh)
       dtime = (timerad - timel)/(timeh - timel)

!********************
!  Temperature      *
!********************
                                                                                         
       TestType: if ( BlackBody ) then 
                                                                                         
!  Find the interpolated temperature and then generate a
!  Planckian energy spectrum at that temperature

         Te  =  ProfID% Values(nl) + dtime*  &
               (ProfID% Values(nh) - ProfID% Values(nl))

         if (Te > zero) then

!  Compute the fraction of the total emission in each energy group
!  The input for RTPLNK is (h*nu)/(k*Te).

           if (ngr == 1) then
  
             ProfID% InterpValues(1) = Mult*ac*Te*Te*Te*Te

           else

             call integrateBlackBodyGroups(Te,kb,ac,ngr,gnu,planck)

             do ig=1,ngr
               ProfID% InterpValues(ig) = Mult*planck(ig)
             enddo

           endif

         else

           ProfID% InterpValues(:) = zero
 
         endif

!**********************
! Frequency-Dependent * 
!**********************

       else

!  For frequency-dependent sources there are NGR values for
!  each time. The code is expecting the user to input values with
!  units of Energy/(Volume*PhotonEnergy)

         nlow  = (nl - 1)*ngr
         nhigh = (nh - 1)*ngr

         do ig=1,ngr
           ProfID% InterpValues(ig) = speed_light*Mult*         &
                                     ( gnu(ig+1) - gnu(ig) )*   &
                     (ProfID% Values(nlow+ig) + dtime*          &
                     (ProfID% Values(nhigh+ig) - ProfID% Values(nlow+ig)))
         enddo


       endif TestType

     endif TestTime

   enddo ProfileLoop


   return
   end subroutine profint

