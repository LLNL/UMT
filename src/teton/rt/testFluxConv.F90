!***********************************************************************
!                        Version 1:  03/2009, PFN                      *
!                                                                      *
!   testFluxConv - Monitors convergence of the incident current on     *
!                  shared boundaries.                                  *
!                                                                      *
!***********************************************************************
   subroutine testFluxConv(cSetID, FluxConverged) 

   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use CommSet_mod
   use iter_control_list_mod
   use iter_control_mod
   use radconstant_mod

   implicit none

!  Arguments

   integer,          intent(in)    :: cSetID
   logical (kind=1), intent(inout) :: FluxConverged

!  Local

   type(IterControl) , pointer :: incidentFluxControl => NULL()
   type(IterControl) , pointer :: temperatureControl  => NULL()
   type(CommSet),      pointer :: CSet

   integer    :: bin
   integer    :: nConv
   integer    :: NangBin

   real(adqt) :: fluxTolerance
   real(adqt) :: tempError
   real(adqt) :: minTolerance
   real(adqt) :: totalFlux
   real(adqt) :: weight
   real(adqt) :: threshold

!  Constants

   parameter (minTolerance=0.01d0)
   parameter (threshold=0.001d0)

   CSet                => getCommSetData(Quad, cSetID)
   incidentFluxControl => getIterationControl(IterControls,"incidentFlux")
   temperatureControl  => getIterationControl(IterControls,"temperature")

!  Use the error in the nonlinear iteration as a guide for
!  converging the incident flux iteration

   tempError     = getGlobalError(temperatureControl)
   fluxTolerance = getEpsilonPoint(incidentFluxControl)

   fluxTolerance = max(fluxTolerance, tempError/twenty)
   fluxTolerance = min(fluxTolerance, minTolerance)

!  Tally the number of transport sweeps performed for this group set

   do bin=1,CSet% NumBin0
     if ( .not. CSet% Converged(bin) ) then
       NangBin          = CSet% NangBinList(bin) 
       CSet% fluxSweeps = CSet% fluxSweeps + NangBin 
     endif
   enddo

!  Compute errors

   totalFlux = CSet% totalIncFlux

   if (abs(totalFlux) < adqtEpsilon*speed_light*rad_constant*Size%tr4floor) then

     CSet% relError(:) = zero

   else

     do bin=1,CSet% NumBin0
       weight = CSet% IncFlux(bin)/totalFlux

       if (weight > threshold) then
         CSet% relError(bin) = abs( (CSet% IncFlux(bin)  -      &
                                     CSet% IncFluxOld(bin)) )/  &
                                     CSet% IncFlux(bin)
       else
         CSet% relError(bin) = zero
       endif
     enddo

   endif

!  Find how many bins have converged 

   nConv = 0
   do bin=1,CSet% NumBin0
     if (CSet% relError(bin) <= fluxTolerance) then

! Do not set Converged to TRUE here as it messes up the MPI communication
! when repeating all of the angles
!       CSet% Converged(bin) = .TRUE.

       nConv = nConv + 1
     endif
   enddo

   if ( nConv == CSet% NumBin0 ) then 
     FluxConverged = .TRUE.
   endif


   return
   end subroutine testFluxConv 

