#include "macros.h"
#include "omp_wrappers.h"

!***********************************************************************
!                       Last Update:  10/2016, PFN                     *
!                                                                      *
!   GTASolver                                                          *
!     Controls acceleration in multidimensions. The acceleration       *
!   is grey and S2 and is solved using source iteration                *
!   pre-conditioned by "stretched" TSA. The system is solved with      *
!   a bi-conjugate gradient "stable" solver (BiCG-Stab).               *
!                                                                      *
!     The first "step" is the usual transport sweep:                   *
!                                                                      *
!       [OMEGA*DEL + SIGT] PSI(l+1/2) = SIGS*PHI(l) + Q                *
!                                                                      *
!   where l is the iteration index.  The second step is the solution   *
!   for the "corrections":                                             *
!                                                                      *
!       [OMEGA*DEL + (SIGT/EPS)] f(l+1) = (SIGT/EPS - EPS*SIGA)*F(l) + *
!                                       EPS*SIGS*(PHI(l+1/2) - PHI(l)) *
!                                                                      *
!       F(l+1) = Integral[f(l+1) d(OMEGA)]                             *
!                                                                      *
!   where F(l+1) is the additive correction and EPS is the             *
!  "stretching" parameter. Note that we choose the stretching          *
!   parameter so that the scattering term vanishes:                    *
!                                                                      *
!        SIGT/EPS - EPS*SIGA = 0 => EPS = 1/SQRT(1-C)                  *
!                                                                      *
!   where "C" is the scattering ratio.                                 *
!                                                                      *
!   The new scalar flux is given by:                                   *
!                                                                      *
!       PHI(l+1) = PHI(l+1/2) + F(l+1)                                 *
!                                                                      *
!                                                                      *
!   Units:   E/e/T/m/L/A/V/t -                                         *
!        energy/photon energy/temperature/mass/length/area/volume/time *
!***********************************************************************

   subroutine GTASolver_GPU

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use mpi_param_mod
   use mpif90_mod
   use iter_control_list_mod
   use iter_control_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use GreyAcceleration_mod
   use QuadratureList_mod
   use ZoneSet_mod
   use ieee_arithmetic

   implicit none

!  Local

   type(IterControl), pointer  :: greyControl => NULL()

   integer    :: c
   integer    :: c0
   integer    :: zone
   integer    :: nCorner
   integer    :: alloc_stat
   integer    :: nGreyIter
   integer    :: izRelErrPoint
   integer    :: ngdart
   integer    :: nzones
   integer    :: zSetID
   integer    :: nZoneSets

   real(adqt) :: errL2
   real(adqt) :: errZone
   real(adqt) :: relErrPoint
   real(adqt) :: maxRelErrPoint
   real(adqt) :: relErrL2
   real(adqt) :: phiL2
   real(adqt) :: phiNew
   real(adqt) :: pz
   real(adqt) :: maxRelErrGrey
   real(adqt) :: maxRelErrGreyLocal
   real(adqt) :: rrproduct
   real(adqt) :: betaCG
   real(adqt) :: alphaCG
   real(adqt) :: omegaCG
   real(adqt) :: rrproductold
   real(adqt) :: dadproduct
   real(adqt) :: omegaNum
   real(adqt) :: omegaDen

!  Note that there is some logic in these functions to floor to zero to avoid
!     underflow errors.  So scat_prod(..) == zero checks aren't as bad as they
!     seem.  Still, we'll use scat_prod(..) < adqtSmall in place of those checks
!     just to be safe.
   real(adqt), external :: scat_prod
   real(adqt), external :: scat_prod1

   logical(kind=1)      :: withSource

   character(len=512)   :: descriptor

!  Dynamic

   real(adqt), allocatable :: pzOld(:)
   real(adqt), allocatable :: CGResidual(:)
   real(adqt), allocatable :: CGDirection(:)
   real(adqt), allocatable :: CGAction(:)
   real(adqt), allocatable :: CGActionS(:)

!  Constants

   greyControl => getIterationControl(IterControls, "grey")

   nzones    = Size%nzones
   nZoneSets = getNumberOfZoneSets(Quad)

!  Allocate memory for BiConjugate Gradient

   allocate( pzOld(nzones) )



   allocate( CGResidual(Size% ncornr) )
   allocate( CGDirection(Size% ncornr) )
   allocate( CGAction(Size% ncornr) )
   allocate( CGActionS(Size% ncornr) )

!  Initialize index of zone with maximum error:
   izRelErrPoint  = -1

!  Sum current solution over groups for convergence test
!  Compute grey source

   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
   TOMPC(shared(nZoneSets, ZSet, Geom, Rad))

   do zSetID=1,nZoneSets

!$omp parallel do default(none) schedule(dynamic)  &
!$omp& shared(zSetID, ZSet, Geom, Rad)

     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       ZSet% sumT(c) = Geom% Volume(c)*sum( Rad% PhiTotal(:,c) )
     enddo
!$omp end parallel do

   enddo

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, Geom, Rad, ZSet)&)
   TOMPC(private(c0, nCorner))

   do zSetID=1,nZoneSets

!$omp parallel do default(none) schedule(dynamic)  &
!$omp& shared(zSetID, Geom, Rad, ZSet) &
!$omp& private(c0, nCorner)

     do zone=Geom% zone1(zSetID),Geom% zone2(zSetID)
       nCorner              = Geom% numCorner(zone)
       c0                   = Geom% cOffSet(zone)
       Rad% radEnergy(zone) = zero

       do c=1,nCorner
         Rad% radEnergy(zone) = Rad% radEnergy(zone) + ZSet% sumT(c0+c)
       enddo

       Rad% radEnergy(zone) = Rad% radEnergy(zone)/Geom% VolumeZone(zone)
     enddo
!$omp end parallel do

   enddo

   TOMP(end target teams distribute)

   TOMP(target update from(Rad% radEnergy))

!  Initialize Transport Matrices

   if (Size% ndim == 2) then
     call InitGreySweepUCBrz_GPU
   else
     call InitGreySweepUCBxyz_GPU
   endif

!  Initialize the additive grey corrections, P, and CG
!  residual

   GTA%GreyCorrection(:) = zero
   pzOld(:)              = zero
   CGResidual(:)         = zero
   GTA%CGResidualB(:,:)  = zero

!  Initialize the CG residual using an extraneous source

   nGreyIter             =  1
   withSource            = .TRUE.
   GTA% nGreySweepIters  =  2

   call GreySweepNEW(GTA%CGResidualB, CGResidual, withSource)

!  Initialize the CG iteration.  Remove entries with zero scattering --
!  they live in the null space of M, where A := [I-M].

   CGDirection(:)        = CGResidual(:)
   GTA%CGDirectionB(:,:) = GTA%CGResidualB(:,:)

   rrProductOld   = scat_prod1(CGResidual)

!  All CG sweeps are performed with zero extraneous source

   withSource = .FALSE.

!  Begin CG loop, iterating on grey corrections

   ngdart = getNumberOfIterations(greyControl)

   GreyIteration: do

     ! This only does something if mod(verbose_level,10) > 2
     write(descriptor,'(A15,I5)') "GTASolver, GreyIteration number ", nGreyIter
     call PrintEnergies(trim(descriptor))

!    Exit CG if the residual is below the minimum. This used to test against zero,
!    but due to differences in rounding errors some platforms would return
!    very small numbers and not zero.

     if (abs(rrProductOld) < adqtSmall) then
       if (nGreyIter <= 2) then
         GTA%GreyCorrection(:) = CGResidual(:)
       endif
       exit GreyIteration
     endif

!    increment the grey iteration counter
     nGreyIter = nGreyIter + 2

!    Perform a transport sweep to compute the action of M on the
!    conjugate direction (stored in CGAction)

     CGAction(:)        = CGDirection(:)
     GTA%CGActionB(:,:) = GTA%CGDirectionB(:,:)

     call GreySweepNEW(GTA%CGActionB, CGAction, withSource)

!    Compute the action of the transport matrix, A, on the conjugate
!    direction.  Recall:  A := [I-M]

     CGAction(:)        = CGDirection(:)        - CGAction(:)
     GTA%CGActionB(:,:) = GTA%CGDirectionB(:,:) - GTA%CGActionB(:,:)

!    Compute the inner product, <d,Ad>

     dAdProduct = scat_prod1(CGAction)

!    Exit CG if the conjugate direction or the action of A on the
!    conjugate direction is zero

     if (abs(dAdProduct) < adqtSmall) then
       exit GreyIteration
     endif

     alphaCG = rrProductOld/dAdProduct

!    Update the residual
     CGResidual(:)        = CGResidual(:)        - alphaCG*CGAction(:)
     GTA%CGResidualB(:,:) = GTA%CGResidualB(:,:) - alphaCG*GTA%CGActionB(:,:)

     CGActionS(:)        = CGResidual(:)
     GTA%CGActionSB(:,:) = GTA%CGResidualB(:,:)

     call GreySweepNEW(GTA%CGActionSB, CGActionS, withSource)

!    Compute the action of the transport matrix, A, on the conjugate
!    direction.  Recall:  A := [I-M]

     CGActionS(:)        = CGResidual(:)        - CGActionS(:)
     GTA%CGActionSB(:,:) = GTA%CGResidualB(:,:) - GTA%CGActionSB(:,:)

     omegaNum = scat_prod(CGActionS,CGResidual)
     omegaDen = scat_prod(CGActionS,CGActionS)

     if (abs(omegaDen) < adqtSmall .or. abs(omegaNum) < adqtSmall) then
       GTA%GreyCorrection(:) = GTA%GreyCorrection(:) + alphaCG*CGDirection(:)

       exit GreyIteration
     endif

     omegaCG = omegaNum/omegaDen

!    Update the Grey additive correction
     GTA%GreyCorrection(:) = GTA%GreyCorrection(:) +   &
                             alphaCG*CGDirection(:) + omegaCG*CGResidual(:)

     CGResidual(:)        = CGResidual(:)        - omegaCG*CGActionS(:)
     GTA%CGResidualB(:,:) = GTA%CGResidualB(:,:) - omegaCG*GTA%CGActionSB(:,:)

!    Compute the inner product, <r,r0>
     rrProduct = scat_prod1(CGResidual)

     betaCG = (rrProduct*alphaCG)/(rrProductOld*omegaCG)

!    update the conjugate direction
     CGDirection(:)        = CGResidual(:)  + betaCG*  &
                            (CGDirection(:) - omegaCG*CGAction(:))

     GTA%CGDirectionB(:,:) = GTA%CGResidualB(:,:)  + betaCG*  &
                            (GTA%CGDirectionB(:,:) - omegaCG*GTA%CGActionB(:,:))

!    Compute the additive grey corrections on zones for convergence tests

     errL2          = zero
     phiL2          = zero
     maxRelErrPoint = zero

     CorrectionZoneLoop: do zone=1,nzones
       nCorner = Geom% numCorner(zone)
       c0      = Geom% cOffSet(zone)

!      Calculate the new zonal correction PZ

       pz = zero
       do c=1,nCorner
         pz = pz + Geom% Volume(c0+c)*GTA%GreyCorrection(c0+c)
       enddo
       pz = pz/Geom% VolumeZone(zone)

       errZone = pz - pzOld(zone)
       errL2   = errL2 + Geom% VolumeZone(zone)*(errZone*errZone)

       phiNew  = Rad% radEnergy(zone) + pz
       phiL2   = phiL2 + Geom% VolumeZone(zone)*(phiNew*phiNew)

       ! Is this too cumbersome? a NaN check on every zone on every GTA iteration?
       if (ieee_is_nan(phiNew) .or. ieee_is_nan(errZone)) then
         izRelErrPoint  = zone  ! The zone where we first see a nan
         print *, "Teton's GTASolver encountered a NaN on iteration", nGreyIter, " on rank ", Size% myRankInGroup, " in zone ", izRelErrPoint
         call sleep(15)
         TETON_FATAL("Grey solver encountered a NaN!")
       else if (phiNew /= zero) then
         relErrPoint = abs(errZone/phiNew)
         if (relErrPoint > maxRelErrPoint) then
           maxRelErrPoint = relErrPoint
           izRelErrPoint  = zone
         endif
       endif

       pzOld(zone) = pz
     enddo CorrectionZoneLoop

     if (phiL2 /= zero) then
       relErrL2 = sqrt( abs(errL2/phiL2) )
     else
       relErrL2 = zero
     endif

     maxRelErrGreyLocal  = max(maxRelErrPoint,relErrL2)
     maxRelErrGrey       = maxRelErrGreyLocal

     call MPIAllReduce(maxRelErrGrey, "max", MY_COMM_GROUP)

!    Check convergence of the Grey Iteration

     if ( GTA% enforceHardGTAIterMax .and. nGreyIter >= getMaxNumberOfIterations(greyControl) ) then

       exit GreyIteration

     else if ( (maxRelErrGrey < getEpsilonPoint(greyControl) .or. &
           nGreyIter >= getMaxNumberOfIterations(greyControl)) .and. &
           maxRelErrGrey < GTA%epsGrey ) then

       exit GreyIteration

     else if ( nGreyIter >= 100*getMaxNumberOfIterations(greyControl)) then

       ! Only print on offending ranks:
       if (maxRelErrGreyLocal >= GTA%epsGrey) then
          print *, "Teton's GTASolver is not converging despite nGreyIter ", nGreyIter, " >= 100*getNumberOfMaxIterations! Maximum error on rank ", Size% myRankInGroup, " is ", maxRelErrPoint, " in zone ", izRelErrPoint
       endif

       ! Provide enough time for the above statement to get printed on every rank
       call sleep(15)

       TETON_FATAL("Grey solver is not converging, has exceeded iteration control's max # iterations * 100")

     else

       rrProductOld = rrProduct
       cycle GreyIteration

     endif

   enddo GreyIteration

   call PrintEnergies("GTASolver, after end of GreyIteration")

   ngdart = ngdart + nGreyIter

   call setNumberOfIterations(greyControl,ngdart)

!  Free memory

   deallocate(pzOld,         stat=alloc_stat)

   deallocate(CGResidual,    stat=alloc_stat)
   deallocate(CGDirection,   stat=alloc_stat)
   deallocate(CGAction,      stat=alloc_stat)
   deallocate(CGActionS,     stat=alloc_stat)


   return
   end subroutine GTASolver_GPU

