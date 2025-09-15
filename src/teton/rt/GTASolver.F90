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
!   20240305 - Comments from BCY:                                      *
!     I've added annotations marking where                             *
!     each term of the BiCGSTAB algorithm on Wikipedia is computed.    *
!     The matrix system being solved looks something like              *
!         (I-Linv*S)*PHI = Linv*(GTA%GreySource)                       *
!     where S is the grey scattering and Linv is the grey sweep        *
!     The action of the preconditioner is given by                     *
!          K^{-1} = I+Linv_{eps}S_{eps}                                *
!     where Linv_{eps} and S_{eps} are the stretched sweep and         *
!     scattering operators.                                            *
!     See the unpreconditioned BICGSTAB algorithm,
!        https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method#Unpreconditioned_BiCGSTAB
!     with A --> K^{-1}(I-Linv*S) and b --> K^{-1}b
!                                                                      *
!     Absent an extenral source, GreySweep applies the operator        *
!        I - (I+Linv_{eps}S_{eps})(I-Linv*S)                           *
!     which is equivalent to                                           *
!        Linv S + Linv_{eps}S_{eps} Linv S - Linv_{eps} S_{eps}        *
!     See rt/GreySweep.F90 for more details.                           *
!                                                                      *
!     The variables ending in "B" seem to                              *
!     be boundary information, though I'm not exactly sure.            *
!                                                                      *
!                                                                      *
!   Units:   E/e/T/m/L/A/V/t -                                         *
!        energy/photon energy/temperature/mass/length/area/volume/time *
!***********************************************************************

   subroutine GTASolver

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
   use ieee_arithmetic
   use, intrinsic :: iso_fortran_env, only : stdout=>output_unit, &
                                             stderr=>error_unit

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
   real(adqt) :: sumRad
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
   real(adqt), allocatable :: CGDirectionB(:,:)
   real(adqt), allocatable :: CGResidualB(:,:)
   real(adqt), allocatable :: CGActionB(:,:)
   real(adqt), allocatable :: CGActionSB(:,:)

!  Constants

   greyControl => getIterationControl(IterControls, "grey")
   nzones      =  Size%nzones

!  Allocate memory for the SI solution of the grey equations

   allocate( pzOld(nzones) )

!  Allocate memory for BiConjugate Gradient

   allocate( CGResidual(Size% ncornr) )
   allocate( CGDirection(Size% ncornr) )
   allocate( CGAction(Size% ncornr) )
   allocate( CGActionS(Size% ncornr) )
   allocate( CGDirectionB(Size% nbelem,Size% nangGTA) )
   allocate( CGResidualB(Size% nbelem,Size% nangGTA) )
   allocate( CGActionB(Size% nbelem,Size% nangGTA) )
   allocate( CGActionSB(Size% nbelem,Size% nangGTA) )

!  Initialize index of zone with maximum error:
   izRelErrPoint  = -1

!  Map data from the GPU

   if ( Size% useGPU ) then
     TOMP(target update from(Rad% PhiTotal))
     TOMP(target update from(GTA% GreySource))
     TOMP(target update to(GTA% Chi))
   endif

!  Sum current solution over groups for convergence test
!  Compute grey source

   ZoneLoop: do zone=1,nzones

     nCorner              = Geom% numCorner(zone)
     c0                   = Geom% cOffSet(zone)
     Rad% radEnergy(zone) = zero

     do c=1,nCorner
       sumRad               = sum( Rad% PhiTotal(:,c0+c) )
       Rad% radEnergy(zone) = Rad% radEnergy(zone) + Geom% Volume(c0+c)*sumRad
     enddo

     Rad% radEnergy(zone) = Rad% radEnergy(zone)/Geom% VolumeZone(zone)

   enddo ZoneLoop

!  Initialize Transport Matrices

   if (Size% useNewGTASolver) then

     if (Size% ndim == 2) then

!$omp parallel do default(none) schedule(static) &
!$omp& shared(nzones)
       do zone=1,nzones
         call InitGreySweepUCBrz(zone)
       enddo
!$omp end parallel do

     else

!$omp parallel do default(none) schedule(static) &
!$omp& shared(nzones)
       do zone=1,nzones
         call InitGreySweepUCBxyz(zone)
       enddo
!$omp end parallel do

     endif

   endif

!  Initialize the additive grey corrections, P, and CG
!  residual

   GTA%GreyCorrection(:) = zero
   ! x_0 = b
   pzOld(:)              = zero
   CGResidual(:)         = zero
   CGResidualB(:,:)      = zero
   ! \vec{\hat{r}}_0 = 1

!  Initialize the CG residual using an extraneous source

   nGreyIter  =  1
   withSource = .TRUE.

   if (Size% useNewGTASolver) then
     call GreySweepNEW(CGResidualB, CGResidual, withSource)
   else
     ! This computes CGResidual = \vec{r}_0 = \vec{\tilde{b}} = K_1^{-1} \vec{b}
     !   i.e., the action of the preconditioner on the initial right hand side
     ! We're starting with an initial guess of \vec{x}_0 = 0
     ! This is the only time that the GreySweep is called with a source.
     call GreySweep(CGResidualB, CGResidual)
   endif

!  Initialize the CG iteration.  Remove entries with zero scattering --
!  they live in the null space of M, where A := [I-M].

   CGDirection(:)    = CGResidual(:) ! \vec{p}_0 = \vec{r}_0
   CGDirectionB(:,:) = CGResidualB(:,:)

   rrProductOld   = scat_prod1(CGResidual) ! \rho_0 = \vec{r}_0 \cdot \vec{\hat{r}}_0

!  All CG sweeps are performed with zero extraneous source

   GTA%GreySource(:) = zero
   withSource        = .FALSE.

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

     CGAction(:)    = CGDirection(:)
     CGActionB(:,:) = CGDirectionB(:,:)

     if (Size% useNewGTASolver) then
       call GreySweepNEW(CGActionB, CGAction, withSource)
     else
       call GreySweep(CGActionB, CGAction)
     endif

!  Compute the action of the transport matrix, A, on the conjugate
!  direction.  Recall:  A := [I-M]

     ! CGAction = \vec{\nu} =  K_1^{-1} A \vec{p}_{i-1}
     CGAction(:)    = CGDirection(:)        - CGAction(:)
     CGActionB(:,:) = CGDirectionB(:,:) - CGActionB(:,:)

!    Compute the inner product, <d,Ad>

     dAdProduct = scat_prod1(CGAction)

!    Exit CG if the conjugate direction or the action of A on the
!    conjugate direction is zero

     if (abs(dAdProduct) < adqtSmall) then
       exit GreyIteration
     endif

     ! \alpha
     alphaCG = rrProductOld/dAdProduct

!    Update the residual
     ! \vec{s} = \vec{r}_{i-1} - \alpha \vec{\nu}
     CGResidual(:)    = CGResidual(:)    - alphaCG*CGAction(:)
     CGResidualB(:,:) = CGResidualB(:,:) - alphaCG*CGActionB(:,:)

     CGActionS(:)    = CGResidual(:)
     CGActionSB(:,:) = CGResidualB(:,:)

     if (Size% useNewGTASolver) then
       call GreySweepNEW(CGActionSB, CGActionS, withSource)
     else
       call GreySweep(CGActionSB, CGActionS)
     endif

!    Compute the action of the transport matrix, A, on the conjugate
!    direction.  Recall:  A := [I-M]

     ! \vec{t}
     CGActionS(:)    = CGResidual(:)    - CGActionS(:)
     CGActionSB(:,:) = CGResidualB(:,:) - CGActionSB(:,:)

     omegaNum = scat_prod(CGActionS,CGResidual)
     omegaDen = scat_prod(CGActionS,CGActionS)

     if (abs(omegaDen) < adqtSmall .or. abs(omegaNum) < adqtSmall) then
       GTA%GreyCorrection(:) = GTA%GreyCorrection(:) + alphaCG*CGDirection(:)

       exit GreyIteration
     endif

     ! \omega
     omegaCG = omegaNum/omegaDen

!    Update the Grey additive correction
     ! \vec{x}_i
     GTA%GreyCorrection(:) = GTA%GreyCorrection(:) +   &
                             alphaCG*CGDirection(:) + omegaCG*CGResidual(:)

     ! \vec{r}_i
     CGResidual(:)    = CGResidual(:)    - omegaCG*CGActionS(:)
     CGResidualB(:,:) = CGResidualB(:,:) - omegaCG*CGActionSB(:,:)

!    Compute the inner product, <r,r0>
     ! \rho_i
     rrProduct = scat_prod1(CGResidual)

     ! \beta
     betaCG = (rrProduct*alphaCG)/(rrProductOld*omegaCG)

!    update the conjugate direction
     ! \vec{p}_i
     CGDirection(:)    = CGResidual(:)  + betaCG*  &
                        (CGDirection(:) - omegaCG*CGAction(:))

     CGDirectionB(:,:) = CGResidualB(:,:)  + betaCG*  &
                        (CGDirectionB(:,:) - omegaCG*CGActionB(:,:))

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

       if (.not. ieee_is_finite(phiNew)) then
         izRelErrPoint  = zone  ! The zone where we first see a nan
         print *, "Teton's GTASolver encountered a bad value of ", phiNew, " on iteration", nGreyIter, " on rank ", Size% myRankInGroup, " in zone ", izRelErrPoint
         flush(stdout)
         TETON_FATAL("Grey solver encountered a NaN!")
       else if (abs(phiNew) > zero) then
         relErrPoint = abs(errZone/phiNew)
         if (relErrPoint > maxRelErrPoint) then
           maxRelErrPoint = relErrPoint
           izRelErrPoint  = zone
         endif
       endif

       pzOld(zone) = pz
     enddo CorrectionZoneLoop

     if (abs(phiL2) > zero) then
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

   if (Size% useNewGTASolver) then
     ngdart = ngdart + nGreyIter
   else
     ngdart = ngdart + 2*nGreyIter
   endif

   call setNumberOfIterations(greyControl,ngdart,.FALSE.)

!  Free memory

   deallocate(pzOld,        stat=alloc_stat)

   deallocate(CGResidual,   stat=alloc_stat)
   deallocate(CGDirection,  stat=alloc_stat)
   deallocate(CGAction,     stat=alloc_stat)
   deallocate(CGActionS,    stat=alloc_stat)
   deallocate(CGDirectionB, stat=alloc_stat)
   deallocate(CGResidualB,  stat=alloc_stat)
   deallocate(CGActionB,    stat=alloc_stat)
   deallocate(CGActionSB,   stat=alloc_stat)


   return
   end subroutine GTASolver

