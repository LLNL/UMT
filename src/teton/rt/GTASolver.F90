#include "macros.h"
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

   subroutine GTASolver 

   use kind_mod
   use constant_mod
   use mpi_param_mod
   use mpif90_mod
   use iter_control_list_mod
   use iter_control_mod
   use Size_mod
   use Geometry_mod
   use GreyAcceleration_mod
   use ieee_arithmetic

   implicit none

!  Local

   type(IterControl),      pointer  :: greyControl => NULL()

   integer    :: c, c0, zone, nCorner, alloc_stat
   integer    :: g, Groups

   integer    :: ncornr, nzones 

   integer    :: nGreyIter, izRelErrPoint, ngdart 

   real(adqt) :: errL2, errZone, relErrPoint, maxRelErrPoint, &
                 relErrL2, phiL2, phiNew, sumRad, &
                 maxRelErrGreyLocal, maxRelErrGrey

   real(adqt) :: rrproduct, betaCG, alphaCG, omegaCG, &
                 rrproductold, dadproduct

   real(adqt) :: omegaNum, omegaDen, pz

   real(adqt), external :: scat_prod, scat_prod1

   logical(kind=1)      :: useGPU
   logical(kind=1)      :: withSource

   character(len=512)   :: descriptor

!  Dynamic

   real(adqt), allocatable :: pzOld(:)
   real(adqt), allocatable :: phitot(:)
   real(adqt), allocatable :: CGResidual(:)
   real(adqt), allocatable :: CGDirection(:)
   real(adqt), allocatable :: CGAction(:)
   real(adqt), allocatable :: CGActionS(:)

!  Set some constants

   greyControl => getIterationControl(IterControls, "grey")

   useGPU      =  getGPUStatus(Size)
   ncornr      =  Size%ncornr
   nzones      =  Size%nzones
   Groups      =  Size% ngr

!  Allocate memory for the SI solution of the grey equations

   allocate( pzOld(nzones) )
   allocate( phitot(nzones) )

!  Allocate memory for BiConjugate Gradient 

   allocate( CGResidual(ncornr) )
   allocate( CGDirection(ncornr) )
   allocate( CGAction(ncornr) )
   allocate( CGActionS(ncornr) )

!  Initialize index of zone with maximum error:
   izRelErrPoint  = -1

!  Sum current solution over groups for convergence test
!  Compute grey source

   phitot(:) = zero

   ZoneLoop: do zone=1,nzones

     nCorner = Geom% numCorner(zone) 
     c0      = Geom% cOffSet(zone) 

     do c=1,nCorner
       sumRad = zero
       do g=1,Groups
         sumRad = sumRad + Geom% PhiTotal(g,c0+c)
       enddo
       phitot(zone) = phitot(zone) + Geom% Volume(c0+c)*sumRad
     enddo

     phitot(zone) = phitot(zone)/Geom% VolumeZone(zone)

   enddo ZoneLoop

!  Initialize Transport Matrices

   if (Size% useNewGTASolver) then

     if ( useGPU ) then
       if (Size% ndim == 2) then
         call InitGreySweepUCBrz_GPU
       else
         call InitGreySweepUCBxyz_GPU
       endif
     else

       if (Size% ndim == 2) then

!$omp parallel do private(zone) schedule(static)
         do zone=1,nzones
           call InitGreySweepUCBrz(zone)
         enddo
!$omp end parallel do

       else

!$omp parallel do private(zone) schedule(static)
         do zone=1,nzones
           call InitGreySweepUCBxyz(zone)
         enddo
!$omp end parallel do

       endif

     endif

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

   if (Size% useNewGTASolver) then
     call GreySweepNEW(GTA%CGResidualB, CGResidual, withSource) 
   else
     call GreySweep(GTA%CGResidualB, CGResidual)
   endif

!  Initialize the CG iteration.  Remove entries with zero scattering --
!  they live in the null space of M, where A := [I-M].

   CGDirection(:)        = CGResidual(:)
   GTA%CGDirectionB(:,:) = GTA%CGResidualB(:,:)

   rrProductOld   = scat_prod1(CGResidual)

!  All CG sweeps are performed with zero extraneous source

   GTA%GreySource(:)    = zero
   withSource           = .FALSE.

!  Begin CG loop, iterating on grey corrections

   ngdart = getNumberOfIterations(greyControl) 

   GreyIteration: do

     ! This only does something if mod(verbose_level,10) > 2
     write(descriptor,'(A15,I5)') "GTASolver, GreyIteration number ", nGreyIter
     call PrintEnergies(trim(descriptor))

!  Exit CG if the residual is identically zero

     if (rrProductOld == zero) then
       if (nGreyIter <= 2) then
         GTA%GreyCorrection(:) = CGResidual(:)
       endif
       exit GreyIteration
     endif

!    increment the grey iteration counter
     nGreyIter = nGreyIter + 2 

!  Perform a transport sweep to compute the action of M on the 
!  conjugate direction (stored in CGAction)

     CGAction(:)        = CGDirection(:)
     GTA%CGActionB(:,:) = GTA%CGDirectionB(:,:) 

     if (Size% useNewGTASolver) then
       call GreySweepNEW(GTA%CGActionB, CGAction, withSource)
     else
       call GreySweep(GTA%CGActionB, CGAction) 
     endif

!  Compute the action of the transport matrix, A, on the conjugate
!  direction.  Recall:  A := [I-M]

     CGAction(:)        = CGDirection(:)        - CGAction(:)
     GTA%CGActionB(:,:) = GTA%CGDirectionB(:,:) - GTA%CGActionB(:,:)

!    Compute the inner product, <d,Ad>

     dAdProduct = scat_prod1(CGAction)

!    Exit CG if the conjugate direction or the action of A on the
!    conjugate direction is zero

     if (dAdProduct == zero) then
       exit GreyIteration
     endif

     alphaCG = rrProductOld/dAdProduct

!    Update the residual
     CGResidual(:)        = CGResidual(:)        - alphaCG*CGAction(:)
     GTA%CGResidualB(:,:) = GTA%CGResidualB(:,:) - alphaCG*GTA%CGActionB(:,:)

     CGActionS(:)        = CGResidual(:)
     GTA%CGActionSB(:,:) = GTA%CGResidualB(:,:)

     if (Size% useNewGTASolver) then
       call GreySweepNEW(GTA%CGActionSB, CGActionS, withSource) 
     else
       call GreySweep(GTA%CGActionSB, CGActionS)
     endif

!  Compute the action of the transport matrix, A, on the conjugate
!  direction.  Recall:  A := [I-M]

     CGActionS(:)        = CGResidual(:)        - CGActionS(:)
     GTA%CGActionSB(:,:) = GTA%CGResidualB(:,:) - GTA%CGActionSB(:,:)

     omegaNum = scat_prod(CGActionS,CGResidual)
     omegaDen = scat_prod(CGActionS,CGActionS)

     if (omegaDen == zero .or. omegaNum == zero) then
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

!  Compute the additive grey corrections on zones for convergence tests

     errL2          = zero
     phiL2          = zero
     maxRelErrPoint = zero

     do zone=1,nzones
       nCorner = Geom% numCorner(zone) 
       c0      = Geom% cOffSet(zone) 

!  Calculate the new zonal correction PZ

       pz = zero
       do c=1,nCorner
         pz = pz + Geom% Volume(c0+c)*GTA%GreyCorrection(c0+c)
       enddo
       pz = pz/Geom% VolumeZone(zone)

       errZone = pz - pzOld(zone)
       errL2   = errL2 + Geom% VolumeZone(zone)*(errZone*errZone)

       phiNew  = phitot(zone) + pz
       phiL2   = phiL2 + Geom% VolumeZone(zone)*(phiNew*phiNew)

       if (phiNew /= zero) then
         relErrPoint = abs(errZone/phiNew)
         if (relErrPoint > maxRelErrPoint) then
           maxRelErrPoint = relErrPoint
           izRelErrPoint  = zone 
         endif
       endif

       pzOld(zone) = pz
     enddo

     if (phiL2 /= zero) then
       relErrL2 = sqrt( abs(errL2/phiL2) )
     else
       relErrL2 = zero
     endif

     maxRelErrGreyLocal  = max(maxRelErrPoint,relErrL2)
     maxRelErrGrey       = maxRelErrGreyLocal

     call MPIAllReduce(maxRelErrGrey, "max", MY_COMM_GROUP)

!  Check convergence of the Grey Iteration

     if ( (maxRelErrGrey < getEpsilonPoint(greyControl) .or.         &
           nGreyIter >= getMaxNumberOfIterations(greyControl)) .and. &
           maxRelErrGrey < GTA%epsGrey ) then

       exit GreyIteration

     else if ( nGreyIter >= 100*getMaxNumberOfIterations(greyControl)) then

       ! Only print on offending ranks:
       if (maxRelErrGreyLocal >= GTA%epsGrey .or. ieee_is_nan(maxRelErrGreyLocal)) then
          print *, "Teton's GTASolver is not converging despite nGreyIter >= 100*getNumberOfMaxIterations! Maximum error on rank ", Size% myRankInGroup, " is ", maxRelErrPoint, " in zone ", izRelErrPoint
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

   call setNumberOfIterations(greyControl,ngdart)

!  Free memory

   deallocate(pzOld,         stat=alloc_stat)
   deallocate(phitot,        stat=alloc_stat)

   deallocate(CGResidual,    stat=alloc_stat)
   deallocate(CGDirection,   stat=alloc_stat)
   deallocate(CGAction,      stat=alloc_stat)
   deallocate(CGActionS,     stat=alloc_stat)


   return
   end subroutine GTASolver 

