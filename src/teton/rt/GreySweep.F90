!***********************************************************************
!                       Last Update:  10/2016, PFN                     *
!                                                                      *
!   GREYSWEEP - Performs a grey sweep with a fixed residual with       *
!               stretched TSA preconditioner. This gives the "action"  *
!               of the grey-transport operator on the input vector P.  *
!                                                                      *
!***********************************************************************

   subroutine GreySweepNEW(PsiB, P, withSource)

   use kind_mod
   use constant_mod
   use Size_mod
   use GreyAcceleration_mod

   implicit none

!  Arguments

   real(adqt), intent(inout)    :: PsiB(Size%nbelem,Size%nangGTA) 
   real(adqt), intent(inout)    :: p(Size%ncornr)

   logical (kind=1), intent(in) :: withSource

!  Local

   integer                     :: zone 
   real(adqt)                  :: wtiso 

   logical (kind=1)            :: useGPU
   logical (kind=1), parameter :: SnSweep=.FALSE.
   logical (kind=1), parameter :: savePsi=.FALSE.

!  Set some constants

   useGPU = getGPUStatus(Size)
   wtiso  = Size% wtiso

!  Perform a transport sweep to update the grey corrections, P

   GTA%TsaSource(:) = wtiso*( GTA%GreySource(:) + GTA%GreySigScat(:)*p(:) )

   call ControlSweep(SnSweep, savePsi, PsiB)

   if ( useGPU ) then
     if ( withSource ) then
       call ScalarIntensityDecompose_GPU(P)
     else
       call ScalarIntensitySolve_GPU(P)
     endif

   else

!$omp parallel do private(zone) shared(P, withSource) schedule(static)
     do zone=1,Size% nzones
       call ScalarIntensityDecompose(zone, P, withSource)
     enddo
!$omp end parallel do

   endif


   return
   end subroutine GreySweepNEW 


   subroutine GreySweep(apsib, P)

   use kind_mod
   use constant_mod
   use Size_mod
   use GreyAcceleration_mod

   implicit none

!  Arguments

   real(adqt), intent(inout) :: apsib(Size%nbelem,Size%nangGTA)
   real(adqt), intent(inout) :: p(Size%ncornr)

!  Local

   real(adqt)                  :: wtiso
   logical (kind=1), parameter :: SnSweep=.FALSE.
   logical (kind=1), parameter :: savePsi=.FALSE.

!  Set some constants

   wtiso = Size% wtiso

!  Perform a transport sweep to update the grey corrections, P

   GTA%TsaSource(:)         = wtiso*( GTA%GreySource(:) +  &
                              GTA%GreySigScat(:)*p(:) )

   GTA%OldGreyCorrection(:) = p(:)

   GTA%ID = 1
   call ControlSweep(SnSweep, savePsi, APSIB, P)

!  Now use "stretched" TSA as a preconditioner. 
!  The source to the TSA equations is the residual
!  from the grey transport sweep. The stretching parameter epsilon
!  has been chosen such that the scattering source vanishes.

   GTA%TsaCorrection(:) = zero
   GTA%TsaPsib(:,:)     = zero

   GTA%TsaSource(:)     = wtiso*GTA%eps(:)*GTA%GreySigScat(:)* &
                        ( p(:) - GTA%OldGreyCorrection(:) )

!  Perform a "stretched" transport sweep

   GTA%ID = 2
   call ControlSweep(SnSweep, savePsi, GTA%TsaPsib, GTA%TsaCorrection)

!  Update the grey corrections

   p(:) = p(:) + GTA%TsaCorrection(:)


   return
   end subroutine GreySweep

