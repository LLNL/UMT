!***********************************************************************
!                        Last Update:  04/2017, PFN                    *
!                                                                      *
!   ConstructSize - Builds the F90 module containing mesh-dependent    * 
!                   parameters.                                        *
!                                                                      *
!***********************************************************************

   subroutine ConstructSize(myRankInGroup,                           &
                            nzones, ncornr, nSides, nbelem,          &
                            maxcf, maxCorner, ncomm,                 &
                            ndim, ngr,                               &
                            functionRNLTE, tfloor,                   &
                            radForceMultiplier, betaNLTE, gammaNLTE, &
                            DopplerShiftOn,                          &
                            useNewNonLinearSolver, useNewGTASolver,  &
                            usePWLD, useSurfaceMassLumping,          & 
                            useGPU, useCUDASweep, useCUDASolver, zoneBatchSize,    &
                            nConcurrentBatches, igeom)               &  
                            BIND(C,NAME="teton_constructsize")

!  Include
   USE ISO_C_BINDING
   use kind_mod
   use Size_mod


   implicit none

!  Arguments

   integer(C_INT),  intent(in) :: myRankInGroup 
   integer(C_INT),  intent(in) :: nzones
   integer(C_INT),  intent(in) :: ncornr
   integer(C_INT),  intent(in) :: nSides
   integer(C_INT),  intent(in) :: nbelem
   integer(C_INT),  intent(in) :: maxcf
   integer(C_INT),  intent(in) :: maxCorner
   integer(C_INT),  intent(in) :: ncomm
   integer(C_INT),  intent(in) :: ndim
   integer(C_INT),  intent(in) :: ngr
   integer(C_INT),  intent(in) :: functionRNLTE

   real(C_DOUBLE),       intent(in) :: tfloor
   real(C_DOUBLE),       intent(in) :: radForceMultiplier
   real(C_DOUBLE),       intent(in) :: betaNLTE
   real(C_DOUBLE),       intent(in) :: gammaNLTE

   logical(C_BOOL), intent(in) :: DopplerShiftOn
   logical(C_BOOL), intent(in) :: useNewNonLinearSolver
   logical(C_BOOL), intent(in) :: useNewGTASolver
   logical(C_BOOL), intent(in) :: usePWLD
   logical(C_BOOL), intent(in) :: useSurfaceMassLumping
   logical(C_BOOL), intent(in) :: useGPU
   logical(C_BOOL), intent(in) :: useCUDASweep
   logical(C_BOOL), intent(in) :: useCUDASolver

   integer(C_INT),  intent(in) :: zoneBatchSize
   integer(C_INT),  intent(in) :: nConcurrentBatches

   integer(C_INT),  intent(in) :: igeom
    
!  Construct Run Parameters

   allocate (Size)

   call construct(Size, myRankInGroup=myRankInGroup,                 &
                        nzones=nzones,                               &
                        ncornr=ncornr,                               &
                        nSides=nSides,                               &
                        nbelem=nbelem,                               &
                        maxcf=maxcf,                                 &
                        maxCorner=maxCorner,                         &
                        ncomm=ncomm,                                 &
                        ndim=ndim,                                   &
                        ngr=ngr,                                     &
                        functionRNLTE=functionRNLTE,                 &
                        tfloor=tfloor,                               &
                        radForceMultiplier=radForceMultiplier,       &
                        betaNLTE=betaNLTE,                           &
                        gammaNLTE=gammaNLTE,                         & 
                        DopplerShiftOn=DopplerShiftOn,               &
                        useNewNonLinearSolver=useNewNonLinearSolver, &
                        useNewGTASolver=useNewGTASolver,             &
                        usePWLD=usePWLD,                             &
                        useSurfaceMassLumping=useSurfaceMassLumping, &
                        useGPU=useGPU,                               &
                        useCUDASweep=useCUDASweep,                   &
                        useCUDASolver=useCUDASolver,                 &
                        zoneBatchSize = zoneBatchSize,               &
                        nConcurrentBatches = nConcurrentBatches,     &
                        igeom=igeom)


   return
   end subroutine ConstructSize

