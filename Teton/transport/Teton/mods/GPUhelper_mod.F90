! GPUhelper module: stores the arrays of unknowns used on the GPU

! these arrays of unknowns may fit on the GPU and persist accross iterations,
! Or they may be batched in. Either way the same name is used here, and they are
! declared allocatable here.

module GPUhelper_mod 

   use kind_mod
   use constant_mod
   use Size_mod
   use Quadrature_mod
   use Geometry_mod
   use ZoneData_mod

   use cudafor
   use nvtx_mod


   ! Cuda streams overlapping stuff
   ! (HtoD and DtoH combined stream, and kernel stream)
   ! integer :: Nbatches =  8 !QuadSet%NumBin
   ! integer(kind=cuda_stream_kind), save :: transfer_stream, kernel_stream
   ! type(cudaEvent) :: Psi_OnDevice(Nbatches), Psi_OnHost(Nbatches)
   ! type(cudaEvent) :: Psib_OnDevice(Nbatches), Psib_OnHost(Nbatches)
   ! type(cudaEvent) :: SweepFinished(Nbatches), STimeFinished(Nbatches)
   ! integer :: s, batch, istat, current, next

   ! zero copy pointers for phi and psib
   !type(C_DEVPTR)                    :: d_phi_p
   type(C_DEVPTR)                    :: d_psib_p
   type(C_DEVPTR)                    :: d_STime_p



   !real(adqt), device, allocatable :: d_psib(:,:,:)
   !real(adqt), device, allocatable :: d_STime(:,:,:)


   real(adqt), device, allocatable :: d_psi(:,:,:,:)
   real(adqt), device, allocatable :: d_STimeBatch(:,:,:,:)
   real(adqt), device, allocatable :: d_psibBatch(:,:,:,:)
   ! d_phi is full size, and persists on the device.
   real(adqt), device, allocatable :: d_phi(:,:)

   ! real(adqt), device :: d_psi(QuadSet%Groups,Size%ncornr,BATCHSIZE,2)
   ! real(adqt), device :: d_STimeBatch(QuadSet%Groups,Size%ncornr,BATCHSIZE,2)
   ! real(adqt), device :: d_psibBatch(QuadSet%Groups,Size%nbelem,BATCHSIZE,2)
   ! ! d_phi is full size, and persists on the device.
   ! real(adqt), device :: d_phi(QuadSet%Groups,Size%ncornr)


contains

end module GPUhelper_mod
