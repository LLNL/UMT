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

   ! Flag to control if problem fits in GPU memory or has to be batched in.
   logical(kind=1) :: FitsOnGPU = .false. ! default is false

   ! Batchsize in number of angles (currently best to set to NangBin as some parts assume this for convenience)
   integer, parameter :: batchsize=32

   ! Cuda streams overlapping stuff
   ! (HtoD and DtoH combined stream, and kernel stream)
   integer(kind=cuda_stream_kind), save :: transfer_stream, kernel_stream
   integer :: Nbatches =  8 !QuadSet%NumBin
   type(cudaEvent), allocatable :: Psi_OnDevice(:), Psi_OnHost(:)
   type(cudaEvent), allocatable :: Psib_OnDevice(:), Psib_OnHost(:)
   type(cudaEvent), allocatable :: SweepFinished(:), STimeFinished(:)
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

  subroutine InitDeviceBuffers()
    implicit none

    integer :: NangBin

    ! EVENTS

    ! allocate space for event recordings (necessary to overlap gpu movement with kernel compute)
    allocate( Psi_OnDevice(Nbatches), Psi_OnHost(Nbatches) )
    allocate( Psib_OnDevice(Nbatches), Psib_OnHost(Nbatches) )
    allocate( SweepFinished(Nbatches), STimeFinished(Nbatches) )


    ! DEVICE BUFFERS

    ! if the problem fits entirely in GPU memory, the buffers can be made the same size as the host arrays

    ! But if the problem does not fit in GPU memory, a double buffer strategy is used and the buffers have
    ! to be sized to fit chunks of the problem. (currently chunked by angle bins)

    if( FitsOnGPU ) then


    else ! Does not fit on GPU, set up double buffer batching strategy:

       ! allocate omega_A_fp sections for batchsize (hardwired to NangBin here)
       NangBin = maxval(QuadSet%NangBinList(:))
       allocate( Geom%ZDataSoA % omega_A_fp(Size% nzones,Size% maxCorner,Size% maxcf, NangBin) )
       allocate( Geom%ZDataSoA % omega_A_ez(Size% nzones,Size% maxCorner,Size% maxcf, NangBin) )
    
    
       allocate(d_psi(QuadSet%Groups,Size%ncornr,BATCHSIZE,2))
       allocate(d_STimeBatch(QuadSet%Groups,Size%ncornr,BATCHSIZE,2))
       allocate(d_psibBatch(QuadSet%Groups,Size%nbelem,BATCHSIZE,2))
       allocate(d_phi(QuadSet%Groups,Size%ncornr))

    endif

  end subroutine InitDeviceBuffers

  subroutine CreateEvents()
    implicit none

    integer :: batch, istat

    do batch = 1, Nbatches
       istat = cudaEventCreate(Psi_OnDevice(batch))
       istat = cudaEventCreate(STimeFinished(batch))
       istat = cudaEventCreate(Psib_OnDevice(batch))
       istat = cudaEventCreate(SweepFinished(batch))
       istat = cudaEventCreate(Psi_OnHost(batch))
    enddo

  end subroutine CreateEvents


end module GPUhelper_mod
