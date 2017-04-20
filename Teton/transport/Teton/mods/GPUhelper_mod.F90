! GPUhelper module: stores the arrays of unknowns used on the GPU

! these arrays of unknowns may fit on the GPU and persist accross iterations,
! Or they may be batched in. Either way the same name is used here, and they are
! declared allocatable here.

module GPUhelper_mod 

   use, intrinsic :: iso_c_binding
   use kind_mod
   use constant_mod
   use Size_mod
   use Quadrature_mod
   use Geometry_mod
   use ZoneData_mod

   use cudafor
   use nvtx_mod

   ! Flag to control if problem fits in GPU memory or has to be batched in.
   logical(kind=1) :: fitsOnGPU = .false. ! default is false

   integer :: numGPUbuffers = 2 ! will be deterimined based on if it fits.

   integer, allocatable :: binSend(:), NangBin(:), anglebatch(:)

   ! Batchsize in number of angles (currently best to set to NangBin as some parts assume this for convenience)
   integer, parameter :: batchsize=32

   ! flag to determine if STime needs to be computed from tau*psi
   logical(kind=1) :: calcSTime

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

   ! these buffers are allocated to fit on the device, either in double buffer batches, or full size if fits.
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


   ! Fortran to C interface for fp_ez_c and snswp3d_c subroutines.
   interface 

    subroutine fp_ez_c ( &
          anglebatch, &
          numzones, &
          numgroups, &
          ncornr, &
          numAngles, &
          AngleOrder, &
          maxcorners, &
          maxfaces, &
          nangbin, &
          nbelem, &
          omega, &
          numCorners, &
          numCFaces, &
          c0, &
          A_fp , &
          omega_A_fp , &
          A_ez , &
          omega_A_ez , &
          Connect , &
          Connect_reorder, &
          next, &
          nextZ, &   
          passZ, &
          streamid &
          ) &
          bind ( c ) 

         use iso_c_binding
         use cudafor
         integer ( c_int ) :: anglebatch
         integer ( c_int ) :: numzones
         integer ( c_int ) :: numgroups
         integer ( c_int ) :: ncornr
         integer ( c_int ) :: numAngles
         integer ( c_int ), device :: AngleOrder(*)
         integer ( c_int ) :: maxcorners
         integer ( c_int ) :: maxfaces
         integer ( c_int ) :: nangbin
         integer ( c_int ) :: nbelem
         real ( c_double ),device :: omega(*)
         integer ( c_int ),device :: numCorners(*) 
         integer ( c_int ),device :: numCFaces(*) 
         integer ( c_int ),device :: c0(*)
         real ( c_double ),device :: A_fp(*) 
         real ( c_double ),device :: omega_A_fp(*) 
         real ( c_double ),device :: A_ez(*) 
         real ( c_double ),device :: omega_A_ez(*) 
         integer ( c_int ),device :: Connect(*) 
         integer ( c_int ),device :: Connect_reorder(*) 
         integer ( c_int ),device :: next(*)
         integer ( c_int ),device :: nextZ(*)
         integer ( c_int ),device :: passZ(*)
         integer ( kind=cuda_stream_kind ),value :: streamid
       end subroutine fp_ez_c

       subroutine snswp3d_c ( &
          anglebatch, &
          numzones, &
          numgroups, &
          ncornr, &
          numAngles, &
          AngleOrder, &
          maxcorners, &
          maxfaces, &
          binRecv, &
          nangbin, &
          nbelem, &
          omega, &
          numCorners, &
          numCFaces, &
          c0, &
          A_fp , &
          omega_A_fp , &
          A_ez , &
          omega_A_ez , &
          Connect , &
          Connect_reorder, &
          STotal , &
          STimeBatch , &
          Volume , &
          psic, &
          psib, &
          next, &
          nextZ, &   
          Sigt, &
          SigtInv, &
          passZ, &
          calcSTime, &
          tau, &
          streamid &
          ) &
          bind ( c ) 

         use iso_c_binding
         use cudafor
         integer ( c_int ) :: anglebatch
         integer ( c_int ) :: numzones
         integer ( c_int ) :: numgroups
         integer ( c_int ) :: ncornr
         integer ( c_int ) :: numAngles
         integer ( c_int ), device :: AngleOrder(*)
         integer ( c_int ) :: maxcorners
         integer ( c_int ) :: maxfaces
         integer ( c_int ) :: binRecv
         integer ( c_int ) :: nangbin
         integer ( c_int ) :: nbelem
         real ( c_double ),device :: omega(*)
         integer ( c_int ),device :: numCorners(*) 
         integer ( c_int ),device :: numCFaces(*) 
         integer ( c_int ),device :: c0(*)
         real ( c_double ),device :: A_fp(*) 
         real ( c_double ),device :: omega_A_fp(*) 
         real ( c_double ),device :: A_ez(*) 
         real ( c_double ),device :: omega_A_ez(*) 
         integer ( c_int ),device :: Connect(*) 
         integer ( c_int ),device :: Connect_reorder(*) 
         real ( c_double ),device :: STotal(*) 
         real ( c_double ),device :: STimeBatch(*)
         real ( c_double ),device :: Volume(*) 
         real ( c_double ),device :: psic(*) 
         real ( c_double ),device :: psib(*) 
         integer ( c_int ),device :: next(*)
         integer ( c_int ),device :: nextZ(*)
         real ( c_double ),device :: Sigt(*) 
         real ( c_double ),device :: SigtInv(*) 
         integer ( c_int ),device :: passZ(*)
         logical ( c_bool ) :: calcSTime
         real ( c_double ) :: tau
         integer ( kind=cuda_stream_kind ),value :: streamid
       end subroutine snswp3d_c


    end interface



contains

  subroutine InitDeviceBuffers()
    implicit none

    integer :: NangBin_max


    ! sweep helper variables
    allocate(binSend(numGPUbuffers), NangBin(numGPUbuffers), anglebatch(numGPUbuffers))


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
       NangBin_max = maxval(QuadSet%NangBinList(:))
       allocate( Geom%ZDataSoA % omega_A_fp(Size% nzones,Size% maxCorner,Size% maxcf, NangBin_max) )
       allocate( Geom%ZDataSoA % omega_A_ez(Size% nzones,Size% maxCorner,Size% maxcf, NangBin_max) )
   
    endif

    ! THese are allocated regardless of size, but their size depends on numGPUbuffers.
    allocate(d_psi(QuadSet%Groups,Size%ncornr,BATCHSIZE,numGPUbuffers))
    allocate(d_STimeBatch(QuadSet%Groups,Size%ncornr,BATCHSIZE,numGPUbuffers))
    allocate(d_psibBatch(QuadSet%Groups,Size%nbelem,BATCHSIZE,numGPUbuffers))
    allocate(d_phi(QuadSet%Groups,Size%ncornr))



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


  subroutine MoveHtoD(d_buffer, h_buffer, bin, mm1, elements, streamid, event)

    implicit none
    
    !  Arguments

    real(adqt), device, intent(in) :: d_buffer(:,:,:,:) ! d_psi, d_psib, or d_STime (may or not fit entire host buffer)
    real(adqt), intent(in) :: h_buffer(:,:,:) ! host buffer
    integer, intent(in) :: bin !bin used to select which section of the buffer is used. What is sent in will be current or next
    integer, intent(in) :: mm1 ! starting angle index within a bin. (will be 1 when batches are sized the same as bins)
    integer, intent(in) :: elements ! number of array elements to be moved
    integer(kind=cuda_stream_kind), intent(in) :: streamid
    type(cudaEvent), intent(in) :: event
    ! local variables
    integer :: istat

    istat=cudaMemcpyAsync(d_buffer(1,1,1,bin),                 &
         h_buffer(1,1,QuadSet%AngleOrder(mm1,binSend(bin))), &
         elements, streamid )
    
    ! Record when movement event finishes (for example, psib on device)
    istat=cudaEventRecord(event, streamid )

  end subroutine MoveHtoD
    


  ! subroutine MoveDtoH(d_buffer, h_buffer, devicebin, hostbin, mm1, elements, streamid, event)
  !   ! NOT USED YET
  !   implicit none
    
  !   !  Arguments

  !   real(adqt), device, intent(in) :: d_buffer(:,:,:,:) ! d_psi, d_psib, or d_STime (may or not fit entire host buffer)
  !   real(adqt), intent(in) :: h_buffer(:,:,:) ! host buffer
  !   integer, intent(in) :: devicebin, hostbin !device bin might be current or next, hostbin will be an anglebin.
  !   integer, intent(in) :: mm1 ! starting angle index within a bin. (will be 1 when batches are sized the same as bins)
  !   integer, intent(in) :: elements ! number of array elements to be moved
  !   integer(kind=cuda_stream_kind), intent(in) :: streamid
  !   type(cudaEvent), intent(in) :: event
  !   ! local variables
  !   integer :: istat
    
  !   ! cuda memcpy syntax is cpy(to, from,...)
  !   istat=cudaMemcpyAsync(h_buffer(1,1,QuadSet%AngleOrder(mm1,hostbin)), &
  !        d_buffer(1,1,1,devicebin), &
  !        elements, streamid )
    
  !   ! Record when movement event finishes (for example, psib on device)
  !   istat=cudaEventRecord(event, streamid )

  ! end subroutine MoveDtoH


  subroutine stageGPUData(bin,batch,mm1)
    implicit none

    integer, intent(in) :: bin, batch, mm1

    ! local variables

    integer istat

    ! bin can be either current or next. Just a way of putting the data movement staging that happens 
    ! before and after the sweep into one reusable function. Bin determines if you are staging it for the
    ! current bin or the next bin.

    !batch+1 should be sent in if using next

               if (calcSTime == .true.) then
              ! Wait for STime to be computed:
              istat = cudaStreamWaitEvent(transfer_stream, STimeFinished(batch), 0)
              ! Update STime to host
              istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend(bin))), &
                   d_STimeBatch(1,1,1,bin), &
                   QuadSet%Groups*Size%ncornr*anglebatch(bin), transfer_stream ) 
           else
              ! STime already computed, just need to move section of STime to device
              call MoveHtoD(d_STimeBatch, Geom%ZDataSoA%STime, bin, mm1, &
                   QuadSet%Groups*Size%ncornr*anglebatch(bin), transfer_stream, STimeFinished(batch))

           endif


           call fp_ez_c(     anglebatch(bin),                     &
                Size%nzones,               &
                QuadSet%Groups,            &
                Size%ncornr,               &
                QuadSet%NumAngles,         &
                QuadSet%d_AngleOrder(mm1,binSend(bin)),        & ! only need angle batch portion
                Size%maxCorner,            &
                Size%maxcf,                &
                NangBin(bin),                   &
                Size%nbelem,                &
                QuadSet%d_omega,             &
                Geom%ZDataSoA%nCorner,                &
                Geom%ZDataSoA%nCFaces,                &
                Geom%ZDataSoA%c0,                &
                Geom%ZDataSoA%A_fp,                &
                Geom%ZDataSoA%omega_A_fp,                &
                Geom%ZDataSoA%A_ez,                &
                Geom%ZDataSoA%omega_A_ez,                &
                Geom%ZDataSoA%Connect,             &
                Geom%ZDataSoA%Connect_reorder,             &
                QuadSet%d_next,              &
                QuadSet%d_nextZ,             &
                QuadSet%d_passZstart,        &
                kernel_stream           &
                )

         end subroutine stageGPUData


end module GPUhelper_mod
