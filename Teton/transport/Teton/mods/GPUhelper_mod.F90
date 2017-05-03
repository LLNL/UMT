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
   !logical(kind=1) :: fitsOnGPU = .true. ! default is false

   integer :: numGPUbuffers = 2 ! will be deterimined based on if it fits.
   !integer :: numGPUbuffers = 8 ! will be deterimined based on if it fits.

   ! create a for the GPU buffers
   type :: gpuStorage
      integer :: owner ! keeps track of which batch's data is currently held in this buffer
      real(adqt), device, allocatable :: data(:,:,:)
   end type gpuStorage


   integer, allocatable :: binSend(:), NangBin(:), anglebatch(:)

   ! Batchsize in number of angles (currently best to set to NangBin as some parts assume this for convenience)
   integer, parameter :: batchsize=32

   ! flag to determine if STime needs to be computed from tau*psi
   logical(kind=1) :: calcSTime
   ! flag to determine if scaleVolume needs to be done
   logical(kind=1) :: scaleVolume

   ! Cuda streams overlapping stuff
   ! (HtoD and DtoH combined stream, and kernel stream)
   integer(kind=cuda_stream_kind), save :: transfer_stream, kernel_stream
   integer :: Nbatches =  8 !QuadSet%NumBin
   type(cudaEvent), allocatable :: Psi_OnDevice(:), Psi_OnHost(:)
   type(cudaEvent), allocatable :: Psib_OnDevice(:), Psib_OnHost(:)
   type(cudaEvent), allocatable :: SweepFinished(:), STimeFinished(:)
   type(cudaEvent), allocatable :: ExitFluxDFinished(:)
   ! integer :: s, batch, istat, current, next

   ! zero copy pointers for phi and psib
   !type(C_DEVPTR)                    :: d_phi_p
   type(C_DEVPTR)                    :: d_psib_p
   type(C_DEVPTR)                    :: d_STime_p

   ! these buffers are allocated to fit on the device, either in double buffer batches, or full size if fits.
   type(gpuStorage), allocatable :: d_psi(:), d_STime(:)

   
   real(adqt), device, allocatable :: d_psibBatch(:,:,:,:)
   ! d_phi is full size, and persists on the device.
   real(adqt), device, allocatable :: d_phi(:,:)

   ! flags for marking if data is already present on the device:
   !   logical(kind=1), allocatable :: psi_present(:),STime_present(:)
   !   logical(kind=1), allocatable :: psib_present(:) ! this one is more of a dummy


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

    ! local variables
    integer :: NangBin_max
    integer :: buffer

    ! sweep helper variables
    allocate(binSend(numGPUbuffers), NangBin(numGPUbuffers), anglebatch(numGPUbuffers))

    ! flags for determining if data for a batch is already on the GPU:
    

    ! psib_present not needed now since psib always must be moved into the device. 
    ! may be useful later if GPU direct communication is used so psib is not moved in.
    !allocate(psib_present(Nbatches))


    ! events:

    ! allocate space for event recordings (necessary to overlap gpu movement with kernel compute)
    allocate( Psi_OnDevice(Nbatches), Psi_OnHost(Nbatches) )
    allocate( Psib_OnDevice(Nbatches), Psib_OnHost(Nbatches) )
    allocate( SweepFinished(Nbatches), STimeFinished(Nbatches) )
    allocate( ExitFluxDFinished(Nbatches) )


    ! device buffers:

    ! if the problem fits entirely in GPU memory, the buffers can be made the same size as the host arrays

    ! But if the problem does not fit in GPU memory, a double buffer strategy is used and the buffers have
    ! to be sized to fit chunks of the problem. (currently chunked by angle bins)

    ! allocate the needed number of copies of the GPU buffers
    allocate(d_psi(numGPUbuffers),d_STime(numGPUbuffers))
    

    do buffer=1, numGPUbuffers

       ! allocate the data in each buffer:
       allocate(d_psi(buffer)% data(QuadSet%Groups,Size%ncornr,BATCHSIZE) )
       allocate(d_STime(buffer)% data(QuadSet%Groups,Size%ncornr,BATCHSIZE) )
       
       ! and set the batch who owns the data in it to null (0)

       d_psi(buffer)% owner = 0
       d_STime(buffer)% owner = 0

    enddo

    ! allocate the other GPU arrays:
    allocate(d_psibBatch(QuadSet%Groups,Size%nbelem,BATCHSIZE,numGPUbuffers))
    allocate(d_phi(QuadSet%Groups,Size%ncornr))



    ! THERE ARE 2 STRATEGIES FOR OMEGA_A_FP, WE COULD CALCULATE EVERY TIME AND USE BUFFER SIZE OF 1 BIN, 
    ! OR WE COULD CALCULATE ONCE BUT HAVE TO STORE ALL 8 BINS. WILL KEEP 1 BIN STRATEGY FOR NOW, LATER CAN EXPERIMENT

    ! CURRENTLY DISABLED:
    if( .false. .and. FitsOnGPU ) then

       ! allocate omega_A_fp for all the angles--this will stay on the GPU and only needs to be calculated once.
       allocate( Geom%ZDataSoA % omega_A_fp(Size% nzones,Size% maxCorner,Size% maxcf, Size%nangSN) )
       allocate( Geom%ZDataSoA % omega_A_ez(Size% nzones,Size% maxCorner,Size% maxcf, Size%nangSN) )
   
    else ! Use just one bin for omega_A_fp and it will be recalculated every time.

       ! allocate omega_A_fp sections for batchsize (hardwired to NangBin here, only single buffer of omega_A_fp needed)
       NangBin_max = maxval(QuadSet%NangBinList(:))
       allocate( Geom%ZDataSoA % omega_A_fp(Size% nzones,Size% maxCorner,Size% maxcf, NangBin_max) )
       allocate( Geom%ZDataSoA % omega_A_ez(Size% nzones,Size% maxCorner,Size% maxcf, NangBin_max) )
   
    endif


  end subroutine InitDeviceBuffers




  subroutine CreateEvents()
    implicit none

    integer :: batch, istat

    do batch = 1, Nbatches
       istat = cudaEventCreate(Psi_OnDevice(batch))
       istat = cudaEventCreate(STimeFinished(batch))
       istat = cudaEventCreate(Psib_OnDevice(batch))
       istat = cudaEventCreate(Psib_OnHost(batch))
       istat = cudaEventCreate(SweepFinished(batch))
       istat = cudaEventCreate(Psi_OnHost(batch))
       istat = cudaEventCreate(ExitFluxDFinished(batch))
    enddo

  end subroutine CreateEvents


  subroutine checkDataOnDevice(d_storage, h_data, batch, buffer, mm1, numelements, streamid, event)

    implicit none
    
    !  Arguments

    type(gpuStorage), intent(inout) :: d_storage(:) ! d_psi, d_psib, or d_STime (may or not fit entire host buffer)
    real(adqt), intent(in) :: h_data(:,:,:) ! host buffer
    integer, intent(in) :: batch ! which data movement batch to check.
    integer, intent(in) :: buffer !used to select which of the spair buffers. Variable sent in will be current or next
    integer, intent(in) :: mm1 ! starting angle index within a bin. (will be 1 when batches are sized the same as angle bins)
    integer, intent(in) :: numelements ! number of array elements to be moved
    integer(kind=cuda_stream_kind), intent(in) :: streamid
    type(cudaEvent), intent(in) :: event
    ! local variables
    integer :: istat

    ! this routine checks whether the given batch already has its data on the device to determine if the move should occur. 

    ! if the data is already on the device, do not actually do a move.
    if(d_storage(buffer)%owner .eq. batch) then

       ! no ops
      
    else ! the data needs to be moved and ownership changed

       istat=cudaMemcpyAsync(d_storage(buffer)%data(1,1,1),                 &
            h_data(1,1,QuadSet%AngleOrder(mm1,binSend(buffer))), &
            numelements, streamid )

       ! change ownership: mark this batch's data as residing in this storage buffer.
       d_storage(buffer)%owner = batch

    endif

    ! Record when movement event finishes (for example, psi on device)
    istat=cudaEventRecord(event, streamid )

  end subroutine CheckDataOnDevice
    


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


  subroutine stageGPUData(buffer,batch,mm1)
    implicit none

    integer, intent(in) :: buffer, batch, mm1

    ! local variables

    integer istat

    ! buffer can be either current or next. Just a way of putting the data movement staging that happens 
    ! before and after the sweep into one reusable function. Buffer determines if you are staging it for the
    ! current buffer or the next buffer.

    !dummybatch = batch+1 should be sent in if using buffer=next

    if( FitsOnGPU ) then 
       ! If data fits on GPU do not worry about data movement. Just record that STime is ready:
       !istat=cudaEventRecord(STimeFinished(batch), transfer_stream )
    else !data does not fit on GPU and has to be streamed

       ! If this is first temp and intensity iteration, STime would have been calculated, needs update to host:
       if (calcSTime == .true.) then
          ! Wait for STime to be computed:
          istat = cudaStreamWaitEvent(transfer_stream, STimeFinished(batch), 0)
          ! Update STime to host
          istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend(buffer))), &
               d_STime(buffer)%data(1,1,1), &
               QuadSet%Groups*Size%ncornr*anglebatch(buffer), transfer_stream ) 

       else ! STime already computed,
          ! just need to move section of STime to device (Never called when FitsOnGPU)
          call checkDataOnDevice(d_STime, Geom%ZDataSoA%STime, batch, buffer, mm1, &
               QuadSet%Groups*Size%ncornr*anglebatch(buffer), transfer_stream, &
               STimeFinished(batch))
                   
       endif

    endif


    ! could select whether to recalculate fp each time. For now always do it.

    call fp_ez_c(     anglebatch(buffer),                     &
         Size%nzones,               &
         QuadSet%Groups,            &
         Size%ncornr,               &
         QuadSet%NumAngles,         &
         QuadSet%d_AngleOrder(mm1,binSend(buffer)),        & ! only need angle batch portion
         Size%maxCorner,            &
         Size%maxcf,                &
         Nangbin(buffer),                   &
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
