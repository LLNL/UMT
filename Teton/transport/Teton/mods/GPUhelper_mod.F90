! GPUhelper module: stores the arrays of unknowns used on the GPU

! these arrays of unknowns may fit on the GPU and persist accross iterations,
! Or they may be batched in. Either way the same name is used here, and they are
! declared allocatable here.

module GPUhelper_mod 

   use, intrinsic :: iso_c_binding
   use kind_mod
   use constant_mod
   !use Size_mod
   use Quadrature_mod
   use Geometry_mod
   use ZoneData_mod

   use cudafor
   use nvtx_mod

   ! Flag to control if problem fits in GPU memory or has to be batched in.
   logical(kind=1) :: fitsOnGPU = .false. ! default is false
   !logical(kind=1) :: fitsOnGPU = .true. ! default is false

   !integer :: numPsi_buffers = 8 ! will be deterimined based on if it fits.
   integer :: numPsi_buffers = 2 ! double buffers for psi
   integer :: numSTime_buffers = 1 ! only need 1 STime buffer the way the code is written
   integer :: numPsib_buffers = 1 !should only need 1 psib buffer
   integer :: numDot_buffers = 1 ! for omega A ez and fp precomputed dot products.


   real(adqt) :: Stime_temp
   real(adqt) :: volumeRatio_temp

   ! create a type for the GPU buffers (3 dimensional arrays)
   type :: gpuStorage
      integer :: slot         ! which storage slot is this
      integer :: owner ! keeps track of which batch's data is currently held in this buffer
      real(adqt), device, allocatable :: data(:,:,:)
   end type gpuStorage

   ! type for buffer tracking of buffers with 4 dimensions
   type :: dotStorage
      integer :: slot         ! which storage slot is this
      integer :: owner        ! which bin is in this storage slot
      !integer :: numelements  ! size of data in array elements
      real(adqt), device, allocatable :: data(:,:,:,:) ! size: nZ*mC*mF*anglebatch
   end type dotStorage

   type :: bufferdata
      integer :: bin             ! bin being swept
      integer :: batch           ! sequentual ordering, 1,2,3... (may not be needed)
      integer :: NangBin         ! number of angles in this bin
      integer :: anglebatch      ! number of angles in this batch (these are the same for now, later handle funny size anglebatchs)
      !integer :: Angles
      type(gpuStorage), pointer  :: psi
      type(gpuStorage), pointer  :: STime
      type(gpuStorage), pointer  :: psib
      type(dotStorage), pointer  :: omega_A_fp
      type(dotStorage), pointer  :: omega_A_ez
   end type bufferdata

   ! buffer managment stuff, we have a current buffer and an other buffer (previous or next)
   type(bufferdata) :: previous, current, next ! I could have just current and other, but I think this will be confusing.

   type(gpuStorage), target, allocatable :: psi_storage(:), psib_storage(:), STime_storage(:)
   type(dotStorage), target, allocatable :: omega_A_ez_storage(:), omega_A_fp_storage(:)

   ! Cuda streams 
   integer :: s, istat

   !integer, allocatable :: batch(:), binSend(:), NangBin(:), anglebatch(:)

   !integer :: previous_batch, previous_binSend, previous_buffer

   ! Batchsize in number of angles (currently best to set to NangBin as some parts assume this for convenience)
   !integer, parameter :: batchsize=32
   integer :: batchsize

   ! flag to determine if STime needs to be computed from tau*psi
   logical(kind=1) :: calcSTime

   ! Cuda streams overlapping stuff
   ! (HtoD and DtoH combined stream, and kernel stream) stream1 is primary stream
   integer(kind=cuda_stream_kind), save :: transfer_stream, kernel_stream, transfer_stream1 
   integer :: Nbatches =  8 !This is the max number of batches (Total number of angle bins that might be swept)
   ! movement events:
   type(cudaEvent), allocatable :: Psi_OnDevice(:), Psi_OnHost(:)
   type(cudaEvent), allocatable :: Psib_OnDevice(:), Psib_OnHost(:)
   type(cudaEvent) :: phi_OnHost
   ! kernel events
   type(cudaEvent), allocatable :: SweepFinished(:), STimeFinished(:)
   type(cudaEvent), allocatable :: ExitFluxDFinished(:), AfpFinished(:)
   type(cudaEvent), allocatable :: snmomentsFinished(:)

   ! these buffers are allocated to fit on the device, either in double buffer batches, or full size if fits.
   !type(gpuStorage), allocatable :: d_psi(:), d_STime(:)

   ! zero copy pointers for psib
   type(C_DEVPTR)                    :: d_psib_p
   !real(adqt), device, allocatable :: pinned_psib(:,:,:)
   real(adqt), device, pointer :: pinned_psib(:,:,:)
   
   ! device resident batch size of psib
   !real(adqt), device, allocatable :: d_psibBatch(:,:,:)
   ! d_phi is full size, and persists on the device.
   real(adqt), device, allocatable :: d_phi(:,:)

   ! create device versions
   !real(adqt), device, allocatable :: d_omega_A_fp(:,:,:,:) ! size: nZ*mC*mF*anglebatch
   !real(adqt), device, allocatable :: d_omega_A_ez(:,:,:,:) ! size: nZ*mC*mF*anglebatch


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



! Interface to CUDA-C api for setting 

INTERFACE

    INTEGER FUNCTION CUDAFUNCSETATTRIBUTE(func, attr, avalue) &

            BIND(C, NAME='cudaFuncSetAttribute')

      EXTERNAL func

      INTEGER, value :: attr

      INTEGER, value :: avalue

    END FUNCTION

END INTERFACE

 

enum, bind(c)

  enumerator :: cudaFuncAttributeMaxDynamicSharedMemorySize = 8

  !enumerator cudaFuncAttributePreferredSharedMemoryCarveout

  !enumerator cudaFuncAttributeMax

end enum

 

enum, bind(c)

  enumerator :: cudaSharedmemCarveoutDefault   = -1

  enumerator :: cudaSharedmemCarveoutMaxShared = 100

end enum

contains


  attributes(global) subroutine GPU_fp_ez_hplane_f( &
                              anglebatch,                     &
                              nzones,               &
                              ncornr,               &
                              NumAngles,         &
                              AngleOrder,        & ! only angle batch portion
                              maxCorner,            &
                              maxcf,                &
                              NangBin,                   &
                              nbelem,                &
                              ZData,       &
                              omega,             &
                              omega_A_fp,                &
                              omega_A_ez,                &
                              next,              &
                              nextZ,             &
                              passZstart             )
    implicit none

    !  Arguments

   integer, value,    intent(in)    :: anglebatch
   integer, value,    intent(in)    :: nzones
   integer, value,    intent(in)    :: ncornr
   integer, value,    intent(in)    :: NumAngles

   integer,    device, intent(in) :: AngleOrder(anglebatch)

   integer, value,    intent(in)    :: maxCorner
   integer, value,    intent(in)    :: maxcf
   integer, value,    intent(in)    :: NangBin
   integer, value,    intent(in)    :: nbelem

   type(GPU_ZoneData), device, intent(in) :: ZData(nzones)

   real(adqt), device, intent(in)    :: omega(3,NumAngles)

   real(adqt), device, intent(out)    :: omega_A_fp(maxcf,maxCorner, nzones, anglebatch) 
   real(adqt), device, intent(out)    :: omega_A_ez(maxcf,maxCorner,nzones, anglebatch) 

   integer,    device, intent(in) :: next(ncornr+1,NumAngles)
   integer,    device, intent(in) :: nextZ(nzones,NumAngles)
   integer,    device, intent(in) :: passZstart(nzones,NumAngles)   

    !  Local Variables

    integer    :: Angle, i, ib, ic, icfp, icface
    integer    :: zone, c, cez, ii, mm, ndone
    integer    :: p, ndoneZ, passZcount
    integer    :: nCorner, nCFaces, c0

    !  Constants

    mm = blockIdx%x
    Angle = AngleOrder(mm)

    p = 0
    ndoneZ = 0
    PassLoop: do while (ndoneZ < nzones)
       p = p + 1
       ! number of zones in this hyperplane:
       passZcount = passZstart(p+1,Angle) - passZstart(p,Angle)
       

       ZoneLoop: do ii=threadIdx%x,passZcount,blockDim%x

          !!FIXME: simplifying assumption that all zones have same nCorner values
          !! (they're all 8 from what we've seen). If this isn't true in general,
          !! just convert this into a table lookup
          ndone = (ndoneZ+ii-1) * maxCorner

          zone = nextZ(ndoneZ+ii,Angle)

          !nCorner = nCornerArray(zone)
          !nCFaces =   nCFacesArray(zone)
          !c0      =   c0Array(zone)

          !write(0,*) zone
          !print *, "zone = ", zone
          nCorner =   ZData(zone)%nCorner
          nCFaces =   ZData(zone)%nCFaces
          c0      =   ZData(zone)%c0


          CornerLoop: do i=1,nCorner

             ic      = next(ndone+i,Angle)
             c       = ic - c0

             !  Calculate Area_CornerFace dot Omega to determine the 
             !  contributions from incident fluxes across external 
             !  corner faces (FP faces)

             do icface=1,ncfaces

                omega_A_fp(icface,c,zone,blockIdx%x) = omega(1,Angle)*ZData(zone)%A_fp(1,icface,c) + &
                     omega(2,Angle)* ZData(zone)%A_fp(2,icface,c) + &
                     omega(3,Angle)* ZData(zone)%A_fp(3,icface,c)
                
                ! Reorder is now done when GPU_ZData is initialized, so following is not needed.
                !icfp    =  Connect(1,icface,c,zone)
                !ib      =  Connect(2,icface,c,zone)
                !cez     =  Connect(3,icface,c,zone)
                
                !Connect_ro(icface,c,1,zone) = icfp
                !Connect_ro(icface,c,2,zone) = ib
                !Connect_ro(icface,c,3,zone) = cez
                
             enddo

             !  Contributions from interior corner faces (EZ faces)

             do icface=1,nCFaces

                omega_A_ez(icface,c,zone,blockIdx%x) = omega(1,Angle)* ZData(zone)%A_ez(1,icface,c) + &
                     omega(2,Angle)* ZData(zone)%A_ez(2,icface,c) + &
                     omega(3,Angle)* ZData(zone)%A_ez(3,icface,c) 

             enddo

          enddo CornerLoop

       enddo ZoneLoop

       ndoneZ = ndoneZ + passZcount

       call syncthreads

    enddo PassLoop

  end subroutine GPU_fp_ez_hplane_f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! fp_ez Caller
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine fp_ez_f ( &
          anglebatch, &
          nzones, &
          ncornr, &
          numAngles, &
          d_AngleOrder, &
          maxCorner, &
          maxcf, &
          NangBin, &
          nbelem, &
          d_omega, &
          d_omega_A_fp , &
          d_omega_A_ez , &
          d_next, &
          d_nextZ, &   
          d_passZstart, &
          streamid &
          ) 


   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use Quadrature_mod
   use Material_mod
   use ZoneData_mod
   use cudafor

   implicit none

!  Arguments
   integer, parameter :: ndim=3

   integer,    intent(in)    :: anglebatch
   integer,    intent(in)    :: nzones
   integer,    intent(in)    :: ncornr
   integer,    intent(in)    :: NumAngles

   integer,    device, intent(in) :: d_AngleOrder(anglebatch)

   integer,    intent(in)    :: maxCorner
   integer,    intent(in)    :: maxcf
   integer,    intent(in)    :: NangBin
   integer,    intent(in)    :: nbelem
   real(adqt), device, intent(in)    :: d_omega(3,NumAngles)

   real(adqt), device, intent(out)    :: d_omega_A_fp(maxcf ,maxCorner, nzones, anglebatch) 
   real(adqt), device, intent(out)    :: d_omega_A_ez(maxcf ,maxCorner, nzones, anglebatch)  


   integer,    device, intent(in) :: d_next(Size%ncornr+1,QuadSet%NumAngles)
   integer,    device, intent(in) :: d_nextZ(Size%nzones,QuadSet%NumAngles)
   integer,    device, intent(in) :: d_passZstart(Size%nzones,QuadSet%NumAngles)   
   integer(kind=cuda_stream_kind), intent(in) :: streamid   

!  Local Variables

   integer    :: mm, Angle,istat,i,ib,ic

   type(dim3) :: threads,blocks
   
   !integer    :: shmem !amount of shared memory need by GPU sweep kernel.

   
   !threads=dim3(QuadSet%Groups,NZONEPAR,1) 
   threads=dim3(128,1,1) 
   blocks=dim3(anglebatch,1,1)

   call GPU_fp_ez_hplane_f<<<blocks,threads,0,streamid>>>( anglebatch,                     &
                              nzones,               &
                              ncornr,               &
                              NumAngles,         &
                              d_AngleOrder,        & ! only angle batch portion
                              maxCorner,            &
                              maxcf,                &
                              NangBin,                   &
                              nbelem,                &
                              Geom%d_GPU_ZData,    &
                              d_omega,             &
                              d_omega_A_fp,                &
                              d_omega_A_ez,                &
                              d_next,              &
                              d_nextZ,             &
                              d_passZstart             )



   return
 end subroutine fp_ez_f



  subroutine GPUmemRequirements(psib,psi,phi,STime,next,omega_a_fp,sweep_mem)
    implicit none
    
    real(adqt), intent(in) :: psib(:,:,:), psi(:,:,:), phi(:,:), STime(:,:,:), omega_a_fp(:,:,:,:)
    integer, intent(in) :: next(:,:)

    integer(kind=8), intent(out) :: sweep_mem 
    
    
    integer(kind=8) :: psib_mem, psi_mem, phi_mem, stime_mem, next_mem, omega_a_fp_mem

    ! calculate sizes of large arrays that will go on the GPU.
      psib_mem = size(psib,kind=8)*8
      psi_mem = size(psi,kind=8)*8
      phi_mem = size(phi,kind=8)*8
      stime_mem = size(Geom%ZDataSoA%STime, kind=8)*8
      next_mem = size(QuadSet%next,kind=8)*4
      omega_a_fp_mem = size(Geom%ZDataSoA%omega_a_fp,kind=8)*8


      ! figure out if this problem fits on the GPU.
      !print *,"size(psib) =", psib_mem
      !print *,"size(psi) =", psi_mem
      !print *,"size(phi) =", phi_mem
      !print *,"size(stime) = ", stime_mem
      !print *,"size(next_mem) = ", next_mem
      !print *,"size(omega_A_fp) =", omega_a_fp_mem

      ! estimate for the total memory requirements of the sweep:
      sweep_mem = &
           psib_mem + &
           psi_mem + &
           phi_mem + &
           stime_mem + &
           next_mem + &
           2*omega_a_fp_mem

    
  end subroutine GPUmemRequirements


  subroutine InitDeviceBuffers()
    use Size_mod
    implicit none


    ! events:

    ! allocate space for event recordings (necessary to overlap gpu movement with kernel compute)
    allocate( Psi_OnDevice(Nbatches), Psi_OnHost(Nbatches) )
    allocate( Psib_OnDevice(Nbatches), Psib_OnHost(Nbatches) )
    allocate( SweepFinished(Nbatches), STimeFinished(Nbatches) )
    allocate( ExitFluxDFinished(Nbatches), AfpFinished(Nbatches) )
    allocate( snmomentsFinished(Nbatches) )


    ! device buffers:

    ! if the problem fits entirely in GPU memory, the buffers can be made the same size as the host arrays

    ! But if the problem does not fit in GPU memory, a double buffer strategy is used and the buffers have
    ! to be sized to fit chunks of the problem. (currently chunked by angle bins)

    batchsize = maxval(QuadSet%NangBinList(:))

    ! allocate the storage space available on the GPU
    call gpuStorage_ctor(psi_storage,numPsi_buffers,QuadSet%Groups,Size%ncornr,BATCHSIZE)
    call gpuStorage_ctor(STime_storage,numSTime_buffers,QuadSet%Groups,Size%ncornr,BATCHSIZE)
    call gpuStorage_ctor(psib_storage,numPsib_buffers,QuadSet%Groups,Size%nbelem,BATCHSIZE)

    call dotStorage_ctor(omega_A_fp_storage,numDot_buffers,Size% nzones,Size% maxCorner,Size% maxcf, BATCHSIZE)
    call dotStorage_ctor(omega_A_ez_storage,numDot_buffers,Size% nzones,Size% maxCorner,Size% maxcf, BATCHSIZE)

    
    ! allocate the other GPU arrays:
    allocate(d_phi(QuadSet%Groups,Size%ncornr))
   

  end subroutine InitDeviceBuffers


  subroutine gpuStorage_ctor(storage,numslots,size1,size2,size3)
    implicit none
    ! allocates storage and storage data for data that has 3 indices.
    type(gpuStorage), allocatable, intent(inout) :: storage(:) ! will be psi_storage, psib_s, or STime_s
    integer, intent(in) :: numslots, size1, size2, size3

    ! local
    integer :: slot
    integer :: numelements

    numelements = size1*size2*size3

    allocate( storage( numslots ) )
    
    do slot = 1, numslots
       allocate( storage(slot)%data(size1,size2,size3) )
       ! initiallly no one ownes this storage slot
       storage(slot)%owner = 0
       ! record the slot number
       storage(slot)%slot = slot
       ! handy to keep the size around
       !storage%numelements=numelements
    enddo



  end subroutine gpuStorage_ctor



  subroutine dotStorage_ctor(storage,numslots,size1,size2,size3,size4)
    implicit none
    ! allocates storage and storage data for data that has 4 indices.
    type(dotStorage), allocatable, intent(inout) :: storage(:) ! will be A_fp or ez
    integer, intent(in) :: numslots, size1, size2, size3, size4

    ! local
    integer :: slot
    integer :: numelements

    numelements = size1*size2*size3*size4

    allocate( storage( numslots ) )
    do slot = 1, numslots
       allocate( storage(slot)%data(size1,size2,size3,size4) )
       ! initiallly no one ownes this storage slot
       storage(slot)%owner = 0
       ! record the slot number
       storage(slot)%slot = slot
       ! handy to keep the size around
       !storage%numelements=numelements
    enddo

  end subroutine dotStorage_ctor



  subroutine CreateEvents()
    implicit none

    ! create an event for each batch.

    integer :: batch, istat

    do batch = 1, Nbatches
       istat = cudaEventCreate(Psi_OnDevice(batch))
       istat = cudaEventCreate(STimeFinished(batch))
       istat = cudaEventCreate(Psib_OnDevice(batch))
       istat = cudaEventCreate(Psib_OnHost(batch))
       istat = cudaEventCreate(SweepFinished(batch))
       istat = cudaEventCreate(Psi_OnHost(batch))
       istat = cudaEventCreate(ExitFluxDFinished(batch))
       istat = cudaEventCreate(AfpFinished(batch))
       istat = cudaEventCreate(snmomentsFinished(batch))

    enddo

    istat = cudaEventCreate(phi_OnHost)       

  end subroutine CreateEvents

  subroutine CreateEventsWithFlags()
    implicit none

    ! create an event for each batch.

    integer :: batch, istat

    do batch = 1, Nbatches
       istat = cudaEventCreateWithFlags(Psi_OnDevice(batch),cudaEventBlockingSync)
       istat = cudaEventCreateWithFlags(STimeFinished(batch),cudaEventBlockingSync)
       istat = cudaEventCreateWithFlags(Psib_OnDevice(batch),cudaEventBlockingSync)
       istat = cudaEventCreateWithFlags(Psib_OnHost(batch),cudaEventBlockingSync)
       istat = cudaEventCreateWithFlags(SweepFinished(batch),cudaEventBlockingSync)
       istat = cudaEventCreateWithFlags(Psi_OnHost(batch),cudaEventBlockingSync)
       istat = cudaEventCreateWithFlags(ExitFluxDFinished(batch),cudaEventBlockingSync)
       istat = cudaEventCreateWithFlags(AfpFinished(batch),cudaEventBlockingSync)
       istat = cudaEventCreateWithFlags(snmomentsFinished(batch),cudaEventBlockingSync)

    enddo

    istat = cudaEventCreateWithFlags(phi_OnHost, cudaEventBlockingSync)

    if(istat .ne. cudaSuccess) write( 0,*) "Error in cuda event creation"

  end subroutine CreateEventsWithFlags



!  subroutine checkDataOnDevice(p_storage, storage, h_data, bin, batch, prevslot, mm1, numelements, streamid, event)
  subroutine checkDataOnDevice(p_storage, storage, bin, prevslot)
    ! This routine checks if the bin is already being stored on the device. It returns a pointer to either where it 
    ! exists on the device, or where it should exist (i.e. that is where it should be moved). The move
    ! routine does the check to see if data really needs to be moved. If data is computed instead of moved,
    ! slot owner should be manually changed after the data is computed.
    implicit none
    
    !  Arguments

    type(gpuStorage), intent(inout), pointer :: p_storage ! pointer that will point into storage
    type(gpuStorage), intent(inout), target :: storage(:) ! storage container with slots for storing bins worth of data
    !real(adqt), intent(in) :: h_data(:,:,:) ! host buffer
    integer, intent(in) :: bin ! bin that is being swept
    !integer, intent(in) :: batch ! the batch number being processed
    integer, intent(in) :: prevslot ! previous slot used so I can use a different one (if possible)
    !integer, intent(in) :: mm1 ! starting angle index within a bin. (will be 1 when batches are sized the same as angle bins)
    !integer,intent(in) :: numelements ! number of array elements to be moved
    !integer(kind=cuda_stream_kind), intent(in) :: streamid
    !type(cudaEvent), intent(in) :: event(:)
    
    ! local variables
    integer :: istat

    integer :: slot, numslots

    !print *, "inside check data"

    ! get the number of slots in storage container
    numslots = size(storage) 
    !print *, "numslots = ", numslots

    ! check each slot to see if the bin is already stored on the device
    CheckBuffer: do slot = 1, numslots
       if(storage(slot)%owner == bin) then
          ! found a slot with the data on it
          ! point at this slot
          p_storage => storage(slot)
          !print *, "found bin ", bin, "in slot ", slot
          return 
       endif
    enddo CheckBuffer

    ! there was not a slot already assigned for this bin, so pick a slot:
    ! use a different slot than was used last time (because the data in that slot may still be in use)
    slot = 1+modulo(prevslot,numslots)
    !print *, "slot selected inside checkdata = ", slot

    ! point p_storage to the storage spot
    p_storage => storage(slot)


  end subroutine checkDataOnDevice
    
  subroutine checkDataOnDeviceDot(p_storage, storage, bin, prevslot)
    ! This routine checks if the bin is already being stored on the device. It returns a pointer to either where it 
    ! exists on the device, or where it should exist (i.e. that is where it should be moved). The move
    ! routine does the check to see if data really needs to be moved. If data is computed instead of moved,
    ! slot owner should be manually changed after the data is computed.
    implicit none
    
    !  Arguments

    type(dotStorage), intent(inout), pointer :: p_storage ! pointer that will point into storage
    type(dotStorage), intent(inout), target :: storage(:) ! storage container with slots for storing bins worth of data
    !real(adqt), intent(in) :: h_data(:,:,:) ! host buffer
    integer, intent(in) :: bin ! bin that is being swept
    !integer, intent(in) :: batch ! the batch number being processed
    integer, intent(in) :: prevslot ! previous slot used so I can use a different one (if possible)
    !integer, intent(in) :: mm1 ! starting angle index within a bin. (will be 1 when batches are sized the same as angle bins)
    !integer,intent(in) :: numelements ! number of array elements to be moved
    !integer(kind=cuda_stream_kind), intent(in) :: streamid
    !type(cudaEvent), intent(in) :: event(:)
    
    ! local variables
    integer :: istat

    integer :: slot, numslots

    ! get the number of slots in storage container
    numslots = size(storage) 

    ! check each slot to see if the bin is already stored on the device
    CheckBuffer: do slot = 1, numslots
       if(storage(slot)%owner == bin) then
          ! found a slot with the data on it
          ! point at this slot
          p_storage => storage(slot)
          return 
       endif
    enddo CheckBuffer

    ! there was not a slot already assigned for this bin, so pick a slot:
    ! use a different slot than was used last time (because the data in that slot may still be in use)
    slot = 1+modulo(prevslot,numslots)

    ! point p_storage to the storage spot
    p_storage => storage(slot)


  end subroutine checkDataOnDeviceDot
    



  subroutine MoveDataOnDevice(p_storage, storage, h_data, bin, batch, prevslot, mm1, numelements, streamid, event)

    implicit none
    
    !  Arguments

    type(gpuStorage), intent(inout), pointer :: p_storage ! pointer that will point into storage
    type(gpuStorage), intent(inout), target :: storage(:) ! storage container with slots for storing bins worth of data
    real(adqt), intent(in) :: h_data(:,:,:) ! host buffer
    integer, intent(in) :: bin ! bin that is being swept
    integer, intent(in) :: batch ! the batch number being processed
    integer, intent(in) :: prevslot ! previous slot used so I can use a different one (if possible)
    integer, intent(in) :: mm1 ! starting angle index within a bin. (will be 1 when batches are sized the same as angle bins)
    integer,intent(in) :: numelements ! number of array elements to be moved
    integer(kind=cuda_stream_kind), intent(in) :: streamid
    type(cudaEvent), intent(in) :: event(:)
    ! local variables
    integer :: istat

    !integer :: slot, numslots

    ! if the data is already on the device do not actually do the move
    ! I searched all storage slots in previous routine--if bin was already on the device,
    ! p_storage would be pointing at that bin and this would evaluate to false.
    if(p_storage%owner /= bin) then

       ! move the data from host into this slot
       istat=cudaMemcpyAsync(p_storage%data(1,1,1),                 &
            h_data(1,1,QuadSet%AngleOrder(mm1,bin)), &
            numelements, streamid )

       ! change ownership: mark this bins data as residing in this storage slot.
       p_storage%owner = bin

    endif

    ! Record when movement event finishes (for example, psi on device) for this batch
    istat=cudaEventRecord( event(batch), streamid )

  end subroutine MoveDataOnDevice
    

  subroutine MoveDataOnDeviceDot(p_storage, storage, h_data, bin, batch, prevslot, mm1, numelements, streamid, event)

    implicit none
    
    !  Arguments

    type(dotStorage), intent(inout), pointer :: p_storage ! pointer that will point into storage
    type(dotStorage), intent(inout), target :: storage(:) ! storage container with slots for storing bins worth of data
    real(adqt), intent(in) :: h_data(:,:,:,:) ! host buffer
    integer, intent(in) :: bin ! bin that is being swept
    integer, intent(in) :: batch ! the batch number being processed
    integer, intent(in) :: prevslot ! previous slot used so I can use a different one (if possible)
    integer, intent(in) :: mm1 ! starting angle index within a bin. (will be 1 when batches are sized the same as angle bins)
    integer,intent(in) :: numelements ! number of array elements to be moved
    integer(kind=cuda_stream_kind), intent(in) :: streamid
    type(cudaEvent), intent(in) :: event(:)
    ! local variables
    integer :: istat

    !integer :: slot, numslots

    ! if the data is already on the device do not actually do the move
    ! I searched all storage slots in previous routine--if bin was already on the device,
    ! p_storage would be pointing at that bin and this would evaluate to false.
    if(p_storage%owner /= bin) then

       ! move the data from host into this slot
       istat=cudaMemcpyAsync(p_storage%data(1,1,1,1),                 &
            h_data(1,1,1,QuadSet%AngleOrder(mm1,bin)), &
            numelements, streamid )

       ! change ownership: mark this bins data as residing in this storage slot.
       p_storage%owner = bin

    endif

    ! Record when movement event finishes (for example, psi on device) for this batch
    istat=cudaEventRecord( event(batch), streamid )

  end subroutine MoveDataOnDeviceDot




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


  ! subroutine stageGPUData(buffer,batch,mm1)
  !   use Size_mod
  !   implicit none

  !   integer, intent(in) :: buffer, batch(:), mm1

  !   ! local variables

  !   integer istat
  !   integer thisbatch

  !   ! select that batch used in this routine.
  !   thisbatch = batch(buffer)

  !   ! buffer can be either current or next. Just a way of putting the data movement staging that happens 
  !   ! before and after the sweep into one reusable function. Buffer determines if you are staging it for the
  !   ! current buffer or the next buffer.



  !   ! ! need to move fp stuff to device
  !   ! istat = cudaMemcpyAsync(d_omega_A_fp(1,1,1,1), &
  !   !          Geom%ZDataSoA%omega_A_fp(1,1,1,QuadSet%AngleOrder(mm1,binSend(buffer))), &
  !   !          Size% nzones*Size% maxCorner*Size% maxcf*anglebatch(buffer), transfer_stream)

  !   ! ! need to move ez stuff to device
  !   ! istat = cudaMemcpyAsync(d_omega_A_ez(1,1,1,1), &
  !   !          Geom%ZDataSoA%omega_A_ez(1,1,1,QuadSet%AngleOrder(mm1,binSend(buffer))), &
  !   !          Size% nzones*Size% maxCorner*Size% maxcf*anglebatch(buffer), transfer_stream)

  !   ! ! fp and ez stuff will be finished if STimeFinished, since it is run in same stream below.
  !   ! ! This means fp and ez are ready on the GPU if STime is ready.



  !   if( FitsOnGPU ) then 
  !      ! If data fits on GPU do not worry about data movement. Just record that STime is ready:
  !      !istat=cudaEventRecord(STimeFinished(batch), transfer_stream )
  !   else !data does not fit on GPU and has to be streamed

  !      ! If this is first temp and intensity iteration, STime would have been calculated, needs update to host:
  !      if (calcSTime == .true.) then

  !      else ! STime already computed,
  !         ! just need to move section of STime to device (Never called when FitsOnGPU)
  !         call checkDataOnDevice(d_STime, Geom%ZDataSoA%STime, batch, buffer, mm1, &
  !              QuadSet%Groups*Size%ncornr*anglebatch(buffer), transfer_stream, &
  !              STimeFinished)
                   
  !      endif

  !   endif


  !   ! could select whether to recalculate fp each time. For now always do it.

  !   ! call fp_ez_c(     anglebatch(buffer),                     &
  !   !      Size%nzones,               &
  !   !      QuadSet%Groups,            &
  !   !      Size%ncornr,               &
  !   !      QuadSet%NumAngles,         &
  !   !      QuadSet%d_AngleOrder(mm1,binSend(buffer)),        & ! only need angle batch portion
  !   !      Size%maxCorner,            &
  !   !      Size%maxcf,                &
  !   !      Nangbin(buffer),                   &
  !   !      Size%nbelem,                &
  !   !      QuadSet%d_omega,             &
  !   !      Geom%ZDataSoA%nCorner,                &
  !   !      Geom%ZDataSoA%nCFaces,                &
  !   !      Geom%ZDataSoA%c0,                &
  !   !      Geom%ZDataSoA%A_fp,                &
  !   !      Geom%ZDataSoA%omega_A_fp,                &
  !   !      Geom%ZDataSoA%A_ez,                &
  !   !      Geom%ZDataSoA%omega_A_ez,                &
  !   !      Geom%ZDataSoA%Connect,             &
  !   !      Geom%ZDataSoA%Connect_reorder,             &
  !   !      QuadSet%d_next,              &
  !   !      QuadSet%d_nextZ,             &
  !   !      QuadSet%d_passZstart,        &
  !   !      kernel_stream           &
  !   !      )

  ! end subroutine stageGPUData




   subroutine scalePsibyVolume(psir, ZData, anglebatch, streamid)
     ! scaling by change in mesh volume that used to be done in advanceRT is done here
     use kind_mod
     use constant_mod
     use Quadrature_mod
     use ZoneData_mod
     use Size_mod
     use cudafor
     
     implicit none
     
     !  Arguments

     real(adqt), device, intent(inout)  :: psir(QuadSet%Groups,Size%ncornr,anglebatch) 
     type(GPU_ZoneData), device, intent(in) :: ZData(Size%nzones)
     !real(adqt), device, intent(in) :: volumeRatio(Size%ncornr)
     integer, intent(in)  :: anglebatch 
     integer(kind=cuda_stream_kind), intent(in) :: streamid

     !  Local

     integer    :: ia, zone,c, ig, nzones,maxCorner,nCorner, Groups
     integer    :: c0, ic
     real(adqt) :: tau


     nzones = Size%nzones
     !ncornr = Size%ncornr
     Groups = QuadSet% Groups   

     maxCorner = Size% maxCorner

     !$cuf kernel do(4) <<< *, *, stream=streamid >>>
     do ia=1,anglebatch
        do zone=1, nzones
           do c=1,maxCorner 
              do ig=1, Groups
                 nCorner = ZData(zone)%nCorner
                 if(c<nCorner) then
                    c0 = ZData(zone)%c0
                    ic = c0+c
                    psir(ig,ic,ia) = psir(ig,ic,ia)*ZData(zone)%volumeRatio(c)
                 endif
              enddo
          enddo
       enddo
     enddo


     return
   end subroutine scalePsibyVolume




   subroutine computeSTime(psiccache, STimeBatch, anglebatch, streamid)
     ! Multiply by tau to get STime
     use kind_mod
     use constant_mod
     use Quadrature_mod
     use Size_mod
     use cudafor
     
     implicit none
     
     !  Arguments

     real(adqt), device, intent(in)  :: psiccache(QuadSet%Groups,Size%ncornr,anglebatch) 
     real(adqt), device, intent(out) :: STimeBatch(QuadSet%Groups,Size%ncornr,anglebatch)
     integer, intent(in)  :: anglebatch 
     integer(kind=cuda_stream_kind), intent(in) :: streamid

     !  Local

     integer    :: ia, ic, ig, ncornr, Groups
     real(adqt) :: tau


     ncornr = Size%ncornr
     Groups = QuadSet% Groups   
     tau    = Size%tau


     !$cuf kernel do(3) <<< *, *, stream=streamid >>>
     do ia=1,anglebatch
        do ic=1,ncornr
           do ig=1, Groups
              STimeBatch(ig,ic,ia) = tau*psiccache(ig,ic,ia)
           enddo
        enddo
     enddo


     return
   end subroutine computeSTime



  attributes(global)   subroutine computeSTimeD(psiccache, STimeBatch, anglebatch, groups, ncornr, tau)
     ! Multiply by tau to get STime
     !use kind_mod
     !use constant_mod
     !use Quadrature_mod
     !use Size_mod
     !use cudafor
     
     implicit none
     
     !  Arguments

     real(adqt), device, intent(in)  :: psiccache(Groups,ncornr,anglebatch) 
     real(adqt), device, intent(out) :: STimeBatch(Groups,ncornr,anglebatch)
     integer, value, intent(in)  :: anglebatch, Groups, ncornr
     real(adqt), value, intent(in) :: tau

     !  Local

     integer    :: ia, ic, ig

     do ia=1,anglebatch !blockIdx.x
        do ic=1,ncornr !threadIdx.y
           do ig=1, Groups !threadIdx.x
              STimeBatch(ig,ic,ia) = tau*psiccache(ig,ic,ia)
           enddo
        enddo
     enddo


     return
   end subroutine computeSTimeD



end module GPUhelper_mod
