!***********************************************************************
!                        Version 1:  09/96, PFN                        *
!                                                                      *
!   SNFLWXYZ - This routine, called by RSWPMD and RTACCELMD, solves    *
!              the fixed-source transport problem on an arbitrary      *
!              grid in either xyz-geometry or rz-geometry.             *
!              An upstream corner-balance spatial discretization is    *
!              used.                                                   *
!                                                                      *
!   Input:                                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!***********************************************************************

! if batches are not currently size of bins, data is staged wrong.
!!!#define BATCHSIZE 32
   subroutine snflwxyz(ipath, PSIB, PSI, PHI, angleLoopTime, intensityIter, tempIter)

   use, intrinsic :: iso_c_binding
   use kind_mod
   use constant_mod
   use Size_mod
   use Quadrature_mod
   use Geometry_mod
   use ZoneData_mod
   use snswp3d_mod
   use cudafor
   use nvtx_mod
   use GPUhelper_mod

#include "assert.h"
!  Assertion checking include file for TETON

   implicit none
   include 'mpif.h'

   ! ! Fortran to C interface
   ! interface 

   !  subroutine fp_ez_c ( &
   !        anglebatch, &
   !        numzones, &
   !        numgroups, &
   !        ncornr, &
   !        numAngles, &
   !        AngleOrder, &
   !        maxcorners, &
   !        maxfaces, &
   !        nangbin, &
   !        nbelem, &
   !        omega, &
   !        numCorners, &
   !        numCFaces, &
   !        c0, &
   !        A_fp , &
   !        omega_A_fp , &
   !        A_ez , &
   !        omega_A_ez , &
   !        Connect , &
   !        Connect_reorder, &
   !        next, &
   !        nextZ, &   
   !        passZ, &
   !        streamid &
   !        ) &
   !        bind ( c ) 

   !       use iso_c_binding
   !       use cudafor
   !       integer ( c_int ) :: anglebatch
   !       integer ( c_int ) :: numzones
   !       integer ( c_int ) :: numgroups
   !       integer ( c_int ) :: ncornr
   !       integer ( c_int ) :: numAngles
   !       integer ( c_int ), device :: AngleOrder(*)
   !       integer ( c_int ) :: maxcorners
   !       integer ( c_int ) :: maxfaces
   !       integer ( c_int ) :: nangbin
   !       integer ( c_int ) :: nbelem
   !       real ( c_double ),device :: omega(*)
   !       integer ( c_int ),device :: numCorners(*) 
   !       integer ( c_int ),device :: numCFaces(*) 
   !       integer ( c_int ),device :: c0(*)
   !       real ( c_double ),device :: A_fp(*) 
   !       real ( c_double ),device :: omega_A_fp(*) 
   !       real ( c_double ),device :: A_ez(*) 
   !       real ( c_double ),device :: omega_A_ez(*) 
   !       integer ( c_int ),device :: Connect(*) 
   !       integer ( c_int ),device :: Connect_reorder(*) 
   !       integer ( c_int ),device :: next(*)
   !       integer ( c_int ),device :: nextZ(*)
   !       integer ( c_int ),device :: passZ(*)
   !       integer ( kind=cuda_stream_kind ),value :: streamid
   !     end subroutine fp_ez_c

   !     subroutine snswp3d_c ( &
   !        anglebatch, &
   !        numzones, &
   !        numgroups, &
   !        ncornr, &
   !        numAngles, &
   !        AngleOrder, &
   !        maxcorners, &
   !        maxfaces, &
   !        binRecv, &
   !        nangbin, &
   !        nbelem, &
   !        omega, &
   !        numCorners, &
   !        numCFaces, &
   !        c0, &
   !        A_fp , &
   !        omega_A_fp , &
   !        A_ez , &
   !        omega_A_ez , &
   !        Connect , &
   !        Connect_reorder, &
   !        STotal , &
   !        STimeBatch , &
   !        Volume , &
   !        psic, &
   !        psib, &
   !        next, &
   !        nextZ, &   
   !        Sigt, &
   !        SigtInv, &
   !        passZ, &
   !        calcSTime, &
   !        tau, &
   !        streamid &
   !        ) &
   !        bind ( c ) 

   !       use iso_c_binding
   !       use cudafor
   !       integer ( c_int ) :: anglebatch
   !       integer ( c_int ) :: numzones
   !       integer ( c_int ) :: numgroups
   !       integer ( c_int ) :: ncornr
   !       integer ( c_int ) :: numAngles
   !       integer ( c_int ), device :: AngleOrder(*)
   !       integer ( c_int ) :: maxcorners
   !       integer ( c_int ) :: maxfaces
   !       integer ( c_int ) :: binRecv
   !       integer ( c_int ) :: nangbin
   !       integer ( c_int ) :: nbelem
   !       real ( c_double ),device :: omega(*)
   !       integer ( c_int ),device :: numCorners(*) 
   !       integer ( c_int ),device :: numCFaces(*) 
   !       integer ( c_int ),device :: c0(*)
   !       real ( c_double ),device :: A_fp(*) 
   !       real ( c_double ),device :: omega_A_fp(*) 
   !       real ( c_double ),device :: A_ez(*) 
   !       real ( c_double ),device :: omega_A_ez(*) 
   !       integer ( c_int ),device :: Connect(*) 
   !       integer ( c_int ),device :: Connect_reorder(*) 
   !       real ( c_double ),device :: STotal(*) 
   !       real ( c_double ),device :: STimeBatch(*)
   !       real ( c_double ),device :: Volume(*) 
   !       real ( c_double ),device :: psic(*) 
   !       real ( c_double ),device :: psib(*) 
   !       integer ( c_int ),device :: next(*)
   !       integer ( c_int ),device :: nextZ(*)
   !       real ( c_double ),device :: Sigt(*) 
   !       real ( c_double ),device :: SigtInv(*) 
   !       integer ( c_int ),device :: passZ(*)
   !       logical ( c_bool ) :: calcSTime
   !       real ( c_double ) :: tau
   !       integer ( kind=cuda_stream_kind ),value :: streamid
   !     end subroutine snswp3d_c


   !  end interface

!  Arguments

   real(adqt), intent(inout) :: psib(QuadSet%Groups,Size%nbelem,QuadSet%NumAngles)
   real(adqt), intent(inout) :: psi(QuadSet%Groups,Size%ncornr,QuadSet%NumAngles)
   real(adqt), intent(inout) :: Phi(QuadSet%Groups,Size%ncornr),angleLoopTime

   character(len=8), intent(in) :: ipath

   integer, intent(in) :: intensityIter, tempIter ! current flux and temperature iteration from rtmainsn

!  Local
   
   logical(kind=1), save :: first_time = .true.

   integer          :: Angle, mm,mm1,mm2,  binRecv 
   integer          :: Groups, fluxIter, ishared


   logical (kind=1) :: FluxConverged

   real(adqt)       :: maxFluxError
   real(adqt)       :: startOMPLoopTime, endOMPLoopTime, theOMPLoopTime

   ! Cuda streams and double buffer managment stuff
   integer :: s, batch, istat, current, next   
   
   integer :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
   integer NumAngles, nbelem, ncornr, NumBin, myrank, info

   !integer :: devnum, cacheconfig

!  Convenient Mesh Constants

   Groups = QuadSet%Groups
   ! NumAngles = QuadSet%NumAngles
   nbelem = Size%nbelem
   ncornr = Size%ncornr
   ! print *, ncornr
   ! NangBin = maxval(QuadSet%NangBinList(:))
   ! NumBin = QuadSet%NumBin
   !call mpi_comm_rank(mpi_comm_world, myrank, info)


   
   ! This sets up to allow zero copy use of phi directly on the device:
   ! Get a device pointer for phi, put it to d_phi_p
   !istat = cudaHostGetDevicePointer(d_phi_p, C_LOC(phi(1,1)), 0)
   ! Translate that C pointer to the fortran array with given dimensions
   !call c_f_pointer(d_phi_p, d_phi, [QuadSet%Groups,Size%ncornr] )
   
   ! if(1) then
   !    ! This sets up to allow zero copy use of STime directly on the device:
   !    ! Get a device pointer for STime, put it to d_STime_p
   !    istat = cudaHostGetDevicePointer(d_STime_p, C_LOC(Geom%ZDataSoA%STime(1,1,1)), 0)
   !    ! Translate that C pointer to the fortran array with given dimensions
   !    call c_f_pointer(d_STime_p, d_STime, [QuadSet%Groups,Size%ncornr, Size%nangSN] )
   ! endif



   theOMPLoopTime=0.0

   !call mpi_comm_rank(mpi_comm_world, myrank, info)
   
   
   ! ! Set the Cache configuration for the GPU (use more L1, less shared)
   ! istat = cudaDeviceGetCacheConfig(cacheconfig)
   ! print *, "cacheconfig =", cacheconfig
   ! if (cacheconfig .eq. cudaFuncCachePreferShared) then
   !    print *, "L1 set for shared memory usage"
   ! elseif (cacheconfig .eq. cudaFuncCachePreferL1) then
   !    print *, "L1 set for hardware caching (prefer L1)."
   ! else
   !    print *, "other L1 configuration present."
   ! endif

   ! cacheconfig = cudaFuncCachePreferL1
   !istat = cudaDeviceSetCacheConfig(cacheconfig)
   !istat = cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte)


   if (first_time) then
      ! Create streams that can overlap 
      istat = cudaStreamCreate(transfer_stream)
      istat = cudaStreamCreate(kernel_stream)

      call InitDeviceBuffers()

      first_time = .false.
   endif

   call nvtxStartRange("createEvents")
   ! Create events to synchronize among different streams
   call CreateEvents()
   call nvtxEndRange

!  Loop over angle bins

   if (ipath == 'sweep') then
     call timer_beg('_setflux')
     call nvtxStartRange("InitExchange")
     call setIncidentFlux(psib)
     call nvtxEndRange
     call timer_end('_setflux')
   endif
                                                                                         
   FluxConverged = .FALSE.
   fluxIter      =  0

   call restoreCommOrder(QuadSet)


   FluxIteration: do

!    Post receives for all data
                                                                                                  
     
     call timer_beg('_initexch')
     call nvtxStartRange("InitExchange")
     call InitExchange
     call nvtxEndRange
     call timer_end('_initexch')

     fluxIter = fluxIter + 1

     if (ipath == 'sweep') then
       !phi(:,:) = zero
        d_phi=0 ! could try something async memset here.
     endif

!    Loop over angles, solving for each in turn:
     startOMPLoopTime = MPI_WTIME()
     
     ! By default STime does not need to be computed.
     calcSTime = .false.

     ! If this is first temp and intensity iteration, need to calculate STime 
     if (intensityIter == 1 .and. tempIter == 1 .and. fluxIter==1) then
        ! compute STime from initial d_psi
        calcSTime = .true.
     else
        ! STime already computed,
        calcSTime = .false.
     endif


     call timer_beg('_anglebins')     
     AngleBin: do binRecv=1,QuadSet% NumBin


        ! loop over batches of angles within the angle bin. 
        ! BatchLoop: do batch=1,Nbatches

        batch = binRecv
        !batch(next) = binRecv+1

        ! cycle the buffers used to hold the batches.

        ! there will be either 2 buffers if the problem does not fit in GPU, or NangBin if it fits.
        current = 1 + mod(batch,numGPUbuffers) ! gives 2,1,2,1,2,1,2
        next = 1+mod(current,numGPUbuffers)      ! gives 1,2,1,2, etc.


        binSend(current) = QuadSet% SendOrder(binRecv)
        binSend(next) = QuadSet% SendOrder(binRecv+1) ! binSend on next iteration
        NangBin(current) = QuadSet% NangBinList(binSend(current))
        NangBin(next) = QuadSet% NangBinList(binSend(next))


        !   lower angle index for this batch
        mm1=1
        !   get the angles upper index for this batch
        mm2=NangBin(current)
        !   ! True size of the batch (not always BATCHSIZE for last batch if not nicely divisible)
        anglebatch(current)=mm2-mm1+1
        anglebatch(next)=NangBin(next)
        
        FirstOctant: if (binRecv == 1) then
           !for other bins, will begin staging in the data at the end of prev
           !iteration of the loop


           call MoveHtoD(d_psi, psi, current, mm1, &
                QuadSet%Groups*Size%ncornr*anglebatch(current), transfer_stream, Psi_OnDevice(batch))

           ! If this is first temp and intensity iteration, need to calculate STime
           if (calcSTime == .true.) then
              ! have kernel stream wait until transfer of psi to device
              istat = cudaStreamWaitEvent(kernel_stream, Psi_OnDevice(batch), 0)
              ! compute STime from initial d_psi
              call computeSTime(d_psi(1,1,1,current), d_STimeBatch(1,1,1,current), anglebatch(current), kernel_stream )

              istat=cudaEventRecord(STimeFinished(batch), kernel_stream )
           endif

        endif FirstOctant

        ! relfected angles is done on CPU, while above is taking place on GPU
        call nvtxStartRange("snreflect")
        ! Set angular fluxes for reflected angles
        do mm=mm1,mm2
           call snreflect(QuadSet%AngleOrder(mm,binSend(current)), PSIB)
        enddo
        call nvtxEndRange

        !Stage batch of psib into GPU (after reflected angles is completed on host)
        ! happens regardless of problem size.
        call MoveHtoD(d_psibBatch, psib, current, mm1, &
             QuadSet%Groups*Size%nbelem*anglebatch(current), transfer_stream, Psib_OnDevice(batch))

        FirstOctant2: if (binRecv == 1) then

           call stageGPUData(current,batch,mm1)

        endif FirstOctant2


        ! Do not launch sweep kernel until psib is on GPU transfer in transfer stream is done.
        istat = cudaStreamWaitEvent(kernel_stream, Psib_OnDevice(batch), 0)

        ! would put snreflect here if done on GPU.

        ! Also make sure STime is ready before launching sweep
        istat = cudaStreamWaitEvent(kernel_stream, STimeFinished(batch), 0)


        !        Sweep the mesh, calculating PSI for each corner; the
        !        boundary flux array PSIB is also updated here.
        !        Mesh cycles are fixed automatically.
        if(0) then ! call CUDA fortran version
           ! call snswp3d(     anglebatch,                     &
           !      Size%nzones,               &
           !      QuadSet%Groups,            &
           !      Size%ncornr,               &
           !      QuadSet%NumAngles,         &
           !      QuadSet%d_AngleOrder(mm1,binSend),        & ! only need angle batch portion
           !      Size%maxCorner,            &
           !      Size%maxcf,                &
           !      NangBin,                   &
           !      Size%nbelem,                &
           !      QuadSet%d_omega,             &
           !      Geom%ZDataSoA%nCorner,                &
           !      Geom%ZDataSoA%nCFaces,                &
           !      Geom%ZDataSoA%c0,                &
           !      Geom%ZDataSoA%A_fp,                &
           !      Geom%ZDataSoA%A_ez,                &
           !      Geom%ZDataSoA%Connect,             &
           !      Geom%ZDataSoA%STotal,              &
           !                      !Geom%ZDataSoA%STime,               &
           !      d_STimeBatch, &
           !      Geom%ZDataSoA%Volume,             &
           !      d_psi,                      &  ! only want angle batch portion
           !      d_psib,                      &
           !      QuadSet%d_next,              &
           !      QuadSet%d_nextZ,             &
           !      Geom%ZDataSoA%Sigt,                &
           !      Geom%ZDataSoA%SigtInv,             &
           !      QuadSet%d_passZstart,              &
           !      kernel_stream)
        else ! Call CUDA c version
           call snswp3d_c(     anglebatch(current),                     &
                Size%nzones,               &
                QuadSet%Groups,            &
                Size%ncornr,               &
                QuadSet%NumAngles,         &
                QuadSet%d_AngleOrder(mm1,binSend(current)),        & ! only need angle batch portion
                Size%maxCorner,            &
                Size%maxcf,                &
                binRecv,                   &
                NangBin(current),                   &
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
                Geom%ZDataSoA%STotal,              &
                                !Geom%ZDataSoA%STime,               &
                d_STimeBatch(1,1,1,current),          &
                Geom%ZDataSoA%Volume,             &
                d_psi(1,1,1,current),                      &  ! only want angle batch portion
                d_psibBatch(1,1,1,current),                      &
                QuadSet%d_next,              &
                QuadSet%d_nextZ,             &
                Geom%ZDataSoA%Sigt,                &
                Geom%ZDataSoA%SigtInv,             &
                QuadSet%d_passZstart,        &
                calcSTime,                  &
                Size%tau,             &
                kernel_stream           &
                )
        endif

        ! record when sweep is finished for this batch
        istat=cudaEventRecord(SweepFinished(batch), kernel_stream )



        NotLastOctants: if (binRecv < QuadSet% NumBin) then
           ! If not the last bin (octant), pre-stage data for next angle bin 
           call MoveHtoD(d_psi, psi, next, mm1, &
                QuadSet%Groups*Size%ncornr*anglebatch(next), transfer_stream, Psi_OnDevice(batch+1))

        endif NotLastOctants
        
        
        if (ipath == 'sweep') then
           call timer_beg('__snmoments')
           ! snmoments only reads d_psi, produces d_phi
           call snmomentsD(d_psi(1,1,1,current), d_phi, QuadSet%d_Weight,     &
                QuadSet%d_AngleOrder(mm1,binSend(current)),      &
                anglebatch(current), kernel_stream) ! GPU version, one batch at a time
           call timer_end('__snmoments')
        endif


        if ( .not. fitsOnGPU ) then
           ! before moving DtoH psi, sweep needs to complete
           istat=cudaStreamWaitEvent(transfer_stream, SweepFinished(batch), 0)

           ! Copy d_psi to host psi.
           istat=cudaMemcpyAsync(psi(1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
                d_psi(1,1,1,current), &
                QuadSet%Groups*Size%ncornr*anglebatch(current), transfer_stream )

           istat=cudaEventRecord(Psi_OnHost(batch), transfer_stream)

        endif


        NotLastOctants2: if (binRecv < QuadSet% NumBin) then
           ! If not the last bin (octant), pre-stage data for next angle bin 

           ! If this is first temp and intensity iteration, need to calculate STime
           if (calcSTime == .true.) then
              ! have kernel stream wait until transfer of psi to device
              istat = cudaStreamWaitEvent(kernel_stream, Psi_OnDevice(batch+1), 0)
              ! compute STime from initial d_psi
              call computeSTime(d_psi(1,1,1,next), d_STimeBatch(1,1,1,next), anglebatch(next), kernel_stream )

              istat=cudaEventRecord(STimeFinished(batch+1), kernel_stream )
           endif

           call stageGPUData(next,batch+1,mm1)
       
        endif NotLastOctants2



        !call timer_beg('__setExitFlux')

        ! Set the exit flux
        if( fitsOnGPU ) then
           call setExitFluxD<<<batchsize,Groups,0,kernel_stream>>>(anglebatch(current), &
                QuadSet%d_AngleOrder(mm1,binSend(current)),  &
                d_psi(1,1,1,current), d_psibBatch(1,1,1,current),&
                QuadSet%d_iExit, groups, ncornr, nbelem )

           ! When setExitFluxD is done, record it has finished in kernel stream
           istat=cudaEventRecord(ExitFluxD(batch), kernel_stream )

           ! transfer stream should wait for event setExitFluxD to finish
           istat = cudaStreamWaitEvent(transfer_stream, ExitFluxD(batch), 0)

           ! need to move psib to Host (or later exchange from GPU).
           istat=cudaMemcpyAsync(psib(1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
                d_psibBatch(1,1,1,current), &
                QuadSet%Groups*Size%nbelem*anglebatch(current), transfer_stream ) 

           ! CPU code should wait until psib is on the host.
           istat=cudaEventSynchronize( psib_OnHost(batch) )

        else ! Problem does not fit in GPU, and psi is moved back from GPU anyway

           ! if ready then set exit flux and move exchange psib.
           istat=cudaEventSynchronize( Psi_OnHost(batch) )
           ! compute exit flux on host. 
           call nvtxStartRange("setExitFlux")
           call setExitFlux(anglebatch(current), QuadSet%AngleOrder(mm1,binSend(current)), psi, psib)
           call nvtxEndRange

        endif

        !call timer_end('__setExitFlux')


        !      Exchange Boundary Fluxes
        ! these need to become non-blocking

        call timer_beg('__exch')
        call nvtxStartRange("exchange")
        call exchange(PSIB, binSend(current), binRecv) 
        call nvtxEndRange
        call timer_end('__exch')

     enddo AngleBin

     call timer_end('_anglebins')

     endOMPLoopTime = MPI_WTIME()
     theOMPLoopTime = theOMPLoopTime + (endOMPLoopTime-startOMPLoopTime)

     istat = cudaDeviceSynchronize()


     if (ipath == 'sweep') then
        call timer_beg('_setflux')
        call setIncidentFlux(psib)
        call timer_end('_setflux')
        call testFluxConv(FluxConverged, fluxIter, maxFluxError)
     else
        FluxConverged = .TRUE.
     endif

     if ( FluxConverged ) then
        exit FluxIteration
     else
        call setCommOrder(QuadSet)
        cycle FluxIteration
     endif

  enddo FluxIteration

  !  Update the scaler flux

  if (ipath == 'sweep') then
     ! move d_phi data to host:
     istat=cudaMemcpyAsync(phi(1,1), &
                   d_phi(1,1), &
                   QuadSet%Groups*Size%ncornr, transfer_stream )

     call restoreCommOrder(QuadSet)
  endif

  angleLoopTime = angleLoopTime + theOMPLoopTime


  return
end subroutine snflwxyz


