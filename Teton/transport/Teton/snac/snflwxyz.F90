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
   use snreflect

#include "assert.h"
!  Assertion checking include file for TETON

   implicit none
   include 'mpif.h'

!  Arguments

   real(adqt), intent(inout) :: psib(QuadSet%Groups,Size%nbelem,QuadSet%NumAngles)
   real(adqt), intent(inout) :: psi(QuadSet%Groups,Size%ncornr,QuadSet%NumAngles)
   real(adqt), intent(inout) :: Phi(QuadSet%Groups,Size%ncornr),angleLoopTime

   character(len=8), intent(in) :: ipath

   integer, intent(in) :: intensityIter, tempIter ! current flux and temperature iteration from rtmainsn

!  Local

   integer          :: Angle, mm,mm1,mm2,  binRecv 
   integer          :: Groups, fluxIter, ishared


   logical (kind=1) :: FluxConverged

   real(adqt)       :: maxFluxError
   real(adqt)       :: startOMPLoopTime, endOMPLoopTime, theOMPLoopTime

   
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
   NumBin = QuadSet%NumBin
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

   
   ! sets of zero copy of psib needed for snrefelctD on device.
   istat = cudaHostGetDevicePointer(d_psib_p, C_LOC(psib(1,1,1)), 0)
   ! Translate that C pointer to the fortran array with given dimensions
   call c_f_pointer(d_psib_p, pinned_psib, [QuadSet%Groups, Size%nbelem, QuadSet%NumAngles] )


   theOMPLoopTime=0.0

   !call mpi_comm_rank(mpi_comm_world, myrank, info)
   
   

   ! if (first_time) then
   !    ! Create streams that can overlap 
   !    istat = cudaStreamCreate(transfer_stream)
   !    istat = cudaStreamCreate(kernel_stream)

   !    call InitDeviceBuffers()

   !    first_time = .false.
   ! endif

   ! NOT NEEDED?
   call nvtxStartRange("createEvents")
   ! Create events to synchronize among different streams
   call CreateEvents()
   call nvtxEndRange

!  Loop over angle bins

   if (ipath == 'sweep') then
     call timer_beg('_setflux')
     call nvtxStartRange("setIncidentFlux")
     call setIncidentFlux(psib)
     call nvtxEndRange
     call timer_end('_setflux')
   endif
                                                                                         
   FluxConverged = .FALSE.
   fluxIter      =  0

   call restoreCommOrder(QuadSet)


   FluxIteration: do
     
     !  Post receives for all data                                                                                               
     
     call timer_beg('_initexch')
     call nvtxStartRange("InitExchange")
     call InitExchange
     call nvtxEndRange
     call timer_end('_initexch')

     fluxIter = fluxIter + 1

     if (ipath == 'sweep') then
        !d_phi=0 ! could try something async memset here.
        istat = cudaMemset(d_phi,zero,QuadSet%Groups*Size%ncornr)
     endif

!    Loop over angles, solving for each in turn:
     startOMPLoopTime = MPI_WTIME()
     
     ! By default STime does not need to be computed.
     calcSTime = .false.

     ! If this is first temp and intensity iteration, several things should be done that were previously done in other routines:
     ! 1. need to calculate STime (was done in rtstrsn.F90). Now done for each storage buffer as they are loaded in.
     ! 2. need to scale phi and psi by the ratio of the volume change.
     if (intensityIter == 1 .and. tempIter == 1 .and. fluxIter==1) then
        ! compute STime from initial d_psi
        calcSTime = .true.
        scaleVolume = .true.

        ! ! Mark all the psi buffers as stale and needing update (trying to find bug)
        ! do buffer=1,8
        !    d_psi(buffer)% owner = 0
        ! enddo

     else
        ! STime already computed,
        calcSTime = .false.
        scaleVolume = .false.
     endif


     call timer_beg('_anglebins')     
     AngleBin: do binRecv=1,QuadSet% NumBin

        ! cycle the buffers used to hold the batches.
        ! there will be either 2 buffers if the problem does not fit in GPU, or NumBin if it fits.
        current = 1 + mod(binRecv-1,numGPUbuffers) ! gives 1,2,1,2...  or 1,2,3...7,8
        next = 1+mod(binRecv,numGPUbuffers)      ! gives 2,1,2,1... or 2,3,4...8,1 

        ! each batch corresponds to a angle bin.
        batch(current) = binRecv
        batch(next) = 1 + mod(binRecv,QuadSet% NumBin)

        
        binSend(current) = QuadSet% SendOrder(binRecv)
        binSend(next) = QuadSet% SendOrder(batch(next)) ! binSend on next iteration

        
        NangBin(current) = QuadSet% NangBinList(binSend(current))
        NangBin(next) = QuadSet% NangBinList(binSend(next))


        !   lower angle index for this batch
        mm1=1
        !   get the angles upper index for this batch
        mm2=NangBin(current)
        !   ! True size of the batch (not always BATCHSIZE for last batch if not nicely divisible)
        anglebatch(current)=mm2-mm1+1
        anglebatch(next)=NangBin(next)
        
        !print *, "current = ", current
        !print *, "anglebatch(current) = ", anglebatch(current)
        !print *, "binSend(current) = ", binSend(current)

        ! Stage batch of psib into GPU 
        ! happens regardless of problem size.
        istat = cudaMemcpyAsync(d_psibBatch(1,1,1,current), psib(1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
             QuadSet%Groups*Size%nbelem*anglebatch(current), transfer_stream)

        istat=cudaEventRecord(Psib_OnDevice( batch(current) ), transfer_stream )



        ! If this is first temp and intensity iteration, can calculate STime while above is copying in
        if (calcSTime == .true.) then
           ! have kernel stream wait until transfer of psi to device (calc depends on psi)
           istat = cudaStreamWaitEvent(kernel_stream, Psi_OnDevice(batch(current)), 0)

           ! compute current batch STime from current batch d_psi (only possible first sweep of the timestep)
           call computeSTime(d_psi(current)%data(1,1,1), d_STime(current)%data(1,1,1), anglebatch(current), kernel_stream )

           ! need to record that current batch of STime is held in this device buffer (so it does not need to be copied in)
           d_STime(current)% owner = batch(current)

           istat=cudaEventRecord(STimeFinished( batch(current) ), kernel_stream )
        endif


        ! Do not launch snreflect kernel until psib is on GPU.
        istat = cudaStreamWaitEvent(kernel_stream, Psib_OnDevice( batch(current) ), 0)



        ! Set angular fluxes for reflected angles
        ! relfected angles could be done on CPU or GPU, both steal CPU bandwidth needed from copies.
        call snreflectD(anglebatch(current), QuadSet%d_AngleOrder(mm1,binSend(current)), &
             d_psibBatch(1,1,1,current), pinned_psib, kernel_stream)



        ! Expect this will not need to be called, because
        ! advanceRT will set up at least the first bin.
        if(1) then
           call checkDataOnDevice(d_psi, psi, batch, current, mm1, &
                QuadSet%Groups*Size%ncornr*anglebatch(current), transfer_stream, &
                Psi_OnDevice )

           if( d_STime(current)%owner /= batch(current) ) print *, "I was called, batch = ", batch(current)
           ! I think this is not needed either (always a no-op)
           call checkDataOnDevice(d_STime, Geom%ZDataSoA%STime, batch, current, mm1, &
                QuadSet%Groups*Size%ncornr*anglebatch(current), transfer_stream, &
                STimeFinished)


        endif

        ! YOU SHOULD BE ABLE TO MOVE STIME OFF THE DEVICE HERE IF CALCSTIME=TRUE
        ! THIS WILL OVERLAP WITH SWEEP KERNEL

        ! Also make sure STime is ready before launching sweep
        istat = cudaStreamWaitEvent(kernel_stream, STimeFinished( batch(current) ), 0)


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
                d_omega_A_fp,                &
                Geom%ZDataSoA%A_ez,                &
                d_omega_A_ez,                &
                Geom%ZDataSoA%Connect,             &
                Geom%ZDataSoA%Connect_reorder,             &
                Geom%ZDataSoA%STotal,              &
                                !Geom%ZDataSoA%STime,               &
                d_STime(current)%data(1,1,1),          &
                Geom%ZDataSoA%Volume,             &
                d_psi(current)%data(1,1,1),                      &  ! only want angle batch portion
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
        istat=cudaEventRecord(SweepFinished( batch(current) ), kernel_stream )

        !!!!! Start of things that will overlap sweep kernel !!!!!

        ! DO ME:::::
        ! take off previous batch psi,
        ! put on next batch psi
        ! this would allow single buffer of STime, saving memory.

        if ( .not. fitsOnGPU ) then
!           if ( binRecv /= 1 ) then
              ! take the previous batch of psi off the device (have previous, current, and next)
              previous = 1+ modulo(binRecv-2,numGPUbuffers) ! gives 2,1,2,1...  or 8,1,2...6,7

              previous_batch = 1 + modulo(binRecv-2,QuadSet% NumBin)

              previous_binSend = QuadSet% SendOrder(previous_batch) ! binSend on previous iteration

              ! before moving DtoH psi, previous sweep needs to complete
              istat=cudaStreamWaitEvent(transfer_stream, SweepFinished(previous_batch), 0)

              !print *, "previous = ", previous
              !print *, "previous_batch = ", previous_batch
              !print *, "anglebatch(previous) = ", anglebatch(previous)

              ! Copy d_psi to host psi.
              istat=cudaMemcpyAsync(psi(1,1,QuadSet%AngleOrder(mm1,previous_binSend )), &
                   d_psi(previous)%data(1,1,1), &
                   QuadSet%Groups*Size%ncornr*anglebatch(previous), transfer_stream )

              istat=cudaEventRecord( Psi_OnHost( previous_batch ), transfer_stream)



              !  pre-stage next batch of psi into the device
              call checkDataOnDevice(d_psi, psi, batch, next, mm1, &
                   QuadSet%Groups*Size%ncornr*anglebatch(next), transfer_stream, &
                   Psi_OnDevice)

!           endif


        endif

        !!!!! End of things that will overlap sweep kernel  !!!!


        call setExitFluxD2(anglebatch(current), &
             QuadSet%d_AngleOrder(mm1,binSend(current)),  &
             d_psi(current)%data(1,1,1), d_psibBatch(1,1,1,current),&
             QuadSet%d_iExit, groups, ncornr, nbelem, kernel_stream )

        ! call setExitFluxD<<<batchsize,Groups,0,kernel_stream>>>(anglebatch(current), &
        !      QuadSet%d_AngleOrder(mm1,binSend(current)),  &
        !      d_psi(current)%data(1,1,1), d_psibBatch(1,1,1,current),&
        !      QuadSet%d_iExit, groups, ncornr, nbelem)


        ! When setExitFluxD is done, record it has finished in kernel stream
        istat=cudaEventRecord(ExitFluxDFinished( batch(current) ), kernel_stream )
        
        ! transfer stream should wait for event setExitFluxD to finish
        istat = cudaStreamWaitEvent(transfer_stream, ExitFluxDFinished( batch(current) ), 0)
        
        ! need to move psib to Host (or later exchange from GPU).
        istat=cudaMemcpyAsync(psib(1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
             d_psibBatch(1,1,1,current), &
             QuadSet%Groups*Size%nbelem*anglebatch(current), transfer_stream ) 
        
        istat=cudaEventRecord(psib_OnHost( batch(current) ), transfer_stream )


        !!!!! Start of things that will overlap exchange !!!!
        
        ! THIS NEXT BATCH OF STIME SHOULD ACTUALLY GO RIGHT INTO THE CURRENT BATCH STIME 
        ! AS LONG AS SWEEP IS FINISHED. (NO DOUBLE BUFFER OF STIME)

        if ( calcSTime ) then

           ! ! transfer stream waits for STime to be computed:
           ! istat = cudaStreamWaitEvent(transfer_stream, STimeFinished(batch(current)), 0)
           ! ! move current batch STime to host
           ! istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
           !      d_STime(current)%data(1,1,1), &
           !      QuadSet%Groups*Size%ncornr*anglebatch(current), transfer_stream ) 


        else

           ! check/move next batch of STime onto GPU. This should be done even for last bin to prepare for next iteration
           call checkDataOnDevice(d_STime, Geom%ZDataSoA%STime, batch, next, mm1, &
                QuadSet%Groups*Size%ncornr*anglebatch(next), transfer_stream, &
                STimeFinished)

        endif
        
        if (ipath == 'sweep') then
           call timer_beg('__snmoments')
           ! snmoments only reads d_psi, produces d_phi
           call snmomentsD(d_psi(current)%data(1,1,1), d_phi, QuadSet%d_Weight,     &
                QuadSet%d_AngleOrder(mm1,binSend(current)),      &
                anglebatch(current), kernel_stream) ! GPU version, one batch at a time

           istat=cudaEventRecord(snmomentsFinished( batch(current) ), kernel_stream )

           call timer_end('__snmoments')
        endif


        ! stage omega_A_fp into GPU
        istat = cudaMemcpyAsync( d_omega_A_fp(1,1,1,1), &
             Geom%ZDataSoA%omega_A_fp(1,1,1,QuadSet%AngleOrder(mm1,binSend(next))), &
             Size% nzones*Size% maxCorner*Size% maxcf*anglebatch(next), transfer_stream)


        ! stage omega_A_ez into GPU
        istat = cudaMemcpyAsync( d_omega_A_ez(1,1,1,1), &
             Geom%ZDataSoA%omega_A_ez(1,1,1,QuadSet%AngleOrder(mm1,binSend(next))), &
             Size% nzones*Size% maxCorner*Size% maxcf*anglebatch(next), transfer_stream)


        !!!! End of things that will overlap exchange

        ! CPU code should wait until psib is on the host before exchanging.
        istat=cudaEventSynchronize( psib_OnHost( batch(current) ) )

        !      Exchange Boundary Fluxes
        ! these need to become non-blocking

        call timer_beg('__exch')
        call nvtxStartRange("exchange")
        call exchange(PSIB, binSend(current), binRecv) 
        call nvtxEndRange
        call timer_end('__exch')

     enddo AngleBin

     call timer_end('_anglebins')


     ! ELIMINATE BY DOING ABOVE SWAP FOR ALL BINS, AND KEEPING BIN OWNER IN ADVANCERT?
     ! if( .not. fitsOnGPU ) then
     !    ! The last bin needs to be moved to the host as well, since inside the loop 
     !    ! before moving DtoH psi, current sweep needs to complete
     !    istat=cudaStreamWaitEvent(transfer_stream, SweepFinished(batch(current)), 0)

     !    ! Copy d_psi to host psi.
     !    istat=cudaMemcpyAsync(psi(1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
     !         d_psi(current)%data(1,1,1), &
     !         QuadSet%Groups*Size%ncornr*anglebatch(current), transfer_stream )

     !    istat=cudaEventRecord( Psi_OnHost( batch(current) ), transfer_stream)
     ! endif


     endOMPLoopTime = MPI_WTIME()
     theOMPLoopTime = theOMPLoopTime + (endOMPLoopTime-startOMPLoopTime)

     ! needed? I think I can get rid of this.
     !istat = cudaDeviceSynchronize()


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

     
     ! transfer stream1 waits for snmoments calc to be finished (the last current one)
     istat = cudaStreamWaitEvent(transfer_stream1, snmomentsFinished(batch(current)), 0)

     

     ! move d_phi data to host:
     istat=cudaMemcpyAsync(phi(1,1), &
                   d_phi(1,1), &
                   QuadSet%Groups*Size%ncornr, transfer_stream1 )
     
     istat=cudaEventRecord( phi_OnHost, transfer_stream1 )
     
     ! CPU code should wait until phi is on the host before using it
     istat=cudaEventSynchronize( phi_OnHost )


     ! May need device sync here if time between sweeps decreases.
     ! istat = cudaDeviceSynchronize()

     call restoreCommOrder(QuadSet)
  endif



  ! There are still some routines (advanceRT) that expect a host psi.
  ! so for now, once converged, move psi back from the device
  ! if( fitsOnGPU ) then
  !    ! Copy d_psi to host psi.
  !    do buffer=1, QuadSet% NumBin0 
  !       binSend(buffer) = QuadSet% SendOrder0(buffer)
  !       !print *, "QuadSet% NumBin = ", QuadSet% NumBin
  !       !print *, "binSend(buffer) = ", binSend(buffer)
  !       !print *, "mm1 = ", mm1
  !       !print *, "buffer = ", buffer
  !       !print *, "anglebatch(buffer) = ", anglebatch(buffer)
  !       istat=cudaMemcpyAsync(psi(1,1,QuadSet%AngleOrder(mm1,binSend(buffer))), &
  !            d_psi(buffer)%data(1,1,1), &
  !            QuadSet%Groups*Size%ncornr*batchsize, 0 )

  !       ! mark the data as un-owned since host will change it, making device version stale:
  !       d_psi(buffer)% owner = 0
  !       ! CHECKME: STime may be marked as stale more often than necessary.
  !       !d_STime(buffer)% owner = 0

  !    enddo

  ! endif



  angleLoopTime = angleLoopTime + theOMPLoopTime


  return
end subroutine snflwxyz


