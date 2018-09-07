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

   integer          :: Angle, mm,mm1,mm2,  binRecv, slot
   integer          :: Groups, fluxIter, ishared


   logical (kind=1) :: FluxConverged

   real(adqt)       :: maxFluxError
   real(adqt)       :: startOMPLoopTime, endOMPLoopTime, theOMPLoopTime

   
   integer :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
   integer NumAngles, nbelem, ncornr, NumBin, myrank, info

   integer :: t0, t1, timelimit
   integer :: shared_mem

   timelimit = 30 ! assume hung if event is waiting more than timelimit seconds

!  Convenient Mesh Constants

   Groups = QuadSet%Groups
   ! NumAngles = QuadSet%NumAngles
   nbelem = Size%nbelem
   ncornr = Size%ncornr
   ! print *, ncornr
   ! NangBin = maxval(QuadSet%NangBinList(:))
   NumBin = QuadSet%NumBin
   
   call mpi_comm_rank(mpi_comm_world, myrank, info)

   istat = cudaDeviceSynchronize() ! strangely seems to be needed to prevent hangs.
   
   ! attempt to set the shared mem/cache config for GPU_sweep:

   shared_mem = 38*1024 ! 96 Kb
   istat = cudaFuncSetAttribute(GPU_sweep,cudaFuncAttributeMaxDynamicSharedMemorySize, shared_mem)
   if(myrank .eq. 0) then
      if(istat .ne. 0) print *, cudaGetErrorString(istat)
      print *, "set the func attribute with error = ", cudaGetErrorString(istat)
   endif


   !istat = cudafuncsetcacheconfig(GPU_sweep,cudaFuncCachePreferShared)
   !if(myrank .eq. 0) then
   !   if(istat .ne. 0) print *, cudaGetErrorString(istat)
   !endif


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

   
   ! sets up zero copy of psib needed for snrefelctD on device.
   istat = cudaHostGetDevicePointer(d_psib_p, C_LOC(psib(1,1,1)), 0)
   ! Translate that C pointer to the fortran array with given dimensions
   call c_f_pointer(d_psib_p, pinned_psib, [QuadSet%Groups, Size%nbelem, QuadSet%NumAngles] )


   theOMPLoopTime=0.0   

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

     if (intensityIter == 1 .and. tempIter == 1 .and. fluxIter==1) then
        ! compute STime from initial d_psi
        calcSTime = .true.


     else
        ! STime already computed,
        calcSTime = .false.

     endif


!     if(myrank == 0 ) print *, "NumBins = ", QuadSet% NumBin
!     if(myrank == 0 ) print *, "binSend(1) = ", QuadSet% SendOrder(1)
!     if(myrank == 0 ) print *, "binSend(NumBin) = ", QuadSet% SendOrder(QuadSet%NumBin)

!     if(myrank == 0 ) print *, "--------------"
!     if(myrank == 0 ) print *, "SendOrder(:) = ", QuadSet% SendOrder(1:8)
!     if(myrank == 0 ) print *, "--------------"


     !call timer_beg('_anglebins')     
     AngleBin: do binRecv=1,QuadSet% NumBin

        call timer_beg('_sweepoct')

        ! which bin is the current bin being worked on?
        current%bin       = QuadSet%SendOrder(binRecv) ! same as old binSend
        current%batch     = binRecv    ! batch 1,2,3... sequentual ordering
        current%NangBin   = QuadSet% NangBinList( current%bin )
        current%anglebatch =current%NangBin ! later can be subset of angles in a bin.
        !current%Angles=QuadSet%AngleOrder(mm1,current%bin)

        !print *, "--------------------------"

        !print *, "current%bin       =", current%bin
        !print *, "current%batch     =", current%batch

        ! debugging optimized code:                                                                                           

        !print *, "psib before starting sweeps: ", psib(1,1,QuadSet%AngleOrder(mm1,current%bin )) 
   
        !print *, "psi before starting sweeps: ", psi(1,1,QuadSet%AngleOrder(mm1,current%bin )) 

        !print *, "STime before starting sweeps: ", Geom%ZData(1)%STime(1,1,QuadSet%AngleOrder(mm1,current%bin ))


        !   lower angle index for this batch (later can help with batches that are not all the angles in the bin)
        mm1=1


        ! Stage batch of psib into GPU 
        ! always needs to be updated, so set the owned slot to 0 so it will be updated:
        do slot=1,numPsib_buffers
           psib_storage(slot)%owner=0
        enddo


        ! setup pointers to psib storage
        call checkDataOnDevice(current%psib, psib_storage, current%bin, previous%psib%slot)

        ! move psib onto GPU
        call moveDataOnDevice(current%psib, psib_storage, psib, current%bin, current%batch, & 
             previous%psib%slot, mm1, QuadSet%Groups*Size%nbelem*current%anglebatch, transfer_stream, &
             Psib_OnDevice)

        
        ! ! If this is first temp and intensity iteration, can calculate STime while above is copying in
        ! if (calcSTime == .true.) then
        !    ! have kernel stream wait until transfer of psi to device (calc depends on psi)
        !    istat = cudaStreamWaitEvent(kernel_stream, Psi_OnDevice(current%batch), 0)

        !    ! scalePsiby Volume and compute STime is also done in advanceRT
        !    ! so is already done for as many psi slots as fit on the GPU
        !    if( current%batch > numPsi_buffers ) then
        !       ! scale current batch of psi
        !       call scalePsibyVolume(current%psi%data(1,1,1), Geom%ZDataSoA%volumeRatio, current%anglebatch, kernel_stream )  
              
        !    endif

        !    if( current%batch > numSTime_buffers ) then

        !       ! set up STime pointers before computing STime
        !       call checkDataOnDevice(current%STime, STime_storage, current%bin, previous%STime%slot)

        !       ! compute current batch STime from current batch d_psi (only possible first sweep of the timestep)
        !       call computeSTime(current%psi%data(1,1,1), current%STime%data(1,1,1), current%anglebatch, kernel_stream )

        !       ! need to record that current batch of STime is held in this storage slot (so it does not need to be copied in)
        !       current%STime% owner = current%bin

        !    endif

        !    istat=cudaEventRecord(STimeFinished( current%batch ), kernel_stream )
        ! endif


        ! Do not launch snreflect kernel until psib is on GPU.
        istat = cudaStreamWaitEvent(kernel_stream, Psib_OnDevice( current%batch ), 0)


        ! Set angular fluxes for reflected angles
        ! relfected angles could be done on CPU or GPU, both steal CPU bandwidth needed from copies.
        call snreflectD(current%anglebatch, QuadSet%d_AngleOrder(mm1,current%bin), &
             current%psib%data(1,1,1), pinned_psib, kernel_stream)


        ! make sure omega_A stuff is on the GPU: (often will not be called)

        ! setup pointers to omega_A stuff.
        call checkDataOnDeviceDot(current%omega_A_fp,omega_A_fp_storage, current%bin, previous%omega_A_fp%slot)
        call checkDataOnDeviceDot(current%omega_A_ez,omega_A_ez_storage, current%bin, previous%omega_A_ez%slot)


        ! stage omega_A_fp into GPU
        call moveDataOnDeviceDot(current%omega_A_fp, omega_A_fp_storage, &
             Geom%ZDataSoA%omega_A_fp, &
             current%bin, current%batch, previous%omega_A_fp%slot, &
             mm1, Size% nzones*Size% maxCorner*Size% maxcf*current%anglebatch, transfer_stream, &
             AfpFinished )


        ! stage omega_A_ez into GPU
        call moveDataOnDeviceDot(current%omega_A_ez, omega_A_ez_storage, &
             Geom%ZDataSoA%omega_A_ez, &
             current%bin, current%batch, previous%omega_A_ez%slot, &
             mm1, Size% nzones*Size% maxCorner*Size% maxcf*current%anglebatch, transfer_stream, &
             AfpFinished )



        ! this is not often called, but is needed because sometimes the order of bins change and bin that was staged in at
        ! the end of the last fluxiter was not the bin we actually start with.
        ! advanceRT will set up at least the first bin for every timestep, so this can take place after computeSTime.
        if(1) then

           ! check if bin is already held in a slot, and get a pointer to a good slot.
           call checkDataOnDevice(current%psi, psi_storage, current%bin, previous%psi%slot)

           ! move the data (only moves if it was not found on the device)
           call moveDataOnDevice(current%psi, psi_storage, psi, current%bin, current%batch, &
                previous%psi%slot,mm1, QuadSet%Groups*Size%ncornr*current%anglebatch, transfer_stream, &
                Psi_OnDevice )


           call checkDataOnDevice(current%STime, STime_storage, current%bin, previous%STime%slot)
           !if( d_STime(current)%owner /= current%batch ) print *, "I was called, batch = ", current%batch
           call moveDataOnDevice(current%STime, STime_storage, Geom%ZDataSoA%STime, current%bin, current%batch, &
                previous%STime%slot, mm1, QuadSet%Groups*Size%ncornr*current%anglebatch, transfer_stream, &
                STimeFinished)


        endif

        !Stime_temp = current%STime%data(1,33,1)
        !volumeRatio_temp = Geom%ZDataSoA%volumeRatio(33)
        !print *, "STime(1,33,1) = ", Stime_temp
        !print *, "volumeRatio = ", volumeRatio_temp

        ! YOU SHOULD BE ABLE TO MOVE STIME OFF THE DEVICE HERE IF CALCSTIME=TRUE
        ! THIS WILL OVERLAP WITH SWEEP KERNEL

        ! Also make sure STime is ready before launching sweep
        istat = cudaStreamWaitEvent(kernel_stream, STimeFinished( current%batch ), 0)


        !        Sweep the mesh, calculating PSI for each corner; the
        !        boundary flux array PSIB is also updated here.
        !        Mesh cycles are fixed automatically.
        ! call CUDA fortran sweep version
        call snswp3d_f(     current%anglebatch,                     &
                Size%nzones,               &
                QuadSet%Groups,            &
                Size%ncornr,               &
                QuadSet%NumAngles,         &
                QuadSet%d_AngleOrder(mm1,current%bin),        & ! only need a bin of angles at a time
                Size%maxCorner,            &
                Size%maxcf,                &
                !binRecv,                   &
                current%NangBin,                   &
                Size%nbelem,                &
                current%omega_A_fp%data,                &
                current%omega_A_ez%data,                &
                                !Geom%ZDataSoA%STime,               &
                current%STime%data(1,1,1),          &
                current%psi%data(1,1,1),                      &  ! only want angle batch portion
                current%psib%data(1,1,1),                      &
                QuadSet%d_next,              &
                QuadSet%d_nextZ,             &
                QuadSet%d_passZstart,        &
                !calcSTime,                  &
                !Size%tau,             &
                kernel_stream           &
                )


        ! record when sweep is finished for this batch
        istat=cudaEventRecord(SweepFinished( current%batch ), kernel_stream )


        !print *, "Sweep finished", fluxiter

        !!!!! Start of things that will overlap sweep kernel !!!!!

        ! DO ME:::::
        ! take off previous batch psi,
        ! put on next batch psi
        ! this would allow single buffer of STime, saving memory.

        if ( .not. fitsOnGPU ) then

           ! at the end of the flux iter, current should become previous. Once you have used
           ! previous, you can use it for next, and then once you use next you can use it for
           ! the new previous again.

           ! previous should be known coming into this region.

           

           ! before moving DtoH psi, previous sweep needs to complete
           istat=cudaStreamWaitEvent(transfer_stream, SweepFinished( previous%batch ), 0)

           ! Copy d_psi to host psi.
           istat=cudaMemcpyAsync(psi(1,1,QuadSet%AngleOrder(mm1,previous%bin )), &
                   previous%psi%data(1,1,1), &
                   QuadSet%Groups*Size%ncornr*previous%anglebatch, transfer_stream )

           istat=cudaEventRecord( Psi_OnHost( previous%batch ), transfer_stream)
           
           !if(myrank == 0 ) print *, "binRecv        = ", binRecv
           !if(myrank == 0 ) print *, "binSend(previous) = ", QuadSet% SendOrder(previous_batch)
           !if(myrank == 0 ) print *, "binSend(current)  = ", QuadSet% SendOrder(current%batch)
           !if(myrank == 0 ) print *, "binSend(next)     = ", QuadSet% SendOrder(batch(next))
           
           
           ! which bin is the next bin that will be worked on? (or if last bin, best guess)
           ! numbin is used instead of numBatches because actually don't always do all 8 batches.
           next%batch     = 1+modulo(binRecv,QuadSet% NumBin)    ! batch 1,2,3... sequentual ordering
           next%bin       = QuadSet%SendOrder(next%batch) ! same as old binSend
           next%NangBin   = QuadSet% NangBinList( next%bin )
           next%anglebatch = next%NangBin ! later can be subset of angles in a bin.

           !  pre-stage next batch of psi into the device

           call checkDataOnDevice(next%psi, psi_storage, next%bin, current%psi%slot)

           call MoveDataOnDevice(next%psi, psi_storage, psi, next%bin, next%batch, &
                current%psi%slot, mm1, QuadSet%Groups*Size%ncornr*next%anglebatch, transfer_stream, &
                Psi_OnDevice )

        endif


        !print *, "bin: (prev, current, next) = ", previous%bin, current%bin, next%bin



        !!!!! End of things that will overlap sweep kernel  !!!!


        call setExitFluxD2(current%anglebatch, &
             QuadSet%d_AngleOrder(mm1,current%bin),  &
             current%psi%data(1,1,1), current%psib%data(1,1,1),&
             QuadSet%m_iExit, groups, ncornr, nbelem, kernel_stream )

        ! call setExitFluxD<<<batchsize,Groups,0,kernel_stream>>>(anglecurrent%batch, &
        !      QuadSet%d_AngleOrder(mm1,binSend(current)),  &
        !      d_psi(current)%data(1,1,1), d_psibBatch(1,1,1,current),&
        !      QuadSet%d_iExit, groups, ncornr, nbelem)



        ! When setExitFluxD is done, record it has finished in kernel stream
        istat=cudaEventRecord(ExitFluxDFinished( current%batch ), kernel_stream )

        if (istat /= 0) then
           write(0,*) "CUDA event record API error:",istat
           stop
        endif
        
        ! transfer stream should wait for event setExitFluxD to finish
        istat = cudaStreamWaitEvent(transfer_stream1, ExitFluxDFinished( current%batch ), 0)

        if (istat /= 0) then
           write(0,*) "CUDA StreamWaitEven API error:",istat
           stop
        endif
        
        ! need to move psib to Host (or later exchange from GPU).
        istat=cudaMemcpyAsync(psib(1,1,QuadSet%AngleOrder(mm1,current%bin)), &
             current%psib%data(1,1,1), &
             QuadSet%Groups*Size%nbelem*current%anglebatch, transfer_stream1 ) 
        
        istat=cudaEventRecord(psib_OnHost( current%batch ), transfer_stream1 )


        !!!!! Start of things that will overlap exchange !!!!
        
        ! THIS NEXT BATCH OF STIME SHOULD ACTUALLY GO RIGHT INTO THE CURRENT BATCH STIME 
        ! AS LONG AS SWEEP IS FINISHED. (NO DOUBLE BUFFER OF STIME)

        ! don't move STime into device until sweep is finished (only 1 buffer)
        istat=cudaStreamWaitEvent(transfer_stream, SweepFinished( current%batch ), 0)

        if ( calcSTime ) then
           
           ! below is commented out because advance RT saves the copy of STime now.
           ! I think this should be done before the sweep for better overlap too.

           ! ! transfer stream waits for STime to be computed:
           ! istat = cudaStreamWaitEvent(transfer_stream, STimeFinished(current%batch), 0)
           ! ! move current batch STime to host
           ! istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
           !      d_STime(current)%data(1,1,1), &
           !      QuadSet%Groups*Size%ncornr*anglecurrent%batch, transfer_stream ) 

           if ( binRecv == QuadSet% NumBin ) then

              ! move next batch STime anyway since it is for batch 1 of next flux iteration:
              call checkDataOnDevice(next%STime, STime_storage, next%bin, current%STime%slot)

              ! check/move next batch of STime onto GPU. This should be done even for last bin to prepare for next iteration
              call moveDataOnDevice(next%STime, STime_storage, Geom%ZDataSoA%STime, next%bin, next%batch, &
                   current%STime%slot, mm1, QuadSet%Groups*Size%ncornr*next%anglebatch, transfer_stream, &
                   STimeFinished)

              
           endif


        else


           call checkDataOnDevice(next%STime, STime_storage, next%bin, current%STime%slot)

           ! check/move next batch of STime onto GPU. This should be done even for last bin to prepare for next iteration
           call moveDataOnDevice(next%STime, STime_storage, Geom%ZDataSoA%STime, next%bin, next%batch, &
                current%STime%slot, mm1, QuadSet%Groups*Size%ncornr*next%anglebatch, transfer_stream, &
                STimeFinished)

        endif

        
        if (ipath == 'sweep') then
           call timer_beg('__snmoments')
           ! snmoments only reads d_psi, produces d_phi
           call snmomentsD(current%psi%data(1,1,1), d_phi, QuadSet%d_Weight,     &
                QuadSet%d_AngleOrder(mm1,current%bin),      &
                current%anglebatch, kernel_stream) ! GPU version, one batch at a time

           istat=cudaEventRecord(snmomentsFinished( current%batch ), kernel_stream )


           ! if this is the last bin, just move phi to the host 
           ! (only need to move when converged, but it is small and this allows better overlap)
           if (binRecv == QuadSet% NumBin ) then
              ! transfer stream1 waits for snmoments calc to be finished (the last current one)
              istat = cudaStreamWaitEvent(transfer_stream1, snmomentsFinished(current%batch), 0)


              ! move d_phi data to host:
              istat=cudaMemcpyAsync(phi(1,1), &
                   d_phi(1,1), &
                   QuadSet%Groups*Size%ncornr, transfer_stream1 )

              istat=cudaEventRecord( phi_OnHost, transfer_stream1 )
           endif


           call timer_end('__snmoments')
        endif

        ! before moving DtoH omega_A, previous sweep needs to complete
        istat=cudaStreamWaitEvent(transfer_stream, SweepFinished( current%batch ), 0)

        ! setup pointers to omega_A stuff.
        call checkDataOnDeviceDot(next%omega_A_fp,omega_A_fp_storage, next%bin, current%omega_A_fp%slot)
        call checkDataOnDeviceDot(next%omega_A_ez,omega_A_ez_storage, next%bin, current%omega_A_ez%slot)

        ! stage omega_A_fp into GPU
        call moveDataOnDeviceDot(next%omega_A_fp, omega_A_fp_storage, &
             Geom%ZDataSoA%omega_A_fp, &
             next%bin, next%batch, current%omega_A_fp%slot, &
             mm1, Size% nzones*Size% maxCorner*Size% maxcf*next%anglebatch, transfer_stream, &
             AfpFinished )

        ! stage omega_A_ez into GPU
        call moveDataOnDeviceDot(next%omega_A_ez, omega_A_ez_storage, &
             Geom%ZDataSoA%omega_A_ez, &
             next%bin, next%batch, current%omega_A_ez%slot, &
             mm1, Size% nzones*Size% maxCorner*Size% maxcf*next%anglebatch, transfer_stream, &
             AfpFinished )




        !!!! End of things that will overlap exchange
        
        ! I should put a cuda get last error here.

        if(.false.) then ! old way
           ! CPU code should wait until psib is on the host before exchanging.
           istat=cudaEventSynchronize( psib_OnHost( current%batch ) )
        else 
           t0 = MPI_WTIME()
           t1 = MPI_WTIME()
           istat=cudaErrorNotReady
           do while (istat == cudaErrorNotReady)
              ! new debug way:
              istat=cudaEventQuery( psib_OnHost( current%batch ) )
              t1 = MPI_WTIME()
              if( t1-t0 > timelimit ) then
                 ! this mpi rank is hung so end the profiler
                 !call cudaProfilerStop()
                 exit ! break out of while loop to see if later gpu call will show crash
              endif
           enddo
           if (istat /= cudaSuccess) print *, "rank: ", myrank, " Event query never succeeded, istat = ", cudaGetErrorString(istat)
        endif



        !      Exchange Boundary Fluxes
        ! these need to become non-blocking

        call timer_beg('__exch')
        call nvtxStartRange("exchange")
        call exchange(PSIB, current%bin, binRecv) 
        call nvtxEndRange
        call timer_end('__exch')



        ! end of loop, so set previous to the current.
        previous = current

        call timer_end('_sweepoct')

     enddo AngleBin

     !call timer_end('_anglebins')


     endOMPLoopTime = MPI_WTIME()
     theOMPLoopTime = theOMPLoopTime + (endOMPLoopTime-startOMPLoopTime)

     if (ipath == 'sweep') then
        call timer_beg('_setflux')
        call setIncidentFlux(psib)
        call timer_end('_setflux')
        call testFluxConv(FluxConverged, fluxIter, maxFluxError)
     else
        FluxConverged = .TRUE.
     endif

     if ( FluxConverged ) then
        if(myrank == 0) print *, "True flux iterations = ", fluxIter
        exit FluxIteration
     else
        call setCommOrder(QuadSet)
        cycle FluxIteration
     endif

  enddo FluxIteration

  !  Update the scaler flux

  if (ipath == 'sweep') then
     
     call restoreCommOrder(QuadSet)

     ! CPU code should wait until phi is on the host before using it
     istat=cudaEventSynchronize( phi_OnHost )

  endif


  angleLoopTime = angleLoopTime + theOMPLoopTime


  return
end subroutine snflwxyz


