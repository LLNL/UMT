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
#define BATCHSIZE 16
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

#include "assert.h"
!  Assertion checking include file for TETON

   implicit none
   include 'mpif.h'

   ! Fortran to C interface
   interface 

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
          STime , &
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
         real ( c_double ),device :: STime(*) 
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

!  Arguments

   real(adqt), intent(inout) :: psib(QuadSet%Groups,Size%nbelem,QuadSet%NumAngles)
   real(adqt), intent(inout) :: psi(QuadSet%Groups,Size%ncornr,QuadSet%NumAngles)
   real(adqt), intent(inout) :: Phi(QuadSet%Groups,Size%ncornr),angleLoopTime

   character(len=8), intent(in) :: ipath

   integer, intent(in) :: intensityIter, tempIter ! current flux and temperature iteration from rtmainsn

!  Local

   ! Cuda streams overlapping stuff
   integer, parameter :: Nbatches=2
   integer :: nStreams = 2*Nbatches*8 ! 2*Nbatches streams per octant (8 max, should be ok with less)
   integer(kind=cuda_stream_kind) :: stream(nStreams)
   type(cudaEvent) :: HtoDdone(nStreams),SweepFinished(nStreams),PsiOnHost(nStreams)
   type(cudaEvent) :: PsibOnDevice(nStreams), PsibOnHost(nStreams)
   integer :: s, batch, istat
   


   integer          :: Angle, mm,mm1,mm2,anglebatch, anglebatch_next
   integer          :: Groups, fluxIter, ishared
   integer          :: binSend, binSend_next, binRecv, NangBin, NangBin_next

   logical (kind=1) :: FluxConverged

   real(adqt)       :: maxFluxError
   real(adqt)       :: startOMPLoopTime, endOMPLoopTime, theOMPLoopTime

   ! zero copy pointers for phi and psib
   type(C_DEVPTR)                    :: d_phi_p
   type(C_DEVPTR)                    :: d_psib_p
   type(C_DEVPTR)                    :: d_STime_p
   real(adqt), device, allocatable :: d_phi(:,:)
   real(adqt), device, allocatable :: d_psib(:,:,:)
   real(adqt), device, allocatable :: d_STime(:,:,:)
   
   logical(kind=1) :: calcSTime

   ! picking up environment variables
   character(len=255) :: envstring

   real(adqt), device :: d_psi(QuadSet%Groups,Size%ncornr,BATCHSIZE,Nbatches)
   real(adqt), device :: d_STimeBatch(QuadSet%Groups,Size%ncornr,BATCHSIZE,Nbatches)
   real(adqt), device :: d_psibBatch(QuadSet%Groups,Size%nbelem,BATCHSIZE,Nbatches)
   
   type(C_PTR) :: cptr
   type(C_DEVPTR) :: dptr

   integer :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
   integer NumAngles, nbelem, ncornr, NumBin, myrank, info

   integer :: devnum, cacheconfig

!  Convenient Mesh Constants

   Groups = QuadSet%Groups

   NumAngles = QuadSet%NumAngles
   nbelem = Size%nbelem
   ncornr = Size%ncornr
   print *, ncornr
   NangBin = maxval(QuadSet%NangBinList(:))
   NumBin = QuadSet%NumBin
   call mpi_comm_rank(mpi_comm_world, myrank, info)



   ! This sets up to allow zero copy use of phi directly on the device:
   ! Get a device pointer for phi, put it to d_phi_p
   istat = cudaHostGetDevicePointer(d_phi_p, C_LOC(phi(1,1)), 0)
   ! Translate that C pointer to the fortran array with given dimensions
   call c_f_pointer(d_phi_p, d_phi, [QuadSet%Groups,Size%ncornr] )
   
   if(1) then
      ! This sets up to allow zero copy use of STime directly on the device:
      ! Get a device pointer for STime, put it to d_STime_p
      istat = cudaHostGetDevicePointer(d_STime_p, C_LOC(Geom%ZDataSoA%STime(1,1,1)), 0)
      ! Translate that C pointer to the fortran array with given dimensions
      call c_f_pointer(d_STime_p, d_STime, [QuadSet%Groups,Size%ncornr, Size%nangSN] )
   endif

   !call get_environment_variable ("ZEROCOPY", envstring )


   ! if(0) then
   ! !if(TRIM(envstring) .eq. "True") then
   !    print *, "Zero copy setup of d_psib"
   !    ! Get a device pointer for psib, put it to d_psib_p
   !    istat = cudaHostGetDevicePointer(d_psib_p, C_LOC(psib(1,1,1)), 0)
   !    ! Translate that C pointer to the fortran array with given dimensions
   !    call c_f_pointer(d_psib_p, d_psib, [QuadSet%Groups,Size%nbelem,QuadSet%NumAngles] )
   !    ! conversion from c to f not really need anymore with cuda c call.
   ! else !NOT WORKING YET.
   !    ! explicitly batch d_psib onto GPU. 
   !    print *, "explicitly batching d_psib onto GPU."
   !    ! Allocate an array to hold subset of angles:
   !    allocate(d_psib(Groups,nbelem,BATCHSIZE))
   ! endif

   theOMPLoopTime=0.0

   call mpi_comm_rank(mpi_comm_world, myrank, info)
   
   ! Double check the rank and device mapping:
   !istat = cudaGetDevice(devnum)
   !write(0,*) 'Rank = ', myrank, 'gpu device = ', devnum, 'istat = ', istat


   !if (myrank .eq. 0) write(0,*) ' groups, ncornr, nbelem, angles, NangBin, NumBin = ', groups, ncornr, nbelem, angles, NangBin, NumBin
   
   ! Set the Cache configuration for the GPU (use more L1, less shared)
   istat = cudaDeviceGetCacheConfig(cacheconfig)
   print *, "cacheconfig =", cacheconfig
   if (cacheconfig .eq. cudaFuncCachePreferShared) then
      print *, "L1 set for shared memory usage"
   elseif (cacheconfig .eq. cudaFuncCachePreferL1) then
      print *, "L1 set for hardware caching (prefer L1)."
   else
      print *, "other L1 configuration present."
   endif

   cacheconfig = cudaFuncCachePreferL1
   !istat = cudaDeviceSetCacheConfig(cacheconfig)
   !istat = cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte)


   ! Create streams that can overlap 
   do s = 1, nStreams
      istat = cudaStreamCreate(stream(s))
   enddo

   ! Create events to synchronize among different streams (only need even ones)
   do s = 1, nStreams, 2
      istat = cudaEventCreate(HtoDdone(s))
      istat = cudaEventCreate(SweepFinished(s))
      istat = cudaEventCreate(PsiOnHost(s))
      istat = cudaEventCreate(PsibOnDevice(s))
      istat = cudaEventCreate(PsibOnHost(s))
   enddo

!  Loop over angle bins

   if (ipath == 'sweep') then
     call timer_beg('_setflux')
     call setIncidentFlux(psib)
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
       phi(:,:) = zero
     endif

!    Loop over angles, solving for each in turn:
     startOMPLoopTime = MPI_WTIME()
     
     ! By default STime does not need to computed.
     calcSTime = .false.

     ! If this is first temp and intensity iteration, need to calculate STime and update to host:
     if (intensityIter == 1 .and. tempIter == 1) then
        ! compute STime from initial d_psi
        calcSTime = .true.
     else
        ! STime already computed,
        calcSTime = .false.
     endif


     call timer_beg('_anglebins')     
     AngleBin: do binRecv=1,QuadSet% NumBin
       binSend = QuadSet% SendOrder(binRecv)
       binSend_next = QuadSet% SendOrder(binRecv+1) ! binSend on next iteration
       NangBin = QuadSet% NangBinList(binSend)
       NangBin_next = QuadSet% NangBinList(binSend_next)
       



       
       !Stage batches of (psi,psib,phi, STime) into GPU
       FirstOctant: if (binRecv == 1) then
          !for other bins, will begin staging in the data at the end of prev
          !iteration of the loop

          ! loop over batches of angles within the angle bin (Later will want loop variable to be batch number.)
          do batch=1,Nbatches
             ! lower angle index for this batch
             mm1=(batch-1)*BATCHSIZE+1
             ! get the angles upper index for this batch
             mm2=min(batch*BATCHSIZE,NangBin)
             ! True size of the batch (not always BATCHSIZE for last batch if not nicely divisible)
             anglebatch=mm2-mm1+1

             ! stream is determined by octant and batch (2 unique streams per batch)
             s = batch*2-1 + 2*Nbatches*(binRecv-1)

             !print *,"first octant:", "mm1 = ", mm1, "mm2 = ", mm2, "anglebatch = ", anglebatch, "s = ", s

             ! stream s wait until data from previous stream is done moving (event HtoDdone(s-2))
             istat=cudaStreamWaitEvent( stream(s), HtoDdone(s-2), 0)

             ! May want to wait until previous psib is moved as well:
             istat = cudaStreamWaitEvent(stream(s), PsibOnDevice(s-2), 0)

             !istat=cudaDeviceSynchronize()

             ! move anglebatch section of psi to d_psi(batch), which has room for BATCHSIZE angles of psi
             istat=cudaMemcpyAsync(d_psi(1,1,1,batch),                 &
                  psi(1,1,QuadSet%AngleOrder(mm1,binSend)), &
                  QuadSet%Groups*Size%ncornr*anglebatch, stream(s) )

             !istat=cudaDeviceSynchronize()

             ! If this is first temp and intensity iteration, need to calculate STime and update to host:
             if (calcSTime == .true.) then
                ! compute STime from initial d_psi
                call computeSTime(d_psi(1,1,1,batch), d_STimeBatch(1,1,1,batch), anglebatch, stream(s) )
                ! Update STime to host
                istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend)), &
                     d_STimeBatch(1,1,1,batch), &
                     QuadSet%Groups*Size%ncornr*anglebatch, stream(s) ) ! can be another stream later?
             else
                ! STime already computed, just need to move section of STime to device
                istat=cudaMemcpyAsync(d_STimeBatch(1,1,1,batch),                 &
                     Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend)), &
                     QuadSet%Groups*Size%ncornr*anglebatch, stream(s) )
             endif

             ! record when HtoD movements are completed for stream s
             istat=cudaEventRecord(HtoDdone(s), stream(s) )

          enddo

          ! Now the data movements have been pipelined for the buffers that fit on GPU.

       endif FirstOctant


       
       ! loop over batches of angles within the angle bin
       do batch=1,Nbatches
          ! lower angle index for this batch
          mm1=(batch-1)*BATCHSIZE+1
          ! get the angles upper index for this batch
          mm2=min(batch*BATCHSIZE,NangBin)
          ! True size of the batch (not always BATCHSIZE for last batch if not nicely divisible)
          anglebatch=mm2-mm1+1

          ! stream is determined by octant and batch (2 unique streams per batch)
          s = batch*2-1 + 2*Nbatches*(binRecv-1)
       
          !print *,"compute section:", "mm1 = ", mm1, "mm2 = ", mm2, "anglebatch = ", anglebatch, "s = ", s

          !istat=cudaDeviceSynchronize()
          call nvtxStartRange("snreflect")
          ! Set angular fluxes for reflected angles
          do mm=mm1,mm2
             call snreflect(QuadSet%AngleOrder(mm,binSend), PSIB)
          enddo
          call nvtxEndRange

          !Stage batch of psib into GPU
          ! move anglebatch section of psib to d_psib, which has room for BATCHSIZE angles of psib
          istat=cudaMemcpyAsync(d_psibBatch(1,1,1,batch),                 &
               psib(1,1,QuadSet%AngleOrder(mm1,binSend)), &
               QuadSet%Groups*Size%nbelem*anglebatch, stream(s) )
          ! record when psib is on device
          istat=cudaEventRecord(PsibOnDevice(s), stream(s) )
          
          
          ! Do not launch sweep kernel in stream s+1 until HtoD transfer in stream s is done.
          istat = cudaStreamWaitEvent(stream(s+1), HtoDdone(s), 0)

          ! Also wait until psib is on device.
          ! Later can run Afpz kernel befor this
          istat = cudaStreamWaitEvent(stream(s+1), PsibOnDevice(s), 0)

          ! currently need this synch to make sure psib updated in snreflect before zero copy to device?
          !istat=cudaDeviceSynchronize()
          
          !        Sweep the mesh, calculating PSI for each corner; the
          !        boundary flux array PSIB is also updated here.
          !        Mesh cycles are fixed automatically.
          if(0) then ! call CUDA fortran version
             call snswp3d(     anglebatch,                     &
                              Size%nzones,               &
                              QuadSet%Groups,            &
                              Size%ncornr,               &
                              QuadSet%NumAngles,         &
                              QuadSet%d_AngleOrder(mm1,binSend),        & ! only need angle batch portion
                              Size%maxCorner,            &
                              Size%maxcf,                &
                              binRecv,                   &
                              NangBin,                   &
                              Size%nbelem,                &
                              QuadSet%d_omega,             &
                              Geom%ZDataSoA%nCorner,                &
                              Geom%ZDataSoA%nCFaces,                &
                              Geom%ZDataSoA%c0,                &
                              Geom%ZDataSoA%A_fp,                &
                              Geom%ZDataSoA%A_ez,                &
                              Geom%ZDataSoA%Connect,             &
                              Geom%ZDataSoA%STotal,              &
                              !Geom%ZDataSoA%STime,               &
                              d_STimeBatch, &
                              Geom%ZDataSoA%Volume,             &
                              d_psi,                      &  ! only want angle batch portion
                              d_psib,                      &
                              QuadSet%d_next,              &
                              QuadSet%d_nextZ,             &
                              Geom%ZDataSoA%Sigt,                &
                              Geom%ZDataSoA%SigtInv,             &
                              QuadSet%d_passZstart,              &
                              stream(1))
          else ! Call CUDA c version
             call snswp3d_c(     anglebatch,                     &
                              Size%nzones,               &
                              QuadSet%Groups,            &
                              Size%ncornr,               &
                              QuadSet%NumAngles,         &
                              QuadSet%d_AngleOrder(mm1,binSend),        & ! only need angle batch portion
                              Size%maxCorner,            &
                              Size%maxcf,                &
                              binRecv,                   &
                              NangBin,                   &
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
                              d_STimeBatch(1,1,1,batch),          &
                              d_STime,                            &
                              Geom%ZDataSoA%Volume,             &
                              d_psi(1,1,1,batch),                      &  ! only want angle batch portion
                              d_psibBatch(1,1,1,batch),                      &
                              QuadSet%d_next,              &
                              QuadSet%d_nextZ,             &
                              Geom%ZDataSoA%Sigt,                &
                              Geom%ZDataSoA%SigtInv,             &
                              QuadSet%d_passZstart,        &
                              calcSTime,                  &
                              Size%tau,             &
                              stream(s+1)           &
                              )
          endif
          
          ! record when stream s+1 finishes s sweep (as opposed to s+2)
          istat=cudaEventRecord(SweepFinished(s), stream(s+1) )

          !istat=cudaDeviceSynchronize()

          if (ipath == 'sweep') then
             call timer_beg('__snmoments')
             ! snmoments only reads d_psi, produces d_phi
             call snmomentsD(d_psi(1,1,1,batch), d_phi, QuadSet%d_Weight,     &
                  QuadSet%d_AngleOrder(mm1,binSend),      &
                  anglebatch, stream(s+1)) ! GPU version, one batch at a time
             call timer_end('__snmoments')
          endif

          ! before moving DtoH psi, sweep needs to complete
          istat=cudaStreamWaitEvent(stream(s), SweepFinished(s) , 0)
          
          ! make sure psi in previous stream done before moving psi in new stream (otherwise they share bandwidth)
          istat=cudaStreamWaitEvent(stream(s), PsiOnHost(s-2) , 0)

          !istat=cudaDeviceSynchronize()

          istat=cudaMemcpyAsync(psi(1,1,QuadSet%AngleOrder(mm1,binSend)), &
               d_psi(1,1,1,batch), &
               QuadSet%Groups*Size%ncornr*anglebatch, stream(s) )

          istat=cudaEventRecord(PsiOnHost(s), stream(s))

       enddo

       NotLastOctants: if (binRecv < QuadSet% NumBin) then
          ! If not the last bin (octant), pre-stage data for next angle bin 

          ! loop over batches of angles within the angle bin (Later will want loop variable to be batch number.)
          do batch=1,Nbatches
             ! lower angle index
             mm1=(batch-1)*BATCHSIZE+1
             ! get the angles upper index for this batch for next octant
             mm2=min(batch*BATCHSIZE,NangBin_next)
             ! True size of the batch (not always BATCHSIZE for last batch if not nicely divisible)
             anglebatch_next=mm2-mm1+1

             ! the stream here corresponds to the next octant streams (binRecv+1-1)
             s = batch*2-1 + Nbatches*2*(binRecv)

             !print *,"NotLastOctants:", "mm1 = ", mm1, "mm2 = ", mm2, "anglebatch_next = ", anglebatch_next, "s = ", s

             ! wait in stream s for the data transfer in the corresponding previous stream to finish.
             istat=cudaStreamWaitEvent(stream(s), HtoDdone(s-2) , 0)

             ! most importantly, cannot overwrite psi before it has been moved off device
             istat=cudaStreamWaitEvent(stream(s), PsiOnHost(s-2*Nbatches) , 0)

             !istat=cudaDeviceSynchronize()

              ! move anglebatch section of psi to d_psi, which has room for BATCHSIZE angles of psi
             istat=cudaMemcpyAsync(d_psi(1,1,1,batch),        &
                  psi(1,1,QuadSet%AngleOrder(mm1,binSend_next)), &
                  QuadSet%Groups*Size%ncornr*anglebatch_next, stream(s) )


             !istat=cudaDeviceSynchronize()

             ! If this is first temp and intensity iteration, need to calculate STime and update to host:
             if (calcSTime == .true.) then
                ! doing this in the kernel with zero-copy of STime for now. Only for calcSTime
                ! compute STime from initial d_psi
                call computeSTime(d_psi(1,1,1,batch), d_STimeBatch(1,1,1,batch), anglebatch, stream(s) )
                ! Update STime to host
                istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend_next)), &
                     d_STimeBatch(1,1,1,batch), &
                     QuadSet%Groups*Size%ncornr*anglebatch_next, stream(s) ) ! can be another stream later?
             else
                ! STime already computed, just need to move section of STime to device
                istat=cudaMemcpyAsync(d_STimeBatch(1,1,1,batch),                 &
                     Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend_next)), &
                     QuadSet%Groups*Size%ncornr*anglebatch_next, stream(s) )
             endif

             ! below is not working yet
             if(0) then
                !if( TRIM(envstring) .ne. "True" ) then
                !Stage batch of psib into GPU
                              
                ! istat=cudaMemcpyAsync(d_psib(1,1,1,batch),                 &
                !      psib(1,1,QuadSet%AngleOrder(mm1,binSend_next)), &
                !      QuadSet%Groups*Size%nbelem*anglebatch_next, stream(s) )
                
             endif


             !istat=cudaDeviceSynchronize()

             ! record when HtoD movements are completed for stream s
             istat = cudaEventRecord(HtoDdone(s), stream(s) )

          enddo
          
       endif NotLastOctants
       
       ! loop over batches of angles within the angle bin
       do batch=1,Nbatches
          ! lower angle index for this batch
          mm1=(batch-1)*BATCHSIZE+1
          ! get the angles upper index for this batch
          mm2=min(batch*BATCHSIZE,NangBin)
          ! True size of the batch (not always BATCHSIZE for last batch if not nicely divisible)
          anglebatch=mm2-mm1+1

          ! stream is determined by octant and batch (2 unique streams per batch)
          s = batch*2-1 + 2*Nbatches*(binRecv-1)

          !print *,"CPU work:", "mm1 = ", mm1, "mm2 = ", mm2, "anglebatch = ", anglebatch, "s = ", s

          ! if ready then set exit flux and move psi.
          istat=cudaEventSynchronize( PsiOnHost(s) )
          
          !istat = cudaDeviceSynchronize()

          !call timer_beg('__setExitFlux')
          call nvtxStartRange("setExitFlux")
          call setExitFlux(anglebatch, QuadSet%AngleOrder(mm1,binSend), psi, psib)
          call nvtxEndRange
          !call timer_end('__setExitFlux')

       enddo

       istat=cudaDeviceSynchronize()

       !      Exchange Boundary Fluxes
       ! these need to become non-blocking

       call timer_beg('_exch')
       call nvtxStartRange("exchange")
       call exchange(PSIB, binSend, binRecv) 
       call nvtxEndRange
       call timer_end('_exch')

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

   ! Destroy streams 
   do s = 1, nStreams
      istat = cudaStreamDestroy(stream(s))
   enddo



!  Update the scaler flux

   if (ipath == 'sweep') then
     call restoreCommOrder(QuadSet)
   endif

   angleLoopTime = angleLoopTime + theOMPLoopTime


   return
   end subroutine snflwxyz


