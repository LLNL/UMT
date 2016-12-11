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
#define BATCHSIZE 32
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
          A_ez , &
          Connect , &
          STotal , &
          STime , &
          Volume , &
          psic, &
          psib, &
          next, &
          nextZ, &   
          Sigt, &
          SigtInv, &
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
         integer ( c_int ) :: binRecv
         integer ( c_int ) :: nangbin
         integer ( c_int ) :: nbelem
         real ( c_double ),device :: omega(*)
         integer ( c_int ),device :: numCorners(*) 
         integer ( c_int ),device :: numCFaces(*) 
         integer ( c_int ),device :: c0(*)
         real ( c_double ),device :: A_fp(*) 
         real ( c_double ),device :: A_ez(*) 
         integer ( c_int ),device :: Connect(*) 
         real ( c_double ),device :: STotal(*) 
         real ( c_double ),device :: STime(*) 
         real ( c_double ),device :: Volume(*) 
         real ( c_double ),device :: psic(*) 
         real ( c_double ),device :: psib(*) 
         integer ( c_int ),device :: next(*)
         integer ( c_int ),device :: nextZ(*)
         real ( c_double ),device :: Sigt(*) 
         real ( c_double ),device :: SigtInv(*) 
         integer ( c_int ),device :: passZ(*)
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

   integer, parameter :: nStreams = BATCHSIZE
   integer(kind=cuda_stream_kind) :: stream(nStreams)
   integer :: i
   

   integer          :: Angle, mm,mm1,mm2,anglebatch
   integer          :: Groups, fluxIter, ishared
   integer          :: binSend, binSend_next, binRecv, NangBin, istat

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
   

   ! picking up environment variables
   character(len=255) :: envstring

   real(adqt), device :: d_psi(QuadSet%Groups,Size%ncornr,BATCHSIZE)
   real(adqt), device :: d_STimeBatch(QuadSet%Groups,Size%ncornr,BATCHSIZE)
   
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


   if(1) then
   !if(TRIM(envstring) .eq. "True") then
      print *, "Zero copy setup of d_psib"
      ! Get a device pointer for psib, put it to d_psib_p
      istat = cudaHostGetDevicePointer(d_psib_p, C_LOC(psib(1,1,1)), 0)
      ! Translate that C pointer to the fortran array with given dimensions
      call c_f_pointer(d_psib_p, d_psib, [QuadSet%Groups,Size%nbelem,QuadSet%NumAngles] )
      ! conversion from c to f not really need anymore with cuda c call.
   else !NOT WORKING YET.
      ! explicitly batch d_psib onto GPU. 
      print *, "explicitly batching d_psib onto GPU."
      ! Allocate an array to hold subset of angles:
      allocate(d_psib(Groups,nbelem,BATCHSIZE))
   endif

   theOMPLoopTime=0.0

   call mpi_comm_rank(mpi_comm_world, myrank, info)
   
   ! Double check the rank and divice mapping:
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
   do i = 1, nStreams
      istat = cudaStreamCreate(stream(i))
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
     call InitExchange
     call timer_end('_initexch')

     fluxIter = fluxIter + 1

     if (ipath == 'sweep') then
       phi(:,:) = zero
     endif

!    Loop over angles, solving for each in turn:
     startOMPLoopTime = MPI_WTIME()

     
     AngleBin: do binRecv=1,QuadSet% NumBin
       binSend = QuadSet% SendOrder(binRecv)
       binSend_next = QuadSet% SendOrder(binRecv+1) ! binSend on next iteration
       NangBin = QuadSet% NangBinList(binSend)
       

       call timer_beg('_angleloop')
       ! loop over batches within the angle bin (Later will want loop variable to be batch number.)
       AngleLoop: do mm1=1,NangBin,BATCHSIZE

         ! get the angleloop upper bound for this batch
         mm2=min(mm1+BATCHSIZE-1,NangBin)
         ! true size of the batch (not always batchsize for last iteration)
         anglebatch=mm2-mm1+1

         !Stage batches of (psi,psib,phi, STime) into GPU
         if (binRecv == 1) then
            !for other bins, will begin staging in the data at the end of prev
            !iteration of the loop

            ! move anglebatch section of psi to d_psi, which has room for BATCHSIZE angles of psi
            istat=cudaMemcpyAsync(d_psi(1,1,1),                 &
                 psi(1,1,QuadSet%AngleOrder(mm1,binSend)), &
                 QuadSet%Groups*Size%ncornr*anglebatch, stream(1) )


            ! If this is first temp and intensity iteration, need to calculate STime and update to host:
            if (intensityIter == 1 .and. tempIter == 1) then
               ! compute STime from initial d_psi
               call computeSTime(d_psi, d_STimeBatch, anglebatch, stream(1) )
               ! Update STime to host
               istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend)), &
                    d_STimeBatch(1,1,1), &
                    QuadSet%Groups*Size%ncornr*anglebatch, stream(1) )

            else
               ! STime already computed, just need to move section of STime to device
               istat=cudaMemcpyAsync(d_STimeBatch(1,1,1),                 &
                    Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend)), &
                    QuadSet%Groups*Size%ncornr*anglebatch, stream(1) )

               ! could just use pinned memory version instead here.
            endif

         endif

!        Set angular fluxes for reflected angles

         do mm=mm1,mm2
            call snreflect(QuadSet%AngleOrder(mm,binSend), PSIB)
         enddo
            !Stage batch of psib into GPU
            if (binRecv == 1) then
               !for other bins, will begin staging in the data at the end of prev
               !iteration of the loop

               ! This is not working yet.
               if(0) then
               !if( TRIM(envstring) .ne. "True" ) then               
                  ! move anglebatch section of psib to d_psib, which has room for BATCHSIZE angles of psib
                  istat=cudaMemcpyAsync(d_psib(1,1,1),                 &
                       psib(1,1,QuadSet%AngleOrder(mm1,binSend)), &
                       QuadSet%Groups*Size%nbelem*anglebatch, stream(1) )
               endif

         endif

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
                              Geom%ZDataSoA%A_ez,                &
                              Geom%ZDataSoA%Connect,             &
                              Geom%ZDataSoA%STotal,              &
                              !Geom%ZDataSoA%STime,               &
                              d_STimeBatch,                          &
                              Geom%ZDataSoA%Volume,             &
                              d_psi,                      &  ! only want angle batch portion
                              d_psib,                      &
                              QuadSet%d_next,              &
                              QuadSet%d_nextZ,             &
                              Geom%ZDataSoA%Sigt,                &
                              Geom%ZDataSoA%SigtInv,             &
                              QuadSet%d_passZstart,        &
                              stream(1)           &
                              )
        endif

        istat=cudaDeviceSynchronize()

        ! need different stream for every angle...
        ! update psi on the host an angle at a time:
        do mm=mm1,mm2
           istat=cudaMemcpyAsync(psi(1,1,QuadSet%AngleOrder(mm,binSend)), &
                d_psi(1,1,mm-mm1+1), &
                QuadSet%Groups*Size%ncornr, stream(mm-mm1+1) )
        enddo
        
        if (ipath == 'sweep') then
           call timer_beg('__snmoments')
           ! snmoments only reads d_psi, produces d_phi
           call snmomentsD(d_psi, d_phi, QuadSet%d_Weight,     &
                QuadSet%d_AngleOrder(mm1,binSend),      &
                anglebatch, stream(1)) ! GPU version, one batch at a time
           call timer_end('__snmoments')
        endif

        do mm=mm1, mm2
           ! synchronize with stream mm to be sure psi(mm) is on the host
           istat = cudaStreamSynchronize( stream(mm-mm1+1) )
           if (istat /= 0) then
              write(*,*) "CUDA stream sync API error:",istat
              stop
           endif

           call timer_beg('__setExitFlux')
           call setExitFlux_singleangle(QuadSet%AngleOrder(mm,binSend), psi, psib)
           call timer_end('__setExitFlux')

           if (binRecv < QuadSet% NumBin) then
              ! If not the last bin, pre-stage data for next angle bin (an angle at a time)

              istat=cudaMemcpyAsync(d_psi(1,1,mm-mm1+1),        &
                   psi(1,1,QuadSet%AngleOrder(mm,binSend_next)), &
                   QuadSet%Groups*Size%ncornr, stream(mm-mm1+1) )

           endif

        enddo
        
        if (binRecv < QuadSet% NumBin) then
           ! If ! not the last bin, pre-stage data for next angle bin
           ! istat=cudaMemcpyAsync(d_psi(1,1,1),        &
           !      psi(1,1,QuadSet%AngleOrder(mm1,QuadSet%SendOrder(binRecv+1))), &
           !      QuadSet%Groups*Size%ncornr*anglebatch, stream(binRecv+1) )
           ! ! this won't work for when all anglebatch sizes are not the same...need next_anglebatch

           ! If this is first temp and intensity iteration, need to calculate STime and update to host:
           if (intensityIter == 1 .and. tempIter == 1) then
              ! sync all streams (d_psi must be on device to compute STime).
              istat=cudaDeviceSynchronize()
              ! compute STime from initial d_psi
              call computeSTime(d_psi, d_STimeBatch, anglebatch, stream(1) )
              ! Update STime to host
              istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend_next)), &
                   d_STimeBatch(1,1,1), &
                   QuadSet%Groups*Size%ncornr*anglebatch, stream(1) )

           else

              ! STime already computed, just need to move section of STime to device
              istat=cudaMemcpyAsync(d_STimeBatch(1,1,1),                 &
                   Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend_next)), &
                   QuadSet%Groups*Size%ncornr*anglebatch, stream(1) )
           endif

           ! below is not working yet
           if(0) then
              !if( TRIM(envstring) .ne. "True" ) then
              !Stage batch of psib into GPU

              ! move anglebatch section of psib to d_psib, which has room for BATCHSIZE angles of psib
              istat=cudaMemcpyAsync(d_psib(1,1,1),                 &
                   psib(1,1,QuadSet%AngleOrder(mm1,QuadSet%SendOrder(binRecv+1))), &
                   QuadSet%Groups*Size%nbelem*anglebatch, stream(1) )

           endif

        endif

       enddo AngleLoop
       call timer_end('_angleloop')


!      Exchange Boundary Fluxes

       call timer_beg('_exch')
       call exchange(PSIB, binSend, binRecv) 
       call timer_end('_exch')

     enddo AngleBin

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
   do i = 1, nStreams
      istat = cudaStreamDestroy(stream(i))
   enddo



!  Update the scaler flux

   if (ipath == 'sweep') then
     call restoreCommOrder(QuadSet)
   endif

   angleLoopTime = angleLoopTime + theOMPLoopTime


   return
   end subroutine snflwxyz


