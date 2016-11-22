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
#define BATCHSIZE 96
   subroutine snflwxyz(ipath, PSIB, PSI, PHI, angleLoopTime)

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

!  Arguments

   real(adqt), intent(inout) :: psib(QuadSet%Groups,Size%nbelem,QuadSet%NumAngles)
   real(adqt), intent(inout) :: psi(QuadSet%Groups,Size%ncornr,QuadSet%NumAngles)
   real(adqt), intent(inout) :: Phi(QuadSet%Groups,Size%ncornr),angleLoopTime

   character(len=8), intent(in) :: ipath

!  Local

   integer, parameter :: nStreams = 8
   integer(kind=cuda_stream_kind) :: stream(nStreams)
   integer :: i
   

   integer          :: Angle, mm,mm1,mm2,anglebatch
   integer          :: Groups, fluxIter, ishared
   integer          :: binSend, binRecv, NangBin, istat

   logical (kind=1) :: FluxConverged

   real(adqt)       :: maxFluxError
   real(adqt)       :: startOMPLoopTime, endOMPLoopTime, theOMPLoopTime

   ! zero copy pointers for phi and psib
   type(C_DEVPTR)                    :: d_phi_p
   type(C_DEVPTR)                    :: d_psib_p
   real(adqt), device, allocatable :: d_phi(:,:)
   real(adqt), device, allocatable :: d_psib(:,:,:)

   ! picking up environment variables
   character(len=255) :: envstring

   real(adqt), device :: d_psi(QuadSet%Groups,Size%ncornr,BATCHSIZE)
   type(C_PTR) :: cptr
   type(C_DEVPTR) :: dptr

   integer :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
   integer angles, nbelem, ncornr, NumBin, myrank, info

   integer :: devnum, cacheconfig

!  Convenient Mesh Constants

   Groups = QuadSet%Groups

   angles = QuadSet%NumAngles
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
   

   ! This sets up to allow zero copy use of phi directly on the device:
   ! Get a device pointer for phi, put it to d_phi_p
   istat = cudaHostGetDevicePointer(d_phi_p, C_LOC(phi(1,1)), 0)
   ! Translate that C pointer to the fortran array with given dimensions
   call c_f_pointer(d_phi_p, d_phi, [QuadSet%Groups,Size%ncornr] )

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
   print *, cacheconfig
   if (cacheconfig .eq. cudaFuncCachePreferShared) then
      print *, "L1 set for shared memory usage"
   elseif (cacheconfig .eq. cudaFuncCachePreferL1) then
      print *, "L1 set for hardware caching (prefer L1)."
   else
      print *, "other L1 configuration present."
   endif

   cacheconfig = cudaFuncCachePreferL1
   istat = cudaDeviceSetCacheConfig(cacheconfig)

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
       NangBin = QuadSet% NangBinList(binSend)
       

       call timer_beg('_angleloop')
       ! loop over batches within the angle bin
       AngleLoop: do mm1=1,NangBin,BATCHSIZE

         mm2=min(mm1+BATCHSIZE-1,NangBin)
         anglebatch=mm2-mm1+1

         !Stage batches of (psi,psib,phi) into GPU
         if (binRecv == 1) then
            !for other bins, will begin staging in the data at the end of prev
            !iteration of the loop

            ! move d_psi, which holds BATCHSIZE angles of psi
            istat=cudaMemcpyAsync(d_psi(1,1,1),                 &
                                   psi(1,1,QuadSet%AngleOrder(mm1,binSend)), &
                                   QuadSet%Groups*Size%ncornr*mm2, stream(binRecv) )

         endif

!        Set angular fluxes for reflected angles

         do mm=mm1,mm2
            call snreflect(QuadSet%AngleOrder(mm,binSend), PSIB)
         enddo

!        Sweep the mesh, calculating PSI for each corner; the
!        boundary flux array PSIB is also updated here.
!        Mesh cycles are fixed automatically.

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
                              Geom%ZDataSoA%STime,               &
                              Geom%ZDataSoA%Volume,             &
                              d_psi,                      &  ! only want angle batch portion
                              d_psib,                      &
                              QuadSet%d_next,              &
                              QuadSet%d_nextZ,             &
                              Geom%ZDataSoA%Sigt,                &
                              Geom%ZDataSoA%SigtInv,             &
                              QuadSet%d_passZstart,              &
                              stream(binRecv))


         ! call snswp3d(anglebatch, QuadSet%AngleOrder(mm1,binSend), &
         !              QuadSet%d_AngleOrder(mm1,binSend),           &
         !              QuadSet%d_next,QuadSet%d_nextZ,              &
         !              QuadSet%d_passZstart, d_psi, psib, stream(binRecv))


         ! update psi on the host:
         istat=cudaMemcpyAsync(psi(1,1,QuadSet%AngleOrder(mm1,binSend)), &
              d_psi(1,1,1), &
              QuadSet%Groups*Size%ncornr*mm2, stream(binRecv) )

         if (ipath == 'sweep') then
            call timer_beg('__snmoments')
            ! snmoments only reads d_psi, produces d_phi
            call snmomentsD(d_psi, d_phi, QuadSet%d_Weight,     &
                           QuadSet%d_AngleOrder(mm1,binSend),      &
                           anglebatch, stream(binRecv)) ! GPU version, one slice at a time
            call timer_end('__snmoments')
         endif

         
         ! synchronize to be sure psi is updated on host. (Later should sync to StreamEvent)
         istat = cudaStreamSynchronize(stream(binRecv) )
         !istat = cudaDeviceSynchronize()
         if (istat /= 0) then
            write(*,*) "CUDA stream sync API error:",istat
            stop
         endif

         call setExitFlux(anglebatch, QuadSet%AngleOrder(mm1,binSend), psi, psib)


         if (binRecv < QuadSet% NumBin) then
           ! pre-stage data for next angle bin while exchange is happening
             istat=cudaMemcpyAsync(d_psi(1,1,1),        &
                                   psi(1,1,QuadSet%AngleOrder(mm1,QuadSet%SendOrder(binRecv+1))), &
                                   QuadSet%Groups*Size%ncornr*mm2, stream(binRecv+1) )
         endif


       enddo AngleLoop
       call timer_end('_angleloop')
       endOMPLoopTime = MPI_WTIME()
       theOMPLoopTime = theOMPLoopTime + (endOMPLoopTime-startOMPLoopTime)

!      Exchange Boundary Fluxes

       call timer_beg('_exch')
       call exchange(PSIB, binSend, binRecv) 
       call timer_end('_exch')

     enddo AngleBin

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
     call restoreCommOrder(QuadSet)
   endif

   angleLoopTime = angleLoopTime + theOMPLoopTime


   return
   end subroutine snflwxyz


