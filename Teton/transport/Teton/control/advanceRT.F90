!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   ADVANCERT - Save zone-average quantities from previous cycle for   *
!               delta-t calculation.  Convert specific radiation       *
!               intensity (i.e. per unit mass) to intensity (per       *
!               unit volume) before the transport calculation.         *
!                                                                      *
!   Input:   tez,trz                                                   *
!                                                                      *
!   Output:  TEZN,TRZN                                                 *
!                                                                      *
!***********************************************************************
 
   subroutine advanceRT(dtrad, psi, PHI, psib) 


   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use Material_mod
   use QuadratureList_mod
   use Quadrature_mod
   use Editor_mod
   use ZoneData_mod

   use GPUhelper_mod
   use snswp3d_mod

   implicit none
   include 'mpif.h'

!  Arguments

   real(adqt), intent(in)    :: dtrad

   real(adqt), intent(inout) :: psib(Size%ngr,Size%nbelem,Size%nangSN), &
                                psi(Size%ngr,Size%ncornr,Size%nangSN), &
                                Phi(Size%ngr,Size%ncornr)

!  Local

   logical(kind=1), save :: first_time = .true.
   integer    :: slot, binRecv, mm1
   integer    :: myrank, info
   integer(kind=8) :: sweep_mem 

   integer    :: ia, ic, ig
   integer    :: c, c0, nCorner, zone
   integer    :: nzones, ngr, numAngles, NumBin

   real(adqt) :: factor, tfloor, tr4min, ratio, eradBOC 
   real(adqt) :: deltaVolume, aveVolume, volumeRatio
   real(adqt) :: sumRad, PhiAve, Tstar

!  Constants

   nzones     = Size%nzones
   ngr        = Size%ngr
   tfloor     = Size%tfloor

   QuadSet    => getQuadrature(Quad,1)
   numAngles  = QuadSet%NumAngles
   NumBin     = QuadSet% NumBin

!***********************************************************************
!  Make corner temperatures consistent with zone averages obtained     *
!  from the host code.                                                 *
!                                                                      *
!  Advance zone temperatures [set old = new]                           *
!***********************************************************************

   call timer_beg('_ZoneLoop0')

!$omp parallel do private(zone,Z,nCorner,c0,Tstar,c,ratio)
   ZoneLoop: do zone=1,nzones

     Z       => getZoneData(Geom, zone)
                                                                                                   
     nCorner = Z% nCorner
     c0      = Z% c0
                                                                                                   
     Tstar = zero
     do c=1,nCorner
       Tstar = Tstar + Z% Volume(c)*Mat%tec(c0+c)
     enddo
      Tstar = Tstar/Z% VolumeZone

     Mat%trz(zone)  = max( Mat%trz(zone), tfloor )
     Mat%tez(zone)  = max( Mat%tez(zone), tfloor )

     Mat%tezn(zone) = Mat%tez(zone)
     ratio          = Mat%tez(zone)/Tstar

     do c=1,nCorner
       Mat%tec(c0+c) = max( ratio*Mat%tec(c0+c), tfloor )
     enddo

   enddo ZoneLoop

   call timer_end('_ZoneLoop0')

!***********************************************************************
!  Set the scaler intensity                                            *
!***********************************************************************

   call mpi_comm_rank(mpi_comm_world, myrank, info)

   ! GPU tranfer related stuff:
   if (first_time) then
      ! Create streams that can overlap 
      istat = cudaStreamCreate(transfer_stream1)
      istat = cudaStreamCreate(transfer_stream)
      istat = cudaStreamCreate(kernel_stream)
      
      call GPUmemRequirements(&
           psib, &
           psi, &
           phi, &
           Geom%ZDataSoA%STime, &
           QuadSet%next, &
           Geom%ZDataSoA%omega_a_fp,&
           sweep_mem)

      if(myrank == 0) print *,"sweep_mem =",sweep_mem

      call InitDeviceBuffers()

      ! Set up dummy previous to prime the loop for first time.
      binRecv = 1 !loop goes from 8..1, so when current is 8, previous is 1.
      previous%bin       = QuadSet%SendOrder(binRecv) ! same as old binSend
      previous%batch     = binRecv    ! batch 1,2,3... sequentual ordering
      previous%NangBin   = QuadSet% NangBinList( previous%bin )
      previous%anglebatch =previous%NangBin ! later can be subset of angles in a bin.
      !previous%Angles=QuadSet%AngleOrder(mm1,previous%bin)
      
      ! set previous pointers to fist slot, but slots remain unowned because there is no data in them yet.
      slot=1
      previous%psi => psi_storage(slot)
      previous%psib => psib_storage(slot)
      previous%STime => STime_storage(slot)
      previous%omega_A_fp => omega_A_fp_storage(slot)
      previous%omega_A_ez => omega_A_ez_storage(slot)

      first_time = .false.
   endif


      ! ! Set up dummy previous to prime the loop for first time.
      ! binRecv = 1 !loop goes from 8..1, so when current is 8, previous is 1.
      ! previous%bin       = QuadSet%SendOrder(binRecv) ! same as old binSend
      ! previous%batch     = binRecv    ! batch 1,2,3... sequentual ordering
      ! previous%NangBin   = QuadSet% NangBinList( previous%bin )
      ! previous%anglebatch =previous%NangBin ! later can be subset of angles in a bin.
      ! !previous%Angles=QuadSet%AngleOrder(mm1,previous%bin)
      
      ! ! set previous pointers to fist slot, but slots remain unowned because there is no data in them yet.
      ! slot=1
      ! previous%psi => psi_storage(slot)
      ! previous%psib => psib_storage(slot)
      ! previous%STime => STime_storage(slot)
      ! previous%omega_A_fp => omega_A_fp_storage(slot)
      ! previous%omega_A_ez => omega_A_ez_storage(slot)



   ! phi has to be set to zero before being used to accumultate. This
   ! used to happen in snmoments in the original code, but since we
   ! do a bin at a time now, we need to set to zero outside of bin loop.

   d_phi=0 ! could try something async memset here.

   call nvtxStartRange("createEvents")
   ! Create events to synchronize among different streams
   call CreateEventsWithFlags()
   call nvtxEndRange

   do slot=1, numSTime_buffers
      ! mark the data as un-owned since a new timestep has made device version stale:
      !d_psi(buffer)% owner = 0
      STime_storage(slot)% owner = 0
   enddo


   ! print *, "--------------"
   ! print *, "SendOrder(:) = ", QuadSet% SendOrder(1:8)
   ! print *, "--------------"

   ! debug1 REMOVE ME?
   !istat = cudaDeviceSynchronize()

   call timer_beg('_snmoments1')

   ! step through the bins from top to bottom, so when sweep starts at bottom data will be on GPU.
   do binRecv=NumBin, 1, -1

      ! which bin is the current bin being worked on?
      current%bin       = QuadSet%SendOrder(binRecv) ! The bin being worked on
      current%batch     = binRecv    ! batch 1,2,3... sequentual ordering
      current%NangBin   = QuadSet% NangBinList( current%bin )
      current%anglebatch =current%NangBin ! later can be subset of angles in a bin.
      
      !print *, "current%bin       =", current%bin
      !print *, "current%batch     =", current%batch

      !   lower angle index for this batch
      mm1=1

      ! get a pointer to where current bin of psi should be moved (or already exists) on GPU
      call checkDataOnDevice(current%psi, psi_storage, current%bin, previous%psi%slot)

      ! check/move current batch of psi to GPU
      call MoveDataOnDevice(current%psi, psi_storage, psi, current%bin, current%batch, &
           previous%psi%slot,mm1, QuadSet%Groups*Size%ncornr*current%anglebatch, transfer_stream, &
           Psi_OnDevice )


      ! setup pointers to omega_A stuff.
      call checkDataOnDeviceDot(current%omega_A_fp,omega_A_fp_storage, current%bin, previous%omega_A_fp%slot)
      call checkDataOnDeviceDot(current%omega_A_ez,omega_A_ez_storage, current%bin, previous%omega_A_ez%slot)


      ! while this is transferring in, compute mesh normal stuff
      if(1) then
         ! CUDA fortran version
         call fp_ez_f(     current%anglebatch,                     &
              Size%nzones,               &
              Size%ncornr,               &
              QuadSet%NumAngles,         &
              QuadSet%d_AngleOrder(mm1,current%bin),        & ! only need angle batch portion
              Size%maxCorner,            &
              Size%maxcf,                &
              current%NangBin,                   &
              Size%nbelem,                &
              QuadSet%d_omega,             &
              current%omega_A_fp%data,                &
              current%omega_A_ez%data,                &
              QuadSet%d_next,              &
              QuadSet%d_nextZ,             &
              QuadSet%d_passZstart,        &
              kernel_stream           &
              )
      else
         ! ! CUDA C version
         ! call fp_ez_c(     current%anglebatch,                     &
         !      Size%nzones,               &
         !      QuadSet%Groups,            &
         !      Size%ncornr,               &
         !      QuadSet%NumAngles,         &
         !      QuadSet%d_AngleOrder(mm1,current%bin),        & ! only need angle batch portion
         !      Size%maxCorner,            &
         !      Size%maxcf,                &
         !      current%NangBin,                   &
         !      Size%nbelem,                &
         !      QuadSet%d_omega,             &
         !      Geom%ZDataSoA%nCorner,                &
         !      Geom%ZDataSoA%nCFaces,                &
         !      Geom%ZDataSoA%c0,                &
         !      Geom%ZDataSoA%A_fp,                &
         !      current%omega_A_fp%data,                &
         !      Geom%ZDataSoA%A_ez,                &
         !      current%omega_A_ez%data,                &
         !      Geom%ZDataSoA%Connect,             &
         !      Geom%ZDataSoA%Connect_reorder,             &
         !      QuadSet%d_next,              &
         !      QuadSet%d_nextZ,             &
         !      QuadSet%d_passZstart,        &
         !      kernel_stream           &
         !      )
      endif

      ! Better mark that these fp are owned
      current%omega_A_fp%owner = current%bin
      current%omega_A_ez%owner = current%bin

      istat=cudaEventRecord(AfpFinished( current%batch ), kernel_stream )

      ! transfer stream should wait until Afp is done calculating.
      istat = cudaStreamWaitEvent(transfer_stream, AfpFinished( current%batch ), 0)

      if(numDot_buffers /= 8) then
         ! not enough buffers to keep A_fp on device throughout computation.
         ! need to move fp stuff back to host.
         istat = cudaMemcpyAsync(Geom%ZDataSoA%omega_A_fp(1,1,1,QuadSet%AngleOrder(mm1,current%bin)), &
              current%omega_A_fp%data(1,1,1,1), &
              Size% nzones*Size% maxCorner*Size% maxcf*current%anglebatch, transfer_stream)



         ! need to move ez stuff back to host here.
         istat = cudaMemcpyAsync(Geom%ZDataSoA%omega_A_ez(1,1,1,QuadSet%AngleOrder(mm1,current%bin)), &
              current%omega_A_ez%data(1,1,1,1), &
              Size% nzones*Size% maxCorner*Size% maxcf*current%anglebatch, transfer_stream)

      endif


      ! have kernel stream wait until transfer of psi to device
      istat = cudaStreamWaitEvent(kernel_stream, Psi_OnDevice( current%batch ), 0)

      call timer_beg('__snmoments')
      ! snmoments only reads d_psi, produces d_phi
      call snmomentsD(current%psi%data(1,1,1), d_phi, QuadSet%d_Weight,     &
             QuadSet%d_AngleOrder(mm1,current%bin),      &
             current%anglebatch, kernel_stream) ! GPU version, one batch at a time
      call timer_end('__snmoments')

      ! Record when snmoments finishes
      istat=cudaEventRecord(snmomentsFinished( current%batch ), kernel_stream )

      ! scale current batch of psi
      call scalePsibyVolume(current%psi%data(1,1,1), Geom%d_GPU_ZData, current%anglebatch, kernel_stream )  

      ! print *, "called scalebyvolume"

      ! don't start Stime until previous Stime is on host. (not needed because transfer stream blocks)


      !print *, "setting up pointers for STime bin ", current%bin

      ! set up STime pointers before computing STime
      call checkDataOnDevice(current%STime, STime_storage, current%bin, previous%STime%slot)

      !print *, "current%STime points to slot", current%STime%slot

      ! compute STime from initial d_psi
      ! (this is done again in snflw, but does not hurt now as transfer dominates here anyway, and better when data fits)
      call computeSTime(current%psi%data(1,1,1), current%STime%data(1,1,1), current%anglebatch, kernel_stream )

      istat=cudaEventRecord(STimeFinished( current%batch ), kernel_stream )
        
      ! need to record who's batch of STime is held in this device buffer:
      current%STime% owner = current%bin

      ! better not set the boundary until previous batch psib is on host,
      ! but I think this is enforced automatically.
      !istat = cudaStreamWaitEvent(kernel_stream, psib_OnHost( previous%batch ), 0)        

      ! setup pointer for psib
      call checkDataOnDevice(current%psib, psib_storage, current%bin, previous%psib%slot)

      ! set the boundary here
      call setbdyD(current%anglebatch, QuadSet%d_AngleOrder(mm1,current%bin), &
             current%psi%data(1,1,1), &
             current%psib%data(1,1,1), kernel_stream)


      ! only prestage next batch psi if not the last bin, 
      ! because we want psi 2 and 1 on device, not 1 and 8
      if(binRecv /= 1) then

         ! numbin is used instead of numBatches because actually don't always do all 8 batches.
         next%batch     = 1+modulo(binRecv-2,QuadSet% NumBin)  ! next bin is really the previous since index counts down.
         next%bin       = QuadSet%SendOrder(next%batch) ! same as old binSend
         next%NangBin   = QuadSet% NangBinList( next%bin )
         next%anglebatch = next%NangBin ! later can be subset of angles in a bin.

         !print *, "next%bin       =", next%bin
         !print *, "next%batch     =", next%batch


         call checkDataOnDevice(next%psi, psi_storage, next%bin, current%psi%slot)

         ! move next batch of psi to GPU (will happen while above kernels are running)
         call MoveDataOnDevice(next%psi, psi_storage, psi, next%bin, next%batch, &
              current%psi%slot, mm1, QuadSet%Groups*Size%ncornr*next%anglebatch, transfer_stream, &
              Psi_OnDevice )


      endif

      ! Do not move psib out until setbdyD has finished
      istat = cudaDeviceSynchronize()

      ! move current batch psib to host
      istat=cudaMemcpyAsync(psib(1,1,QuadSet%AngleOrder(mm1,current%bin)), &
             current%psib%data(1,1,1), &
             QuadSet%Groups*Size%nbelem*current%anglebatch, transfer_stream ) 

      istat=cudaEventRecord(psib_OnHost( current%batch ), transfer_stream )

        
      ! transfer stream waits for STime to be computed:
      istat = cudaStreamWaitEvent(transfer_stream, STimeFinished(current%batch), 0)

      if ( .not. fitsOnGPU ) then

         ! move current batch STime to host
         istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,current%bin)), &
              current%STime%data(1,1,1), &
              QuadSet%Groups*Size%ncornr*current%anglebatch, transfer_stream ) 
         
         !print *, "moving STime to host for bin ", current%bin
         !print *, "the owner was ", current%STime%owner
         
      endif

      if(binRecv /= 1) then
         ! when loop cycles, current becomes previous
         previous=current
         ! did this even work?
         !print *, "previous%psi%slot = ", previous%psi%slot
      else
         ! on the last loop, current stays what it is, and previous is batch 2
      endif
         
      ! debug1
      !istat = cudaDeviceSynchronize()
      

   enddo

   
   
   
   ! transfer stream waits for snmoments calc to be finished (the last one)
   istat = cudaStreamWaitEvent(transfer_stream1, snmomentsFinished(current%batch), 0)
   
   !print *, "moving d_phi to host"

   ! move d_phi data to host:
   istat=cudaMemcpyAsync(phi(1,1), &
        d_phi(1,1), &
        QuadSet%Groups*Size%ncornr, transfer_stream1 )

   istat=cudaEventRecord( phi_OnHost, transfer_stream1 )
   ! Actually could just do a synchronize on event phi on host...
   !istat = cudaDeviceSynchronize()

   ! CPU code should wait until phi is on the host before using it
   istat=cudaEventSynchronize( phi_OnHost )


   !print *, "d_phi finished move"

   call timer_end('_snmoments1')

!***********************************************************************
!  Compute the work done on radiation field due to volume changes.     *
!  This is an external source rate in the radiation transport equation.*
!***********************************************************************

   ! this is never executed in CORAL UMT run, but you would need phi unscaled here,
   ! which is why I do not scale phi on the GPU above. (could scale psi before snmoments)
   if (Size%radForceMultiplier > zero) then

     Mat%qext(:,:) = zero

     factor = -third*Size%radForceMultiplier/(dtrad*speed_light)

   call timer_beg('_ZoneLoop1')

     ZoneLoop1: do zone=1,nzones

       Z    => getZoneData(Geom, zone)

       nCorner = Z% nCorner
       c0      = Z% c0

       do c=1,nCorner
         deltaVolume      =       Z% Volume(c) - Z% VolumeOld(c) 
         aveVolume        = half*(Z% Volume(c) + Z% VolumeOld(c))
         volumeRatio      = factor*deltaVolume/aveVolume

         Mat%qext(:,c0+c) = volumeRatio*Phi(:,c0+c) 
       enddo

     enddo ZoneLoop1

     call timer_end('_ZoneLoop1')

   endif

!***********************************************************************
!  Scale the radiation field to account for volume changes and         *
!  tally beginning-of-cycle radiation energy                           *
!***********************************************************************
 
   eradBOC = zero
   tr4min  = tfloor*tfloor*tfloor*tfloor

   call timer_beg('_ZoneLoop2')

!$omp parallel do private(zone,PhiAve,Z,nCorner,c0,c,volumeRatio,ia,sumRad) reduction(+:eradBOC)
   ZoneLoop2: do zone=1,nzones

     PhiAve  = zero

     Z       => getZoneData(Geom, zone)
     nCorner = Z% nCorner
     c0      = Z% c0

     do c=1,nCorner
       volumeRatio = Z% VolumeOld(c)/Z% Volume(c)

       Phi(:,c0+c) = Phi(:,c0+c)*volumeRatio

       ! This scaling of psi is now done with psi on the GPU in snflwxyz.F90.
       ! or comment that out and do here:
       !do ia=1,numAngles
       !   psi(:,c0+c,ia) = psi(:,c0+c,ia)*volumeRatio
       !enddo

       sumRad = zero
       do ig=1,ngr
         sumRad = sumRad + Phi(ig,c0+c)
       enddo

       PhiAve = PhiAve + Z% Volume(c)*sumRad

     enddo


     eradBOC = eradBOC + PhiAve
     PhiAve  = PhiAve/(Z% VolumeZone*rad_constant*speed_light)

     Mat%trzn(zone) = sqrt( sqrt( max(PhiAve, tr4min) ) )

   enddo ZoneLoop2

   call timer_end('_ZoneLoop2')

   if (Size%igeom == 'rz') then
     RadEdit% EnergyRadBOC = two*pi*eradBOC/speed_light
   else
     RadEdit% EnergyRadBOC = eradBOC/speed_light
   endif


   return
   end subroutine advanceRT



