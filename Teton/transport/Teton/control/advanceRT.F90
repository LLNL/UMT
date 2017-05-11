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

   implicit none

!  Arguments

   real(adqt), intent(in)    :: dtrad

   real(adqt), intent(inout) :: psib(Size%ngr,Size%nbelem,Size%nangSN), &
                                psi(Size%ngr,Size%ncornr,Size%nangSN), &
                                Phi(Size%ngr,Size%ncornr)

!  Local

   logical(kind=1), save :: first_time = .true.
   integer    :: buffer, binRecv, mm1, mm2


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

   ! GPU tranfer related stuff:
   if (first_time) then
      ! Create streams that can overlap 
      istat = cudaStreamCreate(transfer_stream1)
      istat = cudaStreamCreate(transfer_stream)
      istat = cudaStreamCreate(kernel_stream)

      call InitDeviceBuffers()

      first_time = .false.
   endif

   ! phi has to be set to zero before being used to accumultate. This
   ! used to happen in snmoments in the original code, but since we
   ! do a bin at a time now, we need to set to zero outside of bin loop.

   d_phi=0 ! could try something async memset here.

   call nvtxStartRange("createEvents")
   ! Create events to synchronize among different streams
   call CreateEvents()
   call nvtxEndRange

   do buffer=1, numGPUbuffers
      ! mark the data as un-owned since a new timestep has made device version stale:
      !d_psi(buffer)% owner = 0
      d_STime(buffer)% owner = 0
   enddo

   ! but actually buffer 8


   call timer_beg('_snmoments1')

   ! step through the bins from top to bottom, so when sweep starts at bottom data will be on GPU.
   do binRecv=NumBin, 1, -1
        ! cycle the buffers used to hold the batches.
        ! there will be either 2 buffers if the problem does not fit in GPU, or NumBin if it fits.
        current = 1 + mod(binRecv-1,numGPUbuffers) ! gives 1,2,1,2...  or 1,2,3...7,8
        !next = 1+mod(binRecv,numGPUbuffers)        ! gives 2,1,2,1...  or 2,3,4...8,1 
        next = 1+mod(binRecv-2,numGPUbuffers)  ! next bin is really the previous since index counts down.

        ! each batch corresponds to an angle bin.
        batch(current) = binRecv
        !batch(next) = 1 + mod(binRecv,QuadSet% NumBin)
        batch(next) = 1 + mod(binRecv-2,QuadSet% NumBin) ! next is previous

        
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


        ! check/move current batch of psi to GPU
        call checkDataOnDevice(d_psi, psi, batch, current, mm1, &
             QuadSet%Groups*Size%ncornr*anglebatch(current), transfer_stream, &
             Psi_OnDevice )


        ! while this is transferring in, compute mesh normal stuff
        
        call fp_ez_c(     anglebatch(current),                     &
             Size%nzones,               &
             QuadSet%Groups,            &
             Size%ncornr,               &
             QuadSet%NumAngles,         &
             QuadSet%d_AngleOrder(mm1,binSend(current)),        & ! only need angle batch portion
             Size%maxCorner,            &
             Size%maxcf,                &
             Nangbin(current),                   &
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
             QuadSet%d_next,              &
             QuadSet%d_nextZ,             &
             QuadSet%d_passZstart,        &
             kernel_stream           &
             )

        istat=cudaEventRecord(AfpFinished( batch(current) ), kernel_stream )

        ! transfer stream should wait until Afp is done calculating.
        istat = cudaStreamWaitEvent(transfer_stream, AfpFinished( batch(current) ), 0)


        ! need to move fp stuff back to host.
        istat = cudaMemcpyAsync(Geom%ZDataSoA%omega_A_fp(1,1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
             d_omega_A_fp(1,1,1,1), &
             Size% nzones*Size% maxCorner*Size% maxcf*anglebatch(current), transfer_stream)

        ! need to move ez stuff back to host here.
        istat = cudaMemcpyAsync(Geom%ZDataSoA%omega_A_ez(1,1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
             d_omega_A_ez(1,1,1,1), &
             Size% nzones*Size% maxCorner*Size% maxcf*anglebatch(current), transfer_stream)


        ! have kernel stream wait until transfer of psi to device
        istat = cudaStreamWaitEvent(kernel_stream, Psi_OnDevice( batch(current) ), 0)

        call timer_beg('__snmoments')
        ! snmoments only reads d_psi, produces d_phi
        call snmomentsD(d_psi(current)%data(1,1,1), d_phi, QuadSet%d_Weight,     &
             QuadSet%d_AngleOrder(mm1,binSend(current)),      &
             anglebatch(current), kernel_stream) ! GPU version, one batch at a time
        call timer_end('__snmoments')

        ! Record when snmoments finishes
        istat=cudaEventRecord(snmomentsFinished( batch(current) ), kernel_stream )

        ! scale current batch of psi
        call scalePsibyVolume(d_psi(current)%data(1,1,1), Geom%ZDataSoA%volumeRatio, anglebatch(current), kernel_stream )  

        ! don't start Stime until previous Stime is on host. (not needed because transfer stream blocks)


        ! compute STime from initial d_psi
        ! (this is done again in snflw, but does not hurt now as transfer dominates here anyway, and better when data fits)
        call computeSTime(d_psi(current)%data(1,1,1), d_STime(current)%data(1,1,1), anglebatch(current), kernel_stream )

        istat=cudaEventRecord(STimeFinished( batch(current) ), kernel_stream )
        
        ! need to record who's batch of STime is held in this device buffer:
        d_STime(current)% owner = batch(current)

        ! better not set the boundary until previous batch psib is on host,
        ! but I think this is enforced automatically.
        
        ! set the boundary here
        call setbdyD(anglebatch(current), QuadSet%d_AngleOrder(mm1,binSend(current)), &
             d_psi(current)%data(1,1,1), &
             d_psibBatch(1,1,1,current), kernel_stream)


        ! only prestage next batch psi if not the last bin, 
        ! because we want psi 2 and 1 on device, not 1 and 8
        if(binRecv /= 1) then
           ! check/move next batch of psi to GPU (will happen while above kernels are running)
           call checkDataOnDevice(d_psi, psi, batch, next, mm1, &
                QuadSet%Groups*Size%ncornr*anglebatch(next), transfer_stream, &
                Psi_OnDevice )
        endif


        ! move current batch psib to host
        istat=cudaMemcpyAsync(psib(1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
             d_psibBatch(1,1,1,current), &
             QuadSet%Groups*Size%nbelem*anglebatch(current), transfer_stream ) 

        istat=cudaEventRecord(psib_OnHost( batch(current) ), transfer_stream )
        
        ! transfer stream waits for STime to be computed:
        istat = cudaStreamWaitEvent(transfer_stream, STimeFinished(batch(current)), 0)
        ! move current batch STime to host
        istat=cudaMemcpyAsync(Geom%ZDataSoA%STime(1,1,QuadSet%AngleOrder(mm1,binSend(current))), &
             d_STime(current)%data(1,1,1), &
             QuadSet%Groups*Size%ncornr*anglebatch(current), transfer_stream ) 



   enddo
   
   
   ! transfer stream waits for snmoments calc to be finished (the last one)
   istat = cudaStreamWaitEvent(transfer_stream1, snmomentsFinished(batch(current)), 0)
   
   ! move d_phi data to host:
   istat=cudaMemcpyAsync(phi(1,1), &
        d_phi(1,1), &
        QuadSet%Groups*Size%ncornr, transfer_stream1 )

   istat=cudaEventRecord( phi_OnHost, transfer_stream1 )
   ! Actually could just do a synchronize on event phi on host...
   !istat = cudaDeviceSynchronize()

   ! CPU code should wait until phi is on the host before using it
   istat=cudaEventSynchronize( phi_OnHost )


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



