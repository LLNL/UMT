!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   RTSTRTSN - Computes an angle-dependent source using the previous   *
!              time step angular intensity and saves copies of the     *
!              corner temperatures.                                    *
!                                                                      *
!   Input:   psir,tec                                                  *
!                                                                      *
!   Output:  ASRCC,TEZOLD,TECN                                         *
!                                                                      *
!***********************************************************************
 
   subroutine rtstrtsn( Phi )

   use kind_mod
   use constant_mod
   use Size_mod
   use Material_mod
   use Geometry_mod
   use ZoneData_mod
   use cudafor

   implicit none

!  Arguments

   !real(adqt), intent(in)    :: psir(Size%ngr,Size%ncornr,Size%nangSN)
   real(adqt), intent(in)    :: Phi(Size%ngr,Size%ncornr)

   !real(adqt), intent(inout) :: psib(Size%ngr,Size%nbelem,Size%nangSN)

!  Local Variables

   integer    :: ia, ic, ig, ncornr, nzones, ngr, nangSN
   integer    :: c, c0, nCorner, zone
   integer :: tid, nth, iabeg, iaend, istat
   integer, external :: omp_get_thread_num, omp_get_num_threads

   real(adqt) :: tau, mysum, quadwt

!  Mesh Constants

   ncornr = Size%ncornr
   nzones = Size%nzones
   ngr    = Size%ngr
   tau    = Size%tau
   nangSN = Size%nangSN

!  If this is a time-dependent problem compute the angle-dependent
!  source using the old time-step intensity

! This setting of STime could get moved into snswp3d_GPU and only done once for each time step
! psir will already be on the GPU in batches, and STime could be computed in batches and moved back 
! to the CPU. This means STime has to be brought back the first time and batched in after that.

! Alternative is to store entire STime on GPU--then you never pay to move it ever. If I keep this 
! here I'm just paying to move it once--which is pretty close to never.

!if (.not. allocated(Geom%ZDataSoA%STime)) then
!  print *, "Error STime should have been allocated already"
!  call abort
!endif

if (Size%itimsrc == 'exact') then
  ! Move tau*psir to the GPU STime array.
  !Geom%ZDataSoA%STime(:,:,:) = tau*psir(:,:,:)
  !Geom%ZDataSoA%STime = tau*psir
  !Geom%ZDataSoA%STime(:,:,:) = psir(:,:,:)
  ! Move psir into STime array for processing
!   istat=cudaMemcpy(Geom%ZDataSoA%STime(1,1,1),                 &
!                                    psir(1,1,1), &
!                                    int(ngr*ncornr*nangSN,kind=8))
!   print *, "istat = ",istat
!   istat = cudaDeviceSynchronize()

!   ! Multiply by tau to get STime
!   !$cuf kernel do(3) <<< (*,*), (16,16) >>>
   ! do ia=1,nangSN
   !    do ic=1,ncornr
   !       do ig=1, ngr
   !          Geom%ZDataSoA%STime(ig,ic,ia) = tau*psir(ig,ic,ia)
   !       enddo
   !    enddo
   ! enddo
else
!  Geom%ZDataSoA%STime = zero
endif



!    if (Size%itimsrc == 'exact') then
! !$omp parallel private(nCorner,c0,tid,nth,iabeg,iaend)
!      tid = omp_get_thread_num()
!      nth = omp_get_num_threads()
!      call omp_block_partition(tid,nth,1,Size%nangSN,iabeg,iaend)
!      ZoneLoop: do zone=1,nzones
!        Z => getZoneData(Geom, zone)
!        nCorner = Z% nCorner
!        c0      = Z% c0
!        do ia = iabeg, iaend  !YKT  was : ia=1,Size%nangSN
!          do c=1,nCorner
!            Z% STime(:,c,ia) = tau*psir(:,c0+c,ia)
!          enddo
!        enddo
!      enddo ZoneLoop
! !$omp end parallel
!    else
!      do zone=1,nzones
!        Z => getZoneData(Geom, zone)
!        Z% STime(:,:,:) = zero 
!      enddo
!    endif

   call timer_beg('_initialize')
!  Initialize arrays

   do ic=1,ncornr
     Mat%tecn(ic)  = Mat%tec(ic)
     Mat%denec(ic) = zero
   enddo

   call timer_end('_initialize')

!  Initialize the boundary flux array (PSIB)
   !call timer_beg('_setbdy')
   !call setbdy(psir, PSIB)
   !call timer_end('_setbdy')

!  Calculate zone-average energy density for convergence test
   call timer_beg('_zoneaverage')
!$omp parallel do private(zone,Z,nCorner,c0)
   do zone=1,nzones
     Z => getZoneData(Geom, zone)

     nCorner             = Z% nCorner
     c0                  = Z% c0
     Z% EnergyDensityOld = zero

     do c=1,nCorner
       mysum = zero
       do ig=1,ngr
         mysum = mysum + Phi(ig,c0+c)
       enddo
       Z% EnergyDensityOld = Z% EnergyDensityOld + Z%Volume(c)*mysum
     enddo
     Z% EnergyDensityOld = Z% EnergyDensityOld/Z%VolumeZone
   enddo
   call timer_end('_zoneaverage')

!  Compute zone-average temperature for convergence test

   call rtave(Mat%tec, Mat%TEZOLD)


   return
   end subroutine rtstrtsn
 
! block partition do-loop indices i1 to i2 
! for worker : myrank = 0, ..., nranks - 1
! using a method that is as balanced as possible

subroutine omp_block_partition(myrank, nranks, i1, i2, ibeg, iend)
  implicit none
  integer, intent(in) :: myrank, nranks, i1, i2
  integer, intent(out) :: ibeg, iend
  integer nwork, chunk, extra, ntcut  !local variables
  nwork = i2 - i1 + 1
  chunk = nwork/nranks
  extra = nwork - nranks*chunk
  ntcut = nranks - extra
  if (myrank .lt. ntcut) then
    ibeg = i1 + myrank*chunk
    iend = ibeg + chunk - 1
  else
    ibeg = i1 + ntcut*chunk + (myrank - ntcut)*(chunk + 1)
    iend = ibeg + chunk
  end if
  if (iend .gt. i2) iend = i2
end
