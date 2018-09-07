!***********************************************************************
!                        Version 1:  03/01, PFN                        *
!                                                                      *
!   SNSWP3D  - This routine calculates angular fluxes for a single     *
!              direction for an upstream corner-balance spatial        *
!              discretization in 3D.                                   *
!                                                                      *
!   Input:                                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!***********************************************************************

! number of zones that will be processed in parallel from each batch
#define NZONEPAR 48
! number of threads available for groups (must be >= groups)
#define THREADX 16
! the number of corners that can be processed at a time in shared memory
#define CORNERCHUNK 8

module snswp3d_mod
  use kind_mod
  use constant_mod
  use Quadrature_mod
  use ZoneData_mod
  use cudafor

contains






  attributes(global) subroutine GPU_sweep(anglebatch,                     &
                              nzones,               &
                              Groups,            &
                              ncornr,               &
                              NumAngles,         &
                              AngleOrder,        & ! only angle batch portion
                              maxCorner,            &
                              maxcf,                &
                              NangBin,                   &
                              nbelem,                &
                              ZData,      &
                              omega_A_fp,                &
                              omega_A_ez,                &
                              STimeBatch,               & ! only angle batch portion
                              psicbatch,                      &  ! only angle batch portion
                              psib,                      &
                              next,              &
                              nextZ,             &
                              passZstart             )
    implicit none

    !  Arguments

   integer, value,    intent(in)    :: anglebatch
   integer, value,    intent(in)    :: nzones
   integer, value,    intent(in)    :: Groups
   integer, value,    intent(in)    :: ncornr
   integer, value,    intent(in)    :: NumAngles

   integer,    device, intent(in) :: AngleOrder(anglebatch)

   integer, value,    intent(in)    :: maxCorner
   integer, value,    intent(in)    :: maxcf
   integer, value,    intent(in)    :: NangBin
   integer, value,    intent(in)    :: nbelem

   type(GPU_ZoneData), device, intent(in) :: ZData(nzones)

   real(adqt), device, intent(in)    :: omega_A_fp(maxcf,maxCorner, nzones, anglebatch) 
   real(adqt), device, intent(in)    :: omega_A_ez(maxcf,maxCorner,nzones, anglebatch) 

   real(adqt), device, intent(in)    :: STimeBatch(Groups, ncornr, anglebatch)
   real(adqt), device, intent(inout) :: psicbatch(Groups,ncornr,anglebatch)
   real(adqt), device, intent(inout) :: psib(Groups,nbelem,anglebatch)
   integer,    device, intent(in) :: next(ncornr+1,NumAngles)
   integer,    device, intent(in) :: nextZ(nzones,NumAngles)
   integer,    device, intent(in) :: passZstart(nzones,NumAngles)   


    !  Local Variables

    integer    :: Angle, i, ib, ic, icfp, icface, ifp, ig, k, nxez
    integer    :: zone, c, cez, ii, mm, ndone
    integer    :: p, ndoneZ, passZcount
    integer    :: nCorner, nCFaces, c0

    !!FIXME: sizes are hardcoded at present due to a CUDA Fortran limitation
    !real(adqt), shared :: Volume(8,NZONEPAR) ! (maxCorner)
    !real(adqt), shared :: A_fp(3,3,NZONEPAR) ! (ndim,maxcf)
    !real(adqt), shared :: A_ez(3,3,NZONEPAR) ! (ndim,maxcf)
    !integer,    shared :: Connect(3,3,NZONEPAR) ! (3,maxcf)

    !integer    :: ez_exit(3) ! (maxcf)

    real(adqt) :: fouralpha, fouralpha4, aez, aez2, area_opp, psi_opp
    real(adqt) :: source, sigv, sigv2, gnum, gtau, sez, sumArea
    real(adqt) :: Sigt!, SigtInv
    
    real(adqt) :: r_afpm, temp

    real(adqt) :: src(8)      ! (maxCorner)
    real(adqt) :: Q(8)        ! (maxCorner)
    !real(adqt) :: SigtVol(8)  ! (maxCorner)

    !real(adqt) :: afpm(3)     ! (maxcf)
    !real(adqt) :: coefpsic(3) ! (maxcf)

    real(adqt) :: psifp(3)    ! (maxcf)
    !real(adqt) :: tpsic(8)    ! (maxCorner)
    real(adqt) :: tpsic

    ! shared memory:
    real(adqt), shared :: coefpsic(maxcf,blockDim%z) ! (maxcf) 3*32 *8 bytes
    integer, shared    :: ez_exit(maxcf,blockDim%z) ! (maxcf) 3*32 *4 bytes
    real(adqt), shared :: sh_omega_A_fp(maxcf,maxCorner, blockDim%z)     ! (maxcf*maxCorner) 3*8*32 *8 bytes
    real(adqt), shared :: sh_omega_A_ez(maxcf,maxCorner, blockDim%z)     ! (maxcf*maxCorner) 3*8*32 *8 bytes
    integer, shared    :: Connect_ro(maxcf,maxCorner,3,blockDim%z) !(maxcf*maxCorner*ndim) 3*8*3*32 *4 bytes
    real(adqt), shared :: Volume(maxCorner, blockDim%z) !maxCorner 8*32 *8 bytes


    !  Constants

    parameter (fouralpha=1.82d0)
    parameter (fouralpha4=5.82d0)

    ! best if threadIdx is 16 or 32, and multiples of these groups are done by block y.
    ig = threadIdx%x + blockDim%x*(blockIdx%y-1)
    mm = blockIdx%x
    Angle = AngleOrder(mm)

    p = 0
    ndoneZ = 0
    PassLoop: do while (ndoneZ < nzones)
       p = p + 1
       passZcount = passZstart(p+1,Angle) - passZstart(p,Angle)
       
       !if(ig .le. Groups) then
          ZoneLoop: do ii=threadIdx%z, ((passZcount/blockDim%z)+1)*blockDim%z,blockDim%z

             call syncthreads

             if ( ii .le. passZcount ) then

             !!FIXME: simplifying assumption that all zones have same nCorner values
             !! (they're all 8 from what we've seen). If this isn't true in general,
             !! just convert this into a table lookup
             ndone = (ndoneZ+ii-1) * maxCorner

             zone = nextZ(ndoneZ+ii,Angle)

             nCorner =   ZData(zone)%nCorner
             nCFaces =   ZData(zone)%nCFaces
             c0      =   ZData(zone)%c0


             Sigt    =  ZData(zone)%Sigt(ig)
             !SigtInv = one/Sigt !need to thread?
             !SigtInv = SigtInvArray(ig,zone)


             ! coalesced load into volume shared memory
             do c = threadIdx%x, nCorner, blockDim%x
                Volume(c,threadIdx%z) = ZData(zone)%Volume(c)
             enddo

             ! could loop over chunks of ncfaces and nCorner that fit with threadidx.
             ! Could use more shared memory, or could work on set of ncfaces*nCorner at a time?
             !if(threadIdx%x <= ncfaces*nCorner) then
             do k = threadIdx%x, ncfaces*nCorner, blockDim%x

                c = (k-1)/maxcf + 1 ! split thread block x-dimension loop into two dims
                icface = k - ((c-1)*maxcf) ! remainder
                ! this gives c = 1,1,1, 2,2,2, 3,3,3
                !       icface = 1,2,3, 1,2,3, 1,2,3

                ! coalesced loads into shared memory array
                sh_omega_A_fp(icface,c,threadIdx%z) = omega_A_fp(icface,c,zone,mm)
                sh_omega_A_ez(icface,c,threadIdx%z) = omega_A_ez(icface,c,zone,mm)                

                !do i=1,3 !ndim
                Connect_ro(icface,c,1,threadIdx%z) = ZData(zone)%connect_reorder(icface,c,1)
                Connect_ro(icface,c,2,threadIdx%z) = ZData(zone)%connect_reorder(icface,c,2)
                Connect_ro(icface,c,3,threadIdx%z) = ZData(zone)%connect_reorder(icface,c,3)
                !enddo

             enddo

             endif

             call syncthreads

             if ( ii .le. passZcount ) then

             !  Contributions from volume terms

             do c=1,nCorner
                source     =  ZData(zone)%STotal(ig,c) +  STimeBatch(ig,c0+c,mm)
                !Q(c)       = SigtInv*source 
                Q(c)       = source/Sigt 
                src(c)     =  Volume(c,threadIdx%z)*source
                !SigtVol(c) = Sigt*Volume(c,zone)
             enddo


             CornerLoop: do i=1,nCorner

                ic      = next(ndone+i,Angle)
                c       = ic - c0

                !  Calculate Area_CornerFace dot Omega to determine the 
                !  contributions from incident fluxes across external 
                !  corner faces (FP faces)

                sumArea = zero


                do icface=1,ncfaces

                   ! afpm(icface) = omega(1,Angle)* A_fp(1,icface,c,zone) + &
                   !      omega(2,Angle)* A_fp(2,icface,c,zone) + &
                   !      omega(3,Angle)* A_fp(3,icface,c,zone)
                   
                   !afpm(icface) = omega_A_fp(icface,c,zone,mm)
                   r_afpm = sh_omega_A_fp(icface,c,threadIdx%z)
                   
                   !icfp    =  Connect(1,icface,c,zone)
                   !ib      =  Connect(2,icface,c,zone)

                   icfp    =  Connect_ro(icface,c,1,threadIdx%z)
                   ib      =  Connect_ro(icface,c,2,threadIdx%z)

                   if ( r_afpm >= zero ) then
                      sumArea = sumArea + r_afpm
                   else
                      if (icfp == 0) then
                         psifp(icface) = psib(ig,ib,mm)
                      else
                         psifp(icface) = psicbatch(ig,icfp,mm)
                      endif

                      src(c) = src(c) - r_afpm*psifp(icface)
                   endif
                enddo

                !  Contributions from interior corner faces (EZ faces)

                nxez = 0

                do icface=1,nCFaces

                   ! aez = omega(1,Angle)* A_ez(1,icface,c,zone) + &
                   !      omega(2,Angle)* A_ez(2,icface,c,zone) + &
                   !      omega(3,Angle)* A_ez(3,icface,c,zone) 
                   
                   aez = sh_omega_A_ez(icface,c,threadIdx%z)
                   
                   if (aez > zero ) then

                      sumArea        = sumArea + aez
                      area_opp       = zero
                      nxez           = nxez + 1
                      !cez            = Connect(3,icface,c,zone)
                      cez            = Connect_ro(icface,c,3,threadIdx%z)
                      ez_exit(nxez,threadIdx%z)  = cez
                      coefpsic(nxez,threadIdx%z) = aez

                      if (nCFaces == 3) then

                         ifp = mod(icface,nCFaces) + 1                         

                         ! need to do r_afpm = r_afpm(ifp)
                         r_afpm = sh_omega_A_fp(ifp,c,threadIdx%z)

                         if ( r_afpm < zero ) then
                            area_opp   = -r_afpm
                            psi_opp    =  psifp(ifp)
                         endif

                      else

                         ifp        = icface
                         area_opp   = zero
                         psi_opp    = zero

                         do k=1,nCFaces-2
                            ifp = mod(ifp,nCFaces) + 1
                            r_afpm = sh_omega_A_fp(ifp,c,threadIdx%z)
                            if ( r_afpm < zero ) then
                               area_opp   = area_opp   - r_afpm
                               psi_opp    = psi_opp    - r_afpm*psifp(ifp)
                            endif
                         enddo

                         psi_opp = psi_opp/area_opp

                      endif

                      temp = half*aez*( Q(c) - Q(cez) )
                      
                      TestOppositeFace: if (area_opp > zero) then

                         aez2 = aez*aez
                         !sigv    = SigtVol(c)
                         sigv    = Sigt*Volume(c,threadIdx%z)
                         sigv2        = sigv*sigv
                         gnum         = aez2*( fouralpha*sigv2 +              &
                              aez*(four*sigv + three*aez) )

                         gtau         = gnum/                                    &
                              ( gnum + four*sigv2*sigv2 + aez*sigv*(six*sigv2 + &
                              two*aez*(two*sigv + aez)) ) 

                         sez          = gtau*sigv*( psi_opp - Q(c) ) +   &
                              (one - gtau)*temp

                      else

                         sez          = temp

                      endif TestOppositeFace

                      src(c)       = src(c)   + sez
                      src(cez)     = src(cez) - sez


                   endif

                enddo

                !  Corner angular flux

                tpsic = src(c)/(sumArea + Sigt*Volume(c,threadIdx%z) )

                !  Calculate the angular flux exiting all "FP" surfaces
                !  and the current exiting all "EZ" surfaces.
                !  The downstream corner index is "ez_exit."

                !  Zone Interior or "EZ" Faces

                do icface=1,nxez
                   cez      = ez_exit(icface,threadIdx%z)
                   !r_coefpsic = coefpsic(icface,thread
                   src(cez) = src(cez) + coefpsic(icface,threadIdx%z)*tpsic
                enddo

                psicbatch(ig,c0+c,mm) = tpsic

             enddo CornerLoop

             !  Copy temporary flux into the global array

             !print *, "ig, c0, Angle", ig, c0, Angle
             !do c=1,nCorner

                !!!!psicbatch(ig,c0+c,mm) = tpsic(c)

                !if(ig>Groups .or. c0+c > ncornr .or. Angle > NumAngles) then
                !   print *, "ig, c0, c, Angle", ig, c0, c, Angle
                !endif
                !psic(ig,c0+c,Angle) = tpsic(c)
                !psic(ig,c0+c,Angle) = tpsic(c)
             !enddo
             
             endif

          enddo ZoneLoop
       !endif ! ig .le. groups

       ndoneZ = ndoneZ + passZcount

       call syncthreads

    enddo PassLoop

  end subroutine GPU_sweep






  attributes(global) subroutine setExitFluxD(  anglebatch, Angles, psicache, psibBatch, d_iExit, groups,ncornr, nbelem)

    use kind_mod
    use constant_mod
    use Size_mod
    use Geometry_mod
    use Quadrature_mod

    implicit none

    integer, value,  intent(in)    :: anglebatch, groups, ncornr, nbelem

    integer, device, intent(in)    :: Angles(anglebatch)

    real(adqt), device, intent(in)  :: psicache(groups,ncornr,anglebatch) 

    ! ALERT--SHOULDN'T THIS BE DEVICE ARRAY:
    real(adqt), device, intent(out) :: psibBatch(Groups,nbelem,anglebatch)

    type(GPUExit), device, intent(in) :: d_iExit(:) 

    ! Local variables

    integer :: mm, Angle, ig, i, ib, ic


    !  Set exiting boundary fluxes

    do mm=blockIdx%x, anglebatch, gridDim%x
       do ig=threadIdx%x, Groups, blockDim%x
          Angle = Angles(mm)
          do i=1,d_iExit(Angle)% nExit
             ib = d_iExit(Angle)% d_ListExit(1,i)
             ic = d_iExit(Angle)% d_ListExit(2,i)

             psibBatch(ig,ib,mm) = psicache(ig,ic,mm)
          enddo
       enddo

    enddo

  end subroutine setExitFluxD



subroutine setExitFluxD2(  anglebatch, Angles, psicache, psibBatch, d_iExit, groups,ncornr, nbelem, streamid)
   
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use Quadrature_mod
   
    implicit none

   integer,  intent(in)    :: anglebatch, groups, ncornr, nbelem

   integer, device, intent(in)    :: Angles(anglebatch)

   real(adqt), device, intent(in)  :: psicache(groups,ncornr,anglebatch) 

   real(adqt), device, intent(out) :: psibBatch(Groups,nbelem,anglebatch)

   type(GPUExit), device, intent(in) :: d_iExit(:) 

   integer(kind=cuda_stream_kind), intent(in) :: streamid

   ! Local variables

   integer :: mm, Angle, ig, i, ib, ic
   

   !  Set exiting boundary fluxes
   !$cuf kernel do(2) <<< (*,*), (*,*), stream=streamid >>>            
   do mm=1, anglebatch
      do ig=1, Groups
         Angle = Angles(mm)
         do i=1,d_iExit(Angle)% nExit
            ib = d_iExit(Angle)% d_ListExit(1,i)
            ic = d_iExit(Angle)% d_ListExit(2,i)
            
            psibBatch(ig,ib,mm) = psicache(ig,ic,mm)
         enddo
      enddo
      
   enddo

end subroutine setExitFluxD2






  subroutine setExitFlux(  anglebatch, Angles, psic, psib)
   
    use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use Quadrature_mod
   
    implicit none

   integer,   intent(in)    :: anglebatch

   integer,    intent(in) :: Angles(anglebatch)

   !real(adqt), device, intent(in)  :: psicache(QuadSet%Groups,Size%ncornr,anglebatch) 
   real(adqt), intent(in) :: psic(QuadSet%Groups,Size%ncornr,QuadSet%NumAngles)

   real(adqt), intent(out) :: psib(QuadSet%Groups,Size%nbelem,QuadSet%NumAngles)

   integer :: mm, Angle, i, ib, ic


!  Set exiting boundary fluxes

!$OMP parallel do private(Angle,mm,ExitBdy,i,ib,ic)
  do mm=1,anglebatch
   Angle = Angles(mm)

   ExitBdy => getExitList(QuadSet, Angle)
   !iExit => QuadSet%iExit(Angle)

   do i=1,ExitBdy% nExit
     ib = ExitBdy% ListExit(1,i)
     ic = ExitBdy% ListExit(2,i)
     !do ig=1,Groups
        psib(:,ib,Angle) = psic(:,ic,Angle)
     !enddo
   enddo

  enddo

end subroutine setExitFlux






 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! snswp3d Caller
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine snswp3d_f(     anglebatch,                     &
                              nzones,               &
                              Groups,            &
                              ncornr,               &
                              NumAngles,         &
                              d_AngleOrder,        & ! only angle batch portion
                              maxCorner,            &
                              maxcf,                &
                              NangBin,                   &
                              nbelem,                &
                              d_omega_A_fp,                &
                              d_omega_A_ez,                &
                              d_STimeBatch,               & ! only angle batch portion
                              d_psicbatch,                      &  ! only angle batch portion
                              d_psib,                      &
                              d_next,              &
                              d_nextZ,             &
                              d_passZstart,              &
                              streamid)

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
   integer,    intent(in)    :: Groups
   integer,    intent(in)    :: ncornr
   integer,    intent(in)    :: NumAngles

   integer,    device, intent(in) :: d_AngleOrder(anglebatch)

   integer,    intent(in)    :: maxCorner
   integer,    intent(in)    :: maxcf
   integer,    intent(in)    :: NangBin
   integer,    intent(in)    :: nbelem

   real(adqt), device, intent(in)    :: d_omega_A_fp(maxcf ,maxCorner, nzones, anglebatch) 
   real(adqt), device, intent(in)    :: d_omega_A_ez(maxcf ,maxCorner, nzones, anglebatch)  


   real(adqt), device, intent(in)    :: d_STimeBatch(Groups, ncornr, anglebatch)
   real(adqt), device, intent(inout) :: d_psicbatch(QuadSet%Groups,Size%ncornr,anglebatch)
   real(adqt), device, intent(inout) :: d_psib(QuadSet%Groups,Size%nbelem,anglebatch)
   integer,    device, intent(in) :: d_next(Size%ncornr+1,QuadSet%NumAngles)
   integer,    device, intent(in) :: d_nextZ(Size%nzones,QuadSet%NumAngles)
   integer,    device, intent(in) :: d_passZstart(Size%nzones,QuadSet%NumAngles)   
   integer(kind=cuda_stream_kind), intent(in) :: streamid   

!  Local Variables

   integer    :: mm, Angle,istat,i,ib,ic, groupset

   type(dim3) :: threads,blocks
   
   integer    :: shmem !amount of shared memory need by GPU sweep kernel.

   ! Groups are done in chunks of THREADX size, and groupset is number of those chunks.
   groupset = QuadSet%Groups/THREADX

   ! THREADX*NZONEPAR must be .le. 1024 on P100 hardware
   !threads=dim3(QuadSet%Groups,NZONEPAR,1) 
   threads=dim3(THREADX,1,NZONEPAR) 
   blocks=dim3(anglebatch,groupset,1)

   ! shared memory needs:
   shmem = &
        maxcf*NZONEPAR*8 + & ! coefpsic
        maxcf*NZONEPAR*4 + & ! ez_exit
        maxcf*maxCorner*NZONEPAR*8 + & ! afpm with cf and c
        maxcf*maxCorner*NZONEPAR*8 + & ! aez with cf and c
        maxcf*maxCorner*3*NZONEPAR*4 + & ! connect_ro
        maxCorner*NZONEPAR*8           ! volume

   !print *, "shared memory size = ", real(shmem)/real(1024), " kb"
   !print *, "blocks = ", blocks
   !print *, "threads = ", threads

   call GPU_sweep<<<blocks,threads,shmem,streamid>>>( anglebatch,                     &
                              nzones,               &
                              Groups,            &
                              ncornr,               &
                              NumAngles,         &
                              d_AngleOrder,        & ! only angle batch portion
                              maxCorner,            &
                              maxcf,                &
                              NangBin,                   &
                              nbelem,                &
                              Geom%d_GPU_ZData,    &
                              d_omega_A_fp,                &
                              d_omega_A_ez,                &
                              d_STimeBatch,               &  ! only angle batch portion
                              d_psicbatch,                      &  ! only angle batch portion
                              d_psib,                      &
                              d_next,              &
                              d_nextZ,             &
                              d_passZstart             )



   return
   end subroutine snswp3d_f






 end module snswp3d_mod
