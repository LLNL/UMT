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
#define NZONEPAR 5
! number of threads available for groups (must be >= groups)
#define THREADX 192

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
                              Angles,        & ! only angle batch portion
                              maxCorner,            &
                              maxcf,                &
                              binRecv,                   &
                              NangBin,                   &
                              nbelem,                &
                              omega,             &
                              nCornerArray,                &
                              nCFacesArray,                &
                              c0Array,                &
                              A_fp,                &
                              A_ez,                &
                              Connect,             &
                              STotal,              &
                              STimeBatch,               & ! only angle batch portion
                              Volume,             &
                              psicbatch,                      &  ! only angle batch portion
                              psib,                      &
                              next,              &
                              nextZ,             &
                              SigtArray,                &
 !                             SigtInv,             &
                              passZstart             )
    implicit none

    !  Arguments

   integer, value,    intent(in)    :: anglebatch
   integer, value,    intent(in)    :: nzones
   integer, value,    intent(in)    :: Groups
   integer, value,    intent(in)    :: ncornr
   integer, value,    intent(in)    :: NumAngles

   integer,    device, intent(in) :: Angles(anglebatch)

   integer, value,    intent(in)    :: maxCorner
   integer, value,    intent(in)    :: maxcf
   integer, value,    intent(in)    :: binRecv
   integer, value,    intent(in)    :: NangBin
   integer, value,    intent(in)    :: nbelem
   real(adqt), device, intent(in)    :: omega(3,NumAngles)
   integer, device, intent(in):: nCornerArray(nzones)
   integer, device, intent(in):: nCFacesArray(nzones)
   integer, device, intent(in):: c0Array(nzones)
   real(adqt), device, intent(in)    :: A_fp(3,maxcf,maxCorner,nzones) !ndim,maxcf,maxCorner,nzones
   real(adqt), device, intent(in)    :: A_ez(3,maxcf,maxCorner,nzones) !ndim,maxcf,maxCorner,nzones
   integer, device, intent(in):: Connect(3,maxcf,maxCorner,nzones) 
   real(adqt), device, intent(in)    :: STotal(Groups, maxCorner, nzones)
   real(adqt), device, intent(in)    :: STimeBatch(Groups, ncornr, anglebatch)
   real(adqt), device, intent(in)    :: Volume(maxCorner, nzones)
   real(adqt), device, intent(inout) :: psicbatch(Groups,ncornr,anglebatch)
   real(adqt), device, intent(inout) :: psib(Groups,nbelem,NumAngles)
   integer,    device, intent(in) :: next(ncornr+1,NumAngles)
   integer,    device, intent(in) :: nextZ(nzones,NumAngles)
   real(adqt), device, intent(in)    :: SigtArray(Groups, nzones)
!   real(adqt), device, intent(in)    :: SigtInv(Groups, nzones)
   integer,    device, intent(in) :: passZstart(nzones,NumAngles)   


    !  Local Variables

    integer    :: Angle, i, ib, ic, icfp, icface, id, ifp, ig, k, nxez
    integer    :: zone, c, cez, ii, mm, ndone
    integer    :: p, ndoneZ, passZcount
    integer    :: nCorner, nCFaces, c0

    !!FIXME: sizes are hardcoded at present due to a CUDA Fortran limitation
    !real(adqt), shared :: Volume(8,NZONEPAR) ! (maxCorner)
    !real(adqt), shared :: A_fp(3,3,NZONEPAR) ! (ndim,maxcf)
    !real(adqt), shared :: A_ez(3,3,NZONEPAR) ! (ndim,maxcf)
    !integer,    shared :: Connect(3,3,NZONEPAR) ! (3,maxcf)

    integer    :: ez_exit(3) ! (maxcf)

    real(adqt) :: fouralpha, fouralpha4, aez, aez2, area_opp, psi_opp
    real(adqt) :: source, sigv, sigv2, gnum, gtau, sez, sumArea
    real(adqt) :: Sigt, SigtInv

    real(adqt) :: src(8)      ! (maxCorner)
    real(adqt) :: Q(8)        ! (maxCorner)
    real(adqt) :: afpm(3)     ! (maxcf)
    real(adqt) :: coefpsic(3) ! (maxcf)
    real(adqt) :: psifp(3)    ! (maxCorner)
    real(adqt) :: tpsic(8)    ! (maxCorner)

    !  Constants

    parameter (fouralpha=1.82d0)
    parameter (fouralpha4=5.82d0)


    ig = threadIdx%x
    mm = blockIdx%z
    Angle = Angles(mm)

    p = 0
    ndoneZ = 0
    PassLoop: do while (ndoneZ < nzones)
       p = p + 1
       passZcount = passZstart(p+1,Angle) - passZstart(p,Angle)
       
       !if(ig .le. Groups) then
          ZoneLoop: do ii=threadIdx%z,passZcount,blockDim%z

             !!FIXME: simplifying assumption that all zones have same nCorner values
             !! (they're all 8 from what we've seen). If this isn't true in general,
             !! just convert this into a table lookup
             ndone = (ndoneZ+ii-1) * maxCorner

             zone = nextZ(ndoneZ+ii,Angle)

             nCorner = nCornerArray(zone)
             nCFaces =   nCFacesArray(zone)
             c0      =   c0Array(zone)
             Sigt    =  SigtArray(ig,zone)
             SigtInv = one/Sigt !need to thread?

             !  Contributions from volume terms

             do c=1,nCorner
                !do ig= threadIdx%x, Groups, blockDim%x
                source     =  STotal(ig,c,zone) +  STimeBatch(ig,c0+c,mm)
                !enddo
                Q(c)       = SigtInv*source 
                src(c)     =  Volume(c,zone)*source
             enddo

             CornerLoop: do i=1,nCorner

                ic      = next(ndone+i,Angle)
                c       = ic - c0

                sigv    = Sigt* Volume(c,zone)

                !  Calculate Area_CornerFace dot Omega to determine the 
                !  contributions from incident fluxes across external 
                !  corner faces (FP faces)

                sumArea = zero

                do icface=1,ncfaces

                   afpm(icface) = omega(1,Angle)* A_fp(1,icface,c,zone) + &
                        omega(2,Angle)* A_fp(2,icface,c,zone) + &
                        omega(3,Angle)* A_fp(3,icface,c,zone)

                   icfp    =  Connect(1,icface,c,zone)
                   ib      =  Connect(2,icface,c,zone)

                   if ( afpm(icface) >= zero ) then
                      sumArea = sumArea + afpm(icface)
                   else
                      if (icfp == 0) then
                         psifp(icface) = psib(ig,ib,Angle)
                      else
                         psifp(icface) = psicbatch(ig,icfp,mm)
                      endif

                      src(c) = src(c) - afpm(icface)*psifp(icface)
                   endif
                enddo

                !  Contributions from interior corner faces (EZ faces)

                nxez = 0

                do icface=1,nCFaces

                   aez = omega(1,Angle)* A_ez(1,icface,c,zone) + &
                        omega(2,Angle)* A_ez(2,icface,c,zone) + &
                        omega(3,Angle)* A_ez(3,icface,c,zone) 

                   if (aez > zero ) then

                      sumArea        = sumArea + aez
                      area_opp       = zero
                      nxez           = nxez + 1
                      cez            = Connect(3,icface,c,zone)
                      ez_exit(nxez)  = cez
                      coefpsic(nxez) = aez

                      if (nCFaces == 3) then

                         ifp = mod(icface,nCFaces) + 1

                         if ( afpm(ifp) < zero ) then
                            area_opp   = -afpm(ifp)
                            psi_opp    =  psifp(ifp)
                         endif

                      else

                         ifp        = icface
                         area_opp   = zero
                         psi_opp    = zero

                         do k=1,nCFaces-2
                            ifp = mod(ifp,nCFaces) + 1
                            if ( afpm(ifp) < zero ) then
                               area_opp   = area_opp   - afpm(ifp)
                               psi_opp    = psi_opp    - afpm(ifp)*psifp(ifp)
                            endif
                         enddo

                         psi_opp = psi_opp/area_opp

                      endif

                      TestOppositeFace: if (area_opp > zero) then

                         aez2 = aez*aez

                         sigv2        = sigv*sigv
                         gnum         = aez2*( fouralpha*sigv2 +              &
                              aez*(four*sigv + three*aez) )

                         gtau         = gnum/                                    &
                              ( gnum + four*sigv2*sigv2 + aez*sigv*(six*sigv2 + &
                              two*aez*(two*sigv + aez)) ) 

                         sez          = gtau*sigv*( psi_opp - Q(c) ) +   &
                              half*aez*(one - gtau)*( Q(c) - Q(cez) )
                         src(c)       = src(c)   + sez
                         src(cez)     = src(cez) - sez

                      else

                         sez          = half*aez*( Q(c) - Q(cez) )
                         src(c)       = src(c)   + sez
                         src(cez)     = src(cez) - sez

                      endif TestOppositeFace

                   endif

                enddo

                !  Corner angular flux

                tpsic(c) = src(c)/(sumArea + sigv)

                !  Calculate the angular flux exiting all "FP" surfaces
                !  and the current exiting all "EZ" surfaces.
                !  The downstream corner index is "ez_exit."

                !  Zone Interior or "EZ" Faces

                do icface=1,nxez
                   cez      = ez_exit(icface)
                   src(cez) = src(cez) + coefpsic(icface)*tpsic(c)
                enddo

             enddo CornerLoop

             !  Copy temporary flux into the global array

             !print *, "ig, c0, Angle", ig, c0, Angle
             do c=1,nCorner
                psicbatch(ig,c0+c,mm) = tpsic(c)
                !if(ig>Groups .or. c0+c > ncornr .or. Angle > NumAngles) then
                !   print *, "ig, c0, c, Angle", ig, c0, c, Angle
                !endif
                !psic(ig,c0+c,Angle) = tpsic(c)
                !psic(ig,c0+c,Angle) = tpsic(c)
             enddo

          enddo ZoneLoop
       !endif ! ig .le. groups

       ndoneZ = ndoneZ + passZcount

       call syncthreads

    enddo PassLoop

  end subroutine GPU_sweep






  attributes(global) subroutine GPU_sweep_old(Groups, NumAngles, anglebatch, &
       ncornr, nzones, nbelem, &
       maxcf, maxCorner, passZstart, Angles, omega, &
       ZDataSoA, next, nextZ, psicbatch, psib)
    implicit none

    !  Arguments

    integer,    value, intent(in)     :: Groups, NumAngles, anglebatch, ncornr, &
         nzones, nbelem, maxcf, maxCorner

    integer,    device, intent(in)    :: passZstart(nzones,NumAngles)
    integer,    device, intent(in)    :: Angles(anglebatch)
    real(adqt), device, intent(in)    :: omega(3,NumAngles)

    !type(ZoneData),     device, intent(in) :: ZData(nzones)
    type(ZoneData_SoA), device, intent(in) :: ZDataSoA

    integer,    device, intent(in)    :: next(ncornr+1,NumAngles)
    integer,    device, intent(in)    :: nextZ(nzones,NumAngles)

    real(adqt), device, intent(inout) :: psicbatch(Groups,ncornr,anglebatch)
    real(adqt), device, intent(inout) :: psib(Groups,nbelem,NumAngles)

    !  Local Variables

    integer    :: Angle, i, ib, ic, icfp, icface, id, ifp, ig, k, nxez
    integer    :: zone, c, cez, ii, mm, ndone
    integer    :: p, ndoneZ, passZcount
    integer    :: nCorner, nCFaces, c0

    !!FIXME: sizes are hardcoded at present due to a CUDA Fortran limitation
    real(adqt), shared :: Volume(8,NZONEPAR) ! (maxCorner)
    real(adqt), shared :: A_fp(3,3,NZONEPAR) ! (ndim,maxcf)
    real(adqt), shared :: A_ez(3,3,NZONEPAR) ! (ndim,maxcf)
    integer,    shared :: Connect(3,3,NZONEPAR) ! (3,maxcf)

    integer    :: ez_exit(3) ! (maxcf)

    real(adqt) :: fouralpha, fouralpha4, aez, aez2, area_opp, psi_opp
    real(adqt) :: source, sigv, sigv2, gnum, gtau, sez, sumArea
    real(adqt) :: Sigt, SigtInv

    real(adqt) :: src(8)      ! (maxCorner)
    real(adqt) :: Q(8)        ! (maxCorner)
    real(adqt) :: afpm(3)     ! (maxcf)
    real(adqt) :: coefpsic(3) ! (maxcf)
    real(adqt) :: psifp(3)    ! (maxCorner)
    real(adqt) :: tpsic(8)    ! (maxCorner)

    !  Constants

    parameter (fouralpha=1.82d0)
    parameter (fouralpha4=5.82d0)


    ig = threadIdx%x
    mm = blockIdx%z
    Angle = Angles(mm)

    p = 0
    ndoneZ = 0
    PassLoop: do while (ndoneZ < nzones)
       p = p + 1
       passZcount = passZstart(p+1,Angle) - passZstart(p,Angle)
       
       !if(ig .le. Groups) then
          ZoneLoop: do ii=threadIdx%z,passZcount,blockDim%z

             !!FIXME: simplifying assumption that all zones have same nCorner values
             !! (they're all 8 from what we've seen). If this isn't true in general,
             !! just convert this into a table lookup
             ndone = (ndoneZ+ii-1) * maxCorner

             zone = nextZ(ndoneZ+ii,Angle)

             nCorner = ZDataSoA% nCorner(zone)
             nCFaces = ZDataSoA% nCFaces(zone)
             c0      = ZDataSoA% c0(zone)
             Sigt    = ZDataSoA%Sigt(ig,zone)
             SigtInv = one/Sigt !need to thread?

             !  Contributions from volume terms

             do c=1,nCorner
                !do ig= threadIdx%x, Groups, blockDim%x
                source     = ZDataSoA%STotal(ig,c,zone) + ZDataSoA%STime(ig,c0+c,Angle)
                !enddo
                Q(c)       = SigtInv*source 
                src(c)     = ZDataSoA%Volume(c,zone)*source
             enddo

             CornerLoop: do i=1,nCorner

                ic      = next(ndone+i,Angle)
                c       = ic - c0

                sigv    = Sigt*ZDataSoA%Volume(c,zone)

                !  Calculate Area_CornerFace dot Omega to determine the 
                !  contributions from incident fluxes across external 
                !  corner faces (FP faces)

                sumArea = zero

                do icface=1,ncfaces

                   afpm(icface) = omega(1,Angle)*ZDataSoA%A_fp(1,icface,c,zone) + &
                        omega(2,Angle)*ZDataSoA%A_fp(2,icface,c,zone) + &
                        omega(3,Angle)*ZDataSoA%A_fp(3,icface,c,zone)

                   icfp    = ZDataSoA%Connect(1,icface,c,zone)
                   ib      = ZDataSoA%Connect(2,icface,c,zone)

                   if ( afpm(icface) >= zero ) then
                      sumArea = sumArea + afpm(icface)
                   else
                      if (icfp == 0) then
                         psifp(icface) = psib(ig,ib,Angle)
                      else
                         psifp(icface) = psicbatch(ig,icfp,mm)
                      endif

                      src(c) = src(c) - afpm(icface)*psifp(icface)
                   endif
                enddo

                !  Contributions from interior corner faces (EZ faces)

                nxez = 0

                do icface=1,nCFaces

                   aez = omega(1,Angle)*ZDataSoA%A_ez(1,icface,c,zone) + &
                        omega(2,Angle)*ZDataSoA%A_ez(2,icface,c,zone) + &
                        omega(3,Angle)*ZDataSoA%A_ez(3,icface,c,zone) 

                   if (aez > zero ) then

                      sumArea        = sumArea + aez
                      area_opp       = zero
                      nxez           = nxez + 1
                      cez            = ZDataSoA%Connect(3,icface,c,zone)
                      ez_exit(nxez)  = cez
                      coefpsic(nxez) = aez

                      if (nCFaces == 3) then

                         ifp = mod(icface,nCFaces) + 1

                         if ( afpm(ifp) < zero ) then
                            area_opp   = -afpm(ifp)
                            psi_opp    =  psifp(ifp)
                         endif

                      else

                         ifp        = icface
                         area_opp   = zero
                         psi_opp    = zero

                         do k=1,nCFaces-2
                            ifp = mod(ifp,nCFaces) + 1
                            if ( afpm(ifp) < zero ) then
                               area_opp   = area_opp   - afpm(ifp)
                               psi_opp    = psi_opp    - afpm(ifp)*psifp(ifp)
                            endif
                         enddo

                         psi_opp = psi_opp/area_opp

                      endif

                      TestOppositeFace: if (area_opp > zero) then

                         aez2 = aez*aez

                         sigv2        = sigv*sigv
                         gnum         = aez2*( fouralpha*sigv2 +              &
                              aez*(four*sigv + three*aez) )

                         gtau         = gnum/                                    &
                              ( gnum + four*sigv2*sigv2 + aez*sigv*(six*sigv2 + &
                              two*aez*(two*sigv + aez)) ) 

                         sez          = gtau*sigv*( psi_opp - Q(c) ) +   &
                              half*aez*(one - gtau)*( Q(c) - Q(cez) )
                         src(c)       = src(c)   + sez
                         src(cez)     = src(cez) - sez

                      else

                         sez          = half*aez*( Q(c) - Q(cez) )
                         src(c)       = src(c)   + sez
                         src(cez)     = src(cez) - sez

                      endif TestOppositeFace

                   endif

                enddo

                !  Corner angular flux

                tpsic(c) = src(c)/(sumArea + sigv)

                !  Calculate the angular flux exiting all "FP" surfaces
                !  and the current exiting all "EZ" surfaces.
                !  The downstream corner index is "ez_exit."

                !  Zone Interior or "EZ" Faces

                do icface=1,nxez
                   cez      = ez_exit(icface)
                   src(cez) = src(cez) + coefpsic(icface)*tpsic(c)
                enddo

             enddo CornerLoop

             !  Copy temporary flux into the global array

             !print *, "ig, c0, Angle", ig, c0, Angle
             do c=1,nCorner
                psicbatch(ig,c0+c,mm) = tpsic(c)
                !if(ig>Groups .or. c0+c > ncornr .or. Angle > NumAngles) then
                !   print *, "ig, c0, c, Angle", ig, c0, c, Angle
                !endif
                !psic(ig,c0+c,Angle) = tpsic(c)
                !psic(ig,c0+c,Angle) = tpsic(c)
             enddo

          enddo ZoneLoop
       !endif ! ig .le. groups

       ndoneZ = ndoneZ + passZcount

       call syncthreads

    enddo PassLoop

  end subroutine GPU_sweep_old





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

   real(adqt), intent(out) :: psibBatch(Groups,nbelem,anglebatch)

   type(Exit), device, target, intent(in) :: d_iExit(:) 

   ! Local variables

   integer :: mm, Angle, ig, i, ib, ic

   !type(Exit), device, pointer :: d_ExitBdy
   


!  Set exiting boundary fluxes

  do mm=blockIdx%x, anglebatch, gridDim%x
   Angle = Angles(mm)

   !d_iExit(Angle)

     do ig=threadIdx%x, Groups, blockDim%x
        do i=1,d_iExit(Angle)% nExit
           ib = d_iExit(Angle)% d_ListExit(1,i)
           ic = d_iExit(Angle)% d_ListExit(2,i)

           psibBatch(ig,ib,mm) = psicache(ig,ic,mm)
        enddo
     enddo

  enddo

  ! do mm=blockIdx%x, anglebatch, gridDim%x
  !  Angle = Angles(mm)

  !  d_ExitBdy => d_iExit(Angle)

  !    do ig=threadIdx%x, Groups, blockDim%x
  !       do i=1,d_ExitBdy% nExit
  !          ib = d_ExitBdy% d_ListExit(1,i)
  !          ic = d_ExitBdy% d_ListExit(2,i)

  !          psibBatch(ig,ib,mm) = psicache(ig,ic,mm)
  !       enddo
  !    enddo

  ! enddo

end subroutine setExitFluxD






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
! Caller
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine snswp3d(     anglebatch,                     &
                              nzones,               &
                              Groups,            &
                              ncornr,               &
                              NumAngles,         &
                              d_Angles,        & ! only angle batch portion
                              maxCorner,            &
                              maxcf,                &
                              binRecv,                   &
                              NangBin,                   &
                              nbelem,                &
                              d_omega,             &
                              d_nCorner,                &
                              d_nCFaces,                &
                              d_c0,                &
                              d_A_fp,                &
                              d_A_ez,                &
                              d_Connect,             &
                              d_STotal,              &
                              d_STimeBatch,               & ! only angle batch portion
                              d_Volume,             &
                              d_psicbatch,                      &  ! only angle batch portion
                              d_psib,                      &
                              d_next,              &
                              d_nextZ,             &
                              d_Sigt,                &
                              d_SigtInv,             &
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

   integer,    device, intent(in) :: d_Angles(anglebatch)

   integer,    intent(in)    :: maxCorner
   integer,    intent(in)    :: maxcf
   integer,    intent(in)    :: binRecv
   integer,    intent(in)    :: NangBin
   integer,    intent(in)    :: nbelem
   real(adqt), device, intent(in)    :: d_omega(3,NumAngles)
   integer, device, intent(in):: d_nCorner(nzones)
   integer, device, intent(in):: d_nCFaces(nzones)
   integer, device, intent(in):: d_c0(nzones)
   real(adqt), device, intent(in)    :: d_A_fp(ndim,maxcf,maxCorner,nzones) !ndim,maxcf,maxCorner,nzones
   real(adqt), device, intent(in)    :: d_A_ez(ndim,maxcf,maxCorner,nzones) !ndim,maxcf,maxCorner,nzones
   integer, device, intent(in):: d_Connect(3,maxcf,maxCorner,nzones) 
   real(adqt), device, intent(in)    :: d_STotal(Groups, maxCorner, nzones)
   real(adqt), device, intent(in)    :: d_STimeBatch(Groups, ncornr, anglebatch)
   real(adqt), device, intent(in)    :: d_Volume(maxCorner, nzones)
   real(adqt), device, intent(inout) :: d_psicbatch(QuadSet%Groups,Size%ncornr,anglebatch)
   real(adqt), device, intent(inout) :: d_psib(QuadSet%Groups,Size%nbelem,QuadSet%NumAngles)
   integer,    device, intent(in) :: d_next(Size%ncornr+1,QuadSet%NumAngles)
   integer,    device, intent(in) :: d_nextZ(Size%nzones,QuadSet%NumAngles)
   real(adqt), device, intent(in)    :: d_Sigt(Groups, nzones)
   real(adqt), device, intent(in)    :: d_SigtInv(Groups, nzones)
   integer,    device, intent(in) :: d_passZstart(Size%nzones,QuadSet%NumAngles)   
   integer(kind=cuda_stream_kind), intent(in) :: streamid   

!  Local Variables

   integer    :: mm, Angle,istat,i,ib,ic

   type(dim3) :: threads,blocks


#ifdef PROFILING_ON
   integer profiler(2) / 0, 0 /
   save profiler
#endif


#ifdef PROFILING_ON
   call TAU_PROFILE_TIMER(profiler, 'snswp3d')
   call TAU_PROFILE_START(profiler)
#endif

   ! groups*NZONEPAR must be .le. 1024 on K80 hardware
   !threads=dim3(QuadSet%Groups,NZONEPAR,1) 
   threads=dim3(THREADX,1,NZONEPAR) 
   blocks=dim3(1,1,anglebatch)

   call GPU_sweep<<<blocks,threads,0,streamid>>>( anglebatch,                     &
                              nzones,               &
                              Groups,            &
                              ncornr,               &
                              NumAngles,         &
                              d_Angles,        & ! only angle batch portion
                              maxCorner,            &
                              maxcf,                &
                              binRecv,                   &
                              NangBin,                   &
                              nbelem,                &
                              d_omega,             &
                              d_nCorner,                &
                              d_nCFaces,                &
                              d_c0,                &
                              d_A_fp,                &
                              d_A_ez,                &
                              d_Connect,             &
                              d_STotal,              &
                              d_STimeBatch,               &  ! only angle batch portion
                              d_Volume,             &
                              d_psicbatch,                      &  ! only angle batch portion
                              d_psib,                      &
                              d_next,              &
                              d_nextZ,             &
                              d_Sigt,                &
!                              d_SigtInv,             &
                              d_passZstart             )




#ifdef PROFILING_ON
   call TAU_PROFILE_STOP(profiler)
#endif


   return
   end subroutine snswp3d



   subroutine scalePsibyVolume(psir, volumeRatio, anglebatch, streamid)
     ! scaling by change in mesh volume that used to be done in advanceRT is done here
     use kind_mod
     use constant_mod
     use Quadrature_mod
     use Size_mod
     use cudafor
     
     implicit none
     
     !  Arguments

     real(adqt), device, intent(inout)  :: psir(QuadSet%Groups,Size%ncornr,anglebatch) 
     real(adqt), device, intent(in) :: volumeRatio(Size%ncornr)
     integer, intent(in)  :: anglebatch 
     integer(kind=cuda_stream_kind), intent(in) :: streamid

     !  Local

     integer    :: ia, ic, ig, ncornr, Groups
     real(adqt) :: tau


     ncornr = Size%ncornr
     Groups = QuadSet% Groups   

     !$cuf kernel do(3) <<< *, *, stream=streamid >>>
     do ia=1,anglebatch
        do ic=1,ncornr
           do ig=1, Groups
             psir(ig,ic,ia) = psir(ig,ic,ia)*volumeRatio(ic)
          enddo
       enddo
     enddo


     return
   end subroutine scalePsibyVolume




   subroutine computeSTime(psiccache, STimeBatch, anglebatch, streamid)
     ! Multiply by tau to get STime
     use kind_mod
     use constant_mod
     use Quadrature_mod
     use Size_mod
     use cudafor
     
     implicit none
     
     !  Arguments

     real(adqt), device, intent(in)  :: psiccache(QuadSet%Groups,Size%ncornr,anglebatch) 
     real(adqt), device, intent(out) :: STimeBatch(QuadSet%Groups,Size%ncornr,anglebatch)
     integer, intent(in)  :: anglebatch 
     integer(kind=cuda_stream_kind), intent(in) :: streamid

     !  Local

     integer    :: ia, ic, ig, ncornr, Groups
     real(adqt) :: tau


     ncornr = Size%ncornr
     Groups = QuadSet% Groups   
     tau    = Size%tau


     !$cuf kernel do(3) <<< *, *, stream=streamid >>>
     do ia=1,anglebatch
        do ic=1,ncornr
           do ig=1, Groups
              STimeBatch(ig,ic,ia) = tau*psiccache(ig,ic,ia)
           enddo
        enddo
     enddo


     return
   end subroutine computeSTime



  attributes(global)   subroutine computeSTimeD(psiccache, STimeBatch, anglebatch, groups, ncornr, tau)
     ! Multiply by tau to get STime
     !use kind_mod
     !use constant_mod
     !use Quadrature_mod
     !use Size_mod
     !use cudafor
     
     implicit none
     
     !  Arguments

     real(adqt), device, intent(in)  :: psiccache(Groups,ncornr,anglebatch) 
     real(adqt), device, intent(out) :: STimeBatch(Groups,ncornr,anglebatch)
     integer, value, intent(in)  :: anglebatch, Groups, ncornr
     real(adqt), value, intent(in) :: tau

     !  Local

     integer    :: ia, ic, ig

     do ia=1,anglebatch !blockIdx.x
        do ic=1,ncornr !threadIdx.y
           do ig=1, Groups !threadIdx.x
              STimeBatch(ig,ic,ia) = tau*psiccache(ig,ic,ia)
           enddo
        enddo
     enddo


     return
   end subroutine computeSTimeD




 end module snswp3d_mod
