!***********************************************************************
!                        Version 0:  01/97, PFN                        *
!                                                                      *
!   RADTR  - Control program for radiation transport. It initializes   *
!            arrays, controls timestep, calls the transport package    *
!            and performs edits.                                       *
!                                                                      *
!                                                                      *
!***********************************************************************

   subroutine radtr(PSIR, PHI, RADENERGYDENSITY, angleLoopTime)

!  Include

   use, intrinsic :: iso_c_binding
   use kind_mod
   use Size_mod
   use Material_mod
   use Geometry_mod
   use TimeStepControls_mod
   use ProfileList_mod
   use ZoneData_mod
   use constant_mod
   use radconstant_mod

   use cudafor
   use nvtx_mod
   !use GPUhelper_mod


   implicit none

!  Arguments

!  Types
                                                                                         
   real(adqt), intent(inout) :: psir(Size%ngr,Size%ncornr,Size%nangSN),  &
                                Phi(Size%ngr,Size%ncornr),               &
                                RadEnergyDensity(Size%nzones,Size%ngr), angleLoopTime

 
!  Photon Intensities on the problem boundary
   real(adqt), pinned, allocatable, save :: psib(:,:,:)


!  Local

   integer    :: nbelem, ngr, nangSN

   integer    :: zone, istat

   real(adqt) :: dtrad

!  Allocate Memory 

   call constructPsiInc(SourceProfiles)

!***********************************************************************
!                                                                      *
!     ALLOCATE MEMORY                                                  *
!                                                                      *
!***********************************************************************
 
!  Set some scalars used for dimensioning

   nbelem   = Size%nbelem
   ngr      = Size%ngr
   nangSN   = Size%nangSN


   call nvtxStartRange("Allocate psib",6)

!  Photon Intensities on the problem boundary

   if (.not. allocated(psib) ) then
     allocate( psib(ngr,nbelem,nangSN) )
!    print *, "sizeof(psib): ", sizeof(psib)
   endif

   call nvtxEndRange

!***********************************************************************
!     UPDATE DEVICE ZONE DATA (once only)                              *
!***********************************************************************

   call nvtxStartRange("Setup ZData",5)

   ! this gives a device valid way to reference ZData%d_member where 
   ! d_member is a device memory location. It copies the addresses of device 
   ! memory locations (already stored in ZData) to the device structure d_ZData. 

   if (Geom%d_ZData_uptodate == .false.) then
     istat = cudaMemcpyAsync(C_DEVLOC(Geom%d_ZData), C_LOC(Geom%ZData), sizeof(Geom%ZData), 0)
     Geom%d_ZData_uptodate = .true.
   endif

   call nvtxEndRange

!***********************************************************************
!     ADD TIME-ABSORPTION TO THE TOTAL CROSS SECTION ARRAY             *
!***********************************************************************

   dtrad = getRadTimeStep(DtControls)

   if (Size%ittyp == 'timedep') then
     Size%tau = one/(speed_light*dtrad)
   else
     Size%tau = zero
   endif

   call nvtxStartRange("Sigt and Inv",5)

!$omp parallel do
   do zone=1,Size%nzones
     Z => getZoneData(Geom, zone)

     Z% Sigt(:)    = Mat%siga(:,zone) + Mat%sigs(:,zone) + Size%tau
     Z% SigtInv(:) = one/Z% Sigt(:)
   enddo

   call nvtxEndRange

!***********************************************************************
!     INTERPOLATE SOURCE PROFILES                                      *
!***********************************************************************

   call profint
 
!***********************************************************************
!     BOUNDARY CONDITIONS                                              *
!***********************************************************************
 
   call rtbdry

!***********************************************************************
!     VOLUME RADIATION SOURCES                                         *
!***********************************************************************

   call rtvsrc


!***********************************************************************
!     RADIATION TRANSPORT MODULE                                       *
!***********************************************************************

   call rtmainsn(dtrad, PSIR, PHI, psib,  angleLoopTime) 

!***********************************************************************
!     ENERGY UPDATE                                                    *
!***********************************************************************
         
   ! call newenergy
   call rtave(Mat%denec, Mat%DENEZ)

!***********************************************************************
!     EDITS                                                            *
!***********************************************************************

   call rtedit(dtrad, Phi)

!***********************************************************************
!     TIME STEP CONTROL                                                *
!***********************************************************************

   call dtnew

!***********************************************************************
!     UPDATE RADIATION MOMENTS                                         *
!***********************************************************************

   !in: Phi
   !out: RadEnergyDensity(zone,groups)
   call RadMoments(Phi, RadEnergyDensity)

!***********************************************************************
!     RELEASE MEMORY                                                   *
!***********************************************************************

   call destructPsiInc(SourceProfiles)


   return
   end subroutine radtr



   
   subroutine pinmem(PSIR, PHI)
     
     !  Include
     
     use, intrinsic :: iso_c_binding
     use kind_mod
     use Size_mod

     use cudafor
     use nvtx_mod

     implicit none

     !  Arguments

     real(adqt), intent(inout) :: psir(Size%ngr,Size%ncornr,Size%nangSN),  &
          Phi(Size%ngr,Size%ncornr)

     !  Local

     integer    :: istat


     ! DA Nov 2017: no longer pin psi and phi in Fortran,
     ! instead pin when allocated in C++ stl vector.

     call nvtxStartRange("Pin Arrays",6)

     ! pin psi and phi

     istat = cudaHostRegister(C_LOC(psir(1,1,1)), int(Size%ngr,KIND=8)&
          *int(Size%ncornr,KIND=8)&
          *int(Size%nangSN,KIND=8)*8, cudaHostRegisterMapped)
!    print *, "size of psir: ", sizeof(psir)
!    print *, "dimensions of psir: ", Size%ngr,Size%ncornr,Size%nangSN
!    print *, "Correct size used is:", int(Size%ngr,KIND=8)&
!         *int(Size%ncornr,KIND=8)&
!         *int(Size%nangSN,KIND=8)*8
     if(istat .ne. 0) then
        print *, "pinning error, istat = ", istat , LOC(psir(1,1,1))
        !print *, cudaGetErrorString(istat)
     endif


!    print *, "pinning phi, sizeof(phi) = ", sizeof(phi)
     istat = cudaHostRegister(C_LOC(phi(1,1)), sizeof(phi), cudaHostRegisterMapped)


     call nvtxEndRange

   end subroutine pinmem
