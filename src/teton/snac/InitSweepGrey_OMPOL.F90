!***********************************************************************
!                        Last Update:  09/2018, PFN                    *
!                                                                      *
!   InitGreySweepUCB  - This routine calculates the transfer           *
!                       matrices used for the direct solve for the     *
!                       scalar corrections.                            *
!                                                                      *
!***********************************************************************

   subroutine InitGreySweep

   use kind_mod
   use Size_mod
   use QuadratureList_mod

   implicit none

!  Local

   integer         :: setID
   integer         :: nAngleSets
   integer         :: zone
   integer         :: nzones

   logical(kind=1) :: useGPU

!  Constants

   nAngleSets = getNumberOfAngleSets(Quad)
   setID      = nAngleSets + 1
   nzones     = Size%nzones
   useGPU     = getGPUStatus(Size)

   if ( useGPU ) then

     if (Size% ndim == 2) then
!       call InitGreySweepUCBrz_GPU
     else
!       call InitGreySweepUCBxyz_GPU
     endif

   else

     if (Size% ndim == 2) then

!$omp parallel do default(none) schedule(static) &
!$omp& shared(setID, nzones)
       do zone=1,nzones
         call InitGreySweepUCBrz(setID, zone)
       enddo
!$omp end parallel do

     else

!$omp parallel do default(none) schedule(static) &
!$omp& shared(setID, nzones)
       do zone=1,nzones
         call InitGreySweepUCBxyz(setID, zone)
       enddo
!$omp end parallel do

     endif

   endif


   return
   end subroutine InitGreySweep

