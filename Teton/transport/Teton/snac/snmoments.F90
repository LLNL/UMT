!***********************************************************************
!                        Version 1:  05/95, PFN                        *
!                                                                      *
!   SNMOMENTS - This routine, called by SNFLWRZA and SNFLWXYZ          *
!               calculates the required spherical harmonic moments     *
!               [phic] of the angular flux [psic]. It uses the array   *
!               ynm(n,m), whose definition is:                         *
!                                                                      *
!               ynm(n,m) = real part of (l,k)th spherical harmonic,    *
!                          evaluated at the mth direction, where       *
!                                                                      *
!                             n = 1 + l*(l+1)/2 + k                    *
!                                                                      *
!            This routine is designed to accumulate moments as an      *
!            angle is calculated and does not require storage of the   *
!            full angular flux array.  It is assumed that the moment   *
!            array has been initialized before the loop over angles    *
!            in SNFLWRZA or SNFLWXYZ.                                  *
!                                                                      *
!                                                                      *
!   Input:   psic     - angular flux                   (E/A/t/ster)    *
!            quadwt   - quadrature weights                      (0)    *
!            ynm      - spherical harmonics                     (0)    *
!                                                                      *
!   Output:  PHIC     - flux moments                        (E/A/t)    *
!                                                                      *
!***********************************************************************

   subroutine snmoments(psic, PHI)

   use kind_mod
   use constant_mod
   use Quadrature_mod
   use Size_mod

   implicit none

!  Arguments

   real(adqt), intent(in)    :: psic(QuadSet%Groups,Size%ncornr,QuadSet%NumAngles) 

   ! Why is this inout, phi is set to zero inside routine. Looks like it is only out.
   real(adqt), intent(inout) :: Phi(QuadSet%NumMoments*QuadSet%Groups,Size%ncornr)

!  Local

   integer    :: ic, ig, Angle, Groups, ncornr
   integer :: tid,nth,icbeg,icend
   integer, external :: omp_get_thread_num, omp_get_num_threads

   real(adqt) :: quadwt 

!  Add this angles contribution to the flux moments

   Groups = QuadSet% Groups
   ncornr = Size% ncornr

   Phi(:,:) = zero

!$omp parallel private(quadwt,tid,nth,icbeg,icend)
   tid = omp_get_thread_num()
   nth = omp_get_num_threads()
   !do c0=1,ncornr,512 ! cache blocking loop
   !call omp_block_partition(tid,nth,c0,min(ncornr,c0+511),icbeg,icend)
   call omp_block_partition(tid,nth,1,ncornr,icbeg,icend)

   AngleLoop: do Angle=1,QuadSet%NumAngles

     quadwt = QuadSet% Weight(Angle)

     if (quadwt /= zero) then

       do ic = icbeg, icend     ! YKT, was :  ic=1,ncornr
         do ig=1,Groups
           Phi(ig,ic) = Phi(ig,ic) + quadwt*psic(ig,ic,Angle)
         enddo
       enddo

     endif

   enddo AngleLoop

   !enddo
!$omp end parallel


   print *, "host version of snmoments still called here!!!!"
 
   return
   end subroutine snmoments

!!!!!! Device version of the same thing


   subroutine snmomentsD(psiccache, PHI, weight, angleorder, angles, streamid)

   use kind_mod
   use constant_mod
   use Quadrature_mod
   use Size_mod
   use cudafor

   implicit none

!  Arguments

   integer, intent(in) :: angles

   real(adqt), device, intent(in)  :: psiccache(QuadSet%Groups,Size%ncornr,angles) 
   real(adqt), device, intent(out) :: Phi(QuadSet%NumMoments*QuadSet%Groups,Size%ncornr)
   real(adqt), device, intent(in)  :: weight(QuadSet%NumAngles)
   integer, device, intent(in)  :: angleorder(angles)

   integer(kind=cuda_stream_kind), intent(in) :: streamid

!  Local

   integer    :: ic, ig, Groups, ncornr, istat

!  Add this angles contribution to the flux moments

   Groups = QuadSet% Groups
   ncornr = Size% ncornr

   !$cuf kernel do(2) <<< (*,*), (16,16), stream=streamid >>>
   do ic=1,ncornr
     do ig=1,Groups
       Phi(ig,ic) = Phi(ig,ic) + SUM(weight(angleorder(:))*psiccache(ig,ic,:))
     enddo
   enddo
 
   return
   end subroutine snmomentsD

