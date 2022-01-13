!***********************************************************************
!                        Last Update:  09/2018, PFN                    *
!                                                                      *
!   UpdateScalarIntensity:                                             *
!                                                                      *
!   Invert the transport operator to solve for an updated value of     *
!   the scalar intensity.
!                                                                      *
!                                                                      *
!***********************************************************************
   subroutine ScalarIntensityDecompose(zone, P, withSource)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use GreyAcceleration_mod


   implicit none

!  Arguments

   integer,          intent(in)    :: zone

   real(adqt),       intent(inout) :: P(Size% ncornr)
   logical (kind=1), intent(in)    :: withSource

!  Local

   integer    :: c
   integer    :: cc 
   integer    :: c0
   integer    :: nCorner
   integer    :: i, j, k

   real(adqt) :: wtiso
   real(adqt) :: t
   real(adqt) :: v
   real(adqt) :: diagInv

   real(adqt) :: Phi(Size% maxCorner)

   wtiso = Size% wtiso

!  Update Scalar Intensity (solve A*Phi = S)

   nCorner =  Geom% numCorner(zone)
   c0      =  Geom% cOffSet(zone)

   if ( withSource ) then

     do c=1,nCorner
       Phi(c) = GTA% PhiInc(c0+c)
       do cc=1,nCorner
         Phi(c)           =  Phi(c) + GTA% TT(cc,c0+c)*  &
                             wtiso*GTA%GreySource(c0+cc)
         GTA% TT(cc,c0+c) = -wtiso*GTA%GreySigScat(c0+cc)*GTA% TT(cc,c0+c)
       enddo
       GTA% TT(c,c0+c)    = one + GTA% TT(c,c0+c)
     enddo

!    Decompose:  A = LU 

     do i=1,nCorner

       t = zero

       do k=1,i-1
         t = t + GTA% TT(k,c0+i)*GTA% TT(i,c0+k)
       enddo

       GTA% TT(i,c0+i) = GTA% TT(i,c0+i) - t
       diagInv         = one/GTA% TT(i,c0+i)


       do j=i+1,nCorner

         t = zero
         v = zero

         do k=1,i-1
           t = t + GTA% TT(k,c0+i)*GTA% TT(j,c0+k)
           v = v + GTA% TT(k,c0+j)*GTA% TT(i,c0+k)
         enddo

         GTA% TT(j,c0+i) = GTA% TT(j,c0+i) - t
         GTA% TT(i,c0+j) = diagInv*(GTA% TT(i,c0+j) - v)

       enddo

     enddo

   else

     do c=1,nCorner
       Phi(c) = GTA% PhiInc(c0+c)
     enddo

   endif

!  Solve

   call ScalarIntensitySolve(c0, nCorner, Phi)

   do c=1,nCorner
     P(c0+c) = Phi(c)
   enddo


   return
   end subroutine ScalarIntensityDecompose

!********************************************************************
!                     Last Update:  06/2018, PFN                    *
!                                                                   *
!   solveLinearSystem - Solves a small and dense matrix using       *
!                       Gaussian elimination with back              *
!                       substitution.                               *
!                                                                   *
!********************************************************************

   subroutine ScalarIntensitySolve(c0, nCorner, Phi)

   use kind_mod
   use Size_mod
   use constant_mod
   use GreyAcceleration_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: c0 
   integer,    intent(in)    :: nCorner

   real(adqt), intent(inout) :: Phi(Size% maxCorner)

! Local

   integer                   :: i, j, k

   real(adqt)                :: t

!  Solve Ly = S

   do j=2,nCorner
     t = zero
     do i=1,j-1
       t = t - GTA% TT(i,c0+j)*Phi(i)
     enddo
     Phi(j) = Phi(j) + t
   enddo

!  Solve Ux = y

   Phi(nCorner) = Phi(nCorner)/GTA% TT(nCorner,c0+nCorner)

   do k=nCorner-1,1,-1
     t = zero

     do i=k+1,nCorner
       t = t + Phi(i)*GTA% TT(i,c0+k)
     enddo

     Phi(k) = (Phi(k) - t)/GTA% TT(k,c0+k)
   enddo

   return
   end subroutine ScalarIntensitySolve 
                                   
