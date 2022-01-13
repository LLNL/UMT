subroutine harmonics(YLK,omega,num_harm,num_degree,num_angle,num_dim)


!  This subroutine assigns the values of an orthonormal set of surface
!  harmonics functions for all directions in the quadrature set.  The
!  number of surface harmonics required depends on the number of space
!  dimensions:

!     1D:  num_harm == num_degree + 1
!     2D:  num_harm == (num_degree+1)*(num_degree+2)/2
!     3D:  num_harm == (num_degree+1)^2

!  where the function expansion is over l=0 ... num_degree.

!  Upon return, YLK(:,:) contains each degree and order of the surface
!  harmonics evaluated at each discrete ordinate; ier<0 if an error has
!  occured.

!  Currently, the 2D option assumes RZ geometry; this could easily be
!  generalized to XY by passing in a geometry flag.

use kind_mod
use io_mod

!  Variable declarations
implicit none

!  Arguments
integer,    intent(in)     :: num_harm,num_degree,num_angle,num_dim

real(adqt), intent(in)     :: omega(num_dim,num_angle)
real(adqt), intent(inout)  :: ylk(num_harm,num_angle)

!  Local Variables
integer  :: num_harm_required

!  Assign the surface harmonics for the various geometries

      select case (num_dim)

!  One-dimensional geometries
      case (1)

!       Error check the number of moments

        num_harm_required = num_degree + 1
        if (num_harm /= num_harm_required) then
          write(nout,900) num_degree,num_harm_required,num_harm
          call f90fatal("Number of harmonics is inconsistent with degree")
        endif

!       X geometry:
!         mu[m] == direction cosine with respect to the X-axis
!         Y[l,k](OMEGA[m]) = Y[l,0](mu[m])

        call harmonics_1d(YLK,omega(1,:),num_harm,num_degree,num_angle)

!  Two-dimensional geometries
      case (2)

!       Error check the number of moments

        num_harm_required = (num_degree + 1)*(num_degree + 2)/2
        if (num_harm /= num_harm_required) then
          write(nout,900) num_degree,num_harm_required,num_harm
          call f90fatal("Number of harmonics is inconsistent with degree")
        endif

!       RZ geometry:
!          mu[m] == direction cosine with respect to the R-axis
!         eta[m] == direction cosine with respect to the Z-axis
!           w[m] == angle between R-axis and the projection of
!                   OMEGA[m] onto the RTheta-plane
!         Y[l,k](OMEGA[m]) = Y[l,k](eta[m],w[m])

        call harmonics_2d(YLK,omega(2,:),omega(1,:),num_harm, &
                          num_degree,num_angle)

!       XY geometry:
!          mu[m] == direction cosine with respect to the X-axis
!         eta[m] == direction cosine with respect to the Y-axis
!           w[m] == angle between Y-axis and the projection of
!                   OMEGA[m] onto the YZ-plane
!         Y[l,k](OMEGA[m]) = Y[l,k](mu[m],w[m])

        call harmonics_2d(YLK,omega(1,:),omega(2,:),num_harm, &
                          num_degree,num_angle)

!  Three-dimensional geometries
      case(3)

!       Error check the number of moments

        num_harm_required = (num_degree + 1)**2
        if (num_harm /= num_harm_required) then
          write(nout,900) num_degree,num_harm_required,num_harm
          call f90fatal("Number of harmonics is inconsistent with degree")
        endif

!       XYZ geometry:
!          mu[m] == direction cosine with respect to the X-axis
!         eta[m] == direction cosine with respect to the Y-axis
!          xi[m] == direction cosine with respect to the Z-axis
!           w[m] == angle between Y-axis and the projection of
!                   OMEGA[m] onto the YZ-plane
!         Y[l,k](OMEGA[m]) = Y[l,k](mu[m],w[m])

        call harmonics_3d(YLK,omega(1,:),omega(2,:),omega(3,:), &
                          num_harm,num_degree,num_angle)

      end select

!  Format statements

900  format(1x,"** Error **"/3x,"Solver Error:"/6x,"Harmonics of ", &
            "degree ",i4," requested, which requires storage of "/i4, &
            "harmonics; storage for ",i4," harmonics was allocated.")


return

end subroutine harmonics











