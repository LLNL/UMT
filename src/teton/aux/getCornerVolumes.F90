
!***********************************************************************
!                         Version 0: 03/11/21 TSH                      *
!                                                                      *
!         getCornerVolumes -  Called from host to get the corner       *
!                             volumes in the zone.                     *
!                                                                      *
! 02/07/2022 - The value returned by this function was updated to      *
!              be scaled by the geometry factor.                       *
!***********************************************************************

  subroutine getCornerVolumes(zone, cornerVolumes) &
                              BIND(C,NAME="teton_getcornervolumes")

    USE ISO_C_BINDING
    use kind_mod
    use Size_mod
    use Geometry_mod

    implicit none

!   Arguments

    integer(C_INT),   intent(in)     :: zone
    real(C_DOUBLE),   intent(inout)  :: cornerVolumes(Size%maxCorner) 

!   Local

    integer    :: c 
    integer    :: c0 
    integer    :: nCorner 
    real(adqt) :: geometryFactor

!   Set the corner volumes
    c0 = Geom% cOffSet(zone)
    geometryFactor = getGeometryFactor(Size)
    nCorner = Geom% numCorner(zone)
    do c=1,nCorner
      cornerVolumes(c) = geometryFactor * Geom% Volume(c0+c)
    enddo

  end subroutine getCornerVolumes
