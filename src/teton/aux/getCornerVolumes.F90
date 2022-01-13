
!***********************************************************************
!                         Version 0: 03/11/21 TSH                      *
!                                                                      *
!         getCornerVolumes -  Called from host to get the              *
!                             corner volumes in the zone               *
!                                                                      *
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

!   Set the corner volumes
    c0 = Geom% cOffSet(zone)
    nCorner = Geom% numCorner(zone)
    do c=1,nCorner
      cornerVolumes(c) = Geom% Volume(c0+c)
    enddo

  end subroutine getCornerVolumes
