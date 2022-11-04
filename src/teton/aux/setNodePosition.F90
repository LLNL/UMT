!***********************************************************************
!                        Last Update:  03/2012, PFN                    *
!                                                                      *
!    SetNodePosition -  Called from host to set node positions         *
!                       in the Geometry module.                        *
!                                                                      *
!***********************************************************************

   subroutine setNodePosition(zoneID, nodePosition) &
        BIND(C,NAME="teton_setnodeposition")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Arguments

   integer(C_INT), intent(in) :: zoneID 
   real(C_DOUBLE), intent(in) :: nodePosition(Size% ndim,Size% maxCorner) 

!  Local

   integer  :: c, c0, nCorner

   nCorner = Geom% numCorner(zoneID)
   c0      = Geom% cOffSet(zoneID) 

   do c=1,nCorner
     Geom% px(:,c0+c) = nodePosition(:,c)
   enddo


   return
   end subroutine setNodePosition 
