!***********************************************************************
!                        Last Update:  03/2012, PFN                    *
!                                                                      *
!    SetNodePosition -  Called from host to set node positions         *
!                       in the MeshData module.                        *
!                                                                      *
!***********************************************************************

   subroutine setNodePosition(zoneID, nodePosition) &
        BIND(C,NAME="teton_setnodeposition")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use Geometry_mod
   use MeshData_mod
   use radconstant_mod

   implicit none

!  Arguments

   integer(C_INT),    intent(in) :: zoneID 
   real(C_DOUBLE), intent(in) :: nodePosition(Size% ndim,Size% maxCorner) 

!  Local

   integer  :: c, nCorner

   M => getMesh(Geom, zoneID)

   nCorner = M% nCorner

   do c=1,nCorner
     M% px(:,c)  = nodePosition(:,c)
   enddo


   return
   end subroutine setNodePosition 
