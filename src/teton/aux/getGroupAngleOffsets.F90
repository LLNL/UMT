!***********************************************************************
!                       Last Update:  08/2017, PGM                     *
!                                                                      *
!   getGroupAngleOffsets - Return a Group/Angle phase space (Set)'s    *
!                          Global group and angle number offsets       *
!                                                                      *
!***********************************************************************
 
   subroutine getGroupAngleOffsets(setID, g0, angle0) &
        BIND(C,NAME="teton_getgroupangleoffsets")

   USE ISO_C_BINDING
   use QuadratureList_mod
   use SetData_mod

   implicit none 

!  Arguments

   integer(C_INT), intent(in)    :: setID
   integer(C_INT), intent(out)   :: g0
   integer(C_INT), intent(out)   :: angle0

!  Local
   type(SetData), pointer  :: Set

   Set        => getSetData(Quad, setID+1)

   g0     =  Set% g0
   angle0 =  Set% angle0

   return
   end subroutine getGroupAngleOffsets 

