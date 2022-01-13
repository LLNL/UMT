!***********************************************************************
!                       Last Update:  02/2017, PFN                     *
!                                                                      *
!   queryPsi      - Returns a value of the  group/angle dependent      *
!                   radiation intensity (Psi) given three indices      *
!                   (group,corner,angle).                              *
!                                                                      *
!***********************************************************************
 
   subroutine queryPsi(PsiIndices, PsiValue) BIND(C,NAME="teton_querypsi")

   USE ISO_C_BINDING
   use kind_mod
   use QuadratureList_mod
   use SetData_mod

   implicit none 

!  Arguments

   integer(C_INT),    intent(in)    :: PsiIndices(3) 
   real(C_DOUBLE), intent(inout) :: PsiValue

!  Local

   type(SetData),   pointer  :: Set

   integer                   :: setID
   integer                   :: groupID
   integer                   :: cornerID
   integer                   :: angleID
   integer                   :: group
   integer                   :: angle

!  Find the correct set

   groupID    =  PsiIndices(1) + 1
   cornerID   =  PsiIndices(2) + 1
   angleID    =  PsiIndices(3) + 1

   setID      =  getSetIDfromGroupAngle(Quad, groupID, angleID)
   Set        => getSetData(Quad, setID)

   group      =  groupID - Set% g0
   angle      =  angleID - Set% angle0


   PsiValue   =  Set% Psi(group,cornerID,angle)


   return
   end subroutine queryPsi 


