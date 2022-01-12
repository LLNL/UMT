!***********************************************************************
!                       Last Update:  07/2017, TSH                     *
!                                                                      *
!   getPsiPointer - Returns a pointer to the group/angle dependent     *
!                   radiation intensity (Psi) in a phase-space "set"   *
!                   and it's associated dimensions.                    *
!                                                                      *
!***********************************************************************
 
   subroutine getPsiPointer(setID, PsiDims, PsiPtr) &
        BIND(C,NAME="teton_getpsipointer")

   USE ISO_C_BINDING
   use kind_mod
   use QuadratureList_mod
   use SetData_mod

   implicit none 

!  Arguments

   integer(C_INT), intent(in)    :: setID
   integer(C_INT), intent(out)   :: PsiDims(3)
   TYPE(C_PTR),    intent(out)   :: PsiPtr

!  Local
   type(SetData), pointer  :: Set

   Set        => getSetData(Quad, setID+1)

   PsiDims(1) =  Set% Groups
   PsiDims(2) =  Set% nCorner 
   PsiDims(3) =  Set% NumAngles

   PsiPtr =  C_LOC( Set% Psi )


   return
   end subroutine getPsiPointer 

