!***********************************************************************
!                       Last Update:  07/2017, TSH                     *
!                                                                      *
!   getNumSNSets - Returns # of SN sets          .                     *
!                                                                      *
!***********************************************************************
 
   subroutine getNumSNSets(numSNSets) &
        BIND(C,NAME="teton_getnumsnsets")

   USE ISO_C_BINDING
   use QuadratureList_mod, only : Quad

   implicit none

!  Arguments

   integer(C_INT), intent(out)    :: numSNSets

   numSNSets = 1 

   return
   end subroutine getNumSNSets








