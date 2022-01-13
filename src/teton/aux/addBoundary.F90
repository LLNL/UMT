!***********************************************************************
!                         LastUpdate: 02/2012 PFN                      *
!                                                                      *
!    AddBoundary  -  Called from host to add a boundary                *
!                    to the boundary list.                             *
!                                                                      *
!***********************************************************************
   subroutine addBoundary(numBCTotal, BCType, NumBdyElem,    & 
                          NeighborID) & 
                          BIND(C,NAME="teton_addboundary")


   use kind_mod
   use BoundaryList_mod
   use Boundary_mod
   use flags_mod

   use ISO_C_BINDING

   implicit none


!  Arguments

   integer(C_INT), intent(in)          :: numBCTotal

 ! array of boundary types corresponding to the parameters given flags_mod.F90
   integer(C_INT), intent(in)          :: BCType(numBCTotal) 

   integer(C_INT), intent(in)          :: NumBdyElem(numBCTotal)
   integer(C_INT), intent(in)          :: NeighborID(numBCTotal) 

!  Local

   integer          :: BdyID
   integer          :: nReflecting
   integer          :: BdyElem1
              
!  Add this profile to the list 

   BdyElem1 = 1

   do BdyID=1,numBCTotal

     if (BCType(BdyID) == bcType_none ) then 
       call f90fatal("addBoundary: Radiation boundary condition is not set")
     else

       call setBoundary(RadBoundary,       &
                        BdyID,             &
                        BCType(BdyID),     &
                        NumBdyElem(BdyID), &
                        BdyElem1,          &
                        NeighborID(BdyID))

       BdyElem1 = BdyElem1 + NumBdyElem(BdyID)

     endif

   enddo

!  Allocate space for reflecting boundary data

   nReflecting = getNumberOfReflecting(RadBoundary)

   do BdyID=1,nReflecting
     Bdy => getReflecting(RadBoundary,BdyID)

     call constructReflectedAngle(Bdy)
   enddo

   

   return
   end subroutine addBoundary 

