!***********************************************************************
!                        Version 0:  02/2014, PFN                      *
!                                                                      *
!   ConstructRadIntensity - Sets F90 module pointers for the angle-    *
!                           dependent (Psi) and scalar (Phi) radiation *
!                           intensity to memory allocated by the       *
!                           host code.                                 *
!                                                                      *
!***********************************************************************


   subroutine ConstructRadIntensity() &
        BIND(C,NAME="teton_constructradintensity")

!  Include

   USE ISO_C_BINDING
   use kind_mod
   use QuadratureList_mod
   use RadIntensity_mod


   implicit none

!  Local

   type(RadIntensity), pointer  :: RadT

   integer :: setID
   integer :: nSets
   integer :: Groups

!  Construct Material Module 

   nSets = getNumberOfSets(Quad)

   do setID=1,nSets

     RadT      => getRadIntensity(Quad, setID)
     Groups    =  getNumberOfGroups(Quad, setID)
 
     call construct(RadT, Groups)

   enddo


   return
   end subroutine ConstructRadIntensity

