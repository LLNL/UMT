!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   SETTOTALOPACITY - Constructs a total opacity for transport         *
!                     comprised of absorption, time-absorption and     *
!                     scattering.                                      *
!                                                                      *
!***********************************************************************
 
   subroutine setTotalOpacity(gSetID) 


   use kind_mod
   use Size_mod
   use constant_mod
   use Material_mod
   use QuadratureList_mod
   use GroupSet_mod

   implicit none

!  Arguments

   integer, intent(in)     :: gSetID

!  Local

   type(GroupSet), pointer :: GSet

   integer    :: zone
   integer    :: nzones
   integer    :: g
   integer    :: g0
   integer    :: Groups


!  Constants

   GSet   => getGroupSetData(Quad, gSetID)

   nzones =  Size%nzones
   g0     =  GSet% g0 
   Groups =  GSet% Groups

!***********************************************************************
!     ADD TIME-ABSORPTION TO THE TOTAL CROSS SECTION                   *
!***********************************************************************

   ZoneLoop: do zone=1,nzones

     do g=1,Groups
       GSet% Sigt(g,zone) = Mat%SigA(g0+g,zone) + Mat%SigS(g0+g,zone) + Size% tau
     enddo

   enddo ZoneLoop


   return
   end subroutine setTotalOpacity 

