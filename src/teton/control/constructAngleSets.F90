!***********************************************************************
!                        Last Update:  11/2022, PFN                    *
!                                                                      *
!   constructAngleSets - Creates the quadratures used for high-order   *
!                        and GTA transport sweeps.                     *
!                                                                      *
!***********************************************************************
 
   subroutine constructAngleSets(QuadDef, gnu) 

   use kind_mod
   use Size_mod
   use QuadratureList_mod

   implicit none

!  Arguments

   integer,    intent(in) :: QuadDef(6,2)
   real(adqt), intent(in) :: gnu(Size%ngr+1)

!  Local

   integer    :: quadID
   integer    :: quadType
   integer    :: Order
   integer    :: Groups
   integer    :: NumAngles
   integer    :: NumMoments
   integer    :: NPolar
   integer    :: NAzimuthal
   integer    :: PolarAxis
   real(adqt) :: GrpBnds(Size%ngr+1)

!  Find the number of angle sets and the number of groups in them 

   Groups     = 0 
   NumMoments = 1 

   do quadID=1,2 

!  Set the current quadrature definition

     quadType   = QuadDef(1,quadID)
     Order      = QuadDef(2,quadID)
     NPolar     = QuadDef(3,quadID)
     NAzimuthal = QuadDef(4,quadID)
     PolarAxis  = QuadDef(5,quadID)
     NumAngles  = QuadDef(6,quadID)

     if (quadID == 1) then
       Groups     = Size% ngr
       GrpBnds(:) = Gnu(:)
     elseif (quadID == 2) then
       Groups     = 1
       GrpBnds(1) = Gnu(1) 
       GrpBnds(2) = Gnu(Size% ngr+1)
     endif

     call setQuadrature(Quad,          &
                        quadID,        &
                        Groups,        &
                        NumAngles,     &
                        NumMoments,    &
                        Order,         &
                        NPolar,        &
                        NAzimuthal,    &
                        PolarAxis,     &
                        quadType,      &
                        GrpBnds)

   enddo 


   return
   end subroutine constructAngleSets 

