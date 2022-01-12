!***********************************************************************
!                                                                      *
!   ReConstructSetsFromdV - Reconstructs the angle-dependent           * 
!   and angle-integrated intensities after Lagrange mesh motion        *
!   from on changes to the pre/post remap volumes                      *
!                                                                      *
!***********************************************************************
 
   subroutine ReConstructSetsFromdV(setID, Groups)

   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use QuadratureList_mod
   use Geometry_mod
   use SetData_mod

   implicit none 

!  Arguments

   integer,    intent(in)    :: setID
   integer,    intent(in)    :: Groups

!  Local

   type(SetData),   pointer  :: Set

   integer    :: angle 
   integer    :: g 
   integer    :: g0
   integer    :: zone 
   integer    :: NumAngles 
   integer    :: nZones 
   integer    :: nCorner 
   integer    :: c 
   integer    :: c0

   real(adqt) :: factor(Groups) 

!  Parameters 

   Set  => getSetData(Quad, setID)

   nZones    = Size% nzones
   NumAngles = Set% NumAngles
   g0        = Set% g0

!***********************************************************************
!  Adjust the discrete angular intensities to conserve energy          *
!***********************************************************************

   ZoneLoop: do zone=1,nZones

     nCorner = Geom% numCorner(zone)
     c0      = Geom% cOffSet(zone)

!  Compute multiplicative correction

!  Apply correction uniformly based on corner volume changes from the remap 

     do c=1,nCorner

       do angle=1,NumAngles
         Set% Psi(:,c0+c,angle) = &
             (Geom% VolumeOld(c0+c)/Geom% Volume(c0+c)) * Set% Psi(:,c0+c,angle)
       enddo

     enddo

   enddo ZoneLoop



   return
   end subroutine ReConstructSetsFromdV 


