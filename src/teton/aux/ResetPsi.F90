!***********************************************************************
!                        Last Update:  08/2019, TAB                    *
!                                                                      *
!   ResetPsi - Set the intensity in a zone to zero and tally energy as *
!              escaped.                                                *
!                                                                      *
!***********************************************************************

   subroutine ResetPsi(nVoidZones, voidZoneList)  &
         BIND(C,NAME="teton_resetpsi")

   use ISO_C_BINDING
   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use Geometry_mod
   use RadIntensity_mod
   use Material_mod
   use SetData_mod

   implicit none

!  Arguments

   integer(C_INT), intent(in) :: nVoidZones
   integer(C_INT), intent(in) :: voidZoneList(nVoidZones) 

!  Local

   type(SetData),  pointer    :: Set

   integer    :: n
   integer    :: zone
   integer    :: setID
   integer    :: nSets
   integer    :: nCorner
   integer    :: c0
   integer    :: c
   integer    :: angle
   integer    :: g

!  Constants

   nSets = getNumberOfSets(Quad)

!  Reset the void flags

   Mat% isVoid(:) = .FALSE.

!  Loop over the void zone list

   VoidZoneLoop: do n=1,nVoidZones

     zone    = voidZoneList(n)
     nCorner = Geom% numCorner(zone)
     c0      = Geom% cOffSet(zone)

! TODO: Tally energy lost as escaped.

     SetLoop: do setID=1,nSets

       Set => getSetData(Quad, setID)

       do angle=1,Set% NumAngles
         do c=1,nCorner
           Set% Psi(:,c0+c,angle) = zero
         enddo
       enddo

     enddo SetLoop

!    Zero energy densities

     do g=1,Size% ngr
       Rad% RadEnergyDensity(zone,g) = zero
     enddo

     do c=1,nCorner
       Rad% PhiTotal(:,c0+c) = zero
     enddo

!    Need to remove this zone from time step control

     Mat% isVoid(zone) = .TRUE.

   enddo VoidZoneLoop


   return
   end subroutine ResetPsi
