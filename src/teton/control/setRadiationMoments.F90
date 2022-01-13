!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   setRadiationMoments - Calculates the radiation force and           *
!                         radiation flux in a corner.                  *
!                                                                      *
!***********************************************************************
 
   subroutine setRadiationMoments(setID, Force) 

   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use Quadrature_mod
   use ZoneData_mod
   use Material_mod
   use SetData_mod
   use AngleSet_mod
   use RadIntensity_mod

   implicit none 

!  Arguments

   integer,         intent(in) :: setID
   logical(kind=1), intent(in) :: Force

!  Local

   type(SetData),          pointer  :: Set
   type(AngleSet),         pointer  :: ASet
   type(RadIntensity),     pointer  :: RadT

   integer    :: angle 
   integer    :: c 
   integer    :: c0 
   integer    :: nCorner 
   integer    :: g 
   integer    :: g0
   integer    :: Groups 
   integer    :: zone 
   integer    :: nzones
 
   real(adqt) :: depRate 
   real(adqt) :: geometryFactor 
   real(adqt) :: constant
   real(adqt) :: VolFrac
   real(adqt) :: wOmega(Size%ndim) 
   real(adqt) :: wOmegaOmega(Size%ndim)

!  Constants

   Set  => getSetData(Quad, setID)
   ASet => getAngleSetFromSetID(Quad, setID) 
   RadT => getRadIntensity(Quad, setID)

   Groups = getNumberOfGroups(Quad, setID)
   g0     = Set% g0
   nzones = Size% nzones

   geometryFactor = getGeometryFactor(Size)
   constant       = Size% radForceMultiplier*geometryFactor/speed_light

!***********************************************************************
!  Compute the radiation force on the matter and total radiation flux  *
!***********************************************************************

   if ( Force ) then

     ZoneLoop1: do zone=1,nzones

       nCorner = Geom% numCorner(zone)
       c0      = Geom% cOffSet(zone)

       do c=1,nCorner
         RadT% RadiationForce(:,c0+c) = zero
       enddo

       AngleLoop1: do angle=1,Set% NumAngles

         wOmega(:) =  ASet% Weight(angle)*ASet% omega(:,angle)

         CornerLoop1: do c=1,nCorner 

           depRate = zero

           do g=1,Groups
             depRate = depRate + Set% Psi(g,c0+c,angle)*  &
                                (Mat% Siga(g0+g,zone) + Mat% Sigs(g0+g,zone))
           enddo

           RadT% RadiationForce(:,c0+c) = RadT% RadiationForce(:,c0+c) + constant*  &
                                          wOmega(:)*depRate*Geom% Volume(c0+c) 


         enddo CornerLoop1

       enddo AngleLoop1

     enddo ZoneLoop1

   else

     ZoneLoop2: do zone=1,nzones

       nCorner = Geom% numCorner(zone)
       c0      = Geom% cOffSet(zone)

       RadT% RadiationFlux(:,:,zone) = zero
       RadT% radEnergy(zone) = zero
       RadT% EddingtonTensorDiag(:,zone) = zero

       AngleLoop2: do angle=1,Set% NumAngles

         wOmega(:) =  ASet% Weight(angle)*ASet% omega(:,angle)
         wOmegaOmega(:) = wOmega(:)*ASet% omega(:,angle)

         CornerLoop2: do c=1,nCorner

           VolFrac = Geom% Volume(c0+c)/Geom% VolumeZone(zone)

           do g=1,Groups
             RadT% radEnergy(zone) = RadT% radEnergy(zone) + &
                                             ASet% Weight(angle)*VolFrac* &
                                             Set% Psi(g,c0+c,angle)
             RadT% RadiationFlux(:,g,zone) = RadT% RadiationFlux(:,g,zone) +  &
                                             wOmega(:)*VolFrac*               &
                                             Set% Psi(g,c0+c,angle)
             RadT% EddingtonTensorDiag(:,zone) = RadT% EddingtonTensorDiag(:,zone) + &
                                             wOmegaOmega(:)*VolFrac* &
                                             Set% Psi(g,c0+c,angle)
           enddo

         enddo CornerLoop2

       enddo AngleLoop2

     enddo ZoneLoop2

   endif


   return
   end subroutine setRadiationMoments 


