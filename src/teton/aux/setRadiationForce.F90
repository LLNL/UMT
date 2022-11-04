!***********************************************************************
!                        Last Update:  03/2012, PFN                    *
!                                                                      *
!   setRadiationForce   - Calculates the radiation force.              *
!                                                                      *
!***********************************************************************
 
   subroutine setRadiationForce() BIND(C,NAME="teton_setradiationforce")

   USE ISO_C_BINDING
   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use QuadratureList_mod
   use Geometry_mod
   use Material_mod
   use SetData_mod
   use AngleSet_mod
   use RadIntensity_mod

   implicit none 

!  Local

   type(SetData),  pointer  :: Set
   type(AngleSet), pointer  :: ASet

   integer    :: zSetID
   integer    :: nZoneSets
   integer    :: zone
   integer    :: setID
   integer    :: nSets
   integer    :: angle
   integer    :: c
   integer    :: g
   integer    :: g0
   integer    :: Groups

   real(adqt) :: constant
   real(adqt) :: depRate 
   real(adqt) :: wOmega(Size%ndim)

!***********************************************************************
!  Compute the radiation force on the matter                           *
!***********************************************************************

   nZoneSets = getNumberOfZoneSets(Quad)
   nSets     = getNumberOfSets(Quad)

   constant  = Size% radForceMultiplier*Size% geometryFactor/speed_light


!$omp  parallel do default(none) schedule(dynamic)                   &
!$omp& shared(Geom, Mat, Rad, Quad, nZoneSets, nSets, constant)      &
!$omp& private(zone, Set, ASet, g0, Groups, wOmega, depRate)

   ZoneSetLoop: do zSetID=1,nZoneSets

     do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
       Rad% RadiationForce(:,c) = zero
     enddo

     do setID=1,nSets

       Set  => getSetData(Quad, setID)
       ASet => getAngleSetFromSetID(Quad, setID)

       Groups = getNumberOfGroups(Quad, setID)
       g0     = Set% g0

       do angle=1,Set% NumAngles

         wOmega(:) = ASet% Weight(angle)*ASet% omega(:,angle)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           zone = Geom% CToZone(c)

           depRate = zero

           do g=1,Groups
             depRate = depRate + Set% Psi(g,c,angle)*  &
                                (Mat% Siga(g0+g,zone) + Mat% Sigs(g0+g,zone))
           enddo

           Rad% RadiationForce(:,c) = Rad% RadiationForce(:,c) + constant* &
                                      wOmega(:)*depRate*Geom% Volume(c)
         enddo

       enddo

     enddo

   enddo ZoneSetLoop

!$omp end parallel do


   return
   end subroutine setRadiationForce 


