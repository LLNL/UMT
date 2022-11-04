!***********************************************************************
!                       Last Update:  03/2012, PFN                     *
!                                                                      *
!   setRadiationFlux    - Calculates the zone-average radiation        *
!                         flux vector.                                 * 
!                                                                      *
!***********************************************************************
 
   subroutine setRadiationFlux() BIND(C,NAME="teton_setradiationflux")

   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use Geometry_mod
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

   real(adqt) :: VolFrac
   real(adqt) :: wOmega(Size%ndim)
   real(adqt) :: wOmegaOmega(Size%ndim)

!***********************************************************************
!  Compute the radiation flux                                          *
!***********************************************************************

   nZoneSets = getNumberOfZoneSets(Quad)
   nSets     = getNumberOfSets(Quad)


!$omp  parallel do default(none)  &
!$omp& shared(Geom, Rad, Quad, nZoneSets, nSets) &
!$omp& private(Set, ASet, g0, Groups, wOmega, wOmegaOmega, VolFrac) &
!$omp& schedule(dynamic)

   ZoneSetLoop: do zSetID=1,nZoneSets

     do zone=Geom% zone1(zSetID),Geom% zone2(zSetID)
       Rad% radEnergy(zone)             = zero
       Rad% RadiationFlux(:,:,zone)     = zero
       Rad% EddingtonTensorDiag(:,zone) = zero
     enddo

     do setID=1,nSets

       Set  => getSetData(Quad, setID)
       ASet => getAngleSetFromSetID(Quad, setID)

       Groups = getNumberOfGroups(Quad, setID)
       g0     = Set% g0

       do angle=1,Set% NumAngles

         wOmega(:)      = ASet% Weight(angle)*ASet% omega(:,angle)
         wOmegaOmega(:) = wOmega(:)*ASet% omega(:,angle)

         do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           zone = Geom% CToZone(c)

           VolFrac = Geom% Volume(c)/Geom% VolumeZone(zone)

           do g=1,Groups
             Rad% radEnergy(zone)             = Rad% radEnergy(zone) +             &
                                                ASet% Weight(angle)*VolFrac*       &
                                                Set% Psi(g,c,angle)
             Rad% RadiationFlux(:,g0+g,zone)  = Rad% RadiationFlux(:,g0+g,zone) +  &
                                                wOmega(:)*VolFrac*                 &
                                                Set% Psi(g,c,angle)
             Rad% EddingtonTensorDiag(:,zone) = Rad% EddingtonTensorDiag(:,zone) + &
                                                wOmegaOmega(:)*VolFrac*            &
                                                Set% Psi(g,c,angle)
           enddo

         enddo

       enddo


     enddo

   enddo ZoneSetLoop

!$omp end parallel do


   return
   end subroutine setRadiationFlux 


