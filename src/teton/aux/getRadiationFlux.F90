!***********************************************************************
!                       Last Update:  03/2012, PFN                     *
!                                                                      *
!   getRadiationFlux    - Calculates the zone-average radiation        *
!                         flux vector.                                 * 
!                                                                      *
!***********************************************************************
 
   subroutine getRadiationFlux(zoneID, RadiationFlux) &
        BIND(C,NAME="teton_getradiationflux")

   USE ISO_C_BINDING
   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use RadIntensity_mod
   use SetData_mod

   implicit none 

!  Arguments

   integer(C_INT),    intent(in)    :: zoneID

   real(C_DOUBLE), intent(inout) :: RadiationFlux(Size%ndim,Size%ngr)

!  Local

   type(RadIntensity),     pointer  :: RadT
   type(SetData),          pointer  :: Set

   integer    :: setID
   integer    :: nSets
   integer    :: g
   integer    :: g0
   integer    :: Groups

!  Constants

!***********************************************************************
!  Compute the radiation flux                                          *
!***********************************************************************

   nSets   =  getNumberOfSets(Quad)

   RadiationFlux(:,:) = zero


   SetLoop: do setID=1,nSets

     RadT   => getRadIntensity(Quad, setID)
     Set    => getSetData(Quad, setID)

     g0     =  Set% g0
     Groups =  Set% Groups

     do g=1,Groups
       RadiationFlux(:,g0+g) = RadiationFlux(:,g0+g) + &
                               RadT% RadiationFlux(:,g,zoneID) 
     enddo

   enddo SetLoop



   return
   end subroutine getRadiationFlux 


