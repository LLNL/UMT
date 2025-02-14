!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   InitTeton - Calculates an initial radiation distribution.  The     *
!            energy density is a*T^4 where T is the initial radiation  *
!            temperature.  The energy is placed in a Planck spectrum   *
!            with an isotropic angular distribution.  The electron     *
!            temperature is also set.                                  *
!                                                                      *
!***********************************************************************
 
   subroutine InitTeton(erad) BIND(C,NAME="teton_initteton_new")

   USE ISO_C_BINDING
   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use RadIntensity_mod
   use Material_mod
   use QuadratureList_mod
   use SetData_mod

!  From PhysicsUtils library
   use emissionutils

   implicit none

!  Arguments

   real(C_DOUBLE), intent(inout) :: erad

!  Local

   type(SetData),       pointer  :: Set

   integer    :: angle
   integer    :: NumAngles
   integer    :: ngr
   integer    :: zone
   integer    :: nzones
   integer    :: c
   integer    :: c0
   integer    :: nCorner
   integer    :: setID
   integer    :: nSets
   integer    :: g
   integer    :: g0
   integer    :: Groups

   real(adqt) :: efloor
   real(adqt) :: wtiso
   real(adqt) :: geometryFactor
   real(adqt) :: Tr
   real(adqt) :: ac
   real(adqt) :: kb

!  Stack Arrays

   real(adqt) :: gnu(Size%ngr+1)
   real(adqt) :: planck(Size%ngr)

!  Constants

   efloor         = zero
   erad           = zero
   ngr            = Size% ngr
   nzones         = Size% nzones
   wtiso          = Size% wtiso
   geometryFactor = getGeometryFactor(Size)
   nSets          = getNumberOfSets(Quad)

   gnu(:)  = getEnergyGroups(Quad,ngr)

!  Compute the fraction of the total emission in each energy group
!  The input for RTPLNK is (h*nu)/(k*Te).
 
   ZoneLoop: do zone=1,nzones

     nCorner =  Geom% numCorner(zone)
     c0      =  Geom% cOffSet(zone)

!  Initialize the corner temperature

     do c=1,nCorner
       Mat% Tec(c0+c) = Mat% Tez(zone)
     enddo 

!  T4 has units of energy/area/time 

     Tr    = Mat% Trz(zone)
     ac    = speed_light*rad_constant
     kb    = one

!  Compute hnu/kt at upper energy boundary

     call integrateBlackBodyGroups(Tr,kb,ac,ngr,gnu,planck)

     SetLoop: do setID=1,nSets

       Set       => getSetData(Quad, setID)
       Groups    =  getNumberOfGroups(Quad, setID)
       NumAngles =  getNumberOfAngles(Quad, setID)
       g0        =  Set% g0

!  Loop over all angles in group

       do angle=1,NumAngles
         do c=1,nCorner
           do g=1,Groups
             Set% Psi(g,c0+c,angle) = max(wtiso*planck(g0+g),efloor)
           enddo
         enddo
       enddo

     enddo SetLoop

!    Initialize the radiation energy denisty

     planck = planck/speed_light

     do g=1,ngr
       Rad% RadEnergyDensity(zone,g)      = max(planck(g),efloor)
       erad = erad + Geom% VolumeZone(zone)*max(planck(g),efloor)
     enddo

   enddo ZoneLoop
 

   erad = geometryFactor*erad

   call MPIAllReduce(ERAD, "sum", MY_COMM_GROUP)

 
   return
   end subroutine InitTeton 

