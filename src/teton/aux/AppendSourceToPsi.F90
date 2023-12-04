#include "macros.h"

!***********************************************************************
!                       Last Update:  2023/01/10 BCY                   *
!                                                                      *
!   AppendSourceToPsi - Updates Psi from the previous time step to     *
!               to account for external fixed radiation sources.       *
!                                                                      *
!                                                                      *
!   We are basically adding q*c*dt to \psi_old to piggy-back off       *
!   the existing \psi_old/(c*dt) term in the sweep                     *
!                                                                      *
!                                                                      *
!***********************************************************************
 
   subroutine AppendSourceToPsi(numSourceZones, &
                                numSourceAngles, &
                                numSourceGroups, &
                                sourceZoneList,  &
                                sourceValues) BIND(C,NAME="teton_appendsourcetopsi")

   USE ISO_C_BINDING
   use kind_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use TimeStepControls_mod
   use Geometry_mod
   use QuadratureList_mod

   implicit none 

!  Arguments

   integer(C_INT), intent(in)    :: numSourceZones
   integer(C_INT), intent(in)    :: numSourceAngles
   integer(C_INT), intent(in)    :: numSourceGroups

   integer(C_INT), intent(in)    :: sourceZoneList(numSourceZones)
   ! This must be in units of energy per time
   ! That is, sourceValues*c*dt/zoneVolume must have the same units as \psi
   real(C_DOUBLE), intent(in)    :: sourceValues(numSourceGroups, numSourceAngles, numSourceZones)

!  Local

   type(SetData), pointer  :: Set
   type(AngleSet), pointer :: ASet

   integer    :: angle
   integer    :: angle0
   integer    :: polarAngle
   integer    :: c
   integer    :: c0
   integer    :: g0
   integer    :: nCornerInZone
   integer    :: z
   integer    :: zone
   integer    :: NumGroupsInSet
   integer    :: NumAnglesInSet
   integer    :: setID
   integer    :: nSets
   real(adqt) :: cdt
   real(adqt) :: volumeZone
   real(adqt) :: src_value(numSourceGroups)

   integer    :: angleSetID
   integer    :: numAngleSets
   integer    :: numAnglesInternal
   integer    :: numPolarAnglesInternal

!  Constants

   nSets   = getNumberOfSets(Quad)
   numAngleSets = getNumberOfAngleSets(Quad)
   cdt = speed_light*getRadTimeStep(DtControls)

! Some checks on quadrature/angle sizes:
   numAnglesInternal = 0
   do angleSetID = 1,numAngleSets
      ASet => getAngleSetData(Quad,angleSetID)
      numAnglesInternal = numAnglesInternal + ASet%NumAngles
   enddo
   TETON_VERIFY(Quad%QuadPtr(1)%TypeName /= 'lobatto', 'Lobatto quadratures not supported for general sources yet.')
   if (Quad%QuadPtr(1)%TypeName == 'levelsym') then
      TETON_VERIFY(numSourceAngles == 1 .or. numSourceAngles == numAnglesInternal, "For level symmetric, numSourceAngles must be 1 or # of total angles.")
   elseif (Quad%QuadPtr(1)%TypeName == 'product') then
      numPolarAnglesInternal = Quad%QuadPtr(1)%nPolarAngles
      TETON_VERIFY(numSourceAngles == 1 .or. numSourceAngles == numAnglesInternal .or. numSourceAngles == numPolarAnglesInternal, "numSourceAngles must be 1, # of polar angles, or # of total angles.")
   endif

   TETON_VERIFY(numSourceGroups == Size%ngr, "number of source groups must match ngr in Size")

   SetLoop: do setID=1,nSets

     Set       => getSetData(Quad, setID)

     NumGroupsInSet =  Set% Groups
     NumAnglesInSet =  Set% NumAngles
     g0             =  Set% g0    
     angle0         =  Set% angle0    

     AngleLoop: do angle=1,Set%NumAngles

        SourceZoneLoop: do z = 1,numSourceZones

           zone          = sourceZoneList(z)
           c0            = Geom%cOffset(zone)
           nCornerInZone = Geom%numCorner(zone)
           volumeZone    = Geom%VolumeZone(zone)*getGeometryFactor(Size)

           if (numSourceAngles == 1) then
              src_value(1:NumGroupsInSet) = sourceValues(g0+1:g0+NumGroupsInSet, 1, z)*Size%wtiso/volumeZone
           else if (numSourceAngles == numPolarAnglesInternal) then
              polarAngle = Set%polarAngle(angle)
              if (polarAngle < 1) then ! skip starting/finishing directions
                 cycle
              endif
              src_value(1:NumGroupsInSet) = sourceValues(g0+1:g0+NumGroupsInSet, polarAngle, z)*Size%wtiso*two/volumeZone
           else !if (numSourceAngles == numAnglesInternal) then
              src_value(1:NumGroupsInSet) = sourceValues(g0+1:g0+NumGroupsInSet, angle0+angle, z)/volumeZone
           endif
           src_value = cdt*src_value

           do c=c0+1,c0+nCornerInZone
              Set% Psi(:,c,angle) = Set% Psi(:,c,angle) + src_value(1:NumGroupsInSet)
           enddo

        enddo SourceZoneLoop
     enddo AngleLoop
   enddo SetLoop

   return
   end subroutine AppendSourceToPsi
