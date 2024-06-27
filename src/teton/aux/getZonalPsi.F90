#include "macros.h"
!***********************************************************************
!                        Version 1:  10/2016, PFN                      *
!                                                                      *
!   getZonalPsi                                                        *
!                                                                      *
!***********************************************************************

   subroutine getZonalPsi(numAngles, Psi) &
        BIND(C,NAME="teton_getzonalpsi")

   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use constant_mod
   use Geometry_mod
   use QuadratureList_mod
   use AngleSet_mod
   use SetData_mod

   implicit none

!  Arguments

   integer(C_INT), intent(in)     :: numAngles
! TODO is this the ordering we want? Psi is (group, zone, angle) in SetData
!  Balance between ease for host code and Teton's layout
! Note that index ordering in Fortran is backwards compared to C++!
   real(C_DOUBLE), intent(inout)  :: Psi(Size%nzones,Size%ngr,numAngles)

!  Local

   type(AngleSet), pointer  :: ASet
   type(SetData), pointer  :: Set

   integer :: angleSetID, setID
   integer :: numAnglesInternal
   integer :: numAngleSets
   integer :: numSets

   integer :: angle0, group0
   integer :: angle, group
   integer :: numAnglesInSet
   integer :: numGroupsInSet

   integer :: c0, nCorner, zone

   numAngleSets = getNumberOfAngleSets(Quad)
! Sanity check:
   numAnglesInternal = 0
   do angleSetID = 1,numAngleSets
      ASet => getAngleSetData(Quad,angleSetID)
      numAnglesInternal = numAnglesInternal + ASet%NumAngles
   enddo
   TETON_ASSERT(numAngles == numAnglesInternal, "numAngles given to getZonalPsi does not match Teton's internal total number of angles")

   numSets = getNumberOfSets(Quad)
   do setID = 1,numSets
      Set => getSetData(Quad, setID)

      angle0 = Set%angle0
      group0 = Set%g0
      numAnglesInSet = Set%NumAngles
      numGroupsInSet = Set%Groups

      ! Copy over data:
      do group = 1,numGroupsInSet
        do angle = 1,numAnglesInSet
           do zone = 1,Size%nzones
             nCorner = Geom% numCorner(zone)
             c0      = Geom% cOffSet(zone)

             Psi(zone,group0+group,angle0+angle) = &
               dot_product(Set%Psi(group,c0+1:c0+nCorner,angle), Geom%Volume(c0+1:c0+nCorner))/ &
                 Geom%VolumeZone(zone)
           enddo
         enddo
      enddo

   enddo ! Set Loop

   return
   end subroutine getZonalPsi
