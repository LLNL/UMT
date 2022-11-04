#include "macros.h"
!***********************************************************************
!                        Version 2:  01/98, PFN                        *
!                                                                      *
!   findReflectedAngles - This routine perfroms two functions:         *
!                                                                      *
!            1. Checks that all reflecting faces are in the same       *
!               plane (a requirement)                                  *
!            2. For each reflecting boundary, it finds the             *
!               reflected angle for each incident angle                *
!                                                                      *
!***********************************************************************

   subroutine findReflectedAngles(aSetID) 

   use kind_mod
   use constant_mod
   use Size_mod
   use BoundaryList_mod
   use Boundary_mod
   use QuadratureList_mod
   use Geometry_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,         intent(in) :: aSetID

!  Local Variables

   type(AngleSet),  pointer    :: ASet
   type(Boundary),  pointer    :: BdyT

   integer    :: d, ndim 
   integer    :: reflID, ib, nReflecting, nBdyElem
   integer    :: Angle
   integer    :: ReflAngle
   integer    :: c, c0, nCorner
   integer    :: zone
   integer    :: NumAngles

   real(adqt) :: eps, tol, A_mag, delta_A
   real(adqt) :: OmegaDotA
   real(adqt) :: A_set(Size%ndim)
   real(adqt) :: Area(Size%ndim)

!  Constants

   parameter (tol=1.0d-6)
   parameter (eps=1.d-15)

   ndim        =  Size%ndim
   nReflecting =  getNumberOfReflecting(RadBoundary)

   ASet        => getAngleSetData(Quad, aSetID)
   NumAngles   =  ASet% NumAngles

!  Check if all reflecting faces are in the same plane 

   ReflectingLoop: do reflID=1,nReflecting

     BdyT     => getReflecting(RadBoundary, reflID)
     nBdyElem =  getNumberOfBdyElements(BdyT)

     A_set(:) =  BdyT% A_bdy(:,1)
     A_mag    =  DOT_PRODUCT( A_set(:),A_set(:) ) 
     A_set(:) =  A_set(:)/sqrt(A_mag)

     do ib=1,nBdyElem

       Area(:) = BdyT% A_bdy(:,ib)
       A_mag   = DOT_PRODUCT( Area(:),Area(:) ) 
       Area(:) = Area(:)/sqrt(A_mag)

       delta_A = zero

       do d=1,ndim
          delta_A = delta_A + abs(A_set(d) - Area(d))
       enddo

       if ( delta_A > tol ) then
         zone    = BdyT% BdyToZone(ib)
         nCorner = Geom% numCorner(zone)
         c0      = Geom% cOffSet(zone)

         write(6,100) Size% myRankInGroup,zone

         if (ndim == 2) then
           do c=1,nCorner
             write(6,200) Geom% px(1,c0+c), Geom% px(2,c0+c)
           enddo
         else if (ndim == 3) then
           do c=1,nCorner
             write(6,300) Geom% px(1,c0+c), Geom% px(2,c0+c), Geom% px(3,c0+c)
           enddo
         endif

 100     format("Error occured on process ",i6,2x,"local zoneID = ",i8,  &
                " with coordinates:")

 200     format("r = ",1pe14.6,2x,"z = ",1pe14.6)
 300     format("x = ",1pe14.6,2x,"y = ",1pe14.6,2x,"z = ",1pe14.6)

         write(6,*) "findReflectedAngles: Not all faces in reflecting plane. Check for"
         write(6,*) "an inverted or degenerate zone (e.g., a zone in which two or more"
         write(6,*) "vertices are the same) on a reflecting boundary.  Teton does not"
         write(6,*) "support degenerate zones. If you get this message at startup, you"
         write(6,*) "likely have more than one plane in a reflecting BC.  Each plane of"
         write(6,*) "reflection needs its own BC."

         TETON_FATAL("findReflectedAngles: not all faces in BC set are planar.")

       endif

     enddo

!  Find the reflected angles

     Area(:) = BdyT% A_bdy(:,1)

     AngleLoop: do Angle=1,NumAngles

       OmegaDotA =  DOT_PRODUCT( ASet% omega(:,Angle),Area(:) )

       TestIncident: if (OmegaDotA < -eps) then

!  If OmegaDotA<0, Angle is incoming on this set. Here we test OmegaDotA<-eps
!  to account for roundoff errors if the direction is parallel to
!  the reflecting surface (e.g. in Lobatto quadratures). The routine
!  SNMREF computes the angle, MREF, that reflects onto Angle.
!  It also computes a multiplier, cosrat, that makes our
!  reflection algorithm conservative (i.e., net flux=0).
!  AREA contains the components of the outward normal for this boundary set.

         call reflectAxis(ndim, Angle, NumAngles, ReflAngle,  &
                          ASet% omega, Area)

       else

         ReflAngle = -1

       endif TestIncident

       call setReflectedAngle(ASet, reflID, Angle, ReflAngle)

     enddo AngleLoop

   enddo ReflectingLoop


   return
   end subroutine findReflectedAngles 


