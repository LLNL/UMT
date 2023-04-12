!***********************************************************************
!                        Version 1:  02/06, PFN                        *
!                                                                      *
!   setTetonZone1D - This routine sets geometry information in a zone  *
!                    data structure for a 1D spatial mesh.             *
!                                                                      *
!***********************************************************************

   subroutine setTetonZone1D(zoneID, numBCTotal, BCZoneID) &
        BIND(C,NAME="teton_setzone1d")

!  Include
   use ISO_C_BINDING
   use constant_mod
   use flags_mod
   use kind_mod
   use Size_mod
   use Geometry_mod
   use BoundaryList_mod
   use Boundary_mod


   implicit none

!  Arguments

   integer(C_INT),    intent(in)    :: zoneID
   integer(C_INT),    intent(in)    :: numBCTotal
   integer(C_INT),    intent(in)    :: BCZoneID(numBCTotal)

!  Local 

   integer          :: c,c0,b0
   integer          :: bcID

!  Set the zone data structures

   c0 = 2*(zoneID - 1) 

   Geom% cOffSet(zoneID)   = c0
   Geom% numCorner(zoneID) = 2
   Geom% zoneFaces(zoneID) = 2

   ! Assuming only 2 zone faces per zone,
   ! indexed consistently with local corners

   do c=1,2
     Geom% CToFace(1,c0+c) = c
   enddo

   if (zoneID == 1 .or. zoneID == Size% nzones) then
     Geom% BoundaryZone(zoneID) = .TRUE.

     ! BCZoneID is a global boundary condition array
     ! Loop through and figure out which bcID corresponds to this zone
     do bcID=1,numBCTotal
       if (BCZoneID(bcID) == zoneID) then
         Bdy => getBoundary(RadBoundary, bcID)
         Bdy% BdyToZone(1) = zoneID
         b0 = getFirstBdyElement(Bdy)- 1

         if (zoneID == 1) then
           RadBoundary% innerBdyID = bcID

           if ( getBCType(Bdy) == bcType_shared) then
             RadBoundary% inner1DBdyShared = .TRUE.
             RadBoundary% innerSharedID    = min(bcID, Size%ncomm) 
           else
             if ( getBCType(Bdy) == bcType_refl) then
               RadBoundary% inner1DBdyReflect = .TRUE.
             endif
             ! Fill out Geom% cFP unless it is a shared boundary:
             ! I am assuming there is only one corner for each boundary:
             Geom% cFP(1,c0+1) = Size% ncornr+b0+1
           endif

         elseif (zoneID == Size% nzones) then
           RadBoundary% outerBdyID = bcID

           if ( getBCType(Bdy) == bcType_shared) then
             RadBoundary% outer1DBdyShared = .TRUE.
             RadBoundary% outerSharedID    = min(bcID, Size%ncomm)
           else
             if ( getBCType(Bdy) == bcType_refl) then
               RadBoundary% outer1DBdyReflect = .TRUE.
             endif
             ! Fill out Geom% cFP unless it is a shared boundary:
             ! I am assuming there are only two possible boundary corners
             Geom% cFP(1,c0+2) = Size% ncornr+b0+1
           endif

         endif
       endif
     enddo

   else
     Geom% BoundaryZone(zoneID) = .FALSE.
   endif

   !! Surface edit needs these:
   do c = 1,Geom% numCorner(zoneID)
     Geom% cToZone(c0+c)      = zoneID
     Geom% nCFacesArray(c0+c) = 1
   enddo
   ! Outward pointing normals:
   Geom% A_fp(1,1,c0+1) = -one
   Geom% A_fp(1,1,c0+2) = one
   ! Fill out Geom% cFP for non-boundary zones:
   if ( zoneID /= 1 ) then
     Geom% cFP(1,c0+1)  = c0
   endif
   if ( zoneID /= Size% nzones ) then
     Geom% cFP(1,c0+2)  = c0+3
   endif


   return
   end subroutine setTetonZone1D 


