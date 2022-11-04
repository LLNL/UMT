!***********************************************************************
!                        Version 1:  10/2017, PFN                      *
!                                                                      *
!   setOppositeFace - This routine creates a list of opposite faces    * 
!                     for each zone. The faceIDs are local, such that: *
!                     1 <= faceID <= nFaces.                           *
!                                                                      *
!***********************************************************************

   subroutine setOppositeFace() BIND(C,NAME="teton_setoppositeface") 

!  Include
   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Local 

   integer :: face
   integer :: face2
   integer :: nFaces
   integer :: nFaces2
   integer :: zone
   integer :: nZones 
   integer :: zoneOpp

   nZones = Size% nZones 

!  First set boundary element numbers

   ZoneLoop: do zone=1,nZones

     Geom% faceOpp(:,zone) = 0
     nFaces = Geom% zoneFaces(zone)

     FaceLoop: do face=1,nFaces
       zoneOpp = Geom% zoneOpp(face,zone)

       if (zoneOpp < 0) then
         Geom% faceOpp(face,zone) = -1
         Geom% BoundaryZone(zone) = .TRUE.
       else
         nFaces2 = Geom% zoneFaces(zoneOpp)
         do face2=1,nFaces2
           if (Geom% zoneOpp(face2,zoneOpp) == zone) then
             Geom% faceOpp(face,zone) = face2
           endif
         enddo
       endif

       if ( Geom% faceOpp(face,zone) == 0 ) then
         call f90fatal("Inconsistent definition of opposite zone/face")
       endif

     enddo FaceLoop

   enddo ZoneLoop


   return
   end subroutine setOppositeFace 


