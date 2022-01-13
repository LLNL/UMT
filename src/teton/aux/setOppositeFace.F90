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
   use MeshData_mod
   use ZoneData_mod

   implicit none

!  Local 

   type(MeshData), pointer   :: M2

   integer :: face
   integer :: face2
   integer :: zone
   integer :: nZones 
   integer :: zoneOpp

   nZones = Size% nZones 

!  First set boundary element numbers

   ZoneLoop: do zone=1,nZones
     Z => getZoneData(Geom, zone)
     M => getMesh(Geom, zone)
     M% faceOpp(:) = 0

     FaceLoop: do face=1,M% nFaces
       zoneOpp = M% zoneOpp(face)

       if (zoneOpp < 0) then
         M% faceOpp(face) = -1
         Z% BoundaryZone  = .TRUE.
       else
         M2 => getMesh(Geom, zoneOpp)

         do face2=1,M2% nFaces
           if (M2% zoneOpp(face2) == zone) then
             M% faceOpp(face) = face2
           endif
         enddo
       endif

       if ( M% faceOpp(face) == 0 ) then
         call f90fatal("Inconsistent definition of opposite zone/face")
       endif

     enddo FaceLoop

   enddo ZoneLoop


   return
   end subroutine setOppositeFace 


