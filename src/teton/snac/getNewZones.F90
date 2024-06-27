!***********************************************************************
!                        Last Update:  05/2023, PFN                    *
!                                                                      *
!   getNewZones   - This routine creates a list of starting zones.     *
!                   These zones are on the boundary of the grid        *
!                   and require no incident fluxes except from         *
!                   boundary conditions.  There may be situations      *
!                   where no starting zones can be found; in this      *
!                   situation, we are forced to use some old           *
!                   information to get started.                        * 
!                                                                      *
!***********************************************************************
   subroutine getNewZones(newZones, meshCycles, needZ, listZone,  &
                          cycleList, exitFace, onCycleList)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Arguments

   integer,    intent(inout)       :: newZones 
   integer,    intent(inout)       :: meshcycles 

   integer,    intent(inout)       :: needZ(Size%nzones)
   integer,    intent(inout)       :: listZone(Size%nzones) 
   integer,    intent(inout)       :: cycleList(Size%ncornr)

   logical (kind=1), intent(inout) :: exitFace(Size%maxFaces,Size%nzones)
   logical (kind=1), intent(inout) :: onCycleList(Size%nzones)

!  Local Variables

   integer :: c
   integer :: zone
   integer :: nzones 
   integer :: zoneID 
   integer :: zoneOpp
   integer :: minNeed
   integer :: face
   integer :: faceOpp
   integer :: nFaces

!  Mesh Constants

   nzones = Size% nzones

!  Create a list of zones to begin the sweep 

   newZones = 0

   ZoneLoop: do zone=1,nzones
     if (needZ(zone) == 0) then
       newZones           = newZones + 1
       listZone(newZones) = zone 
     endif
   enddo ZoneLoop


   if (newZones == 0) then

!  If no seeds were found, find a zone on the boundary that requires
!  the fewest incident fluxes 

     minNeed = nzones
     zoneID  = 0

     BoundaryZoneLoop: do zone=1,nzones
       if ( Geom% BoundaryZone(zone) ) then
         if (needZ(zone) < minNeed) then
           zoneID  = zone
           minNeed = needZ(zone)
         endif
       endif
     enddo BoundaryZoneLoop

     newZones      = 1
     listZone(1)   = zoneID
     needZ(zoneID) = 0
     nFaces        = Geom% zoneFaces(zoneID)

     do face=1,nFaces
       if ( .not. exitFace(face,zoneID) ) then
         zoneOpp = Geom% zoneOpp(face,zoneID)
         faceOpp = Geom% faceOpp(face,zoneID)

         if (zoneOpp > 0) then

           do c=1,Geom% numCorner(zoneOpp)
             meshCycles            = meshCycles + 1
             cycleList(meshCycles) = Geom% cOffSet(zoneOpp) + c
           enddo

           exitFace(faceOpp,zoneOpp) = .FALSE.
           onCycleList(zoneOpp)      = .TRUE.
         endif
       endif
     enddo

   endif

!  Error Check

   if (newZones == 0) then 
     call f90fatal("No starting zones found in getNewZones!")
   endif



   return
   end subroutine getNewZones 

