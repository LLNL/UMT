!***********************************************************************
!                        Version 1:  04/02, PFN                        *
!                                                                      *
!   FINDSEEDS    - This routine creates a list of starting points or   *
!                  "seeds" for the grid sweep.  The seeds are on the   *
!                  boundary of the grid and require no incident        *
!                  fluxes except from boundary conditions.  There may  *
!                  be situations where no seeds can be found; this     *
!                  will occur if there is a mesh cycle right at the    *
!                  boundary.  In this situation, we are forced to use  *
!                  some old information to get started.                * 
!                                                                      *
!   Input:                                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!***********************************************************************
   subroutine findseeds(NSEED, MESHCYCLES, needZ, listZone,  &
                        cycleList, exitFace, onCycleList)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Arguments

   integer,    intent(inout)       :: nseed
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

!  Create a list of zone "seeds"

   nseed = 0

   ZoneLoop: do zone=1,nzones
     if (needZ(zone) == 0) then
       nseed           = nseed + 1
       listZone(nseed) = zone 
     endif
   enddo ZoneLoop


   if (nseed == 0) then

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

     nseed         = 1
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

   if (nseed == 0) then 
     call f90fatal("No seeds found in FINDSEEDS!")
   endif



   return
   end subroutine findseeds 

