!***********************************************************************
!                        Last Update:  04/2024, PFN                    *
!                                                                      *
!   sccSearchZone- This recursive routine search the dependency graph  *
!                  for strongly-connected components (SCC).            *
!                                                                      *
!***********************************************************************
   recursive subroutine sccSearchZone(zone, ngraph, ncount, stackindex, &
                                      nBreaks, meshCycles, dfnum,  &
                                      lowlink, needZ, stack, new,  &
                                      onstack, exitFace, tempList, &
                                      cycleList, zoneBreakList,   &
                                      onCycleList)

   use kind_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Arguments

   integer,    intent(in)          :: zone, ngraph

   integer,    intent(inout)       :: ncount
   integer,    intent(inout)       :: stackindex 
   integer,    intent(inout)       :: nBreaks
   integer,    intent(inout)       :: meshCycles 

   integer,    intent(inout)       :: dfnum(Size%nzones)
   integer,    intent(inout)       :: lowlink(Size%nzones)
   integer,    intent(inout)       :: needZ(Size%nzones)
   integer,    intent(inout)       :: stack(ngraph)
   integer,    intent(inout)       :: zoneBreakList(ngraph)
   integer,    intent(inout)       :: tempList(ngraph)
   integer,    intent(inout)       :: cycleList(Size%ncornr)

   logical (kind=1), intent(inout) :: new(Size%nzones)
   logical (kind=1), intent(inout) :: onstack(Size%nzones)
   logical (kind=1), intent(inout) :: exitFace(Size%maxFaces,Size%nzones)
   logical (kind=1), intent(inout) :: onCycleList(Size%nzones)

!  Local Variables

   integer    :: i
   integer    :: c
   integer    :: zone2
   integer    :: cyclesize
   integer    :: faceBreak
   integer    :: zoneBreak
   integer    :: lowlinkZ
   integer    :: face
   integer    :: nFaces

!  Start the search procedure

   ncount        = ncount + 1
   dfnum(zone)   = ncount
   lowlink(zone) = ncount
   new(zone)     = .FALSE.

!  Put current "zone" on the stack

   stackindex        = stackindex + 1
   stack(stackindex) = zone 
   onstack(zone)     = .TRUE. 

!  Loop over all downstream zones that have not been completed 

   nFaces = Geom% zoneFaces(zone)

   FaceLoop: do face=1,nFaces 

     if ( exitFace(face,zone) ) then

       zone2 = Geom% zoneOpp(face,zone) 

       if (zone2 > 0) then

         if ( new(zone2) ) then

           call sccSearchZone(zone2, ngraph, ncount, stackindex,    &
                              nBreaks, meshCycles, dfnum, lowlink,  &
                              needZ, stack, new, onstack, exitFace, &
                              tempList, cycleList, zoneBreakList,   &
                              onCycleList)

           if (lowlink(zone2) < lowlink(zone)) then
             lowlink(zone) = lowlink(zone2)
           endif

         else

           if (dfnum(zone2) < dfnum(zone) .and.  &
               onstack(zone2)             .and.  &
               lowlink(zone2) < lowlink(zone)) then

             lowlink(zone) = lowlink(zone2)
           endif
 
         endif

       endif

     endif

   enddo FaceLoop

!  Cycle Check

   CheckCycle: if (lowlink(zone) == dfnum(zone)) then

     zone2          = stack(stackindex)
     stackindex     = stackindex - 1
     onstack(zone2) = .FALSE. 

     DetectCycle: if (zone2 /= zone) then

       cyclesize  = 0

       do while (zone2 /= zone)
         cyclesize           = cyclesize + 1
         tempList(cyclesize) = zone2 

         zone2               = stack(stackindex)
         stackindex          = stackindex - 1
       enddo

       cyclesize             = cyclesize + 1
       tempList(cyclesize)   = zone2
       onstack(tempList(1))  = .TRUE.

!***********************************************************************
!  Now break all connections of zones on the stack to the lowest       *
!  link.                                                               *
!***********************************************************************

       lowlinkZ = tempList(cyclesize)

!  Loop over all neighbors for this zone and find the ones on the stack

       nFaces = Geom% zoneFaces(lowlinkZ)

       FaceLoop2: do face=1,nFaces

         zoneBreak = Geom% zoneOpp(face,lowlinkZ)
         faceBreak = Geom% faceOpp(face,lowlinkZ)

         if (zoneBreak > 0) then

           if ( onstack(zoneBreak) .and. exitFace(faceBreak,zoneBreak) ) then

             if ( .not. onCycleList(zoneBreak) ) then
               do c=1,Geom% numCorner(zoneBreak)
                 meshCycles            = meshCycles + 1
                 cycleList(meshCycles) = Geom% cOffSet(zoneBreak) + c
               enddo
               onCycleList(zoneBreak)  = .TRUE.
             endif

             needZ(lowlinkZ)               = needZ(lowlinkZ) - 1
             exitFace(faceBreak,zoneBreak) = .FALSE.

             if (needZ(lowlinkZ) == 0) then
               nBreaks                = nBreaks    + 1
               zoneBreakList(nBreaks) = lowlinkZ
             endif
               
           endif

         endif

       enddo FaceLoop2

!  Reset the stack

       do i=1,cyclesize
         onstack( tempList(i) ) = .FALSE.
       enddo

     endif DetectCycle

   endif CheckCycle



   return
   end subroutine sccSearchZone 
 
