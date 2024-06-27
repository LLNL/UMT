!***********************************************************************
!                        Last Update:  04/2024, PFN                    *
!                                                                      *
!  sccSearchCorner- This recursive routine search the dependency graph *
!                   for strongly-connected components (SCC).           *
!                                                                      *
!***********************************************************************
   recursive subroutine sccSearchCorner(c, ngraph, ncount, stackindex,       &
                                        nBreaks, meshCycles, dfnum, lowlink, &
                                        need, stack, new, onstack, tempList, &
                                        cycleList, cBreakList, nDSC, DSC,    &
                                        onCycleList)

   use kind_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Arguments

   integer,    intent(in)          :: c 
   integer,    intent(in)          :: ngraph
   integer,    intent(inout)       :: ncount
   integer,    intent(inout)       :: stackindex 
   integer,    intent(inout)       :: nBreaks
   integer,    intent(inout)       :: meshCycles 

   integer,    intent(inout)       :: dfnum(Size%ncornr)
   integer,    intent(inout)       :: lowlink(Size%ncornr)
   integer,    intent(inout)       :: need(Size%ncornr)
   integer,    intent(inout)       :: stack(ngraph)
   integer,    intent(inout)       :: cBreakList(ngraph)
   integer,    intent(inout)       :: tempList(ngraph)
   integer,    intent(inout)       :: cycleList(Size%ncornr)
   integer,    intent(inout)       :: nDSC(Size%ncornr)
   integer,    intent(inout)       :: DSC(2*Size%maxcf,Size%ncornr)

   logical (kind=1), intent(inout) :: new(Size%ncornr)
   logical (kind=1), intent(inout) :: onstack(Size%ncornr)
   logical (kind=1), intent(inout) :: onCycleList(Size%ncornr)

!  Local Variables

   integer    :: i
   integer    :: cExit
   integer    :: cyclesize
   integer    :: cBreak
   integer    :: lowlinkC
   integer    :: cface
   integer    :: nCFaces
   integer    :: c0
   integer    :: zone

!  Start the search procedure

   ncount     = ncount + 1
   dfnum(c)   = ncount
   lowlink(c) = ncount
   new(c)     = .FALSE.

!  Put current "corner" on the stack

   stackindex        = stackindex + 1
   stack(stackindex) = c 
   onstack(c)        = .TRUE. 

!  Loop over all downstream corners that have not been completed 

   DownStreamC: do i=1,nDSC(c) 

     cExit = DSC(i,c)

     if ( new(cExit) ) then

       call sccSearchCorner(cExit, ngraph, ncount, stackindex,    &
                            nBreaks, meshCycles, dfnum, lowlink,  &
                            need, stack, new, onstack, tempList,  &
                            cycleList, cBreakList, nDSC, DSC,     &
                            onCycleList)

       if (lowlink(cExit) < lowlink(c)) then
         lowlink(c) = lowlink(cExit)
       endif

     else

       if (dfnum(cExit) < dfnum(c) .and.  &
           onstack(cExit)          .and.  &
           lowlink(cExit) < lowlink(c)) then

         lowlink(c) = lowlink(cExit)
       endif
 
     endif

   enddo DownStreamC

!  Cycle Check

   CheckCycle: if (lowlink(c) == dfnum(c)) then

     cExit          = stack(stackindex)
     stackindex     = stackindex - 1
     onstack(cExit) = .FALSE. 

     DetectCycle: if (cExit /= c) then

       cyclesize  = 0

       do while (cExit /= c)
         cyclesize           = cyclesize + 1
         tempList(cyclesize) = cExit 

         cExit               = stack(stackindex)
         stackindex          = stackindex - 1
       enddo

       cyclesize             = cyclesize + 1
       tempList(cyclesize)   = cExit 
       onstack(tempList(1))  = .TRUE.

!***********************************************************************
!  Now break all connections of corners on the stack to the lowest     *
!  link.                                                               *
!***********************************************************************

       lowlinkC = tempList(cyclesize)
       zone     = Geom% CToZone(lowlinkC)
       c0       = Geom% cOffSet(zone)

!  Loop over all neighbors for this corner and find the ones on the stack

       if (Size% ndim == 2) then
         nCFaces = 2 
       else
         nCFaces = Geom% nCFacesArray(lowlinkC)
       endif

       CornerFaceLoop: do cface=1,nCFaces

!  Corners in the same zone

         cBreak = Geom% cEZ(cface,lowlinkC) + c0 

         if ( onstack(cBreak) .and. (.not. onCycleList(cBreak))) then

           do i=1,nDSC(cBreak)
             cExit = DSC(i,cBreak)

             if (cExit == lowlinkC) then
               DSC(i,cBreak)         = Size% ncornr + 1
               meshCycles            = meshCycles + 1
               cycleList(meshCycles) = cBreak
               onCycleList(cBreak)   = .TRUE.
               need(lowlinkC)        = need(lowlinkC) - 1
             endif
           enddo

         endif

!  Corners in neighboring zones

         cBreak = Geom% cFP(cface,lowlinkC)

         if ( cBreak <= Size% ncornr ) then

           if ( onstack(cBreak) .and. (.not. onCycleList(cBreak))) then

             do i=1,nDSC(cBreak)
               cExit = DSC(i,cBreak)

               if (cExit == lowlinkC) then
                 DSC(i,cBreak)         = Size% ncornr + 1
                 meshCycles            = meshCycles + 1
                 cycleList(meshCycles) = cBreak
                 onCycleList(cBreak)   = .TRUE.
                 need(lowlinkC)        = need(lowlinkC) - 1
               endif
             enddo

           endif
         endif

       enddo CornerFaceLoop

       if (need(lowlinkC) == 0) then
         nBreaks             = nBreaks + 1
         cBreakList(nBreaks) = lowlinkC
       endif

!  Reset the stack

       do i=1,cyclesize
         onstack( tempList(i) ) = .FALSE.
       enddo

     endif DetectCycle

   endif CheckCycle



   return
   end subroutine sccSearchCorner 
 
