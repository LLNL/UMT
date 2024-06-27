!***********************************************************************
!                        Last Update:  04/2024, PFN                    *
!                                                                      *
!   CycleBreakerZone - This routine breaks cycles in the mesh by       *
!                      selecting a zone that will use some old         *
!                      (i.e. previous iterate) incident fluxes.        *
!                                                                      *
!***********************************************************************
   subroutine cycleBreakerZone(ndoneZ, MESHCYCLES, nextZone,  &
                               addedZones,needZ, listZone, cycleList, &
                               exitFace, onCycleList) 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Arguments

   integer,    intent(in)          :: ndoneZ

   integer,    intent(inout)       :: meshcycles 
   integer,    intent(inout)       :: nextZone
   integer,    intent(inout)       :: addedZones

   integer,    intent(inout)       :: needZ(Size%nzones)
   integer,    intent(inout)       :: listZone(Size%nzones)
   integer,    intent(inout)       :: cycleList(Size%ncornr)

   logical (kind=1), intent(inout) :: exitFace(Size%maxFaces,Size%nzones)
   logical (kind=1), intent(inout) :: onCycleList(Size%nzones)

!  Local Variables

   integer :: i,ngraph,nleft,ncount,stackindex,nzones,zone
   integer :: nBreaks

!  Dynamic

   integer, allocatable :: listZ(:)
   integer, allocatable :: zoneBreakList(:)
   integer, allocatable :: dfnum(:)
   integer, allocatable :: lowlink(:)
   integer, allocatable :: stack(:)
   integer, allocatable :: tempList(:)

   logical (kind=1), allocatable :: new(:)
   logical (kind=1), allocatable :: onstack(:)

!  Mesh Constants

   nzones = Size%nzones

!  Allocate arrays for the number of zones in the graph (= nzones - ndone)

   ngraph = nzones - ndoneZ

   allocate( listZ(ngraph) )
   allocate( zoneBreakList(ngraph) )
   allocate( dfnum(nzones) )
   allocate( lowlink(nzones) )
   allocate( stack(ngraph) )
   allocate( tempList(ngraph) )

   allocate( new(nzones) )
   allocate( onstack(nzones) )

!  Initialize arrays and counters

   new(:)     = .TRUE. 
   onstack(:) = .FALSE. 

   nBreaks    = 0
   ncount     = 0
   stackindex = 0

   stack(:)   = 0

!  Make a list of all remaining zones 

   nleft = 0

   do zone=1,nzones 
     if (needZ(zone) == 0) then
       new(zone)    = .FALSE. 
     else
       nleft        = nleft + 1
       listZ(nleft) = zone 
     endif
   enddo

   if (nleft /= ngraph) then
     call f90fatal("Miscount of remaining zones in CycleBreakerZone")
   endif

!  Loop over the number of zones in the graph

   do i=1,ngraph

     zone = listZ(i)

     if ( new(zone) ) then

       call sccSearchZone(zone, ngraph, ncount, stackindex,     &
                          nBreaks, meshCycles, dfnum, lowlink,  &
                          needZ, stack, new, onstack, exitFace, &
                          tempList, cycleList, zoneBreakList,   &
                          onCycleList)

     endif

   enddo


   if (nBreaks == 0) then 

     call f90fatal("CycleBreakerZone: detection failed, no dependencies broken")

   else

     addedZones   = 0

     do i=1,nBreaks
       zone = zoneBreakList(i)

       if (needZ(zone) == 0) then
         nextZone           = nextZone   + 1
         addedZones         = addedZones + 1
         listZone(nextZone) = zone
       elseif (needZ(zone) < 0) then
         call f90fatal("CycleBreakerZone, needZ < 0")
       endif
     enddo

     if (addedZones == 0) then
       call f90fatal("CycleBreakerZone found, but not broken")
     endif

   endif

!  Release memory

   deallocate( listZ )
   deallocate( zoneBreakList )
   deallocate( dfnum )
   deallocate( lowlink )
   deallocate( stack )
   deallocate( tempList )
   deallocate( new )
   deallocate( onstack )


   return
   end subroutine cycleBreakerZone
 
