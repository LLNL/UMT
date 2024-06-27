!***********************************************************************
!                        Last Update:  04/2024, PFN                    *
!                                                                      *
!   CycleBreakerCorner - This routine breaks cycles in the mesh by     * 
!                        selecting a corner that will use some old     *
!                        (i.e. previous iterate) incident fluxes.      *
!                                                                      *
!***********************************************************************
   subroutine cycleBreakerCorner(ndone, meshCycles, nextCorner,   &
                                 addedCorners, need, cornerList,  &
                                 cycleList, nDSC, DSC, onCycleList) 

   use kind_mod
   use constant_mod
   use Size_mod

   implicit none

!  Arguments

   integer,    intent(in)          :: ndone
   integer,    intent(inout)       :: meshCycles 
   integer,    intent(inout)       :: nextCorner
   integer,    intent(inout)       :: addedCorners

   integer,    intent(inout)       :: need(Size%ncornr)
   integer,    intent(inout)       :: cornerList(Size%ncornr)
   integer,    intent(inout)       :: cycleList(Size%ncornr)
   integer,    intent(inout)       :: nDSC(Size%ncornr)
   integer,    intent(inout)       :: DSC(2*Size%maxcf,Size%ncornr)

   logical (kind=1), intent(inout) :: onCycleList(Size%ncornr)

!  Local Variables

   integer :: i
   integer :: ngraph
   integer :: nleft
   integer :: ncount
   integer :: stackindex
   integer :: ncornr
   integer :: c
   integer :: nBreaks

!  Dynamic

   integer, allocatable :: list(:)
   integer, allocatable :: cBreakList(:)
   integer, allocatable :: dfnum(:)
   integer, allocatable :: lowlink(:)
   integer, allocatable :: stack(:)
   integer, allocatable :: tempList(:)

   logical (kind=1), allocatable :: new(:)
   logical (kind=1), allocatable :: onstack(:)

!  Mesh Constants

   ncornr = Size% ncornr

!  Allocate arrays for the number of corners in the graph (= ncornr - ndone)

   ngraph = ncornr - ndone

   allocate( list(ngraph) )
   allocate( cBreakList(ngraph) )
   allocate( dfnum(ncornr) )
   allocate( lowlink(ncornr) )
   allocate( stack(ngraph) )
   allocate( tempList(ngraph) )

   allocate( new(ncornr) )
   allocate( onstack(ncornr) )

!  Initialize arrays and counters

   new(:)     = .TRUE. 
   onstack(:) = .FALSE. 

   nBreaks    = 0
   ncount     = 0
   stackindex = 0

   stack(:)   = 0

!  Make a list of all remaining corners 

   nleft = 0

   do c=1,ncornr 
     if (need(c) == 0) then
       new(c)    = .FALSE. 
     else
       nleft        = nleft + 1
       list(nleft) = c 
     endif
   enddo

   if (nleft /= ngraph) then
     call f90fatal("Miscount of remaining corners in CycleBreakerCorner")
   endif

!  Loop over the number of corners in the graph

   do i=1,ngraph

     c = list(i)

     if ( new(c) ) then

       call sccSearchCorner(c, ngraph, ncount, stackindex,        &
                            nBreaks, meshCycles, dfnum, lowlink,  &
                            need, stack, new, onstack, tempList,  &
                            cycleList, cBreakList, nDSC, DSC,     &
                            onCycleList)

     endif

   enddo


   if (nBreaks == 0) then 

     call f90fatal("CycleBreakerCorner: detection failed, no dependencies broken")

   else

     addedCorners = 0
     do i=1,nBreaks
       c = cBreakList(i)

       if (need(c) == 0) then
         nextCorner             = nextCorner   + 1
         addedCorners           = addedCorners + 1
         cornerList(nextCorner) = c 
       elseif (need(c) < 0) then
         call f90fatal("CycleBreakerCorner, need < 0")
       endif
     enddo

     if (addedCorners == 0) then
       call f90fatal("Cycles found, but not broken in CycleBreakerCorner")
     endif

   endif

!  Release memory

   deallocate( list )
   deallocate( cBreakList )
   deallocate( dfnum )
   deallocate( lowlink )
   deallocate( stack )
   deallocate( tempList )
   deallocate( new )
   deallocate( onstack )


   return
   end subroutine cycleBreakerCorner
 
