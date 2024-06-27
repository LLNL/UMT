!***********************************************************************
!                        Last Update:  04/2019, PFN                    *
!                                                                      *
!   getNewCorners - This routine creates a list of starting corners    *
!                   These corners are on the boundary of the grid      *
!                   and require no incident fluxes except from         *
!                   boundary conditions.  There may be situations      *
!                   where no starting corners can be found; in this    *
!                   situation, we are forced to use some old           *
!                   information to get started.                        * 
!                                                                      *
!***********************************************************************
   subroutine getNewCorners(newCorners, meshCycles, need,   &
                            cornerList, cycleList, nDSC, DSC)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Arguments

   integer, intent(inout)    :: newCorners 
   integer, intent(inout)    :: meshcycles 

   integer, intent(inout)    :: need(Size%ncornr)
   integer, intent(inout)    :: cornerList(Size%ncornr) 
   integer, intent(inout)    :: cycleList(Size%ncornr)
   integer, intent(inout)    :: nDSC(Size%ncornr)
   integer, intent(inout)    :: DSC(2*Size%maxcf,Size%ncornr)

!  Local Variables

   type(Boundary), pointer   :: BdyT

   integer :: i 
   integer :: c
   integer :: c0
   integer :: cez
   integer :: nCorner
   integer :: cface
   integer :: nCFaces
   integer :: n
   integer :: nBoundary
   integer :: b
   integer :: NumBdyElem
   integer :: zone

!  Create a list of corners to begin the sweep 

   nBoundary  = getNumberOfBoundaries(RadBoundary)
   newCorners = 0

   CornerLoop: do c=1,Size% ncornr
     if (need(c) == 0) then
       newCorners             = newCorners + 1
       cornerList(newCorners) = c 
     endif
   enddo CornerLoop


   if (newCorners == 0) then

!  If no corners were found, find a corner on the boundary that requires
!  only one incident flux. Also, add all corners on "ez" faces to the
!  cycle list to be conservative.

     BoundaryLoop: do n=1,nBoundary

       BdyT       => getBoundary(RadBoundary, n)
       NumBdyElem =  getNumberOfBdyElements(BdyT)

       do b=1,NumBdyElem
         c    = BdyT% BdyToC(b)
         zone = Geom% CToZone(c)
         c0   = Geom% cOffSet(zone)

         if (need(c) == 1) then
           newCorners             = 1
           cornerList(newCorners) = c
           need(c)                = 0

           if (Size% ndim == 3) then
             nCFaces = Geom% nCFacesArray(c)
           elseif (Size% ndim == 2) then
             nCFaces = 2
           endif

           do cface=1,nCFaces
             cez = c0 + Geom% cEZ(cface,c)

!  If corner "c" is downstream of "cez", add cez to the cycle list
!  and break the connection

             do i=1,nDSC(cez)
               if ( DSC(i,cez) == c ) then
                 meshCycles            = meshCycles   + 1
                 cycleList(meshCycles) = cez
                 DSC(i,cez)            = Size% ncornr + 1
               endif
             enddo
           enddo

           exit BoundaryLoop

         endif
       enddo

     enddo BoundaryLoop

   endif

!  Error Check

   if (newCorners == 0) then 
     call f90fatal("No starting corners found in getNewCorners!")
   endif



   return
   end subroutine getNewCorners 

