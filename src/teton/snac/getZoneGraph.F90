!***********************************************************************
!                        Last Update:  04/2024, PFN                    *
!                                                                      *
!   getZoneGraph - This routine builds the sweep ordering array        *
!                  (a directed graph) of ZONES for a single direction. *
!                                                                      *
!***********************************************************************
   subroutine getZoneGraph(aSetID, angle, nHyperDomains) 

   use kind_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: aSetID
   integer,    intent(in)    :: angle 
   integer,    intent(in)    :: nHyperDomains

!  Local Variables

   type(AngleSet), pointer   :: ASet

   integer    :: c
   integer    :: c0
   integer    :: Zexit
   integer    :: newZones
   integer    :: addedZones
   integer    :: lastZone
   integer    :: nextZone
   integer    :: meshCycles
   integer    :: ncornr
   integer    :: nzones

   integer    :: face 
   integer    :: nFaces
   integer    :: zone
   integer    :: zID
   integer    :: ndoneZ
   integer    :: nHyperPlanes

   real(adqt) :: omega(Size%ndim)

!  Dynamic

   integer,          allocatable :: needZ(:)
   integer,          allocatable :: listZone(:)
   integer,          allocatable :: zonesInPlane(:)
   integer,          allocatable :: cycleList(:)
   integer,          allocatable :: CToHypPlane(:)

   logical (kind=1), allocatable :: badZone(:)
   logical (kind=1), allocatable :: onCycleList(:)
   logical (kind=1), allocatable :: doneZ(:)
   logical (kind=1), allocatable :: exitFace(:,:)

!  Constants

   ASet     => getAngleSetData(Quad, aSetID) 

   ncornr   =  Size% ncornr
   nzones   =  Size% nzones
   omega(:) =  ASet% Omega(:,angle)

!  Allocate arrays

   allocate( needZ(nzones) )
   allocate( listZone(nzones) )
   allocate( zonesInPlane(nzones) )
   allocate( cycleList(ncornr) )
   allocate( CToHypPlane(ncornr) )
   allocate( badZone(nzones) )
   allocate( onCycleList(nzones) )
   allocate( doneZ(nzones) )
   allocate( exitFace(Size%maxFaces,nzones) )

   doneZ(:)        = .FALSE.
   onCycleList(:)  = .FALSE.

   meshCycles      = 0

!  Build NEED array by computing Outward_Normal dot Omega(m)

   call getZoneDependency(meshCycles, omega, NEEDZ, cycleList, exitFace, onCycleList)

!  Create a list of zones to start the sweep

   call getNewZones(newZones, meshCycles, needZ, listZone,  &
                    cycleList, exitFace, onCycleList) 

!  Check for zones the have circular dependencies

   call getDownStreamData(omega, aSetID, angle, meshCycles,  &
                          cycleList, badZone)

!  Create the "next" array. 

   ndoneZ       = 0
   nextZone     = 0
   lastZone     = 0
   nHyperPlanes = 0


   OuterIteration: do

!  Advance to a new hyper-plane

     nHyperPlanes               = nHyperPlanes + 1
     zonesInPlane(nHyperPlanes) = newZones
     nextZone                   = lastZone + newZones 
     addedZones                 = 0

!    Loop over all zones in the current list

     ZoneLoop: do zID=1,newZones

       zone    =  listZone(lastZone+zID)
       nFaces  =  Geom% zoneFaces(zone)

       ndoneZ      =  ndoneZ + 1
       doneZ(zone) = .TRUE.

!  Loop over the down-stream zones for the zone just added
!  to the nextZ list, decrementing the needZ array for these
!  neighboring zones 

       FaceLoop: do face=1,nFaces

         if ( exitFace(face,zone) ) then

           Zexit = Geom% zoneOpp(face,zone)

           if ( Zexit > 0) then
             if ( .not. doneZ(Zexit) ) then

               needZ(Zexit) = needZ(Zexit) - 1

               if (needZ(Zexit) == 0) then
                 nextZone           = nextZone   + 1
                 addedZones         = addedZones + 1
                 listZone(nextZone) = Zexit
               elseif (needZ(Zexit) < 0) then
                 write(6,100) Zexit, zone
                 call f90fatal("needZ < 0 in SNNEXT!")
               endif

             endif
           endif

         endif

       enddo FaceLoop

       if ( badZone(zone) ) then
         ASet% nextZ(ndoneZ,angle) = -zone
       else
         ASet% nextZ(ndoneZ,angle) =  zone
       endif

!    For the zone just added, create a map from corner to hyperplane

     c0 = Geom% cOffSet(zone)

     do c=1,Geom% numCorner(zone)
       CToHypPlane(c0+c) = nHyperPlanes
     enddo

     enddo ZoneLoop

     lastZone = lastZone + newZones

     if (lastZone == nzones) then

       exit OuterIteration

     else

       if (addedZones > 0) then

         newZones = addedZones

       elseif (addedZones == 0) then

!        Break a cycle to add a zone to the list

         call cycleBreakerZone(ndoneZ, meshCycles, nextZone, addedZones,  &
                               needZ, listZone, cycleList, exitFace,      &
                               onCycleList) 

         newZones = addedZones 

       endif

       cycle OuterIteration

     endif

   enddo OuterIteration

!  End of Outer Loop, save the number of hyperplanes

   if (meshCycles > ncornr) then
     call f90fatal("MeshCycles exceeds the number of corners in SNNEXT!") 
   endif

   ASet% numCycles(angle) = meshCycles

   call constructHyperPlane( ASet, angle, nHyperPlanes, meshCycles,   &
                             nHyperDomains, zonesInPlane(1:nHyperPlanes),  &
                             CToHypPlane, cycleList(1:meshCycles) )

!  Set the number of hyperplanes in the set module for this angle

   ASet% nHyperPlanes(angle) = nHyperPlanes

!  Final error check

   if (ndoneZ /= nzones) then
     call f90fatal("Wrong number of zones in SNNEXT!")
   endif

 100   format('zone ',i9,' has already been done and is down stream of ',i9)

!  Release memory

   deallocate( needZ )
   deallocate( listZone )
   deallocate( zonesInPlane )
   deallocate( cycleList )
   deallocate( CToHypPlane )
   deallocate( badZone )
   deallocate( onCycleList )
   deallocate( doneZ )
   deallocate( exitFace ) 

 
   return
   end subroutine getZoneGraph 

