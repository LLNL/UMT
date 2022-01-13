!***********************************************************************
!                        Last Update:  11/2017, PFN                    *
!                                                                      *
!   getDownStreamData - Creates a list of downstream corners for the   *
!                       input zone and angle.                          *
!                                                                      *
!***********************************************************************

   subroutine getDownStreamData(omega, aSetID, angle, meshCycles,  &
                                cycleList, badZone)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use AngleSet_mod

   implicit none

!  Arguments

   real(adqt),       intent(in)    :: omega(Size%ndim)

   integer,          intent(in)    :: aSetID
   integer,          intent(in)    :: angle
   integer,          intent(inout) :: meshCycles
   integer,          intent(inout) :: cycleList(Size%ncornr)
   logical (kind=1), intent(inout) :: badZone(Size%nzones)

!  Local Variables

   type(AngleSet), pointer   :: ASet

   integer, dimension (1)    :: cfirst

   integer    :: zone
   integer    :: i
   integer    :: c 
   integer    :: c0 
   integer    :: cc 
   integer    :: cez 
   integer    :: cface 
   integer    :: nCorner 
   integer    :: nCFaces
   integer    :: Cexit
   integer    :: minNeed
   integer    :: ndim

   integer    :: need(Size%maxCorner)
   integer    :: DownStreamC(Size%maxcf,Size%maxCorner)
   integer    :: nDSC(Size%maxCorner)

   real(adqt) :: aez 

!  Constants

   ASet => getAngleSetData(Quad, aSetID)

   ndim = Size% ndim

!  Find the corners in this zone with need=0; if a zone cannot
!  be started or there is a circular dependency with the zone, we
!  add the zone to the cycleList (done in "fixZone")

!  For incoming corner-faces we increment the need array; for outgoing
!  corner-faces we put the downstream corner number into an index list.


   ZoneLoop: do zone=1,Size% nzones

     nDSC(:) = 0 
     need(:) = 0 

     nCorner = Geom% numCorner(zone)
     c0      = Geom% cOffSet(zone)

     CornerLoop: do c=1,nCorner

       cc = c0 + c

       if (ndim == 2) then
         nCFaces = 2
       else
         nCFaces = Geom% nCFacesArray(cc)
       endif

       CornerFaceLoop: do cface=1,nCFaces
 
!        Get downstream corner number

         cez  = Geom% cEZ(cface,cc)

!  Omega dot Outward normal - IMPORTANT: the dot product must be
!  coded this way to be compatible with the coding in SNSWP3D and SNSWP2D.
!  Failure to comply results in wrong answers!

         if (cez > c) then

           aez = DOT_PRODUCT( omega(:),Geom% A_ez(:,cface,cc) )

           if (aez < zero) then
             need(c)                    = need(c)   + 1
             nDSC(cez)                  = nDSC(cez) + 1
             DownStreamC(nDSC(cez),cez) = c
           elseif (aez > zero) then
             need(cez)                  = need(cez) + 1
             nDSC(c)                    = nDSC(c)   + 1
             DownStreamC(nDSC(c),c)     = cez
           endif

         endif

       enddo CornerFaceLoop

     enddo CornerLoop

!    Find the corners in this zone with need=0; if a zone cannot
!    be started or there is a circular dependency with the zone, we
!    add the zone to the cycleList (done in "fixZone")

     badZone(zone) = .FALSE.

     do i=1,nCorner
       cfirst                  = minloc( need(1:nCorner) )
       c                       = cfirst(1)
       minNeed                 = need(c)
       ASet% nextC(c0+i,angle) = c

       if ( minNeed /= 0 ) then
         badZone(zone) = .TRUE.
       endif

       do cface=1,nDSC(c)
         Cexit       = DownStreamC(cface,c)
         need(Cexit) = need(Cexit) - 1
       enddo
       need(c) = 99
     enddo

     if ( badZone(zone) ) then

       do i=1,nCorner
         ASet% nextC(c0+i,angle) = i
       enddo

       call fixZone(zone, meshCycles, cycleList)

     endif

   enddo ZoneLoop

 
   return
   end subroutine getDownStreamData 

