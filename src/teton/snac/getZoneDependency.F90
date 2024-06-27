!***********************************************************************
!                        Last Update:  05/2023, PFN                    *
!                                                                      *
!   getZoneDependency - This routine builds the NEED array which       *
!                       indicates the number of incoming fluxes        *
!                       required to compute the outgoing flux for      *
!                       a particular direction.                        *
!                                                                      *
!***********************************************************************

   subroutine getZoneDependency(MESHCYCLES, omega, NEEDZ, cycleList, &
                                exitFace, onCycleList) 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod

   implicit none

!  Arguments

   integer,    intent(inout)       :: MeshCycles

   integer,    intent(inout)       :: needZ(Size%nzones)
   integer,    intent(inout)       :: cycleList(Size%ncornr)

   logical (kind=1), intent(inout) :: exitFace(Size%maxFaces,Size%nzones)
   logical (kind=1), intent(inout) :: onCycleList(Size%nzones)

   real(adqt), intent(in)          :: omega(Size%ndim)

!  Local Variables

   integer    :: izero, nzones, nFaces, nCorner, nCFaces
   integer    :: c, c0, cc, cface, zone, zoneOpp, face
   integer    :: faceOpp

   real(adqt) :: afpm 

   integer          :: nInc(Size% maxCorner)
   integer          :: nExit(Size% maxCorner)
   real(adqt)       :: afpm_Face(Size% maxCorner)

!  Constants

   parameter (izero=0)

   nzones     = Size%nzones
   MeshCycles = 0

!  For incoming corner-faces we increment the need array; for outgoing
!  corner-faces we put the downstream corner number into an index list.

   needZ(:)        = izero
   exitFace(:,:)   = .FALSE.

   if (Size% ndim == 2) then

     ZoneLoop2D: do zone=1,nzones

       nCorner =  Geom% numCorner(zone)
       c0      =  Geom% cOffSet(zone)

       CornerLoop2D: do c=1,nCorner 

         face    = Geom% CToFace(1,c0+c)
         zoneOpp = Geom% zoneOpp(face,zone)

         if (zone < zoneOpp) then

           faceOpp = Geom% faceOpp(face,zone)
           afpm    = DOT_PRODUCT( omega(:),Geom% A_fp(:,1,c0+c) )

           if (afpm < zero) then
             needZ(zone)               = needZ(zone)   + 1
             exitFace(faceOpp,zoneOpp) = .TRUE.
           elseif (afpm > zero) then
             needZ(zoneOpp)            = needZ(zoneOpp) + 1
             exitFace(face,zone)       = .TRUE.
           endif
         endif

       enddo CornerLoop2D

     enddo ZoneLoop2D

   elseif (Size% ndim == 3) then

     ZoneLoop: do zone=1,nzones

       nCorner =  Geom% numCorner(zone)
       c0      =  Geom% cOffSet(zone)
       nFaces  =  Geom% zoneFaces(zone)

       afpm_Face(:) = zero
       nInc(:)      = izero
       nExit(:)     = izero

       CornerLoop: do c=1,nCorner

         cc      = c0 + c

         nCFaces = Geom% nCFacesArray(cc)

         CornerFaceLoop: do cface=1,nCFaces
 
!          Get downstream zone number

           face    = Geom% CToFace(cface,cc)
           zoneOpp = Geom% zoneOpp(face,zone)

!  Omega dot Outward normal - IMPORTANT: the dot product must be
!  coded this way to be compatible with the coding in SNSWP3D and SNSWP2D.
!  Failure to comply results in wrong answers!

!  Corner Face FP (neighboring corner, neighboring zone)
       
           if (zoneOpp > zone) then

             afpm = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,cc) )

             afpm_Face(face) = afpm_Face(face) + afpm

             if (afpm < zero) then
               nInc(face)  = nInc(face)  + 1
             elseif (afpm > zero) then
               nExit(face) = nExit(face) + 1 
             endif

           endif

         enddo CornerFaceLoop

       enddo CornerLoop

       FaceLoop: do face=1,nFaces

         zoneOpp = Geom% zoneOpp(face,zone)

         if (zoneOpp > zone) then

           faceOpp = Geom% faceOpp(face,zone)

           if (afpm_Face(face) < zero) then
             needZ(zone)               = needZ(zone) + 1
             exitFace(faceOpp,zoneOpp) = .TRUE.

             if ( nExit(face) > 0 ) then
               if ( .not. onCycleList(zone) ) then
!  Only add this zone to the cycle list once
                 do c=1,nCorner
                   meshCycles            = meshCycles + 1
                   cycleList(meshCycles) = c0 + c
                 enddo

                 onCycleList(zone) = .TRUE.

               endif
             endif

           elseif (afpm_Face(face) > zero) then
             needZ(zoneOpp)            = needZ(zoneOpp) + 1
             exitFace(face,zone)       = .TRUE.

             if ( nInc(face) > 0 ) then
               if ( .not. onCycleList(zoneOpp) ) then
!  Add the neighboring zone across this face
                 do c=1,Geom% numCorner(zoneOpp)
                   meshCycles            = meshCycles + 1
                   cycleList(meshCycles) = Geom% cOffSet(zoneOpp) + c
                 enddo

                 onCycleList(zoneOpp) = .TRUE.
               endif
             endif

           endif

         endif

       enddo FaceLoop

     enddo ZoneLoop

   endif


   return
   end subroutine getZoneDependency 

