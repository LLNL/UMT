!***********************************************************************
!                        Version 1:  02/06, PFN                        *
!                                                                      *
!   setTetonZone - This routine sets geometry information in a zone    *
!                  data structure.                                     *
!                                                                      *
!***********************************************************************

   subroutine setTetonZone(zoneID, corner0, zoneFaces, cornerFaces,    &
                           zoneNCorner, zoneOpp, CornerID, CornerOpp,  &
                           nCPerFace, FaceToBCList) &
                           BIND(C,NAME="teton_setzone")

!  Include
   USE ISO_C_BINDING
   use flags_mod
   use kind_mod
   use Size_mod
   use Geometry_mod
   use ZoneData_mod
   use MeshData_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Arguments

   integer(C_INT),    intent(in)    :: zoneID
   integer(C_INT),    intent(in)    :: corner0
   integer(C_INT),    intent(in)    :: zoneFaces
   integer(C_INT),    intent(in)    :: cornerFaces
   integer(C_INT),    intent(in)    :: zoneNCorner 
   integer(C_INT),    intent(in)    :: zoneOpp(zoneFaces)
   integer(C_INT),    intent(in)    :: CornerID(cornerFaces)
   integer(C_INT),    intent(inout) :: CornerOpp(cornerFaces)
   integer(C_INT),    intent(in)    :: nCPerFace(zoneFaces)
   integer(C_INT),    intent(in)    :: FaceToBCList(zoneFaces)

!  Local 

   integer :: nCorner 
   integer :: numC
   integer :: nSides
   integer :: side0
   integer :: face
   integer :: i 
   integer :: ii 
   integer :: faceIndex 
   integer :: iCW, iCCW
   integer :: c, c1, c2, cCWLast
   integer :: bcID, b0, bdyelem

   integer, parameter :: fp=1
   integer, parameter :: ez=2
   integer, parameter :: cc=3

   integer :: cFaceID(Size%maxCorner)
   integer :: cCW(Size%maxcf,Size%maxCorner)
   integer :: cCCW(Size%maxcf,Size%maxCorner)
   integer :: cOppZone(Size%maxcf,Size%maxCorner)
   integer :: cFP(Size%maxcf,Size%maxCorner)
   integer :: cEZ(Size%maxcf,Size%maxCorner)
   integer :: faceID(Size%maxcf,Size%maxCorner)
   integer :: CToFace(Size%maxcf,Size%maxCorner)
   integer :: nCFaces(Size%maxCorner)


   nCorner    = zoneNCorner 
   nSides     = nCorner
   side0      = corner0
   nCFaces(:) = 0
   cFP(:,:)   = 0
   cEZ(:,:)   = 0

!  Find the maximum number of faces per zone - used in several routines
!  to dimension temporary arrays

   Size% maxFaces = max(Size% maxFaces, zoneFaces)
   Size% maxSides = max(Size% maxSides, nSides)

!  First set boundary element numbers

   faceIndex = 0

   do face=1,zoneFaces
     numC = nCPerFace(face)
     if (zoneOpp(face) < 0) then
       bcID =  FaceToBCList(face)
       Bdy  => getBoundary(RadBoundary, bcID)
       b0   =  getFirstBdyElement(Bdy) - 1

       if (Bdy% BCType == bcType_shared) then

!  These faces are set in "setSharedFace" 

       else

!  The rest are ordered by corner ID

         do i=1,numC
           Bdy% BdyElemCtr         = Bdy% BdyElemCtr + 1
           bdyelem                 = Bdy% BdyElemCtr
           Bdy% BdyToC(bdyelem)    = corner0 + CornerID(faceIndex+i)
           Bdy% BdyToZone(bdyelem) = zoneID
           CornerOpp(faceIndex+i)  = Size% ncornr + (b0 + bdyelem)
         enddo

       endif

     endif
     faceIndex = faceIndex + numC
   enddo

!  Make a list of neighboring corners in the zone 

   if (Size% ndim == 2) then

     do face=1,zoneFaces

       c1          = CornerID(2*face-1)
       c2          = CornerID(2*face)

       cFP(2,c1)   = CornerOpp(2*face-1)
       cFP(1,c2)   = CornerOpp(2*face)

       cEZ(1,c1)   = c2
       cEZ(2,c2)   = c1 

       nCFaces(c1) = nCFaces(c1) + 1
       nCFaces(c2) = nCFaces(c2) + 1

       CToFace(2,c1) = face
       CToFace(1,c2) = face

     enddo

   elseif (Size% ndim == 3) then 

     cFaceID(:) = 0 
     faceIndex  = 0

     do face=1,zoneFaces

       numC = nCPerFace(face) 

       do i=1,numC
         iCCW = numC - mod(numC-i+1,numC)
         iCW  = mod(i,numC) + 1
         c    = CornerID(faceIndex+i)

         cFaceID(c) = cFaceID(c) + 1

         cCW(cFaceID(c),c)      = CornerID(faceIndex+iCW) 
         cCCW(cFaceID(c),c)     = CornerID(faceIndex+iCCW) 
         cOppZone(cFaceID(c),c) = CornerOpp(faceIndex+i) 
         faceID(cFaceID(c),c)   = face
         nCFaces(c)             = nCFaces(c) + 1
       enddo

       faceIndex = faceIndex + numC

     enddo

!    Set the zone connectivity

     do c=1,nCorner

       cEZ(1,c)     = cCCW(1,c)
       cFP(1,c)     = cOppZone(1,c)
       cCWLast      = cCW(1,c)
       CToFace(1,c) = faceID(1,c)

       do i=2,nCFaces(c)
         CornerFaceLoop: do ii=2,nCFaces(c)
           if (cCCW(ii,c) == cCWLast) then
             cEZ(i,c)     = cCCW(ii,c)
             cFP(i,c)     = cOppZone(ii,c)
             cCWLast      = cCW(ii,c)
             CToFace(i,c) = faceID(ii,c) 
             exit CornerFaceLoop
           endif
         enddo CornerFaceLoop 
       enddo

     enddo

   endif

!  Set the zone data structures

   Z => getZoneData(Geom, zoneID)
   M => getMesh(Geom, zoneID)

   call constructZone( Z, nCorner, corner0, zoneFaces,   &
                       nSides, side0, cFP )

   call constructMesh( M, nCorner, corner0, zoneFaces,   &
                       zoneOpp, CToFace, nCFaces )

   call setConnectivity(Geom, zoneID, nCorner, corner0,  &
                        cFP, cEZ, nCFaces)


   return
   end subroutine setTetonZone 


