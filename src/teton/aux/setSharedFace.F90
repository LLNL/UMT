!***********************************************************************
!                        Version 1:  11/2010, PFN                      *
!                                                                      *
!   setSharedFace   - This routine sets geometry information for       *
!                     faces on shared boundaries.                      *
!                                                                      *
!***********************************************************************

   subroutine setSharedFace(bcID, zoneID, faceID, cornerID) &
        BIND(C,NAME="teton_setsharedface")

!  Include
   USE ISO_C_BINDING
   use kind_mod
   use Size_mod
   use Geometry_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Arguments

   integer(C_INT),    intent(in)    :: bcID
   integer(C_INT),    intent(in)    :: zoneID
   integer(C_INT),    intent(in)    :: faceID
   integer(C_INT),    intent(in)    :: cornerID

!  Local 

   integer :: ii
   integer :: c0, cID, nCFaces
   integer :: b0, bdyelem, nbelem

   nbelem = Size% nbelem

!  First set boundary element numbers

   Bdy  => getBoundary(RadBoundary, bcID)

   b0   =  getFirstBdyElement(Bdy) - 1
   c0   =  Geom% cOffSet(zoneID+1)

!  Add a boundary element to the list for this corner-face 

   Bdy% BdyElemCtr         = Bdy% BdyElemCtr + 1
   bdyelem                 = Bdy% BdyElemCtr 

   cID                     = cornerID + 1 - c0
   Bdy% BdyToZone(bdyelem) = zoneID   + 1
   Bdy% BdyToC(bdyelem)    = cornerID + 1 

   if (Size% ndim == 2) then

     do ii=1,2
       if (Geom% cFP(ii,c0+cID) == -(nbelem+faceID+1)) then
         Geom% cFP(ii,c0+cID) = Size% ncornr + (b0 + bdyelem)
       endif
     enddo

   elseif (Size% ndim == 3) then

     nCFaces = Geom% nCFacesArray(c0+cID)
     do ii=1,nCFaces
       if (Geom% cFP(ii,c0+cID) == -(nbelem+faceID+1)) then
         Geom% cFP(ii,c0+cID) = Size% ncornr + (b0 + bdyelem)
       endif
     enddo

   elseif (Size% ndim == 1) then

     Geom% cFP(1,c0+cID) = Size% ncornr + (b0 + bdyelem)

   endif


   return
   end subroutine setSharedFace 


