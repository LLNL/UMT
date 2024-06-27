!***********************************************************************
!                        Last Update:  04/2019, PFN                    *
!                                                                      *
!   getCornerDependency - This routine builds the NEED array which     *
!                         indicates the number of incoming fluxes      *
!                         required to compute the outgoing flux for    *
!                         a particular direction.                      *
!                                                                      *
!***********************************************************************

   subroutine getCornerDependency(omega, need, nDSC, DSC)

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod

   implicit none

!  Arguments

   integer,          intent(inout) :: need(Size%ncornr)
   integer,          intent(inout) :: nDSC(Size%ncornr)
   integer,          intent(inout) :: DSC(2*Size%maxcf,Size%ncornr)

   real(adqt),       intent(in)    :: omega(Size%ndim)

!  Local Variables

   integer         :: c
   integer         :: c0
   integer         :: cez
   integer         :: cfp
   integer         :: cface 
   integer         :: nCFaces
   integer         :: nCorner
   integer         :: zone

   real(adqt)      :: aez
   real(adqt)      :: afp 

!  For incoming corner-faces we increment the need array; for outgoing
!  corner-faces we put the downstream corner number into an index list.

   need(:) = 0 
   nDSC(:) = 0
   nCorner = Size% ncornr

   CornerLoopEZ: do c=1,nCorner 

     zone = Geom% CToZone(c)
     c0   = Geom% cOffSet(zone)

     if (Size% ndim == 2) then
       nCFaces = 2
     else
       nCFaces = Geom% nCFacesArray(c)
     endif

     CornerFaceLoopEZ: do cface=1,nCFaces

       cez = c0 + Geom% cEZ(cface,c)

       if (c < cez) then

         aez = DOT_PRODUCT( omega(:),Geom% A_ez(:,cface,c) )

         if ( aez < zero ) then
           need(c)            = need(c)   + 1
           nDSC(cez)          = nDSC(cez) + 1
           DSC(nDSC(cez),cez) = c
         elseif ( aez > zero ) then
           need(cez)          = need(cez) + 1
           nDSC(c)            = nDSC(c)   + 1
           DSC(nDSC(c),c)     = cez
         endif

       endif

     enddo CornerFaceLoopEZ

   enddo CornerLoopEZ

!  Zone Faces

   if (Size% ndim == 2) then

     CornerLoop2D: do c=1,nCorner

       CornerFaceLoop2D: do cface=1,2

         cfp = Geom% cFP(cface,c)

!  If an fp-face is on a boundary (cfp > nCorner) then do not
!  increment "need"

       if (c < cfp .and. cfp <= nCorner) then
         afp = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c) )

         if (afp < zero) then
           need(c)            = need(c)   + 1
           nDSC(cfp)          = nDSC(cfp) + 1
           DSC(nDSC(cfp),cfp) = c

         elseif (afp > zero) then
           need(cfp)          = need(cfp) + 1
           nDSC(c)            = nDSC(c)   + 1
           DSC(nDSC(c),c)     = cfp
         endif
 
       endif

       enddo CornerFaceLoop2D

     enddo CornerLoop2D

   elseif (Size% ndim == 3) then

     CornerLoop3D: do c=1,nCorner

       nCFaces = Geom% nCFacesArray(c)

       CornerFaceLoop3D: do cface=1,nCFaces

         cfp = Geom% cFP(cface,c)

!  If an fp-face is on a boundary (cfp > nCorner) then do not
!  increment "need"

       if (c < cfp .and. cfp <= nCorner) then
         afp = DOT_PRODUCT( omega(:),Geom% A_fp(:,cface,c) )

         if (afp < zero) then
           need(c)            = need(c)   + 1
           nDSC(cfp)          = nDSC(cfp) + 1
           DSC(nDSC(cfp),cfp) = c

         elseif (afp > zero) then
           need(cfp)          = need(cfp) + 1
           nDSC(c)            = nDSC(c)   + 1
           DSC(nDSC(c),c)     = cfp
         endif

       endif

       enddo CornerFaceLoop3D

     enddo CornerLoop3D

   endif


   return
   end subroutine getCornerDependency 

