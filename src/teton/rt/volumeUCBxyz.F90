!***********************************************************************
!                        Last Update:  05/2017, PFN                    *
!                                                                      *
!   volumeUCBxyz   - Calculates corner and zone volumes for the        *
!                    upstream corner-balance (UCB) transport method    *
!                    in 3D cartesian geometry.                         *
!                                                                      *
!***********************************************************************
   subroutine volumeUCBxyz

   use kind_mod
   use Size_mod
   use Geometry_mod
   use constant_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Local Variables

   integer    :: cf
   integer    :: cface1
   integer    :: cface2
   integer    :: cfp
   integer    :: cez1
   integer    :: cez2
   integer    :: face
   integer    :: zone
   integer    :: zoneOpp
   integer    :: nCFaces

   integer    :: nzones 
   integer    :: c 
   integer    :: cc
   integer    :: c0 
   integer    :: cOpp 
   integer    :: nCorner 
   integer    :: nFaces
   integer    :: b
   integer    :: b0
   integer    :: n
   integer    :: nBoundary
   integer    :: nBdyElem

   real(adqt) :: tdl(3)
   real(adqt) :: tfl(3)
   real(adqt) :: tzl(3)
   real(adqt) :: A_fep(3)
   real(adqt) :: pnt0(3)
   real(adqt) :: pnt1(3)
   real(adqt) :: pnt2(3)
   real(adqt) :: zoneCenter(3)
   real(adqt) :: faceCenter(3,Size% maxFaces)

!  Dynamic arrays

   real(adqt), allocatable :: A_bdy(:,:)
   real(adqt), allocatable :: A_fp(:,:,:)

   allocate( A_bdy(3,Size% nbelem) )
   allocate( A_fp(3,Size% maxcf,Size% ncornr) )

!  Mesh Constants

   nzones    = Size%nzones
   nBoundary = getNumberOfBoundaries(RadBoundary)

!  Compute Vectors from edge-center to point (TEL), face-center to
!  point (TFL) and zone-center to point (TZL).  These are used
!  to compute outward normals on corner faces.  The corner-face
!  area is the sum of two half-side areas. 

   ZoneLoop: do zone=1,nzones

     nFaces = Geom% zoneFaces(zone)

!  Calculate the location of the face-centers and zone-centers

     zoneCenter(:)          = getZoneCenter(Geom, zone)
     faceCenter(:,1:nFaces) = getFaceCenter(Geom, zone, nFaces)

     nCorner                =  Geom% numCorner(zone)
     c0                     =  Geom% cOffSet(zone)
     Geom% VolumeZone(zone) = zero

     do c=1,nCorner
       Geom% Volume(c0+c)   = zero
       Geom% A_ez(:,:,c0+c) = zero
     enddo

     CornerLoop: do c=1,nCorner 

       cc      = c0 + c
       nCFaces = Geom% nCFacesArray(cc)

       CornerFaceLoop: do cface1=1,nCFaces

         cface2  = mod(cface1,nCFaces) + 1

         cfp     = Geom% cFP(cface1,cc)
         cez1    = Geom% cEZ(cface1,cc)
         cez2    = Geom% cEZ(cface2,cc)

         face    = Geom% CToFace(cface1,cc)

         pnt0(:) = Geom% px(:,c0+c)
         pnt1(:) = Geom% px(:,c0+cez1)
         pnt2(:) = Geom% px(:,c0+cez2)


         tdl(:)  = half*(pnt1(:)            - pnt2(:))
         tfl(:)  =       faceCenter(:,face) - pnt0(:)
         tzl(:)  =       zoneCenter(:)      - pnt0(:)

!  Calculate the components of the outward normals on
!  "FP" corner faces; this is the sum of two tet area vectors

         A_fep(1) = half*( tfl(3)*tdl(2) - tfl(2)*tdl(3) )
         A_fep(2) = half*( tfl(1)*tdl(3) - tfl(3)*tdl(1) )
         A_fep(3) = half*( tfl(2)*tdl(1) - tfl(1)*tdl(2) )

         zoneOpp  = Geom% zoneOpp(face,zone)

         if (zoneOpp > 0) then

!  Ensure that outward normals on FP faces are equal and opposite

           if ( zone < zoneOpp ) then
             A_fp(:,cface1,c0+c) = A_fep(:)

             do cf=1,Geom% nCFacesArray(cfp)
               if (Geom% cFP(cf,cfp) == cc) then
                 A_fp(:,cf,cfp) = -A_fep(:)
               endif
             enddo

           else
             A_fep(:) = A_fp(:,cface1,c0+c)
           endif

!  Set the outward normal on problem boundary

         elseif (zoneOpp < 0) then
           cOpp          = cfp - Size% ncornr
           A_bdy(:,cOpp) = A_fep(:)
         endif

!  Accumulate corner volumes

         Geom% Volume(cc) = Geom% Volume(cc) + third*abs( tzl(1)*A_fep(1) +  &
                                                          tzl(2)*A_fep(2) +  &
                                                          tzl(3)*A_fep(3) )

       enddo CornerFaceLoop

!  Accumulate zone volumes 

       Geom%VolumeZone(zone) = Geom% VolumeZone(zone) + Geom% Volume(cc) 

     enddo CornerLoop

   enddo ZoneLoop

!  Load components of the area vectors for the half-sides on
!  the problem boundary (only FP corner-faces live on the boundary)

   do n=1,nBoundary
     Bdy      => getBoundary(RadBoundary, n)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     b0       =  getFirstBdyElement(Bdy) - 1
     do b=1,nBdyElem
       Bdy% A_bdy(:,b) = A_bdy(:,b+b0)
     enddo
   enddo

!  Release temporary arrays

   deallocate( A_bdy )
   deallocate( A_fp  )



   return
   end subroutine volumeUCBxyz

