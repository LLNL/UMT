!***********************************************************************
!                        Last Update:  05/2017, PFN                    *
!                                                                      *
!   geometryUCBxyz - Calculates certain geometry factors for the       *
!                    upstream corner-balance (UCB) transport method    *
!                    in 3D cartesian geometry.                         *
!                                                                      *
!***********************************************************************
   subroutine geometryUCBxyz

   use kind_mod
   use Size_mod
   use Geometry_mod
   use constant_mod

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
   integer    :: cez 
   integer    :: nCorner 
   integer    :: nFaces 

   real(adqt) :: tdl(3)
   real(adqt) :: tfl(3)
   real(adqt) :: tzl(3)
   real(adqt) :: tfz(3)
   real(adqt) :: tfe1(3)
   real(adqt) :: tfe2(3)
   real(adqt) :: A_fep(3)
   real(adqt) :: A_fez1(3)
   real(adqt) :: A_fez2(3)
   real(adqt) :: pnt0(3)
   real(adqt) :: pnt1(3)
   real(adqt) :: pnt2(3)
   real(adqt) :: zoneCenter(3)
   real(adqt) :: faceCenter(3,Size% maxFaces)

!  Mesh Constants

   nzones = Size%nzones

!  Compute Vectors from edge-center to point (TEL), face-center to
!  point (TFL) and zone-center to point (TZL).  These are used
!  to compute outward normals on corner faces.  The corner-face
!  area is the sum of two half-side areas. 

   ZoneLoop: do zone=1,nzones

     nFaces = Geom% zoneFaces(zone) 

!  Calculate the location of the face-centers and zone-centers

     zoneCenter(:)          = getZoneCenter(Geom, zone)
     faceCenter(:,1:nFaces) = getFaceCenter(Geom, zone, nFaces)

     nCorner =  Geom% numCorner(zone)
     c0      =  Geom% cOffSet(zone)

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
         tfz(:)  =       faceCenter(:,face) - zoneCenter(:)
         tfe1(:) =       faceCenter(:,face) - half*(pnt1(:) + pnt0(:))
         tfe2(:) =       faceCenter(:,face) - half*(pnt2(:) + pnt0(:))

!  Calculate the components of the outward normals on
!  "FP" corner faces; this is the sum of two tet area vectors

         A_fep(1) = half*( tfl(3)*tdl(2) - tfl(2)*tdl(3) )
         A_fep(2) = half*( tfl(1)*tdl(3) - tfl(3)*tdl(1) )
         A_fep(3) = half*( tfl(2)*tdl(1) - tfl(1)*tdl(2) )

         zoneOpp  = Geom% zoneOpp(face,zone)

         if (zoneOpp > 0) then

!  Ensure that outward normals on FP faces are equal and opposite

           if ( zoneOpp > zone) then
             Geom% A_fp(:,cface1,cc) = A_fep(:)

             do cf=1,Geom% nCFacesArray(cfp)
               if (Geom% cFP(cf,cfp) == cc) then
                 Geom% A_fp(:,cf,cfp) = -A_fep(:)
               endif
             enddo
           else
             A_fep(:) = Geom% A_fp(:,cface1,cc)
           endif

!  Set the outward normal on problem boundary

         elseif (zoneOpp < 0) then
           Geom% A_fp(:,cface1,cc) = A_fep(:)
         endif

!  "EZ" corner faces; here we add the tet area vectors
!  to two different "EZ" faces

         A_fez1(1) = half*( tfz(3)*tfe1(2) - tfz(2)*tfe1(3) )
         A_fez1(2) = half*( tfz(1)*tfe1(3) - tfz(3)*tfe1(1) )
         A_fez1(3) = half*( tfz(2)*tfe1(1) - tfz(1)*tfe1(2) )

         A_fez2(1) = half*( tfz(2)*tfe2(3) - tfz(3)*tfe2(2) )
         A_fez2(2) = half*( tfz(3)*tfe2(1) - tfz(1)*tfe2(3) )
         A_fez2(3) = half*( tfz(1)*tfe2(2) - tfz(2)*tfe2(1) )

!  Accumulate corner surface areas on "EZ" faces

         Geom% A_ez(:,cface1,cc) = Geom% A_ez(:,cface1,cc) + A_fez1(:)
         Geom% A_ez(:,cface2,cc) = Geom% A_ez(:,cface2,cc) + A_fez2(:)

!  Accumulate corner volumes

         Geom% Volume(cc) = Geom% Volume(cc) + third*abs( tzl(1)*A_fep(1) +  &
                                                          tzl(2)*A_fep(2) +  &
                                                          tzl(3)*A_fep(3) )


       enddo CornerFaceLoop

     enddo CornerLoop

!  Ensure that outward normals on "EZ" faces are equal and opposite

     do c=1,nCorner
       cc = c0 + c
       do cface1=1,Geom% nCFacesArray(cc)
         cez = Geom% cEZ(cface1,cc)
         if (cez > c) then
           do cface2=1,Geom% nCFacesArray(c0+cez)
             if (Geom% cEZ(cface2,c0+cez) == c) then
               Geom% A_ez(:,cface2,c0+cez) = -Geom% A_ez(:,cface1,cc)
             endif
           enddo
         endif
       enddo
     enddo

   enddo ZoneLoop


   return
   end subroutine geometryUCBxyz

