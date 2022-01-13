!***********************************************************************
!                        Version 1:  02/06, PFN                        *
!                                                                      *
!   updateGeometry - This routine updates phase-space set geometry.    *
!                                                                      *
!***********************************************************************

   subroutine updateGeometry(setID)

!  Include

   use kind_mod
   use Size_mod
   use Geometry_mod
   use ZoneData_mod
   use QuadratureList_mod
   use SetData_mod

   implicit none

!  Arguments

   integer,  intent(in)    :: setID

!  Local

   type(SetData),  pointer :: Set
   type(ZoneData), pointer :: ZT

   integer                 :: c
   integer                 :: c0
   integer                 :: nCorner
   integer                 :: zone

!  "Set" geometry


   Set => getSetData(Quad, setID)

   if (Size%ndim == 2) then

     if ( Size% usePWLD ) then
       call geometryPWLDrz(setID)
     else
       call geometryUCBrz(setID)
     endif

   else if (Size%ndim == 3) then

     do zone=1,Size%nzones
       ZT => getZoneData(Geom, zone)

       nCorner = ZT% nCorner
       c0      = ZT% c0

       do c=1,nCorner
         Set% Volume(c0+c)   = ZT% Volume(c)
         Set% A_fp(:,:,c0+c) = ZT% A_fp(:,:,c)
         Set% A_ez(:,:,c0+c) = ZT% A_ez(:,:,c)
       enddo
     enddo

   endif


   return
   end subroutine updateGeometry 


