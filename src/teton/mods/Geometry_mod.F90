! Geometry Modules:  geometry and connectivity information used by Sn
                                                                                 
module Geometry_mod 

  use kind_mod

  private

! public interfaces
                                                                                             
  public construct 
  public destruct 
  public getZoneAverage
  public setConnectivity
  public getZoneCenter
  public getFaceCenter
                                                                                 
                                                                                 
  type, public :: Geometry 

!    Connectivity 

     integer,         pointer, contiguous :: cEZ(:,:) => null()
     integer,         pointer, contiguous :: cFP(:,:) => null()
     integer,         pointer, contiguous :: nCFacesArray(:) => null()
     integer,         pointer, contiguous :: cOffSet(:) => null()
     integer,         pointer, contiguous :: numCorner(:) => null()
     integer,         pointer, contiguous :: CToZone(:) => null()
     integer,         pointer, contiguous :: zone1(:) => null()
     integer,         pointer, contiguous :: zone2(:) => null()
     integer,         pointer, contiguous :: corner1(:) => null()
     integer,         pointer, contiguous :: corner2(:) => null()
     integer,         pointer, contiguous :: zoneFaces(:) => null()
     integer,         pointer, contiguous :: zoneOpp(:,:) => null()
     integer,         pointer, contiguous :: faceOpp(:,:) => null()
     integer,         pointer, contiguous :: CToFace(:,:) => null()

!    Geometry Factors

     real(adqt),      pointer, contiguous :: Volume(:) => null()
     real(adqt),      pointer, contiguous :: VolumeOld(:) => null()
     real(adqt),      pointer, contiguous :: VolumeZone(:) => null()
     real(adqt),      pointer, contiguous :: Area(:) => null()
     real(adqt),      pointer, contiguous :: Radius(:) => null()
     real(adqt),      pointer, contiguous :: RadiusEZ(:,:) => null()
     real(adqt),      pointer, contiguous :: RadiusFP(:,:) => null()
     real(adqt),      pointer, contiguous :: A_fp(:,:,:) => null()
     real(adqt),      pointer, contiguous :: A_ez(:,:,:) => null()
     real(adqt),      pointer, contiguous :: VoC(:,:) => null()
     real(adqt),      pointer, contiguous :: px(:,:)  => null()

! PWLD

     real(adqt),      pointer, contiguous :: LL(:,:,:,:) => null()
     real(adqt),      pointer, contiguous :: AD(:,:,:) => null()
     real(adqt),      pointer, contiguous :: MM(:,:) => null()

! Unlumped PWLD

     real(adqt),      pointer, contiguous :: LLu(:,:,:,:) => null()
     real(adqt),      pointer, contiguous :: ADu(:,:,:) => null()
     real(adqt),      pointer, contiguous :: MMu(:,:,:) => null()

! 1D only

     real(adqt), pointer, contiguous :: zoneWidth(:) => null() ! delta radius
     real(adqt), pointer, contiguous :: Rave(:)      => null() ! average radius
     real(adqt), pointer, contiguous :: Rmin(:)      => null() ! inner radius
     real(adqt), pointer, contiguous :: Rmax(:)      => null() ! outer radius

! Misc arrays

     logical(kind=1), pointer, contiguous :: BoundaryZone(:) ! is zone on a boundary

!    Misc
     character(len=8) :: label ! A string descriptor for this set.

  end type Geometry

  type(Geometry), pointer, public :: Geom => null()

  interface construct
    module procedure Geometry_ctor
  end interface

  interface getZoneAverage
    module procedure Geometry_getZoneAverage
  end interface

  interface setConnectivity 
    module procedure Geometry_setConnectivity
  end interface

  interface getZoneCenter
    module procedure Geometry_getZoneCenter
  end interface

  interface getFaceCenter
    module procedure Geometry_getFaceCenter
  end interface

  interface destruct
    module procedure Geometry_dtor
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine Geometry_ctor(self, nZoneSets)

    use Size_mod
    use constant_mod
    use MemoryAllocator_mod

    implicit none

!   Passed variables

    type(Geometry),  intent(inout) :: self

    integer,         intent(in)    :: nZoneSets

!   Local

    integer                        :: zonesPerSet
    integer                        :: zonesTotal
    integer                        :: nSetsP
    integer                        :: setID

    self%label = "geometry"

!   Pointers

    allocate( self% VoC(Size% ndim,Size% ncornr) )
    allocate( self% px(Size% ndim,Size% ncornr) )
    allocate( self% zoneFaces(Size% nzones) )
    allocate( self% zoneOpp(Size%maxFaces,Size% nzones) )
    allocate( self% faceOpp(Size%maxFaces,Size% nzones) )
    allocate( self% CToFace(Size% maxcf,Size% ncornr) )
    allocate( self% BoundaryZone(Size% nzones) )

    call Allocator%allocate(Size%usePinnedMemory,self%label,"Volume",       self% Volume,    Size% ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"VolumeOld",    self% VolumeOld, Size% ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"VolumeZone",   self% VolumeZone,Size% nzones)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"cOffSet",      self% cOffSet,   Size% nzones)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"numCorner",    self% numCorner, Size% nzones)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"CToZone",      self% CToZone,   Size% ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"zone1",        self% zone1,     nZoneSets)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"zone2",        self% zone2,     nZoneSets)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"corner1",      self% corner1,   nZoneSets)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"corner2",      self% corner2,   nZoneSets)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"cFP",          self% cFP,       Size% maxcf,Size% ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"cEZ",          self% cEZ,       Size% maxcf,Size% ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"A_fp",         self% A_fp,      Size% ndim,Size% maxcf,Size% ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"A_ez",         self% A_ez,      Size% ndim,Size% maxcf,Size% ncornr)
    call Allocator%allocate(Size%usePinnedMemory,self%label,"nCFacesArray", self% nCFacesArray, Size% ncornr)

    self% VolumeZone(:)   = zero
    self% VoC(:,:)        = zero
    self% px(:,:)         = zero
    self% zoneFaces(:)    = 0
    self% zoneOpp(:,:)    = 0
    self% faceOpp(:,:)    = 0
    self% CToFace(:,:)    = 0
    self% BoundaryZone(:) = .FALSE.

    if (Size% ndim == 1) then
      allocate( self% Area(Size% ncornr) )
      allocate( self% Radius(Size% ncornr) )
      allocate( self% zoneWidth(Size% nzones) )
      allocate( self% Rave(Size% nzones) )
      allocate( self% Rmin(Size% nzones) )
      allocate( self% Rmax(Size% nzones) )

      self% Area(:)      = zero
      self% Radius(:)    = zero
      self% zoneWidth(:) = zero
      self% Rave(:)      = zero
      self% Rmin(:)      = zero
      self% Rmax(:)      = zero

    elseif (Size% ndim == 2) then

      call Allocator%allocate(Size%usePinnedMemory,self%label,"Area",     self% Area,      Size% ncornr)
      call Allocator%allocate(Size%usePinnedMemory,self%label,"RadiusEZ", self% RadiusEZ,  2,Size% ncornr)
      call Allocator%allocate(Size%usePinnedMemory,self%label,"RadiusFP", self% RadiusFP,  2,Size% ncornr)

    endif

!   PWLD

    if ( Size% usePWLD ) then
      if ( Size% useSurfaceMassLumping ) then
        allocate( self% LL(2,Size% maxSides,Size% maxSides,Size% nZones) )
        allocate( self% AD(Size% maxSides,Size% maxSides,Size% nZones) )
        allocate( self% MM(Size% maxSides,Size% nZones) )
      else
        allocate( self% LLu(2,Size% maxSides,Size% maxSides,Size% nZones) )
        allocate( self% ADu(Size% maxSides,Size% maxSides,Size% nZones) )
        allocate( self% MMu(Size% maxSides,Size% maxSides,Size% nZones) )
      endif
    endif

!   Zone Sets: Assign zone range to each zone set; sets 1->nSetsP get an extra
!   zone if there is a remainder

    zonesPerSet = Size% nzones/nZoneSets
    nSetsP      = Size% nzones - (nZoneSets*zonesPerSet)
    zonesTotal  = 0

    do setID=1,nSetsP
      self% zone1(setID) = zonesTotal + 1
      self% zone2(setID) = zonesTotal + zonesPerSet + 1
      zonesTotal         = zonesTotal + zonesPerSet + 1
    enddo

    do setID=nSetsP+1,nZoneSets
      self% zone1(setID) = zonesTotal + 1
      self% zone2(setID) = zonesTotal + zonesPerSet
      zonesTotal         = zonesTotal + zonesPerSet
    enddo

    if (zonesTotal /= Size% nzones) then
      call f90fatal("Construct Geometry: incorrect zone count in Zone sets")
    endif


    return
                                                                                             
  end subroutine Geometry_ctor

!=======================================================================
! set connectivity interface
!=======================================================================

  subroutine Geometry_setConnectivity(self,      &
                                      zone,      &
                                      nCorner,   &
                                      c0,        &
                                      zoneFaces, &
                                      zoneOpp,   &
                                      CToFace,   &
                                      cFP,       &
                                      cEZ,       &
                                      nCFacesArray)

    use Size_mod

    implicit none

!   Passed variables

    type(Geometry),  intent(inout) :: self

    integer, intent(in)            :: zone
    integer, intent(in)            :: nCorner
    integer, intent(in)            :: c0
    integer, intent(in)            :: zoneFaces
    integer, intent(in)            :: zoneOpp(zoneFaces)
    integer, intent(in)            :: CToFace(Size%maxcf,Size%maxCorner)
    integer, intent(in)            :: cFP(Size%maxcf,Size%maxCorner)
    integer, intent(in)            :: cEZ(Size%maxcf,Size%maxCorner)
    integer, intent(in)            :: nCFacesArray(Size%maxCorner)

!   Local

    integer                        :: c, cface, face

    self% cOffSet(zone)   = c0
    self% numCorner(zone) = nCorner
    self% zoneFaces(zone) = zoneFaces

    do c=1,nCorner
      self% CToZone(c0+c) = zone
    enddo

    do face=1,zoneFaces
      self% zoneOpp(face,zone) = zoneOpp(face)
    enddo

!   Set corner connectivity

    do c=1,nCorner
      self% nCFacesArray(c0+c) = nCFacesArray(c)
      do cface=1,Size%maxcf
        self% cFP(cface,c0+c)     = cFP(cface,c)
        self% cEZ(cface,c0+c)     = cEZ(cface,c)
        self% CToFace(cface,c0+c) = CToFace(cface,c)
      enddo
    enddo


    return

  end subroutine Geometry_setConnectivity

!=======================================================================
! getZoneAverage interface
!=======================================================================
  function Geometry_getZoneAverage(self, zone, cornerVar) &
                                   result(ZoneAverage)

     use constant_mod
     use Size_mod

!   Return the zone average of a corner variable 

!    variable declarations
     implicit none

!    passed variables
     type(Geometry), intent(in) :: self
     integer,        intent(in) :: zone
     real(adqt),     intent(in) :: cornerVar(Size% ncornr)
     real(adqt)                 :: ZoneAverage

!    Local
     integer  :: c
     integer  :: c0
     integer  :: nCorner

     ZoneAverage = zero
     nCorner     = self% numCorner(zone)
     c0          = self% cOffSet(zone)

     do c=1,nCorner
       ZoneAverage = ZoneAverage + cornerVar(c0+c)*self% Volume(c0+c)
     enddo

     ZoneAverage = ZoneAverage/self% VolumeZone(zone)


     return

  end function Geometry_getZoneAverage

!=======================================================================
! getZoneCenter interface
!=======================================================================

  function Geometry_getZoneCenter(self, zone) result(zoneCenter)

!    Return the zone center

     use Size_mod
     use constant_mod

     implicit none

!    passed variables
     type(Geometry), intent(in) :: self
     integer,        intent(in) :: zone
     real(adqt)                 :: zoneCenter(Size% ndim)

!    Local
     integer  :: c
     integer  :: c0
     integer  :: nCorner

!    Accumulate sum of coordinates:

     nCorner = self% numCorner(zone)
     c0      = self% cOffSet(zone)

     zoneCenter(:) = zero

     do c=1,nCorner
       zoneCenter(:) = zoneCenter(:) + self% px(:,c0+c)
     enddo

!    Divide by number of corners to get average coordinate:

     zoneCenter(:) = zoneCenter(:)/real(nCorner,adqt)


     return

  end function Geometry_getZoneCenter

!=======================================================================
! getFaceCenter interface
!=======================================================================

  function Geometry_getFaceCenter(self, zone, nFaces) result(faceCenter)

!    Return the zone center

     use Size_mod
     use constant_mod

     implicit none

!    passed variables
     type(Geometry), intent(in) :: self
     integer,        intent(in) :: zone
     integer,        intent(in) :: nFaces
     real(adqt)                 :: faceCenter(3,nFaces)

!    Local
     integer  :: c
     integer  :: c0
     integer  :: nCorner
     integer  :: cface
     integer  :: face
     integer  :: nCFaces
     integer  :: nc_face(nFaces)

!    Sum all point coordinates associated with a face and
!    the number of corner faces on a zone face

     nCorner         = self% numCorner(zone)
     c0              = self% cOffSet(zone)
     nc_face(:)      = 0
     faceCenter(:,:) = zero

     do c=1,nCorner
       nCFaces = self% nCFacesArray(c0+c)
       do cface=1,nCFaces
         face               = self% CToFace(cface,c0+c)
         nc_face(face)      = nc_face(face) + 1
         faceCenter(:,face) = faceCenter(:,face) + self% px(:,c0+c)
       enddo
     enddo

!    Calculate the face-center coordinates as the average of point
!    coordinates associated with that face

     do face=1,nFaces
       if (nc_face(face) /= 0) then
         faceCenter(:,face) = faceCenter(:,face)/real(nc_face(face),adqt)
       else
         faceCenter(:,face) = zero
       endif
     enddo

     return

  end function Geometry_getFaceCenter
                                                                                             
!=======================================================================
! destruct interface
!=======================================================================

  subroutine Geometry_dtor(self)

    use Size_mod
    use MemoryAllocator_mod

    implicit none

!   Passed variables

    type(Geometry),  intent(inout) :: self

!   Pointers 

    deallocate( self% VoC )
    deallocate( self% px  )
    deallocate( self% zoneFaces )
    deallocate( self% zoneOpp )
    deallocate( self% faceOpp )
    deallocate( self% CToFace )
    deallocate( self% BoundaryZone )

    call Allocator%deallocate(Size%usePinnedMemory,self%label,"Volume",       self% Volume)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"VolumeOld",    self% VolumeOld)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"VolumeZone",   self% VolumeZone)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"cOffSet",      self% cOffSet)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"numCorner",    self% numCorner)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"CToZone",      self% CToZone)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"zone1",        self% zone1)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"zone2",        self% zone2)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"corner1",      self% corner1)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"corner2",      self% corner2)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"cFP",          self% cFP)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"cEZ",          self% cEZ)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"A_fp",         self% A_fp)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"A_ez",         self% A_ez)
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"nCFacesArray", self% nCFacesArray)

    if (Size% ndim == 1) then

      deallocate( self% Area )
      deallocate( self% Radius )
      deallocate( self% zoneWidth )
      deallocate( self% Rave )
      deallocate( self% Rmin )
      deallocate( self% Rmax )

    elseif (Size% ndim == 2) then

      call Allocator%deallocate(Size%usePinnedMemory,self%label,"Area",     self% Area)
      call Allocator%deallocate(Size%usePinnedMemory,self%label,"RadiusEZ", self% RadiusEZ)
      call Allocator%deallocate(Size%usePinnedMemory,self%label,"RadiusFP", self% RadiusFP)

    endif

!   PWLD
    if ( Size% usePWLD ) then
      if ( Size% useSurfaceMassLumping ) then
        deallocate( self% LL )
        deallocate( self% AD )
        deallocate( self% MM )
      else
        deallocate( self% LLu )
        deallocate( self% ADu )
        deallocate( self% MMu )
      endif
    endif


    return

  end subroutine Geometry_dtor

end module Geometry_mod
