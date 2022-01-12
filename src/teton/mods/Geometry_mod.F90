! Geometry Modules:  geometry and connectivity information used by Sn
                                                                                 
module Geometry_mod 

  use kind_mod
  use ZoneData_mod
  use MeshData_mod

  private

! public interfaces
                                                                                             
  public construct 
  public destruct 
  public getZoneData
  public getMesh
  public getZoneAverage
  public setConnectivity
                                                                                 
                                                                                 
  type, public :: Geometry 

!    Connectivity 

     integer,         pointer, contiguous :: cEZ(:,:) => null()
     integer,         pointer, contiguous :: cFP(:,:) => null()
     integer,         pointer, contiguous :: nCFacesArray(:) => null()
     integer,         pointer, contiguous :: cOffSet(:) => null()
     integer,         pointer, contiguous :: numCorner(:) => null()
     integer,         pointer, contiguous :: cToZone(:) => null()
     integer,         pointer, contiguous :: zone1(:) => null()
     integer,         pointer, contiguous :: zone2(:) => null()
     integer,         pointer, contiguous :: corner1(:) => null()
     integer,         pointer, contiguous :: corner2(:) => null()

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

! Placeholder

     real(adqt),      pointer, contiguous :: PhiTotal(:,:) => null()
     real(adqt),      pointer, contiguous :: CollisionRate(:) => null()
     real(adqt),      pointer, contiguous :: radEnergy(:) => null()
     real(adqt),      pointer, contiguous :: RadEnergyDensity(:,:) => null()

! PWLD

     real(adqt),      pointer, contiguous :: LL(:,:,:,:) => null()
     real(adqt),      pointer, contiguous :: AD(:,:,:) => null()
     real(adqt),      pointer, contiguous :: MM(:,:) => null()

! Unlumped PWLD

     real(adqt),      pointer, contiguous :: LLu(:,:,:,:) => null()
     real(adqt),      pointer, contiguous :: ADu(:,:,:) => null()
     real(adqt),      pointer, contiguous :: MMu(:,:,:) => null()

! Pointers to other modules

     type(ZoneData), pointer :: ZData(:)   => null()   ! zone data pointers
     type(MeshData), pointer :: MData(:)   => null()   ! mesh data pointers

!    Misc
     character(len=8) :: label ! A string descriptor for this set.

  end type Geometry

  type(Geometry), pointer, public :: Geom => null()

  interface construct
    module procedure Geometry_ctor
  end interface

  interface getZoneData
    module procedure Geometry_getZone
  end interface

  interface getMesh
    module procedure Geometry_getMesh
  end interface

  interface getZoneAverage
    module procedure Geometry_getZoneAverage
  end interface

  interface setConnectivity 
    module procedure Geometry_setConnectivity
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
    logical                        :: usePinnedMemory

    usePinnedMemory = Size%useGPU

    self%label = "geometry"

!   Pointers

    allocate( self% ZData(Size%nzones) )
    allocate( self% MData(Size%nzones) )

    allocate( self% VolumeZone(Size%nzones) )
    self% VolumeZone(:) = zero
    allocate( self% VoC(Size% ndim,Size% ncornr) )
    self% VoC(:,:)      = zero

    call Allocator%allocate(usePinnedMemory,self%label,"Volume",       self% Volume,    Size% ncornr)
    call Allocator%allocate(usePinnedMemory,self%label,"VolumeOld",    self% VolumeOld, Size% ncornr)
    call Allocator%allocate(usePinnedMemory,self%label,"cOffSet",      self% cOffSet,   Size% nzones)
    call Allocator%allocate(usePinnedMemory,self%label,"numCorner",    self% numCorner, Size% nzones)
    call Allocator%allocate(usePinnedMemory,self%label,"cToZone",      self% cToZone,   Size% ncornr)
    call Allocator%allocate(usePinnedMemory,self%label,"zone1",        self% zone1,     nZoneSets)
    call Allocator%allocate(usePinnedMemory,self%label,"zone2",        self% zone2,     nZoneSets)
    call Allocator%allocate(usePinnedMemory,self%label,"corner1",      self% corner1,   nZoneSets)
    call Allocator%allocate(usePinnedMemory,self%label,"corner2",      self% corner2,   nZoneSets)
    call Allocator%allocate(usePinnedMemory,self%label,"cFP",          self% cFP,       Size% maxcf,Size% ncornr)
    call Allocator%allocate(usePinnedMemory,self%label,"cEZ",          self% cEZ,       Size% maxcf,Size% ncornr)
    call Allocator%allocate(usePinnedMemory,self%label,"A_fp",         self% A_fp,      Size% ndim,Size% maxcf,Size% ncornr)
    call Allocator%allocate(usePinnedMemory,self%label,"A_ez",         self% A_ez,      Size% ndim,Size% maxcf,Size% ncornr)
    call Allocator%allocate(usePinnedMemory,self%label,"nCFacesArray", self% nCFacesArray, Size% ncornr)

    if (Size% ndim == 1) then
      allocate( self% Area(Size% ncornr) )
      allocate( self% Radius(Size% ncornr) )

      self% Area(:)   = zero
      self% Radius(:) = zero

    elseif (Size% ndim == 2) then

      call Allocator%allocate(usePinnedMemory,self%label,"Area",     self% Area,      Size% ncornr)
      call Allocator%allocate(usePinnedMemory,self%label,"RadiusEZ", self% RadiusEZ,  2,Size% ncornr)
      call Allocator%allocate(usePinnedMemory,self%label,"RadiusFP", self% RadiusFP,  2,Size% ncornr)

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

    call Allocator%allocate(usePinnedMemory,self%label,"PhiTotal", self% PhiTotal, Size% ngr, Size% ncornr)

    allocate( self% CollisionRate(Size% ncornr) )
    self% CollisionRate(:)      = zero

    call Allocator%allocate(.FALSE.,self%label,"radEnergy", self% radEnergy, Size%nZones)
    call Allocator%allocate(.FALSE.,self%label,"RadEnergyDensity", self% RadEnergyDensity, Size%nzones, Size%ngr)

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

  subroutine Geometry_setConnectivity(self,     &
                                      zone,     &
                                      nCorner,  &
                                      c0,       &
                                      cFP,      &
                                      cEZ,      &
                                      nCFacesArray)

    use Size_mod

    implicit none

!   Passed variables

    type(Geometry),  intent(inout) :: self

    integer, intent(in)            :: zone
    integer, intent(in)            :: nCorner
    integer, intent(in)            :: c0
    integer, intent(in)            :: cFP(Size%maxcf,Size%maxCorner)
    integer, intent(in)            :: cEZ(Size%maxcf,Size%maxCorner)
    integer, intent(in)            :: nCFacesArray(Size%maxCorner)

!   Local

    integer                        :: c, cface

    self% cOffSet(zone)   = c0
    self% numCorner(zone) = nCorner

    do c=1,nCorner
      self% cToZone(c0+c) = zone
    enddo

!   Set corner connectivity

    do c=1,nCorner
      self% nCFacesArray(c0+c) = nCFacesArray(c)
      do cface=1,Size%maxcf
        self% cFP(cface,c0+c) = cFP(cface,c)
        self% cEZ(cface,c0+c) = cEZ(cface,c)
      enddo
    enddo


    return

  end subroutine Geometry_setConnectivity

!=======================================================================
! get Zone Data interface
!=======================================================================
  function Geometry_getZone(self,zoneID) result(ZData)

!    Return a pointer to a zone definition

!    variable declarations
     implicit none

!    passed variables
     type(Geometry),     intent(in)   :: self
     integer,            intent(in)   :: zoneID
     type(ZoneData),     pointer      :: ZData

     ZData => self % ZData(zoneID)

     return

  end function Geometry_getZone

!=======================================================================
! getMesh interface
!=======================================================================
  function Geometry_getMesh(self,zoneID) result(MData)
 
!    Return a pointer to a zone definition
 
!    variable declarations
     implicit none
 
!    passed variables
     type(Geometry),     intent(in)   :: self
     integer,            intent(in)   :: zoneID
     type(MeshData),     pointer      :: MData
 
     MData => self % MData(zoneID)
 
     return
 
  end function Geometry_getMesh

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
! destruct interface
!=======================================================================

  subroutine Geometry_dtor(self)

    use Size_mod
    use MemoryAllocator_mod

    implicit none

!   Passed variables

    type(Geometry),  intent(inout) :: self

    logical                        :: usePinnedMemory
    usePinnedMemory = Size%useGPU

!   Pointers 

    deallocate( self% ZData )
    deallocate( self% MData )

    deallocate( self% VolumeZone )
    deallocate( self% VoC )

    call Allocator%deallocate(usePinnedMemory,self%label,"Volume",       self% Volume)
    call Allocator%deallocate(usePinnedMemory,self%label,"VolumeOld",    self% VolumeOld)
    call Allocator%deallocate(usePinnedMemory,self%label,"cOffSet",      self% cOffSet)
    call Allocator%deallocate(usePinnedMemory,self%label,"numCorner",    self% numCorner)
    call Allocator%deallocate(usePinnedMemory,self%label,"cToZone",      self% cToZone)
    call Allocator%deallocate(usePinnedMemory,self%label,"zone1",        self% zone1)
    call Allocator%deallocate(usePinnedMemory,self%label,"zone2",        self% zone2)
    call Allocator%deallocate(usePinnedMemory,self%label,"corner1",      self% corner1)
    call Allocator%deallocate(usePinnedMemory,self%label,"corner2",      self% corner2)
    call Allocator%deallocate(usePinnedMemory,self%label,"cFP",          self% cFP)
    call Allocator%deallocate(usePinnedMemory,self%label,"cEZ",          self% cEZ)
    call Allocator%deallocate(usePinnedMemory,self%label,"A_fp",         self% A_fp)
    call Allocator%deallocate(usePinnedMemory,self%label,"A_ez",         self% A_ez)
    call Allocator%deallocate(usePinnedMemory,self%label,"nCFacesArray", self% nCFacesArray)

    if (Size% ndim == 1) then

      deallocate( self% Area )
      deallocate( self% Radius )

    elseif (Size% ndim == 2) then

      call Allocator%deallocate(usePinnedMemory,self%label,"Area",     self% Area)
      call Allocator%deallocate(usePinnedMemory,self%label,"RadiusEZ", self% RadiusEZ)
      call Allocator%deallocate(usePinnedMemory,self%label,"RadiusFP", self% RadiusFP)

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

    call Allocator%deallocate(usePinnedMemory,self%label,"PhiTotal", self% PhiTotal)

    deallocate( self% CollisionRate )

    call Allocator%deallocate(.FALSE.,self%label,"radEnergy", self% radEnergy)
    call Allocator%deallocate(.FALSE.,self%label,"RadEnergyDensity", self% RadEnergyDensity)

    return

  end subroutine Geometry_dtor

end module Geometry_mod
