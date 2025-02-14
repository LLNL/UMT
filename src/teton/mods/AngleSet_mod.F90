#include "macros.h"
! Angle Set Module:  Contains data structures for a group set 

module AngleSet_mod

  use kind_mod
  use constant_mod
  use Quadrature_mod
  use Size_mod

  USE ISO_C_BINDING

  private

! public interfaces

  public construct
  public destruct
  public constructHyperPlane
  public getNumberOfHyperPlanes
  public getZonesInPlane
  public destructHyperPlane
  public getReflectedAngle
  public setReflectedAngle
  public constructBdyExitList
  public destructBdyExitList
  public constructCycleList
  public destructCycleList
  public constructIncidentTest
  public getIncidentTest
  public destructIncidentTest


  type, public :: AngleSet

     integer                              :: NumAngles    ! number of angles 
     integer                              :: nPolarAngles
     integer                              :: Order
     integer                              :: NumBin
     integer                              :: NumBin0
     integer                              :: angle0       ! angle offset
     integer                              :: totalCycles
     integer                              :: maxInterface
     logical (kind=1)                     :: GTASet   

!    Quadrature related
     integer,         pointer, contiguous :: nExit(:) => null()
     integer,         pointer, contiguous :: ExitAngleList(:,:) => null()
     integer,         pointer, contiguous :: ReflectedAngle(:,:) => null()
     integer,         pointer, contiguous :: NangBinList(:) => null()
     integer,         pointer, contiguous :: AngleToBin(:) => null()
     integer,         pointer, contiguous :: PolarAngle(:) => null()

     real(adqt),      pointer, contiguous :: PolarAngleList(:) => null()
     real(kind=C_DOUBLE),      pointer, contiguous :: Omega(:,:) => null()
     real(kind=C_DOUBLE),      pointer, contiguous :: Weight(:) => null()

!    Angular coefficients used in 1D and 2D
     real(adqt),      pointer, contiguous :: Alpha(:) => null()
     real(adqt),      pointer, contiguous :: Tau(:) => null()
     real(adqt),      pointer, contiguous :: quadTauW1(:) => null()
     real(adqt),      pointer, contiguous :: quadTauW2(:) => null()
     real(adqt),      pointer, contiguous :: angDerivFac(:) => null()
     real(adqt),      pointer, contiguous :: Falpha(:) => null()
     real(adqt),      pointer, contiguous :: Adweta(:) => null()

!    Outward Normals

     real(adqt),      pointer, contiguous :: AfpNorm(:,:) => null()
     real(adqt),      pointer, contiguous :: AezNorm(:,:) => null()
     real(adqt),      pointer, contiguous :: ANormSum(:) => null()

     logical(kind=C_BOOL), pointer, contiguous :: StartingDirection(:) => null()
     logical(kind=C_BOOL), pointer, contiguous :: FinishingDirection(:) => null()

!    For transport Sweeps

     integer(kind=C_INT),  pointer, contiguous :: nextZ(:,:) => null()
     integer(kind=C_INT),  pointer, contiguous :: nextC(:,:) => null()
     integer(kind=C_INT),  pointer, contiguous :: nHyperPlanes(:) => null()
     integer(kind=C_INT),  pointer, contiguous :: numCycles(:) => null()
     integer(kind=C_INT),  pointer, contiguous :: cycleOffSet(:) => null()
     integer(kind=C_INT),  pointer, contiguous :: cycleList(:) => null()

!    Pointers 

     type(HypPlane),  pointer             :: HypPlanePtr(:) => null()
     type(BdyExit),   pointer             :: BdyExitPtr(:) => null()
     type(IncidentTest), pointer          :: TestPtr(:) => null()

!    Misc
     character(len=22) :: label ! A string descriptor for this set.

  end type AngleSet 

  type, public :: HypPlane
     integer                              :: maxZones
     integer                              :: maxCorners
     integer                              :: interfaceLen
     integer,         pointer, contiguous :: zonesInPlane(:) => null()
     integer,         pointer, contiguous :: cornersInPlane(:) => null()
     integer,         pointer, contiguous :: badCornerList(:) => null()
     integer,         pointer, contiguous :: interfaceList(:) => null()
     integer,         pointer, contiguous :: hplane1(:) => null()
     integer,         pointer, contiguous :: hplane2(:) => null()
     integer,         pointer, contiguous :: ndone(:) => null()
     integer,         pointer, contiguous :: c1(:) => null()
     integer,         pointer, contiguous :: c2(:) => null()
  end type HypPlane

  type, public :: BdyExit
     integer                         :: nxBdy
     integer,    pointer, contiguous :: bdyList(:,:) => null()
  end type BdyExit

  type, public :: IncidentTest
     ! Sizes of IncTest and IncTestR, respectively
     integer                         :: nSend
     integer                         :: nRecv

     integer                         :: request(2)

! If these two arrays are initialized to null(), IBM XLF will throw FPE's when
! the AngleSet's TestPtr is allocated.
#if defined(__ibmxl__)
     integer,    pointer, contiguous :: IncTest(:)
     integer,    pointer, contiguous :: IncTestR(:)
#else
     integer,    pointer, contiguous :: IncTest(:) => null()
     integer,    pointer, contiguous :: IncTestR(:) => null()
#endif
  end type IncidentTest

  interface construct
    module procedure AngleSet_ctor
  end interface

  interface destruct
    module procedure AngleSet_dtor
  end interface

  interface constructHyperPlane
    module procedure AngleSet_ctorHypPlane
  end interface

  interface getNumberOfHyperPlanes
    module procedure AngleSet_getNumHypPlane
  end interface

  interface getZonesInPlane
    module procedure AngleSet_getZonesInPlane
  end interface

  interface destructHyperPlane
    module procedure AngleSet_dtorHypPlane
  end interface

  interface getReflectedAngle
    module procedure AngleSet_getReflAngle
  end interface

  interface setReflectedAngle
    module procedure AngleSet_setReflAngle
  end interface

  interface constructBdyExitList
    module procedure AngleSet_ctorBdyExitList
  end interface

  interface destructBdyExitList
    module procedure AngleSet_dtorBdyExitList
  end interface

  interface constructCycleList
    module procedure AngleSet_ctorCycleList
  end interface

  interface destructCycleList
    module procedure AngleSet_dtorCycleList
  end interface

  interface constructIncidentTest
    module procedure AngleSet_ctorIncTest
  end interface

  interface getIncidentTest
    module procedure AngleSet_getIncTest
  end interface

  interface destructIncidentTest
    module procedure AngleSet_dtorIncTest
  end interface


contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine AngleSet_ctor(self,                 &
                           NumAngles,            &
                           angle0,               &
                           nZones,               &
                           nReflecting,          &
                           GTASet,               &
                           QuadPtr)

    use Size_mod
    use constant_mod
    use MemoryAllocator_mod

    implicit none

!   Passed variables

    type(AngleSet), intent(inout)        :: self

    integer,         intent(in)          :: NumAngles
    integer,         intent(in)          :: angle0
    integer,         intent(in)          :: nZones
    integer,         intent(in)          :: nReflecting
    logical(kind=1), intent(in)          :: GTASet

    type(Quadrature), target, intent(in) :: QuadPtr

!   Local

    integer :: n
    integer :: ndim
    integer :: angle
    integer :: nLevels

    write(self%label,'(I0.3)') angle0
    self%label = QuadPtr%label//"_angle_set_"//self%label

!   Set Properties
    ndim               = Size% ndim
    self% NumAngles    = NumAngles
    self% nPolarAngles = QuadPtr% nPolarAngles
    self% angle0       = angle0
    self% Order        = QuadPtr% Order
    self% maxInterface = 1
    self% GTASet       = GTASet

!   Allocate Memory

    if (Size%ncomm > 0) then
      allocate( self% TestPtr(Size%ncomm) )
    endif

    allocate( self% nExit(nReflecting) )
    allocate( self% ExitAngleList(self% NumAngles,nReflecting) )
    allocate( self% ReflectedAngle(self% NumAngles, nReflecting) )
    allocate( self% AngleToBin(self% NumAngles) )

    allocate( self% PolarAngle(self% NumAngles) )
    allocate( self% PolarAngleList(self% nPolarAngles) )
    allocate( self% Omega(ndim,self% NumAngles) )
    allocate( self% Weight(self% NumAngles) )
    allocate( self% StartingDirection(self% NumAngles) )
    allocate( self% FinishingDirection(self% NumAngles) )

    self% nExit(:)            =  0
    self% ExitAngleList(:,:)  =  0
    self% ReflectedAngle(:,:) = -1
    self% PolarAngleList(:)   = QuadPtr% PolarAngleList(:)

    if (ndim > 1) then
      call Allocator%allocate(Size%usePinnedMemory,self%label,"nextZ",    self% nextZ,    nzones,      self% NumAngles)
      call Allocator%allocate(Size%usePinnedMemory,self%label,"nextC",    self% nextC,    Size%ncornr, self% NumAngles)

      if ( .not. self% GTASet ) then
        call Allocator%allocate(Size%usePinnedMemory,self%label,"AfpNorm",  self% AfpNorm,  Size%maxcf,  Size%ncornr)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"AezNorm",  self% AezNorm,  Size%maxcf,  Size%ncornr)
        call Allocator%allocate(Size%usePinnedMemory,self%label,"ANormSum", self% ANormSum, Size%ncornr)
      endif
    endif

    nLevels = 0

    do n=1,NumAngles

      angle               = QuadPtr% angleList(angle0+n)
      self% PolarAngle(n) = QuadPtr% PolarAngle(angle)

      self% Omega(:,n)            = QuadPtr% Omega(:,angle)
      self% Weight(n)             = QuadPtr% Weight(angle)

      self% StartingDirection(n)  = QuadPtr% StartingDirection(angle)
      self% FinishingDirection(n) = QuadPtr% FinishingDirection(angle)

      if ( self% StartingDirection(n) ) then
        nLevels = nLevels + 1
      endif

      self% AngleToBin(n) = n

    enddo

!   Angular Coefficents

    if (ndim == 1) then

      self% NumBin  = NumAngles
      self% NumBin0 = NumAngles

      allocate( self% NangBinList(self% NumBin) )
      allocate( self% Alpha(self% NumAngles) )
      allocate( self% Tau(self% NumAngles) )
      allocate( self% Falpha(self% NumAngles) )
      allocate( self% Adweta(self% NumAngles) )

      self% NangBinList(:) = 1

      call AngleCoef1D(self)

    elseif (ndim == 2) then

      self% NumBin  = nLevels
      self% NumBin0 = nLevels

      allocate( self% NangBinList(self% NumBin) )
      allocate( self% Alpha(self% NumAngles) )
      allocate( self% Tau(self% NumAngles) )
      allocate( self% quadTauW1(self% NumAngles) )
      allocate( self% quadTauW2(self% NumAngles) )
      allocate( self% angDerivFac(self% NumAngles) )

      nLevels              = 0
      self% NangBinList(:) = 0

      do n=1,NumAngles
        if (self% StartingDirection(n)) then
          nLevels = nLevels + 1
        endif

        self% AngleToBin(n)        = nLevels
        self% NangBinList(nLevels) = self% NangBinList(nLevels) + 1
      enddo

      call AngleCoef2D(self)

      do n=1,NumAngles
        if (self% StartingDirection(n) .or. self% FinishingDirection(n)) then
          self% angDerivFac(n) = zero
          self% quadTauW1(n)   = one
          self% quadTauW2(n)   = zero
        else
          self% angDerivFac(n) = self% Omega(1,n) +   &
                                 self% Alpha(n)/(self% Weight(n)*self% Tau(n))
          self% quadTauW1(n)   = one/self% Tau(n)
          self% quadTauW2(n)   = (one - self% Tau(n))/self% Tau(n)
        endif
      enddo

    elseif (ndim == 3) then
      self% NumBin  = NumAngles
      self% NumBin0 = NumAngles

      allocate( self% NangBinList(self% NumBin) )

      self% NangBinList(:) = 1

    endif

!   Hyperplanes

    if (Size% ndim > 1) then
      allocate( self% numCycles(self% NumAngles) )
      allocate( self% cycleOffSet(self% NumAngles) )
      allocate( self% nHyperPlanes(self% NumAngles) )
      allocate( self% HypPlanePtr(self% NumAngles) )
      allocate( self% BdyExitPtr(self% NumAngles) )

      self% totalCycles     = 0
      self% numCycles(:)    = 0
      self% cycleOffSet(:)  = 0
      self% nHyperPlanes(:) = 0

      do n=1,NumAngles
        self% HypPlanePtr(n)% interfaceLen = 0
      enddo
    endif


    return

  end subroutine AngleSet_ctor

!=======================================================================
! destruct interface
!=======================================================================

  subroutine AngleSet_dtor(self)

    use Size_mod
    use MemoryAllocator_mod

    implicit none

!   Passed variables

    type(AngleSet), intent(inout)    :: self

    if (Size% ndim > 1) then
      call Allocator%deallocate(Size%usePinnedMemory,self%label,"nextZ",    self% nextZ)
      call Allocator%deallocate(Size%usePinnedMemory,self%label,"nextC",    self% nextC)

      if ( .not. self% GTASet ) then
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"AfpNorm",  self% AfpNorm)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"AezNorm",  self% AezNorm)
        call Allocator%deallocate(Size%usePinnedMemory,self%label,"ANormSum", self% ANormSum)
      endif
    endif

    deallocate( self% nExit )
    deallocate( self% ExitAngleList )
    deallocate( self% ReflectedAngle )
    deallocate( self% NangBinList )
    deallocate( self% AngleToBin )

    if (Size% ncomm > 0) then
      deallocate( self% TestPtr )
    endif

    deallocate( self% PolarAngle )
    deallocate( self% PolarAngleList )
    deallocate( self% Omega )
    deallocate( self% Weight )
    deallocate( self% StartingDirection )
    deallocate( self% FinishingDirection )

!   Space for angular coefficients

    if (Size% ndim == 1) then
      deallocate( self% Alpha )
      deallocate( self% Tau )
      deallocate( self% Falpha )
      deallocate( self% Adweta )
    else if (Size% ndim == 2) then
      deallocate( self% Alpha )
      deallocate( self% Tau )
      deallocate( self% quadTauW1 )
      deallocate( self% quadTauW2 )
      deallocate( self% angDerivFac )
    endif

!   Hyperplanes

    if (Size% ndim > 1) then
      deallocate( self% numCycles )
      deallocate( self% cycleOffSet )
      deallocate( self% nHyperPlanes )
      deallocate( self% HypPlanePtr )
      deallocate( self% BdyExitPtr )
    endif


    return

  end subroutine AngleSet_dtor

!=======================================================================
! construct HyperPlane
!=======================================================================
  subroutine AngleSet_ctorHypPlane(self, angle, nHyperPlanes, meshCycles, &
                                   nHyperDomains, nZoneSets,              &
                                   elementsInPlane, CToHypPlane, cycleList)

    use Size_mod
    use Geometry_mod
    use MemoryAllocator_mod
    use Options_mod

    implicit none

!   Passed variables
    type(AngleSet),  intent(inout) :: self
    integer,         intent(in)    :: angle
    integer,         intent(in)    :: nHyperPlanes
    integer,         intent(in)    :: meshCycles
    integer,         intent(in)    :: nHyperDomains
    integer,         intent(in)    :: nZoneSets

    integer,         intent(in)    :: elementsInPlane(nHyperPlanes)
    integer,         intent(in)    :: CToHypPlane(Size% ncornr) 
    integer,         intent(in)    :: cycleList(meshCycles)

!   Local
    type(HypPlane),  pointer       :: HypPlanePtr

    integer                        :: sweepVersion
    integer                        :: hPlane
    integer                        :: hPlane1
    integer                        :: hPlane2
    integer                        :: mCycle
    integer                        :: elementSum
    integer                        :: elementsPerDomain
    integer                        :: elementTarget
    integer                        :: domID
    integer                        :: numC
    integer                        :: c
    integer                        :: c0
    integer                        :: cfp
    integer                        :: cez
    integer                        :: cface
    integer                        :: nCorner
    integer                        :: nCFaces
    integer                        :: nCFacesEZ
    integer                        :: i
    integer                        :: ii
    integer                        :: zone
    integer                        :: nzones
    integer                        :: ndone
    integer                        :: hp
    integer                        :: setID
    integer                        :: cornersPerSet
    integer                        :: cornersTotal
    integer                        :: nSetsP

    real(adqt)                     :: afp
    real(adqt)                     :: aez

    integer,          allocatable  :: cornerList(:)
    logical (kind=1), allocatable  :: done(:)

    logical(kind=1)                :: notDone

!   Allocate Memory 

! Corner sweep is not yet implemented for grey sweep.
    if ( self% GTASet ) then
      sweepVersion = 1
    else
      sweepVersion = Options% getSweepVersion()
    endif

    HypPlanePtr  => self% HypPlanePtr(angle)

    if ( sweepVersion == 1 ) then
      call Allocator%allocate(Size%usePinnedMemory, self%label, "zonesInPlane", HypPlanePtr% zonesInPlane, nHyperPlanes)

      elementsPerDomain            = Size% nzones/nHyperDomains
      HypPlanePtr% zonesInPlane(:) = elementsInPlane(:)
      HypPlanePtr% maxZones        = maxval( elementsInPlane(1:nHyperPlanes) )
    else if (sweepVersion == 2 ) then
      call Allocator%allocate(Size%usePinnedMemory, self%label, "cornersInPlane", HypPlanePtr% cornersInPlane, nHyperPlanes)

      elementsPerDomain              = Size% ncornr/nHyperDomains
      HypPlanePtr% cornersInPlane(:) = elementsInPlane(:)
      HypPlanePtr% maxCorners        = maxval( elementsInPlane(1:nHyperPlanes) )
    endif

    call Allocator%allocate(Size%usePinnedMemory, self%label, "hplane1", HypPlanePtr% hplane1, nHyperDomains+1)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "hplane2", HypPlanePtr% hplane2, nHyperDomains)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "ndone",   HypPlanePtr% ndone,   nHyperDomains+1)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "c1",      HypPlanePtr% c1,      nZoneSets)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "c2",      HypPlanePtr% c2,      nZoneSets)

    allocate( HypPlanePtr% badCornerList(meshCycles+1) )
    allocate( cornerList(Size% ncornr) )
    allocate( done(Size% ncornr) )

    domID                 = 1
    elementSum            = 0
    elementTarget         = elementsPerDomain
    notDone               = .TRUE.
    HypPlanePtr% ndone(1) = 0

!   If nHyperDomains=1 there is no interface list, but it is given a length of 1
!   because it is mapped. In this case, it uses cornerList(1) which may be
!   unset. We set it here.

    cornerList(1) = 1

    do hPlane=1,nHyperPlanes
      elementSum = elementSum + elementsInPlane(hPlane)

      if (elementSum > elementTarget .and. notDone) then

        HypPlanePtr% ndone(domID+1)   = elementSum
        HypPlanePtr% hplane1(domID+1) = hPlane + 1
        HypPlanePtr% hplane2(domID)   = hPlane
        elementTarget                 = elementTarget + elementsPerDomain
        domID                         = domID + 1

        if (domID == nHyperDomains) then
          notDone = .False.
        endif
      endif
    enddo

    HypPlanePtr% hplane1(1)             = 1
    HypPlanePtr% hplane2(nHyperDomains) = nHyperPlanes

!   We need to store the upstream corner fluxes at the hyperdomain
!   boundaries so compute that size and corner list here, Note
!   that we do not count the first domain as there are no
!   upstream hyperdomains 

    numC    =  0
    done(:) = .FALSE.

    do domID=2,nHyperDomains
      hPlane1 = HypPlanePtr% hplane1(domID)
      hPlane2 = HypPlanePtr% hplane2(domID)
      ndone   = HypPlanePtr% ndone(domID)

      HyperPlaneLoop: do hPlane=hPlane1,hPlane2

        if ( sweepVersion == 1 ) then

          nzones = elementsInPlane(hPlane)

          ZoneLoop: do ii=1,nzones
            zone    = iabs( self% nextZ(ndone+ii,angle) )
            c0      = Geom% cOffSet(zone)
            nCorner = Geom% numCorner(zone)

            CornerLoop: do c=1,nCorner
              nCFaces = Geom% nCFacesArray(c0+c)

              CornerFaceLoop: do cface=1,nCFaces
                cfp = Geom% cFP(cface,c0+c)

!               Eliminate entries outside the mesh (cfp > ncornr)
                if ( cfp <= Size% ncornr ) then

                  hp  = CToHypPlane(cfp)

!                 Consider all corners in hyperplanes upstream from hplane1 
                  if (hp < hplane1) then

                    afp = DOT_PRODUCT( self% omega(:,angle),Geom% A_fp(:,cface,c0+c) )

                    if ( afp < zero ) then
                      if ( .not. done(cfp) ) then
                        numC             =  numC + 1
                        cornerList(numC) =  cfp
                        done(cfp)        = .TRUE.
                      endif
                    endif
                  endif
                endif
              enddo CornerFaceLoop
            enddo CornerLoop
          enddo ZoneLoop

          ndone = ndone + nzones

        else if ( sweepVersion == 2 ) then

          nCorner = elementsInPlane(hPlane)

          CornerLoop1: do ii=1,nCorner
            c       = self% nextC(ndone+ii,angle)
            zone    = Geom% CToZone(c)
            c0      = Geom% cOffSet(zone)
            nCFaces = Geom% nCFacesArray(c)

            CornerFaceLoop1: do cface=1,nCFaces
              cfp = Geom% cFP(cface,c)

!             Eliminate entries outside the mesh (cfp > ncornr)
              if ( cfp <= Size% ncornr ) then

                hp  = CToHypPlane(cfp)

!               Consider all corners in hyperplanes upstream from hplane1 
                if (hp < hplane1) then

                  afp = DOT_PRODUCT( self% omega(:,angle),Geom% A_fp(:,cface,c) )

                  if ( afp < zero ) then
                    if ( .not. done(cfp) ) then
                      numC             =  numC + 1
                      cornerList(numC) =  cfp
                      done(cfp)        = .TRUE.
                    endif
                  endif
                endif
              endif

!             For the corner sweep we also add corners in the same zone
              cez = c0 + Geom% cEZ(cface,c)
              hp  = CToHypPlane(cez)

!             Consider all corners in hyperplanes upstream from hplane1 
              if (hp < hplane1) then 
                aez = DOT_PRODUCT( self% omega(:,angle),Geom% A_ez(:,cface,c) )

                if ( aez < zero ) then
                  if ( .not. done(cez) ) then
                    numC             =  numC + 1
                    cornerList(numC) =  cez
                    done(cez)        = .TRUE.
                  endif

                  nCFacesEZ = Geom% nCFacesArray(cez)

                  do i=1,nCFacesEZ
                    cfp = Geom% cFP(i,cez)

                    if ( cfp <= Size% ncornr ) then

                      if (CToHypPlane(cfp) < hplane1) then

                        afp = DOT_PRODUCT( self% omega(:,angle),Geom% A_fp(:,i,cez) )

                        if ( afp < zero ) then
                          if ( .not. done(cfp) ) then
                            numC             =  numC + 1
                            cornerList(numC) =  cfp
                            done(cfp)        = .TRUE.
                          endif
                        endif
                      endif

                    endif
                  enddo

                endif
              endif
            enddo CornerFaceLoop1
          enddo CornerLoop1

          ndone = ndone + nCorner

        endif

      enddo HyperPlaneLoop
    enddo

!   Use the maximum to avoid mapping a zero-length array when nHyperDomains=1
    HypPlanePtr% interfaceLen = numC
    numC                      = max( numC,1 )

    call Allocator%allocate(Size%usePinnedMemory, self%label, "interfaceList", HypPlanePtr% interfaceList, numC)

    do ii=1,numC
      HypPlanePtr% interfaceList(ii) = cornerList(ii)
    enddo

    do mCycle=1,meshCycles
      HypPlanePtr% badCornerList(mCycle) = cycleList(mCycle)
    enddo

!   The interface list can be very large so divide up the list
!   Assign a range to each zone set; sets 1->nSetsP get an extra
!   element if there is a remainder

    if ( nHyperDomains > 1 ) then

      cornersPerSet = HypPlanePtr% interfaceLen/nZoneSets
      nSetsP        = HypPlanePtr% interfaceLen - (nZoneSets*cornersPerSet)
      cornersTotal  = 0

      do setID=1,nSetsP
        HypPlanePtr% c1(setID) = cornersTotal + 1
        HypPlanePtr% c2(setID) = cornersTotal + cornersPerSet + 1
        cornersTotal           = cornersTotal + cornersPerSet + 1
      enddo

      do setID=nSetsP+1,nZoneSets
        HypPlanePtr% c1(setID) = cornersTotal + 1
        HypPlanePtr% c2(setID) = cornersTotal + cornersPerSet
        cornersTotal           = cornersTotal + cornersPerSet
      enddo

    else

!     If there are no interface elements (nHyperDomains = 1) these are not used
      HypPlanePtr% c1(:) = 0
      HypPlanePtr% c2(:) = 0

    endif


    deallocate( cornerList )
    deallocate( done )


    return

  end subroutine AngleSet_ctorHypPlane

!=======================================================================
! get number of HyperPlanes    
!=======================================================================
  function AngleSet_getNumHypPlane(self, angle) result(nHyperPlanes)

    implicit none

!   Passed variables
    type(AngleSet),   intent(inout) :: self
    integer,          intent(in)    :: angle
    integer                         :: nHyperPlanes

!   Local

    nHyperPlanes =  self% nHyperPlanes(angle)

    return

  end function AngleSet_getNumHypPlane

!=======================================================================
! get zones in theHyperPlane    
!=======================================================================
  function AngleSet_getZonesInPlane(self, angle, hPlane) result(nZonesPlane)

    implicit none

!   Passed variables
    type(AngleSet),   intent(inout) :: self
    integer,          intent(in)    :: angle
    integer,          intent(in)    :: hPlane
    integer                         :: nZonesPlane

!   Local
    type(HypPlane),   pointer       :: HypPlanePtr


    HypPlanePtr => self% HypPlanePtr(angle)

    nZonesPlane =  HypPlanePtr% zonesInPlane(hPlane)

    return

  end function AngleSet_getZonesInPlane

!=======================================================================
! destruct HyperPlane    
!=======================================================================
  subroutine AngleSet_dtorHypPlane(self, sweepVersion)

    use MemoryAllocator_mod
    implicit none

!   Passed variables
    type(AngleSet),   intent(inout) :: self
    integer,          intent(in)    :: sweepVersion

!   Local
    type(HypPlane),   pointer       :: HypPlanePtr

    integer                         :: angle

!   Release Memory 

    do angle=1,self% NumAngles
      if ( .not. self% FinishingDirection(angle) ) then

        HypPlanePtr => self% HypPlanePtr(angle)

        if ( sweepVersion == 1 ) then
          call Allocator%deallocate(Size%usePinnedMemory, self%label, "zonesInPlane",   HypPlanePtr% zonesInPlane)
        else if ( sweepVersion == 2 ) then
          call Allocator%deallocate(Size%usePinnedMemory, self%label, "cornersInPlane", HypPlanePtr% cornersInPlane)
        endif

        call Allocator%deallocate(Size%usePinnedMemory, self%label, "hplane1",       HypPlanePtr% hplane1)
        call Allocator%deallocate(Size%usePinnedMemory, self%label, "hplane2",       HypPlanePtr% hplane2)
        call Allocator%deallocate(Size%usePinnedMemory, self%label, "ndone",         HypPlanePtr% ndone)
        call Allocator%deallocate(Size%usePinnedMemory, self%label, "interfaceList", HypPlanePtr% interfaceList)
        call Allocator%deallocate(Size%usePinnedMemory, self%label, "c1",            HypPlanePtr% c1)
        call Allocator%deallocate(Size%usePinnedMemory, self%label, "c2",            HypPlanePtr% c2)

        deallocate( HypPlanePtr% badCornerList )
      endif
    enddo


    return

  end subroutine AngleSet_dtorHypPlane

!-----------------------------------------------------------------------
  function AngleSet_getReflAngle(self,reflID,IncAngle) result(ReflAngle)

!    Return a pointer to a reflected angle given an incident angle 

!    variable declarations
     implicit none

!    passed variables
     type(AngleSet),    intent(in)   :: self
     integer,           intent(in)   :: reflID
     integer,           intent(in)   :: IncAngle
     integer                         :: ReflAngle


     ReflAngle = self% ReflectedAngle(IncAngle,reflID)

     return

  end function AngleSet_getReflAngle

!-----------------------------------------------------------------------
  subroutine AngleSet_setReflAngle(self,reflID,IncAngle,ReflAngle)


!    variable declarations
     implicit none

!    passed variables
     type(AngleSet),    intent(inout)  :: self
     integer,           intent(in)     :: reflID
     integer,           intent(in)     :: IncAngle
     integer,           intent(in)     :: ReflAngle


     self% ReflectedAngle(IncAngle,reflID) = ReflAngle

     return

  end subroutine AngleSet_setReflAngle

!=======================================================================
! construct Boundary Exit List
!=======================================================================
  subroutine AngleSet_ctorBdyExitList(self, angle, nxBdy, bdyList)

    implicit none

!   Passed variables
    type(AngleSet), intent(inout) :: self
    integer,        intent(in)    :: angle
    integer,        intent(in)    :: nxBdy
    integer,        intent(in)    :: bdyList(2,nxBdy)

!   Local
    type(BdyExit),  pointer       :: BdyExitPtr

    integer                       :: b

!   Allocate Memory 

    BdyExitPtr => self% BdyExitPtr(angle)

    BdyExitPtr% nxBdy = nxBdy

    allocate( BdyExitPtr% bdyList(2,nxBdy) )

    do b=1,nxBdy
      BdyExitPtr% bdyList(:,b) = bdyList(:,b)
    enddo


    return

  end subroutine AngleSet_ctorBdyExitList

!=======================================================================
! destruct Boundary Exit List
!=======================================================================
  subroutine AngleSet_dtorBdyExitList(self)

    implicit none

!   Passed variables
    type(AngleSet), intent(inout) :: self

!   Local
    type(BdyExit),  pointer       :: BdyExitPtr

    integer                       :: angle

!   Deallocate Memory

    do angle=1,self% NumAngles
      BdyExitPtr => self% BdyExitPtr(angle)

      deallocate( BdyExitPtr% bdyList )
    enddo


    return

  end subroutine AngleSet_dtorBdyExitList

!=======================================================================
! construct Cycle List
!=======================================================================
  subroutine AngleSet_ctorCycleList(self, cycleList)
    use MemoryAllocator_mod

    implicit none

!   Passed variables
    type(AngleSet), intent(inout) :: self
    integer,        intent(in)    :: cycleList(self% totalCycles)

!   Local

    integer                       :: mCycle 

!   Allocate Memory 
    call Allocator%allocate(Size%usePinnedMemory,self%label,"cycleList", self% cycleList, self% totalCycles+1)

    self% cycleList(:) = 0

    do mCycle=1,self% totalCycles
      self% cycleList(mCycle) = cycleList(mCycle)
    enddo

    return

  end subroutine AngleSet_ctorCycleList

!=======================================================================
! destruct Cycle List
!=======================================================================
  subroutine AngleSet_dtorCycleList(self)
    use MemoryAllocator_mod

    implicit none

!   Passed variables
    type(AngleSet), intent(inout) :: self

!   Deallocate Memory 
    call Allocator%deallocate(Size%usePinnedMemory,self%label,"cycleList",    self% cycleList)


    return

  end subroutine AngleSet_dtorCycleList

!=======================================================================
! construct incident test 
!=======================================================================
  subroutine AngleSet_ctorIncTest(self, sharedID, NumBdyElem, neighborRank)

    implicit none

!   Passed variables
    type(AngleSet),     intent(inout) :: self
    integer,            intent(in)    :: sharedID
    integer,            intent(in)    :: NumBdyElem
    integer,            intent(in)    :: neighborRank

!   Local
    type(IncidentTest), pointer       :: TestPtr

!   Allocate Memory 

    TestPtr => self% TestPtr(sharedID)

    if (self% NumAngles > 1) then
      TETON_VERIFY( mod(self% NumAngles, 2) == 0, "Number of angles in each set must be 1 or even." )

      ! In general, sets have an even number of angles.  Half of them will be
      !   receiving on this boundary, the other half will be sending.
      TestPtr% nSend = NumBdyElem*self% NumAngles/2
      TestPtr% nRecv = TestPtr% nSend
    else
      ! In the case where NumAngles = 1, half the IncTest/IncTestR objects will have
      !  boundary info but the other half won't.
      if (Size%myRankInGroup < neighborRank) then
        ! This should be 0, but 1 seems safer, with a minimal memory burden:
        TestPtr% nSend = 1
        TestPtr% nRecv = NumBdyElem
      else
        TestPtr% nSend = NumBdyElem
        ! This should be 0, but 1 seems safer, with a minimal memory burden:
        TestPtr% nRecv = 1
      endif
    endif

    allocate( TestPtr% IncTest(TestPtr% nSend) )
    allocate( TestPtr% IncTestR(TestPtr% nRecv) )

    TestPtr% IncTest(:)  = 0
    TestPtr% IncTestR(:) = 0

    return

  end subroutine AngleSet_ctorIncTest

!=======================================================================
! get incident test 
!=======================================================================
  function AngleSet_getIncTest(self, sharedID) result(TestPtr)

    implicit none

!   Passed variables
    type(AngleSet),     intent(inout) :: self
    integer,            intent(in)    :: sharedID
    type(IncidentTest), pointer       :: TestPtr


    TestPtr => self% TestPtr(sharedID)

    return

  end function AngleSet_getIncTest

!=======================================================================
! destruct incident test 
!=======================================================================
  subroutine AngleSet_dtorIncTest(self, sharedID)

    use mpi_param_mod
    use mpif90_mod

    implicit none

!   Passed variables
    type(AngleSet),     intent(inout) :: self
    integer,            intent(in)    :: sharedID

    type(IncidentTest), pointer       :: TestPtr

!   Local

    integer  :: ierr

!   Release Memory 

    TestPtr => self% TestPtr(sharedID)

    call MPI_Request_Free(TestPtr% request(1), ierr)
    call MPI_Request_Free(TestPtr% request(2), ierr)

    deallocate( TestPtr% IncTest )
    deallocate( TestPtr% IncTestR )


    return

  end subroutine AngleSet_dtorIncTest


end module AngleSet_mod
