! Zone Set Module:  Contains temporary storage used in zone set loops 

module ZoneSet_mod

  use kind_mod
  use constant_mod
  USE ISO_C_BINDING

  private

! public interfaces

  public construct
  public destruct

  type, public :: ZoneSet

     integer,    pointer, contiguous :: nCornerSet(:)        => null() ! total corners 
     integer,    pointer, contiguous :: nCornerBatch(:)      => null() ! corners in current batch
     integer,    pointer, contiguous :: offset(:)            => null()

     integer,    pointer, contiguous :: cornerList(:)        => null() ! list of corners in batch 
     integer,    pointer, contiguous :: cornerMap(:)         => null()
     integer,    pointer, contiguous :: zoneList(:)          => null()
     integer,    pointer, contiguous :: cornerConverged(:)   => null() ! equals 1 if converged

     real(adqt), pointer, contiguous :: Te(:)                => null()
     real(adqt), pointer, contiguous :: TeOld(:)             => null()
     real(adqt), pointer, contiguous :: delta(:)             => null()
     real(adqt), pointer, contiguous :: sumT(:)              => null()
     real(adqt), pointer, contiguous :: netRate(:)           => null()
     real(adqt), pointer, contiguous :: dTCompton(:)         => null()
     real(adqt), pointer, contiguous :: B(:,:)               => null()
     real(adqt), pointer, contiguous :: dBdT(:,:)            => null()
     real(adqt), pointer, contiguous :: Snu0(:,:)            => null()
     real(adqt), pointer, contiguous :: dSnu0dT(:,:)         => null()
     real(adqt), pointer, contiguous :: AD(:,:)              => null()
     real(adqt), pointer, contiguous :: AU(:,:)              => null()
     real(adqt), pointer, contiguous :: AL(:,:)              => null()
     real(adqt), pointer, contiguous :: z(:,:)               => null()
     real(adqt), pointer, contiguous :: fk2(:,:)             => null()
     real(adqt), pointer, contiguous :: nI(:,:)              => null()
     real(adqt), pointer, contiguous :: nS(:,:)              => null()
     real(adqt), pointer, contiguous :: ex(:,:)              => null()
     real(adqt), pointer, contiguous :: expPH(:,:)           => null()
     real(adqt), pointer, contiguous :: comptonDeltaEr(:,:)  => null()
     real(adqt), pointer, contiguous :: dComptonDT(:,:)      => null()
     real(adqt), pointer, contiguous :: comptonSe(:)         => null()

!    Misc
     character(len=7) :: label

  end type ZoneSet 

  type(ZoneSet), pointer, public :: ZSet => null()

  interface construct
    module procedure ZoneSet_ctor
  end interface

  interface destruct
    module procedure ZoneSet_dtor
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine ZoneSet_ctor(self,        &
                          nZoneSets)

    use Size_mod
    use constant_mod
    use Geometry_mod
    use MemoryAllocator_mod

    implicit none

!   Passed variables

    type(ZoneSet), intent(inout)    :: self

    integer,        intent(in)      :: nZoneSets

!   Local

    integer  :: zSetID
    integer  :: nCornerSet
    integer  :: setSize

    self% label = "zoneSet"

!   Set Properties

    call Allocator%allocate(Size%usePinnedMemory, self%label, "nCornerSet",   self% nCornerSet,   nZoneSets)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "nCornerBatch", self% nCornerBatch, nZoneSets)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "offset",       self% offset,       nZoneSets)

    do zSetID=1,nZoneSets
      nCornerSet = Geom% corner2(zSetID) - Geom% corner1(zSetID) + 1
      self% nCornerSet(zSetID)   = nCornerSet
      self% nCornerBatch(zSetID) = nCornerSet 
      self% offset(zSetID)       = Geom% corner1(zSetID) - 1
    enddo

!   For debugging purpose, assign the "setSize" to all corners. In the near
!   future, we will limit the setSize and run batches

    setSize = Size% ncornr

!   Allocate Memory

    call Allocator%allocate(Size%usePinnedMemory, self%label, "cornerList",      self% cornerList,      setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "cornerMap",       self% cornerMap,       setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "zoneList",        self% zoneList,        setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "cornerConverged", self% cornerConverged, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "Te",              self% Te,        setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "TeOld",           self% TeOld,     setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "delta",           self% delta,     setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "sumT",            self% sumT,      setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "netRate",         self% netRate,   setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "dTCompton",       self% dTCompton, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "B",               self% B,       Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "dBdT",            self% dBdT,    Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "Snu0",            self% Snu0,    Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "dSnu0dT",         self% dSnu0dT, Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "AD",              self% AD,      Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "AU",              self% AU,      Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "AL",              self% AL,      Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "z",               self% z,       4, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "fk2",             self% fk2,     4, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "nI",              self% nI,      Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "nS",              self% nS,      Size% ngr+1, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "ex",              self% ex,      Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "expPH",           self% expPH,   Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "comptonDeltaEr",  self% comptonDeltaEr, Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "dComptonDT",      self% dComptonDT,     Size% ngr, setSize)
    call Allocator%allocate(Size%usePinnedMemory, self%label, "comptonSe",       self% comptonSe, setSize)

    self% cornerList(:)        = 0 
    self% cornerMap(:)         = 0
    self% zoneList(:)          = 0
    self% cornerConverged(:)   = 0 
    self% Te(:)                = zero 
    self% TeOld(:)             = zero
    self% delta(:)             = zero
    self% sumT(:)              = zero
    self% netRate(:)           = zero
    self% dTCompton(:)         = zero
    self% B(:,:)               = zero
    self% dBdT(:,:)            = zero
    self% Snu0(:,:)            = zero
    self% dSnu0dT(:,:)         = zero
    self% AD(:,:)              = zero
    self% AU(:,:)              = zero
    self% AL(:,:)              = zero
    self% z(:,:)               = zero
    self% fk2(:,:)             = zero
    self% nI(:,:)              = zero
    self% nS(:,:)              = zero
    self% ex(:,:)              = zero
    self% expPH(:,:)           = zero
    self% comptonDeltaEr(:,:)  = zero
    self% dComptonDT(:,:)      = zero
    self% comptonSe(:)         = zero

    return

  end subroutine ZoneSet_ctor

!=======================================================================
! destruct interface
!=======================================================================

  subroutine ZoneSet_dtor(self)

    use Size_mod
    use MemoryAllocator_mod

    implicit none

!   Passed variables

    type(ZoneSet), intent(inout)    :: self

    call Allocator%deallocate(Size%usePinnedMemory, self%label, "nCornerSet",      self% nCornerSet)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "nCornerBatch",    self% nCornerBatch)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "offset",          self% offset)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "cornerList",      self% cornerList)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "cornerMap",       self% cornerMap)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "zoneList",        self% zoneList)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "cornerConverged", self% cornerConverged)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "Te",              self% Te)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "TeOld",           self% TeOld)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "delta",           self% delta)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "sumT",            self% sumT)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "netRate",         self% netRate)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "dTCompton",       self% dTCompton)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "B",               self% B)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "dBdT",            self% dBdT)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "Snu0",            self% Snu0)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "dSnu0dT",         self% dSnu0dT)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "AD",              self% AD)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "AU",              self% AU)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "AL",              self% AL)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "z",               self% z)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "fk2",             self% fk2)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "nI",              self% nI)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "nS",              self% nS)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "ex",              self% ex)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "expPH",           self% expPH)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "comptonDeltaEr",  self% comptonDeltaEr)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "dComptonDT",      self% dComptonDT)
    call Allocator%deallocate(Size%usePinnedMemory, self%label, "comptonSe",       self% comptonSe)

    return

  end subroutine ZoneSet_dtor

end module ZoneSet_mod
