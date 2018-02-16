! Geometry Modules:  geometry and connectivity information used by Sn
                                                                                 
module Geometry_mod 

  use kind_mod
  use ZoneData_mod
  use cudafor

  private

! public interfaces
                                                                                             
  public construct, destruct, getZoneData
                                                                                 
                                                                                 
  type, public :: Geometry 

     real(adqt), pointer :: px(:,:)           ! px(ndim,npnts)  - coordinates of zone vertices
     real(adqt), pointer :: volc(:)           ! volc(ncornr)    - corner volumes

!    1D & 2D Specific Arrays
     real(adqt), pointer :: areac(:)          ! areac(ncornr)
     real(adqt), pointer :: rj2(:)            ! rj2(nzones)
     real(adqt), pointer :: r2(:)             ! r2(nzones+1)
     real(adqt), pointer :: hj(:)             ! hj(nzones)

     integer, pointer :: CToFace(:,:)         ! CToFace(maxcf,ncornr)
     integer, pointer :: CToZone(:)           ! CToZone(ncornr)
     integer, pointer :: CToPoint(:)          ! CToPoint(ncornr)
     integer, pointer :: ZoneToSrc(:)         ! ZoneToSrc(nzones)
     integer, pointer :: nfpc(:)              ! nfpc(ncornr)

     type(ZoneData), contiguous, pointer :: ZData(:)      ! zone data pointers
     type(ZoneData), device, allocatable :: d_ZData(:)

     type(GPU_ZoneData), pinned, allocatable :: GPU_ZData(:)
     type(GPU_ZoneData), device, allocatable :: d_GPU_ZData(:)

     type(ZoneData_SoA), allocatable :: ZDataSoA ! is this host version really needed?
     type(ZoneData_SoA), device, allocatable :: d_ZDataSoA

     logical :: d_ZData_uptodate

  end type Geometry

  type(Geometry), pointer, public :: Geom

  interface construct
    module procedure Geometry_ctor
  end interface

  interface getZoneData
    module procedure Geometry_getZone
  end interface

  interface destruct
    module procedure Geometry_dtor
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine Geometry_ctor(self)

    use, intrinsic :: iso_c_binding
    use Size_mod

    implicit none

!   Passed variables

    type(Geometry),  intent(inout) :: self

!   Local variables

    integer :: istat

!!$    allocate( self % px(Size%ndim,Size%npnts) )
    allocate( self % volc(Size%ncornr) )

!   Geometry Specific Arrays

    if (Size%ndim == 1) then
      allocate( self % rj2(Size%nzones) )
      allocate( self % r2(Size%nzones+1) )
      allocate( self % hj(Size%nzones) )
      allocate( self % areac(Size%ncornr) )
    endif

!   Mesh Connectivity

    allocate( self % CToFace(Size%maxcf,Size%ncornr) )
    allocate( self % CToZone(Size%ncornr) )
    allocate( self % CToPoint(Size%ncornr) )
    allocate( self % ZoneToSrc(Size%nzones) )
    allocate( self % nfpc(Size%ncornr) )

!   Pointers

    allocate( self % ZData(Size%nzones) )
    allocate( self % d_ZData(Size%nzones) ) 

    allocate( self % GPU_ZData(Size%nzones) )
    allocate( self % d_GPU_ZData(Size%nzones) )
    
    allocate( self % ZDataSoA )
    allocate( self % d_ZDataSoA )

    ! Pin the array pointed to by ZData:
    istat = cudaHostRegister(C_LOC(self%ZData(1)), sizeof(self%ZData), cudaHostRegisterMapped)

    ! is host version of ZDataSoA really needed?
    call constructZones_SoA(self % ZDataSoA)
    istat = cudaMemcpyAsync(C_DEVLOC(self%d_ZDataSoA), C_LOC(self%ZDataSoA), sizeof(self%ZDataSoA), 0)
    ! still not up to date, because ZDataSoA now just points to correct places, but those places have not been
    ! popluated with data yet.
    self % d_ZData_uptodate = .false.

    return
                                                                                             
  end subroutine Geometry_ctor

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
! destruct interface
!=======================================================================

  subroutine Geometry_dtor(self)

    use, intrinsic :: iso_c_binding
    use Size_mod

    implicit none

!   Passed variables

    type(Geometry),  intent(inout) :: self

!   Locals

    integer :: istat


!!$    deallocate( self % px )
    deallocate( self % volc )

!   Geometry Specific Arrays

    if (Size%ndim == 1) then
      deallocate( self % rj2 )
      deallocate( self % r2 )
      deallocate( self % hj )
      deallocate( self % areac )
    endif

!   Mesh Connectivity

    deallocate( self % CToFace )
    deallocate( self % CToZone )
    deallocate( self % CToPoint )
    deallocate( self % ZoneToSrc )
    deallocate( self % nfpc )

!   Pointers

    !! should probably deallocate these here; not doing so
    !! because self%ZData was not deallocated in original code

    !istat = cudaHostUnregister(C_LOC(self % ZData(1)))
    !deallocate( self % ZData )
    !deallocate( self % d_ZData )
    !destructZones_SoA(self%ZDataSoA)
    !deallocate( self % ZDataSoA )
    !deallocate( self % d_ZDataSoA )

    return

  end subroutine Geometry_dtor

                                                      
end module Geometry_mod

