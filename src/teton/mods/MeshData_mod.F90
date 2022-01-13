! MeshData Module:  Contains data structures for the mesh 

module MeshData_mod 

  use kind_mod
  use constant_mod

  private

! public interfaces

  public constructMesh
  public destructMesh 
  public getZoneCenter
  public getFaceCenter

  type, public :: MeshData
     integer                          :: nCorner          ! number of corners
     integer                          :: c0
     integer                          :: nFaces           ! number of faces in zone
     integer                          :: nCFaces

     integer,    pointer, contiguous  :: zoneOpp(:)       ! zone ID across this face
     integer,    pointer, contiguous  :: faceOpp(:)
     integer,    pointer, contiguous  :: CToFace(:,:)     ! CToFace(maxcf,ncornr)
     integer,    pointer, contiguous  :: nCFacesArray(:)

     real(adqt), pointer, contiguous  :: px(:,:)          ! coordinates of zone vertices

  end type MeshData
                                                                                 
  type(MeshData),  pointer, public :: M => null() 


  interface constructMesh
    module procedure MeshData_ctor,  &
                     MeshData1D_ctor
  end interface

  interface destructMesh
    module procedure MeshData_dtor
  end interface

  interface getZoneCenter
    module procedure MeshData_getZoneCenter
  end interface

  interface getFaceCenter
    module procedure MeshData_getFaceCenter
  end interface

contains

!=======================================================================
! construct Mesh interface
!=======================================================================
  subroutine MeshData_ctor(self,     &
                           nCorner,  &
                           c0,       &
                           nFaces,   &
                           zoneOpp,  &
                           CToFace,  &
                           nCFacesArray)

    use Size_mod

    implicit none
    
!   Passed variables 
    type(MeshData),    intent(inout) :: self
    integer,           intent(in)    :: nCorner
    integer,           intent(in)    :: c0
    integer,           intent(in)    :: nFaces
    integer,           intent(in)    :: zoneOpp(nFaces)
    integer,           intent(in)    :: CToFace(Size%maxcf,Size%maxCorner)
    integer, intent(in)              :: nCFacesArray(Size%maxCorner)

!   Local

    integer          :: cID, nCFaces1

    self% nCorner    = nCorner
    self% c0         = c0
    self% nFaces     = nFaces
    self% nCFaces    = 2

    allocate( self% CToFace(Size%maxcf,Size%maxCorner) )
    allocate( self% zoneOpp(self% nFaces) )
    allocate( self% faceOpp(self% nFaces) )
    allocate( self% px(Size% ndim,self% nCorner) )

    self% CToFace(:,:) = CToFace(:,:)
    self% faceOpp(:)   = 0
    self% zoneOpp(:)   = zoneOpp(:)

    if (Size% ndim == 3) then
      allocate( self% nCFacesArray(self% nCorner) )

!     Allow for a variable number of corner-faces per corner

      do cID=1,self% nCorner
        self% nCFacesArray(cID) = nCFacesArray(cID)
      enddo

      nCFaces1      = self% nCFacesArray(1)
      self% nCFaces = self% nCFacesArray(1)

      do cID=2,nCorner
        if (self% nCFacesArray(cID) /= nCFaces1) then
          self% nCFaces = 0
        endif
      enddo

    endif

    return 

  end subroutine MeshData_ctor

!=======================================================================
! construct 1D Mesh interface
!=======================================================================
  subroutine MeshData1D_ctor(self, c0)

    use Size_mod

    implicit none

!   Passed variables 
    type(MeshData),    intent(inout) :: self
    integer,           intent(in)    :: c0

!   Local variables
    integer                          :: cLocal

    self% nCorner = 2 
    self% nFaces  = 2 
    self% nCFaces = 1
    self% c0      = c0

    ! Assuming only 2 zone faces per zone,
    ! indexed consistently with local corners
    allocate( self% CToFace(Size% ndim,self% nCorner) )
    do cLocal=1,self%nCorner
      self% CToFace(1,cLocal) = cLocal
    enddo

    allocate( self% px(Size% ndim,self% nCorner) )

    return

  end subroutine MeshData1D_ctor

!=======================================================================
! getZoneCenter interface
!=======================================================================

  function MeshData_getZoneCenter(self) result(zoneCenter)

!    Return the zone center

     use Size_mod

     implicit none

!    passed variables
     type(MeshData), intent(in) :: self
     real(adqt)                 :: zoneCenter(Size% ndim) 

!    Local
     integer  :: c

!    Accumulate sum of coordinates:

     zoneCenter(:) = zero

     do c=1,self% nCorner
       zoneCenter(:) = zoneCenter(:) + self% px(:,c)
     enddo

!    Divide by number of corners to get average coordinate:

     zoneCenter(:) = zoneCenter(:)/real(self% nCorner,adqt)


     return

  end function MeshData_getZoneCenter

!=======================================================================
! getFaceCenter interface
!=======================================================================

  function MeshData_getFaceCenter(self) result(faceCenter)

!    Return the zone center

     use Size_mod

     implicit none

!    passed variables
     type(MeshData), intent(in) :: self
     real(adqt)                 :: faceCenter(3,self% nFaces)

!    Local
     integer  :: c
     integer  :: cface
     integer  :: face
     integer  :: nCFaces
     integer  :: nc_face(self% nFaces)

!    Sum all point coordinates associated with a face and
!    the number of corner faces on a zone face

     nc_face(:)      = 0
     faceCenter(:,:) = zero

     do c=1,self% nCorner
       nCFaces = self% nCFacesArray(c)
       do cface=1,nCFaces
         face               = self% CToFace(cface,c)
         nc_face(face)      = nc_face(face) + 1
         faceCenter(:,face) = faceCenter(:,face) + self% px(:,c)
       enddo
     enddo

!    Calculate the face-center coordinates as the average of point
!    coordinates associated with that face

     do face=1,self% nFaces
       if (nc_face(face) /= 0) then
         faceCenter(:,face) = faceCenter(:,face)/real(nc_face(face),adqt)
       else
         faceCenter(:,face) = zero
       endif
     enddo

     return

  end function MeshData_getFaceCenter

!=======================================================================
! destruct Mesh interface
!=======================================================================
                                                                                    
  subroutine MeshData_dtor(self)

    use Size_mod

    implicit none

!   Passed variables
    type(MeshData),  intent(inout)   :: self


    if (Size% ndim > 1) then
      deallocate( self% CToFace )
      deallocate( self% zoneOpp )
      deallocate( self% faceOpp )
    endif

    if (Size% ndim == 3) then
      deallocate( self% nCFacesArray )
    endif

    deallocate( self% px )

    return

  end subroutine MeshData_dtor


end module MeshData_mod

