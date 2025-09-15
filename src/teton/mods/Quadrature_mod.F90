#include "macros.h"
! Quadrature Module:  Contains data structures describing angular quadrature 

module Quadrature_mod 

  use flags_mod
  use kind_mod
  use constant_mod

  private

! public interfaces

  public construct
  public destruct
  public getNumberOfEnergyGroups
  public getNumberOfAngles
                                                                                 
  type, public :: Quadrature 

     integer              :: QuadID            ! quadrature ID
     integer              :: Groups            ! number of energy groups 
     integer              :: NumAngles         ! number of angles 
     integer              :: nComputeAngles
     integer              :: maxAngleSets
     integer              :: NumMoments        ! number of angular moments
     integer              :: nPolarAngles      ! number of polar angles
     integer              :: Order             ! quadrature order 
     integer              :: NPolar            ! number of polar angles
     integer              :: NAzimuthal        ! number of azimuthal angles
     integer              :: PolarAxis         ! polar axis (1,2 or 3)

     character(len=8)     :: TypeName          ! quadrature type

     logical(kind=1), pointer, contiguous :: StartingDirection(:)
     logical(kind=1), pointer, contiguous :: FinishingDirection(:)

     real(adqt),      pointer, contiguous :: Gnu(:)            ! energy group boundaries
     real(adqt),      pointer, contiguous :: GnuBar(:)         ! average group energy
     real(adqt),      pointer, contiguous :: Omega(:,:)        ! direction cosines 
     real(adqt),      pointer, contiguous :: Weight(:)         ! quadrature weights 
     real(adqt),      pointer, contiguous :: polarAngleList(:)

!    Spherical Harmonics
     real(adqt),      pointer, contiguous :: Ynm(:,:)          ! Ynm(NumMoments,NumAngles)
     real(adqt),      pointer, contiguous :: Pnm(:,:)          ! Pnm(NumMoments,NumAngles)

!    Sweep and communication data structures
     integer                              :: totcycles         ! total number of mesh cycles
     integer                              :: NumBin0          
     integer                              :: NumAnglesDyn
     integer                              :: NumBin            ! number of angle bins for communication

     integer,         pointer, contiguous :: NangBinList(:)    ! NangBinList(NumBin)
     integer,         pointer, contiguous :: AngleToBin(:)     ! AngleToBin(NumAngles) 
     integer,         pointer, contiguous :: PolarAngle(:)     
     integer,         pointer, contiguous :: angleList(:)
     integer,         pointer, contiguous :: angleSetSize(:)

     character(len=8) :: label ! A string descriptor for this quadrature

  end type Quadrature 

  type(Quadrature), pointer, public :: QuadSet => null()

  interface construct
    module procedure Quadrature_ctor
  end interface

  interface destruct
    module procedure Quadrature_dtor
  end interface

  interface getNumberOfEnergyGroups
    module procedure Quadrature_getNumberOfEnergyGroups
  end interface

  interface getNumberOfAngles
    module procedure Quadrature_getNumberOfAngles
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================
                                                                                   
  subroutine Quadrature_ctor(self,       &
                             QuadID,     &
                             Groups,     &
                             NumAngles,  &
                             NumMoments, &
                             Order,      &
                             NPolar,     &
                             NAzimuthal, &
                             PolarAxis,  &
                             TypeName,   &
                             Gnu) 

    use Size_mod
    use constant_mod

    implicit none

!   Passed variables

    type(Quadrature), intent(inout)    :: self

    integer, intent(in)                :: QuadID
    integer, intent(in)                :: Groups       
    integer, intent(in)                :: NumAngles
    integer, intent(in)                :: NumMoments
    integer, intent(in)                :: Order 
    integer, intent(in)                :: NPolar
    integer, intent(in)                :: NAzimuthal
    integer, intent(in)                :: PolarAxis
    character(len=8), intent(in)       :: TypeName
    real(adqt), intent(in)             :: Gnu(Groups+1)

!   Local

    integer, dimension (1) :: imin

    integer          :: i, n, isctp1, Ndim, nlevel, Nang, delNang
    integer          :: bin, NangBin, angle
    integer          :: pAxis, levelNeg, levelPos
    integer          :: g
    integer          :: GroupsP1

    real(adqt)       :: pAngle
    integer          :: igeom

    logical (kind=1), allocatable :: notDone(:)

!   Set Properties

    write(self%label, '(I0.3)') QuadID
    self%label = "quad_"//self%label

    self% QuadID        = QuadID
    self% Groups        = Groups 
    self% NumAngles     = NumAngles 
    self% NumAnglesDyn  = NumAngles
    self% maxAngleSets  = 1

    self% NumMoments    = NumMoments
    self% Order         = Order 
    self% NPolar        = NPolar
    self% NAzimuthal    = NAzimuthal
    self% PolarAxis     = PolarAxis
    self% TypeName      = TypeName

    isctp1              = 1
    Ndim                = Size% ndim
    igeom               = Size% igeom
    GroupsP1            = self% Groups + 1


    allocate( self% Gnu(GroupsP1) )
    allocate( self% gnuBar(Groups) )
    allocate( self% StartingDirection(self% NumAngles) )
    allocate( self% FinishingDirection(self% NumAngles) )
    allocate( self% Omega(Ndim,self% NumAngles) )
    allocate( self% Weight(self% NumAngles) )
    allocate( self% Ynm(self% NumMoments,self% NumAngles) )
    allocate( self% Pnm(self% NumMoments,self% NumAngles) )


    self % Gnu(:) = Gnu(:)

    do g=1,Groups
      self% gnuBar(g) = sqrt( self% gnu(g+1)*self% gnu(g) )
    enddo

    if (self% gnuBar(1) < adqtEpsilon*Size%tfloor) then
      self% gnuBar(1) = self% gnu(2)/sqrt(three)
    endif

!   Space for sweep data structures

    self% totcycles = 0
    delNang         = 0
    pAxis           = 1
    NangBin         = 1

    if (Ndim == 1) then

      self% NumBin         = self% NumAngles
      NangBin              = 1 

      if (igeom == geometry_slab) then
        self% nPolarAngles   = self% NumAngles
        self% nComputeAngles = self% NumAngles
      else
        self% nPolarAngles   = self% NumAngles - 2
        self% nComputeAngles = self% NumAngles - 1
      endif

    elseif (Ndim == 2) then

      if (self% TypeName == 'levelsym') then
        self% NumBin         = self% Order 
        NangBin              = self% Order + 2
        self% nPolarAngles   = self% Order
        self% nComputeAngles = self% NumAngles - self% Order
        delNang              = 2 
      elseif (self% TypeName == 'product') then
        self% NumBin         = 2*self% NPolar
        NangBin              = self% NumAngles/self% NumBin
        self% nPolarAngles   = 2*self% NPolar
        self% nComputeAngles = self% NumAngles - 2*self% NPolar
        delNang              = 0
      else
        TETON_VERIFY(.false., "Invalid quadrature type specified for RZ. Options are level-symmetric and product (1 or 2)")
      endif

    elseif (Ndim == 3) then

      self% nComputeAngles = self% NumAngles
      self% NumBin         = self% NumAngles
      NangBin              = 1

      if (self% TypeName == 'levelsym') then
        pAxis              = 1
        self% nPolarAngles = self% Order
      elseif (self% TypeName == 'product') then
        pAxis              = self% PolarAxis
        self% nPolarAngles = 2*self% NPolar
      elseif (self% TypeName == 'lobatto') then
        pAxis              = self% PolarAxis
        self% nPolarAngles = 2*self% NPolar+2 ! +2 for the two angles on the poles
      else
        TETON_VERIFY(.false., "Invalid quadrature type specified for 3D. Options are level-symmetric, product, and lobatto (1, 2, or 3)");
      endif

    endif

    self% NumBin0 = self% NumBin

    allocate( self% NangBinList(self% NumBin) )
    allocate( self% AngleToBin(self% NumAngles) )
    allocate( self% PolarAngle(self% NumAngles) )
    allocate( self% angleList(self% NumAngles) )
    allocate( self% angleSetSize(self% NumAngles) )
    allocate( self% PolarAngleList(self% nPolarAngles) )

!   Set quadrature points and weights

    call rtquad(self)
    call snynmset(self, Ndim, isctp1)
    call snpnmset(self, Ndim, isctp1)

!   Allow for a variable number of angles per bin

    if (Ndim == 2) then

      nlevel = self% NumBin/2
      Nang   = NangBin

      do i=1,nlevel
        self% NangBinList(2*i-1) = Nang
        self% NangBinList(2*i)   = Nang
        Nang                     = Nang - delNang
      enddo

    else 
      self% NangBinList(:) = NangBin
    endif

!   Create a map from angle to polar angle

    if (Ndim == 1) then

      n = 0
      do angle=1,self% NumAngles
        if (self% weight(angle) < adqtEpsilon) then
          self% PolarAngle(angle) = 0
        else
          n = n + 1
          self% PolarAngle(angle) = n
          self% PolarAngleList(n) = self% omega(1,angle)
        endif
      enddo

    elseif (Ndim == 2) then

      angle    = 0
      levelNeg = self% nPolarAngles/2
      levelPos = levelNeg + 1

      ! This assumes a certain ordering of angles in omega.
      ! Filling this list could be done when we fill in the angles.
      do i=1,self% nPolarAngles,2
        Nang  = self% NangBinList(i)
        angle = angle + Nang
        self% PolarAngleList(levelNeg) = self% omega(2,angle)
        angle = angle + Nang
        self% PolarAngleList(levelPos) = self% omega(2,angle)

        levelNeg = levelNeg - 1
        levelPos = levelPos + 1
      enddo

      do angle=1,self% NumAngles
        if (self% weight(angle) < adqtEpsilon) then
          self% PolarAngle(angle) = 0
        else
          do i=1,self% nPolarAngles
            if (self% omega(2,angle) == self% PolarAngleList(i)) then
              self% PolarAngle(angle) = i
              exit
            endif
          enddo
        endif
      enddo

    elseif (Ndim == 3) then

      allocate( notDone(self% NumAngles) )

      notDone(:) = .TRUE.
      levelNeg   = self% nPolarAngles/2
      levelPos   = self% nPolarAngles

      do i=1,levelNeg

        imin   = minloc( self% omega(pAxis,:) , notDone )
        pAngle = self% omega(pAxis,imin(1))

        self% PolarAngleList(i)        =  pAngle
        self% PolarAngleList(levelPos) = -pAngle

        do angle=1,self% NumAngles

          if (self% omega(pAxis,angle) == pAngle) then
            self% PolarAngle(angle) = i
            notDone(angle)          = .FALSE.
          endif

          if (self% omega(pAxis,angle) == -pAngle) then
            self% PolarAngle(angle) = levelPos 
            notDone(angle)          = .FALSE.
          endif
        enddo

        levelPos = levelPos - 1

      enddo

      deallocate( notDone )

    endif

!   Initialize an angle list used for angle decomposition

    do angle=1,self% NumAngles
      self% angleList(angle)    = angle
      self% angleSetSize(angle) = 0
    enddo

!   Create a map from angle ID to bin ID

    angle = 0
    do bin=1,self% NumBin
      NangBin = self% NangBinList(bin)
      do i=1,NangBin
        self% AngleToBin(angle+i) = bin
      enddo
      angle = angle + NangBin
    enddo
    
  return


  end subroutine Quadrature_ctor

!=======================================================================
! destruct interface
!=======================================================================
                                                                                    
  subroutine Quadrature_dtor(self)

    use Size_mod

    implicit none

!   Passed variables
                                                                                     
    type(Quadrature),  intent(inout) :: self

!   Local

    deallocate( self% Gnu )
    deallocate( self% gnuBar )
    deallocate( self% StartingDirection )
    deallocate( self% FinishingDirection )
    deallocate( self% Omega )
    deallocate( self% Weight )
    deallocate( self% Ynm )
    deallocate( self% Pnm )

!   Space for sweep data structures

    deallocate( self% NangBinList )
    deallocate( self% AngleToBin )
    deallocate( self% PolarAngle )
    deallocate( self% angleList )
    deallocate( self% angleSetSize )
    deallocate( self% PolarAngleList )


    return

  end subroutine Quadrature_dtor

!-----------------------------------------------------------------------
!    Returns the number of groups
!-----------------------------------------------------------------------
  function Quadrature_getNumberOfEnergyGroups(self) result(nGroups)
     type(Quadrature), intent(in) :: self
     integer                      :: nGroups

     nGroups = self% Groups

     return

  end function Quadrature_getNumberOfEnergyGroups

!-----------------------------------------------------------------------
!    Returns the number of energy groups, c callable
!-----------------------------------------------------------------------
  function Teton_Quadrature_getNumberOfEnergyGroups(cptr) bind(c) result(nGroups)

     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_int
     type(c_ptr), value, intent(in) :: cptr
     type(Quadrature), pointer      :: fptr

     integer(kind=c_int)            :: nGroups
     call c_f_pointer(cptr, fptr)
     nGroups = getNumberOfEnergyGroups( fptr )

     return

  end function Teton_Quadrature_getNumberOfEnergyGroups

!-----------------------------------------------------------------------
!    Returns the number of compute angles
!-----------------------------------------------------------------------
  function Quadrature_getNumberOfAngles(self) result(nAngles)
     type(Quadrature), intent(in) :: self
     integer                      :: nAngles

     nAngles = self% NumAngles

     return

  end function Quadrature_getNumberOfAngles

!-----------------------------------------------------------------------
!    Returns the number of compute angles, c callable
!-----------------------------------------------------------------------
  function Teton_Quadrature_getNumberOfAngles(cptr) bind(c) result(nAngles)

     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_int
     type(c_ptr), value, intent(in) :: cptr
     type(Quadrature), pointer      :: fptr

     integer(kind=c_int)            :: nAngles
     call c_f_pointer(cptr, fptr)
     nAngles = getNumberOfAngles( fptr )

     return

  end function Teton_Quadrature_getNumberOfAngles

end module Quadrature_mod

