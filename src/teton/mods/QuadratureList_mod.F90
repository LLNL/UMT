#include "macros.h"

! QuadratureList Module:  Contains structures for quadrature and phase-space
! sets

module QuadratureList_mod
  use kind_mod
  use Quadrature_mod
  use SetData_mod
  use AngleSet_mod
  use GroupSet_mod
  use CommSet_mod
  use Size_mod
  use CodeChecks_mod
  implicit none

  private

! public interfaces

  public construct 
  public destruct  
  public setQuadrature
  public getQuadrature 
  public getNumberOfSets
  public getNumberOfGTASets
  public getNumberOfAngleSets
  public getNumberOfGroupSets
  public getNumberOfCommSets
  public getNumberOfZoneSets
  public getNumberOfHyperDomains
  public getNumberOfHyperElements
  public constructSetPointers
  public getGTAQuadrature 
  public getSNQuadrature 
  public getEnergyGroups
  public getGroupAverageEnergy
  public setSweepCounters
  public getTotalSweeps
  public getSetData
  public getGTASetData
  public getAngleSetData
  public getGroupSetData
  public getCommSetData
  public getAngleSetFromSetID
  public getGroupSetFromSetID
  public getCommSetFromSetID
  public getNumberOfGroups
  public getNumberOfAngles
  public getSetIDfromGroupAngle

  type, public :: QuadratureList

     integer                         :: nGTASets ! Number of computational sets for GTA sweeps
     integer                         :: nSets    ! Number of computational sets for Sn sweeps

     ! nSets is the maximum number of sets we can break the Sn quadrature sets
     !   into, as determined by the user input nSetsMaster and also hardware
     !   constraints.
     ! WARNING: nSets may NOT actually be the # of sets used.  It is generally
     !   an upper bound on the number of sets Teton will use in its GPU sweeps.
     !   At the point that the value of nSets is determined, we have not yet
     !   begun decomposing the SN problem into sets.

     integer                         :: nHyperDomains(2)
     integer                         :: nHyperElements(2)
     integer                         :: nAngleSets
     integer                         :: nGroupSets
     integer                         :: nCommSets
     integer                         :: nZoneSets
     integer, pointer, contiguous    :: groupID(:)
     integer, pointer, contiguous    :: angleID(:)
     integer, pointer, contiguous    :: commID(:)
     integer, pointer, contiguous    :: SetIDList(:,:)

     type(Quadrature),   pointer     :: QuadPtr(:)    => null() ! Pointers to quadrature sets
     type(SetData),      pointer     :: SetDataPtr(:) => null()
     type(AngleSet),     pointer     :: AngSetPtr(:)  => null()
     type(GroupSet),     pointer     :: GrpSetPtr(:)  => null()
     type(CommSet),      pointer     :: CommSetPtr(:) => null()

  end type QuadratureList

  type(QuadratureList), pointer, public :: Quad => null()

  interface construct
    module procedure QuadratureList_ctor
  end interface

  interface constructSetPointers
    module procedure QuadratureList_ctorSetPointers
  end interface

  interface setQuadrature 
    module procedure QuadratureList_set
  end interface

  interface getQuadrature
    module procedure QuadratureList_getQuad
  end interface

  interface getNumberOfSets
    module procedure QuadratureList_getNumberOfSets
  end interface

  interface getNumberOfGTASets
    module procedure QuadratureList_getNumberOfGTASets
  end interface

  interface getNumberOfAngleSets
    module procedure QuadratureList_getNumberOfAngleSets
  end interface

  interface getNumberOfGroupSets
    module procedure QuadratureList_getNumberOfGroupSets
  end interface

  interface getNumberOfCommSets
    module procedure QuadratureList_getNumberOfCommSets
  end interface

  interface getNumberOfZoneSets
    module procedure QuadratureList_getNumberOfZoneSets
  end interface

  interface getNumberOfHyperDomains
    module procedure QuadratureList_getNumberOfHyperDomains
  end interface

  interface getNumberOfHyperElements
    module procedure QuadratureList_getNumberOfHyperElements
  end interface

  interface getGTAQuadrature
    module procedure QuadratureList_getGTAQuad
  end interface

  interface getSNQuadrature
    module procedure QuadratureList_getSNQuad
  end interface

  interface getEnergyGroups
    module procedure QuadratureList_getEnergyGroups
  end interface

  interface getGroupAverageEnergy
    module procedure QuadratureList_getGroupAverageEnergy
  end interface

  interface getNumberOfGroups
    module procedure QuadratureList_getNumberOfGroups
  end interface

  interface getNumberOfAngles
    module procedure QuadratureList_getNumberOfAngles
  end interface

  interface setSweepCounters
    module procedure QuadratureList_setCounters
  end interface

  interface getTotalSweeps
    module procedure QuadratureList_getSweeps
  end interface

  interface getSetData
    module procedure QuadratureList_getSetData
  end interface

  interface getAngleSetData
    module procedure QuadratureList_getAngleSetData
  end interface

  interface getGroupSetData
    module procedure QuadratureList_getGroupSetData
  end interface

  interface getCommSetData
    module procedure QuadratureList_getCommSetData
  end interface

  interface getAngleSetFromSetID
    module procedure QuadratureList_getAngleSetFromSetID
  end interface

  interface getGroupSetFromSetID
    module procedure QuadratureList_getGroupSetFromSetID
  end interface

  interface getCommSetFromSetID
    module procedure QuadratureList_getCommSetFromSetID
  end interface

  interface getGTASetData
    module procedure QuadratureList_getGTASetData
  end interface

  interface getSetIDfromGroupAngle
    module procedure QuadratureList_getSetIDfromGroupAngle
  end interface

  interface destruct
    module procedure QuadratureList_dtor
  end interface

contains

!=======================================================================
! construct interface
!=======================================================================

  subroutine QuadratureList_ctor(self, nAnglesSn, nSetsMaster, nSets)

    use cmake_defines_mod, only : omp_device_team_thread_limit
    use Options_mod, only : Options
    use constant_mod

!   Passed variables

    type(QuadratureList), intent(inout) :: self
    integer,              intent(in)    :: nAnglesSn 
    integer,              intent(in)    :: nSetsMaster
    integer,              intent(inout) :: nSets
!   Warning: nSets may not be the value you think it should be.  See comments
!     above in the QuadratureList type definition

!   Local

#if defined(TETON_ENABLE_OPENMP)
    integer :: nOmpMaxThreads
    integer :: nOmpMaxTeams
#endif


    self% nSets             = 1
    self% nGroupSets        = 1
    self% nZoneSets         = 1
    self% nHyperElements(:) = 0 

#if defined(TETON_ENABLE_OPENMP)
    if (Size%useGPU) then
       ! If running on gpu, set these to use all available gpu processors.
       nOmpMaxTeams    = Options%getNumDeviceProcessors()
       self% nSets     = max(1, nOmpMaxTeams)
       self% nZoneSets = min(Size%nzones, nOmpMaxTeams)
    else
       ! If not running on gpu, set these to use all available CPU threads.
       nOmpMaxThreads  = Options%getNumOmpMaxThreads()
       self% nSets     = max(1, nOmpMaxThreads)
       self% nZoneSets = min(Size%nzones, nOmpMaxThreads)
    endif
#endif

!   For debugging, it may be useful to assign the number of sets
    if (nSetsMaster <= 0) then
      nSets = self% nSets
    else
      self% nSets     = nSetsMaster
      nSets           = nSetsMaster
      self% nZoneSets = min(nSetsMaster, Size%nzones)
    endif

    allocate( self% SetIDList(Size% ngr,nAnglesSn) )
    allocate( self% QuadPtr(2) )

    return

  end subroutine QuadratureList_ctor

!=======================================================================
! construct set pointer interface
!=======================================================================

  subroutine QuadratureList_ctorSetPointers(self,       &
                                            nSetsNew,   &
                                            nAngleSets, &
                                            nGroupSets, &
                                            nCommSets,  &
                                            nGTASets)

!   Passed variables

    type(QuadratureList), intent(inout) :: self
    integer,              intent(in)    :: nSetsNew
    integer,              intent(in)    :: nAngleSets
    integer,              intent(in)    :: nGroupSets
    integer,              intent(in)    :: nCommSets
    integer,              intent(in)    :: nGTASets


    self% nSets      = nSetsNew
    self% nAngleSets = nAngleSets
    self% nGroupSets = nGroupSets
    self% nCommSets  = nCommSets
    self% nGTASets   = nGTASets

    if ( .not. associated( self% SetDataPtr ) ) then
      allocate( self% SetDataPtr(self% nSets+self% nGTASets) )
      allocate( self% AngSetPtr(self% nAngleSets+self% nGTASets) )
      allocate( self% CommSetPtr(self% nCommSets+self% nGTASets) )
      allocate( self% GrpSetPtr(self% nGroupSets) )

      allocate( self% groupID(self% nSets) )
      allocate( self% angleID(self% nSets+self% nGTASets) )
      allocate( self% commID(self% nSets+self% nGTASets) )
    endif


    return

  end subroutine QuadratureList_ctorSetPointers

!=======================================================================
! get C pointer to QuadratureList, c callable
!=======================================================================
  function Teton_QuadratureList_getQuadList() bind(c) result(QuadList)

     use, intrinsic :: iso_c_binding, only : c_ptr, c_loc
   
     type(c_ptr) :: QuadList
     QuadList = c_loc(Quad)

     return

  end function Teton_QuadratureList_getQuadList

!=======================================================================
! set interface
!=======================================================================

  subroutine QuadratureList_set(self,          &
                                quadID,        &
                                Groups,        &
                                NumAngles,     &
                                NumMoments,    &
                                Order,         &
                                NPolar,        &
                                NAzimuthal,    &
                                PolarAxis,     &
                                QuadType,      &
                                Gnu)

!   Passed variables

    type(QuadratureList), intent(inout) :: self

    integer,    intent(in)              :: quadID
    integer,    intent(in)              :: Groups
    integer,    intent(in)              :: NumAngles
    integer,    intent(in)              :: NumMoments
    integer,    intent(in)              :: Order
    integer,    intent(in)              :: NPolar
    integer,    intent(in)              :: NAzimuthal
    integer,    intent(in)              :: PolarAxis
    integer,    intent(in)              :: QuadType
    real(adqt), intent(in)              :: Gnu(Groups+1)

!   Local

    character(len=8) :: TypeName

    select case (QuadType)
      case (1)
        TypeName = 'levelsym'
      case (2)
        TypeName = 'product'
      case (3)
        TypeName = 'lobatto'
    end select 
 
    call construct(self% QuadPtr(quadID), &
                         quadID,          &
                         Groups,          &
                         NumAngles ,      &
                         NumMoments,      &
                         Order,           &
                         NPolar,          &
                         NAzimuthal,      &
                         PolarAxis,       &
                         TypeName,        &
                         Gnu)

    return

  end subroutine QuadratureList_set

!=======================================================================
! getSetIDfromGroupAngle interface
!=======================================================================

  function QuadratureList_getSetIDfromGroupAngle(self, group, angle) result(setID)

!   Passed variables

    type(QuadratureList),  intent(inout) :: self
    integer,               intent(in)    :: group
    integer,               intent(in)    :: angle

    integer                              :: setID

    TETON_CHECK_BOUNDS2(self%SetIDList, group, angle)
    setID = self% SetIDList(group,angle)

    return

  end function QuadratureList_getSetIDfromGroupAngle

!=======================================================================
! destruct interface
!=======================================================================

  subroutine QuadratureList_dtor(self)

!   Passed variables

    type(QuadratureList),  intent(inout) :: self

!   Local


    deallocate( self% SetIDList )
    deallocate( self% QuadPtr )
    deallocate( self% SetDataPtr )
    deallocate( self% GrpSetPtr )
    deallocate( self% AngSetPtr )
    deallocate( self% CommSetPtr )

    deallocate( self% angleID )
    deallocate( self% groupID )
    deallocate( self% commID )

    return

  end subroutine QuadratureList_dtor

!-----------------------------------------------------------------------
!    Returns the number of phase-space sets (nSets)
!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfSets(self) result(nSets)

     type(QuadratureList), intent(in) :: self
     integer                          :: nSets

     nSets = self% nSets

     return

  end function QuadratureList_getNumberOfSets

!-----------------------------------------------------------------------
!    Returns the number of phase-space sets (nSets), c callable
!-----------------------------------------------------------------------
  function Teton_QuadratureList_getNumberOfSets(cptr) bind(c) result(nSets)

     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_int
     type(c_ptr), value, intent(in) :: cptr
     type(QuadratureList), pointer  :: fptr

     integer(kind=c_int)            :: nSets
     call c_f_pointer(cptr, fptr)
     nSets = getNumberOfSets( fptr )

     return

  end function Teton_QuadratureList_getNumberOfSets

!-----------------------------------------------------------------------
!    Returns the number of GTA phase-space sets, c callable
!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfGTASets(self) result(nGTASets)
!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nGTASets

     nGTASets = self% nGTASets

     return

  end function QuadratureList_getNumberOfGTASets

!-----------------------------------------------------------------------
!    Returns the number of GTA phase-space sets, c callable
!-----------------------------------------------------------------------
  function Teton_QuadratureList_getNumberOfGTASets(cptr) bind(c) result(nGTASets)

     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_int
     type(c_ptr), value, intent(in)   :: cptr
     type(QuadratureList), pointer    :: fptr

     integer(kind=c_int)              :: nGTASets
     call c_f_pointer(cptr, fptr)
     nGTASets = getNumberOfGTASets( fptr )

     return

  end function Teton_QuadratureList_getNumberOfGTASets

!-----------------------------------------------------------------------
!    Returns the number of angle sets
!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfAngleSets(self) result(nAngleSets)
!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nAngleSets

     nAngleSets = self% nAngleSets

     return

  end function QuadratureList_getNumberOfAngleSets

!-----------------------------------------------------------------------
!    Returns the number of angle sets, c callable
!-----------------------------------------------------------------------
  function Teton_QuadratureList_getNumberOfAngleSets(cptr) bind(c) result(nAngleSets)

     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_int
     type(c_ptr), value, intent(in) :: cptr
     type(QuadratureList), pointer  :: fptr

     integer(kind=c_int)            :: nAngleSets
     call c_f_pointer(cptr, fptr)
     nAngleSets = getNumberOfAngleSets( fptr )

     return

  end function Teton_QuadratureList_getNumberOfAngleSets


!-----------------------------------------------------------------------
!    Returns the number of energy group sets
!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfGroupSets(self) result(nGroupSets)
!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nGroupSets

     nGroupSets = self% nGroupSets

     return

  end function QuadratureList_getNumberOfGroupSets

!-----------------------------------------------------------------------
!    Returns the number of energy group sets, c callable
!-----------------------------------------------------------------------
  function Teton_QuadratureList_getNumberOfGroupSets(cptr) bind(c) result(nGroupSets)

     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_int
     type(c_ptr), value, intent(in) :: cptr
     type(QuadratureList), pointer  :: fptr

     integer(kind=c_int)            :: nGroupSets
     call c_f_pointer(cptr, fptr)
     nGroupSets = getNumberOfGroupSets( fptr )

     return

  end function Teton_QuadratureList_getNumberOfGroupSets



!-----------------------------------------------------------------------
!    Returns the number of communication sets
!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfCommSets(self) result(nCommSets)

!    Returns the number of communication sets
!      nCommSets   number of communication sets

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nCommSets

     nCommSets = self% nCommSets

     return

  end function QuadratureList_getNumberOfCommSets

!-----------------------------------------------------------------------
!    Returns the number of communication sets, c callable
!-----------------------------------------------------------------------
  function Teton_QuadratureList_getNumberOfCommSets(cptr) bind(c) result(nCommSets)

     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_int
     type(c_ptr), value, intent(in) :: cptr
     type(QuadratureList), pointer  :: fptr

     integer(kind=c_int)            :: nCommSets
     call c_f_pointer(cptr, fptr)
     nCommSets = getNumberOfCommSets( fptr )

     return

  end function Teton_QuadratureList_getNumberOfCommSets

!-----------------------------------------------------------------------
!    Returns the number of zone sets
!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfZoneSets(self) result(nZoneSets)

!    Returns the number of zone sets
!      nZoneSets   number of zone sets

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nZoneSets

     nZoneSets = self% nZoneSets

     return

  end function QuadratureList_getNumberOfZoneSets

!-----------------------------------------------------------------------
!    Returns the number of zone sets, c callable
!-----------------------------------------------------------------------
  function Teton_QuadratureList_getNumberOfZoneSets(cptr) bind(c) result(nZoneSets)

     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_int
     type(c_ptr), value, intent(in) :: cptr
     type(QuadratureList), pointer  :: fptr

     integer(kind=c_int)            :: nZoneSets
     call c_f_pointer(cptr, fptr)
     nZoneSets = getNumberOfZoneSets( fptr )

     return

  end function Teton_QuadratureList_getNumberOfZoneSets

!-----------------------------------------------------------------------
!    Returns the number of hyper domains used for sweeps, c callable
!      nHyperDomains   number of hyper domains 
!      ID = 1   High-order angle set
!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfHyperDomains(self, ID) result(nHyperDomains)
!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: ID
     integer                          :: nHyperDomains 

     TETON_ASSERT(ID <= 2, "Invalid ID for getNumberOfHyperDomains, must be 1 or 2")
     nHyperDomains = self% nHyperDomains(ID) 

     return

  end function QuadratureList_getNumberOfHyperDomains

!-----------------------------------------------------------------------
!    Returns the number of hyper domains used for sweeps, c callable
!      nHyperDomains   number of hyper domains 
!      ID = 1   High-order angle set
!      ID = 2   GTA angle set
!-----------------------------------------------------------------------
  function Teton_QuadratureList_getNumberOfHyperDomains(cptr, ID) bind(c) result(nHyperDomains)

     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_int
     type(c_ptr), value, intent(in) :: cptr
     integer(kind=c_int), value     :: ID
     type(QuadratureList), pointer  :: fptr

     integer(kind=c_int)            :: nHyperDomains
     call c_f_pointer(cptr, fptr)
     nHyperDomains = getNumberOfHyperDomains( fptr, ID )

     return

  end function Teton_QuadratureList_getNumberOfHyperDomains

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfHyperElements(self, ID) result(nHyperElements)

!    Returns the number of interface elements at hyper-domain
!    boundaries  used for sweeps 
!      nHyperElements   number of interface elements 
!      ID = 1   High-order angle set
!      ID = 2   GTA angle set

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: ID
     integer                          :: nHyperElements 

     TETON_ASSERT(ID <= 2, "Invalid ID for getNumberOfHyperElements, must be 1 or 2")
     nHyperElements = self% nHyperElements(ID)

     return

  end function QuadratureList_getNumberOfHyperElements

!-----------------------------------------------------------------------
  function QuadratureList_getQuad(self,QuadID) result(QuadPtr)

!    Return a pointer to a quadrature set 
!      QuadID   quadrature set ID number 
!      QuadPtr  pointer to the quadrature set 

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: QuadID 
     type(Quadrature),     pointer    :: QuadPtr

            
     TETON_CHECK_BOUNDS1(self%QuadPtr, QuadID)
     QuadPtr => self% QuadPtr(QuadID)
            
     return

  end function QuadratureList_getQuad 

!-----------------------------------------------------------------------
!    Return a C pointer to a quadrature set, C callable.
!-----------------------------------------------------------------------
  function Teton_QuadratureList_getQuad(quadlist_cptr, quadID) bind(c) result(quad_cptr)
     use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, c_loc, c_int
     type(c_ptr), value, intent(in)         :: quadlist_cptr
     integer(kind=c_int), value, intent(in) :: quadID
     type(c_ptr)                            :: quad_cptr

     type(QuadratureList), pointer          :: quadlist_fptr
     type(Quadrature), pointer              :: quad_fptr

     call c_f_pointer(quadlist_cptr, quadlist_fptr)
     quad_fptr => getQuadrature(quadlist_fptr, quadID)
     quad_cptr = c_loc(quad_fptr)
     
     return
  end function Teton_QuadratureList_getQuad

!-----------------------------------------------------------------------
!    Return a pointer to the GTA quadrature set
!-----------------------------------------------------------------------
  function QuadratureList_getGTAQuad(self) result(QuadPtr)
!    variable declarations
                                                                                           
!    passed variables
     type(QuadratureList), intent(in) :: self
     type(Quadrature),     pointer    :: QuadPtr
                                                                                           
                                                                                           
     QuadPtr => self% QuadPtr(2)
                                                                                           
     return
                                                                                           
  end function QuadratureList_getGTAQuad

!-----------------------------------------------------------------------
!    Return a pointer to the SN quadrature set
!-----------------------------------------------------------------------
function QuadratureList_getSNQuad(self) result(QuadPtr)
                                                                                           
!      QuadPtr  pointer to the quadrature set
                                                                                           
!    variable declarations
                                                                                           
!    passed variables
     type(QuadratureList), intent(in) :: self
     type(Quadrature),     pointer    :: QuadPtr
                                                                                           
                                                                                           
     QuadPtr => self% QuadPtr(1)

     return
                                                                                           
  end function QuadratureList_getSNQuad

!-----------------------------------------------------------------------
  function QuadratureList_getSetData(self,setID) result(SetDataPtr)

!    Return a pointer to an energy group/angle set 
!      setID       set ID number 
!      SetDataPtr  pointer to the set 

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     type(SetData),        pointer    :: SetDataPtr 


     TETON_CHECK_BOUNDS1(self%SetDataPtr, setID)
     SetDataPtr => self% SetDataPtr(setID)

     return

  end function QuadratureList_getSetData

!-----------------------------------------------------------------------
!    Return a pointer to an angle set
!      angleSetID  angle set ID number
!-----------------------------------------------------------------------
  function QuadratureList_getAngleSetData(self,angleSetID) result(AngSetPtr)

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: angleSetID
     type(AngleSet),       pointer    :: AngSetPtr


     TETON_CHECK_BOUNDS1(self%AngSetPtr, angleSetID)
     AngSetPtr => self% AngSetPtr(angleSetID)

     return

  end function QuadratureList_getAngleSetData

!-----------------------------------------------------------------------
!    Return a pointer to an energy group set 
!      groupSetID  group set ID number    
!      GrpSetPtr   pointer to the group set 
!-----------------------------------------------------------------------
  function QuadratureList_getGroupSetData(self,groupSetID) result(GrpSetPtr)

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in)       :: self
     integer(kind=c_int), value, intent(in) :: groupSetID
     type(GroupSet), pointer                :: GrpSetPtr


     TETON_CHECK_BOUNDS1(self%GrpSetPtr, groupSetID)
     GrpSetPtr => self% GrpSetPtr(groupSetID)

     return

  end function QuadratureList_getGroupSetData

!-----------------------------------------------------------------------
  function QuadratureList_getCommSetData(self,commSetID) result(CommSetPtr)

!    Return a pointer to a communication set 
!      commSetID    comm set ID number   
!      CommSetPtr   pointer to the communication set 

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: commSetID
     type(CommSet),        pointer    :: CommSetPtr


     TETON_CHECK_BOUNDS1(self%CommSetPtr, commSetID)
     CommSetPtr => self% CommSetPtr(commSetID)

     return

  end function QuadratureList_getCommSetData

!-----------------------------------------------------------------------
  function QuadratureList_getAngleSetFromSetID(self,setID) result(AngSetPtr)

!    Return a pointer to an angle set 
!      setID       set ID number 
!      AngSetPtr   pointer to the angle set 

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     type(AngleSet),       pointer    :: AngSetPtr

     integer                          :: localID

     TETON_CHECK_BOUNDS1(self%angleID, setID)
     localID = self% angleID(setID)


     TETON_CHECK_BOUNDS1(self%AngSetPtr, localID)
     AngSetPtr => self% AngSetPtr(localID)

     return

  end function QuadratureList_getAngleSetFromSetID

!-----------------------------------------------------------------------
  function QuadratureList_getGroupSetFromSetID(self,setID) result(GrpSetPtr)

!    Return a pointer to an energy group set 
!      setID       set ID number 
!      GrpSetPtr   pointer to the group set 

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     type(GroupSet),       pointer    :: GrpSetPtr

     integer                          :: localID

     TETON_CHECK_BOUNDS1(self%groupID, setID)
     localID = self% groupID(setID)


     TETON_CHECK_BOUNDS1(self%GrpSetPtr, localID)
     GrpSetPtr => self% GrpSetPtr(localID)

     return

  end function QuadratureList_getGroupSetFromSetID

!-----------------------------------------------------------------------
  function QuadratureList_getCommSetFromSetID(self,setID) result(CommSetPtr)

!    Return a pointer to a communication set 
!      setID        set ID number 
!      CommSetPtr   pointer to the communication set 

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     type(CommSet),        pointer    :: CommSetPtr

     integer                          :: localID

     TETON_CHECK_BOUNDS1(self%commID, setID)
     localID = self% commID(setID)


     TETON_CHECK_BOUNDS1(self%CommSetPtr, localID)
     CommSetPtr => self% CommSetPtr(localID)

     return

  end function QuadratureList_getCommSetFromSetID

!-----------------------------------------------------------------------
  function QuadratureList_getGTASetData(self,setID) result(SetDataPtr)

!    Return a pointer to an energy group/angle set 
!      setID       set ID number 
!      SetDataPtr  pointer to the set 

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     type(SetData),        pointer    :: SetDataPtr


     TETON_CHECK_BOUNDS1(self%SetDataPtr, self% nSets + setID)
     SetDataPtr => self% SetDataPtr(self% nSets + setID)

     return

  end function QuadratureList_getGTASetData

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfGroups(self,setID) result(nGroups)

!    Return the number of groups in this group/angle set 
!      setID       set ID number 
!      Groups      number of energy groups 

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     integer                          :: nGroups
     type(SetData),        pointer    :: SetDataPtr


     TETON_CHECK_BOUNDS1(self%SetDataPtr, setID)
     SetDataPtr => self% SetDataPtr(setID)
     nGroups     =  SetDataPtr% Groups

     return

  end function QuadratureList_getNumberOfGroups

!-----------------------------------------------------------------------
!    Returns the number of angles in a phase space set's angle set.
!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfAngles(self,setID) result(NumAngles)

!    Return the number of angles in this group/angle set 
!      setID       set ID number 
!      NumAngles   number of angles 

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     integer                          :: NumAngles 
     type(SetData),        pointer    :: SetDataPtr


     TETON_CHECK_BOUNDS1(self%SetDataPtr, setID)
     SetDataPtr => self% SetDataPtr(setID)
     NumAngles  =  SetDataPtr% NumAngles 

     return

  end function QuadratureList_getNumberOfAngles

!-----------------------------------------------------------------------
  function QuadratureList_getEnergyGroups(self,numGroups) result(GrpBnds)
                                                                                            
!    Returns all energy group bounds
!      GrpBnds    array of energy group bounds 
                                                                                            
!    variable declarations
                                                                                            
!    passed variables
     type(QuadratureList), intent(in) :: self
     integer, intent(in)              :: numGroups

!    Local
     type(Quadrature), pointer    :: QuadPtr
     integer                      :: g
     real(adqt)                   :: GrpBnds(numGroups+1)
                                                                                            
     QuadPtr => self% QuadPtr(1)
     TETON_CHECK_BOUNDS1(QuadPtr%Gnu, numGroups+1)
     do g=1,numGroups+1
       GrpBnds(g) = QuadPtr% Gnu(g)
     enddo

     return
                                                                                            
  end function QuadratureList_getEnergyGroups

!-----------------------------------------------------------------------
  function QuadratureList_getGroupAverageEnergy(self,numGroups) result(gnuBar)

!    Returns all group-average energies
!      gnuBar    array of group-average energy 

     use constant_mod

!    variable declarations

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer, intent(in)              :: numGroups

!    Local
     type(Quadrature), pointer    :: QuadPtr
     integer                      :: g
     real(adqt)                   :: gnuBar(0:numGroups+1)

     QuadPtr => self% QuadPtr(1)
     TETON_CHECK_BOUNDS1(QuadPtr%gnuBar, numGroups)
     do g=1,numGroups
       gnuBar(g) = QuadPtr% gnuBar(g)
     enddo

     gnuBar(0)           = zero
     gnuBar(numGroups+1) = zero

     return

  end function QuadratureList_getGroupAverageEnergy

!=======================================================================
! setCounters interface
!=======================================================================

  subroutine QuadratureList_setCounters(self)


!   Passed variables

    type(QuadratureList), intent(inout) :: self

!   Local
    type(CommSet),        pointer       :: CommSetPtr
    integer                             :: setID

    TETON_CHECK_BOUNDS1(self%CommSetPtr, self%nCommSets)
    do setID=1,self% nCommSets 
      CommSetPtr => self% CommSetPtr(setID)
      CommSetPtr% fluxSweeps = 0
    enddo


    return

  end subroutine QuadratureList_setCounters

!-----------------------------------------------------------------------
  function QuadratureList_getSweeps(self) result(nTotalSweeps)

!    Return the total number of transport sweeps performed 

     use Size_mod

  
!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nTotalSweeps

!    Local
     type(CommSet), pointer           :: CommSetPtr 
     integer                          :: setID 
     integer                          :: nSweeps 
     integer                          :: Groups
     integer                          :: totalAngles
     real(adqt)                       :: fracSweep

     nSweeps     = 0     
     totalAngles = 0

     TETON_CHECK_BOUNDS1(self%CommSetPtr, self%nCommSets)
     do setID=1,self% nCommSets 
       CommSetPtr     => self% CommSetPtr(setID)
       Groups         =  CommSetPtr% Groups
       nSweeps        =  nSweeps + Groups*CommSetPtr% fluxSweeps 
       totalAngles    =  totalAngles + Groups*CommSetPtr% NumAngles
       fracSweep      =  real(CommSetPtr% fluxSweeps,adqt)/real(CommSetPtr% NumAngles, adqt)
     enddo

     nTotalSweeps = nSweeps/totalAngles


     return

  end function QuadratureList_getSweeps

end module QuadratureList_mod
