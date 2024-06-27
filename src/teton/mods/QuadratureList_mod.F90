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


    use, intrinsic :: iso_c_binding, only : c_int
    use cmake_defines_mod, only : omp_device_num_processors, omp_device_team_thread_limit
    use Options_mod, only : Options
    use constant_mod

    implicit none

!   Passed variables

    type(QuadratureList), intent(inout) :: self
    integer,              intent(in)    :: nAnglesSn 
    integer,              intent(in)    :: nSetsMaster
    integer,              intent(inout) :: nSets
!   Warning: nSets may not be the value you think it should be.  See comments
!     above in the QuadratureList type definition

!   Local

    integer :: nOmpMaxThreads
    integer :: nOmpMaxTeams


    self% nSets             = 1
    self% nGroupSets        = 1
    self% nZoneSets         = 1
    self% nHyperElements(:) = 0 

#if defined(TETON_ENABLE_OPENMP)
    if (Size%useGPU) then
       ! If running on gpu, set these to use all available gpu processors.
       nOmpMaxTeams    = omp_device_num_processors
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

    implicit none

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

    implicit none

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

    implicit none

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

    implicit none

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
  function QuadratureList_getNumberOfSets(self) result(nSets)

!    Returns the number of phase-space sets (nSets)

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nSets

     nSets = self% nSets

     return

  end function QuadratureList_getNumberOfSets

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfGTASets(self) result(nGTASets)

!    Returns the number of GTA quadrature sets
!      nGTASets   number of group/anglesets

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nGTASets

     nGTASets = self% nGTASets

     return

  end function QuadratureList_getNumberOfGTASets

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfAngleSets(self) result(nAngleSets)

!    Returns the number of angle sets
!      nAngleSets   number of angle sets

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nAngleSets

     nAngleSets = self% nAngleSets

     return

  end function QuadratureList_getNumberOfAngleSets

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfGroupSets(self) result(nGroupSets)

!    Returns the number of energy group sets
!      nGroupSets   number of group sets

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nGroupSets

     nGroupSets = self% nGroupSets

     return

  end function QuadratureList_getNumberOfGroupSets

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfCommSets(self) result(nCommSets)

!    Returns the number of communication sets
!      nCommSets   number of communication sets

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nCommSets

     nCommSets = self% nCommSets

     return

  end function QuadratureList_getNumberOfCommSets

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfZoneSets(self) result(nZoneSets)

!    Returns the number of zone sets
!      nZoneSets   number of zone sets

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: nZoneSets

     nZoneSets = self% nZoneSets

     return

  end function QuadratureList_getNumberOfZoneSets

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfHyperDomains(self, ID) result(nHyperDomains)

!    Returns the number of hyper domains used for sweeps 
!      nHyperDomains   number of hyper domains 
!      ID = 1   High-order angle set
!      ID = 2   GTA angle set

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: ID
     integer                          :: nHyperDomains 

     TETON_ASSERT(ID <= 2, "Invalid ID for getNumberOfHyperDomains, must be 1 or 2")
     nHyperDomains = self% nHyperDomains(ID) 

     return

  end function QuadratureList_getNumberOfHyperDomains

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfHyperElements(self, ID) result(nHyperElements)

!    Returns the number of interface elements at hyper-domain
!    boundaries  used for sweeps 
!      nHyperElements   number of interface elements 
!      ID = 1   High-order angle set
!      ID = 2   GTA angle set

!    variable declarations
     implicit none

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
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: QuadID 
     type(Quadrature),     pointer    :: QuadPtr

            
     TETON_CHECK_BOUNDS1(self%QuadPtr, QuadID)
     QuadPtr => self% QuadPtr(QuadID)
            
     return

  end function QuadratureList_getQuad 

!-----------------------------------------------------------------------
  function QuadratureList_getGTAQuad(self) result(QuadPtr)
                                                                                           
!    Return a pointer to the GTA quadrature set
!      QuadPtr  pointer to the quadrature set
                                                                                           
!    variable declarations
     implicit none
                                                                                           
!    passed variables
     type(QuadratureList), intent(in) :: self
     type(Quadrature),     pointer    :: QuadPtr
                                                                                           
                                                                                           
     QuadPtr => self% QuadPtr(2)
                                                                                           
     return
                                                                                           
  end function QuadratureList_getGTAQuad

!-----------------------------------------------------------------------
  function QuadratureList_getSNQuad(self) result(QuadPtr)
                                                                                           
!    Return a pointer to the SN quadrature set
!      QuadPtr  pointer to the quadrature set
                                                                                           
!    variable declarations
     implicit none
                                                                                           
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
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     type(SetData),        pointer    :: SetDataPtr 


     TETON_CHECK_BOUNDS1(self%SetDataPtr, setID)
     SetDataPtr => self% SetDataPtr(setID)

     return

  end function QuadratureList_getSetData

!-----------------------------------------------------------------------
  function QuadratureList_getAngleSetData(self,angleSetID) result(AngSetPtr)

!    Return a pointer to an angle set
!      angleSetID  angle set ID number
!      AngSetPtr   pointer to the angle set

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: angleSetID
     type(AngleSet),       pointer    :: AngSetPtr


     TETON_CHECK_BOUNDS1(self%AngSetPtr, angleSetID)
     AngSetPtr => self% AngSetPtr(angleSetID)

     return

  end function QuadratureList_getAngleSetData

!-----------------------------------------------------------------------
  function QuadratureList_getGroupSetData(self,groupSetID) result(GrpSetPtr)

!    Return a pointer to an energy group set 
!      groupSetID  group set ID number    
!      GrpSetPtr   pointer to the group set 

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: groupSetID
     type(GroupSet),       pointer    :: GrpSetPtr


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
     implicit none

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
     implicit none

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
     implicit none

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
     implicit none

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
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     type(SetData),        pointer    :: SetDataPtr


     TETON_CHECK_BOUNDS1(self%SetDataPtr, self% nSets + setID)
     SetDataPtr => self% SetDataPtr(self% nSets + setID)

     return

  end function QuadratureList_getGTASetData

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfGroups(self,setID) result(Groups)

!    Return the number of groups in this group/angle set 
!      setID       set ID number 
!      Groups      number of energy groups 

!    variable declarations
     implicit none

!    passed variables
     type(QuadratureList), intent(in) :: self
     integer,              intent(in) :: setID
     integer                          :: Groups
     type(SetData),        pointer    :: SetDataPtr


     TETON_CHECK_BOUNDS1(self%SetDataPtr, setID)
     SetDataPtr => self% SetDataPtr(setID)
     Groups     =  SetDataPtr% Groups

     return

  end function QuadratureList_getNumberOfGroups

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfAngles(self,setID) result(NumAngles)

!    Return the number of angles in this group/angle set 
!      setID       set ID number 
!      NumAngles   number of angles 

!    variable declarations
     implicit none

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
     implicit none
                                                                                            
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
     implicit none

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

    implicit none

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

     implicit none
  
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
