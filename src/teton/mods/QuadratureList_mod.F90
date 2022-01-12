! QuadratureList Module:  Contains a list of source profiles
#include "macros.h"
module QuadratureList_mod

  use kind_mod
  use Quadrature_mod
  use SetData_mod
  use RadIntensity_mod
  use AngleSet_mod
  use GroupSet_mod
  use CommSet_mod
  use Size_mod

  private

! public interfaces

  public construct 
  public destruct  
  public setQuadrature
  public getQuadrature 
  public getNumQuadSets 
  public getNumSnSets
  public getNumberOfSets
  public getNumberOfGTASets
  public getNumberOfAngleSets
  public getNumberOfGroupSets
  public getNumberOfCommSets
  public getNumberOfZoneSets
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
  public getRadIntensity
  public getZoneRadiationEnergy
  public getSetIDfromGroupAngle

  type, public :: QuadratureList

     ! These first three variables refer to quadrature sets and NOT phase-space sets.
     integer                         :: NumQuadSets  ! Total number of quadrature sets 
     integer                         :: NumSnSets    ! Number of quadrature sets used for Sn sweeps
     integer                         :: nGTASets     ! Number of quadrature sets used for GTA

     integer                         :: nSets
     ! nSets is the maximum number of sets we can break the Sn quadrature sets
     !   into, as determined by the user input nSetsMaster and also hardware
     !   constraints.
     ! WARNING: nSets may NOT actually be the # of sets used.  It is generally
     !   an upper bound on the number of sets Teton will use in its GPU sweeps.
     !   At the point that the value of nSets is determined, we have not yet
     !   begun decomposing the SN problem into sets.

     integer                         :: nHyperDomains
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
     type(RadIntensity), pointer     :: RadIntPtr(:)  => null()

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

  interface getNumQuadSets
    module procedure QuadratureList_getNumQuadSets
  end interface

  interface getNumSnSets
    module procedure QuadratureList_getNumSnSets
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

  interface getRadIntensity
    module procedure QuadratureList_getRadIntensity
  end interface

  interface getZoneRadiationEnergy
    module procedure QuadratureList_getZoneRadEnergy
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

  subroutine QuadratureList_ctor(self, NumQuadSets, ngr,   &
                                 totalAngles, nSetsMaster, nSets)

    use, intrinsic :: iso_c_binding, only : c_int                                                                                                                                                                  
    use constant_mod
    use Options_mod

    implicit none

!   Passed variables

    type(QuadratureList), intent(inout) :: self
    integer,              intent(in)    :: NumQuadSets 
    integer,              intent(in)    :: ngr
    integer,              intent(in)    :: totalAngles
    integer,              intent(in)    :: nSetsMaster
    integer,              intent(inout) :: nSets
!   Warning: nSets may not be the value you think it should be.  See comments
!     above in the QuadratureList type definition

!   Local

    integer :: nOmpMaxThreads
    integer :: nOmpMaxTeams

    TETON_VERIFY(NumQuadSets == 2, "Teton: Support for multiple Sn quadrature sets is deprecated.")
    ! Support for multiple Sn quad sets is being removed.
    ! See Paul's issue for more information
    ! https://rzlc.llnl.gov/gitlab/deterministic-transport/TRT/Teton/-/issues/43

    self% NumQuadSets   = 2
    self% NumSnSets     = 1
    self% nSets         = 1
    self% nGroupSets    = 1
    self% nZoneSets     = 1
    self% nHyperDomains = 1

#if defined(TETON_ENABLE_OPENMP)
    if (Size%useGPU) then
       nOmpMaxTeams = Options%getNumOmpMaxTeams()
       self% nSets       = max(self% NumSnSets, nOmpMaxTeams)
       self% nZoneSets   = min(Size%nzones, nOmpMaxTeams)
    else
       nOmpMaxThreads = Options%getNumOmpMaxThreads()
       self% nSets       = max(self% NumSnSets, nOmpMaxThreads)
       self% nZoneSets   = min(Size%nzones, nOmpMaxThreads)
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

    allocate( self% SetIDList(ngr,totalAngles) )
    allocate( self% QuadPtr(self% NumQuadSets) )

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
      allocate( self% RadIntPtr(self% nSets) )
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
                                QuadID,        &
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

    integer,    intent(in)              :: QuadID
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
 
    call construct(self% QuadPtr(QuadID), &
                         QuadID,          &
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
! getSetIDfromGroup interface
!=======================================================================

  function QuadratureList_getSetIDfromGroupAngle(self, group, angle) result(setID)

    implicit none

!   Passed variables

    type(QuadratureList),  intent(inout) :: self
    integer,               intent(in)    :: group
    integer,               intent(in)    :: angle

    integer                              :: setID

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
    deallocate( self% RadIntPtr )

    deallocate( self% angleID )
    deallocate( self% groupID )
    deallocate( self% commID )

    return

  end subroutine QuadratureList_dtor

!-----------------------------------------------------------------------
  function QuadratureList_getNumQuadSets(self) result(NumQuadSets)
                                                                                            
!    Returns the number of unique quadrature sets
!      NumQuadSets   number of quadrature sets

!    variable declarations
     implicit none
                                                                                            
!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: NumQuadSets

     NumQuadSets = self% NumQuadSets
                                                                                            
     return
                                                                                            
  end function QuadratureList_getNumQuadSets

!-----------------------------------------------------------------------
  function QuadratureList_getNumSnSets(self) result(NumSnSets)
                                                                                            
!    Returns the number of non-acceleration quadrature sets
!      NumSnSets   number of Sn sets
                                                                                            
!    variable declarations
     implicit none
                                                                                            
!    passed variables
     type(QuadratureList), intent(in) :: self
     integer                          :: NumSnSets
                                                                                            
     NumSnSets = self% NumSnSets
                                                                                            
     return
                                                                                            
  end function QuadratureList_getNumSnSets

!-----------------------------------------------------------------------
  function QuadratureList_getNumberOfSets(self) result(nSets)

!    Returns the number of non-acceleration quadrature sets
!      NumSnSets   number of group/anglesets

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
                                                                                           
                                                                                           
     QuadPtr => self% QuadPtr(self% NumSnSets + 1)
                                                                                           
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

     localID = self% angleID(setID)


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

     localID = self% groupID(setID)


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

     localID = self% commID(setID)


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


     SetDataPtr => self% SetDataPtr(self% nSets + setID)

     return

  end function QuadratureList_getGTASetData

!=======================================================================
! access the radiation intensity variables for a Set 
!=======================================================================

  function QuadratureList_getRadIntensity(self,setID) result(RadIntPtr)

    implicit none

!   Passed variables
    type(QuadratureList), intent(in) :: self
    integer,              intent(in) :: setID
    type(RadIntensity),   pointer    :: RadIntPtr


    RadIntPtr => self% RadIntPtr(setID)


    return

  end function QuadratureList_getRadIntensity

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
     integer                      :: i, ig, ng
     real(adqt)                   :: GrpBnds(numGroups+1)
                                                                                            
     ng = 0
     do i=1,self % NumSnSets
       QuadPtr => self% QuadPtr(i)
       do ig=1,QuadPtr% Groups+1
         GrpBnds(ng+ig) = QuadPtr% Gnu(ig)
       enddo
       ng = ng + QuadPtr% Groups
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
     integer                      :: qset
     integer                      :: g
     integer                      :: ng
     real(adqt)                   :: gnuBar(0:numGroups+1)

     ng = 0
     do qset=1,self% NumSnSets
       QuadPtr => self% QuadPtr(qset)
       do g=1,QuadPtr% Groups
         gnuBar(ng+g) = QuadPtr% gnuBar(g)
       enddo
       ng = ng + QuadPtr% Groups
     enddo

     gnuBar(0)           = zero
     gnuBar(numGroups+1) = zero

     return

  end function QuadratureList_getGroupAverageEnergy

!=======================================================================
! get the zone radiation energy, summed over sets 
!=======================================================================

  function QuadratureList_getZoneRadEnergy(self, zoneID) result(radEnergy)
   
    use constant_mod
    use RadIntensity_mod

    implicit none

!   Passed variables
    type(QuadratureList),     intent(in) :: self
    integer,                  intent(in) :: zoneID
    real(adqt)                           :: radEnergy

!   Local
    integer                          :: setID
    type(RadIntensity),     pointer  :: RadIntPtr

!   Return the zonal radiation energy 

    radEnergy = zero
    do setID=1,self% nSets
      RadIntPtr => self% RadIntPtr(setID)

      radEnergy = radEnergy + RadIntPtr% radEnergy(zoneID) 
    enddo


    return

  end function QuadratureList_getZoneRadEnergy 

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
