!*****************************************************************************
! loadTetonVariables                                                         *
!                                                                            *
! Called from the host code to read the angle-dependent radiation intensity  *
! (Psi) and corner electron temperature, along with runtime controls,        *
! from a SILO restart file.                                                  *
!                                                                            *
! Teton also requires the last cycle's Volume field on a restart load, as    *
! Teton needs a delta Volume across cycles when advancing cycles.            *
!                                                                            *
! This field is not being explictly saved/loaded on this SILO restart.       *
! Instead, it requires that the host code update Teton's mesh positions and  *
! call teton_getvolume() to populate the Volume field with the values from   *
! the restarted cycle.                                                       *
!*****************************************************************************

   subroutine loadTetonVariables(filePtrID, lenTe,  &
                                 cornerTe, nSetsMaster, success) &
                                 BIND(C,NAME="teton_loadvariables")

   USE ISO_C_BINDING
   use kind_mod
   use SetData_mod
   use QuadratureList_mod
   use MemoryAllocator_mod
   use iter_control_mod
   use iter_control_list_mod
   use Datastore_mod

   implicit none

#include "silo_f9x.inc"

!  Arguments

   integer(C_INT),          intent(in)    :: filePtrID
   integer(C_INT),          intent(in)    :: lenTe
   real(C_DOUBLE),          intent(inout) :: cornerTe(lenTe)
   integer(C_INT),          intent(inout) :: nSetsMaster
   integer(C_INT),          intent(inout) :: success

!  Local

   type(SetData),      pointer :: Set

   integer                     :: setID
   integer                     :: nSets  
   integer                     :: nAngleSets
   integer                     :: nGroupSets
   integer                     :: nCommSets
   integer                     :: nGTASets
   integer                     :: lenName
   integer                     :: setDesc(6)
   integer                     :: maxIters(5)
   integer                     :: d1
   integer                     :: d10
   integer                     :: d100

   character (len=15)          :: nSetName = "nPhaseSpaceSets"
   character (len=10)          :: nAngleSetName = "nAngleSets"
   character (len=10)          :: nGroupSetName = "nGroupSets"
   character (len=9)           :: nCommSetName = "nCommSets"
   character (len=8)           :: nGTASetName = "nGTASets"
   character (len=8)           :: TeName   = "cornerTe"
   character (len=8)           :: TetonDirName = "/_Teton/"
   character (len=12)          :: BaseName = "/_Teton/Set_"
   character (len=10)          :: setNum   = "1234567890"
   character (len=16)          :: directoryName
   character (len=14)          :: setDescName = "setDescription"
   character (len=3)           :: PsiName = "Psi"
   character (len=7)           :: relTolsName = "relTols"
   character (len=8)           :: maxItersName = "maxIters"
   character (len=1)           :: endName = "/"
   character (len=19)          :: setLabel 

   logical                     :: usePinnedMemory

   real(C_DOUBLE)              :: relTols(5)

   type(IterControl), pointer  :: itCon => NULL()

   ! Assume that if host pinned memory allocator was provided, then, we should use
   ! it.  We can't check Size%useGPU, as that module isn't constructed until
   ! later in the restart load process.
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
   usePinnedMemory = .TRUE.
#else
   usePinnedMemory = .FALSE.
#endif

!  Loop over sets and read the angle-dependent intensity and
!  set descriptors from the restart file 

   success = -1

   lenName = 8
   success = dbsetdir(filePtrID, TetonDirName, lenName)

   !! NOTE: The failure to find a directory "_TETON" is ignored, since
   !!       we might be turned off and on in the middle of a simulation.
   if (success < 0) then 
       success = 4
       return
   end if

   lenName = 15
   success = dbrdvar(filePtrID, nSetName, lenName, nSets)

   lenName = 10
   success = dbrdvar(filePtrID, nAngleSetName, lenName, nAngleSets)

   lenName = 10
   success = dbrdvar(filePtrID, nGroupSetName, lenName, nGroupSets)

   lenName = 9 
   success = dbrdvar(filePtrID, nCommSetName, lenName, nCommSets)

   lenName = 8
   success = dbrdvar(filePtrID, nGTASetName, lenName, nGTASets)

   if (success < 0) return

   nSetsMaster = nSets

!  Allocate Memory if the Teton object has not been created

   if ( .not. associated(Quad) ) then
     allocate (Quad)
     call constructSetPointers(Quad, nSets, nAngleSets, nGroupSets,  &
                               nCommSets, nGTASets)
   endif

   SetLoop: do setID=1,nSets

     Set  => getSetData(Quad, SetID)

     d1   =  mod(setID, 10)
     d10  =  setID/10

     if (d10 >= 10) then
       d100 = d10/10
       d10  = d10 - 10
     else
       d100 = 10
     endif

     if (d10 == 0) then
       d10 = 10
     endif

     if (d1 == 0) then
       d1 = 10
     endif

     directoryName =  &
        BaseName//setNum(d100:d100)//setNum(d10:d10)//setNum(d1:d1)//endName

     lenName       = 16

!    Set the directory to read from
     success = dbsetdir(filePtrID, directoryName, lenName)

     if (success < 0) return

!    Read the set descriptors

     lenName = 14

     success = dbrdvar(filePtrID, setDescName, lenName, setDesc)

     if (success < 0) return

     Set% Groups    = setDesc(1)
     Set% nCorner   = setDesc(2)
     Set% NumAngles = setDesc(3)
     Set% g0        = setDesc(4)
     Set% angle0    = setDesc(5)
     Set% QuadID    = setDesc(6)

     if ( .not. associated( Set% Psi ) ) then
        ! Make sure this stays consistent with how the set string label is
        ! defined in the Set constructor.
        write(setLabel, '(I0.3)') setID
        setLabel = "phase_space_set_"//setLabel

        call Allocator%allocate(usePinnedMemory, setLabel,"Psi", Set% Psi, Set% Groups,Set% nCorner,Set% NumAngles)
     endif

!    Read Radiation Intensity
     lenName = 3 

     success = dbrdvar(filePtrID, PsiName, lenName, Set% Psi)

     if (success < 0) return

   enddo SetLoop

!  Read corner electron temperature
   lenName = 8
   success = dbsetdir(filePtrID, TetonDirName, lenName)

   lenName = 8

   success = dbrdvar(filePtrID, TeName, lenName, cornerTe)

   lenName = 7
   success = dbrdvar(filePtrID, relTolsName, lenName, relTols)
   if (success < 0) return

   lenName = 8
   success = dbrdvar(filePtrID, maxItersName, lenName, maxIters)
   if (success < 0) return

   itCon => getIterationControl(IterControls,"temperature")
   call setControls(itCon, epsilonPoint=relTols(1), maxNumberOfIterations=maxIters(1))
   itCon => getIterationControl(IterControls,"intensity")
   ! Intensity control max iterations ignored, do not set it.
   call setControls(itCon, epsilonPoint=relTols(2))
   itCon => getIterationControl(IterControls,"grey")
   call setControls(itCon, epsilonPoint=relTols(3), maxNumberOfIterations=maxIters(3))
   itCon => getIterationControl(IterControls,"incidentFlux")
   call setControls(itCon, epsilonPoint=relTols(4), maxNumberOfIterations=maxIters(4))
   itCon => getIterationControl(IterControls,"nonLinear")
   call setControls(itCon, epsilonPoint=relTols(5), maxNumberOfIterations=maxIters(5))

   return
   end subroutine loadTetonVariables

