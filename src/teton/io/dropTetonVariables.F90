!*****************************************************************************
! dropTetonVariables                                                         *
!                                                                            *
! Called from the host code to write the angle-dependent radiation intensity *
! (Psi) and corner electron temperature, along with run-time controls,       *
! to a SILO restart file.                                                    *
!                                                                            *
!*****************************************************************************

   subroutine dropTetonVariables(filePtrID, lenTe, cornerTe, success) &
        BIND(C,NAME="teton_dropvariables")

   USE ISO_C_BINDING
   use kind_mod
   use SetData_mod
   use QuadratureList_mod
   use iter_control_mod
   use iter_control_list_mod

   implicit none

#include "silo_f9x.inc"

!  Arguments

   integer(C_INT),    intent(in)      :: filePtrID
   integer(C_INT),    intent(in)      :: lenTe
   real(C_DOUBLE),    intent(in)      :: cornerTe(lenTe)
   integer(C_INT),    intent(inout)   :: success

!  Local

   type(SetData),      pointer :: Set

   integer                     :: setID
   integer                     :: nSets  
   integer                     :: nAngleSets
   integer                     :: nGroupSets
   integer                     :: nCommSets
   integer                     :: nGTASets
   integer                     :: nDims
   integer                     :: dims(3)
   integer                     :: setDesc(6)
   integer                     :: maxIters(5)
   integer                     :: lenName
   integer                     :: d1
   integer                     :: d10
   integer                     :: d100
   integer                     :: dirstatus

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

   real(C_DOUBLE)              :: relTols(5)

   type(IterControl), pointer  :: itCon => NULL()

!  Loop over sets and write the angle-dependent intensity and 
!  set descriptors to the restart file

   success = -1

!  Give Teton it's own directory
   lenName = 8
   success =  dbmkdir(filePtrID, TetonDirName, lenName, dirstatus)
   success = dbsetdir(filePtrID, TetonDirName, lenName)

   if (success < 0) return

   nSets   = getNumberOfSets(Quad)
   lenName = 15
   dims(1) = 1
   nDims   = 1

   success = dbwrite(filePtrID, nSetName, lenName, nSets, dims(1), &
                     nDims, DB_INT)

   nAngleSets = getNumberOfAngleSets(Quad)
   lenName    = 10
   dims(1)    = 1
   nDims      = 1

   success = dbwrite(filePtrID, nAngleSetName, lenName, nAngleSets, dims(1), &
                     nDims, DB_INT)

   nGroupSets = getNumberOfGroupSets(Quad)
   lenName    = 10 
   dims(1)    = 1
   nDims      = 1

   success = dbwrite(filePtrID, nGroupSetName, lenName, nGroupSets, dims(1), &
                     nDims, DB_INT)

   nCommSets  = getNumberOfCommSets(Quad)
   lenName    = 9 
   dims(1)    = 1
   nDims      = 1

   success = dbwrite(filePtrID, nCommSetName, lenName, nCommSets, dims(1), &
                     nDims, DB_INT)

   nGTASets = getNumberOfGTASets(Quad)
   lenName  = 8 
   dims(1)  = 1
   nDims    = 1

   success = dbwrite(filePtrID, nGTASetName, lenName, nGTASets, dims(1), &
                     nDims, DB_INT)

   if (success < 0) return

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

!    Make a new directory for each phase-space set
     success =  dbmkdir(filePtrID, directoryName, lenName, dirstatus)
     success = dbsetdir(filePtrID, directoryName, lenName)

     if (success < 0) return

!    Write Radiation Intensity
     lenName  = 3 
     nDims    = 3
     dims(1)  = Set% Groups
     dims(2)  = Set% nCorner
     dims(3)  = Set% NumAngles 

     success  = dbwrite(filePtrID, PsiName, lenName, Set% Psi, dims, &
                        nDims, DB_DOUBLE)

     if (success < 0) return

!    Write the set descriptors 

     nDims      = 1 
     dims(1)    = 6
     lenName    = 14
     setDesc(1) = Set% Groups 
     setDesc(2) = Set% nCorner 
     setDesc(3) = Set% NumAngles
     setDesc(4) = Set% g0
     setDesc(5) = Set% angle0
     setDesc(6) = Set% QuadID

     success    = dbwrite(filePtrID, setDescName, lenName, setDesc, dims, &
                          nDims, DB_INT)

     if (success < 0) return

   enddo SetLoop

!  Write corner electron temperature

   lenName = 8
   success = dbsetdir(filePtrID, TetonDirName, lenName)

   nDims   = 1
   dims(1) = lenTe
   lenName = 8

   success = dbwrite(filePtrID, TeName, lenName, cornerTe, dims(1), & 
                     nDims, DB_DOUBLE )

   if (success < 0) return

!  Write iteration controls.

   itCon => getIterationControl(IterControls,"temperature")
   relTols(1)  =  getEpsilonPoint(itCon)
   maxIters(1) =  getMaxNumberOfIterations(itCon)
   itCon => getIterationControl(IterControls,"intensity")
   relTols(2)  =  getEpsilonPoint(itCon)
   ! Ignored, but setting it to something so memory is initialized.
   maxIters(2) = 1
   itCon => getIterationControl(IterControls,"grey")
   relTols(3)  =  getEpsilonPoint(itCon)
   maxIters(3) =  getMaxNumberOfIterations(itCon)
   itCon => getIterationControl(IterControls,"incidentFlux")
   relTols(4)  =  getEpsilonPoint(itCon)
   maxIters(4) =  getMaxNumberOfIterations(itCon)
   itCon => getIterationControl(IterControls,"nonLinear")
   relTols(5)  =  getEpsilonPoint(itCon)
   maxIters(5) =  getMaxNumberOfIterations(itCon)
   nDims   = 1
   dims(1) = 5
   lenName = 7

   success = dbwrite(filePtrID, relTolsName, lenName, relTols, dims(1), &
                     nDims, DB_DOUBLE )

   if (success < 0) return

   lenName = 8
   success = dbwrite(filePtrID, maxItersName, lenName, maxIters, dims(1), &
                     nDims, DB_INT )

   return
   end subroutine dropTetonVariables

