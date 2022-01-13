!***********************************************************************
!                        Version 1:  09/96, PFN                        *
!                                                                      *
!   SNFLWXYZ - This routine, called by RSWPMD and RTACCELMD, solves    *
!              the fixed-source transport problem on an arbitrary      *
!              grid in either xyz-geometry or rz-geometry.             *
!              An upstream corner-balance spatial discretization is    *
!              used.                                                   *
!                                                                      *
!***********************************************************************

   subroutine GTASweep(GTAsetID, PsiB)


   use kind_mod
   use constant_mod
   use Size_mod
   use Quadrature_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use SetData_mod
   use CommSet_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,    intent(in)    :: GTAsetID
   real(adqt), intent(inout) :: PsiB(Size%nbelem,Size%nangGTA)

!  Local

   type(SetData),    pointer :: Set
   type(CommSet),    pointer :: CSet
   type(AngleSet),   pointer :: ASet

   integer                   :: nSets
   integer                   :: nCommSets
   integer                   :: nAngleSets
   integer                   :: setID
   integer                   :: cSetID
   integer                   :: aSetID
   integer                   :: Angle
   integer                   :: ndim 
   integer                   :: sendIndex

   logical (kind=1)          :: SnSweep 

!  Constants

   nSets      = getNumberOfSets(Quad)
   nCommSets  = getNumberOfCommSets(Quad)
   nAngleSets = getNumberOfAngleSets(Quad)

   setID      = nSets      + GTAsetID
   cSetID     = nCommSets  + GTAsetID
   aSetID     = nAngleSets + GTAsetID

   Set     => getSetData(Quad, setID)
   CSet    => getCommSetData(Quad, cSetID)
   ASet    => getAngleSetData(Quad, aSetID)

   ndim    =  Size% ndim
   SnSweep = .FALSE.

!  Restore the initial communication order
   call restoreCommOrder(CSet)

!  Post receives for all data
   call InitExchange(cSetID)

!  Initialize partial scalar correction

   if (.not. Size% useNewGTASolver) then
     Set% tPhi(:) = zero
   endif

   if (ndim == 2) then
     Set% tPsiM(:) = zero
     Set% tInc(:)  = zero
   endif

!  Loop over angles, solving for each in turn:

   AngleLoop: do sendIndex=1,Set% NumAnglesDyn

     Angle = CSet% AngleOrder(sendIndex)

!    Send the boundary information needed by my neighbors
     call SendFlux(SnSweep, cSetID, sendIndex, PsiB)

!    Test for completion of the sends needed by my neighbors
     call TestSend(cSetID, sendIndex)

!    Receive the boundary information needed to compute this angle
     call RecvFlux(SnSweep, cSetID, Angle, PsiB)


!    Sweep the mesh, calculating PSI for each corner; the 
!    boundary flux array PSIB is also updated here. 
!    Mesh cycles are fixed automatically.

     AngleType: if ( .not. ASet% FinishingDirection(Angle) ) then

       call snreflect(SnSweep, setID, Angle, PsiB)

       if (Size% useNewGTASolver) then

         if (ndim == 3) then
           call SweepGreyUCBxyzNEW(setID, Angle, PsiB)
         elseif (ndim == 2) then
           call SweepGreyUCBrzNEW(setID, Angle, PsiB)
         endif

       else

         if (ndim == 3) then
           call SweepGreyUCBxyz(setID, Angle, PsiB)
         elseif (ndim == 2) then
           call SweepGreyUCBrz(setID, Angle, PsiB)
         endif

       endif

     endif AngleType

   enddo AngleLoop


   return
   end subroutine GTASweep 


