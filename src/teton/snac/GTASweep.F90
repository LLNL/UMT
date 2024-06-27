!***********************************************************************
!                        Version 1:  09/2017, PFN                      *
!                                                                      *
!   GTASweep      - This routine controls the GTA set sweeps           *
!                   communication when running on a CPU.               *
!                                                                      *
!***********************************************************************

   subroutine GTASweep(P, PsiB)


   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use Geometry_mod
   use SetData_mod
   use CommSet_mod
   use AngleSet_mod

   implicit none

!  Arguments

   real(adqt), intent(in)       :: P(Size%ncornr)
   real(adqt), intent(inout)    :: PsiB(Size%nbelem,Size%nangGTA)

!  Local

   type(SetData),    pointer :: Set
   type(CommSet),    pointer :: CSet
   type(AngleSet),   pointer :: ASet

   integer                   :: nGTASets
   integer                   :: nCommSets
   integer                   :: nZoneSets
   integer                   :: setID
   integer                   :: cSetID
   integer                   :: zSetID
   integer                   :: Angle
   integer                   :: ndim 
   integer                   :: sendIndex
   integer                   :: NumAnglesDyn
   integer                   :: GTASetID
   integer                   :: c

   real(adqt)                :: wtiso

   logical (kind=1)          :: SnSweep 

!  Constants

   nGTASets     =  getNumberOfGTASets(Quad)
   nCommSets    =  getNumberOfCommSets(Quad)
   nZoneSets    =  getNumberOfZoneSets(Quad)

   Set          => getGTASetData(Quad, 1) 
   NumAnglesDyn =  Set% NumAnglesDyn

   ndim         =  Size% ndim
   SnSweep      = .FALSE.
   wtiso        =  Size% wtiso


!  Initialize Communication

   do cSetID=nCommSets+1,nCommSets+nGTASets
     CSet => getCommSetData(Quad, cSetID)

!    Restore the initial communication order
     call restoreCommOrder(CSet)

!    Post receives for all data
     call InitExchange(cSetID)
   enddo

!  Set the TSA source before the sweeps

   if ( GTA% ID == 1) then

     !$omp  parallel do default(none)  &
     !$omp& shared(nZoneSets, Geom, GTA, wtiso, P)
     ZoneSetLoop1: do zSetID=1,nZoneSets
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         GTA% TsaSource(c) = wtiso*(GTA% GreySigScat(c)*P(c) + &
                             GTA% GreySource(c)) 
       enddo
     enddo ZoneSetLoop1
     !$omp end parallel do

   else

     !$omp  parallel do default(none)  &
     !$omp& shared(nZoneSets,Geom, GTA, wtiso, P)
     ZoneSetLoop2: do zSetID=1,nZoneSets
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         GTA% TsaSource(c) = wtiso*GTA% GreySigScat(c)*GTA%eps(c)*  &
                             (P(c) - GTA%OldGreyCorrection(c))
       enddo
     enddo ZoneSetLoop2
     !$omp end parallel do

   endif

!  Initialize partial scalar correction

   if (Size% useNewGTASolver) then
     GTA% PhiInc(:) = zero
   else
     do GTAsetID=1,nGTASets
       Set  => getGTASetData(Quad, GTAsetID)
       Set% tPhi(:) = zero
     enddo
   endif

   if (ndim == 2) then
     do GTAsetID=1,nGTASets
       Set  => getGTASetData(Quad, GTAsetID)
       Set% tPsiM(:) = zero
       Set% tInc(:)  = zero
     enddo
   endif

!  Loop over angles, solving for each in turn:

   AngleLoop: do sendIndex=1,Set% NumAnglesDyn

     do cSetID=nCommSets+1,nCommSets+nGTASets

       CSet  => getCommSetData(Quad, cSetID)
       Angle =  CSet% AngleOrder(sendIndex)
       setID =  CSet% set1

!      Send the boundary information needed by my neighbors
       call SendFlux(SnSweep, cSetID, sendIndex, PsiB)

!      Test for completion of the sends needed by my neighbors
       call TestSend(cSetID, sendIndex)

!      Receive the boundary information needed to compute this angle
       call RecvFlux(SnSweep, cSetID, Angle, PsiB)

!      Sweep the mesh, calculating PSI for each corner; the 
!      boundary flux array PSIB is also updated here. 
!      Mesh cycles are fixed automatically.

       ASet => getAngleSetData(Quad, cSetID)

       AngleType: if ( .not. ASet% FinishingDirection(Angle) ) then

         call snreflect(SnSweep, setID, Angle, PsiB)

         if (ndim == 3) then
           call SweepGreyUCBxyz(setID, Angle, PsiB)
         elseif (ndim == 2) then
           call SweepGreyUCBrz(setID, Angle, PsiB)
         endif

       endif AngleType

     enddo

   enddo AngleLoop


   return
   end subroutine GTASweep 


