#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Version 1:  09/2017, PFN                      *
!                                                                      *
!   GTASweep_GPU  - This routine controls the GTA set sweeps           *
!                   communication when running on a GPU.               *
!                                                                      *
!***********************************************************************

   subroutine GTASweep_GPU(P, PsiB, withSource)

   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use GreyAcceleration_mod
   use SetData_mod
   use CommSet_mod

   implicit none

!  Arguments

   real(adqt), intent(in)       :: P(Size%ncornr)
   real(adqt), intent(inout)    :: PsiB(Size% nSurfElem,Size%nangGTA)

   logical (kind=1), intent(in) :: withSource

!  Local

   type(SetData),    pointer :: Set
   type(CommSet),    pointer :: CSet

   integer                   :: c
   integer                   :: setID
   integer                   :: cSetID
   integer                   :: zSetID
   integer                   :: nSets
   integer                   :: nGTASets
   integer                   :: nZoneSets
   integer                   :: nCommSets

   integer                   :: Angle
   integer                   :: ndim 
   integer                   :: sendIndex

   real(adqt)                :: wtiso

   logical (kind=1)          :: SnSweep 

!  Constants
!
   nSets     = getNumberOfSets(Quad)
   nGTASets  = getNumberOfGTASets(Quad)
   nZoneSets = getNumberOfZoneSets(Quad)
   nCommSets = getNumberOfCommSets(Quad)
   ndim      = Size% ndim
   SnSweep   = .FALSE.
   wtiso     = Size% wtiso

   Set       => getSetData(Quad, nSets+1)

!  Initialize Communication 

   TOMP_MAP(target enter data map(to: wtiso))

   do cSetID=nCommSets+1,nCommSets+nGTASets

     CSet => getCommSetData(Quad, cSetID)

!    Restore the initial communication order
     call restoreCommOrder(CSet)

!    Post receives for all data
     call InitExchange(cSetID)

   enddo

!  Initialize angle-independent variables before the sweeps

   if ( withSource ) then

#ifdef TETON_ENABLE_OPENACC
     !$acc parallel loop gang num_gangs(nZoneSets) &
     !$acc& vector_length(omp_device_team_thread_limit)
#else
     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
     TOMPC(shared(nZoneSets, Geom, GTA, wtiso, P))
#endif

     ZoneSetLoop: do zSetID=1,nZoneSets

#ifdef TETON_ENABLE_OPENACC
       !$acc loop vector
#else
       !$omp  parallel do default(none)  &
       !$omp& shared(Geom, GTA, zSetID, wtiso, P)
#endif
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         GTA% PhiInc(c)    = GTA% Sscat(c) 
         GTA% Q(c)         = wtiso*GTA%GreySigtInv(c)*( GTA% GreySource(c) +  &
                                   GTA% GreySigScat(c)*P(c) )
         GTA% TsaSource(c) = wtiso*Geom% Volume(c)*( GTA% GreySource(c) +  &
                                   GTA% GreySigScat(c)*P(c) )
       enddo
#ifndef TETON_ENABLE_OPENACC
       !$omp end parallel do
#endif

     enddo ZoneSetLoop

#ifdef TETON_ENABLE_OPENACC
     !$acc end parallel loop
#else
     TOMP(end target teams distribute)
#endif

   else

#ifdef TETON_ENABLE_OPENACC
     !$acc parallel loop gang num_gangs(nZoneSets) &
     !$acc& vector_length(omp_device_team_thread_limit)
#else
     TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
     TOMPC(shared(nZoneSets, Geom, GTA, wtiso, P))
#endif

     ZoneSetLoop2: do zSetID=1,nZoneSets

#ifdef TETON_ENABLE_OPENACC
       !$acc loop vector
#else
       !$omp  parallel do default(none)  &
       !$omp& shared(Geom, GTA, zSetID, wtiso, P)
#endif
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         GTA% PhiInc(c)    = zero
         GTA% Q(c)         = wtiso*GTA%GreySigtInv(c)*GTA% GreySigScat(c)*P(c)
         GTA% TsaSource(c) = wtiso*Geom% Volume(c)*GTA% GreySigScat(c)*P(c)
       enddo
#ifndef TETON_ENABLE_OPENACC
       !$omp end parallel do
#endif

     enddo ZoneSetLoop2

#ifdef TETON_ENABLE_OPENACC
     !$acc end parallel loop
#else
     TOMP(end target teams distribute)
#endif

   endif

TOMP_MAP(target exit data map(release: wtiso))

!  Loop over angles, solving for each in turn:

   AngleLoop: do sendIndex=1,Set% NumAnglesDyn

!$omp  parallel do default(none) schedule(static)  &
!$omp& private(CSet, Angle)  &
!$omp& shared(Quad, sendIndex, nCommSets, nGTASets, SnSweep,  PsiB) 
     do cSetID=nCommSets+1,nCommSets+nGTASets

       CSet  => getCommSetData(Quad, cSetID)
       Angle =  CSet% AngleOrder(sendIndex)

!      Send the boundary information needed by my neighbors
       call SendFlux(SnSweep, cSetID, sendIndex, PsiB)

!      Test for completion of the sends needed by my neighbors
       call TestSend(cSetID, sendIndex)

!      Receive the boundary information needed to compute this angle
       call RecvFlux(SnSweep, cSetID, Angle, PsiB)

!      Update incident fluxes on reflecting boundaries
       do setID=CSet% set1,CSet%set2
         call snreflect(SnSweep, setID, Angle, PsiB)
       enddo

     enddo
!$omp end parallel do


!    Sweep the mesh, calculating PSI for each corner; the 
!    boundary flux array PSIB is also updated here. 
!    Mesh cycles are fixed automatically.

     if (ndim == 3) then
       call SweepGreyUCBxyzNEW_GPU(sendIndex, PsiB)
     elseif (ndim == 2) then
       call SweepGreyUCBrzNEW_GPU(sendIndex, PsiB)
     endif

   enddo AngleLoop


   return
   end subroutine GTASweep_GPU 


