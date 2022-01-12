#include "macros.h"
#include "omp_wrappers.h"
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

   subroutine GTASweep_GPU(PsiB)

   use, intrinsic :: iso_c_binding, only : c_int
   use Options_mod
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

   real(adqt),           intent(inout) :: PsiB(Size%nbelem,Size%nangGTA)

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

   integer(kind=c_int) :: nOmpMaxTeamThreads

   logical (kind=1)          :: SnSweep 

!  Constants
!
   nOmpMaxTeamThreads = Options%getNumOmpMaxTeamThreads()
   nSets     = getNumberOfSets(Quad)
   nGTASets  = getNumberOfGTASets(Quad)
   nZoneSets = getNumberOfZoneSets(Quad)
   nCommSets = getNumberOfCommSets(Quad)
   ndim      =  Size% ndim
   SnSweep   = .FALSE.

   Set       => getSetData(Quad, nSets+1)

!  Update TSA source on the GPU

   TOMP(target update to(GTA%TsaSource) )

   do cSetID=nCommSets+1,nCommSets+nGTASets

     CSet => getCommSetData(Quad, cSetID)

!    Restore the initial communication order
     call restoreCommOrder(CSet)

!    Post receives for all data
     call InitExchange(cSetID)

   enddo

!  Initialize angle-independent variables before the sweeps

TOMP(target teams distribute num_teams(nZoneSets) thread_limit(nOmpMaxTeamThreads) private(zSetID))

     ZoneSetLoop: do zSetID=1,nZoneSets

!$omp  parallel do default(none)  &
!$omp& shared(Geom, GTA, zSetID) private(c)
       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         GTA% PhiInc(c)    = zero
         GTA% Q(c)         = GTA%GreySigtInv(c)*GTA% TsaSource(c)
         GTA% TsaSource(c) = Geom% Volume(c)*GTA% TsaSource(c)
       enddo
!$omp end parallel do

     enddo ZoneSetLoop

TOMP(end target teams distribute)

!  Loop over angles, solving for each in turn:

   AngleLoop: do sendIndex=1,Set% NumAnglesDyn

!$omp parallel do private(cSetID, setID, CSet, Angle)  &
!$omp& shared(PsiB) schedule(static)
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


