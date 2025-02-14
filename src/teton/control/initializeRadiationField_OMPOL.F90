#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   ADVANCERT - Save zone-average quantities from previous cycle for   *
!               delta-t calculation.  Convert specific radiation       *
!               intensity (i.e. per unit mass) to intensity (per       *
!               unit volume) before the transport calculation.         *
!                                                                      *
!***********************************************************************
 
   subroutine initializeRadiationField 

   use, intrinsic :: iso_c_binding, only : c_int
   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod

   implicit none

!  Local

   type(SetData),  pointer  :: Set
   type(AngleSet), pointer  :: ASet
   type(BdyExit),  pointer  :: BdyExitPtr
   type(HypPlane), pointer  :: HypPlanePtr

   integer    :: setID
   integer    :: nSets
   integer    :: angle
   integer    :: b
   integer    :: c
   integer    :: NumAngles
   integer    :: i
   integer    :: Groups
   integer    :: g
   integer    :: mCycle
   integer    :: offSet

!  Constants

   nSets = getNumberOfSets(Quad)

   if (Size% ndim == 1) then
     return
   endif

!***********************************************************************
!  Scale the radiation field to account for volume changes and         *
!  tally beginning-of-cycle radiation energy.                          *
!                                                                      *
!  Compute the work done on radiation field due to volume changes.     *
!  This is an external source rate in the radiation transport equation.*
!***********************************************************************

   if (Size%useGPU) then

!  Update Boundary data

     TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nSets, Quad)&)
     TOMPC(private(Set, ASet, BdyExitPtr, HypPlanePtr, offSet, Groups, NumAngles, c, b, angle))

     SetLoop: do setID=1,nSets

       Set        => Quad% SetDataPtr(setID)
       ASet       => Quad% AngSetPtr(Set% angleSetID)
       Groups     =  Set% Groups
       NumAngles  =  Set% NumAngles

       do angle=1,NumAngles
         BdyExitPtr => ASet% BdyExitPtr(angle)

         !$omp  parallel do collapse(2) default(none) &
         !$omp& shared(Set, BdyExitPtr, angle, Groups) &
         !$omp& private(b,c)

         do i=1,BdyExitPtr% nxBdy
           do g=1,Groups
             b = BdyExitPtr% bdyList(1,i)
             c = BdyExitPtr% bdyList(2,i)

             Set% PsiB(g,b,angle) = Set% Psi(g,c,angle)
           enddo
         enddo

        !$omp end parallel do

       enddo

!    Update Psi in the cycle list

       do angle=1,NumAngles
         offSet = ASet% cycleOffSet(angle)

         !$omp  parallel do collapse(2) default(none) &
         !$omp& shared(Set, ASet, angle, offSet, Groups) &
         !$omp& private(c)

         do mCycle=1,ASet% numCycles(angle)
           do g=1,Groups
             c                              = ASet% cycleList(offSet+mCycle)
             Set% cyclePsi(g,offSet+mCycle) = Set% Psi(g,c,angle)
           enddo
         enddo

         !$omp end parallel do
         
       enddo

!      Update the Psi saved at hyper-domain interfaces

       do angle=1,NumAngles
         HypPlanePtr => ASet% HypPlanePtr(angle)

         !$omp  parallel do collapse(2) default(none) &
         !$omp& shared(Set, HypPlanePtr, angle, Groups) &
         !$omp& private(c)

         do i=1,HypPlanePtr% interfaceLen
           do g=1,Groups
             c = HypPlanePtr% interfaceList(i)
             Set% PsiInt(g,i,angle) = Set% Psi(g,c,angle)
           enddo
         enddo

         !$omp end parallel do
       enddo

     enddo SetLoop

     TOMP(end target teams distribute)

   else

     !$omp parallel do default(none) schedule(static) &
     !$omp& shared(Quad, nSets) &
     !$omp& private(Set, ASet, BdyExitPtr, NumAngles, b, c)
     do setID=1,nSets

       Set  => getSetData(Quad, setID)
       ASet => getAngleSetFromSetID(Quad, setID)

       NumAngles = Set% NumAngles

!    Initialize exiting boundary fluxes

       do angle=1,NumAngles
         BdyExitPtr => ASet% BdyExitPtr(angle)

         do i=1,BdyExitPtr% nxBdy
           b = BdyExitPtr% bdyList(1,i)
           c = BdyExitPtr% bdyList(2,i)

           Set% PsiB(:,b,angle) = Set% Psi(:,c,angle)
         enddo
       enddo

!    Initialize Psi in the cycle list
       call initCyclePsi(setID)

     enddo
     !$omp end parallel do

   endif ! if Size%useGPU
 
   return
   end subroutine initializeRadiationField 
