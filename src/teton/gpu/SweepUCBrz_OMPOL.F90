#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   SweepUCBrz_GPU - This routine calculates angular fluxes for a      *
!                    single direction and multiple energy groups for   *
!                    for an upstream corner-balance (UCB) spatial      *
!                    in rz-geometry.                                   *
!                                                                      *
!                    The work is offloaded to a GPU in which each      *
!                    computational "set" is a GPU block. The threads   *
!                    assigned to the block compute one group in one    *
!                    zone in a hyperplane (by definition, all of the   *
!                    zones in a hyperplane are independent).           *
!                                                                      *
!***********************************************************************

   subroutine SweepUCBrz_GPU(nSets, sendIndex, savePsi)

   use, intrinsic :: iso_c_binding, only : c_int
   use cmake_defines_mod, only : omp_device_team_thread_limit
   use Options_mod
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use GroupSet_mod
   use ArrayChecks_mod

   implicit none

!  Arguments

   integer,          intent(in)    :: nSets
   integer,          intent(in)    :: sendIndex
   logical (kind=1), intent(in)    :: savePsi

!  Local

   type(SetData),    pointer       :: Set
   type(AngleSet),   pointer       :: ASet
   type(GroupSet),   pointer       :: GSet
   type(HypPlane),   pointer       :: HypPlanePtr
   type(BdyExit),    pointer       :: BdyExitPtr

   integer            :: setID
   integer            :: zSetID
   integer            :: Angle
   integer            :: g
   integer            :: Groups

   integer            :: mCycle
   integer            :: offset
   integer            :: nAngleSets
   integer            :: nZoneSets
   integer            :: nzones
   integer            :: c
   integer            :: ii
   integer            :: ndoneZ
   integer            :: hyperPlane
   integer            :: nHyperplanes

   real(adqt)         :: tau

!  Local

   integer    :: b
   integer    :: i
   integer    :: cface
   integer    :: cez
   integer    :: cfp
   integer    :: c0

   integer    :: zone
   integer    :: nCorner

   real(adqt), parameter :: fouralpha=1.82d0
   real(adqt) :: fac
   real(adqt) :: sigA
   real(adqt) :: sigA2
   real(adqt) :: source

   real(adqt) :: area
   real(adqt) :: sig
   real(adqt) :: sez
   real(adqt) :: gnum
   real(adqt) :: gden
   real(adqt) :: quadTauW1
   real(adqt) :: quadTauW2

   real(adqt) :: afp
   real(adqt) :: aez
   real(adqt) :: R
   real(adqt) :: R2
   real(adqt) :: R_afp
   real(adqt) :: R_afp2
   real(adqt) :: denom

!  Dynamic

   integer, allocatable :: angleList(:)

!  Constants
   tau        = Size% tau
   nAngleSets = getNumberOfAngleSets(Quad)
   nZoneSets  = getNumberOfZoneSets(Quad)

   allocate( angleList(nAngleSets) )

   do setID=1,nSets
     Set => Quad% SetDataPtr(setID)
     angleList(Set% angleSetID) = Set% AngleOrder(sendIndex)
   enddo

!  Here the maximum block size is the product of the maximum
!  number of zones in a hyperplane and the number of groups;
!  The maximum value is the over all teams  

!  Note: num_blocks = nSets and the number of threads per 
!  team (a.k.a. "block") <= block_threads

   ! Verify we won't get out-of-bounds accesses below.
   TETON_CHECK_BOUNDS1(Quad%AngSetPtr, nAngleSets)
   TETON_CHECK_BOUNDS1(Geom%corner1, nZoneSets)
   TETON_CHECK_BOUNDS1(Geom%corner2, nZoneSets)

   TOMP(target data map(to: tau, sendIndex, angleList))


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(nZoneSets, nAngleSets,Geom, angleList, Quad)&)
   TOMPC(private(ASet, Angle))

   ZoneSetLoop: do zSetID=1,nZoneSets

!    Loop over angle sets

     do setID=1,nAngleSets

       ASet  => Quad% AngSetPtr(setID)
       angle =  angleList(setID)

!$omp  parallel do collapse(2) default(none)  &
!$omp& shared(Geom, ASet, Angle, zSetID) private(c,cface)

       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         do cface=1,2
           ASet% AfpNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_fp(:,cface,c) )
           ASet% AezNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_ez(:,cface,c) )
         enddo
       enddo

!$omp end parallel do

     enddo

   enddo ZoneSetLoop

TOMP(end target teams distribute)

TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none)&)
TOMPC(shared(nZoneSets, nAngleSets, Geom, angleList, Quad)&)
TOMPC(private(ASet, angle, fac, R_afp,R_afp2,R,R2))

   ZoneSetLoop2: do zSetID=1,nZoneSets

!    Loop over angle sets

     do setID=1,nAngleSets

       ASet  => Quad% AngSetPtr(setID)
       angle =  angleList(setID)
       fac   =  ASet% angDerivFac(Angle)

!$omp  parallel do default(none)  &
!$omp& shared(Geom, ASet, angle, fac, zSetID) private(R_afp,R_afp2,R,R2)

       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         R_afp  = Geom% RadiusFP(1,c)
         R_afp2 = Geom% RadiusFP(2,c)
         R      = Geom% RadiusEZ(1,c)
         R2     = Geom% RadiusEZ(2,c)

         ASet% ANormSum(c) = fac*Geom% Area(c) - half*(  &
            R_afp *(ASet% AfpNorm(1,c) - abs(ASet% AfpNorm(1,c))) +    &
            R_afp2*(ASet% AfpNorm(2,c) - abs(ASet% AfpNorm(2,c))) +    &
            R     *(ASet% AezNorm(1,c) - abs(ASet% AezNorm(1,c))) +    &
            R2    *(ASet% AezNorm(2,c) - abs(ASet% AezNorm(2,c))) ) 
       enddo

!$omp end parallel do

     enddo

   enddo ZoneSetLoop2

TOMP(end target teams distribute)

! TODO:
! IBM XLF segfaults if 'mCycle', 'b', and 'g' are not scoped to private below.  This should not
! be necessary, as these are loop control variables which the runtime should automatically scope to private.
!
! Relevant portions of OpenMP spec:
! `The loop iteration variable in any associated loop of a for, parallel for,
! taskloop, or distribute construct is private.`
!
! `A loop iteration variable for a sequential loop in a parallel or task
! generating construct is private in the innermost such construct that encloses
! the loop.`
! 
! Look into reporting this bug to IBM, using UMT as a reproducer.
TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
TOMPC(shared(nSets, Quad, Geom, sendIndex, tau)&)
TOMPC(private(c, cfp, Set, ASet, GSet, HypPlanePtr, Angle) &)
TOMPC(private(Groups, nHyperPlanes, ndoneZ, mCycle, b, g, offset, hyperPlane, nzones, fac)&)
TOMPC(private(c0,cez,zone,nCorner, sigA,sigA2,source,area,sig,sez,gnum,gden, aez,afp,R,R_afp,denom))

   SetLoop: do setID=1,nSets

     Set          => Quad% SetDataPtr(setID)
     ASet         => Quad% AngSetPtr(Set% angleSetID)
     GSet         => Quad% GrpSetPtr(Set% groupSetID)

     Groups       =  Set% Groups
     Angle        =  Set% AngleOrder(sendIndex)
     nHyperPlanes =  ASet% nHyperPlanes(Angle)
     ndoneZ       =  0

!    Angle Constants

     fac          = ASet% angDerivFac(Angle)

!    Initialize variable.
     cfp          = -1

!  Initialize boundary values in Psi1 and interior values on the cycle
!  list

     HypPlanePtr => ASet% HypPlanePtr(Angle)
     offSet      =  ASet% cycleOffSet(angle)

     do mCycle=1,ASet% numCycles(Angle)
       c = ASet% cycleList(offSet+mCycle)

       do g=1,Groups
         Set% Psi1(g,c) = Set% cyclePsi(g,offSet+mCycle)
       enddo
     enddo

!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, Groups, Angle)
     do b=1,Set%nbelem
       do g=1,Groups
         Set% Psi1(g,Set%nCorner+b) = Set% PsiB(g,b,Angle)
       enddo
     enddo
!$omp end parallel do

     HyperPlaneLoop: do hyperPlane=1,nHyperPlanes

       nzones = HypPlanePtr% zonesInPlane(hyperPlane)

!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, Geom, ASet, GSet, Angle, nzones, Groups, ndoneZ, tau, fac) &
!$omp& private(c0,cez,cfp,zone,nCorner,sigA,sigA2,source,area,sig,sez,gnum,gden) &
!$omp& private(aez,afp,R,R_afp,denom)

       ZoneLoop: do ii=1,nzones
         GroupLoop: do g=1,Groups

!  Loop through all of the corners using the NEXT list

           zone    = ASet% nextZ(ndoneZ+ii,Angle)
           nCorner = Geom% numCorner(zone)
           c0      = Geom% cOffSet(zone)
           sig     = GSet% Sigt(g,zone)

!  Contributions from volume terms (if a starting direction add angular
!  derivative)

           do c=1,nCorner
             source         = GSet% STotal(g,c0+c) + tau*Set% Psi(g,c0+c,Angle)
             Set% Q(g,c,ii) = source
             Set% S(g,c,ii) = Geom% Volume(c0+c)*source +    &
                              fac*Geom% Area(c0+c)*Set% PsiM(g,c0+c)
           enddo

           CornerLoop: do c=1,nCorner

             CornerFaceLoop: do cface=1,2

               afp = ASet% AfpNorm(cface,c0+c) 
               aez = ASet% AezNorm(cface,c0+c) 

               if ( afp < zero ) then
                 cfp            = Geom% cFP(cface,c0+c)
                 R_afp          = Geom% RadiusFP(cface,c0+c)*afp
                 Set% S(g,c,ii) = Set% S(g,c,ii) - R_afp*Set% Psi1(g,cfp)
               endif

               if ( aez > zero ) then

                 R   = Geom% RadiusEZ(cface,c0+c)
                 cez = Geom% cEZ(cface,c0+c)

                 if ( afp < zero ) then

                   area     = Geom% Area(c0+c)
                   sigA     = sig*area
                   sigA2    = sigA*sigA

                   gnum     = aez*aez*( fouralpha*sigA2 +  &
                              aez*(four*sigA + three*aez) )

                   gden     = area*(four*sigA*sigA2 + aez*(six*sigA2 + &
                              two*aez*(two*sigA + aez)))

                   sez      = R*( area*gnum*  &
                              ( sig*Set% Psi1(g,cfp) - Set% Q(g,c,ii) ) + &
                              half*aez*gden*( Set% Q(g,c,ii) - Set% Q(g,cez,ii) ) )/  &
                              ( gnum + gden*sig )

                   Set% S(g,c,ii)   = Set% S(g,c,ii)   + sez
                   Set% S(g,cez,ii) = Set% S(g,cez,ii) - sez

                 else

                   sez              = half*R*aez*( Set% Q(g,c,ii) - Set% Q(g,cez,ii) )/sig
                   Set% S(g,c,ii)   = Set% S(g,c,ii)   + sez
                   Set% S(g,cez,ii) = Set% S(g,cez,ii) - sez

                 endif

               endif

             enddo CornerFaceLoop

           enddo CornerLoop

!  Solve the corners in the correct order

           do i=1,nCorner

             c = ASet% nextC(c0+i,angle)

!            Corner angular flux
             denom             = ASet% ANormSum(c0+c) + sig*Geom% Volume(c0+c)
             Set% Psi1(g,c0+c) = Set% S(g,c,ii)/denom

!            Calculate the contribution of this flux to the sources of
!            downstream corners in this zone. The downstream corner index is
!            "ez_exit."

             do cface=1,2
               aez = ASet% AezNorm(cface,c0+c)

               if (aez > zero) then
                 R                = Geom% RadiusEZ(cface,c0+c)
                 cez              = Geom% cEZ(cface,c0+c) 
                 Set% S(g,cez,ii) = Set% S(g,cez,ii) + R*aez*Set% Psi1(g,c0+c)
               endif
             enddo

           enddo

         enddo GroupLoop
       enddo ZoneLoop

!$omp end parallel do 

       ndoneZ = ndoneZ + nzones

     enddo HyperPlaneLoop

!    Update Psi in the cycle list

     do mCycle=1,ASet% numCycles(angle)
       c = ASet% cycleList(offSet+mCycle)

       do g=1,Groups
         Set% cyclePsi(g,offSet+mCycle) = Set% Psi1(g,c)
       enddo
     enddo
   enddo SetLoop

TOMP(end target teams distribute)
TOMP(end target data)

TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
TOMPC(shared(nSets, Quad, sendIndex)&)
TOMPC(private(Set, ASet, Angle, Groups, quadTauW1, quadTauW2))

     SetLoop2: do setID=1,nSets

       Set    => Quad% SetDataPtr(setID)
       ASet   => Quad% AngSetPtr(Set% angleSetID)

       Groups =  Set% Groups
       Angle  =  Set% AngleOrder(sendIndex)

!      Set the "half-angle" angular intensity (PsiM) for the next angle

       if ( ASet% StartingDirection(Angle) ) then

!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, Groups)

         do c=1,Set% nCorner
           do g=1,Groups
             Set% PsiM(g,c) = Set% Psi1(g,c)
           enddo 
         enddo 

!$omp end parallel do

       else

         quadTauW1 = ASet% quadTauW1(Angle)
         quadTauW2 = ASet% quadTauW2(Angle)

!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, Groups, quadTauW1, quadTauW2)

         do c=1,Set% nCorner
           do g=1,Groups
             Set% PsiM(g,c) = quadTauW1*Set% Psi1(g,c) -  &
                              quadTauW2*Set% PsiM(g,c)
           enddo 
         enddo

!$omp end parallel do

       endif

     enddo SetLoop2

TOMP(end target teams distribute)
TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none)&)
TOMPC(shared(nSets, Quad, sendIndex)&)
TOMPC(private(Set, ASet, BdyExitPtr, Angle, Groups, b, c))

!!TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) private(setID, Set, ASet, BdyExitPtr, Angle, Groups))

     SetLoop3: do setID=1,nSets

       Set        => Quad% SetDataPtr(setID)
       ASet       => Quad% AngSetPtr(Set% angleSetID)
       Groups     =  Set% Groups
       Angle      =  Set% AngleOrder(sendIndex)
       BdyExitPtr => ASet% BdyExitPtr(Angle)

!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, BdyExitPtr, Angle, Groups) private(b,c)

       do i=1,BdyExitPtr% nxBdy
         do g=1,Groups
           b = BdyExitPtr% bdyList(1,i)
           c = BdyExitPtr% bdyList(2,i)

           Set% PsiB(g,b,Angle) = Set% Psi1(g,c)
         enddo
       enddo

!$omp end parallel do

       if ( ASet% FinishingDirection(Angle+1) ) then

!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, BdyExitPtr, Angle, Groups) private(b,c)

         do i=1,BdyExitPtr% nxBdy
           do g=1,Groups
             b = BdyExitPtr% bdyList(1,i)
             c = BdyExitPtr% bdyList(2,i)

             Set% PsiB(g,b,Angle+1) = Set% PsiM(g,c)
           enddo
         enddo

!$omp end parallel do

       endif 

     enddo SetLoop3

TOMP(end target teams distribute)


!  We only store Psi if this is the last transport sweep in the time step

   if ( savePsi ) then

TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
TOMPC(shared(nSets, sendIndex, Quad)&)
TOMPC(private(Set, ASet, Angle, Groups))

     SetLoop4: do setID=1,nSets

       Set    => Quad% SetDataPtr(setID)
       ASet   => Quad% AngSetPtr(Set% angleSetID)

       Groups =  Set% Groups
       Angle  =  Set% AngleOrder(sendIndex)

!$omp  parallel do collapse(2) default(none) &
!$omp& shared(Set, ASet, Angle, Groups)

       CornerLoop4: do c=1,Set% nCorner
         GroupLoop4: do g=1,Groups

           Set% Psi(g,c,Angle) = Set% Psi1(g,c)

           if ( ASet% FinishingDirection(Angle+1) ) then
             Set% Psi(g,c,Angle+1) = Set% PsiM(g,c)
           endif

         enddo GroupLoop4
       enddo CornerLoop4

!$omp end parallel do

     enddo SetLoop4

TOMP(end target teams distribute)

   endif

   deallocate( angleList )


   return
   end subroutine SweepUCBrz_GPU

!***********************************************************************
!                        Last Update:  01/2018, PFN                    *
!                                                                      *
!   SweepUCBrz_kernel  - This routine calculates angular fluxes for a  *
!                        single direction, energy group and zone for   *
!                        for an upstream corner-balance (UCB) spatial  *
!                        in rz-geometry.                               *
!                                                                      *
!***********************************************************************
