#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  03/2023, PFN                    *
!                                                                      *
!  CornerSweepUCBxyz_GPU:                                              *
!                                                                      * 
!                    This routine calculates angular fluxes for a      *
!                    single direction and multiple energy groups for   *
!                    for an upstream corner-balance (UCB) spatial      *
!                    in xyz-geometry.                                  *
!                                                                      *
!                    The work is offloaded to a GPU in which each      *
!                    computational "set" is a GPU block. The threads   *
!                    assigned to the block compute one group in one    *
!                    corner in a hyperplane (by definition, all of     *
!                    the corners in a hyperplane are independent).     *
!                                                                      *
!***********************************************************************

   subroutine CornerSweepUCBxyz_GPU(nSets, sendIndex, savePsi)

   use, intrinsic :: iso_c_binding, only : c_int
   use cmake_defines_mod, only : omp_device_team_thread_limit
   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use SetData_mod
   use AngleSet_mod
   use GroupSet_mod
   use CodeChecks_mod

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

   integer    :: setID
   integer    :: zSetID
   integer    :: Angle
   integer    :: g
   integer    :: Groups

   integer    :: mCycle
   integer    :: offSet
   integer    :: nAngleSets
   integer    :: nZoneSets
   integer    :: nHyperDomains

   integer    :: nzones
   integer    :: ii
   integer    :: ndone
   integer    :: hyperPlane
   integer    :: domID
   integer    :: hplane1
   integer    :: hplane2

   real(adqt) :: tau

!  Local

   integer    :: b 
   integer    :: i 
   integer    :: cfp 
   integer    :: cface 
   integer    :: ifp
   integer    :: c
   integer    :: c0 
   integer    :: cez
   integer    :: nCFaces

   integer    :: zone
   integer    :: nCorner
   integer    :: c1
   integer    :: c2

   real(adqt), parameter :: fouralpha=1.82d0

   real(adqt) :: aez
   real(adqt) :: area_opp

   real(adqt) :: source
   real(adqt) :: sig
   real(adqt) :: vol
   real(adqt) :: sigv
   real(adqt) :: sigv2
   real(adqt) :: gnum
   real(adqt) :: gden
   real(adqt) :: sez
   real(adqt) :: psi_opp
   real(adqt) :: afp
   real(adqt) :: Qez
   real(adqt) :: QQ
   real(adqt) :: SS
   real(adqt) :: mult

!  Dynamic

   integer, allocatable :: angleList(:)

!  Constants

   tau           = Size% tau
   nAngleSets    = getNumberOfAngleSets(Quad)
   nZoneSets     = getNumberOfZoneSets(Quad)
   nHyperDomains = getNumberOfHyperDomains(Quad,1)

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

!  Initialize

   ! Verify we won't get out-of-bounds accesses below.
   TETON_CHECK_BOUNDS1(Quad%AngSetPtr, nAngleSets)
   TETON_CHECK_BOUNDS1(Geom%corner1, nZoneSets)
   TETON_CHECK_BOUNDS1(Geom%corner2, nZoneSets)


   TOMP_MAP(target enter data map(to: tau, sendIndex, angleList))

   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(private(ASet, setID, Angle) &)
   TOMPC(shared(nZoneSets, angleList, Quad, Geom, nAngleSets) )

   ZoneSetLoop0: do zSetID=1,nZoneSets

!    Loop over angle sets

     do setID=1,nAngleSets

       ASet  => Quad% AngSetPtr(setID)
       angle =  angleList(setID)

! NOTE: This loop doesn't support a collapse(2), as its not a canonical form
! loop ( the inner loop bounds can not be predetermined ); it's significantly
! faster to split into two loops as below

       !$omp  parallel do collapse(2) default(none) &
       !$omp& shared(Geom, ASet, Angle, zSetID)

       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         do cface=1,3
           ASet% AfpNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_fp(:,cface,c) )
           ASet% AezNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_ez(:,cface,c) )
         enddo
       enddo

       !$omp end parallel do

       !$omp  parallel do default(none)  &
       !$omp& shared(Geom, ASet, Angle, zSetID)

       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         do cface=4,Geom% nCFacesArray(c)
           ASet% AfpNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_fp(:,cface,c) )
           ASet% AezNorm(cface,c) = DOT_PRODUCT( ASet% omega(:,angle),Geom% A_ez(:,cface,c) )
         enddo
       enddo

       !$omp end parallel do

     enddo

   enddo ZoneSetLoop0

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nZoneSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(private(ASet, setID) &)
   TOMPC(shared(nZoneSets, nAngleSets, Quad, Geom))

   ZoneSetLoop1: do zSetID=1,nZoneSets

!    Loop over angle sets

     do setID=1,nAngleSets

       ASet  => Quad% AngSetPtr(setID)

       !$omp  parallel do default(none)  &
       !$omp& shared(Geom, ASet, zSetID)

       do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
         ASet% ANormSum(c) = zero
         do cface=1,Geom% nCFacesArray(c)
           ASet% ANormSum(c) = ASet% ANormSum(c) + half*   &
                              (ASet% AfpNorm(cface,c) + abs( ASet% AfpNorm(cface,c) ) +  & 
                               ASet% AezNorm(cface,c) + abs( ASet% AezNorm(cface,c) ) )
         enddo
       enddo

       !$omp end parallel do

     enddo

   enddo ZoneSetLoop1

   TOMP(end target teams distribute)


   TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
   TOMPC(shared(sendIndex, Quad, nSets) &)
   TOMPC(private(Set, ASet, HypPlanePtr, angle, Groups, offSet, c))

   SetLoop0: do setID=1,nSets

     Set    => Quad% SetDataPtr(setID)
     ASet   => Quad% AngSetPtr(Set% angleSetID)

     Groups =  Set% Groups
     angle  =  Set% AngleOrder(sendIndex)
     offSet =  ASet% cycleOffSet(angle)

!  Initialize boundary values in Psi1 and interior values on the cycle
!  list

     !$omp  parallel do collapse(2) default(none) &
     !$omp& shared(angle, Set, ASet, offSet, Groups) private(c)

     do mCycle=1,ASet% numCycles(angle)
       do g=1,Groups
         c              = ASet% cycleList(offSet+mCycle)
         Set% Psi1(g,c) = Set% cyclePsi(g,offSet+mCycle)
       enddo
     enddo

     !$omp end parallel do


     !$omp  parallel do collapse(2) default(none) &
     !$omp& shared(Set, Groups, angle)

     do b=1,Set%nbelem
       do g=1,Groups
         Set% Psi1(g,Set%nCorner+b) = Set% PsiB(g,b,angle)
       enddo
     enddo

     !$omp end parallel do

   enddo SetLoop0

   TOMP(end target teams distribute)

   if ( nHyperDomains > 1 ) then

     TOMP(target teams distribute collapse(2) num_teams(nZoneSets*nSets) &)
     TOMPC(thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(sendIndex, Quad, nZoneSets, nSets) &)
     TOMPC(private(Set, ASet, HypPlanePtr, angle, Groups, c))

     ZoneSetLoop2: do zSetID=1,nZoneSets
       SetLoop2: do setID=1,nSets

         Set         => Quad% SetDataPtr(setID)
         ASet        => Quad% AngSetPtr(Set% angleSetID)
         Groups      =  Set% Groups
         angle       =  Set% AngleOrder(sendIndex)
         HypPlanePtr => ASet% HypPlanePtr(angle)

!        Initialize values at hyper-domain interfaces

         !$omp  parallel do collapse(2) default(none) &
         !$omp& shared(Set, HypPlanePtr, zSetID, Groups, angle) private(c)

         do i=HypPlanePtr% c1(zSetID),HypPlanePtr% c2(zSetID)
           do g=1,Groups
             c              = HypPlanePtr% interfaceList(i)
             Set% Psi1(g,c) = Set% PsiInt(g,i,angle)
           enddo
         enddo
         !$omp end parallel do

       enddo SetLoop2
     enddo ZoneSetLoop2

     TOMP(end target teams distribute)

   endif


! TODO:
! IBM XLF segfaults if 'hyperPlane' is not scoped to private below.
! This should not be necessary, as this is a loop control variables which the runtime should automatically scope to
! private.
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

   TOMP(target teams distribute collapse(2) num_teams(nSets*nHyperDomains) &)
   TOMPC(thread_limit(omp_device_team_thread_limit) &)
   TOMPC(shared(Quad, sendIndex, nSets, nHyperDomains) &)
   TOMPC(private(Set, ASet, GSet, HypPlanePtr, Angle, Groups) &)
   TOMPC(private(hplane1, hplane2, ndone, nCorner, hyperPlane)) 

   SetLoop: do setID=1,nSets
     DomainLoop: do domID=1,nHyperDomains

     Set          => Quad% SetDataPtr(setID)
     ASet         => Quad% AngSetPtr(Set% angleSetID) 
     GSet         => Quad% GrpSetPtr(Set% groupSetID) 

     Groups       =  Set% Groups
     Angle        =  Set% AngleOrder(sendIndex)
     HypPlanePtr  => ASet% HypPlanePtr(Angle)
     hplane1      =  HypPlanePtr% hplane1(domID)
     hplane2      =  HypPlanePtr% hplane2(domID)
     ndone        =  HypPlanePtr% ndone(domID) 

     HyperPlaneLoop: do hyperPlane=hplane1,hplane2

       nCorner = HypPlanePtr% cornersInPlane(hyperPlane)

       !$omp  parallel do collapse(2) default(none) &
       !$omp& shared(Set, Geom, ASet, GSet, Angle, nCorner, Groups, ndone, tau) &
       !$omp& private(c, c0, cfp, ifp, cez, zone, cface, i, nCFaces) &
       !$omp& private(aez, area_opp, sig, vol, Qez, QQ, SS, source) &
       !$omp& private(sigv, sigv2, sez, gnum, gden, psi_opp, afp, c1, c2, mult)

       CornerLoop1: do ii=1,nCorner
         GroupLoop1: do g=1,Groups

!          Loop through the corners using the NEXT list

           c    = ASet% nextC(ndone+ii,Angle)
           zone = Geom% CToZone(c)
           c0   = Geom% cOffSet(zone)
           sig  = GSet% Sigt(g,zone)

!          Calculate Area_CornerFace dot Omega to determine the
!          contributions from incident fluxes across external
!          corner faces (FP faces)

           nCFaces = Geom% nCFacesArray(c)

!          Contributions from volume terms 

           source = GSet% STotal(g,c) + tau*Set% Psi(g,c,Angle)
           SS     = Geom% Volume(c)*source

           do cface=1,nCFaces

             afp = ASet% AfpNorm(cface,c)
             cfp = Geom% cFP(cface,c)

             if ( afp < zero ) then
               SS = SS - afp*Set% Psi1(g,cfp)
             endif
           enddo

!          Contributions from interior corner faces (EZ faces)

           do cface=1,nCFaces

             aez  = ASet% AezNorm(cface,c)
             mult = zero 

             if (aez > zero ) then

               c1   = c
               c2   = c0 + Geom% cEZ(cface,c)
               ifp  = mod(cface,nCFaces) + 1
               QQ   = source
               Qez  = GSet% STotal(g,c2) + tau*Set% Psi(g,c2,Angle)
               mult = one

             elseif (aez < zero ) then

               c2   = c
               c1   = c0 + Geom% cEZ(cface,c)

!              Contributions from upsteam fluxes in the same zone
               SS = SS - aez*Set% Psi1(g,c1)

               do i=1,nCFaces
                 if (Geom% cEZ(i,c1) + c0 == c2) then
                   ifp  = mod(i,nCFaces) + 1
                 endif
               enddo

               Qez  = source
               QQ   = GSet% STotal(g,c1) + tau*Set% Psi(g,c1,Angle)
               mult = -one
               aez  = -aez

             endif

             psi_opp  = zero
             area_opp = zero

             afp = ASet% AfpNorm(ifp,c1)

             if ( afp < zero ) then
               cfp      =  Geom% cFP(ifp,c1)
               area_opp = -afp
               psi_opp  = -afp*Set% Psi1(g,cfp)
             endif

             do i=2,nCFaces-2
               ifp = mod(ifp,nCFaces) + 1
               afp = ASet% AfpNorm(ifp,c1)
               if ( afp < zero ) then
                 cfp      = Geom% cFP(ifp,c1)
                 area_opp = area_opp - afp
                 psi_opp  = psi_opp  - afp*Set% Psi1(g,cfp)
               endif
             enddo

             TestOppositeFace: if ( area_opp > zero ) then

               psi_opp  = psi_opp/area_opp

               vol      = Geom% Volume(c1)

               sigv     = sig*vol
               sigv2    = sigv*sigv

               gnum     = aez*aez*( fouralpha*sigv2 +              &
                          aez*(four*sigv + three*aez) )

               gden     = vol*( four*sigv*sigv2 + aez*(six*sigv2 + &
                          two*aez*(two*sigv + aez)) )

               sez      = ( vol*gnum*( sig*psi_opp - QQ ) +   &
                            half*aez*gden*( QQ - Qez ) )/ &
                          ( gnum + gden*sig)

               SS       = SS   + mult*sez

             else

               sez  = half*aez*( QQ - Qez )/sig
               SS   = SS   + mult*sez

             endif TestOppositeFace

           enddo

!          Corner angular flux
           Set% Psi1(g,c) = SS/(ASet% ANormSum(c) + sig*Geom% Volume(c))

         enddo GroupLoop1
       enddo CornerLoop1

       !$omp end parallel do

       ndone = ndone + nCorner

     enddo HyperPlaneLoop

     enddo DomainLoop
   enddo SetLoop

   TOMP(end target teams distribute)

!  Update Boundary data

     TOMP(target teams distribute num_teams(nSets) thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nSets, Quad, sendIndex)&)
     TOMPC(private(Set, ASet, BdyExitPtr, HypPlanePtr, offSet, angle, Groups, b, c))

     SetLoop3: do setID=1,nSets

       Set         => Quad% SetDataPtr(setID)
       ASet        => Quad% AngSetPtr(Set% angleSetID)
       Groups      =  Set% Groups
       angle       =  Set% AngleOrder(sendIndex)
       offSet      =  ASet% cycleOffSet(angle)
       BdyExitPtr  => ASet% BdyExitPtr(angle)
       HypPlanePtr => ASet% HypPlanePtr(angle)

       !$omp  parallel do collapse(2) default(none) &
       !$omp& shared(Set, BdyExitPtr, Groups, Angle) private(b,c)

       do i=1,BdyExitPtr% nxBdy
         do g=1,Groups
           b = BdyExitPtr% bdyList(1,i)
           c = BdyExitPtr% bdyList(2,i)

           Set% PsiB(g,b,Angle) = Set% Psi1(g,c)
         enddo
       enddo

       !$omp end parallel do

!      Update Psi in the cycle list

       !$omp  parallel do collapse(2) default(none) &
       !$omp& shared(angle, Set, ASet, offSet, Groups) private(c)

       do mCycle=1,ASet% numCycles(angle)
         do g=1,Groups
           c                              = ASet% cycleList(offSet+mCycle)
           Set% cyclePsi(g,offSet+mCycle) = Set% Psi1(g,c)
         enddo
       enddo

       !$omp end parallel do

     enddo SetLoop3

     TOMP(end target teams distribute)

   if ( nHyperDomains > 1 ) then

     TOMP(target teams distribute collapse(2) num_teams(nZoneSets*nSets) &)
     TOMPC(thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(sendIndex, Quad, nZoneSets, nSets) &)
     TOMPC(private(Set, ASet, HypPlanePtr, angle, Groups, c))

     ZoneSetLoop4: do zSetID=1,nZoneSets
       SetLoop4: do setID=1,nSets

         Set         => Quad% SetDataPtr(setID)
         ASet        => Quad% AngSetPtr(Set% angleSetID)
         Groups      =  Set% Groups
         angle       =  Set% AngleOrder(sendIndex)
         HypPlanePtr => ASet% HypPlanePtr(angle)

!        Update values at hyper-domain interfaces

         !$omp  parallel do collapse(2) default(none) &
         !$omp& shared(Set, HypPlanePtr, Groups, zSetID, angle) private(c)

         do i=HypPlanePtr% c1(zSetID),HypPlanePtr% c2(zSetID)
           do g=1,Groups
             c                      = HypPlanePtr% interfaceList(i)
             Set% PsiInt(g,i,angle) = Set% Psi1(g,c)
           enddo
         enddo

         !$omp end parallel do

       enddo SetLoop4
     enddo ZoneSetLoop4

     TOMP(end target teams distribute)

   endif


!  We only store Psi if this is the last transport sweep in the time step

   if ( savePsi ) then

     TOMP(target teams distribute collapse(2) num_teams(nZoneSets*nSets) &)
     TOMPC(thread_limit(omp_device_team_thread_limit) default(none) &)
     TOMPC(shared(nZoneSets, nSets, Quad, Geom, sendIndex) &)
     TOMPC(private(Set, setID, Angle, Groups))

     ZoneSetLoop5: do zSetID=1,nZoneSets
       SetLoop5: do setID=1,nSets

         Set    => Quad% SetDataPtr(setID)
         Groups =  Set% Groups
         Angle  =  Set% AngleOrder(sendIndex)

         !$omp  parallel do collapse(2) default(none) &
         !$omp& shared(Set, Angle, Groups, zSetID, Geom)

         CornerLoop5: do c=Geom% corner1(zSetID),Geom% corner2(zSetID)
           GroupLoop5: do g=1,Groups

             Set% Psi(g,c,Angle) = Set% Psi1(g,c)

           enddo GroupLoop5
         enddo CornerLoop5

         !$omp end parallel do

       enddo SetLoop5
     enddo ZoneSetLoop5

     TOMP(end target teams distribute)

   endif

   TOMP_MAP(target exit data map(release: tau, sendIndex, angleList))


   deallocate( angleList )


   return
   end subroutine CornerSweepUCBxyz_GPU

