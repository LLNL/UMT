!***********************************************************************
!                        Last Update:  10/2016, PFN                    *
!                                                                      *
!   SweepUCBxyz  - TODO:  Fill in expl.
!
!                 This routine calculates angular fluxes for a         *
!                 single direction and multiple energy groups for      *
!                 for an upstream corner-balance (UCB) spatial         *
!                 in xyz-geometry.                                     *
!                                                                      *
!***********************************************************************

subroutine SweepUCBxyzToGPU (Set, setID, Groups, Angle, savePsi, Psi1) 

  use kind_mod
  use constant_mod
  use Size_mod
  use Quadrature_mod
  use QuadratureList_mod
  use BoundaryList_mod
  use Boundary_mod
  use SetData_mod
  use AngleSet_mod
  use GroupSet_mod
  use Geometry_mod
  use omp_lib

  implicit none

  interface 
     subroutine gpu_streamSynchronize( streamId ) bind (c)
       use iso_c_binding
       integer (c_int) :: streamId
     end subroutine gpu_streamSynchronize
  end interface

  interface 
     subroutine gpu_deviceSynchronize() bind (c)
       use iso_c_binding
     end subroutine gpu_deviceSynchronize
  end interface

  interface
     subroutine gpu_sweepucbxyz ( Angle, nHyperPlanes, nZonesInPlane, nextZ, nextC, STotal, tau, Psi, &
          Groups, Volume, Sigt, nCFacesArray, ndim, maxcf, &
          ncorner, A_fp, omega, cFP, Psi1, nbelem, A_ez, cEZ, NumAngles, quadwt, Phi, PsiB, maxCorner, mem0solve1, threadNum, numThreads, &
          savePsi, numCycles, cycleOffset, cyclePsi, cycleList, b0, nBdyElem, PsiBMref, Mref, Geo_numCorner, Geo_cOffSet) bind (c)
       use iso_c_binding
       integer (c_int) :: Angle
       integer (c_int) :: nHyperPlanes
       integer (c_int) :: nZonesInPlane           !! seems this is the proper way to do it
       integer (c_int) :: nextZ                  !! 
       integer (c_int) :: nextC
       real (c_double) :: STotal                 !! 
       real (c_double) :: tau
       real (c_double) :: Psi
       integer (c_int) :: Groups                    !!
       real (c_double) :: Volume
       real (c_double) :: Sigt
       integer (c_int) :: nCFacesArray
       integer (c_int) :: ndim
       integer (c_int) :: maxcf
       integer (c_int) :: ncorner
       real (c_double) :: A_fp
       real (c_double) :: omega
       integer (c_int) :: cFP
       real (c_double) :: Psi1
       integer (c_int) :: nbelem
       real (c_double) :: A_ez
       integer (c_int) :: cEZ
       integer (c_int) :: NumAngles
       real (c_double) :: quadwt
       real (c_double) :: Phi
       real (c_double) :: PsiB
       integer (c_int) :: maxCorner
       integer (c_int) :: mem0solve1
       integer (c_int) :: threadNum
       integer (c_int) :: numThreads
       integer (c_int) :: savePsi
       integer (c_int) :: numCycles
       integer (c_int) :: cycleOffset
       real (c_double) :: cyclePsi
       integer (c_int) :: cycleList
       integer (c_int) :: b0 
       integer (c_int) :: nBdyElem
       real (c_double) :: PsiBMref       
       integer (c_int) :: Mref
       integer (c_int) :: Geo_numCorner
       integer (c_int) :: Geo_cOffSet
     end subroutine gpu_sweepucbxyz
  end interface

  !  Arguments
  type(SetData),    intent(inout) :: Set
 
  integer,          intent(in)    :: setID 
  integer,          intent(in)    :: Groups 
  integer,          intent(in)    :: Angle 

  logical (kind=1), intent(in)    :: savePsi

  real(adqt),       intent(inout) :: Psi1(Groups,Set%nCorner+Set%nbelem)

  !  Local Variables

  type(AngleSet),   pointer       :: ASet
  type(GroupSet),   pointer       :: GSet

  integer    :: nHyperplanes
  integer    :: intSavePsi
  integer    :: nReflecting
  integer    :: nBdyElem
  integer    :: b0
  integer    :: Mref

  type(HypPlane), pointer :: HypPlanePtr
  integer :: firstAngle

  type(Boundary), pointer :: BdyT

  ASet         => getAngleSetFromSetID(Quad, setID)
  GSet         => getGroupSetFromSetID(Quad, setID)

  nHyperPlanes = getNumberOfHyperPlanes(ASet,Angle)

  !  Initialize boundary values in Psi1 and interior values on the cycle list
  ! moved to GPU
  ! call initFromCycleList(Set, angle, Psi1)

  ! moved to GPU
  !  do ib=1,Set%nbelem
  !     Psi1(:,Set%nCorner+ib) = Set% PsiB(:,ib,Angle)
  !  enddo

  if ( savePsi ) then
     intSavePsi = 1
  else
     intSavePsi = 0
  endif

  if ( Set% AngleOrder(1) .eq. Angle ) then
     firstAngle = 1
  else
     firstAngle = 0
  endif

  !! DIVERT TO GPU HERE
  HypPlanePtr => ASet%HypPlanePtr(Angle)

  ! Flag '0' indicates this will just copy the data to the GPU.
  ! point being, we need to copy the input data before the cpu overwrites it (if we are doing the cpu computation alongside) 

  nReflecting = getNumberOfReflecting(RadBoundary)

  ! RCC - Steve says can only hand 1 reflecting boundary, check snreflect.F90
  ! for actual implementation
  if ( nReflecting > 1 ) then
     print *, "ERROR - as yet we can only handle nReflecting == 1"
  end if
  
  if (nReflecting .eq. 1 ) then
     BdyT => getReflecting(RadBoundary, 1)
     nBdyElem = getNumberOfBdyElements(BdyT)
     b0 = getFirstBdyElement(BdyT)-1
     Mref = getReflectedAngle(ASet, 1, Angle)
     if ( Mref < 1 ) then
        nBdyElem = 0;
!     else 
!        print *, "snreflect (device) inputs nReflecgting, Minc, Mref, b0, nBdyElem = ", nReflecting, Angle, Mref, b0, nBdyElem
     end if
  else
     nBdyElem = 0
     b0 = 0
     Mref = -1
  end if

! TODO - Need to add ASet%cycleOffSet(angle) to the parameters, so we can update
! kernel to do this loop correctly:
! do m=1,ASet% numCycles(angle)
!  mCycle  = offSet + m
!  c0      = ASet% cycleList(1,mCycle)
!  nCorner = ASet% cycleList(2,mCycle)
!  do c=1,nCorner
!     Psi1(:,c0+c) = Set% cyclePsi(:,c,mCycle)
! enddo
!


! TODO - Failing to compile - these args changed from allocatables to pointers.

#if defined(TETON_ENABLE_CUDA)
!$omp critical  
! RCC - As long as kernels operate on independent data, omp critical section
!       is not necessary. However, it is kept to aid better ordering in the
!       NVIDIA profiler visualization.

! Do the iso c bindings need to change?
  CALL gpu_sweepucbxyz ( &
       Angle, &
       nHyperPlanes, &
       HypPlanePtr%zonesInPlane(1), &
       ASet%nextZ(1,Angle), &
       ASet%nextC(1,Angle), &
       GSet%STotal(1,1), &
       Size%tau,                      &
       Set%Psi(1,1,Angle), &
       Groups, &
       Geom%Volume(1),                                                                                               &
       GSet%Sigt(1,1), &
       Geom%nCFacesArray(1), &
       Size%ndim, &
       Size%maxcf, &
       Set%nCorner, &
       Geom%A_fp(1,1,1), &
       ASet%omega(1,Angle), &
       Geom%cFP(1,1), &
       Psi1(1,1),    &
       Set%nbelem, &
       Geom%A_ez(1,1,1), &
       Geom%cEZ(1,1), &
       Set%NumAngles, &
       ASet%weight(Angle), &
       Set%Phi(1,1), &
       Set%PsiB(1,1,Angle), &
       Size%maxCorner, &
       firstAngle, &
       setID, &
       omp_get_num_threads(), &
       intSavePsi, &
       ASet%numCycles(Angle), &
       ASet%cycleOffSet(angle), &
       Set%cyclePsi(1,1),      &
       ASet%cycleList(1), &
       b0, &
       nBdyElem, &
       Set%PsiB(1,1,Mref), &
       Mref, &
       Geom%numCorner(1), &
       Geom%cOffSet(1) &
       )

!$omp end critical

  call gpu_streamSynchronize(mod(setId,80))
  
! TODO - should we trust Set%cyclePsi?  this used to be hyperplane%cyclePsi.
  ! moved to GPU
  ! call updateCycleList(Set, angle, Psi1)
#endif

  return
end subroutine SweepUCBxyzToGPU
