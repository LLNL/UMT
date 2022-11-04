!***********************************************************************
!                        Version 1:  09/96, PFN                        *
!                                                                      *
!   SWEEPSCHEDULER - This routine determines the order in which the    *
!                    discrete ordinates will be solved.                *
!                                                                      *
!                                                                      *
! The sweep scheduler looks at three criteria for ordering the angles: *
! 1. Reflecting Boundaries                                             *
! 2. Angular Derivatives (RZ only)                                     *
! 3. The net flux (incident - outgoing) on the shared boundaries for   *
!    each angle                                                        *
!                                                                      *
! Scheduling the order of angles to solve first takes into #1 and #2   *
! to account for dependencies.                                         *
!                                                                      *
! Then, a weight for each angle is calculated based on the net flux    *
! (incident - outgoing ) and angle dependencies                        *
!                                                                      *
! The scheduler will order the angles from smallest to largest weight  *
!                                                                      *
! This prioritizes solving angles with the largest outgoing flux       * 
! first, and defers solving angles with a large incident flux from     *
! neighbors.  To intent is to allow neighboring domains to solve those *
! angles first and minimize the amount of lagged angle data.           *
!                                                                      *
! Important note:                                                      *
! The other important mechanism is the ability to dynamically iterate  *
! and repeat sweeps on angles with lagged data from neighbors to       *
! improve convergence.                                                 *
!***********************************************************************
   subroutine SweepScheduler(cSetID)

   use kind_mod
   use constant_mod
   use mpif90_mod
   use mpi
   use Size_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use Boundary_mod
   use SetData_mod
   use CommSet_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,          intent(in) :: cSetID

!  Local Variables

   type(SetData),      pointer  :: Set
   type(CommSet),      pointer  :: CSet
   type(AngleSet),     pointer  :: ASet
   type(Boundary),     pointer  :: BdyT

   integer :: setID
   integer :: i 
   integer :: ndone 
   integer :: neighbor
   integer :: ierr
   integer :: ishared
   integer :: nsend
   integer :: myRankInGroup 
   integer :: NumBin0
   integer :: NumAngles
   integer :: NangBin
   integer :: nShared
   integer :: nReflecting
   integer :: angle
   integer :: angleRef
   integer :: bin
   integer :: newbin
   integer :: binRef
   integer :: binRecv 
   integer :: n
   integer :: recvID
   integer :: angleID
   integer :: COMM_GROUP
   integer :: tag
   integer :: tagBase

   integer, dimension (1) :: imin
   integer, dimension (1) :: imax

!  Dynamic

   real(adqt),       allocatable :: depend(:)
   real(adqt),       allocatable :: weightComm(:,:)
   integer,          allocatable :: nRefl(:)
   integer,          allocatable :: BinOffSet(:)
   integer,          allocatable :: depAngle(:,:)
   integer,          allocatable :: binDone(:,:)
   integer,          allocatable :: nRecv(:)
   logical (kind=1), allocatable :: notDone(:)
   logical (kind=1), allocatable :: Ready(:)

!  Mesh Constants

   CSet          => getCommSetData(Quad, cSetID)
   ASet          => CSet% AngleSetPtr
   myRankInGroup =  Size% myRankInGroup 
   nShared       =  getNumberOfShared(RadBoundary)
   nReflecting   =  getNumberOfReflecting(RadBoundary)
   COMM_GROUP    =  getCommunicationGroup(CSet)
   tagBase       =  8000 + 100*cSetID

   NumBin0       =  CSet% NumBin0
   NumAngles     =  CSet% NumAngles 

   call MPIBarrier(COMM_GROUP)

!  Allocate temporaries

   allocate( depend(NumBin0) )
   allocate( weightComm(nShared,NumBin0) )
   allocate( nRefl(NumBin0) )
   allocate( BinOffSet(NumBin0) )
   allocate( depAngle(nReflecting,NumAngles) )
   allocate( binDone(nShared,NumBin0) )
   allocate( nRecv(nShared) )
   allocate( notDone(NumBin0) )
   allocate( Ready(NumBin0) )

!  Tally message size by angle bin; "angle-bin" is a xi-level in 2d
!  and an angle in 3D

   DecompTest: if (Size% ncomm > 0) then
     weightComm(:,:) = CSet% NetFlux(:,:)
   else
     weightComm(:,:) = one 
   endif DecompTest

   BinOffSet(1) = 0
   do bin=2,NumBin0
     BinOffSet(bin) = BinOffSet(bin-1) + CSet% NangBinList(bin-1)
   enddo

!  Initialize the dependency-weight to the message weight.
!  In addition, add the message weight from angles that
!  reflect into each angle. 

   depend(:) = zero  

   do bin=1,NumBin0
     do ishared=1,nShared
       depend(bin) = depend(bin) + weightComm(ishared,bin)
     enddo
   enddo

!  Set the reflecting boundary dependencies 

   nRefl(:)      = 0 
   depAngle(:,:) = 0

   do n=1,nReflecting
     do angle=1,NumAngles
       angleRef = getReflectedAngle(ASet, n, angle)
       bin      = CSet% AngleToBin(angle)
       if (angleRef > 0) then
         depAngle(n,angleRef) = angle
         nRefl(bin)           = nRefl(bin) + 1 
       endif
     enddo
   enddo

!  Order the angles starting with those that have the
!  smallest (or no) dependencies 

   ndone        = 0
   nsend        = 1
   binDone(:,:) = -1
   nRecv(:)     = 0
   notDone(:)   = .TRUE.

   OuterIteration: do i=1,NumBin0

!    Post Receives

     do ishared=1,nShared
       BdyT     => getShared(RadBoundary, ishared)
       neighbor =  getNeighborID(BdyT)
       tag      =  tagBase + myRankInGroup + neighbor

       call MPI_Irecv(binDone(ishared,i), nsend, MPI_INTEGER, neighbor,  &
                      tag, COMM_GROUP, CSet% request(2*ishared), ierr)
     enddo

!  Select a new angle bin

     imin = minloc( nRefl, notDone )
     if (nRefl(imin(1)) /= 0) then
       nRefl(imin(1)) = 0
     endif

     do bin=1,NumBin0
       if ( notDone(bin) ) then
         if (nRefl(bin) == 0) then
           Ready(bin) = .TRUE.
         else
           Ready(bin) = .FALSE.
         endif
       else
         Ready(bin) = .FALSE.
       endif
     enddo

     imax    = maxloc( depend, Ready )
     newbin  = imax(1)

     NangBin = CSet% NangBinList(newbin)

!    Send my bin to all neighbors

     do ishared=1,nShared
       BdyT     => getShared(RadBoundary, ishared)
       neighbor =  getNeighborID(BdyT)
       tag      =  tagBase + myRankInGroup + neighbor

       call MPI_Isend(newbin, nsend, MPI_INTEGER, neighbor, &
                      tag, COMM_GROUP, CSet% request(2*ishared-1), ierr)

       call MPIWait(CSet% request(2*ishared-1))
     enddo

!    Decrement for reflecting boundaries
 
     do angleID=1,NangBin
       ndone = ndone + 1
       angle = BinOffSet(newbin) + angleID 

       CSet% AngleOrder(ndone)  = angle 
       CSet% AngleOrder0(ndone) = angle

       do n=1,nReflecting
         angleRef = depAngle(n,angle)
         if (angleRef > 0) then
           binRef = CSet% AngleToBin(angleRef)
           if ( notDone(binRef) ) then
             nRefl(binRef) = nRefl(binRef) - 1 
           endif
         endif
       enddo
     enddo

     if ( notDone(newbin) ) then
       notDone(newbin)   = .FALSE.
     else
       write(6,888)
 888   format("ATTEMPTING to REPEAT an ANGLE!!!!!!!!!!!!!")
     endif

!    Decrement depend for shared boundaries 

     do ishared=1,nShared
       call MPIWait(CSet% request(2*ishared))
       binRecv = binDone(ishared,i)
       NangBin = CSet% NangBinList(binRecv)
       recvID  = nRecv(ishared)

       do angleID=1,NangBin
         CSet% RecvOrder0(recvID+angleID,ishared) = BinOffSet(binRecv) + angleID 
         CSet% RecvOrder(recvID+angleID,ishared)  = BinOffSet(binRecv) + angleID 
       enddo

       nRecv(ishared) = nRecv(ishared) + NangBin

       if ( notDone(binRecv) ) then
         depend(binRecv) = depend(binRecv) - weightComm(ishared,binRecv)
       endif
     enddo

     call MPIBarrier(COMM_GROUP)

   enddo OuterIteration

   do setID=CSet% set1,CSet% set2
     Set => getSetData(Quad, setID)

     do angle=1,NumAngles
       Set% AngleOrder(angle) = CSet% AngleOrder(angle)
     enddo
   enddo

   do bin=1,NumBin0
     if ( notDone(bin) ) then
       write(6,250) myRankInGroup, bin
     endif
   enddo 

 250 format("meshDomain = ",i2,2x,"angle bin = ",i2," was missed")


!  Release Temporaries

   deallocate( depend )
   deallocate( weightComm )
   deallocate( nRefl )
   deallocate( BinOffSet )
   deallocate( depAngle )
   deallocate( binDone )
   deallocate( nRecv )
   deallocate( notDone )
   deallocate( Ready )

   call MPIBarrier(COMM_GROUP)



   return
   end subroutine SweepScheduler 

