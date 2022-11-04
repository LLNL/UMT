!***********************************************************************
!                        Version 1:  09/2017, PFN                      *
!                                                                      *
!   SetSweep      - This routine controls the computational set        *
!                   sweeps and communication when running sweeps       *
!                   on a CPU.                                          *
!                                                                      *
!***********************************************************************

   subroutine SetSweep(savePsi)


   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use SetData_mod
   use CommSet_mod
   use AngleSet_mod
   use mpi_param_mod
   use mpif90_mod
   use iter_control_list_mod
   use iter_control_mod

   implicit none

!  Arguments

   logical (kind=1), intent(in)  :: savePsi

!  Local

   type(SetData),     pointer    :: Set
   type(CommSet),     pointer    :: CSet
   type(AngleSet),    pointer    :: ASet
   type(IterControl), pointer    :: incidentFluxControl => NULL()

   integer                       :: setID
   integer                       :: cSetID
   integer                       :: Angle
   integer                       :: NumAnglesDyn
   integer                       :: ndim 
   integer                       :: Groups
   integer                       :: sendIndex
   integer                       :: fluxIter
   integer                       :: nCommSets
   integer                       :: nConv
   integer                       :: nNotConv
   integer                       :: maxIters

   logical (kind=1)              :: SnSweep 

   logical (kind=1), allocatable :: FluxConverged(:)

!  Constants

   incidentFluxControl => getIterationControl(IterControls,"incidentFlux")
   maxIters            =  getMaxNumberOfIterations(incidentFluxControl)
   nCommSets           =  getNumberOfCommSets(Quad)
   ndim                =  Size% ndim
   SnSweep             = .TRUE.

   allocate( FluxConverged(nCommSets) )

!  If this is the first flux iteration, initialize the communication
!  order and incident flux on shared boundaries

   do cSetID=1,nCommSets

     CSet => getCommSetData(Quad, cSetID)

     call restoreCommOrder(CSet)
     call setIncidentFlux(cSetID)
   enddo

!  Begin Flux Iteration

   FluxConverged(:) = .FALSE.
   fluxIter         =  0

   FluxIteration: do

     fluxIter = fluxIter + 1

!  Initialize

!$omp parallel do default(none) schedule(dynamic) &
!$omp& private(CSet,Set)  &
!$omp& shared(nCommSets,Quad,ndim)
     do cSetID=1,nCommSets
       CSet => getCommSetData(Quad, cSetID)

       do setID=CSet% set1,CSet% set2
         Set => getSetData(Quad, setID)

         if (ndim == 2) then
           Set% PsiM(:,:) = zero
         endif

         Set% Phi(:,:) = zero
       enddo
     enddo
!$omp end parallel do

!  Post receives for all data

     do cSetID=1,nCommSets
       call InitExchange(cSetID)
     enddo

!  Loop over angles, solving for each in turn:

!$omp parallel do default(none) schedule(dynamic) &
!$omp& private(CSet,ASet,NumAnglesDyn,Angle,Set,Groups) &
!$omp& shared(nCommSets,SnSweep,savePsi,ndim,Quad)
     CommLoop: do cSetID=1,nCommSets

       CSet         => getCommSetData(Quad, cSetID)
       ASet         => CSet% AngleSetPtr
       NumAnglesDyn =  CSet% NumAnglesDyn

       AngleLoop: do sendIndex=1,NumAnglesDyn

         Angle =  CSet% AngleOrder(sendIndex)

!        Send the boundary information needed by my neighbors
         call SendFlux(SnSweep, cSetID, sendIndex)

!        Test for completion of the sends needed by my neighbors
         call TestSend(cSetID, sendIndex)

!        Receive the boundary information needed to compute this angle
         call RecvFlux(SnSweep, cSetID, Angle)


!        Sweep the mesh, calculating PSI for each corner; the 
!        boundary flux array PSIB is also updated here. 
!        Mesh cycles are fixed automatically.

         SetLoop: do setID=CSet% set1,CSet% set2

           Set    => getSetData(Quad, setID)
           Groups =  Set% Groups

           AngleType: if ( .not. ASet% FinishingDirection(Angle) ) then

             call snreflect(SnSweep, setID, Angle)

             call initFromCycleList(setID, Angle, Groups, Set% Psi1)

             if (ndim == 3) then
               call SweepUCBxyz(Set, setID, Groups, Angle, savePsi)
             elseif (ndim == 2) then
!              if ( Size% usePWLD ) then
!                call SweepPWLDrz(Set, setID, Groups, Angle, savePsi)
!              else
               call SweepUCBrz(Set, setID, Groups, Angle, savePsi)
!              endif
             endif

             call updateCycleList(setID, Angle, Groups, Set% Psi1)

           endif AngleType

         enddo SetLoop

       enddo AngleLoop

     enddo CommLoop
!$omp end parallel do

!    Test convergence of incident fluxes

!$omp parallel do default(none) schedule(static) &
!$omp& shared(nCommSets,FluxConverged)
     do cSetID=1,nCommSets
       call setIncidentFlux(cSetID)
       call testFluxConv(cSetID, FluxConverged(cSetID))
     enddo
!$omp end parallel do

!    If this is the end of the radiation step and we are saving Psi do
!    not perform additional sweeps

     if (savePsi) then
       exit FluxIteration
     endif

     nConv = 0
     do cSetID=1,nCommSets
       if ( FluxConverged(cSetID) ) then
         nConv = nConv + 1
       endif
     enddo

     nNotConv = nCommSets - nConv

!    Make sure all processes are in sync
     call MPIAllReduce(nNotConv, "max", MY_COMM_GROUP)

     if ( nNotConv == 0 .or. fluxIter >= maxIters ) then
       exit FluxIteration
     else
       cycle FluxIteration
     endif

   enddo FluxIteration


   deallocate( FluxConverged )


   return
   end subroutine SetSweep 


