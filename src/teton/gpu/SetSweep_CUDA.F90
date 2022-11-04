!***********************************************************************
!                        Version 1:  01/2021, RCC                      *
!                                                                      *
!   SetSweep      - This routine controls the computational set        *
!                   sweeps and communication when running CUDA sweeps  *
!                   on a GPU.                                          *
!                                                                      *
!***********************************************************************

subroutine SetSweep_CUDA(savePsi)


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

  interface
     subroutine gpu_streamSynchronize ( streamId ) bind (c)
       use iso_c_binding
       integer (c_int) :: streamId
     end subroutine gpu_streamSynchronize
  end interface

  !  Arguments

  logical (kind=1), intent(in) :: savePsi


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
!$omp& shared(nCommSets,Quad,ndim) &
!$omp& private(CSet,Set)
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

     ! Post receives for all data

     do cSetID=1,nCommSets
       call InitExchange(cSetID)
     enddo

     !  Loop over angles, solving for each in turn:

     ! STRATEGY: Eventually, this loop should be on the GPU
     ! This is a loop over all angles in a set of angles with equivalent hyperplanes, etc.
     ! For _now_ just get the kernel working per single angle.

!$omp parallel do default(none) schedule(dynamic) &
!$omp& shared(nCommSets,SnSweep,savePsi,ndim,Quad) &
!$omp& private(CSet,ASet,NumAnglesDyn,Angle,Set,Groups)
     CommLoop: do cSetID=1,nCommSets

       CSet         => getCommSetData(Quad, cSetID)
       ASet         => CSet% AngleSetPtr
       NumAnglesDyn =  CSet% NumAnglesDyn

       AngleLoop: do sendIndex=1,NumAnglesDyn

          Angle = CSet% AngleOrder(sendIndex)

          ! STRATEGY: MPI communication should all be pipeines w/ packing/memcopies/GPUsweep.
          ! again, that can be implemented later - after the kernel is verified to work.

          !      Send the boundary information needed by my neighbors
          call SendFlux(SnSweep, cSetID, sendIndex)

          !      Test for completion of the sends needed by my neighbors
          call TestSend(cSetID, sendIndex)

          !      Receive the boundary information needed to compute this angle
          call RecvFlux(SnSweep, cSetID, Angle)


          !      Sweep the mesh, calculating PSI for each corner; the 
          !      boundary flux array PSIB is also updated here. 
          !      Mesh cycles are fixed automatically.

          SetLoop: do setID=CSet% set1,CSet% set2

            Set    => getSetData(Quad, setID)
            Groups =  Set% Groups

            AngleType: if ( .not. ASet% FinishingDirection(Angle) ) then

               ! RCC done within GPU_sweepucbxyz.cu call snreflect(SnSweep, setID, Angle)

               ! RCC done within GPU_sweepucbxyz.cu call initFromCycleList(setID, Angle, Groups, Set% Psi1)

               if (ndim == 3) then 

  ! snreflect and initFromCycleList are handled within SweepUCBxyzToGPU. (check
  ! this)

  ! Steve is not using Set Psi1 anymore.  He creates Psi1 solely on the GPU.
  ! Can probably remove this from args ( test later ).

                  call SweepUCBxyzToGPU (Set, setID, Groups, Angle, savePsi, Set%Psi1(:,:))

               elseif (ndim == 2) then
                  call f90fatal("CUDA Sweep is not available for 2d.")
               endif

             ! RCC done within GPU_sweepucbxyz.cu call updateCycleList(setID, Angle, Groups, Set% Psi1)

             endif AngleType

         enddo SetLoop

       enddo AngleLoop

     enddo CommLoop
!$omp end parallel do

#if defined(TETON_ENABLE_CUDA)
     CALL gpu_streamSynchronize ( mod(cSetID,20) ) 
#endif

!$omp parallel do default(none) schedule(static) &
!$omp& shared(nCommSets,FluxConverged)
    do cSetID=1,nCommSets
      call setIncidentFlux(cSetID)
      call testFluxConv(cSetID, FluxConverged(cSetID))
    enddo
!$omp end parallel do

!   If this is the end of the radiation step and we are saving Psi do
!   not perform additional sweeps

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

!   Make sure all processes are in sync
    call MPIAllReduce(nNotConv, "max", MY_COMM_GROUP)

    if ( nNotConv == 0 .or. fluxIter >= maxIters ) then
      exit FluxIteration
    else
      cycle FluxIteration
    endif

  enddo FluxIteration


  deallocate( FluxConverged )


  return
end subroutine SetSweep_CUDA


