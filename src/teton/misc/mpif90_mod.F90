module mpif90_mod

#include "macros.h"
use kind_mod
use mpi_param_mod

!=======================================================================
!                       Version 1.1: 02/99, MRZ
!                       Version 1.0: 05/92, PFN
!-----------------------------------------------------------------------
! MPI
!   This class wraps MPI functions.  This facilitates turning off MPI.
!
!-----------------------------------------------------------------------
! v1.0: Original implementation
! v1.1: MPI functions wrapped in a Fortran90 class
!=======================================================================

private

#include "mpif90.if"

!   Public Interfaces

    integer, public :: MY_COMM_GROUP

contains

!=======================================================================
! MPIAllReduce interface
!=======================================================================

  subroutine mpi_MPIAllReduce_r(recvBuf,mpiOp,comm)

!    This subroutine performs an MPI All Reduce operation with a
!    barrier.  The data to be broadcast, a floating point scalar, is
!    passed in as recvBuf; the reduced data is returned in recvBuf.
!
!      recvBuf  send and receive buffer
!      mpiOp    MPI operation
!      comm     MPI communicator

!    variable declarations
     implicit none

!    passed variables
     real(double),   intent(inout) :: recvBuf
     character(*), intent(in)    :: mpiOp
     integer,      intent(in)    :: comm

!    local variables
     integer    :: length, ierror
     real(double) :: sendBuf

!     character(4), dimension(4) :: mpiOps = &
!                                   (/"min ","max ","prod","sum "/)

!    assertions
!     TETON_ASSERT(any(mpiOp==mpiOps(:)),"Invalid MPI Reduction Operation")

!      copy the send buffer into temporary storage
       sendBuf = recvBuf

!      MPI Barrier is implicit for MPI_Allreduce

!      MPI Reduction
       length = 1
       select case (mpiOp)
       case ("min")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_REAL8, &
            MPI_MIN, comm, ierror)
       case ("max")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_REAL8, &
            MPI_MAX, comm, ierror)
       case ("prod")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_REAL8, &
            MPI_PROD, comm, ierror)
       case ("sum")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_REAL8, &
            MPI_SUM, comm, ierror)
       end select

       if (ierror /= MPI_SUCCESS) then
          call f90fatal("MPI Reduction Failed")
       endif

     return
  end subroutine mpi_MPIAllReduce_r

!-----------------------------------------------------------------------
  subroutine mpi_MPIAllReduce_r_(recvBuf,mpiOp,comm)

!    This subroutine performs an MPI All Reduce operation with a
!    barrier.  The data to be broadcast, a floating point array, is
!    passed in as recvBuf; the reduced data is returned in recvBuf.
!
!      recvBuf  send and receive buffer
!      mpiOp    MPI operation
!      comm     MPI communicator

!    variable declarations
     implicit none

!    passed variables
     real(double),   contiguous, intent(inout) :: recvBuf(:)
     character(*),             intent(in)    :: mpiOp
     integer,                  intent(in)    :: comm

!    local variables
     integer                 :: length, ierror, alloc_stat
     real(double), allocatable :: sendBuf(:)

!     character(4), dimension(4) :: mpiOps = &
!                                   (/"min ","max ","prod","sum "/)

!    assertions
!     TETON_ASSERT(any(mpiOp==mpiOps(:)),"Invalid MPI Reduction Operation")

!      allocate memory for the send buffer
       allocate(sendBuf(size(recvBuf)))

!      copy the send buffer into temporary storage
       sendBuf(:) = recvBuf(:)

!      MPI Barrier is implicit for MPI_Allreduce

!      MPI Reduction
       length = size(recvBuf(:))
       select case (mpiOp)
       case ("min")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_REAL8, &
            MPI_MIN, comm, ierror)
       case ("max")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_REAL8, &
            MPI_MAX, comm, ierror)
       case ("prod")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_REAL8, &
            MPI_PROD, comm, ierror)
       case ("sum")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_REAL8, &
            MPI_SUM, comm, ierror)
       end select

       if (ierror /= MPI_SUCCESS) then
          call f90fatal("MPI Reduction Failed")
       endif

!      free memory
       deallocate(sendBuf, stat=alloc_stat)

     return
  end subroutine mpi_MPIAllReduce_r_

!-----------------------------------------------------------------------
  subroutine mpi_MPIAllReduce_i(recvBuf,mpiOp,comm)

!    This subroutine performs an MPI All Reduce operation with a
!    barrier.  The data to be broadcast, an integer scalar, is
!    passed in as recvBuf; the reduced data is returned in recvBuf.
!
!      recvBuf  send and receive buffer
!      mpiOp    MPI operation
!      comm     MPI communicator

!    variable declarations
     implicit none

!    passed variables
     integer,      intent(inout) :: recvBuf
     character(*), intent(in)    :: mpiOp
     integer,      intent(in)    :: comm

!    local variables
     integer    :: length, ierror
     integer    :: sendBuf

!     character(4), dimension(4) :: mpiOps = &
!                                   (/"min ","max ","prod","sum "/)

!    assertions
!     TETON_ASSERT(any(mpiOp==mpiOps(:)),"Invalid MPI Reduction Operation")

!      copy the send buffer into temporary storage
       sendBuf = recvBuf

!      MPI Barrier is implicit for MPI_Allreduce

!      MPI Reduction
       length = 1
       select case (mpiOp)
       case ("min")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_INTEGER, &
            MPI_MIN, comm, ierror)
       case ("max")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_INTEGER, &
            MPI_MAX, comm, ierror)
       case ("prod")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_INTEGER, &
            MPI_PROD, comm, ierror)
       case ("sum")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_INTEGER, &
            MPI_SUM, comm, ierror)
       end select

       if (ierror /= MPI_SUCCESS) then
          call f90fatal("MPI Reduction Failed")
       endif

     return
  end subroutine mpi_MPIAllReduce_i

!-----------------------------------------------------------------------
  subroutine mpi_MPIAllReduce_i_(recvBuf,mpiOp,comm)

!    This subroutine performs an MPI All Reduce operation with a
!    barrier.  The data to be broadcast, an integer array, is
!    passed in as recvBuf; the reduced data is returned in recvBuf.
!
!      recvBuf  send and receive buffer
!      mpiOp    MPI operation
!      comm     MPI communicator

!    variable declarations
     implicit none

!    passed variables
     integer,      contiguous, intent(inout) :: recvBuf(:)
     character(*),             intent(in)    :: mpiOp
     integer,                  intent(in)    :: comm

!    local variables
     integer              :: length, ierror, alloc_stat
     integer, allocatable :: sendBuf(:)

!     character(4), dimension(4) :: mpiOps = &
!                                   (/"min ","max ","prod","sum "/)

!    assertions
!     TETON_ASSERT(any(mpiOp==mpiOps(:)),"Invalid MPI Reduction Operation")

!      allocate memory for the send buffer
       allocate(sendBuf(size(recvBuf)))

!      copy the send buffer into temporary storage
       sendBuf(:) = recvBuf(:)

!      MPI Barrier is implicit for MPI_Allreduce

!      MPI Reduction
       length = size(recvBuf(:))
       select case (mpiOp)
       case ("min")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_INTEGER, &
            MPI_MIN, comm, ierror)
       case ("max")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_INTEGER, &
            MPI_MAX, comm, ierror)
       case ("prod")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_INTEGER, &
            MPI_PROD, comm, ierror)
       case ("sum")
          call MPI_Allreduce(sendBuf, recvBuf, length, MPI_INTEGER, &
            MPI_SUM, comm, ierror)
       end select

       if (ierror /= MPI_SUCCESS) then
          call f90fatal("MPI Reduction Failed")
       endif

!      free memory
       deallocate(sendBuf, stat=alloc_stat)

     return
  end subroutine mpi_MPIAllReduce_i_

!=======================================================================
! MPIBarrier interface
!=======================================================================

  subroutine mpi_MPIBarrier(comm)

!    This subroutine performs an MPI Barrier operation.
!      comm   MPI communicator

!    variable declarations
     implicit none

!    passed variables
     integer, intent(in) :: comm

!    local variables
     integer :: ierror

!    MPI Barrier
     call MPI_Barrier(comm, ierror)
     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Barrier Failed")
     endif

     return
  end subroutine mpi_MPIBarrier

!=======================================================================
! MPIAbort interface
!=======================================================================

  subroutine mpi_MPIAbort(comm)

!    This subroutine performs an MPI Abort operation.
!      comm   MPI communicator

!    variable declarations
     implicit none

!    passed variables
     integer, intent(in) :: comm

!    local variables
     integer :: errorcode = 1
     integer :: ierror

!    MPI Abort 
     call MPI_Abort(comm, errorcode, ierror)
     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Abort Failed")
     endif

     return
  end subroutine mpi_MPIAbort

!=======================================================================
! MPICommRank interface
!=======================================================================

  subroutine mpi_MPICommRank(comm,rank)

!    This subroutine determines the rank of the calling process in the
!    communicator.
!      rank   processor rank
!      comm   MPI communicator

!    variable declarations
     implicit none

!    passed variables
     integer, intent(in)  :: comm
     integer, intent(out) :: rank

!    local variables
     integer :: ierror

!    MPI Communicator Rank
     call MPI_Comm_rank(comm, rank, ierror)
     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Barrier Failed")
     endif

     return
  end subroutine mpi_MPICommRank

!=======================================================================
! MPICommSize interface
!=======================================================================

  subroutine mpi_MPICommSize(comm,commSize)

!    This subroutine determines the size of the group associated with
!    the communicator.
!      comm       MPI communicator
!      commSize   communicator size

!    variable declarations
     implicit none

!    passed variables
     integer, intent(in)  :: comm
     integer, intent(out) :: commSize

!    local variables
     integer :: ierror

!    MPI Communicator Size
     call MPI_Comm_size(comm, commSize, ierror)
     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Barrier Failed")
     endif

     return
  end subroutine mpi_MPICommSize

!=======================================================================
! MPIGather interface
!=======================================================================

  subroutine mpi_MPIGather_r_(sendBuf,recvBuf,gatherNode,comm)

!    This subroutine performs an MPI Gather operation.  The data to
!    be broadcast, a floating point array, is passed in as send Buf;
!    the gathered data is returned in recvBuf.
!
!      sendBuf      send buffer
!      recvBuf      receive buffer
!      gatherNode   node to which all data is gathered
!      comm         MPI communicator

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(in)    :: sendBuf(:)
     real(double), contiguous, intent(inout) :: recvBuf(:,:)
     integer,    intent(in)    :: gatherNode
     integer,    intent(in)    :: comm

!    local variables
     logical :: tf

     integer    :: commSize, myNode, sendCount, recvCount, ierror
     real(double) :: recvBufDum(1,1)

!    determine size and rank
     commSize = getMPISize(comm)
     myNode   = getMPIRank(comm)

!    assertions
     if (myNode == gatherNode) then
        tf = size(recvBuf,1)==size(sendBuf,1)
        TETON_ASSERT(tf,"Invalid MPI Gather")
        tf = (size(recvBuf,2)==commSize)
        TETON_ASSERT(tf,"Invalid MPI Gather")
     endif

!    MPI Barrier before performing the gather
     call MPIBarrier(comm)

!    MPI Gather
     sendCount = size(sendBuf,1)
     recvCount = sendCount

     if (myNode == gatherNode) then
!       on the gather node, perform the gather operation into the
!       allocated receive buffer

        call MPI_Gather(sendBuf, sendCount, MPI_REAL8, &
                        recvBuf, recvCount, MPI_REAL8, &
                        gatherNode, comm, ierror)
     else
!       on non-gather nodes, the receiver buffer is dereferenced due
!       to a Fortran90 copy-in/copy-out operation.  To avoid
!       dereferencing a null pointer, pass a dummy (allocated) receive
!       buffer, which MPI ignores.

        call MPI_Gather(sendBuf, sendCount, MPI_REAL8, &
                        recvBufDum, recvCount, MPI_REAL8, &
                        gatherNode, comm, ierror)
     endif

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Reduction Failed")
     endif

     return
  end subroutine mpi_MPIGather_r_

!-----------------------------------------------------------------------
  subroutine mpi_MPIGather_r__(sendBuf,recvBuf,gatherNode,comm)

!    This subroutine performs an MPI Gather operation.  The data to
!    be broadcast, a floating point array, is passed in as send Buf;
!    the gathered data is returned in recvBuf.
!
!      sendBuf      send buffer
!      recvBuf      receive buffer
!      gatherNode   node to which all data is gathered
!      comm         MPI communicator

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(in)    :: sendBuf(:,:)
     real(double), contiguous, intent(inout) :: recvBuf(:,:,:)
     integer,    intent(in)    :: gatherNode
     integer,    intent(in)    :: comm

!    local variables
     logical(kind=1) :: tf

     integer    :: commSize, myNode, sendCount, recvCount, ierror
     real(double) :: recvBufDum(1,1,1)

!    determine size and rank
     commSize = getMPISize(comm)
     myNode   = getMPIRank(comm)

!    assertions
     if (myNode == gatherNode) then
        tf = size(recvBuf,1)==size(sendBuf,1)
        TETON_ASSERT(tf,"Invalid MPI Gather")
        tf = size(recvBuf,2)==size(sendBuf,2)
        TETON_ASSERT(tf,"Invalid MPI Gather")
        tf = size(recvBuf,3)==commSize
        TETON_ASSERT(tf,"Invalid MPI Gather")
     endif

!    MPI Barrier before performing the gather
     call MPIBarrier(comm)

!    MPI Gather
     sendCount = size(sendBuf,1)*size(sendBuf,2)
     recvCount = sendCount

     if (myNode == gatherNode) then
!       on the gather node, perform the gather operation into the
!       allocated receive buffer

        call MPI_Gather(sendBuf, sendCount, MPI_REAL8, &
                        recvBuf, recvCount, MPI_REAL8, &
                        gatherNode, comm, ierror)
     else
!       on non-gather nodes, the receiver buffer is dereferenced due
!       to a Fortran90 copy-in/copy-out operation.  To avoid
!       dereferencing a null pointer, pass a dummy (allocated) receive
!       buffer, which MPI ignores.

        call MPI_Gather(sendBuf, sendCount, MPI_REAL8, &
                        recvBufDum, recvCount, MPI_REAL8, &
                        gatherNode, comm, ierror)
     endif

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Reduction Failed")
     endif

     return
  end subroutine mpi_MPIGather_r__

!-----------------------------------------------------------------------
  subroutine mpi_MPIBcast_r(sendBuf,nsend,root,comm)

!    This subroutine performs an MPI broadcast operation.  The
!    data to be sent, a one-dimensional floating point array,
!    is passed in as sendBuf.
!
!      sendBuf      send buffer
!      nsend        length of send buffer
!      root         node that broadcasts sendBuf
!      comm         MPI communicator

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(inout) :: sendBuf(:)
     integer,                intent(in)    :: nsend
     integer,                intent(in)    :: root
     integer,                intent(in)    :: comm

!    local variables
     integer :: ierror

     call MPI_Bcast(sendBuf, nsend, MPI_REAL8, root, &
                    comm, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Bcast Failed")
     endif

     return
  end subroutine mpi_MPIBcast_r

!-----------------------------------------------------------------------
  subroutine mpi_MPIBcast_i(sendBuf,nsend,root,comm)

!    This subroutine performs an MPI broadcast operation.  The
!    data to be sent, a one-dimensional integer array,
!    is passed in as sendBuf.
!
!      sendBuf      send buffer
!      nsend        length of send buffer
!      root         node that broadcasts sendBuf
!      comm         MPI communicator

!    variable declarations
     implicit none

!    passed variables
     integer, contiguous, intent(inout) :: sendBuf(:)
     integer,             intent(in)    :: nsend
     integer,             intent(in)    :: root
     integer,             intent(in)    :: comm

!    local variables
     integer :: ierror

     call MPI_Bcast(sendBuf, nsend, MPI_INTEGER, root, &
                    comm, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Bcast Failed")
     endif

     return
  end subroutine mpi_MPIBcast_i

!-----------------------------------------------------------------------
  subroutine mpi_MPISendInit(sendBuf,nsend,isend,tag,comm,request)

!    This subroutine initalizes an MPI nonblocking send operation.  The
!    data to be sent is passed in as sendBuf.
!
!      sendBuf      send buffer
!      nsend        length of send buffer
!      isend        node to send to
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(in) :: sendBuf(:,:)
     integer,                intent(in) :: nsend
     integer,                intent(in) :: isend
     integer,                intent(in) :: tag
     integer,                intent(in) :: comm
     integer,                intent(out) :: request

!    local variables
     integer :: ierror

     call MPI_Send_Init(sendBuf, nsend, MPI_REAL8, isend, &
                        tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Send_init Failed")
     endif

     return
  end subroutine mpi_MPISendInit

!-----------------------------------------------------------------------
  subroutine mpi_MPISendInit1(sendBuf,nsend,isend,tag,comm,request)

!    This subroutine initalizes an MPI nonblocking send operation.  The
!    data to be sent is passed in as sendBuf.
!
!      sendBuf      send buffer
!      nsend        length of send buffer
!      isend        node to send to
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(in) :: sendBuf(:)
     integer,                intent(in) :: nsend
     integer,                intent(in) :: isend
     integer,                intent(in) :: tag
     integer,                intent(in) :: comm
     integer,                intent(out) :: request

!    local variables
     integer :: ierror

     call MPI_Send_Init(sendBuf, nsend, MPI_REAL8, isend, &
                        tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Send_init Failed")
     endif

     return
  end subroutine mpi_MPISendInit1

!-----------------------------------------------------------------------
  subroutine mpi_MPISendInit1_i(sendBuf,nsend,isend,tag,comm,request)

!    This subroutine initalizes an MPI nonblocking send operation.  The
!    data to be sent is passed in as sendBuf.
!
!      sendBuf      send buffer
!      nsend        length of send buffer
!      isend        node to send to
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     integer, contiguous, intent(in) :: sendBuf(:)
     integer,             intent(in) :: nsend
     integer,             intent(in) :: isend
     integer,             intent(in) :: tag
     integer,             intent(in) :: comm
     integer,             intent(out) :: request

!    local variables
     integer :: ierror

     call MPI_Send_Init(sendBuf, nsend, MPI_INTEGER, isend, &
                        tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Send_init Failed")
     endif

     return
  end subroutine mpi_MPISendInit1_i

!-----------------------------------------------------------------------
  subroutine mpi_MPIRecvInit(recvBuf,nrecv,irecv,tag,comm,request)

!    This subroutine initalizes an MPI nonblocking receive operation.
!    The received data is returned in as recvBuf.
!
!      recvBuf      recv buffer
!      nrecv        length of receive buffer
!      irecv        node to receive from
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(inout) :: recvBuf(:,:)
     integer,                intent(inout) :: nrecv
     integer,                intent(in)    :: irecv
     integer,                intent(in)    :: tag
     integer,                intent(in)    :: comm
     integer,                intent(out)   :: request

!    local variables
     integer :: ierror

     call MPI_Recv_Init(recvBuf, nrecv, MPI_REAL8, irecv, &
                        tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Recv_init Failed")
     endif

     return
  end subroutine mpi_MPIRecvInit

!-----------------------------------------------------------------------
  subroutine mpi_MPIRecvInit1(recvBuf,nrecv,irecv,tag,comm,request)

!    This subroutine initalizes an MPI nonblocking receive operation.
!    The received data is returned in as recvBuf.
!
!      recvBuf      recv buffer
!      nrecv        length of receive buffer
!      irecv        node to receive from
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(inout) :: recvBuf(:)
     integer,                intent(inout) :: nrecv
     integer,                intent(in)    :: irecv
     integer,                intent(in)    :: tag
     integer,                intent(in)    :: comm
     integer,                intent(out)   :: request

!    local variables
     integer :: ierror

     call MPI_Recv_Init(recvBuf, nrecv, MPI_REAL8, irecv, &
                        tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Recv_init Failed")
     endif

     return
  end subroutine mpi_MPIRecvInit1

!-----------------------------------------------------------------------
  subroutine mpi_MPIRecvInit1_i(recvBuf,nrecv,irecv,tag,comm,request)

!    This subroutine initalizes an MPI nonblocking receive operation.
!    The received data is returned in as recvBuf.
!
!      recvBuf      recv buffer
!      nrecv        length of receive buffer
!      irecv        node to receive from
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     integer, contiguous, intent(inout) :: recvBuf(:)
     integer,             intent(inout) :: nrecv
     integer,             intent(in)    :: irecv
     integer,             intent(in)    :: tag
     integer,             intent(in)    :: comm
     integer,             intent(out)   :: request

!    local variables
     integer :: ierror

     call MPI_Recv_Init(recvBuf, nrecv, MPI_INTEGER, irecv, &
                        tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Recv_init Failed")
     endif

     return
  end subroutine mpi_MPIRecvInit1_i

!-----------------------------------------------------------------------
  subroutine mpi_MPIIsend_r(sendBuf,nsend,isend,tag,comm,request)

!    This subroutine performs an MPI nonblocking send operation.  The
!    data to be sent, a one-dimensional floating point array, 
!    is passed in as sendBuf.
!
!      sendBuf      send buffer
!      nsend        length of send buffer
!      isend        node to send to
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(in) :: sendBuf(:)
     integer,                intent(in) :: nsend
     integer,                intent(in) :: isend
     integer,                intent(in) :: tag
     integer,                intent(in) :: comm
     integer,                intent(out) :: request

!    local variables
     integer :: ierror

     call MPI_Isend(sendBuf, nsend, MPI_REAL8, isend, &
                    tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Isend Failed")
     endif

     return
  end subroutine mpi_MPIIsend_r

!-----------------------------------------------------------------------
  subroutine mpi_MPIIsend_r2(sendBuf,nsend,isend,tag,comm,request)

!    This subroutine performs an MPI nonblocking send operation.  The
!    data to be sent, a two-dimensional floating point array, 
!    is passed in as sendBuf.
!
!      sendBuf      send buffer
!      nsend        length of send buffer
!      isend        node to send to
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(in) :: sendBuf(:,:)
     integer,                intent(in) :: nsend
     integer,                intent(in) :: isend
     integer,                intent(in) :: tag
     integer,                intent(in) :: comm
     integer,                intent(out) :: request

!    local variables
     integer :: ierror

     call MPI_Isend(sendBuf, nsend, MPI_REAL8, isend, &
                    tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Isend Failed")
     endif

     return
  end subroutine mpi_MPIIsend_r2

!-----------------------------------------------------------------------
  subroutine mpi_MPIIsend_i(sendBuf,nsend,isend,tag,comm,request)


!    This subroutine performs an MPI nonblocking send operation.  The
!    data to be sent, a one-dimensional integer array, 
!    is passed in as sendBuf.
!
!      sendBuf      send buffer
!      nsend        length of send buffer
!      isend        node to send to
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     integer, contiguous, intent(in) :: sendBuf(:)
     integer,             intent(in) :: nsend
     integer,             intent(in) :: isend
     integer,             intent(in) :: tag
     integer,             intent(in) :: comm
     integer,             intent(out) :: request

!    local variables
     integer :: ierror

     call MPI_Isend(sendBuf, nsend, MPI_INTEGER, isend, &
                    tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Isend Failed")
     endif

     return
  end subroutine mpi_MPIIsend_i

!-----------------------------------------------------------------------
  subroutine mpi_MPIIsend_i2(sendBuf,nsend,isend,tag,comm,request)


!    This subroutine performs an MPI nonblocking send operation.  The
!    data to be sent, a two-dimensional integer array, 
!    is passed in as sendBuf.
!
!      sendBuf      send buffer
!      nsend        length of send buffer
!      isend        node to send to
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     integer, contiguous, intent(in) :: sendBuf(:,:)
     integer,             intent(in) :: nsend
     integer,             intent(in) :: isend
     integer,             intent(in) :: tag
     integer,             intent(in) :: comm
     integer,             intent(out) :: request

!    local variables
     integer :: ierror

     call MPI_Isend(sendBuf, nsend, MPI_INTEGER, isend, &
                    tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Isend Failed")
     endif

     return
  end subroutine mpi_MPIIsend_i2

!-----------------------------------------------------------------------
  subroutine mpi_MPIIrecv_r(recvBuf,nrecv,irecv,tag,comm,request)

!    This subroutine performs an MPI nonblocking receive operation.
!    The received data, a one-dimensional floating point array, 
!    is returned in recvBuf.
!
!      recvBuf      recv buffer
!      nrecv        length of receive buffer
!      irecv        node to receive from
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(inout) :: recvBuf(:)
     integer,                intent(inout) :: nrecv
     integer,                intent(in)    :: irecv
     integer,                intent(in)    :: tag
     integer,                intent(in)    :: comm
     integer,                intent(out)   :: request

!    local variables
     integer :: ierror

     call MPI_Irecv(recvBuf, nrecv, MPI_REAL8, irecv, &
                    tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Irecv Failed")
     endif

     return
  end subroutine mpi_MPIIrecv_r

!-----------------------------------------------------------------------
  subroutine mpi_MPIIrecv_r2(recvBuf,nrecv,irecv,tag,comm,request)

!    This subroutine performs an MPI nonblocking receive operation.
!    The received data, a two-dimensional floating point array, 
!    is returned in recvBuf.
!
!      recvBuf      recv buffer
!      nrecv        length of receive buffer
!      irecv        node to receive from
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     real(double), contiguous, intent(inout) :: recvBuf(:,:)
     integer,                intent(inout) :: nrecv
     integer,                intent(in)    :: irecv
     integer,                intent(in)    :: tag
     integer,                intent(in)    :: comm
     integer,                intent(out)   :: request

!    local variables
     integer :: ierror

     call MPI_Irecv(recvBuf, nrecv, MPI_REAL8, irecv, &
                    tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Irecv Failed")
     endif

     return
  end subroutine mpi_MPIIrecv_r2

!-----------------------------------------------------------------------
  subroutine mpi_MPIIrecv_i(recvBuf,nrecv,irecv,tag,comm,request)


!    This subroutine performs an MPI nonblocking receive operation.
!    The received data, a one-dimensional integer array, 
!    is returned in recvBuf.
!
!      recvBuf      recv buffer
!      nrecv        length of receive buffer
!      irecv        node to receive from
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     integer, contiguous, intent(inout) :: recvBuf(:)
     integer,             intent(inout) :: nrecv
     integer,             intent(in)    :: irecv
     integer,             intent(in)    :: tag
     integer,             intent(in)    :: comm
     integer,             intent(out)   :: request

!    local variables
     integer :: ierror

     call MPI_Irecv(recvBuf, nrecv, MPI_INTEGER, irecv, &
                    tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Irecv Failed")
     endif

     return
  end subroutine mpi_MPIIrecv_i

!-----------------------------------------------------------------------
  subroutine mpi_MPIIrecv_i2(recvBuf,nrecv,irecv,tag,comm,request)


!    This subroutine performs an MPI nonblocking receive operation.
!    The received data, a two-dimensional integer array, 
!    is returned in recvBuf.
!
!      recvBuf      recv buffer
!      nrecv        length of receive buffer
!      irecv        node to receive from
!      tag          message tag
!      comm         MPI communicator
!      request      handle

!    variable declarations
     implicit none

!    passed variables
     integer, contiguous, intent(inout) :: recvBuf(:,:)
     integer,             intent(inout) :: nrecv
     integer,             intent(in)    :: irecv
     integer,             intent(in)    :: tag
     integer,             intent(in)    :: comm
     integer,             intent(out)   :: request

!    local variables
     integer :: ierror

     call MPI_Irecv(recvBuf, nrecv, MPI_INTEGER, irecv, &
                    tag, comm, request, ierror)

     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Irecv Failed")
     endif

     return
  end subroutine mpi_MPIIrecv_i2

!=======================================================================
! MPIStart interface
!=======================================================================

  subroutine mpi_MPIStart(comm)

!    This subroutine starts a send for a communicator

!    variable declarations
     implicit none

!    passed variables
     integer, intent(inout) :: comm

!    local variables
     integer :: ierror

!    MPI Start
     call MPI_Start(comm,ierror)
     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Start Failed")
     endif

     return
  end subroutine mpi_MPIStart

!=======================================================================
! MPIWait interface
!=======================================================================

  subroutine mpi_MPIWait(request)

!    This subroutine waits for a requested send or receive to complete

!    variable declarations
     implicit none

!    passed variables
     integer, intent(inout) :: request

!    local variables
     integer :: status(MPI_STATUS_SIZE)
     integer :: ierror

!    MPI Wait
     call MPI_Wait(request,status,ierror)
     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Wait Failed")
     endif

     return
  end subroutine mpi_MPIWait

!=======================================================================
! MPIWaitall interface
!=======================================================================

  subroutine mpi_MPIWaitall(count, request)

!    This subroutine waits for all requested sends or receives to complete

!    variable declarations
     implicit none

!    passed variables
     integer, intent(in) :: count
     integer, intent(inout) :: request(count)

!    local variables
     integer :: status(MPI_STATUS_SIZE,count)
     integer :: ierror

!    MPI Waitall
     call MPI_Waitall(count,request,status,ierror)
     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Waitall Failed")
     endif

     return
  end subroutine mpi_MPIWaitall

!=======================================================================
! MPIWtime interface
!=======================================================================

  double precision function mpi_MPIWtime() result(MPIWtime)

!    This subroutine returns time in seconds from some arbitrary time

!    variable declarations
     implicit none

!    MPI Wtime
     MPIWtime = MPI_Wtime()

     return
  end function mpi_MPIWtime

!=======================================================================
! getMPIRank interface
!=======================================================================

  function mpi_getMPIRank(comm) result(MPIRank)

!    This subroutine determines the rank of the calling process in the
!    communicator.
!
!      comm      MPI communicator
!      MPIRank   processor rank

!    variable declarations
     implicit none

!    passed variables
     integer, intent(in) :: comm
     integer             :: MPIRank

!    local variables
     integer :: ierror

!    MPI Communicator Rank
     call MPI_Comm_rank(comm, MPIrank, ierror)
     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Barrier Failed")
     endif

     return
  end function mpi_getMPIRank

!=======================================================================
! getMPISize interface
!=======================================================================

  function mpi_getMPISize(comm) result(MPISize)

!    This subroutine determines the size of the group associated with
!    the communicator.
!      comm      MPI communicator
!      MPISize   communicator size

!    variable declarations
     implicit none

!    passed variables
     integer, intent(in) :: comm
     integer             :: MPISize

!    local variables
     integer :: ierror

!    MPI Communicator Size
     call MPI_Comm_size(comm, MPISize, ierror)
     if (ierror /= MPI_SUCCESS) then
        call f90fatal("MPI Barrier Failed")
     endif

     return
  end function mpi_getMPISize

end module mpif90_mod
