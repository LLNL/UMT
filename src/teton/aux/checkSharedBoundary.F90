!***********************************************************************
!                        Version 1:  02/2011, PFN                      *
!                                                                      *
!   checkSharedBoundary - This routine checks that faces on shared     *
!                         boundaries are properly aligned.             *
!                                                                      *
!***********************************************************************

   subroutine checkSharedBoundary() BIND(C,NAME="teton_checksharedboundary")

   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use Geometry_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Local 

   integer    :: request(2*Size%ncomm)
   integer    :: d, i, ib
   integer    :: c, ndim
   integer    :: nBdyElem
   integer    :: nShared
   integer    :: neighbor
   integer    :: nSend
   integer    :: nRecv

   real(adqt) :: error(Size% ndim)
   real(adqt) :: maxError

   real(adqt), parameter :: tolerance = 1.0d-12

!  Check all of the shared boundaries to insure that all faces
!  and nodes are properly aligned

   nShared = getNumberOfShared(RadBoundary)
   ndim    = Size% ndim

   do i=1,nShared
     Bdy      => getShared(RadBoundary, i)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     neighbor =  getNeighborID(Bdy)

     allocate( Bdy% nodePosition(ndim,nBdyElem) )
     allocate( Bdy% nodePosition2(ndim,nBdyElem) )

     nRecv = ndim*nBdyElem

!  Post Receives

     call MPIIrecv(Bdy% nodePosition2(:,1:nBdyElem), nRecv,   &
                   neighbor, 888, MY_COMM_GROUP, request(2*i))

   enddo

   do i=1,nShared
     Bdy      => getShared(RadBoundary, i)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     neighbor =  getNeighborID(Bdy)

     do ib=1,nBdyElem
       c  = Bdy% BdyToC(ib)

       Bdy% nodePosition(:,ib) = Geom% px(:,c)
     enddo

     nSend = ndim*nBdyElem

!  Send my information to all neighbors

     call MPIIsend(Bdy% nodePosition(:,1:nBdyElem), nSend,      &
                   neighbor, 888, MY_COMM_GROUP, request(2*i-1))

   enddo

!  Verify completion of sends

   do i=1,nShared
     call MPIWait(request(2*i-1))
   enddo

!  Alignment must be perfect!

   do i=1,nShared
     Bdy      => getShared(RadBoundary, i)
     nBdyElem =  getNumberOfBdyElements(Bdy)
     neighbor =  getNeighborID(Bdy)
     nSend    =  ndim*nBdyElem

     call MPIWait(request(2*i))

     do ib=1,nBdyElem

       error(:) = abs( Bdy% nodePosition(:,ib) - Bdy% nodePosition2(:,ib) )
       maxError = maxval( error(1:ndim) )

       if ( maxError > tolerance ) then 

         write(6,100) Size% myRankInGroup,neighbor,maxError,tolerance

         do d=1,ndim
           write(6,200) d,Bdy% nodePosition(d,ib),Bdy% nodePosition2(d,ib)
         enddo

         call f90fatal("Boundary Elements are misaligned!!!")

       endif
     enddo

     deallocate( Bdy% nodePosition  )
     deallocate( Bdy% nodePosition2 )

   enddo


 100 format("myMeshDomainID = ",i8,2x,"neighborMeshDomainID = ",i8,2x, &
            "maxErr = ",1pe15.8,2x,"tol = ",1pe15.8)
 200 format("dim = ",i2,2x,"my position = ",1pe22.15,2x,  &
            "neighbor position = ",1pe22.15)


   return
   end subroutine checkSharedBoundary 


