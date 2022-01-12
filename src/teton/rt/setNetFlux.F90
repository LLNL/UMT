!***********************************************************************
!                        Version 1:  03/2009, PFN                      *
!                                                                      *
!   setNetFlux - Calculates the net flux on shared boundaries which    *
!                is used to schedule the transport sweeps (see         *
!                SweepScheduler)                                       *
!                                                                      *
!***********************************************************************
   subroutine setNetFlux(cSetID) 

   use kind_mod
   use constant_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use QuadratureList_mod
   use Communicator_mod
   use BoundaryList_mod
   use Boundary_mod
   use SetData_mod
   use CommSet_mod
   use AngleSet_mod

   implicit none

!  Arguments

   integer,   intent(in)    :: cSetID

!  Local

   type(SetData),         pointer :: Set
   type(CommSet),         pointer :: CSet
   type(AngleSet),        pointer :: ASet
   type(Communicator),    pointer :: CommT
   type(CommunicateFlux), pointer :: CommX
   type(Boundary),        pointer :: BdyT

   integer    :: bin
   integer    :: i
   integer    :: b
   integer    :: b0
   integer    :: g
   integer    :: setID
   integer    :: sharedID
   integer    :: nShared
   integer    :: Angle
   integer    :: NumBin
   integer    :: NumAngles
   integer    :: COMM_GROUP

   real(adqt) :: dot
   real(adqt) :: sumExitFlux

!  Constants

   CSet       => getCommSetData(Quad, cSetID)
   ASet       => CSet% AngleSetPtr

   nShared    =  getNumberOfShared(RadBoundary)
   COMM_GROUP =  getCommunicationGroup(CSet)
   NumAngles  =  CSet% NumAngles
   NumBin     =  CSet% NumBin0 

!  Post Receives

   do sharedID=1,nShared
     CommX => getSharedFlux(CSet, sharedID)

     call MPIStart(CommX% irequest(2))
   enddo

!  Compute exiting flux

   SharedBoundary: do sharedID=1,nShared
     CommX => getSharedFlux(CSet, sharedID)
     BdyT  => getShared(RadBoundary, sharedID)
     b0    =  getFirstBdyElement(BdyT) - 1

     CommX% ExitFlux(:) = zero

     AngleLoop: do Angle=1,NumAngles
       CommT  => getMessage(CSet, sharedID, Angle)
       bin    =  CSet% AngleToBin(Angle)

       sumExitFlux = zero

       if (CommT% nSend > 0) then

         do setID=CSet% set1,CSet% set2
           Set => getSetData(Quad, setID)

           do i=1,CommT% nSend
             b   = CommT% ListSend(1,i)
             dot = DOT_PRODUCT( ASet%omega(:,Angle),BdyT%A_bdy(:,b-b0) )

             do g=1,Set% Groups
               sumExitFlux = sumExitFlux + dot*Set% Psib(g,b,Angle)
             enddo
           enddo

         enddo

       endif

       CommX% ExitFlux(bin) = CommX% ExitFlux(bin) +  &
                              ASet%weight(Angle)*sumExitFlux

     enddo AngleLoop

   enddo SharedBoundary

!  Send exiting flux to all neighbors

   do sharedID=1,nShared
     CommX => getSharedFlux(CSet, sharedID)

     call MPIStart(CommX% irequest(1))
   enddo

!  Check for completion of Sends

   do sharedID=1,nShared
     CommX => getSharedFlux(CSet, sharedID)

     call MPIWait(CommX% irequest(1))
   enddo

!  Update the net flux across a shared boundary

   do sharedID=1,nShared
     CommX => getSharedFlux(CSet, sharedID)
     call MPIWait(CommX% irequest(2))

     do bin=1,NumBin
       CSet% NetFlux(sharedID,bin) = CommX% ExitFlux(bin) - CommX% IncFlux(bin)
     enddo
   enddo



   return
   end subroutine setNetFlux 

