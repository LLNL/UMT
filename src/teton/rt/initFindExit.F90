!***********************************************************************
!                        Version 1:  11/98, PFN                        *
!                                                                      *
!   initFindExit - Initializes persistent communication objects for    *
!                  creating exit lists on shared boundaries.           *
!                                                                      *
!***********************************************************************

   subroutine initFindExit(nAngleSets, nGTASets)

   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use Boundary_mod
   use AngleSet_mod
   use CommSet_mod

   implicit none

!  Arguments

   integer,               intent(in) :: nAngleSets
   integer,               intent(in) :: nGTASets 

!  Local

   type(CommSet),         pointer    :: CSet
   Type(AngleSet),        pointer    :: ASet
   type(IncidentTest),    pointer    :: TSet
   type(Boundary),        pointer    :: BdyT

   integer :: aSetID
   integer :: sharedID
   integer :: nShared
   integer :: neighborRank
   integer :: tag
   integer :: tagBase
   integer :: numBdyElem
   integer :: numAngles

   integer :: COMM_GROUP
   integer :: myRankInGroup

!  Constants

   myRankInGroup = Size% myRankInGroup
   nShared       = getNumberOfShared(RadBoundary) 

   SharedLoop: do sharedID=1,nShared

     Bdy          => getShared(RadBoundary, sharedID)
     numBdyElem   =  getNumberOfBdyElements(Bdy)
     neighborRank =  getNeighborID(Bdy)

     AngleSetLoop: do aSetID=1,nAngleSets+nGTASets
       ASet => getAngleSetData(Quad, aSetID)

       call constructIncidentTest(ASet, sharedID, numBdyElem, neighborRank)
     enddo AngleSetLoop

   enddo SharedLoop

!  Each message depends on angleSet and shared boundary 

   AngleSetLoop2: do aSetID=1,nAngleSets+nGTASets

     ASet       => getAngleSetData(Quad, aSetID)
     CSet       => getCommSetData(Quad, aSetID) 
     COMM_GROUP =  CSet% COMM_GROUP
     tagBase    =  4000 + 100*aSetID
     numAngles  =  ASet% numAngles

!  Loop over the number of shared surfaces

     SharedLoop2: do sharedID=1,nShared

       TSet         => getIncidentTest(ASet, sharedID)
       BdyT         => getShared(RadBoundary, sharedID)
       neighborRank =  getNeighborID(BdyT)
       numBdyElem   =  getNumberOfBdyElements(BdyT)

!  Initialize Persistant Communication Objects (receives are "even"
!  requests and sends are "odd" requests)

       tag    = tagBase + myRankInGroup + neighborRank 

       call MPIRecvInit1(TSet% IncTestR, TSet% nRecv, neighborRank,  &
                         tag, COMM_GROUP, TSet% request(2))

       call MPISendInit1(TSet% IncTest, TSet% nSend, neighborRank,   &
                         tag, COMM_GROUP, TSet% request(1)) 

     enddo SharedLoop2

   enddo AngleSetLoop2

   return
   end subroutine initFindExit

