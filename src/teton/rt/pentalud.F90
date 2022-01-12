!***********************************************************************

   subroutine pentalud

   use kind_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod
   use GreyAcceleration_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none 

!  Local Variables

   integer         :: j
   integer         :: nCorner
   integer         :: neighborID
   integer         :: outerNeighborID
   integer         :: nMU
   integer         :: request(2)
   integer         :: innerBdyID
   integer         :: outerBdyID

   real(adqt)      :: MUrecv(3,2)
   real(adqt)      :: MUsend(3,2)

!  Constants

   nMU        = 6
   nCorner    = Size% ncornr
   innerBdyID = getInnerBdyID(RadBoundary)
   outerBdyID = getOuterBdyID(RadBoundary)


   if ( innerBdyShared(RadBoundary) ) then
     Bdy        => getBoundary(RadBoundary, innerBdyID)
     neighborID =  getNeighborID(Bdy)

     call MPIIrecv(MUrecv, nMU, neighborID, 800, MY_COMM_GROUP, request(2))
   endif

   if ( outerBdyShared(RadBoundary) ) then
     Bdy             => getBoundary(RadBoundary, outerBdyID)
     outerNeighborID =  getNeighborID(Bdy)
   endif


   if ( innerBdyShared(RadBoundary) ) then

     call MPIWait(request(2))

       GTA% MatrixL(1,1) =  GTA% MatrixA(1,1)/MUrecv(1,1)
       GTA% MatrixL(2,1) = (GTA% MatrixA(2,1) - GTA% MatrixL(1,1)*    &
                            MUrecv(2,1)) / MUrecv(1,2)
       GTA% MatrixU(1,1) =  GTA% MatrixA(3,1) -   &
                            GTA% MatrixL(1,1)*MUrecv(3,1) - &
                            GTA% MatrixL(2,1)*MUrecv(2,2)
       GTA% MatrixU(2,1) =  GTA% MatrixA(4,1) -   &
                            GTA% MatrixL(2,1)*MUrecv(3,2)


       GTA% MatrixL(1,2) =  GTA% MatrixA(1,2)/MUrecv(1,2)
       GTA% MatrixL(2,2) = (GTA% MatrixA(2,2) - GTA% MatrixL(1,2)*    &
                            MUrecv(2,2)) / GTA% MatrixU(1,1) 
       GTA% MatrixU(1,2) =  GTA% MatrixA(3,2) -   &
                            GTA% MatrixL(1,2)*MUrecv(3,2) - &
                            GTA% MatrixL(2,2)*GTA% MatrixU(2,1)
       GTA% MatrixU(2,2) =  GTA% MatrixA(4,2) -   &
                            GTA% MatrixL(2,2)*GTA% MatrixU(3,1)


     do j=3,nCorner
       GTA% MatrixL(1,j) =  GTA% MatrixA(1,j)/GTA% MatrixU(1,j-2)
       GTA% MatrixL(2,j) = (GTA% MatrixA(2,j) - GTA% MatrixL(1,j)*    &
                            GTA% MatrixU(2,j-2)) / GTA% MatrixU(1,j-1)
       GTA% MatrixU(1,j) =  GTA% MatrixA(3,j) -   &
                            GTA% MatrixL(1,j)*GTA% MatrixU(3,j-2) - &
                            GTA% MatrixL(2,j)*GTA% MatrixU(2,j-1)
       GTA% MatrixU(2,j) =  GTA% MatrixA(4,j) -   &
                            GTA% MatrixL(2,j)*GTA% MatrixU(3,j-1)

     enddo

   else

     GTA% MatrixU(1,1) = GTA% MatrixA(3,1)
     GTA% MatrixU(2,1) = GTA% MatrixA(4,1)

     GTA% MatrixL(2,2) = GTA% MatrixA(2,2)/GTA% MatrixU(1,1)
     GTA% MatrixU(1,2) = GTA% MatrixA(3,2) - GTA% MatrixL(2,2)*GTA% MatrixU(2,1)
     GTA% MatrixU(2,2) = GTA% MatrixA(4,2) - GTA% MatrixL(2,2)*GTA% MatrixU(3,1)

     do j=3,nCorner
       GTA% MatrixL(1,j) =  GTA% MatrixA(1,j)/GTA% MatrixU(1,j-2)
       GTA% MatrixL(2,j) = (GTA% MatrixA(2,j) - GTA% MatrixL(1,j)*    &
                            GTA% MatrixU(2,j-2)) / GTA% MatrixU(1,j-1)
       GTA% MatrixU(1,j) =  GTA% MatrixA(3,j) -   &
                            GTA% MatrixL(1,j)*GTA% MatrixU(3,j-2) - &
                            GTA% MatrixL(2,j)*GTA% MatrixU(2,j-1)
       GTA% MatrixU(2,j) =  GTA% MatrixA(4,j) -   &
                            GTA% MatrixL(2,j)*GTA% MatrixU(3,j-1)

     enddo


   endif

   if ( outerBdyShared(RadBoundary) ) then
     MUsend(:,1) = GTA% MatrixU(:,nCorner-1)
     MUsend(:,2) = GTA% MatrixU(:,nCorner)

     call MPIIsend(MUsend, nMU, outerNeighborID, 800, MY_COMM_GROUP, request(1))

     call MPIWait(request(1))
   endif


 
   return
   end subroutine pentalud



