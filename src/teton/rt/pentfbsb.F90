!***********************************************************************

   subroutine pentfbsb

   use kind_mod
   use constant_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod
   use BoundaryList_mod
   use Boundary_mod
   use GreyAcceleration_mod
   use iter_control_list_mod
   use iter_control_mod
   use Geometry_mod

   implicit none

!  Local Variables

   type(IterControl),      pointer  :: greyControl => NULL()
   
   integer     :: j
   integer     :: n
   integer     :: nZones
   integer     :: innerNeighborID
   integer     :: outerNeighborID
   integer     :: request(4)
   integer     :: nGS
   integer     :: nGC
   integer     :: nGDAiters
   integer     :: innerBdyID
   integer     :: outerBdyID

   real(adqt)  :: GSrecv(2)
   real(adqt)  :: GSsend(2)
   real(adqt)  :: GCrecv(2)
   real(adqt)  :: GCsend(2)

!***********************************************************************
!
!     This routine performs forward & back substitution to solve the
!     matrix equation LUp = s.  Here L is a lower-triangular matrix of
!     dimension n, U is upper triangular dimension n, s is a given
!     source vector and p is the solution.  L is assumed to have 2
!     nonzero diagonals just below its main diagonal;  U is assumed to
!     have 2 nonzero diagonal just above its main diagonal.
!
!     It is assumed that the main diagonal of L is 1.
!
!     n = dimension of system
!     v = given lowermost diagonal of L
!     w = given next diagonal of L
!     x = given main diagonal of U
!     y = given next diagonal of U
!     z = given uppermost diagonal of U
!     s = given source
!     p = solution
!
!***********************************************************************

   nGS        = 2 
   nGC        = 2
   n          = Size% ncornr
   nZones     = Size% nzones

   innerBdyID = getInnerBdyID(RadBoundary)
   outerBdyID = getOuterBdyID(RadBoundary)

   if ( innerBdyShared(RadBoundary) ) then
     Bdy        => getBoundary(RadBoundary, innerBdyID)
     innerNeighborID =  getNeighborID(Bdy)

     call MPIIrecv(GSrecv, nGS, innerNeighborID, 850, MY_COMM_GROUP, request(2))
   endif

   if ( outerBdyShared(RadBoundary) ) then
     Bdy             => getBoundary(RadBoundary, outerBdyID)
     outerNeighborID =  getNeighborID(Bdy)

     call MPIIrecv(GCrecv, nGC, outerNeighborID, 860, MY_COMM_GROUP, request(4))
   endif



   greyControl => getIterationControl(IterControls, "grey")
   nGDAiters   =  getNumberOfIterations(greyControl)

!  Modify the residual source for 1D

   GTA%GreySource(:) = two*Geom% Volume(:)*GTA%GreySource(:)

!  First forward substitution to get the solution to Lt=s.
!  We use the memory of S to store T.

   if ( innerBdyShared(RadBoundary) ) then

     call MPIWait(request(2))

     GTA% GreySource(1) = GTA% GreySource(1) -   &
                          GSrecv(2)*GTA% MatrixL(2,1) -  &
                          GSrecv(1)*GTA% MatrixL(1,1)

     GTA% GreySource(2) = GTA% GreySource(2) -   &
                          GTA% GreySource(1)*GTA% MatrixL(2,2) -  &
                          GSrecv(2)*GTA% MatrixL(1,2)

     do j=3,n
       GTA% GreySource(j) = GTA% GreySource(j) -   &
                            GTA% GreySource(j-1)*GTA% MatrixL(2,j) -  &
                            GTA% GreySource(j-2)*GTA% MatrixL(1,j)
     enddo

   else
 
     GTA% GreySource(2) = GTA% GreySource(2) -   &
                          GTA% GreySource(1)*GTA% MatrixL(2,2)
 
     do j=3,n
       GTA% GreySource(j) = GTA% GreySource(j) -   &
                            GTA% GreySource(j-1)*GTA% MatrixL(2,j) -  &
                            GTA% GreySource(j-2)*GTA% MatrixL(1,j)
     enddo
 
   endif


   if ( outerBdyShared(RadBoundary) ) then
     GSsend(1) = GTA% GreySource(n-1)
     GSsend(2) = GTA% GreySource(n)

     call MPIIsend(GSsend, nGS, outerNeighborID, 850, MY_COMM_GROUP, request(1))

     call MPIWait(request(1))
   endif
 
!  Now back substitution to get the solution of Up=s:
 
   if ( outerBdyShared(RadBoundary) ) then

     call MPIWait(request(4))

     GTA% GreyCorrection(n) = (GTA% GreySource(n) -   &
                               GTA% MatrixU(2,n)*GCrecv(1) - &
                               GTA% MatrixU(3,n)*GCrecv(2))/ &
                               GTA% MatrixU(1,n)

     GTA% GreyCorrection(n-1) = (GTA% GreySource(n-1) -   &
                               GTA% MatrixU(2,n-1)*GTA% GreyCorrection(n) - &
                               GTA% MatrixU(3,n-1)*GCrecv(1))/ &
                               GTA% MatrixU(1,n-1)

     do j=n-2,1,-1
       GTA% GreyCorrection(j) = (GTA% GreySource(j) -   &
                                 GTA% MatrixU(2,j)*GTA% GreyCorrection(j+1) - &
                                 GTA% MatrixU(3,j)*GTA% GreyCorrection(j+2))/ &
                                 GTA% MatrixU(1,j)
     enddo

   else

     GTA% GreyCorrection(n)   = GTA% GreySource(n) / GTA% MatrixU(1,n)
 
     GTA% GreyCorrection(n-1) = (GTA% GreySource(n-1) -   &
                                 GTA% MatrixU(2,n-1)*GTA% GreyCorrection(n)) /  &
                                 GTA% MatrixU(1,n-1)
 
     do j=n-2,1,-1
       GTA% GreyCorrection(j) = (GTA% GreySource(j) -   &
                                 GTA% MatrixU(2,j)*GTA% GreyCorrection(j+1) -   &
                                 GTA% MatrixU(3,j)*GTA% GreyCorrection(j+2))/ &
                                 GTA% MatrixU(1,j)
     enddo
 
   endif

   if ( innerBdyShared(RadBoundary) ) then
     GCsend(1) = GTA% GreyCorrection(1)
     GCsend(2) = GTA% GreyCorrection(2)

     call MPIIsend(GCsend, nGC, innerNeighborID, 860, MY_COMM_GROUP, request(3))

     call MPIWait(request(3))
   endif


   nGDAiters = nGDAiters + 1
   call setNumberOfIterations(greyControl, nGDAiters, .FALSE.)
 

   return
   end subroutine pentfbsb



