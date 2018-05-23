!***********************************************************************
!                        Version 1:  09/96, PFN                        *
!                                                                      *
!   RTORDER - This routine builds an ordered list of corners for each  *
!             unique direction.                                        *
!                                                                      *
!   Input:                                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!***********************************************************************
   subroutine rtorder 

   use kind_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use Quadrature_mod
   use BoundaryList_mod
   use Boundary_mod
   use cudafor
   use nvtx_mod

   implicit none

!  Local Variables

   integer    :: nBoundary, nBdyElem, b0

   integer    :: i,ii,ia,ib,ic,iz,istat
   integer    :: ncornr,nbelem,ndim,nExit
   integer    :: NumQuadSets

   real(adqt) :: dot

!  Dynamic

   integer,  allocatable :: ListExit(:,:)

!  Constants

   ncornr      = Size%ncornr
   ndim        = Size%ndim
   nbelem      = Size%nbelem
   NumQuadSets = getNumQuadSets(Quad)
   nBoundary   = getNumberOfBoundaries(RadBoundary)


   allocate( ListExit(2,nbelem) )

!  Build NEXT for all angle sets 

   AngleSetLoop: do i=1,NumQuadSets

     QuadSet => getQuadrature(Quad, i)

!  Initialize total number of cycles, cycle counter and offsets
!  per angle

     QuadSet% totcycles = 0

!$OMP PARALLEL DO  PRIVATE(ia)
     OMPAngleLoop: do ia=1,QuadSet% NumAngles
        call nvtxStartRange("snnext",3)
        call snnext(ia) 
        call nvtxEndRange

     enddo OMPAngleLoop

     !!! synchronous, but easy syntax
     !QuadSet%d_next = QuadSet%next
     !QuadSet%d_nextZ = QuadSet%nextZ
     !QuadSet%d_passZstart = QuadSet%passZstart

     !!! asynchronous version
     istat = cudaMemcpyAsync(QuadSet%d_next, QuadSet%next, (Size%ncornr+1)*QuadSet%NumAngles, 0)
     istat = cudaMemcpyAsync(QuadSet%d_nextZ, QuadSet%nextZ,Size%nzones*QuadSet%NumAngles, 0)
     istat = cudaMemcpyAsync(QuadSet%d_passZstart, QuadSet%passZstart, Size%nzones*QuadSet%NumAngles, 0)

     AngleLoop: do ia=1,QuadSet% NumAngles

!  Make a list of exiting angles on all boundaries

       nExit = 0

       BoundaryLoop: do ii=1,nBoundary
         Bdy      => getBoundary(RadBoundary, ii)
         nBdyElem =  getNumberOfBdyElements(Bdy)
         b0       =  getFirstBdyElement(Bdy) - 1
         do ib=1,nBdyElem
           ic = Bdy% BdyToC(ib)

           dot = DOT_PRODUCT( QuadSet%omega(:,ia),Bdy%A_bdy(:,ib) )

           if (dot > zero) then
             nExit             = nExit + 1
             ListExit(1,nExit) = ib + b0 
             ListExit(2,nExit) = ic
           endif
         enddo

       enddo BoundaryLoop

       call constructExitList(QuadSet, ia, nExit, ListExit)


     enddo AngleLoop


   enddo AngleSetLoop


   deallocate( ListExit )


   return
   end subroutine rtorder


