!***********************************************************************
!                        Version 1:  03/2009, PFN                      *
!                                                                      *
!   setIncidentFlux - Calculates the incident flux on shared           *
!                     boundaries.                                      *
!                                                                      *
!***********************************************************************
   subroutine setIncidentFlux(psib) 

   use kind_mod
   use constant_mod
   use Size_mod
   use QuadratureList_mod
   use Quadrature_mod
   use Communicator_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Arguments

   real(adqt), intent(in) :: psib(Size%ngr,Size%nbelem,Size%nangSN) 

!  Local

   integer    :: bin,i,ia,ib,ig
   integer    :: b0,ishared,nShared,Groups,nrec,nrecv0,Angle,NangBin,maxNangBin
   integer, allocatable :: message_base(:)
   type(Boundary), pointer :: Bdary

   real(adqt) :: dot, sumdot, mysum

!  Constants

   nShared =  getNumberOfShared(RadBoundary)
   Groups  =  QuadSet% Groups
   maxNangBin = maxval(QuadSet% NangBinList(:))
   allocate(message_base(maxNangBin))

!  Edit

   QuadSet% IncFluxOld(:) = QuadSet% IncFlux(:)
   QuadSet% IncFlux(:)    = zero
   QuadSet% Flux(:,:)     = zero

   SharedBoundary: do ishared=1,nShared
! YKT: use a local pointer Bdary instead of Bdy because Bdy is thread-private
     Bdary => getShared(RadBoundary, ishared)
     b0  =  getFirstBdyElement(Bdy) - 1

     AngleBin: do bin=1,QuadSet% NumBin0
       Comm => getMessage(QuadSet, bin, ishared)
       NangBin = QuadSet% NangBinList(bin)

       nrecv0 = 0
       do ia=1,NangBin
         message_base(ia) = nrecv0
         nrec = Comm% nrecv(ia)
         nrecv0 = nrecv0 + nrec
       end do

       sumdot = zero
       if (Comm% lenrecv > 0) then

         sumdot = zero
!$omp parallel do private(ia,Angle,nrec,mysum,i,ib,dot,ig) reduction(+:sumdot)
         do ia=1,NangBin
           Angle = QuadSet% AngleOrder(ia,bin)
           nrec = Comm% nrecv(ia)
           mysum = zero
           do i=1,nrec
             ib = Comm% ListRecv(message_base(ia)+i)
             dot = QuadSet%omega(1,Angle)*Bdary%A_bdy(1,ib-b0)  &
                 + QuadSet%omega(2,Angle)*Bdary%A_bdy(2,ib-b0)  &
                 + QuadSet%omega(3,Angle)*Bdary%A_bdy(3,ib-b0)  
             do ig=1,Groups
               mysum = mysum - dot*psib(ig,ib,Angle)
             enddo
           enddo
           sumdot = sumdot + mysum
         enddo

       endif

       QuadSet% IncFlux(bin)      = sumdot + QuadSet% IncFlux(bin) 
       QuadSet% Flux(bin,ishared) = sumdot 

     enddo AngleBin

   enddo SharedBoundary 

   deallocate(message_base)
 
   return
   end subroutine setIncidentFlux 

