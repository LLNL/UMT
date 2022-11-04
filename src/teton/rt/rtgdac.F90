!***********************************************************************
!                        Version 1:  05/92, PFN                        *
!                                                                      *
!   RTGDAC - calculates coefficients for the cell-centered grey-       *
!            diffusion acceleration equations and solves the resulting *
!            pentadiagonal matrix for the grey corrections             *
!                                                                      *
!   Output:  MATRIXL   - diagonals of the lower triangular matrix      *
!            MATRIXU   - diagonals of the upper triangular matrix      *
!                                                                      *
!***********************************************************************
   subroutine rtgdac

   use kind_mod
   use constant_mod
   use mpi_param_mod
   use mpif90_mod
   use Size_mod
   use Geometry_mod
   use GreyAcceleration_mod
   use BoundaryList_mod
   use Boundary_mod

   implicit none

!  Local

   integer          :: zone, nzones, ncornr, nGeom
   integer          :: innerBdyID
   integer          :: outerBdyID 
   integer          :: neighborID

   integer          :: request(4)

   real(adqt)       :: D, Dm1, Dp1, Dp2
   real(adqt)       :: sigvol1, sigvol2
   real(adqt)       :: hinv, hinvm1, hinvp1
   real(adqt)       :: recvGeom(2,2)
   real(adqt)       :: sendGeom(2,2)

!  Constants

   real(adqt), parameter :: rho = 0.2886751_adqt


   ncornr     = Size%ncornr
   nzones     = Size%nzones
   nGeom      = 2
   innerBdyID = getInnerBdyID(RadBoundary)
   outerBdyID = getOuterBdyID(RadBoundary)

!  Inner Boundary

   if ( innerBdyShared(RadBoundary) ) then
     Bdy        => getBoundary(RadBoundary, innerBdyID)
     neighborID =  getNeighborID(Bdy)

     call MPIIrecv(recvGeom(1:nGeom,innerBdyID), nGeom, neighborID, 800,  &
                   MY_COMM_GROUP, request(2*innerBdyID))


     sendGeom(1,innerBdyID) = GTA% GreyDiffCoef(1)
     sendGeom(2,innerBdyID) = one/Geom% zoneWidth(1)

     call MPIIsend(sendGeom(1:nGeom,innerBdyID), nGeom, neighborID, 800,   &
                   MY_COMM_GROUP, request(2*innerBdyID-1))

     call MPIWait(request(2*innerBdyID-1))

   endif

!  Outer Boundary

   if ( outerBdyShared(RadBoundary) ) then

     Bdy        => getBoundary(RadBoundary, outerBdyID)
     neighborID =  getNeighborID(Bdy)

     call MPIIrecv(recvGeom(1:nGeom,outerBdyID), nGeom, neighborID, 800,  &
                   MY_COMM_GROUP, request(2*outerBdyID))


     sendGeom(1,outerBdyID) = GTA% GreyDiffCoef(ncornr)
     sendGeom(2,outerBdyID) = one/Geom% zoneWidth(nzones)

     call MPIIsend(sendGeom(1:nGeom,outerBdyID), nGeom, neighborID, 800,   &
                   MY_COMM_GROUP, request(2*outerBdyID-1))

     call MPIWait(request(2*outerBdyID-1))

   endif

!    The left half-cell equation for first zone is modified depending 
!    on the boundary condition.

   D       = GTA% GreyDiffCoef(1)
   Dp1     = GTA% GreyDiffCoef(2)
   Dp2     = GTA% GreyDiffCoef(3)
   sigvol1 = GTA% GreySigEff(1)*Geom% Volume(1)
   sigvol2 = GTA% GreySigEff(2)*Geom% Volume(2)
   hinv    = one/Geom% zoneWidth(1)
   hinvp1  = one/Geom% zoneWidth(2)

   if ( innerBdyReflect(RadBoundary) ) then

     Dm1    = D
     hinvm1 = hinv

     GTA% MatrixA(1,1) = -Dm1*Geom% Rmin(1)*hinvm1
     GTA% MatrixA(2,1) = -two*rho*Geom% Rmin(1) + (Dm1*Geom% Rmin(1))*hinvm1
     GTA% MatrixA(3,1) = (Geom% Rave(1)*(D + Dp1) - D*Geom% Rmin(1))*hinv +  &
                                   two*rho*Geom% Rmin(1) + two*sigvol1
     GTA% MatrixA(4,1) = (Geom% Rmin(1)*D - Geom% Rave(1)*(D + Dp1))*hinv
     GTA% MatrixU(3,1) =  zero

   elseif ( innerBdyShared(RadBoundary) ) then

     call MPIWait( request(2*innerBdyID) )

     Dm1    = recvGeom(1,innerBdyID)
     hinvm1 = recvGeom(2,innerBdyID)

     GTA% MatrixA(1,1) = -Dm1*Geom% Rmin(1)*hinvm1
     GTA% MatrixA(2,1) = -two*rho*Geom% Rmin(1) + (Dm1*Geom% Rmin(1))*hinvm1
     GTA% MatrixA(3,1) = (Geom% Rave(1)*(D + Dp1) - D*Geom% Rmin(1))*hinv +  &
                                 two*rho*Geom% Rmin(1) + two*sigvol1
     GTA% MatrixA(4,1) = (Geom% Rmin(1)*D - Geom% Rave(1)*(D + Dp1))*hinv
     GTA% MatrixU(3,1) =  zero
 
   else 
 
     GTA% MatrixA(1,1) =  zero
     GTA% MatrixA(2,1) =  zero
     GTA% MatrixA(3,1) =  Geom% Rave(1)*(D + Dp1)*hinv + two*sigvol1 
     GTA% MatrixA(4,1) = -Geom% Rave(1)*(Dp1 + D)*hinv
     GTA% MatrixU(3,1) =  zero
 
   endif
 
   GTA% MatrixA(1,2) =  zero
   GTA% MatrixA(2,2) = (Dp1*Geom% Rmax(1) - Geom% Rave(1)*(D + Dp1))*hinv
   GTA% MatrixA(3,2) = (Geom% Rave(1)*(D + Dp1) - Dp1*Geom% Rmax(1))*hinv +  &
                        two*rho*Geom% Rmax(1) + two*sigvol2
   GTA% MatrixA(4,2) =  Geom% Rmax(1)*(Dp2*hinvp1 - two*rho)
   GTA% MatrixU(3,2) = -Dp2*Geom% Rmax(1)*hinvp1
 
!  Loop over interior zones
 
   hinvm1 = hinv

   do zone=2,nzones-1

     Dm1     = GTA%GreyDiffCoef(2*zone-2)
     D       = GTA%GreyDiffCoef(2*zone-1)
     Dp1     = GTA%GreyDiffCoef(2*zone)
     Dp2     = GTA%GreyDiffCoef(2*zone+1)
     sigvol1 = GTA%GreySigEff(2*zone-1)*Geom% Volume(2*zone-1)
     sigvol2 = GTA%GreySigEff(2*zone)*Geom% Volume(2*zone)

     hinv    = one/Geom% zoneWidth(zone)
     hinvp1  = one/Geom% zoneWidth(zone+1)
 
     GTA% MatrixA(1,2*zone-1) = -Dm1*Geom% Rmin(zone)*hinvm1
     GTA% MatrixA(2,2*zone-1) = -two*rho*Geom% Rmin(zone) + (Dm1*Geom% Rmin(zone))*hinvm1
     GTA% MatrixA(3,2*zone-1) = (Geom% Rave(zone)*(D + Dp1) - D*Geom% Rmin(zone))*hinv +  &
                                 two*rho*Geom% Rmin(zone) + two*sigvol1 
     GTA% MatrixA(4,2*zone-1) = (Geom% Rmin(zone)*D - Geom% Rave(zone)*(D + Dp1))*hinv 
     GTA% MatrixU(3,2*zone-1) =  zero
 
     GTA% MatrixA(1,2*zone)   =  zero
     GTA% MatrixA(2,2*zone)   = (-Geom% Rave(zone)*(D + Dp1) + Dp1*Geom% Rmax(zone))*hinv
     GTA% MatrixA(3,2*zone)   = ( Geom% Rave(zone)*(D + Dp1) - Dp1*Geom% Rmax(zone))*hinv +  & 
                                  two*rho*Geom% Rmax(zone) + two*sigvol2 
     GTA% MatrixA(4,2*zone)   = -two*rho*Geom% Rmax(zone) + (Dp2*Geom% Rmax(zone))*hinvp1
     GTA% MatrixU(3,2*zone)   = -Dp2*Geom% Rmax(zone)*hinvp1

     hinvm1 = hinv
 
   enddo
 
!  Outer Boundary

   Dm1     = GTA%GreyDiffCoef(ncornr-2)
   D       = GTA%GreyDiffCoef(ncornr-1)
   Dp1     = GTA%GreyDiffCoef(ncornr)
   sigvol1 = GTA%GreySigEff(ncornr-1)*Geom% Volume(ncornr-1)
   sigvol2 = GTA%GreySigEff(ncornr)*Geom% Volume(ncornr)
   hinv    = one/Geom% zoneWidth(nzones)
 
   GTA% MatrixA(1,ncornr-1) = -Dm1*Geom% Rmin(nzones)*hinvm1
   GTA% MatrixA(2,ncornr-1) =  Geom% Rmin(nzones)*(-two*rho + (Dm1*hinvm1))
   GTA% MatrixA(3,ncornr-1) = (Geom% Rave(nzones)*(D + Dp1) -  &
                               D*Geom% Rmin(nzones))*hinv +             &
                               two*rho*Geom% Rmin(nzones) + two*sigvol1 
   GTA% MatrixA(4,ncornr-1) = (-Geom% Rave(nzones)*(D + Dp1) + &
                                D*Geom% Rmin(nzones))*hinv
   GTA% MatrixU(3,ncornr-1) = zero


   if ( outerBdyShared(RadBoundary) ) then

     call MPIWait( request(2*outerBdyID) )

     Dp2    = recvGeom(1,outerBdyID)
     hinvp1 = recvGeom(2,outerBdyID)

     GTA% MatrixA(1,ncornr)   =  zero
     GTA% MatrixA(2,ncornr)   = (-Geom% Rave(nzones)*(D + Dp1) +  &
                                 Dp1*Geom%Rmax(nzones))*hinv
     GTA% MatrixA(3,ncornr)   = ( Geom% Rave(nzones)*(D + Dp1) -  &
                                 Dp1*Geom% Rmax(nzones))*hinv +  &
                                   two*rho*Geom% Rmax(nzones) + two*sigvol2
     GTA% MatrixA(4,ncornr)   = -two*rho*Geom% Rmax(nzones) +  &
                                (Dp2*Geom% Rmax(nzones))*hinvp1
     GTA% MatrixU(3,ncornr)   = -Dp2*Geom% Rmax(nzones)*hinvp1

   elseif ( outerBdyReflect(RadBoundary) ) then

     Dp2    = Dp1
     hinvp1 = hinv

     GTA% MatrixA(1,ncornr)   =  zero
     GTA% MatrixA(2,ncornr)   = (-Geom% Rave(nzones)*(D + Dp1) +  &
                                  Dp1*Geom% Rmax(nzones))*hinv
     GTA% MatrixA(3,ncornr)   = ( Geom% Rave(nzones)*(D + Dp1) -  &
                                  Dp1*Geom% Rmax(nzones))*hinv +  &
                                  two*rho*Geom% Rmax(nzones) + two*sigvol2
     GTA% MatrixA(4,ncornr)   = -two*rho*Geom% Rmax(nzones) +  &
                                (Dp2*Geom% Rmax(nzones))*hinvp1
     GTA% MatrixU(3,ncornr)   = -Dp2*Geom% Rmax(nzones)*hinvp1

   else

     GTA% MatrixA(1,ncornr)   = zero
     GTA% MatrixA(2,ncornr)   = (-Geom% Rave(nzones)*(D + Dp1) +  &
                                  Dp1*Geom% Rmax(nzones))*hinv
     GTA% MatrixA(3,ncornr)   = ( Geom% Rave(nzones)*(D + Dp1) -  &
                                  Dp1*Geom% Rmax(nzones))*hinv + &
                                  two*rho*Geom% Rmax(nzones) + two*sigvol2
     GTA% MatrixA(4,ncornr)   = zero
     GTA% MatrixU(3,ncornr)   = zero

   endif
 
 
!  LU decompose pentadiagonal matrix
 
   call pentalud

 
   return
   end subroutine rtgdac

