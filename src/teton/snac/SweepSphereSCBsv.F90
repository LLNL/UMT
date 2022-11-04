!***********************************************************************
!                       Last Update:  10/2016, PFN                     *
!                                                                      *
!   SweepSphereSCB - Transport sweep for 1D spherical or cylindrical   *
!                    geometry with the Simple Corner Balance spatial   *
!                    discretization.                                   *
!                                                                      *
!***********************************************************************
   subroutine SweepSphereSCBsv(setID, Groups, NumAngles) 

   use flags_mod
   use kind_mod
   use constant_mod
   use Size_mod
   use mpi_param_mod
   use mpif90_mod
   use Communicator_mod
   use Geometry_mod
   use QuadratureList_mod
   use Quadrature_mod
   use BoundaryList_mod
   use Boundary_mod
   use SetData_mod
   use CommSet_mod
   use AngleSet_mod
   use GroupSet_mod

   implicit none

!  Arguments

   integer, intent(in)       :: setID
   integer, intent(in)       :: Groups
   integer, intent(in)       :: NumAngles

!  Local

   type(SetData),          pointer  :: Set
   type(CommSet),          pointer  :: CSet
   type(AngleSet),         pointer  :: ASet
   type(GroupSet),         pointer  :: GSet
   type(Boundary),         pointer  :: BdyT
   type(Communicator),     pointer  :: CommT

   integer    :: Angle,g,zone,nzones
   integer    :: c0
   integer    :: xiLevel
   integer    :: mstdr
   integer    :: nAngLevel
   integer    :: nStartDir
   integer    :: nAngLevelNeg
   integer    :: mNeg1, mNeg2, mPos1, mPos2
   integer    :: angleRef
   integer    :: bIn, bOut, innerBdyID, outerBdyID
   integer    :: QuadID
   integer    :: innerSharedID, outerSharedID

   real(adqt) :: a0, a1, b0, den, Q1, Q2
   real(adqt) :: angTerm1, angTerm2
   real(adqt) :: omega, weight, falpha, adweta, quadtau
   real(adqt) :: tau

   real(adqt) :: psi_inc(Groups, NumAngles)
   real(adqt) :: PsiM1(Groups,2)

   real(adqt), allocatable :: PsiM(:,:)

!  Constants

   Set  => getSetData(Quad, setID)
   CSet => getCommSetFromSetID(Quad, setID)
   ASet => getAngleSetFromSetID(Quad, setID)
   GSet => getGroupSetFromSetID(Quad, setID)

   allocate( PsiM(Groups, Set% nCorner) )

   nzones        = Size% nzones
   tau           = Size% tau
   QuadID        = Set% QuadID

   innerBdyID    = getInnerBdyID(RadBoundary)
   outerBdyID    = getOuterBdyID(RadBoundary)
   innerSharedID = getInnerSharedID(RadBoundary)
   outerSharedID = getOuterSharedID(RadBoundary)

   BdyT => getBoundary(RadBoundary, innerBdyID)
   bIn  =  getFirstBdyElement(BdyT)

   BdyT => getBoundary(RadBoundary, outerBdyID)
   bOut =  getFirstBdyElement(BdyT)

   nStartDir = 0

   if (Size%igeom == geometry_sphere) then
     nStartDir = 1
   elseif (Size%igeom == geometry_cylinder) then
     nStartDir = Set%Order/2
   endif

   mstdr = 1

   XiLevelLoop: do xiLevel=1,nStartDir

   nAngLevel     = Set%Order - 2*(xiLevel-1)
   nAngLevelNeg  = nAngLevel/2

   mNeg1 = mstdr + 1
   mNeg2 = mstdr + nAngLevelNeg
   mPos1 = mNeg2 + 1
   mPos2 = mstdr + nAngLevel

!***************
!    mu < 0    *
!***************

!  Set the boundary fluxes on the outer boundary

   do Angle=1,mNeg2

     if ( outerBdyShared(RadBoundary) ) then
       CommT => getMessage(CSet, outerSharedID, Angle)

       call MPIWait(CommT% irequest(2))

       Set% PsiB(:,bOut,Angle) = CommT% psibrecv(:,1)

     elseif ( outerBdyReflect(RadBoundary) ) then
       BdyT => getBoundary(RadBoundary, outerBdyID)

       angleRef = getReflectedAngle(ASet, outerBdyID, Angle)
       if (angleRef > 0) then
         Set% PsiB(:,bOut,Angle) = Set% PsiB(:,bOut,angleRef)
       endif
     endif

     psi_inc(:,Angle) = Set% PsiB(:,bOut,Angle)

   enddo

!  Perform sweep starting at outer boundary

   ZoneLoop1: do zone=nzones,1,-1

     c0 = 2*(zone - 1)

!  Initialize the scalar intensity

     Set% Phi(:,c0+1) = zero
     Set% Phi(:,c0+2) = zero

!  Starting Direction

     a0 = -half*Geom% Rave(zone)

     do g=1,Groups

       Q1  = Geom% Volume(c0+1)*(GSet% STotal(g,c0+1) + tau*Set% Psi(g,c0+1,mstdr))
       Q2  = Geom% Volume(c0+2)*(GSet% STotal(g,c0+2) + tau*Set% Psi(g,c0+2,mstdr)) +  &
                                 Geom% Rmax(zone)*psi_inc(g,1)

       b0  = -a0 + GSet% Sigt(g,zone)*Geom% Volume(c0+1)
       a1  =  a0 + GSet% Sigt(g,zone)*Geom% Volume(c0+2) + Geom% Rmax(zone)
       den =  one/(a0*a0 + b0*a1)

       Set% Psi(g,c0+1,mstdr) = ( a1*Q1 - a0*Q2 )*den
       Set% Psi(g,c0+2,mstdr) = ( a0*Q1 + b0*Q2 )*den

       psi_inc(g,1)           = Set% Psi(g,c0+1,mstdr)

       PsiM1(g,1)             = Set% Psi(g,c0+1,mstdr)
       PsiM1(g,2)             = Set% Psi(g,c0+2,mstdr)

     enddo

     AngleLoop1: do Angle=mNeg1,mNeg2

!  Angular coefficients

       omega   = ASet% omega(1,Angle)
       weight  = ASet% weight(Angle)
       falpha  = ASet% falpha(Angle)
       adweta  = ASet% adweta(Angle)
       quadtau = ASet% tau(Angle)

       a0       = half*omega*Geom% Rave(zone)
       angTerm1 = adweta*Geom% Area(c0+1) + a0 - Geom% Rmin(zone)*omega
       angTerm2 = adweta*Geom% Area(c0+2) - a0

       do g=1,Groups

         Q1  = Geom% Volume(c0+1)*(GSet% STotal(g,c0+1) + tau*Set% Psi(g,c0+1,Angle)) + &
               falpha*Geom% Area(c0+1)*PsiM1(g,1)
         Q2  = Geom% Volume(c0+2)*(GSet% STotal(g,c0+2) + tau*Set% Psi(g,c0+2,Angle)) + &
               falpha*Geom% Area(c0+2)*PsiM1(g,2)   - &
               omega*Geom% Rmax(zone)*psi_inc(g,Angle)

         b0  =  GSet% Sigt(g,zone)*Geom% Volume(c0+1) + angTerm1 
         a1  =  GSet% Sigt(g,zone)*Geom% Volume(c0+2) + angTerm2 
         den =  one/(a0*a0 + b0*a1)

         Set% Psi(g,c0+1,Angle) = ( a1*Q1 - a0*Q2 )*den
         Set% Psi(g,c0+2,Angle) = ( a0*Q1 + b0*Q2 )*den
         psi_inc(g,Angle)       = Set% Psi(g,c0+1,Angle)

!        Compute (m+1/2) and scalar intensities

         PsiM1(g,1)  = quadtau*Set% Psi(g,c0+1,Angle) +  &
                          (one - quadtau)*PsiM1(g,1)
         PsiM1(g,2)  = quadtau*Set% Psi(g,c0+2,Angle) +  &
                          (one - quadtau)*PsiM1(g,2)

         Set% Phi(g,c0+1) = Set% Phi(g,c0+1) + weight*Set% Psi(g,c0+1,Angle)
         Set% Phi(g,c0+2) = Set% Phi(g,c0+2) + weight*Set% Psi(g,c0+2,Angle)

       enddo

     enddo AngleLoop1

     PsiM(:,c0+1) = PsiM1(:,1)
     PsiM(:,c0+2) = PsiM1(:,2)

   enddo ZoneLoop1

!  Inner Boundary

   do Angle=1,mNeg2
     Set% PsiB(:,bIn,Angle) = psi_inc(:,Angle)
   enddo

   if ( innerBdyShared(RadBoundary) ) then
     BdyT  => getBoundary(RadBoundary, innerBdyID)

     do Angle=1,mNeg2
       CommT => getMessage(CSet, innerSharedID, Angle)

       CommT% psibsend(:,1) = Set% PsiB(:,bIn,Angle) 

       call MPIStart(CommT% irequest(1))
       call MPIWait(CommT% irequest(1))
     enddo

   endif

!***************
!    mu > 0    *
!***************

   BdyT => getBoundary(RadBoundary, innerBdyID)

   if ( innerBdyShared(RadBoundary) ) then

     do Angle=mPos1,mPos2
       CommT => getMessage(CSet, innerSharedID, Angle)

       call MPIWait(CommT% irequest(2))

       Set% PsiB(:,bIn,Angle) = CommT% psibrecv(:,1)
     enddo

   elseif ( innerBdyReflect(RadBoundary) ) then

     do Angle=mPos1,mPos2
       angleRef = getReflectedAngle(ASet, innerBdyID, Angle)
       if (angleRef > 0) then
         Set% PsiB(:,bIn,Angle) = Set% PsiB(:,bIn,angleRef)
       endif
     enddo

   endif

   do Angle=mPos1,mPos2
     psi_inc(:,Angle) = Set% PsiB(:,bIn,Angle)
   enddo


!  Perform sweep starting at inner boundary

   ZoneLoop2: do zone=1,nzones

     c0 = 2*(zone - 1)

     PsiM1(:,1) = PsiM(:,c0+1)
     PsiM1(:,2) = PsiM(:,c0+2)

     AngleLoop2: do Angle=mPos1,mPos2

!  Angular coefficients

       omega   = ASet% omega(1,Angle)
       weight  = ASet% weight(Angle)
       falpha  = ASet% falpha(Angle)
       adweta  = ASet% adweta(Angle)
       quadtau = ASet% tau(Angle)

       a0       = half*omega*Geom% Rave(zone)
       angTerm1 = adweta*geom% Area(c0+1) + a0
       angTerm2 = adweta*Geom% Area(c0+2) + Geom% Rmax(zone)*omega - a0

       do g=1,Groups

         Q1  = Geom% Volume(c0+1)*(GSet% STotal(g,c0+1) + tau*Set% Psi(g,c0+1,Angle)) + &
               falpha*Geom% Area(c0+1)*PsiM1(g,1)   + &
               omega*Geom% Rmin(zone)*psi_inc(g,Angle)
         Q2  = Geom% Volume(c0+2)*(GSet% STotal(g,c0+2) + tau*Set% Psi(g,c0+2,Angle)) + &
               falpha*Geom% Area(c0+2)*PsiM1(g,2)

         b0  =  GSet% Sigt(g,zone)*Geom% Volume(c0+1) + angTerm1 
         a1  =  GSet% Sigt(g,zone)*Geom% Volume(c0+2) + angTerm2 
         den =  one/(a0*a0 + b0*a1)

         Set% Psi(g,c0+1,Angle) = ( a1*Q1 - a0*Q2 )*den
         Set% Psi(g,c0+2,Angle) = ( b0*Q2 + a0*Q1 )*den
         psi_inc(g,Angle)       = Set% Psi(g,c0+2,Angle)

!        Compute (m+1/2) and scalar intensities

         PsiM1(g,1)  = quadtau*Set% Psi(g,c0+1,Angle) +  &
                      (one - quadtau)*PsiM1(g,1)
         PsiM1(g,2)  = quadtau*Set% Psi(g,c0+2,Angle) +  &
                      (one - quadtau)*PsiM1(g,2)

         Set% Phi(g,c0+1) = Set% Phi(g,c0+1) + weight*Set% Psi(g,c0+1,Angle)
         Set% Phi(g,c0+2) = Set% Phi(g,c0+2) + weight*Set% Psi(g,c0+2,Angle)

       enddo

     enddo AngleLoop2

!  Finishing Direction

     Set% Psi(:,c0+1,NumAngles) = PsiM1(:,1)
     Set% Psi(:,c0+2,NumAngles) = PsiM1(:,2)

   enddo ZoneLoop2

   do Angle=mPos1,mPos2
     Set% PsiB(:,bOut,Angle) = psi_inc(:,Angle)
   enddo


   if ( outerBdyShared(RadBoundary) ) then
     do Angle=mPos1,mPos2
       CommT => getMessage(CSet, outerSharedID, Angle)

       CommT% psibsend(:,1) = Set% PsiB(:,bOut,Angle)

       call MPIStart(CommT% irequest(1))
       call MPIWait(CommT% irequest(1))
     enddo
   endif

   mstdr = mstdr + nAngLevel + 1

   enddo XiLevelLoop

   deallocate( PsiM )

!  Tally the number of intensity iterations performed

   CSet% fluxSweeps = CSet% fluxSweeps + NumAngles


   return
   end subroutine SweepSphereSCBsv 

