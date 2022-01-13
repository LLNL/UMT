!***********************************************************************
!                       Last Update:  10/2016, PFN                     *
!                                                                      *
!   SweepPlanarSCB - Transport sweep for 1D planar geometry with the   *
!                    Simple Corner Balance spatial discretization.     *
!                                                                      *
!***********************************************************************

   subroutine SweepPlanarSCB(setID, Groups, NumAngles, savePsi)

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

   integer, intent(in)         :: setID
   integer, intent(in)         :: Groups
   integer, intent(in)         :: NumAngles
   logical(kind=1), intent(in) :: savePsi

!  Local

   type(SetData),          pointer  :: Set
   type(CommSet),          pointer  :: CSet
   type(AngleSet),         pointer  :: ASet
   type(GroupSet),         pointer  :: GSet
   type(Boundary),         pointer  :: BdyT
   type(Communicator),     pointer  :: CommT

   integer    :: c0
   integer    :: Angle,g,zone,nzones
   integer    :: mNeg1, mNeg2, mPos1, mPos2
   integer    :: angleRef
   integer    :: bIn, bOut, innerBdyID, outerBdyID
   integer    :: QuadID
   integer    :: innerSharedID, outerSharedID

   real(adqt) :: a0, a1, b0, den, Q1, Q2
   real(adqt) :: omega
   real(adqt) :: weight
   real(adqt) :: Psi1
   real(adqt) :: Psi2
   real(adqt) :: tau

   real(adqt) :: psi_inc(Groups, NumAngles)

!  Constants

   Set  => getSetData(Quad, setID)
   CSet => getCommSetFromSetID(Quad, setID)
   ASet => getAngleSetFromSetID(Quad, setID) 
   GSet => getGroupSetFromSetID(Quad, setID)  

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

   mNeg1 = 1
   mNeg2 = NumAngles/2 
   mPos1 = mNeg2 + 1
   mPos2 = NumAngles 

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

     AngleLoop1: do Angle=mNeg1,mNeg2

!  Angular coefficients

       omega  = ASet% omega(1,Angle)
       weight = ASet% weight(Angle)
       a0     = half*omega

       do g=1,Groups

         Q1  = Geom% Volume(c0+1)*(GSet% STotal(g,c0+1) + tau*Set% Psi(g,c0+1,Angle)) 
         Q2  = Geom% Volume(c0+2)*(GSet% STotal(g,c0+2) + tau*Set% Psi(g,c0+2,Angle)) - &
               omega*psi_inc(g,Angle)

         b0  =  GSet% Sigt(g,zone)*Geom% Volume(c0+1) - a0 
         a1  =  GSet% Sigt(g,zone)*Geom% Volume(c0+2) - a0 
         den =  one/(a0*a0 + b0*a1)

         Psi1 = ( a1*Q1 - a0*Q2 )*den
         Psi2 = ( a0*Q1 + b0*Q2 )*den
         psi_inc(g,Angle)     = Psi1

         if ( savePsi ) then
           Set% Psi(g,c0+1,Angle) = Psi1
           Set% Psi(g,c0+2,Angle) = Psi2
         endif

         Set% Phi(g,c0+1) = Set% Phi(g,c0+1) + weight*Psi1
         Set% Phi(g,c0+2) = Set% Phi(g,c0+2) + weight*Psi2

       enddo

     enddo AngleLoop1

   enddo ZoneLoop1

!  Inner Boundary

   do Angle=1,mNeg2
     Set% PsiB(:,bIn,Angle) = psi_inc(:,Angle)
   enddo

   if ( innerBdyShared(RadBoundary) ) then

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

   if ( innerBdyShared(RadBoundary) ) then

     do Angle=mPos1,mPos2
       CommT => getMessage(CSet, innerSharedID, Angle)

       call MPIWait(CommT% irequest(2))

       Set% PsiB(:,bIn,Angle) = CommT% psibrecv(:,1)
     enddo

   elseif ( innerBdyReflect(RadBoundary) ) then

     BdyT => getBoundary(RadBoundary, innerBdyID)

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

     AngleLoop2: do Angle=mPos1,mPos2

!  Angular coefficients

       omega  = ASet% omega(1,Angle)
       weight = ASet% weight(Angle)
       a0     = half*omega

       do g=1,Groups

         Q1  = Geom%Volume(c0+1)*(GSet% STotal(g,c0+1) + tau*Set% Psi(g,c0+1,Angle)) + &
               omega*psi_inc(g,Angle)
         Q2  = Geom%Volume(c0+2)*(GSet% STotal(g,c0+2) + tau*Set% Psi(g,c0+2,Angle)) 

         b0  =  GSet% Sigt(g,zone)*Geom%Volume(c0+1) + a0 
         a1  =  GSet% Sigt(g,zone)*Geom%Volume(c0+2) + a0 
         den =  one/(a0*a0 + b0*a1)

         Psi1 = ( a1*Q1 - a0*Q2 )*den
         Psi2 = ( b0*Q2 + a0*Q1 )*den
         psi_inc(g,Angle)     = Psi2

         if ( savePsi ) then
           Set% Psi(g,c0+1,Angle) = Psi1
           Set% Psi(g,c0+2,Angle) = Psi2
         endif

         Set% Phi(g,c0+1) = Set% Phi(g,c0+1) + weight*Psi1
         Set% Phi(g,c0+2) = Set% Phi(g,c0+2) + weight*Psi2

       enddo

     enddo AngleLoop2

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

!  Tally the number of intensity iterations performed

   CSet% fluxSweeps = CSet% fluxSweeps + NumAngles


   return
   end subroutine SweepPlanarSCB 

