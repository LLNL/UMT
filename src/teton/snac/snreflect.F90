!***********************************************************************
!                        Version 1:  09/96, PFN                        *
!                                                                      *
!   SNREFLECT - This routine, called by SNFLWXYZ, computes the         *
!               boundary flux (PSIB) for angle m, on reflecting        *
!               boundaries for which angle m is incoming.              *
!                                                                      *
!***********************************************************************

   subroutine snreflect(SnSweep, setID, Minc, PsiB) 

   use kind_mod
   use constant_mod
   use Size_mod
   use Quadrature_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use Boundary_mod
   use SetData_mod
   use AngleSet_mod

   implicit none

!  Arguments

   logical (kind=1),     intent(in)    :: SnSweep
   integer,              intent(in)    :: setID
   integer,              intent(in)    :: Minc
   real(adqt), optional, intent(inout) :: PsiB(Size% nSurfElem,Size% nangGTA)

!  Local Variables

   type(SetData),    pointer  :: Set
   type(AngleSet),   pointer  :: ASet
   type(Boundary),   pointer  :: BdyT

   integer    :: reflID 
   integer    :: angle0
   integer    :: b 
   integer    :: b0
   integer    :: nBdyElem
   integer    :: Mref
   integer    :: nReflecting

!  Constants

   nReflecting =  getNumberOfReflecting(RadBoundary)
   Set         => getSetData(Quad, setID)
   ASet        => getAngleSetFromSetID(Quad, setID)
   angle0      =  Set% angle0

!  Loop over reflecting-boundary sets:
 
   ReflectingLoop: do reflID=1,nReflecting

     BdyT      => getReflecting(RadBoundary, reflID)
     nBdyElem  =  getNumberOfBdyElements(BdyT)
     b0        =  getFirstBdyElement(BdyT) - 1

     Mref      =  getReflectedAngle(ASet, reflID, Minc)

     if (Mref > 0) then

       if ( SnSweep ) then

         do b=1,nBdyElem
           Set% PsiB(:,b0+b,Minc) = Set% PsiB(:,b0+b,Mref)
         enddo

       else

         do b=1,nBdyElem
           PsiB(b0+b,angle0+Minc) = PsiB(b0+b,angle0+Mref)
         enddo
       endif

     endif
 
   enddo ReflectingLoop 

 
   return
   end subroutine snreflect

