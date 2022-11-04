!***********************************************************************
!                        Version 1:  08/94, PFN                        *
!                                                                      *
!   BoundaryEdit - Computes edits of escaping and incident power on    *
!                  external boundaries. Note that edits are computed   *
!                  in the "Lab" frame.                                 * 
!                                                                      *
!***********************************************************************
   subroutine BoundaryEdit(setID) 

   use flags_mod
   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use BoundaryList_mod
   use Boundary_mod
   use SetData_mod
   use AngleSet_mod
   use CommSet_mod

   implicit none

!  Arguments

   integer,  intent(in)         :: setID

!  Local Variables

   type(SetData),      pointer  :: Set
   type(AngleSet),     pointer  :: ASet
   type(CommSet),      pointer  :: CSet
   type(Boundary),     pointer  :: BdyT

   integer    :: angle
   integer    :: NumAngles
   integer    :: b
   integer    :: b0
   integer    :: nBdyElem
   integer    :: c
   integer    :: g
   integer    :: g0
   integer    :: ngr
   integer    :: Groups
   integer    :: n
   integer    :: nVacuum
   integer    :: nSource
   integer    :: polarAngle
   integer    :: COMM_GROUP

   real(adqt) :: angdota
   real(adqt) :: current
   real(adqt) :: geometryFactor
   real(adqt) :: factor
   real(adqt) :: weight
   real(adqt) :: VoCmag2
   real(adqt) :: D
   real(adqt) :: lambdaD
   real(adqt) :: lambdaD3

!  Constants

   Set  => getSetData(Quad, setID)
   CSet => getCommSetFromSetID(Quad, setID)
   ASet => getAngleSetFromSetID(Quad, setID)

   geometryFactor = getGeometryFactor(Size)
   nVacuum        = getNumberOfVacuum(RadBoundary)
   nSource        = getNumberOfSource(RadBoundary)

   COMM_GROUP     = CSet% COMM_GROUP
   NumAngles      = ASet% NumAngles
   ngr            = Size% ngr
   Groups         = Set% Groups
   g0             = Set% g0

!  We accumulate the escaping power using the quadrature set polar bins.
!  Mapping to other tally bins is performed by the host code as needed.
!  and then interpolate to populate the tallies.

!  Initialize partial current edits - Note that these must be full-size
!  arrays (not just the Set dimensions) so that the sums and all reduces
!  work in all cases

   Set% RadPowerEscape(:)         = zero
   Set% RadPowerIncident(:)       = zero
   Set% PolarSectorPowerEscape(:,:) = zero

!  Accumulate exiting and incoming partial currents for all
!  non-reflecting boundaries.  Partial currents are accumulated
!  by group and group/polar sector.

   VacuumLoop: do n=1,nVacuum

     BdyT     => getVacuum(RadBoundary, n)
     nBdyElem =  getNumberOfBdyElements(BdyT)
     b0       =  getFirstBdyElement(BdyT) - 1

!  Compute (unit normal) dot (omega)*area and incident/exiting currents 

     AngleLoop1: do angle=1,NumAngles

       weight     = ASet% weight(angle)
       polarAngle = ASet% polarAngle(angle)

       if (weight > zero) then

         BELoop1: do b=1,nBdyElem

           c = BdyT% BdyToC(b)

!  Compute the transformation to the Lab frame
!  Reference: G. Pomraning, "Radiation Hydsrodynamics", p. 274

           VoCmag2  = DOT_PRODUCT( Geom% VoC(:,c), Geom% VoC(:,c) )
           D        = one + DOT_PRODUCT( ASet%omega(:,angle), Geom% VoC(:,c) )
           lambdaD  = D/sqrt(one - VoCmag2) 
           lambdaD3 = lambdaD*lambdaD*lambdaD
 
           if (Size% igeom == geometry_rz) then
             factor = weight*geometryFactor*BdyT% Radius(b)
           elseif (Size% igeom == geometry_sphere) then
             factor = weight*geometryFactor*BdyT% Radius(b)*BdyT% Radius(b)
           else
             factor = weight*geometryFactor
           endif

           angdota = DOT_PRODUCT( ASet%omega(:,angle),BdyT%A_bdy(:,b) )

           if (angdota > zero) then

             do g=1,Set% Groups
               current = factor*lambdaD3*angdota*Set% Psi(g,c,angle)
               Set% RadPowerEscape(g)                     = Set% RadPowerEscape(g)                     + current
               Set% PolarSectorPowerEscape(g, polarAngle) = Set% PolarSectorPowerEscape(g, polarAngle) + current
             enddo
           endif

         enddo BELoop1

       endif

     enddo AngleLoop1

   enddo VacuumLoop

!  Source Boundaries

   SourceLoop: do n=1,nSource

     BdyT     => getSource(RadBoundary, n)
     nBdyElem =  getNumberOfBdyElements(BdyT)
     b0       =  getFirstBdyElement(BdyT) - 1

!  Compute (unit normal) dot (omega)*area and incident/exiting currents
                                                                                                    
     AngleLoop2: do angle=1,NumAngles

       weight     = ASet% weight(angle)
       polarAngle = ASet% polarAngle(angle)

       if (weight > zero) then

         BELoop2: do b=1,nBdyElem

           c = BdyT% BdyToC(b)

!  Compute the transformation to the Lab frame
!  Reference: G. Pomraning, "Radiation Hydrodynamics", p. 274

           VoCmag2  = DOT_PRODUCT( Geom% VoC(:,c), Geom% VoC(:,c) )
           D        = one + DOT_PRODUCT( ASet%omega(:,angle), Geom% VoC(:,c) )
           lambdaD  = D/sqrt(one - VoCmag2)
           lambdaD3 = lambdaD*lambdaD*lambdaD

           if (Size% igeom == geometry_rz) then
             factor = weight*geometryFactor*BdyT% Radius(b)
           elseif (Size% igeom == geometry_sphere) then
             factor = weight*geometryFactor*BdyT% Radius(b)*BdyT% Radius(b)
           else
             factor = weight*geometryFactor
           endif

           angdota = DOT_PRODUCT( ASet%omega(:,angle),BdyT%A_bdy(:,b) )
 
           if (angdota > zero) then

             do g=1,Set% Groups
               current = factor*lambdaD3*angdota*Set% Psi(g,c,angle)
               Set% RadPowerEscape(g)                     = Set% RadPowerEscape(g)                     + current
               Set% PolarSectorPowerEscape(g, polarAngle) = Set% PolarSectorPowerEscape(g, polarAngle) + current
             enddo

           else

             do g=1,Set% Groups
               current = factor*lambdaD3*angdota*Set% PsiB(g,b0+b,angle)
               Set% RadPowerIncident(g) = Set% RadPowerIncident(g) - current
             enddo

           endif
                                                                                                    
         enddo BELoop2

       endif
 
     enddo AngleLoop2
          
   enddo SourceLoop

   return
   end subroutine BoundaryEdit 


