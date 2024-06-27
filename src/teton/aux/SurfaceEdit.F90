#include "macros.h"

!***********************************************************************
!                        Version 1:  08/94, PFN                        *
!                        Version 2:  04/20, BCY                        *
!                                                                      *
!   SurfaceEdit - Computes edits of escaping and incident              *
!        power on an arbitrary surface, binned by angle, group, and    *
!        time shift.                                                   *
!                                                                      *
!***********************************************************************
   subroutine SurfaceEdit(nCornerFaces, labFrame,    &
                          cornerList, zoneFaceList,  &
                          timeShift, centerPoint,    &
                          numAngleBins,              &
                          numGroups,                 &
                          numTimeBins,               &
                          timeBinBoundaries,         &
                          computeIncident,           &
                          scaleTally,                &
                          calcErrorMetrics,          &
                          tally,                     &
                          tallyIncident,             &
                          errEstShift,               &
                          errEstSrcSize) BIND(C,NAME="teton_surfaceedit_internal")

   use ISO_C_BINDING
   use flags_mod
   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use constant_mod
   use Size_mod
   use Geometry_mod
   use QuadratureList_mod
   use Editor_mod
   use SetData_mod
   use AngleSet_mod
   use radconstant_mod
   use TimeStepControls_mod

   implicit none

!  Arguments

   integer(C_INT),  intent(in)    :: nCornerFaces
   logical(C_BOOL), intent(in)    :: labFrame
   integer(C_INT),  intent(in)    :: cornerList(nCornerFaces)
   integer(C_INT),  intent(in)    :: zoneFaceList(nCornerFaces)
   logical(C_BOOL), intent(in)    :: timeShift
   real(C_DOUBLE),  intent(in)    :: centerPoint(Size% nDim)
   integer(C_INT),  intent(in)    :: numAngleBins
   integer(C_INT),  intent(in)    :: numGroups
   integer(C_INT),  intent(in)    :: numTimeBins
   real(C_DOUBLE),  intent(in)    :: timeBinBoundaries(numTimeBins+1)
   logical(C_BOOL), intent(in)    :: computeIncident
   real(C_DOUBLE),  intent(in)    :: scaleTally
   logical(C_BOOL), intent(in)    :: calcErrorMetrics
   real(C_DOUBLE),  intent(inout) :: tally(numTimeBins*numAngleBins*numGroups)
   real(C_DOUBLE),  intent(inout) :: tallyIncident(numTimeBins*numAngleBins &
                                                    *numGroups)
   real(C_DOUBLE),  intent(inout) :: errEstShift(numTimeBins*2)
   real(C_DOUBLE),  intent(inout) :: errEstSrcSize(numTimeBins*2)

!  Local Variables

   real(C_DOUBLE)                 :: tempTally(numTimeBins*numAngleBins*numGroups)
   real(C_DOUBLE), allocatable    :: tempTallyIncident(:)

   real(C_DOUBLE), allocatable    :: tempErrEstShift(:)
   real(C_DOUBLE), allocatable    :: tempErrEstSrcSize(:)
   real(C_DOUBLE)                 :: timeBinDistribution(numTimeBins)

   type(SetData),  pointer  :: Set
   type(AngleSet), pointer  :: ASet

   logical    :: calcErrorMetricsConfirmed ! check whether time shifting is on

   integer    :: angle
   integer    :: polarAngle
   integer    :: polarAngle0
   integer    :: timeBin
   integer    :: timeBinFinal
   integer    :: timeBin0
   integer    :: timeBinFinal0
   integer    :: timeBinIterator
   integer    :: timeBinPlusNBins
   integer    :: timeBinFinalPlusNBins
   integer    :: NumAngles
   integer    :: numGroupAngleBins ! numAngleBins*numGroups
   integer    :: g
   integer    :: g0
   integer    :: iCornerFace
   integer    :: setID
   integer    :: nSets
   integer    :: c
   integer    :: cOpp
   integer    :: zface
   integer    :: cface
   integer    :: zone
   integer    :: nCFaces

   real(adqt) :: angdota
   real(adqt) :: current
   real(adqt) :: geometryFactorTimesDt
   real(adqt) :: factor
   real(adqt) :: weight
   real(adqt) :: VoCmag2
   real(adqt) :: D
   real(adqt) :: lambdaD
   real(adqt) :: lambdaD3
   real(adqt) :: distParallel
   real(adqt) :: distPerp
   real(adqt) :: oldRadTime
   real(adqt) :: currentRadTime
   real(adqt) :: shiftedRadTimes(2)
   real(adqt) :: radTimeShift
   real(adqt) :: dtrad

   real(adqt), allocatable :: deltasFromCenter(:,:)
   real(adqt), allocatable :: sqDistsFromCenter(:)

! Local corner face indices that make up the surface:
   integer,    allocatable :: cFaceList(:)

! Quiet compiler 'variable may not be initialized' warnings on a few variables.
   timebin0 = 0
   factor   = zero

   ! Check some inputs:
   if (numGroups /= 1) then
     TETON_ASSERT(Size% ngr == numGroups, "numGroups in teton_surfaceedit must be either 1 or the # of Teton groups")
   endif
   TETON_ASSERT(numTimeBins > 0, "Number of time bins must be positive.")
   TETON_ASSERT(minval(cornerList) > 0, "corner indices must be greater than 0")

!  Constants

   numGroupAngleBins = numGroups*numAngleBins
   nSets = getNumberOfSets(Quad)
   dtrad = getRadTimeStep(DtControls)
   geometryFactorTimesDt = getGeometryFactor(Size)*dtrad

   if (nCornerFaces > 0) then
     allocate( cFaceList(nCornerFaces) )
     cFaceList(:) = -1
   endif

!  Get the starting time

   currentRadTime = getRadTime(DtControls)
   oldRadTime = currentRadTime - dtrad

!  For the input corner list and opposite corners, find the correct
!  corner-face index

   if (timeShift .and. nCornerFaces > 0) then
     allocate(deltasFromCenter(Size% nDim, nCornerFaces))
     allocate(sqDistsFromCenter(nCornerFaces))
   endif

   do iCornerFace=1,nCornerFaces
     c     = cornerList(iCornerFace)
     zface = zoneFaceList(iCornerFace)
     zone  = Geom% cToZone(c)

     ! This spatial vector = the corner face center minus the centerPoint

     if (timeShift) then
       deltasFromCenter(:, iCornerFace) = Geom% px(:,c) - centerPoint
       sqDistsFromCenter(iCornerFace) = DOT_PRODUCT(deltasFromCenter(:,iCornerFace),deltasFromCenter(:,iCornerFace))
     endif

!    Find the corner-face index associated with this (local) corner, c,
!    and this (local) zone face

     nCFaces = Geom% nCFacesArray(c)

     CFaceLoop : do cface=1,nCFaces
       if (Geom% CToFace(cface,c) == zface) then
         cFaceList(iCornerFace) = cface
         exit CFaceLoop
       endif
     enddo CFaceLoop

     TETON_ASSERT(cFaceList(iCornerFace) > 0, "Could not find corner face from corner "//char(c)//" and zface"//char(zface))
   enddo

!  Initialize temporary arrays:

   tempTally(:) = zero

   if (computeIncident) then
     allocate( tempTallyIncident(numTimeBins*numAngleBins*numGroups) )
     tempTallyIncident(:) = zero
   endif

   calcErrorMetricsConfirmed = timeShift .and. calcErrorMetrics
   if (calcErrorMetricsConfirmed) then
     if (computeIncident) then
       allocate( tempErrEstShift(numTimeBins*2)   )
       allocate( tempErrEstSrcSize(numTimeBins*2) )
     else
       allocate( tempErrEstShift(numTimeBins)   )
       allocate( tempErrEstSrcSize(numTimeBins) )
     endif

     tempErrEstShift(:)   = zero
     tempErrEstSrcSize(:) = zero
   endif

   SetLoop: do setID=1,nSets

     Set    => getSetData(Quad, setID)
     ASet   => getAngleSetFromSetID(Quad, setID)

     NumAngles = ASet% NumAngles
     if (numGroups /= 1) then
       g0 = Set% g0
     else
       g0 = 0
     endif

     if (numAngleBins /= 1) then
       TETON_ASSERT(ASet% nPolarAngles == numAngleBins, "numAngleBins must either be 1 or # of Teton polar angles")
     else
       polarAngle = 1
     endif

     ! Accumulate exiting and incoming partial currents this surface

     ! Compute (unit normal) dot (omega)*area and incident/exiting currents

     AngleLoop: do angle=1,NumAngles

       weight = ASet% weight(angle)
       if (numAngleBins > 1) then
         polarAngle = ASet% polarAngle(angle)
         polarAngle0 = g0 + (polarAngle - 1)*numGroups
       else
         polarAngle0 = g0
       endif

       if (weight > zero) then

         CFLoop: do iCornerFace=1,nCornerFaces
           c      = cornerList(iCornerFace) ! corner for this corner-face
           cface  = cFaceList(iCornerFace)
           cOpp   = Geom% cFP(cface,c)

!          Compute |\Omega \cdot n| * (area of corner face),
!          where n is normal to the corner face

           angdota = DOT_PRODUCT( ASet%omega(:,angle), Geom% A_fp(:,cface,c))

           if (.not. computeIncident .and. angdota <= zero) then
             cycle
           endif

           ! Compute time shift
           if (timeShift) then
             ! sign(1., angdota)*ASet%omega(:,angle) flips the angle if it's
             !    incoming, does nothing if it's outgoing.
             distParallel = DOT_PRODUCT( deltasFromCenter(:,iCornerFace), &
                                        sign(one,angdota)*ASet%omega(:,angle)    )
             radTimeShift = distParallel / speed_light
             ! Compute error metric quantities:
             distPerp = sqrt(sqDistsFromCenter(iCornerFace) - &
                             distParallel*distParallel )
           else
             distParallel = 0.0
             distPerp = 0.0
             radTimeShift = zero
           endif

           timeBinDistribution(:) = zero
           ! Time step spans (oldRadTime,currentRadTime)
           ! Shifted time step spans (oldRadTime-radTimeShift,
           !                          currentRadTime-radTimeShift)

           ! Find the first bin corresponding to the shifted time step:
           shiftedRadTimes(1) = oldRadTime-radTimeShift
           TimeBinLoop: do timeBin=1,numTimeBins-1
             if (shiftedRadTimes(1) < timeBinBoundaries(timeBin+1)) then
               exit TimeBinLoop
             endif
           enddo TimeBinLoop

           ! Find the last bin corresponding to the shifted time step:
           shiftedRadTimes(2) = currentRadTime-radTimeShift
           TimeBinFinalLoop: do timeBinFinal=timeBin,numTimeBins-1
             if (shiftedRadTimes(2) < timeBinBoundaries(timeBinFinal+1)) then
               exit TimeBinFinalLoop
             endif
           enddo TimeBinFinalLoop

           ! Now compute the fraction of the tally to be distributed in each
           !   time bin:
           if (timeBin == timeBinFinal) then
             ! If the step falls within one bin:
             timeBinDistribution(timeBin) = one
           else
             ! Else, distribute over the time bins:
             timeBinDistribution(timeBin) = (timeBinBoundaries(timeBin+1)-shiftedRadTimes(1))/dtrad
             do timeBinIterator=timeBin+1,timeBinFinal-1
               timeBinDistribution(timeBinIterator) = &
                 (timeBinBoundaries(timeBinIterator+1)-timeBinBoundaries(timeBinIterator))/dtrad
             enddo
             timeBinDistribution(timeBinFinal) = (shiftedRadTimes(2)-timeBinBoundaries(timeBinFinal))/dtrad
           endif
           TETON_ASSERT(abs(sum(timeBinDistribution) - one) < 1.e-12_adqt, "timeBinDistribution must sum to one.")

           ! offsets for timeBin and timeBinFinal:
           timeBin0 = (timeBin-1)*numGroupAngleBins+polarAngle0
           timeBinFinal0 = (timeBinFinal-1)*numGroupAngleBins+polarAngle0

           ! Compute the transformation to the Lab frame
           ! Reference: G. Pomraning, "Radiation Hydrodynamics", p. 274

           if ( labFrame ) then
             VoCmag2  = DOT_PRODUCT( Geom% VoC(:,c), Geom% VoC(:,c) )
             D        = one + DOT_PRODUCT( ASet%omega(:,angle), Geom% VoC(:,c) )
             lambdaD  = D/sqrt(one - VoCmag2)
             lambdaD3 = lambdaD*lambdaD*lambdaD
           else
             lambdaD3 = one
           endif

           if (Size% igeom == geometry_rz) then
             factor = weight*geometryFactorTimesDt*lambdaD3*Geom% RadiusFP(cface,c)
           elseif (Size% igeom == geometry_xyz) then
             factor = weight*geometryFactorTimesDt*lambdaD3
           elseif (Size% igeom == geometry_sphere) then
             factor = weight*geometryFactorTimesDt*lambdaD3*Geom% Radius(c)*Geom% Radius(c)
           else
             TETON_ASSERT(.false., "Unknown geometry type in Teton's SurfaceEdit.F90")
           endif

           if (angdota > zero) then

             do g=1,Set% Groups
               current = factor*angdota*Set% Psi(g,c,angle)
               if (numGroups > 1) then
                 tempTally(timeBin0+g:timeBinFinal0+g:numGroupAngleBins) = &
                   tempTally(timeBin0+g:timeBinFinal0+g:numGroupAngleBins) &
                     + current*timeBinDistribution(timeBin:timeBinFinal)
               else
                 tempTally(timeBin0+1:timeBinFinal0+1:numGroupAngleBins) = &
                   tempTally(timeBin0+1:timeBinFinal0+1:numGroupAngleBins) &
                     + current*timeBinDistribution(timeBin:timeBinFinal)
               endif

               if (calcErrorMetricsConfirmed) then
                 tempErrEstShift(timeBin:timeBinFinal) &
                   = tempErrEstShift(timeBin:timeBinFinal) &
                     + current*abs(distParallel)*timeBinDistribution(timeBin:timeBinFinal)
                 tempErrEstSrcSize(timeBin:timeBinFinal) &
                   = tempErrEstSrcSize(timeBin:timeBinFinal) &
                     + current*distPerp*timeBinDistribution(timeBin:timeBinFinal)
               endif
             enddo

           else if (angdota < zero) then

             ! Note that we only ever reach this part of the code if
             !  computeIncident = true
             TETON_ASSERT(computeIncident, "Should not try to compute incident power in SurfaceEdit if computeIncident is .false.")
             TETON_ASSERT(cOpp > 0, "cOpp must be a positive index")

             timeBinPlusNBins = timeBin + numTimeBins
             timeBinFinalPlusNBins = timeBinFinal + numTimeBins

             do g=1,Set% Groups
               ! Check if it is a boundary corner face:
               if (cOpp > Size% ncornr) then ! This is a boundary corner face
                 current = factor*angdota*Set% PsiB(g,cOpp-Size% ncornr,angle)
               else
                 current = factor*angdota*Set% Psi(g,cOpp,angle)
               endif
               ! currents should be negative if psi are positive!

               if (numGroups > 1) then
                 tempTallyIncident(timeBin0+g:timeBinFinal0+g:numGroupAngleBins) = &
                   tempTallyIncident(timeBin0+g:timeBinFinal0+g:numGroupAngleBins) &
                     - current*timeBinDistribution(timeBin:timeBinFinal)
               else
                 tempTallyIncident(timeBin0+1:timeBinFinal0+1:numGroupAngleBins) = &
                   tempTallyIncident(timeBin0+1:timeBinFinal0+1:numGroupAngleBins) &
                     - current*timeBinDistribution(timeBin:timeBinFinal)
               endif

               if (calcErrorMetricsConfirmed) then
                 tempErrEstShift(timeBinPlusNBins:timeBinFinalPlusNBins) &
                   = tempErrEstShift(timeBinPlusNBins:timeBinFinalPlusNBins) &
                     - current*abs(distParallel)*timeBinDistribution(timeBin:timeBinFinal)
                 tempErrEstSrcSize(timeBinPlusNBins:timeBinFinalPlusNBins) &
                   = tempErrEstSrcSize(timeBinPlusNBins:timeBinFinalPlusNBins) &
                     - current*distPerp*timeBinDistribution(timeBin:timeBinFinal)
               endif

             enddo

           endif ! This if checks whether angle is escape or incident

         enddo CFLoop

       endif ! If weight of this angle bin > 0

     enddo AngleLoop

   enddo SetLoop

   call MPIAllReduce(tempTally, "sum", MY_COMM_GROUP)
   tempTally(:) = tempTally(:)*scaleTally
   tally(:) = tally(:) + tempTally(:)

   if (computeIncident) then

     call MPIAllReduce(tempTallyIncident, "sum", MY_COMM_GROUP)
     tempTallyIncident(:) = tempTallyIncident(:)*scaleTally
     tallyIncident(:) = tallyIncident(:) + tempTallyIncident(:)
     deallocate(tempTallyIncident)

   endif

   if (nCornerFaces > 0) then
     deallocate(cFaceList)
   endif

   if (calcErrorMetricsConfirmed) then

     call MPIAllReduce(tempErrEstShift  , "sum", MY_COMM_GROUP)
     call MPIAllReduce(tempErrEstSrcSize, "sum", MY_COMM_GROUP)

     tempErrEstShift(:)   = tempErrEstShift(:)   * scaleTally
     tempErrEstSrcSize(:) = tempErrEstSrcSize(:) * scaleTally

     if (computeIncident) then

       errEstShift(:)   = errEstShift(:)   + tempErrEstShift(:)
       errEstSrcSize(:) = errEstSrcSize(:) + tempErrEstSrcSize(:)

     else

       errEstShift(1:numTimeBins)   = errEstShift(1:numTimeBins)   + tempErrEstShift(:)
       errEstSrcSize(1:numTimeBins) = errEstSrcSize(1:numTimeBins) + tempErrEstSrcSize(:)

     endif

     deallocate( tempErrEstShift   )
     deallocate( tempErrEstSrcSize )
   endif

   if (timeShift .and. nCornerFaces > 0) then
     deallocate( deltasFromCenter  )
     deallocate( sqDistsFromCenter )
   endif

   return
   end subroutine SurfaceEdit


