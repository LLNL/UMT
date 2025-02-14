!***********************************************************************
!                        Last Update:  01/2012, PFN                    *
!                                                                      *
!   RTEDIT  - Computes various radiation edits at the conclusion       *
!             of the radiation cycle as an edit.                       *
!                                                                      *
!***********************************************************************
   subroutine rtedit() BIND(C,NAME="teton_rtedit_new")

   USE ISO_C_BINDING
   use kind_mod
   use mpi_param_mod
   use mpif90_mod
   use constant_mod
   use radconstant_mod
   use Size_mod
   use Geometry_mod
   use Material_mod
   use RadIntensity_mod
   use QuadratureList_mod
   use TimeStepControls_mod
   use iter_control_list_mod
   use iter_control_mod
   use Editor_mod
   use SetData_mod

   implicit none

!  Local

   type(SetData),      pointer :: Set
   type(IterControl),  pointer :: nonLinearControl => NULL()
   
   integer, dimension (1) :: zoneTrMax 
   integer, dimension (1) :: zoneTeMax 

   integer    :: c, c0, nCorner, zone, nzones
   integer    :: g, gGlobal
   integer    :: bin, bin0
   integer    :: setID
   integer    :: nSets
   integer    :: processTr, processTe
   integer    :: maxNonLinearIters
   integer    :: zoneNonLinearIters
   integer    :: zoneNLItersMax
   integer    :: my_NLItersMax
   integer    :: processNLIters
   integer    :: numSpectrumAngleBins
   integer    :: ngr

   ! Use doubles to avoid overflow errors on very large meshes.
   real(adqt) :: numNonLinearIters
   real(adqt) :: totalNonLinearIters
   real(adqt) :: nZonesGlobal

   real(adqt) :: dtrad
   real(adqt) :: geometryFactor
   real(adqt) :: tr4min, mass, tfloor

   real(adqt) :: my_TrMax, my_TeMax, TrMax, TeMax, PhiAve
   real(adqt) :: deltaEMat, deltaERad
   real(adqt) :: EnergyRadiation
   real(adqt) :: PowerIncident
   real(adqt) :: PowerEscape
   real(adqt) :: PowerAbsorbed
   real(adqt) :: PowerEmitted
   real(adqt) :: PowerExtSources
   real(adqt) :: PowerCompton
   real(adqt) :: EnergyCheck
   real(adqt) :: ERad

!  Iteration Controls

   nonLinearControl => getIterationControl(IterControls, "nonLinear")

!  Parameters 

   numSpectrumAngleBins = getNumberOfSpectrumAngleBins(RadEdit)
   nSets                = getNumberOfSets(Quad)
   dtrad                = getRadTimeStep(DtControls)
   geometryFactor       = getGeometryFactor(Size)
   ngr                  = Size% ngr
   nzones               = Size% nzones
   tfloor               = Size% tfloor

!  Compute zone-average radiation and electron temperatures,
!  end-of-cycle radiation energy and energy change due to
!  external radiation sources (desrcRad)
!  "Phi" has units of energy/area/time 

   EnergyRadiation     = zero
   deltaEMat           = zero
   PowerAbsorbed       = zero
   PowerEmitted        = zero
   PowerCompton        = zero
   PowerExtSources     = zero
   TrMax               = zero
   TeMax               = zero
   tr4min              = Size% tr4floor 
   maxNonLinearIters   = 0
   totalNonLinearIters = 0

!  Initialize global radiation boundary edits

   RadEdit% RadPowerEscape(:)         = zero
   RadEdit% RadPowerIncident(:)       = zero
   RadEdit% PolarSectorPowerEscape(:) = zero

   SetLoop: do setID=1,nSets
     Set => getSetData(Quad, setID)

!  Boundary Edits

     do g=1,Set% Groups
       gGlobal = Set% g0 + g
       RadEdit% RadPowerEscape(gGlobal)   = RadEdit% RadPowerEscape(gGlobal)   + &
                                      Set% RadPowerEscape(g)
       RadEdit% RadPowerIncident(gGlobal) = RadEdit% RadPowerIncident(gGlobal) + &
                                      Set% RadPowerIncident(g)
     enddo

     do bin=1,numSpectrumAngleBins
       bin0 = (bin - 1)*ngr
       do g=1,Set% Groups
         gGlobal = Set% g0 + g
         RadEdit% PolarSectorPowerEscape(bin0+gGlobal) =  &
                  RadEdit% PolarSectorPowerEscape(bin0+gGlobal) +  &
                  Set% PolarSectorPowerEscape(g,bin)
       enddo
     enddo

   enddo SetLoop

   call MPIAllReduce(RadEdit% RadPowerEscape,         "sum", MY_COMM_GROUP)
   call MPIAllReduce(RadEdit% RadPowerIncident,       "sum", MY_COMM_GROUP)
   call MPIAllReduce(RadEdit% PolarSectorPowerEscape, "sum", MY_COMM_GROUP)

   ZoneLoop: do zone=1,nzones

     nCorner = Geom% numCorner(zone)
     c0      = Geom% cOffSet(zone)
     ERad    = zero

     do c=1,nCorner
       do g=1,ngr
         Erad          = Erad          + Geom% Volume(c0+c)*Rad% PhiTotal(g,c0+c)
         PowerAbsorbed = PowerAbsorbed + Geom% Volume(c0+c)*Rad% PhiTotal(g,c0+c)* &
                                         Mat%Siga(g,zone)
       enddo

!  These are computed in UpdateMaterialCoupling
       PowerEmitted     = PowerEmitted     + Mat% PowerEmitted(c0+c)
       PowerCompton     = PowerCompton     + Mat% PowerCompton(c0+c)

     enddo

     EnergyRadiation  = EnergyRadiation + ERad 
     PhiAve           = ERad/(Geom% VolumeZone(zone)*rad_constant*speed_light)

!  Radiation temperature

     Mat%trz(zone) = sqrt( sqrt( max(PhiAve, tr4min) ) )

!  Maximum temperatures

     if (Mat%trz(zone) > TrMax) then
       TrMax     = Mat%trz(zone)
       zoneTrMax = zone
     endif

     if (Mat%tez(zone) > TeMax) then
       TeMax     = Mat%tez(zone)
       zoneTeMax = zone
     endif

!  Non-linear Iterations

     zoneNonLinearIters = 0
     do c=1,nCorner
       totalNonLinearIters = totalNonLinearIters + Mat%nonLinearIterations(c0+c)
       zoneNonLinearIters  = zoneNonLinearIters  + Mat%nonLinearIterations(c0+c)
     enddo

     if (zoneNonLinearIters > maxNonLinearIters) then
       maxNonLinearIters = zoneNonLinearIters 
       zoneNLItersMax    = zone
     endif

!  Energy changes in the matter

     mass            = geometryFactor*Mat%rho(zone)*Geom% VolumeZone(zone)
     deltaEMat       = deltaEMat       + mass*Mat%denez(zone)
     PowerExtSources = PowerExtSources + mass*Mat%SMatEff(zone)

   enddo ZoneLoop

   EnergyRadiation = geometryFactor*EnergyRadiation/speed_light
   PowerAbsorbed   = geometryFactor*PowerAbsorbed
   PowerEmitted    = geometryFactor*PowerEmitted
   PowerCompton    = geometryFactor*PowerCompton

!  Find total radiation energy and the radiation converted
!  to matter energy

   deltaERad = EnergyRadiation - RadEdit% EnergyRadBOC 

!  Find incident/escaping power on problem boundaries. The MPIAllReduce
!  of the group-dependent variables was performed in BoundaryEdit

   PowerEscape   = sum( RadEdit% RadPowerEscape(:) )
   PowerIncident = sum( RadEdit% RadPowerIncident(:) )

!  Energy Balance

   call MPIAllReduce(deltaEMat,       "sum", MY_COMM_GROUP)
   call MPIAllReduce(deltaERad,       "sum", MY_COMM_GROUP)
   call MPIAllReduce(EnergyRadiation, "sum", MY_COMM_GROUP)
   call MPIAllReduce(PowerAbsorbed,   "sum", MY_COMM_GROUP)
   call MPIAllReduce(PowerEmitted,    "sum", MY_COMM_GROUP)
   call MPIAllReduce(PowerExtSources, "sum", MY_COMM_GROUP)
   call MPIAllReduce(PowerCompton,    "sum", MY_COMM_GROUP)

   EnergyCheck = dtrad*(PowerIncident - PowerEscape + PowerExtSources) -  &
                 deltaERad - deltaEMat

!  Maximum temperatures and their locations

   my_TrMax = TrMax
   my_TeMax = TeMax

   call MPIAllReduce(TRMAX,    "max", MY_COMM_GROUP)
   call MPIAllReduce(TEMAX,    "max", MY_COMM_GROUP)

   if (my_TrMax == TrMax) then
     processTr = Size% myRankInGroup 
   else
     processTr = -1
     zoneTrMax = -1
   endif

   if (my_TeMax == TeMax) then
     processTe = Size% myRankInGroup 
   else
     processTe = -1
     zoneTeMax = -1
   endif
                                                                                       
   call MPIAllReduce(processTr, "max", MY_COMM_GROUP)
   call MPIAllReduce(processTe, "max", MY_COMM_GROUP)
   call MPIAllReduce(zoneTrMax, "max", MY_COMM_GROUP)
   call MPIAllReduce(zoneTeMax, "max", MY_COMM_GROUP)

!  Maximum number of non-linear iterations and their locations

   my_NLItersMax = maxNonLinearIters
   nZonesGlobal  = nzones

   ! Get maximum number of non linear iterations across all mesh domains.
   call MPIAllReduce(maxNonLinearIters, "max", MY_COMM_GROUP)
   ! Get total number of zones across all mesh domains
   call MPIAllReduce(nZonesGlobal, "sum", MY_COMM_GROUP)
   ! Get partial average of the non linear iterations ( local domain only )
   numNonLinearIters = totalNonLinearIters/nZonesGlobal
   ! Add partial averages together to get average.
   call MPIAllReduce(numNonLinearIters, "sum", MY_COMM_GROUP)

   if (my_NLItersMax == maxNonLinearIters) then
     processNLIters = Size% myRankInGroup
   else
     processNLIters = -1
   endif 

   call MPIAllReduce(processNLIters, "max", MY_COMM_GROUP)

   call setGlobalMaxIterationsTaken(nonLinearControl, maxNonLinearIters)
   call setNumberOfIterations(nonLinearControl, int(numNonLinearIters))

   call setProcessOfMax(nonLinearControl, processNLIters)  ! of processes in this group

!  Update Edit module

   call setEdits(RadEdit,                          &
                 TrMaxZone=zoneTrMax,              & 
                 TeMaxZone=zoneTeMax,              &
                 TrMaxProcess=processTr,           &
                 TeMaxProcess=processTe,           &
                 TrMax=TrMax,                      &
                 TeMax=TeMax,                      &
                 deltaEMat=deltaEMat,              &
                 EnergyRadiation=EnergyRadiation,  &
                 PowerIncident=PowerIncident,      &
                 PowerEscape=PowerEscape,          &
                 PowerAbsorbed=PowerAbsorbed,      &
                 PowerEmitted=PowerEmitted,        &
                 PowerExtSources=PowerExtSources,  &
                 PowerCompton=PowerCompton,        &
                 EnergyCheck=EnergyCheck)

!  Timing

   call MPIAllReduce(Size% MatCoupTimeCycle, "max", MY_COMM_GROUP)
   call MPIAllReduce(Size% SweepTimeCycle, "max", MY_COMM_GROUP)
   call MPIAllReduce(Size% GPUSweepTimeCycle, "max", MY_COMM_GROUP)
   call MPIAllReduce(Size% GTATimeCycle, "max", MY_COMM_GROUP)
   call MPIAllReduce(Size% RadtrTimeCycle, "max", MY_COMM_GROUP)
   call MPIAllReduce(Size% InitTimeCycle, "max", MY_COMM_GROUP)
   call MPIAllReduce(Size% FinalTimeCycle, "max", MY_COMM_GROUP)

   Size% MatCoupTimeTotal  = Size% MatCoupTimeTotal  + Size% MatCoupTimeCycle 
   Size% SweepTimeTotal    = Size% SweepTimeTotal    + Size% SweepTimeCycle
   Size% GPUSweepTimeTotal = Size% GPUSweepTimeTotal + Size% GPUSweepTimeCycle
   Size% GTATimeTotal      = Size% GTATimeTotal      + Size% GTATimeCYCLE
   Size% RadtrTimeTotal    = Size% RadtrTimeTotal    + Size% RadtrTimeCycle
   Size% InitTimeTotal     = Size% InitTimeTotal     + Size% InitTimeCycle
   Size% FinalTimeTotal    = Size% FinalTimeTotal    + Size% FinalTimeCycle


   return
   end subroutine rtedit

