target_sources( teton PRIVATE
  AdjustFluxExchangeGoals.F90
  AdjustGreyGoals.F90
  AdjustNonlinearGoals.F90
  AdjustRadEnergyDensityGoal.F90
  AdjustTemperatureGoals.F90
  AdjustFluxExchangeRelTol.F90
  AdjustFluxExchangeMaxIts.F90
  AdjustGreyRelTol.F90
  AdjustGreyMaxIts.F90
  AdjustNonlinearRelTol.F90
  AdjustNonlinearMaxIts.F90
  AdjustRadEnergyDensityRelTol.F90
  AdjustTemperatureRelTol.F90
  AdjustTemperatureMaxIts.F90
  AppendSourceToPsi.F90
  ConstructBoundary.F90
  ConstructDtControls.F90
  ConstructEditor.F90
  ConstructGeometry.F90
  ConstructIterControls.F90
  ConstructMaterial.F90
  ConstructPhaseSpaceSets.F90
  ConstructQuadrature.F90
  ConstructRadIntensity.F90
  ConstructSize.F90
  DestructMeshData.F90
  GetDefaultGreyRelTol.F90
  GetDefaultGreyMaxIt.F90
  GetDefaultIncidentFluxRelTol.F90
  GetDefaultIncidentFluxMaxIt.F90
  GetDefaultInnerNlRelTol.F90
  GetDefaultInnerNlMaxIt.F90
  GetDefaultOuterIntensityRelTol.F90
  GetDefaultOuterMaxIt.F90
  GetDefaultOuterTempRelTol.F90
  GetFluxExchangeMaxIts.F90
  GetFluxExchangeRelTol.F90
  GetGreyMaxIts.F90
  GetGreyRelTol.F90
  GetNonlinearMaxIts.F90
  GetNonlinearRelTol.F90
  GetRadEnergyDensityRelTol.F90
  GetTemperatureMaxIts.F90
  GetTemperatureRelTol.F90
  InitTeton.F90
  ResetPsi.F90
  ScalePsir.F90
  SurfaceEdit.F90
  addBoundary.F90
  addProfile.F90
  checkInputSanity.F90
  checkSharedBoundary.F90
  dtnew.F90
  getAngleBins.F90
  getBeamMetric.F90
  getCornerVolumes.F90
  getDtMessage.F90
  getDtControlInfo.F90
  getEdits.F90
  getGroupAngleOffsets.F90
  getMaterialTemperature.F90 
  getNumAngleBins.F90
  getNumSNSets.F90
  getOpacity.F90
  getPsiPointer.F90
  getRadiationDeposited.F90
  getRadiationEnergyDensity.F90
  getRadiationFlux.F90
  getRadiationForce.F90
  getRadiationTemperature.F90
  getRunStats.F90
  getScalarIntensity.F90
  getVolume.F90
  getZonalPsi.F90
  initMaterial.F90
  initOpacity.F90
  normalizeMaterial.F90
  publishEdits.F90
  profint.F90
  queryPsi.F90
  radtr.F90
  resetProfile.F90
  resetTimers.F90
  rtedit.F90
  setCommunicationGroup.F90
  setGTAOptions.F90
  setMaterial.F90
  setMaterialSource.F90
  setNodePosition.F90
  setNodeVelocity.F90
  setOpacity.F90
  setOppositeFace.F90
  setRadiationFlux.F90
  setRadiationForce.F90
  setSharedFace.F90
  setTetonZone.F90
  setTetonZone1D.F90
  setTimeStep.F90
  ResetSize.F90
)
