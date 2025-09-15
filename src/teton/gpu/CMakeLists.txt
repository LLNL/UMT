target_sources( teton PRIVATE
  CornerSweepUCBrz_OMPOL.F90
  CornerSweepUCBxyz_OMPOL.F90
  InitSweepGreyUCBrz_OMPOL.F90
  InitSweepGreyUCBxyz_OMPOL.F90
  SetSweep_OMPOL.F90
  SetSweep_CUDA.F90     # CUDA-C sweep
  SweepUCBxyzToGPU.F90  # CUDA-C sweep
  SweepGreyUCBrz_OMPOL.F90
  SweepGreyUCBxyz_OMPOL.F90
  SweepUCBrz_OMPOL.F90
  SweepUCBxyz_OMPOL.F90
  OMPWrappers_mod.F90
  OMPUtilities_mod.F90
  UpdateScalarIntensity_OMPOL.F90
  cuda_utils_mod.F90
  finalizeGPUMemory.F90
  initializeGPUMemory.F90
)

if (ENABLE_CUDA)
  target_sources( tetoncuda PRIVATE
    GPU_SweepUCBxyz.cu  # CUDA-C sweep
  )

endif()
