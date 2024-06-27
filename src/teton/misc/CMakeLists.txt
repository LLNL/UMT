target_sources( teton PRIVATE
  CodeChecks_mod.F90
  f90errors.F90
  mpi_param_mod.F90
  mpif90_mod.F90
  PrintEnergies.F90
)

# These provide a small replacement for PhysicsUtils and computing the Planck Function.
if (ENABLE_MINIAPP_BUILD)
   target_sources( teton PRIVATE
      EmissionUtils.f90
      NormalizedBlackBody.cc
      )
   target_sources( teton INTERFACE
      NormalizedBlackBody.hh )
endif()

