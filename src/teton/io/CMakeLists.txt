target_sources( teton PRIVATE
  ConduitCheckpoint_mod.F90
  )

if(NOT ENABLE_MINIAPP_BUILD AND ENABLE_SILO)

  target_sources( teton PRIVATE
    dropTetonVariables.F90
    loadTetonVariables.F90
  )

endif()
