if (ENABLE_BLUEPRINT_INTERFACE)
  target_sources( teton PRIVATE
    TetonSources.cc
    TetonBlueprint.cc
    TetonConduitInterface.cc
    TetonSurfaceTallies.cc
  )
endif()
