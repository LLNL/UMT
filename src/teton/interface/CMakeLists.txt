if (ENABLE_BLUEPRINT_INTERFACE)
  target_sources( teton PRIVATE
    TetonBlueprint.cc
    TetonConduitInterface.cc
  )
endif()
