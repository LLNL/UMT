include(../cmake/SubDirList.cmake)

# Quiet the warning that the relative source path for INTERFACE
# will be converted to absolute path.
cmake_policy(SET CMP0076 NEW)

############################
# LIBRARY
############################

# add teton library target
add_library ( teton STATIC )

if( ENABLE_CUDA )
  add_library( tetoncuda OBJECT )
  target_link_libraries( teton tetoncuda ${CMAKE_EXTRA_CUDA_LIBS})
endif()

target_sources( teton INTERFACE include/TetonInterface.hh )

if (ENABLE_UMPIRE)
   target_include_directories( teton PRIVATE include
                               ${UMPIRE_INCLUDE_DIR}
                             )
endif()

if (ENABLE_CALIPER)
   target_include_directories( teton PRIVATE include
                               ${CALIPER_INCLUDE_DIR}
                               ${CALIPER_INCLUDE_DIR}/caliper/fortran
                             )
endif()

target_sources( teton INTERFACE include/TetonConduitInterface.hh )
target_include_directories( teton PRIVATE include
                               ${CONDUIT_INCLUDE_DIR}
                               ${CONDUIT_FORTRAN_MODULES_DIR})

if (NOT ENABLE_MINIAPP_BUILD)
   target_include_directories( teton PRIVATE include
                               ${PHYSICSUTILS_INCLUDE_DIR}
                               ${SILO_INCLUDE_DIR}
                             )
endif()

configure_file (include/TetonVersion.hh.in include/TetonVersion.hh )

target_sources( teton INTERFACE ${PROJECT_BINARY_DIR}/teton/include/TetonVersion.hh )
target_include_directories( teton PUBLIC include ${PROJECT_BINARY_DIR}/teton/include )

#################################################################
# WORKAROUND FOR IBM XLF COMPILER CRUSHING CPU THREAD PERFORMANCE
# Apply the -qoffload flag to only those source files with
# target pragmas.  Applying this flag to all source files reduces
# CPU thread performance by 10%-20%
#################################################################
if (NOT "${TETON_OpenMP_Offload_Fortran_FLAGS}" STREQUAL "" )
  set_source_files_properties(SOURCE
                                     control/initializeSets.F90
                                     control/initPhiTotal_OMPOL.F90
                                     control/initializeRadiationField_OMPOL.F90
                                     control/finalizeSets.F90
                                     control/getPhiTotal_OMPOL.F90
                                     gpu/initializeGPUMemory.F90
                                     gpu/finalizeGPUMemory.F90
                                     gpu/InitSweepGreyUCBrz_OMPOL.F90
                                     gpu/InitSweepGreyUCBxyz_OMPOL.F90
                                     gpu/GTASweep_OMPOL.F90
                                     gpu/OMPWrappers_mod.F90
                                     gpu/SetSweep_OMPOL.F90
                                     gpu/SweepGreyUCBrz_OMPOL.F90
                                     gpu/SweepGreyUCBxyz_OMPOL.F90
                                     gpu/SweepUCBrz_OMPOL.F90
                                     gpu/SweepUCBxyz_OMPOL.F90
                                     gpu/UpdateScalarIntensity_OMPOL.F90
                                     rt/initcomm.F90
                                     rt/RecvFlux.F90
                                     rt/SendFlux.F90

                              PROPERTY COMPILE_FLAGS ${TETON_OpenMP_Offload_Fortran_FLAGS})

  if (NOT ENABLE_MINIAPP_BUILD)
    set_source_files_properties(SOURCE
                                     nlsolver/UpdateMaterialCoupling_OMPOL.F90
                              PROPERTY COMPILE_FLAGS ${TETON_OpenMP_Offload_Fortran_FLAGS})
  endif()

endif()

if (NOT "${TETON_OpenMP_Offload_CXX_FLAGS}" STREQUAL "" )
  set_source_files_properties(SOURCE driver/test_driver.cc
                              PROPERTY COMPILE_FLAGS ${TETON_OpenMP_Offload_Fortran_FLAGS})
endif()

install( FILES "include/TetonInterface.hh" DESTINATION include )
install( FILES "include/TetonConduitInterface.hh" DESTINATION include )
install( FILES "${PROJECT_BINARY_DIR}/teton/include/TetonVersion.hh" DESTINATION include )

install( TARGETS teton ARCHIVE DESTINATION lib)

# Process any sub-dirs with CMakeList.txt files.
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_LIST_DIR})

foreach(subdir ${SUBDIRS})
  message( STATUS "Adding files from ${subdir}")
  add_subdirectory(${subdir})
endforeach()


############################
# TEST DRIVER EXECUTABLE
############################

if(ENABLE_TESTS)

  add_executable( test_driver
                  blueprint/TetonBlueprint.cc
                  driver/test_driver.cc )

  if (DEFINED TETON_OpenMP_Offload_CXX_FLAGS)
    set_source_files_properties(SOURCE
                                driver/test_driver.cc
                                PROPERTY COMPILE_FLAGS ${TETON_OpenMP_Offload_CXX_FLAGS})
  endif()

  target_include_directories( test_driver PRIVATE
                              ../src/include
                              ${CONDUIT_INCLUDE_DIR}
                            )
  target_link_libraries( test_driver
                         ${Fortran_COMPILER_LIB_DIR_LINK}
                         teton
                         ${CONDUIT_LIBRARIES}
                         ${CONDUITRELAY_LIBRARIES}
                         ${CONDUITBLUEPRINT_LIBRARIES}
                         ${SILO_LIBRARIES}
                         ${HDF5_LIBRARIES}
                         ${Z_LIBRARIES}
                         ${SZIP_LIBRARIES}
                       )

  if (NOT ENABLE_MINIAPP_BUILD)
    target_link_libraries( test_driver ${PHYSICSUTILS_LIBRARIES} )
  endif()

  if (ENABLE_MEMUSAGES)
    target_include_directories( test_driver PRIVATE ${MEMUSAGES_INCLUDE_DIR} )
    target_link_libraries( test_driver ${MEMUSAGES_LIBRARIES} )
  endif()

  if (ENABLE_UMPIRE)
    target_include_directories( test_driver PRIVATE ${UMPIRE_INCLUDE_DIR} )
    target_link_libraries( test_driver ${UMPIRE_LIBRARIES} )
  endif()

  if (ENABLE_CALIPER)
    target_include_directories( test_driver PRIVATE ${CALIPER_INCLUDE_DIR} ${ADIAK_INCLUDE_DIR})
    target_link_libraries( test_driver ${CALIPER_LIBRARIES} ${ADIAK_LIBRARIES} )
  endif()

  if (ENABLE_CUDA OR ENABLE_OPENMP_OFFLOAD)
    set_target_properties(test_driver PROPERTIES LINK_FLAGS "${TETON_LINK_FLAGS} ${TETON_OpenMP_Offload_LINK_FLAGS}")
    target_link_libraries(test_driver ${CUBLAS_LIBRARIES}
                                       ${CUPTI_LIBRARIES}
                                       ${NVTOOLS_LIBRARIES}
                                       ${TETON_OpenMP_Offload_LINK_LIBRARIES})
  endif()

  if (ENABLE_MFEM)
    target_include_directories( test_driver PRIVATE ${MFEM_INCLUDE_DIR} ${HYPRE_INCLUDE_DIR})
    target_link_libraries( test_driver ${MFEM_LIBRARIES} ${HYPRE_LIBRARIES} ${METIS_LIBRARIES})
  endif()

  # These are additional link libraries set in the compiler config files (typically system libraries).
  target_link_libraries( test_driver ${TETON_LINK_LIBRARIES} )

  install( TARGETS test_driver
           RUNTIME DESTINATION bin )

endif()
