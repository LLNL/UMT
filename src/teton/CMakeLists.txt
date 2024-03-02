include(../cmake/SubDirList.cmake)

# Quiet the warning that the relative source path for INTERFACE
# will be converted to absolute path.
cmake_policy(SET CMP0076 NEW)

############################
# LIBRARY
############################

# add teton library target
add_library ( teton STATIC )

if (ENABLE_FIND_MPI)
  target_link_libraries( teton PUBLIC MPI::MPI_Fortran MPI::MPI_CXX )
endif()

if( ENABLE_CUDA )
  add_library( tetoncuda OBJECT )
  target_link_libraries( teton PRIVATE tetoncuda )

  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.18.0")
     target_link_libraries( teton PUBLIC CUDA::cublas
                                         CUDA::cupti
                                         CUDA::nvToolsExt
                          )
  endif()
endif()

if( ENABLE_OPENMP )
  target_link_libraries( teton PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries( teton PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Use CMake generator to only apply this flag on Fortran files. - black27
# TODO - Remove OpenACC support as soon as HPE has performant OpenMP in the CCE compiler.
if( ENABLE_OPENACC )
  target_compile_options(teton PUBLIC "$<$<COMPILE_LANGUAGE:Fortran>:${OpenACC_Fortran_FLAGS}>")
endif()

target_sources( teton INTERFACE include/TetonInterface.hh )

if (ENABLE_UMPIRE)
   target_include_directories( teton PRIVATE include
                               ${UMPIRE_INCLUDE_DIR}
                               ${UMPIRE_FORTRAN_MODULES_DIR}
                             )
endif()

# Note:
# Caliper installs its Fortran modules down in include/caliper/fortran,
# so we must add that additional include path
if (ENABLE_CALIPER)
   target_include_directories( teton PRIVATE include
                               ${CALIPER_INCLUDE_DIR}
                               ${CALIPER_INCLUDE_DIR}/caliper/fortran
                             )
endif()

if (ENABLE_BLUEPRINT_INTERFACE)
  target_sources( teton INTERFACE include/TetonSources.hh )
  target_sources( teton INTERFACE include/TetonBlueprint.hh )
  target_sources( teton INTERFACE include/TetonSurfaceTallies.hh )
  target_sources( teton INTERFACE include/TetonConduitInterface.hh )
endif()

# Note:
# The MFEM conduit data collection class expects to have an include
# path to the internal 'include/conduit' directory, so add that
# additional include path until MFEM can fix that behavior.
target_include_directories( teton PRIVATE include
                            ${CONDUIT_INCLUDE_DIR}
                            ${CONDUIT_INCLUDE_DIR}/conduit
                            ${CONDUIT_FORTRAN_MODULES_DIR})

# Note:
# Some installations of physicsutils may have the Fortran modules in
# include/physicsutils instead of the base include directory.
if (NOT ENABLE_MINIAPP_BUILD)
   target_include_directories( teton PRIVATE include
                               ${PHYSICSUTILS_INCLUDE_DIR}
                               ${PHYSICSUTILS_INCLUDE_DIR}/PhysicsUtils
                               ${PHYSICSUTILS_INCLUDE_DIR}/fortran
                               ${SILO_INCLUDE_DIR}
                             )
endif()

configure_file (include/TetonVersion.hh.in include/TetonVersion.hh )

target_sources( teton INTERFACE ${PROJECT_BINARY_DIR}/teton/include/TetonVersion.hh )
target_include_directories( teton PUBLIC include ${PROJECT_BINARY_DIR}/teton/include )

install( FILES "include/TetonInterface.hh" DESTINATION include )

if (ENABLE_BLUEPRINT_INTERFACE)
   install( FILES "include/TetonConduitInterface.hh" DESTINATION include )
   install( FILES "include/TetonBlueprint.hh" DESTINATION include )
   install( FILES "include/TetonSurfaceTallies.hh" DESTINATION include )
   install( FILES "include/TetonSources.hh" DESTINATION include )
endif()

install( FILES "${PROJECT_BINARY_DIR}/teton/include/TetonVersion.hh" DESTINATION include )

install( TARGETS teton ARCHIVE DESTINATION lib)

# Process any sub-dirs with CMakeList.txt files.
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_LIST_DIR})

foreach(subdir ${SUBDIRS})
  message( STATUS "Adding files from ${subdir}")
  add_subdirectory(${subdir})
endforeach()

################################################
# TEST DRIVER EXECUTABLE
################################################

if(ENABLE_TESTS)

  if (NOT ENABLE_BLUEPRINT_INTERFACE)
    message( ERROR "Blueprint interface support must be enabled for the test driver build.")
  endif()

  add_executable( test_driver
                  driver/test_driver.cc )

  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang")
    set_target_properties( test_driver PROPERTIES LINK_FLAGS -fno-fortran-main)
    set(OpenMP_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -fno-fortran-main")
  endif()

  if (ENABLE_FIND_MPI)
    target_link_libraries( test_driver PUBLIC MPI::MPI_Fortran MPI::MPI_CXX )
  endif()

  if( ENABLE_CUDA )
    target_link_libraries( test_driver PRIVATE tetoncuda )
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.18.0")
      target_link_libraries( test_driver PRIVATE CUDA::cudart
                                          CUDA::cublas
                                          CUDA::cupti
                                          CUDA::nvToolsExt
                            )
    endif()
  endif()

  # Generally preferred to link with the Fortran compiler.
  # Intel has trouble with this, fallback on using the C++ compiler to link.
  if ( NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL IntelLLVM)
    set_target_properties(test_driver PROPERTIES LINKER_LANGUAGE Fortran)
  endif ()

  target_include_directories( test_driver PUBLIC
                              ../src/include
                              ${CONDUIT_INCLUDE_DIR}
                              ${CONDUIT_INCLUDE_DIR}/conduit
                            )

  target_link_libraries( test_driver PUBLIC
                         teton
                         ${CONDUIT_LIBRARIES}
                         ${CONDUITRELAY_LIBRARIES}
                         ${CONDUITBLUEPRINT_LIBRARIES}
                         ${CONDUITBLUEPRINTMPI_LIBRARIES}
                         ${CONDUITRELAYMPI_LIBRARIES}
                         ${CONDUITRELAYMPIIO_LIBRARIES}
                       )

# TODO - Remove OpenACC support as soon as HPE has performant OpenMP in the CCE compiler.
  if( ENABLE_OPENACC )
    set_target_properties(test_driver PROPERTIES LINK_FLAGS ${OpenACC_Fortran_FLAGS})
  endif()

  if( ENABLE_OPENMP )
    target_link_libraries( test_driver PUBLIC OpenMP::OpenMP_Fortran OpenMP::OpenMP_CXX)

    # The target_link_options command has trouble with adding flags if there is
    # a space, it keeps putting quotes around the whole line.  Use the older
    # LINK_FLAGS property instead.
    set_target_properties( test_driver PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS})
    #target_link_options(test_driver PRIVATE ${OpenMP_Fortran_FLAGS})

  endif()

  if(ENABLE_SILO)
     target_link_libraries( test_driver PUBLIC
                            ${SILO_LIBRARIES}
                          )
  endif()

  # This is only enabled if one of the TPLs require HDF5, so it can be added to the link line.
  # Our code does not have any direct dependencies on HDF5.
  if (ENABLE_HDF5)
     target_link_libraries( test_driver PUBLIC
                            ${HDF5_LIBRARIES}
                            ${Z_LIBRARIES}
                            dl
                          )
  endif()

  if (NOT ENABLE_MINIAPP_BUILD)
    target_link_libraries( test_driver PUBLIC ${PHYSICSUTILS_LIBRARIES} )
  endif()

  if (ENABLE_UMPIRE)
    target_include_directories( test_driver PUBLIC ${UMPIRE_INCLUDE_DIR} 
                                                    ${UMPIRE_FORTRAN_MODULES_DIR})
    target_link_libraries( test_driver PUBLIC ${UMPIRE_LIBRARIES} )
    if (ENABLE_CAMP)
      target_include_directories( test_driver PUBLIC ${CAMP_INCLUDE_DIR})
      target_link_libraries( test_driver PUBLIC ${CAMP_LIBRARIES} )
    endif()
  endif()

  if (ENABLE_CALIPER)
    target_include_directories( test_driver PUBLIC ${CALIPER_INCLUDE_DIR} ${ADIAK_INCLUDE_DIR})
    target_link_libraries( test_driver PUBLIC ${CALIPER_LIBRARIES} ${ADIAK_LIBRARIES} rt)
  endif()

  if (ENABLE_MFEM)
    target_include_directories( test_driver PUBLIC ${MFEM_INCLUDE_DIR} ${HYPRE_INCLUDE_DIR})
    target_link_libraries( test_driver PUBLIC ${MFEM_LIBRARIES} ${HYPRE_LIBRARIES} ${METIS_LIBRARIES})
  endif()

  target_link_options(test_driver PUBLIC "LINKER:${TETON_LINK_OPTIONS}")

  install( TARGETS test_driver
           RUNTIME DESTINATION bin )

# Executable for creating unstructured box mesh.
  if (ENABLE_MFEM)
    add_executable( makeUnstructuredBox
                    driver/makeUnstructuredBox.cc )
    if (ENABLE_FIND_MPI)
      target_link_libraries( makeUnstructuredBox PUBLIC MPI::MPI_CXX )
    endif()
    if (ENABLE_HDF5)
        target_link_libraries( makeUnstructuredBox PUBLIC
                               ${HDF5_LIBRARIES}
                               ${Z_LIBRARIES}
                               dl
                             )
    endif()

# Note:
# The MFEM conduit data collection class expects to have an include
# path to the internal 'include/conduit' directory, so add that
# additional include path until MFEM can fix that behavior.
    target_include_directories( makeUnstructuredBox PUBLIC ${CONDUIT_INCLUDE_DIR}
                                                           ${CONDUIT_INCLUDE_DIR}/conduit
                                                           ${MFEM_INCLUDE_DIR}
                                                           ${HYPRE_INCLUDE_DIR}
                              )
    target_link_options( makeUnstructuredBox PUBLIC "LINKER:${TETON_LINK_OPTIONS}")
    target_link_libraries( makeUnstructuredBox PUBLIC ${MFEM_LIBRARIES} ${HYPRE_LIBRARIES} ${METIS_LIBRARIES})
    install( TARGETS makeUnstructuredBox
             RUNTIME DESTINATION bin )
  endif()

endif()
