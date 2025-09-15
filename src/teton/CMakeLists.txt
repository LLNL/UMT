include(../cmake/SubDirList.cmake)

# Quiet the warning that the relative source path for INTERFACE
# will be converted to absolute path.
cmake_policy(SET CMP0076 NEW)

############################
# LIBRARY
############################

# add teton library target
add_library ( teton STATIC )

# For AMD C++ and HPE CC Fortran builds that have OpenMP target kernels, the vendor recommendation is to link with the CCE Fortran compiler,
# but linking with CCE Fortran can be difficult to implement on large code bases across multiple libraries.
#
# If the AMD C++ compiler is used to link, some extra steps are needed to manually inject OpenMP runtime symbols before and after the
# library.
#
# Experimental support is provided in the CMake below to perform these steps as a post_build step.
#
# This functionally depends on:
# - the location of the cce_omp_offload_linker being in a known location relative to the Fortran compiler.
# - the linker tool requiring only these three env vars to be set:
# -- CRAY_OMP_TOOL_PATH
# -- CRAY_BINUTILS_BIN
# -- ROCM_PATH
#
# A ticket has been created for AMD and HPE requesting that the above items be added as command line args to
# the linker to reduce dependencies on env vars and the fragility of this functionality.

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")

   if (INJECT_OPENMP_OFFLOAD_SYMBOLS_EXPERIMENTAL AND ENABLE_OPENMP_OFFLOAD)
      message(STATUS Detected Fortran compiler id  ${CMAKE_Fortran_COMPILER_ID})
      message(STATUS Detected Fortran compiler name  ${CMAKE_Fortran_COMPILER_VERSION})
      message(STATUS Detected Fortran compiler path  ${CMAKE_Fortran_COMPILER})

      if (NOT DEFINED HIP_ROOT_DIR)
         message( FATAL_ERROR " Must set HIP_ROOT_DIR when using Cray compiler openmp offload functionality.")
      endif()

      set(CRAY_OMP_TOOL_PATH /opt/cray/pe/cce/${CMAKE_Fortran_COMPILER_VERSION}/cce/x86_64/bin/cce_omp_offload_linker)
      if (NOT EXISTS ${CRAY_OMP_TOOL_PATH})
         message (FATAL_ERROR " Cray OMP offload tool not found at: ${CRAY_OMP_TOOL_PATH}.  You may need to set the CMake var CRAY_OMP_TOOL_PATH manually.")
      endif()

      set(CRAY_BINUTILS_BIN /opt/cray/pe/cce/${CMAKE_Fortran_COMPILER_VERSION}/binutils/x86_64/bin)
      if (NOT EXISTS ${CRAY_BINUTILS_BIN})
         message (FATAL_ERROR " Cray bintuils bin path not found at: ${CRAY_BINUTILS_BIN}.  You may need to set the CMake var CRAY_BINUTILS_BIN manually.")
      endif()

      set(CRAY_PE_USE_CLANG /opt/cray/pe/cce/${CMAKE_Fortran_COMPILER_VERSION}/cce-clang/x86_64/bin/clang)
      if (NOT EXISTS ${CRAY_PE_USE_CLANG})
         message (FATAL_ERROR " Cray clang not found at: ${CRAY_PE_USE_CLANG}.  You may need to set the CMake var CRAY_PE_USE_CLANG manually.")
      endif()

      add_custom_command( TARGET teton POST_BUILD
         COMMAND ROCM_PATH=${ROCM_PATH} CRAY_BINUTILS_BIN=${CRAY_BINUTILS_BIN} CRAY_PE_USE_CLANG=${CRAY_PE_USE_CLANG} ${CRAY_OMP_TOOL_PATH} --arch ${CMAKE_HIP_ARCHITECTURES} --device-link-only --output-begin=begin.o --output-end=end.o -- libteton.a
         COMMAND ${CMAKE_COMMAND} -E make_directory tmp_obj_files
         COMMAND ${CMAKE_COMMAND} -E rename libteton.a tmp_obj_files/libteton.original.a
         COMMAND ${CMAKE_COMMAND} -E chdir tmp_obj_files ${CMAKE_AR} x libteton.original.a
         COMMAND ${CMAKE_AR} rc libteton.a begin.o tmp_obj_files/*.o end.o
         COMMAND ${CMAKE_COMMAND} -E remove_directory tmp_obj_files
         COMMAND ${CMAKE_COMMAND} -E remove begin.o
         COMMAND ${CMAKE_COMMAND} -E remove end.o
         COMMENT "Executed post-build script"
      )
   endif()
endif()

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
  target_link_libraries( teton PUBLIC OpenMP::OpenMP_Fortran)

# For mixed vendor builds ( example: Amd clang and Cray crayftn ) the openmp libraries
# can conflict.  Do not add the omp runtimes from both C++ and Fortran.
# We'll add the needed Fortran omp runtime but for the C++ only add the compiler flags.
  get_target_property(imported_cpp_compile_options OpenMP::OpenMP_CXX INTERFACE_COMPILE_OPTIONS)
  target_compile_options(teton PRIVATE ${imported_cpp_compile_options})
endif()

target_sources( teton INTERFACE include/TetonInterface.hh )

if (ENABLE_UMPIRE)
   target_include_directories( teton PRIVATE include
                               ${UMPIRE_INCLUDE_DIR}
                               ${UMPIRE_FORTRAN_MODULES_DIR}
                             )
   if (FMT_ROOT)
     target_include_directories( teton PUBLIC ${FMT_INCLUDE_DIR})
   endif()

   if (CAMP_ROOT)
     target_include_directories( teton PUBLIC ${CAMP_INCLUDE_DIR})
   endif()
endif()

# Note:
# Caliper installs its Fortran modules down in include/caliper/fortran,
# so we must add that additional include path
if (ENABLE_CALIPER)
   if (ENABLE_FIND_PACKAGE_CONFIG_MODE)
     target_link_libraries( teton PUBLIC caliper )
   else()
     target_include_directories( teton PRIVATE include
                                 ${CALIPER_INCLUDE_DIR}
                                 ${CALIPER_INCLUDE_DIR}/caliper/fortran
                               )
   endif()

   if (ADIAK_ROOT)
      if (ENABLE_FIND_PACKAGE_CONFIG_MODE)
        target_link_libraries( teton PUBLIC adiak::adiak )
      else()
        target_include_directories( teton PRIVATE include
                                    ${ADIAK_INCLUDE_DIR}
                                  )
      endif()
   endif()
endif()


target_sources( teton INTERFACE include/TetonSources.hh )
target_sources( teton INTERFACE include/TetonBlueprint.hh )
target_sources( teton INTERFACE include/TetonSurfaceTallies.hh )
target_sources( teton INTERFACE include/TetonConduitInterface.hh )
target_sources( teton INTERFACE include/TetonNDAccessor.hh )
target_sources( teton INTERFACE include/TetonTesting.hh )
target_sources( teton INTERFACE include/TetonUtilities.hh )
target_sources( teton INTERFACE include/TetonModulesCInterfaces.hh )
target_sources( teton INTERFACE include/processEnvVars.hh )

target_include_directories( teton PRIVATE include
                            ${CONDUIT_INCLUDE_DIR}
                            ${CONDUIT_FORTRAN_MODULES_DIR})

# Note:
# Some installations of physicsutils may have the Fortran modules in
# include/physicsutils instead of the base include directory.
if (NOT ENABLE_MINIAPP_BUILD)
   target_include_directories( teton PRIVATE include
                               ${PHYSICSUTILS_INCLUDE_DIR}
                               ${PHYSICSUTILS_INCLUDE_DIR}/PhysicsUtils
                               ${SILO_INCLUDE_DIR}
                             )

   # There was short lived bug in physicsutils when it was putting the fortran modules into 'fortran'.
   IF(EXISTS ${PHYSICSUTILS_INCLUDE_DIR}/fortran)
      target_include_directories( teton PRIVATE include
                               ${PHYSICSUTILS_INCLUDE_DIR}/fortran
                                )
   endif()

   if (RAJA_ROOT)
      target_include_directories( teton PRIVATE include
                                    ${RAJA_INCLUDE_DIR}
                                  )
   endif()
endif()

configure_file (include/TetonVersion.hh.in include/TetonVersion.hh )

target_sources( teton INTERFACE ${PROJECT_BINARY_DIR}/teton/include/TetonVersion.hh )
target_include_directories( teton PUBLIC include ${PROJECT_BINARY_DIR}/teton/include )

install( FILES "include/TetonInterface.hh" DESTINATION include )

install( FILES "include/TetonConduitInterface.hh" DESTINATION include )
install( FILES "include/TetonBlueprint.hh" DESTINATION include )
install( FILES "include/TetonSurfaceTallies.hh" DESTINATION include )
install( FILES "include/TetonSources.hh" DESTINATION include )

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

  add_executable( test_driver
                  driver/test_driver.cc )

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

  if( ENABLE_OPENMP )
    target_link_libraries( test_driver PUBLIC OpenMP::OpenMP_Fortran)

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
  if (ENABLE_HDF5 OR CONDUIT_REQUIRES_HDF5)
     target_link_libraries( test_driver PUBLIC
                            ${HDF5_LIBRARIES}
                            dl
                          )
  endif()

  if (NOT ENABLE_MINIAPP_BUILD)
    target_link_libraries( test_driver PUBLIC ${PHYSICSUTILS_LIBRARIES} )
    if (RAJA_ROOT)
      target_link_libraries( test_driver PUBLIC ${RAJA_LIBRARIES})
    endif()
  endif()

  if (ENABLE_UMPIRE)
    target_include_directories( test_driver PUBLIC ${UMPIRE_INCLUDE_DIR} 
                                                   ${UMPIRE_FORTRAN_MODULES_DIR})
    target_link_libraries( test_driver PUBLIC ${UMPIRE_LIBRARIES})

    if (FMT_ROOT)
      target_include_directories( test_driver PUBLIC ${FMT_INCLUDE_DIR})
      target_link_libraries( test_driver PUBLIC ${FMT_LIBRARIES} )
    endif()

    if (CAMP_ROOT)
      target_include_directories( test_driver PUBLIC ${CAMP_INCLUDE_DIR})
      target_link_libraries( test_driver PUBLIC ${CAMP_LIBRARIES} )
    endif()
  endif()

  if (ENABLE_CALIPER)
    if (ENABLE_FIND_PACKAGE_CONFIG_MODE)
      target_link_libraries( test_driver PUBLIC caliper )
      if (ADIAK_ROOT)
	      target_link_libraries( test_driver PUBLIC adiak::adiak )
      endif()

    else()
      target_include_directories( test_driver PUBLIC ${CALIPER_INCLUDE_DIR})
      target_link_libraries( test_driver PUBLIC ${CALIPER_LIBRARIES} rt)

      if (ADIAK_ROOT)
	      target_include_directories( test_driver PUBLIC ${ADIAK_INCLUDE_DIR})
	      target_link_libraries( test_driver PUBLIC ${ADIAK_LIBRARIES})
      endif()
    endif()

    if (ENABLE_HIP)
      target_link_libraries( test_driver PUBLIC ${ROCTRACER_LIBRARIES})
    endif()

  endif()


  if (CONDUIT_REQUIRES_PARMETIS)
    target_link_libraries( test_driver PUBLIC ${PARMETIS_LIBRARIES} ${METIS_LIBRARIES})
  endif()

  if (Z_ROOT)
     target_link_libraries( test_driver PUBLIC
                            ${Z_LIBRARIES}
                          )
  endif()

  if (TETON_TEST_DRIVER_LIBRARIES)
    target_link_libraries( test_driver PUBLIC ${TETON_TEST_DRIVER_LIBRARIES})
  endif()

  if (TETON_TEST_DRIVER_LIBRARIES)
    target_link_libraries( test_driver PUBLIC ${TETON_TEST_DRIVER_LIBRARIES})
  endif()

  target_link_options(test_driver PUBLIC "LINKER:${TETON_LINK_OPTIONS}")

  install( TARGETS test_driver
           RUNTIME DESTINATION bin )

endif()
