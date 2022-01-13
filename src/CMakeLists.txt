#<<<<<=====----------------------------------------------------------=====>>>>>
#                      CMake script for Teton library
#<<<<<=====----------------------------------------------------------=====>>>>>
cmake_minimum_required( VERSION 3.13.3 )

set(PROJECT_NAME teton)
# Update version information in this file.

set(TETON_VERSION_MAJOR 4)
set(TETON_VERSION_MINOR 12)
set(TETON_VERSION_PATCH 0)
include (cmake/Version.cmake)

# The CUDA Boltzmann Compton solver source file requires C++11.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

set(TETON_INTEL_VER_REQUIRED 19.0.4.20190416)
set(TETON_GNU_VER_REQUIRED 8.1.0)
set(TETON_XL_VER_REQUIRED 16.01.0001.0009)
set(TETON_CRAY_VER_REQUIRED 11.0.3)

# Tell CMake to not clear out the RPATH when installing an executable.
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if( NOT DEFINED CMAKE_BUILD_TYPE ) # Set a default build type (i.e. Release, Debug)
   set( CMAKE_BUILD_TYPE "Release" CACHE STRING CMAKE_BUILD_TYPE )
endif()

# ----- Set up compiler flags and defines based on user provided options -----

# -----
# NOTE: Do this BEFORE defining the Teton project below.  We want these compiler
# flags included in the CMake compiler tests so the right compiler libraries are
# found.
# -----

# If MPI_lang_COMPILER is set use that, otherwise assume CMAKE_lang_COMPILER is
# an mpi compiler.
if(DEFINED MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_Fortran_COMPILER} CACHE STRING CMAKE_C_COMPILER FORCE)
endif()
if(DEFINED MPI_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE STRING CMAKE_CXX_COMPILER FORCE)
endif()
if(DEFINED MPI_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} CACHE STRING CMAKE_Fortran_COMPILER FORCE)
endif()

# Limits use of preprocessor macros to conform with stricter FPP implementations.
option( STRICT_FPP_MODE "Use preprocessor macros that strictly conform to FPP." NO)
if(STRICT_FPP_MODE)
   message( STATUS "Use preprocessor macros that strictly conform to FPP.")
   add_compile_definitions( "TETON_STRICT_FPP_MODE")
endif()

# -----
# Enable mini-app code version that focuses on the transport sweep algorithm.
# -----
option( ENABLE_MINIAPP_BUILD "Enable transport sweep mini-app." YES )
if(ENABLE_MINIAPP_BUILD)
   message( STATUS "Enabling transport sweep mini-app build." )
   add_compile_definitions( "TETON_ENABLE_MINIAPP_BUILD" )
endif()

# Enable some extra function hooks around the intensity iteration for testing, manufactured solutions, etc.
option( ENABLE_INTENSITY_SOLVE_PREPOST_HOOKS "Enable additional intensity iteration pre/post function hooks" NO )
if(ENABLE_INTENSITY_SOLVE_PREPOST_HOOKS)
   message( STATUS "Enable pre/post function hooks around intensity solve" )
   add_compile_definitions( "TETON_ENABLE_INTENSITY_SOLVE_PREPOST_HOOKS" )
endif()

option( ENABLE_OPENMP "Enable OpenMP pragmas" NO )
option( ENABLE_OPENMP_OFFLOAD "Enable OpenMP target offload pragmas" NO )
option( OPENMP_HAS_USE_DEVICE_ADDR "Compiler OpenMP implementation includes use_device_addr pragma" NO )
option( OPENMP_HAS_FORTRAN_INTERFACE "Compiler OpenMP implementation includes Fortran interface" NO )

if(ENABLE_OPENMP)
   message( STATUS "Enable OpenMP pragmas" )
   add_compile_definitions( "TETON_ENABLE_OPENMP" )

   if(OPENMP_HAS_FORTRAN_INTERFACE)
      message( STATUS "OpenMP includes Fortran interface, disable creation of iso_c_bindings." )
      add_compile_definitions( "TETON_OPENMP_HAS_FORTRAN_INTERFACE" )
   endif()

   if(ENABLE_OPENMP_OFFLOAD)
      message( STATUS "Enable OpenMP target offload pragmas" )
      add_compile_definitions( "TETON_ENABLE_OPENMP_OFFLOAD" )

      if(OPENMP_HAS_USE_DEVICE_ADDR)
         message( STATUS "Enable use of OpenMP use_device_addr pragmas" )
         add_compile_definitions( "TETON_OPENMP_HAS_USE_DEVICE_ADDR" )
      endif()
   endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
   message( STATUS "Configuring Release build type.")
   if(ENABLE_OPENMP)
      set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${TETON_OpenMP_Fortran_FLAGS_RELEASE}" CACHE STRING CMAKE_Fortran_FLAGS FORCE )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TETON_OpenMP_CXX_FLAGS_RELEASE}" CACHE STRING CMAKE_CXX_FLAGS FORCE )
   endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
   message( STATUS "Configuring Debug build type.")
   if(ENABLE_OPENMP)
      set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${TETON_OpenMP_Fortran_FLAGS_DEBUG}" CACHE STRING CMAKE_Fortran_FLAGS FORCE )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TETON_OpenMP_CXX_FLAGS_DEBUG}" CACHE STRING CMAKE_CXX_FLAGS FORCE )
   endif()
   add_compile_definitions( "TETON_CHECK_OUT_OF_BOUNDS_ARRAY_ACCESSES" )
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
   message( STATUS "Configuring Release with Debug info build type.")
   if(ENABLE_OPENMP)
# Use the release flags for OpenMP.  Sierra doesn't support debug symbols on optimized GPU kernels.
      set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${TETON_OpenMP_Fortran_FLAGS_RELEASE}" CACHE STRING CMAKE_Fortran_FLAGS FORCE )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TETON_OpenMP_CXX_FLAGS_RELEASE}" CACHE STRING CMAKE_CXX_FLAGS FORCE )
   endif()
   add_compile_definitions( "TETON_CHECK_OUT_OF_BOUNDS_ARRAY_ACCESSES" )
else()
   if(ENABLE_OPENMP)
      message(FATAL_ERROR "Missing OpenMP flags for build type: ${CMAKE_BUILD_TYPE}, check Teton CMakeLists.txt.")
   endif()
endif()


option( ENABLE_CUDA "Enable compiling CUDA code" NO )
if( ENABLE_CUDA )
   message( STATUS "Enable CUDA kernels" )
   enable_language( CUDA )

   if( NOT DEFINED CUDA_HOME ) # Read the CUDA_HOME environment variable
      execute_process( COMMAND printenv CUDA_HOME OUTPUT_VARIABLE CUDA_HOME )
      string( STRIP ${CUDA_HOME} CUDA_HOME ) # CUDA_HOME not specified, use env var CUDA_HOME
   endif()
   message( STATUS "Using CUDA at ${CUDA_HOME}" )
   set( CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} )
   set( CMAKE_CUDA_COMPILER ${CUDA_HOME}/bin/nvcc CACHE STRING CMAKE_CUDA_COMPILER )
   set( CMAKE_CUDA_FLAGS "-ccbin ${CMAKE_CUDA_HOST_COMPILER} -Xcompiler '${CMAKE_CXX_FLAGS}'" )
   include_directories( ${CUDA_HOME}/include )
   set( CMAKE_Fortran_FLAGS " ${Fortran_CUDA_FLAGS} ${CMAKE_Fortran_FLAGS}" CACHE STRING CMAKE_Fortran_FLAGS FORCE )
   add_compile_definitions( "TETON_ENABLE_CUDA" )
endif()

option( ENABLE_TOMP_TRACES "Enable print traces before any OpenMP target pragmas." NO )
if(ENABLE_TOMP_TRACES)
   message( STATUS "Enabling OpenMP target pragma print traces" )
   add_compile_definitions( "TETON_ENABLE_TOMP_TRACES" )
endif()

option( ENABLE_TESTS "Enable test driver executable." OFF )
if (ENABLE_TESTS)
   message( STATUS "Enable teton test driver." )
   option( ENABLE_MFEM"Enable mfem mesh support." OFF )
   if (ENABLE_MFEM)
      message( STATUS "Enable mfem mesh support.")
      add_compile_definitions( "TETON_ENABLE_MFEM" )
   endif()
else()
   message( STATUS "Disable teton test driver." )
endif()

# ----- Set up Teton project source files and build targets -----
project( ${PROJECT_NAME} LANGUAGES CXX Fortran VERSION ${TETON_VERSION_MAJOR}.${TETON_VERSION_MINOR}.${TETON_VERSION_PATCH})


# ----- Gather gpu hardware details, will be set in Teton build_info module.
include (cmake/GetGPUInfo.cmake)

# ----- Check if any additional Fortran library directories are needed at link time.
#       These are used if compiling the test driver.
include(cmake/FindFortranCompilerLibraries.cmake)


# ----- Check if tested compiler is being used. -----
option(ENABLE_COMPILER_VERSION_CHECK "Check that compiler meets minimum tested version." YES)
if (ENABLE_COMPILER_VERSION_CHECK)
   include(cmake/GetFortranCompilerVersion.cmake)
   message(STATUS "Checking ${CMAKE_Fortran_COMPILER_ID} ${Fortran_COMPILER_VERSION} against minimum supported version...")
   if (CMAKE_Fortran_COMPILER_ID STREQUAL "Gnu" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS ${TETON_GNU_VER_REQUIRED})
      message(FATAL_ERROR "GNU gfortran ${CMAKE_Fortran_COMPILER_VERSION} detected.  This version is untested.  Minimum tested compiler version is ${TETON_GNU_VER_REQUIRED}.")
   endif()

   if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS ${TETON_INTEL_VER_REQUIRED})
      message(FATAL_ERROR "Intel ifort ${CMAKE_Fortran_COMPILER_VERSION} detected.  This version is untested.  Minimum tested compiler version is ${TETON_INTEL_VER_REQUIRED}.")
   endif()
   if (CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
      if (${Fortran_COMPILER_VERSION} VERSION_LESS ${TETON_XL_VER_REQUIRED})
         message(FATAL_ERROR "IBM xlf ${Fortran_COMPILER_VERSION} detected.  This version is untested.  Minimum tested compiler version is ${TETON_XL_VER_REQUIRED}.")
      endif()
   endif()
   if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
      include(cmake/GetFortranCompilerVersion.cmake)
      if (${Fortran_COMPILER_VERSION} VERSION_LESS ${TETON_CRAY_VER_REQUIRED})
         message(FATAL_ERROR "Cray ftn ${Fortran_COMPILER_VERSION} detected.  This version is untested.  Minimum tested compiler version is ${TETON_CRAY_VER_REQUIRED}.")
      endif()
   endif()
else()
    message( STATUS "Minimum compiler version check is disabled." )
endif()

# Add third party FindXXX.cmake modules.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/thirdparty/")

option( ENABLE_UMPIRE "Link and enable the use of the UMPIRE memory library." NO )
if(ENABLE_UMPIRE)
   find_package(Umpire)
   add_compile_definitions( "TETON_ENABLE_UMPIRE" )
endif()

option( ENABLE_MEMUSAGES "Link and enable the use of the MEMUSAGES to track memory usage in test driver." YES )
if(ENABLE_MEMUSAGES)
   find_package(Memusages)
   add_compile_definitions( "TETON_ENABLE_MEMUSAGES" )
endif()

option( ENABLE_CALIPER "Link and enable the use of the CALIPER performance measurement library." NO )
if(ENABLE_CALIPER)
   find_package(Caliper)
   find_package(Adiak)
   add_compile_definitions( "TETON_ENABLE_CALIPER" )
endif()

find_package(Conduit)
find_package(ConduitRelay)
find_package(ConduitBlueprint)

if(NOT ENABLE_MINIAPP_BUILD)
   # Need the source headers for these, regardless of whether we compile the test driver executable.
   find_package(PhysicsUtils)
endif()

option( ENABLE_SILO "Link and enable the use of the SILO file i/o library." YES )
if (ENABLE_SILO)
  find_package(Silo)
endif()

if (ENABLE_TESTS)
   # Only needed if we compile the test driver executable.
   find_package(HDF5)
   find_package(SZ)
   find_package(Z)
   # Only needed if we compile the test driver mfem support.
   if (ENABLE_MFEM)
      find_package(MFEM)
      find_package(Hypre)
      find_package(Metis)
   endif()
endif()

if(ENABLE_CUDA)
  find_package(NVTools)
  find_package(Cublas)
  if(ENABLE_CALIPER)
    find_package(Cupti)
  endif()
endif()

add_subdirectory(teton)
