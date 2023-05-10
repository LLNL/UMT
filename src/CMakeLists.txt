#<<<<<=====----------------------------------------------------------=====>>>>>
#                      CMake script for Teton library
#<<<<<=====----------------------------------------------------------=====>>>>>
cmake_minimum_required( VERSION 3.14 )

set(PROJECT_NAME teton)
# Update version information in this file.

set(TETON_VERSION_MAJOR 5)
set(TETON_VERSION_MINOR 1)
set(TETON_VERSION_PATCH 0)
include (cmake/Version.cmake)

# The CUDA Boltzmann Compton solver source file requires C++11.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# ----- Set up Teton project source files and build targets -----
project( ${PROJECT_NAME} LANGUAGES CXX Fortran VERSION ${TETON_VERSION_MAJOR}.${TETON_VERSION_MINOR}.${TETON_VERSION_PATCH})

# Tell CMake to not clear out the RPATH when installing an executable.
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Add third party FindXXX.cmake modules.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/thirdparty/")

# Set up our default flags for each CMake build type.
include(cmake/InitBuildTypeCompilerFlags.cmake)
# Set correct linker wrapper flags.  This is a workaround for HPE CCE compiler not setting this
# for their Fortran compiler.  ( Ticket submitted to HPE  -black27).
include(cmake/InitLinkerWrapperFlags.cmake)

# ----- Set up compiler flags and defines based on user provided options -----

# Limits use of preprocessor macros to conform with strict FPP implementations.
option( STRICT_FPP_MODE "Use preprocessor macros that strictly conform to FPP." NO)
if(STRICT_FPP_MODE)
   message( STATUS "Use preprocessor macros that strictly conform to FPP.")
   add_compile_definitions( "TETON_STRICT_FPP_MODE")
endif()

# -----
# Enable mini-app code version that focuses on the transport sweep algorithm.
# This will enable the test driver, and MFEM mesh support by default.
# -----
option( ENABLE_MINIAPP_BUILD "Enable transport sweep mini-app." TRUE )
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

# Only allows OpenACC on Cray compiler, as this is only meant to be an interim workaround for OpenMP performance issues.
# TODO - Remove OpenACC support as soon as HPE has performant OpenMP in the CCE compiler.
option( ENABLE_OPENACC "Enable OpenACC pragmas" NO )
if(ENABLE_OPENACC)
   if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
      message( STATUS "Enable OpenACC pragmas" )
      set(OpenACC_Fortran_VERSION 2.0)
      set(OpenACC_Fortran_SPEC_DATE 201306)
      set(OpenACC_Fortran_FLAGS -hacc)

      message( STATUS "Compiler supports OpenACC ${OpenACC_Fortran_VERSION}" )
      message( STATUS "Compiler supports OpenACC spec date ${OpenACC_Fortran_SPEC_DATE}")
      message( STATUS "Compiler OpenACC Fortran flags: ${OpenACC_Fortran_FLAGS}")

      add_compile_definitions( "TETON_ENABLE_OPENACC" )
   else()
      message(ERROR "OpenACC is only supported for Cray compiler.")
   endif()
endif()

option( ENABLE_OPENMP "Enable OpenMP pragmas" NO )
option( ENABLE_OPENMP_OFFLOAD "Enable OpenMP target offload pragmas" NO )
option( OPENMP_HAS_USE_DEVICE_ADDR "Compiler OpenMP implementation includes use_device_addr pragma" NO )
option( OPENMP_HAS_FORTRAN_INTERFACE "Compiler OpenMP implementation includes Fortran interface" NO )

# TODO - add check that OpenMP_Fortran_HAVE_OMPLIB_MODULE is TRUE, since we use the omp_lib module.
if(ENABLE_OPENMP)
   message( STATUS "Enable OpenMP pragmas" )

   if(ENABLE_OPENMP_OFFLOAD)
     message( STATUS "Enable OpenMP target offload pragmas" )

     if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "XL" )
       message(STATUS "Detected XL compiler, manually adding -qoffload to OpenMP flags...")
       set(OpenMP_Fortran_FLAG "-qsmp=omp -qoffload")
     endif()

     if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "XLClang")
       set(OpenMP_CXX_FLAG "-qsmp=omp -qoffload")
     endif()

     add_compile_definitions( "TETON_ENABLE_OPENMP_OFFLOAD" )

     if(OPENMP_HAS_USE_DEVICE_ADDR)
       message( STATUS "Enable use of OpenMP use_device_addr pragmas" )
       add_compile_definitions( "TETON_OPENMP_HAS_USE_DEVICE_ADDR" )
     endif()
   endif()

   find_package(OpenMP REQUIRED COMPONENTS Fortran CXX)

   message( STATUS "Compiler supports OpenMP ${OpenMP_Fortran_VERSION}" )
   message( STATUS "Compiler supports OpenMP spec date ${OpenMP_Fortran_SPEC_DATE}")
   message( STATUS "Compiler OpenMP C++ flags: ${OpenMP_CXX_FLAGS}")
   message( STATUS "Compiler OpenMP Fortran flags: ${OpenMP_Fortran_FLAGS}")

   add_compile_definitions( "TETON_ENABLE_OPENMP" )

   if(OPENMP_HAS_FORTRAN_INTERFACE)
      message( STATUS "OpenMP includes Fortran interface, disable creation of iso_c_bindings." )
      add_compile_definitions( "TETON_OPENMP_HAS_FORTRAN_INTERFACE" )
   endif()
endif()

option( ENABLE_CUDA "Enable compiling CUDA code" NO )
if( ENABLE_CUDA )
   message( STATUS "Enable CUDA kernels" )
   enable_language( CUDA )
   set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
   add_compile_definitions( "TETON_ENABLE_CUDA" )

   if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.18.0") 
      find_package(CUDAToolkit REQUIRED)
   else()
      message( STATUS "Older CMake version detected (pre-3.18.0).  CUDA compiler flags will need to be manually added by cmake config file or on command line.")
   endif()

endif()

# CMake fails to find the MPI libraries and flags on some mixed vendor compiler configurations, for
# example Clang + XLF on BLUEOS.  In these cases, disable the cmake find package for mpi and instead
# use the mpi compiler wrappers.
# Requires that MPI_CXX_COMPILER and MPI_Fortran_COMPILER are set.
option(ENABLE_FIND_MPI "Enables CMake's Find MPI support (Turn off when compiling with the mpi wrapper directly)" ON)
if (ENABLE_FIND_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran CXX)
else()
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

option( ENABLE_TOMP_TRACES "Enable print traces before any OpenMP target pragmas." NO )
if(ENABLE_TOMP_TRACES)
   message( STATUS "Enabling OpenMP target pragma print traces" )
   add_compile_definitions( "TETON_ENABLE_TOMP_TRACES" )
endif()

option( ENABLE_TESTS "Enable test driver executable." ${ENABLE_MINIAPP_BUILD} )
if (ENABLE_TESTS)
   message( STATUS "Enable teton test driver." )
   option( ENABLE_MFEM "Enable mfem mesh support." ${ENABLE_MINIAPP_BUILD} )
   message( STATUS "ENABLE_MFEM is ${ENABLE_MFEM}")
   if (ENABLE_MFEM)
      message( STATUS "Enable mfem mesh support.")
      add_compile_definitions( "TETON_ENABLE_MFEM" )
   endif()
else()
   message( STATUS "Disable teton test driver." )
endif()


# ----- Gather gpu hardware details, will be set in Teton build_info module.
include (cmake/GetGPUInfo.cmake)

option( ENABLE_UMPIRE "Link and enable the use of the UMPIRE memory library." NO )
if(ENABLE_UMPIRE)
   find_package(Umpire)
   if(ENABLE_CAMP)
      find_package(Camp)
   endif()
   add_compile_definitions( "TETON_ENABLE_UMPIRE" )
endif()

option( ENABLE_CALIPER "Link and enable the use of the CALIPER performance measurement library." NO )
if(ENABLE_CALIPER)
   find_package(Caliper)
   find_package(Adiak)
   add_compile_definitions( "TETON_ENABLE_CALIPER" )
endif()

option( ENABLE_BLUEPRINT_INTERFACE "Enable support for configuring teton via conduit blueprint interface." YES )

find_package(Conduit)
find_package(ConduitRelay)
find_package(ConduitRelayMPI)
find_package(ConduitRelayMPIIO)
find_package(ConduitBlueprint)
find_package(ConduitBlueprintMPI)

# Needed for restart SILO file support.
option( ENABLE_SILO "Link and enable the use of the SILO file i/o library." YES)
if (ENABLE_SILO)
   find_package(Silo)
endif()

if(NOT ENABLE_MINIAPP_BUILD)
   # Need the source headers for these, regardless of whether we compile the test driver executable.
   find_package(PhysicsUtils)
endif()

if (ENABLE_TESTS)
   # Only needed if we compile the test driver mfem support.
   if (ENABLE_MFEM)
      find_package(MFEM)
      find_package(Hypre)
      find_package(Metis)
   endif()
endif()

# Needed if conduit was built against HDF5.
option( ENABLE_HDF5 "Link in the HDF5 library." NO)
if (ENABLE_HDF5)
   find_package(HDF5)
   find_package(Z)
endif()

add_subdirectory(teton)

message(STATUS "User specified CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "User specific Fortran Flags: ${CMAKE_Fortran_FLAGS}")
