# This is an example CMake configuration cache file for UMT.

# It expects that you have required 3rd party libraries for UMT installed in a directory structure like:
# $TPL_ROOT/ library_name / version
# Set TPL_ROOT to where your libraries are installed.

set(MPI_C_COMPILER  "mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER  "mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "mpif90" CACHE PATH "")

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall" CACHE PATH "")
set(CMAKE_Fortran_FLAGS "-ffree-line-length-none -Wall" CACHE PATH "")

set(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -g" CACHE PATH "")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g" CACHE PATH "")

set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g" CACHE PATH "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE PATH "")

set(TETON_OpenMP_Fortran_FLAGS_RELEASE "-fopenmp" CACHE STRING "")
set(TETON_OpenMP_Fortran_FLAGS_DEBUG "-fopenmp" CACHE STRING "")
set(TETON_OpenMP_CXX_FLAGS_RELEASE "-fopenmp" CACHE STRING "")
set(TETON_OpenMP_CXX_FLAGS_DEBUG "-fopenmp" CACHE STRING "")

set(ENABLE_OPENMP OFF CACHE BOOL "")
set(ENABLE_OPENMP_OFFLOAD OFF CACHE BOOL "")
set(ENABLE_CUDA OFF CACHE BOOL "")
set(ENABLE_CALIPER OFF CACHE BOOL "")
set(ENABLE_MEMUSAGES OFF CACHE BOOL "")
set(ENABLE_MFEM ON CACHE BOOL "")
set(ENABLE_TESTS ON CACHE BOOL "")
set(ENABLE_UMPIRE OFF CACHE BOOL "")
set(ENABLE_SILO ON CACHE BOOL "")
set(ENABLE_PHYSICSUTILS OFF CACHE BOOL "")
set(ENABLE_ADIAK OFF CACHE BOOL "")

# These libraries are required compile time dependencies
set(SILO_ROOT ${TPL_ROOT}/silo/4.10.3 CACHE PATH "")
set(CONDUIT_ROOT ${TPL_ROOT}/conduit/0.8.0 CACHE PATH "")

# These libraries are link time dependencies.
set(MFEM_ROOT ${TPL_ROOT}/mfem/4.3 CACHE PATH "")
set(HDF5_ROOT ${TPL_ROOT}/hdf5/1.8.10p1 CACHE PATH "")
set(SZ_ROOT ${TPL_ROOT}/bzip2/1.0.8 CACHE PATH "")
set(Z_ROOT ${TPL_ROOT}/zlib/1.2.11 CACHE PATH "")
