# Example cmake cache configuration for UMT.

# Tested with gnu 10.2.1
set(TPL_ROOT "/path/to/libraries" CACHE PATH "")

set(MPI_C_COMPILER  "mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER  "mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "mpif90" CACHE PATH "")

set(CMAKE_CXX_FLAGS "" CACHE PATH "")
set(CMAKE_Fortran_FLAGS "-ffree-line-length-none" CACHE PATH "")

set(TETON_OpenMP_Fortran_FLAGS_RELEASE "-fopenmp" CACHE STRING "")
set(TETON_OpenMP_Fortran_FLAGS_DEBUG "-fopenmp" CACHE STRING "")
set(TETON_OpenMP_CXX_FLAGS_RELEASE "-fopenmp" CACHE STRING "")
set(TETON_OpenMP_CXX_FLAGS_DEBUG "-fopenmp" CACHE STRING "")

set(ENABLE_OPENMP ON CACHE BOOL "")
set(ENABLE_OPENMP_OFFLOAD OFF CACHE BOOL "")

set(CONDUIT_ROOT ${TPL_ROOT}/conduit/0.8.2 CACHE PATH "")
set(MFEM_ROOT ${TPL_ROOT}/mfem/4.4 CACHE PATH "")
set(HYPRE_ROOT ${TPL_ROOT}/hypre/2.24.0 CACHE PATH "")
set(METIS_ROOT ${TPL_ROOT}/metis/5.1.0 CACHE PATH "")
