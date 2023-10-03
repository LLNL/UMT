# Example cmake cache configuration for UMT.

# Tested with gnu 10.2.1
set(TPL_ROOT "/path/to/libraries" CACHE PATH "")

set(CMAKE_C_COMPILER  "nvcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER  "nvcc" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "nvfortran" CACHE PATH "")

set(CMAKE_CXX_FLAGS "" CACHE PATH "")
set(CMAKE_Fortran_FLAGS "" CACHE PATH "")

set(ENABLE_OPENMP ON CACHE BOOL "")
set(ENABLE_OPENMP_OFFLOAD OFF CACHE BOOL "")

set(CONDUIT_ROOT ${TPL_ROOT}/conduit/0.8.2 CACHE PATH "")
set(MFEM_ROOT ${TPL_ROOT}/mfem/4.4 CACHE PATH "")
set(HYPRE_ROOT ${TPL_ROOT}/hypre/2.24.0 CACHE PATH "")
set(METIS_ROOT ${TPL_ROOT}/metis/5.1.0 CACHE PATH "")
