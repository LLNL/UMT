# Example cmake cache configuration for UMT.

# Tested with gnu 10.2.1
set(TPL_ROOT "/path/to/libraries" CACHE PATH "")

set(CMAKE_C_COMPILER  "gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER  "g++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "gfortran" CACHE PATH "")

set(CMAKE_CXX_FLAGS "" CACHE PATH "")
set(CMAKE_Fortran_FLAGS "-ffree-line-length-none" CACHE PATH "")

set(ENABLE_OPENMP ON CACHE BOOL "")
set(ENABLE_OPENMP_OFFLOAD OFF CACHE BOOL "")

set(CONDUIT_ROOT ${TPL_ROOT}/conduit/0.8.2 CACHE PATH "")
