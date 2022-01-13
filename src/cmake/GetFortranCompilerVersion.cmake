# Find compiler version using compiler version output.
# CMake is not providing the full version number for a few compilers:
# - XLF
# - Cray FTN
# This module pulls the version number by querying the compiler and
# getting the string from stdout.
#
# Sets:
# FORTRAN_COMPILER_VERSION - version of the Fortran compiler.

if (CMAKE_Fortran_COMPILER_ID STREQUAL XL)
  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} -qversion OUTPUT_VARIABLE COMPILER_VERSION_STDOUT)
  STRING(REGEX MATCH "([0-9]+)\\.[0-9]+\\.[0-9]+\\.[0-9]+" COMPILER_VERSION "${COMPILER_VERSION_STDOUT}")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} --version OUTPUT_VARIABLE COMPILER_VERSION_STDOUT)
  STRING(REGEX MATCH "([0-9]+)\\.[0-9]+\\.[0-9]+" COMPILER_VERSION "${COMPILER_VERSION_STDOUT}")
else()
  set(COMPILER_VERSION ${CMAKE_Fortran_COMPILER_VERSION})
endif()

set(Fortran_COMPILER_VERSION "${COMPILER_VERSION}")

mark_as_advanced(Fortran_COMPILER_VERSION)
