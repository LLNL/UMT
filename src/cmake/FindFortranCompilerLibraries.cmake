# Find location of Fortran compiler libraries.
# Fortran_COMPILER_LIB_DIR   - path to Fortran compiler libraries
# Fortran_COMPILER_LIB_DIR   - "-L path to Fortran compiler libraries"

if (CMAKE_Fortran_COMPILER_ID STREQUAL XL)
  # Calling the compiler with no input file will result in:
  # /usr/tce/packages/xl/xl-2020.11.12/xlf/16.1.1/bin/.orig/xlf90_r: 1501-294 (S) No input file specified. Please use -qhelp for more information.
  # We parse this line to grab the XL package location.
  EXECUTE_PROCESS( COMMAND ${CMAKE_Fortran_COMPILER} ERROR_VARIABLE OUTPUT_STDERR)
  STRING(REGEX MATCH "/usr/tce/packages/xl/xl-([0-9]+)\\.[0-9]+\\.[0-9]+" XL_PACKAGE_DIR "${OUTPUT_STDERR}")
  set(Fortran_COMPILER_LIB_DIR "${XL_PACKAGE_DIR}/alllibs")
  set(Fortran_COMPILER_LIB_DIR_LINK "-L${Fortran_COMPILER_LIB_DIR}")
else()
  set(Fortran_COMPILER_LIB_DIR "")
  set(Fortran_COMPILER_LIB_DIR_LINK "")
endif()

mark_as_advanced(Fortran_COMPILER_LIB_DIR Fortran_COMPILER_LIB_DIR_LINK)
