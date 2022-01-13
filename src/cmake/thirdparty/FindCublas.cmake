# Try to Umpire
# Once done, this will define
#
# CUBLAS_FOUND         - system has cublas
# CUBLAS_INCLUDE_DIR   - cublas include directory
# CUBLAS_LIBRARIES     - cublas library

include(FindPackageHandleStandardArgs)

find_path(
  CUBLAS_INCLUDE_DIR
  NAMES cublas.h
  PATHS ${CUBLAS_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  CUBLAS_LIBRARIES
  NAMES cublas
  PATHS ${CUBLAS_ROOT}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    Cublas
    DEFAULT_MSG
    CUBLAS_LIBRARIES CUBLAS_INCLUDE_DIR)

mark_as_advanced(CUBLAS_LIBRARIES CUBLAS_INCLUDE_DIR)
