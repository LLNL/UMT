# Try to find CUDA CUPTI library
# Once done, this will define
#
# CUPTI_FOUND         - system has cupti
# CUPTI_INCLUDE_DIR   - cupti include directory
# CUPTI_LIBRARIES     - cupti library

include(FindPackageHandleStandardArgs)

find_path(
  CUPTI_INCLUDE_DIR
  NAMES cupti.h
  PATHS ${CUPTI_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  CUPTI_LIBRARIES
  NAMES cupti
  PATHS ${CUPTI_ROOT}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    Cupti
    DEFAULT_MSG
    CUPTI_LIBRARIES CUPTI_INCLUDE_DIR)

mark_as_advanced(CUPTI_LIBRARIES CUPTI_INCLUDE_DIR)
