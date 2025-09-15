# Try to find RAJA
# Once done, this will define
#
# RAJA_FOUND         - system has raja
# RAJA_INCLUDE_DIR   - raja include directory
# RAJA_LIBRARIES     - raja library

include(FindPackageHandleStandardArgs)

find_path(
  RAJA_INCLUDE_DIR
  NAMES RAJA/RAJA.hpp
  PATHS ${RAJA_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  RAJA_LIBRARIES
  NAMES RAJA
  PATHS ${RAJA_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    RAJA
    DEFAULT_MSG
    RAJA_LIBRARIES RAJA_INCLUDE_DIR)

mark_as_advanced(RAJA_LIBRARIES RAJA_INCLUDE_DIR)
