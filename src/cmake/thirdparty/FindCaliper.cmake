# Try to Umpire
# Once done, this will define
#
# CALIPER_FOUND         - system has caliper
# CALIPER_INCLUDE_DIR   - caliper include directory
# CALIPER_LIBRARIES     - caliper library

include(FindPackageHandleStandardArgs)

find_path(
  CALIPER_INCLUDE_DIR
  NAMES caliper/Caliper.h
  PATHS ${CALIPER_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  CALIPER_LIBRARIES
  NAMES caliper
  PATHS ${CALIPER_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    Caliper
    DEFAULT_MSG
    CALIPER_LIBRARIES CALIPER_INCLUDE_DIR)

mark_as_advanced(CALIPER_LIBRARIES CALIPER_INCLUDE_DIR)
