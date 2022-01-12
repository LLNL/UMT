# Try to Umpire
# Once done, this will define
#
# Z_FOUND         - system has z
# Z_INCLUDE_DIR   - z include directory
# Z_LIBRARIES     - z library

include(FindPackageHandleStandardArgs)

find_path(
  Z_INCLUDE_DIR
  NAMES zlib.h
  PATHS ${Z_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  Z_LIBRARIES
  NAMES z
  PATHS ${Z_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    Z
    DEFAULT_MSG
    Z_LIBRARIES Z_INCLUDE_DIR)

mark_as_advanced(Z_LIBRARIES Z_INCLUDE_DIR)
