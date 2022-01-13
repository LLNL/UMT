# Try to Umpire
# Once done, this will define
#
# SZ_FOUND         - system has sz
# SZ_INCLUDE_DIR   - sz include directory
# SZ_LIBRARIES     - sz library

include(FindPackageHandleStandardArgs)

find_path(
  SZ_INCLUDE_DIR
  NAMES szlib.h
  PATHS ${SZ_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  SZ_LIBRARIES
  NAMES sz
  PATHS ${SZ_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    SZ
    DEFAULT_MSG
    SZ_LIBRARIES SZ_INCLUDE_DIR)

mark_as_advanced(SZ_LIBRARIES SZ_INCLUDE_DIR)
