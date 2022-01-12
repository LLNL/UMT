# Try to find Memusages
# Once done, this will define
#
# MEMUSAGES_FOUND         - system has memusages
# MEMUSAGES_INCLUDE_DIR   - memusages include directory
# MEMUSAGES_LIBRARIES     - memusages library

include(FindPackageHandleStandardArgs)

find_path(
  MEMUSAGES_INCLUDE_DIR
  NAMES MemUsages.h
  PATHS ${MEMUSAGES_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  MEMUSAGES_LIBRARIES
  NAMES MemUsages
  PATHS ${MEMUSAGES_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)

find_package_handle_standard_args(
    Memusages
    DEFAULT_MSG
    MEMUSAGES_LIBRARIES MEMUSAGES_INCLUDE_DIR)

mark_as_advanced(MEMUSAGES_LIBRARIES MEMUSAGES_INCLUDE_DIR)
