# Try to find Camp
# Once done, this will define
#
# CAMP_FOUND         - system has camp
# CAMP_INCLUDE_DIR   - camp include directory
# CAMP_LIBRARIES     - camp library

include(FindPackageHandleStandardArgs)

find_path(
  CAMP_INCLUDE_DIR
  NAMES camp/camp.hpp
  PATHS ${CAMP_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  CAMP_LIBRARIES
  NAMES camp
  PATHS ${CAMP_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)

find_package_handle_standard_args(
    Camp
    DEFAULT_MSG
    CAMP_LIBRARIES CAMP_INCLUDE_DIR)

mark_as_advanced(CAMP_LIBRARIES CAMP_INCLUDE_DIR)
