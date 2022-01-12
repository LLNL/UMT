# Try to find Adiak
# Once done, this will define
#
# ADIAK_FOUND         - system has adiak
# ADIAK_INCLUDE_DIR   - adiak include directory
# ADIAK_LIBRARIES     - adiak library

include(FindPackageHandleStandardArgs)

find_path(
  ADIAK_INCLUDE_DIR
  NAMES adiak.hpp
  PATHS ${ADIAK_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  ADIAK_LIBRARIES
  NAMES adiak
  PATHS ${ADIAK_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    Adiak
    DEFAULT_MSG
    ADIAK_LIBRARIES ADIAK_INCLUDE_DIR)

mark_as_advanced(ADIAK_LIBRARIES ADIAK_INCLUDE_DIR)
