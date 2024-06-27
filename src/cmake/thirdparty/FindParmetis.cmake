# Once done, this will define
#
# PARMETIS_FOUND         - system has parmetis
# PARMETIS_INCLUDE_DIR   - parmetis include directory
# PARMETIS_LIBRARIES     - parmetis library

include(FindPackageHandleStandardArgs)

find_path(
  PARMETIS_INCLUDE_DIR
  NAMES parmetis.h
  PATHS ${PARMETIS_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  PARMETIS_LIBRARIES
  NAMES parmetis
  PATHS ${PARMETIS_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)

find_package_handle_standard_args(
    Parmetis
    DEFAULT_MSG
    PARMETIS_LIBRARIES PARMETIS_INCLUDE_DIR)

mark_as_advanced(PARMETIS_LIBRARIES PARMETIS_INCLUDE_DIR)
mark_as_advanced(PARMETIS_LIBRARIES)
