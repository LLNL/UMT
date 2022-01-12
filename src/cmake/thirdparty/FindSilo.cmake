# Try to Umpire
# Once done, this will define
#
# SILO_FOUND         - system has silo
# SILO_INCLUDE_DIR   - silo include directory
# SILO_LIBRARIES     - silo library

include(FindPackageHandleStandardArgs)

find_path(
  SILO_INCLUDE_DIR
  NAMES silo.h
  PATHS ${SILO_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  SILO_LIBRARIES
  NAMES siloh5
  PATHS ${SILO_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    Silo
    DEFAULT_MSG
    SILO_LIBRARIES SILO_INCLUDE_DIR)

mark_as_advanced(SILO_LIBRARIES SILO_INCLUDE_DIR)
