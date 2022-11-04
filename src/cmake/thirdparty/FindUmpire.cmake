# Try to Umpire
# Once done, this will define
#
# UMPIRE_FOUND         - system has umpire
# UMPIRE_INCLUDE_DIR   - umpire include directory
# UMPIRE_LIBRARIES     - umpire library

include(FindPackageHandleStandardArgs)

find_path(
  UMPIRE_INCLUDE_DIR
  NAMES umpire/Umpire.hpp
  PATHS ${UMPIRE_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  UMPIRE_LIBRARIES
  NAMES umpire
  PATHS ${UMPIRE_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    Umpire
    DEFAULT_MSG
    UMPIRE_LIBRARIES UMPIRE_INCLUDE_DIR)

set(UMPIRE_FORTRAN_MODULES_DIR ${UMPIRE_INCLUDE_DIR}/umpire)

mark_as_advanced(UMPIRE_LIBRARIES UMPIRE_INCLUDE_DIR)
