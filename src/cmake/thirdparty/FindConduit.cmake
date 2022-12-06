# Once done, this will define
#
# CONDUIT_FOUND         - system has conduit
# CONDUIT_INCLUDE_DIR   - conduit include directory
# CONDUIT_LIBRARIES     - conduit library

include(FindPackageHandleStandardArgs)

find_path(
  CONDUIT_INCLUDE_DIR
  NAMES conduit/conduit.hpp
  PATHS ${CONDUIT_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  CONDUIT_LIBRARIES
  NAMES conduit
  PATHS ${CONDUIT_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)

find_package_handle_standard_args(
    Conduit
    DEFAULT_MSG
    CONDUIT_LIBRARIES CONDUIT_INCLUDE_DIR)

set(CONDUIT_FORTRAN_MODULES_DIR ${CONDUIT_INCLUDE_DIR}/conduit)
mark_as_advanced(CONDUIT_LIBRARIES CONDUIT_INCLUDE_DIR CONDUIT_FORTRAN_MODULES_DIR)