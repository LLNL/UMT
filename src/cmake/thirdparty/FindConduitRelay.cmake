# Once done, this will define
#
# CONDUITRELAY_FOUND         - system has conduit relay
# CONDUITRELAY_INCLUDE_DIR   - conduit relay include directory
# CONDUITRELAY_LIBRARIES     - conduit relay library

include(FindPackageHandleStandardArgs)

find_path(
  CONDUITRELAY_INCLUDE_DIR
  NAMES conduit/conduit_relay.hpp
  PATHS ${CONDUIT_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  CONDUITRELAY_LIBRARIES
  NAMES conduit_relay
  PATHS ${CONDUIT_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    ConduitRelay
    DEFAULT_MSG
    CONDUITRELAY_LIBRARIES CONDUITRELAY_INCLUDE_DIR)

mark_as_advanced(CONDUITRELAY_LIBRARIES CONDUITRELAY_INCLUDE_DIR)
