# Try to PhysicsUtils
# Once done, this will define
#
# PHYSICSUTILS_FOUND         - system has physicutils
# PHYSICSUTILS_ROOT          - physicsutils include directory
# PHYSICSUTILS_LIBRARIES     - physicsutils library

include(FindPackageHandleStandardArgs)

find_path(
  PHYSICSUTILS_INCLUDE_DIR
  NAMES PhysicsUtils/PhysicsUtilsVersion.hh
  PATHS ${PHYSICSUTILS_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  PHYSICSUTILS_LIBRARIES
  NAMES PhysicsUtils
  PATHS ${PHYSICSUTILS_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    PhysicsUtils
    DEFAULT_MSG
    PHYSICSUTILS_LIBRARIES PHYSICSUTILS_INCLUDE_DIR)

mark_as_advanced(PHYSICSUTILS_LIBRARIES PHYSICSUTILS_INCLUDE_DIR)
