# Try to Umpire
# Once done, this will define
#
# NVTOOLS_FOUND         - system has nvtools
# NVTOOLS_INCLUDE_DIR   - nvtools include directory
# NVTOOLS_LIBRARIES     - nvtools library

include(FindPackageHandleStandardArgs)

find_path(
  NVTOOLS_INCLUDE_DIR
  NAMES nvToolsExt.h
  PATHS ${NVTOOLS_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  NVTOOLS_LIBRARIES
  NAMES nvToolsExt
  PATHS ${NVTOOLS_ROOT}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)

find_package_handle_standard_args(
    NVTools
    DEFAULT_MSG
    NVTOOLS_LIBRARIES NVTOOLS_INCLUDE_DIR)

mark_as_advanced(NVTOOLS_LIBRARIES NVTOOLS_INCLUDE_DIR)
