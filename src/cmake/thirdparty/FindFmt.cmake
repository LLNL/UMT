# The FMT library.  Currently used by newer Umpire versions.
#
# FMT_FOUND         - system has fmt
# FMT_INCLUDE_DIR   - fmt include directory
# FMT_LIBRARIES     - fmt library

include(FindPackageHandleStandardArgs)

find_path(
  FMT_INCLUDE_DIR
  NAMES fmt/format.h
  PATHS ${FMT_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  FMT_LIBRARIES
  NAMES fmt
  PATHS ${FMT_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    Fmt
    DEFAULT_MSG
    FMT_LIBRARIES FMT_INCLUDE_DIR)

mark_as_advanced(FMT_LIBRARIES FMT_INCLUDE_DIR)
