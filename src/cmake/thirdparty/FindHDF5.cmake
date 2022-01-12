# Try to Umpire
# Once done, this will define
#
# HDF5_FOUND         - system has hdf5
# HDF5_INCLUDE_DIR   - hdf5 include directory
# HDF5_LIBRARIES     - hdf5 library

include(FindPackageHandleStandardArgs)

find_path(
  HDF5_INCLUDE_DIR
  NAMES hdf5.h
  PATHS ${HDF5_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  HDF5_LIBRARIES
  NAMES hdf5
  PATHS ${HDF5_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_package_handle_standard_args(
    HDF5
    DEFAULT_MSG
    HDF5_LIBRARIES HDF5_INCLUDE_DIR)

mark_as_advanced(HDF5_LIBRARIES HDF5_INCLUDE_DIR)
