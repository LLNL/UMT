# Once done, this will define
#
# METIS_FOUND         - system has hypre
# METIS_INCLUDE_DIR   - hypre include directory
# METIS_LIBRARIES     - hypre library

include(FindPackageHandleStandardArgs)

find_path(
  METIS_INCLUDE_DIR
  NAMES metis.h
  PATHS ${METIS_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  METIS_LIBRARIES
  NAMES metis
  PATHS ${METIS_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)

#find_package_handle_standard_args(
#    metis
#    DEFAULT_MSG
#    METIS_LIBRARIES METIS_INCLUDE_DIR)

find_package_handle_standard_args(
    Metis
    DEFAULT_MSG
    METIS_LIBRARIES METIS_INCLUDE_DIR)

mark_as_advanced(METIS_LIBRARIES METIS_INCLUDE_DIR)
mark_as_advanced(METIS_LIBRARIES)
