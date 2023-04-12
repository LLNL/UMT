# Once done, this will define
#
# MFEM_FOUND         - system has hypre
# MFEM_INCLUDE_DIR   - hypre include directory
# MFEM_LIBRARIES     - hypre library

include(FindPackageHandleStandardArgs)

find_path(
  MFEM_INCLUDE_DIR
  NAMES mfem.hpp
  PATHS ${MFEM_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  MFEM_LIBRARIES
  NAMES mfem
  PATHS ${MFEM_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)

find_package_handle_standard_args(
    MFEM
    DEFAULT_MSG
    MFEM_LIBRARIES MFEM_INCLUDE_DIR)

mark_as_advanced(MFEM_LIBRARIES MFEM_INCLUDE_DIR)
