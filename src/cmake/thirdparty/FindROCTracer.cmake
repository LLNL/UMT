# Try to find roctracer library
# Once done, this will define
#
# ROCTRACER_FOUND         - system has roctracer
# ROCTRACER_INCLUDE_DIR   - roctracer include directory
# ROCTRACER_LIBRARIES     - roctracer library
#
# AMD does not provided an exported roctracer target. See issue
# https://rzlc.llnl.gov/jira/browse/ELCAP-578

include(FindPackageHandleStandardArgs)

find_path(
  ROCTRACER_INCLUDE_DIR
  NAMES roctracer/roctracer.h
  PATHS ${HIP_ROOT_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
) 

find_library(
  ROCTRACER_LIBRARY
  NAMES roctracer64
  PATHS ${HIP_ROOT_DIR}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)


find_library(
  ROCTX_LIBRARY
  NAMES roctx64
  PATHS ${HIP_ROOT_DIR}
  PATH_SUFFIXES lib lib64
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
)

if(ROCTRACER_LIBRARY AND ROCTX_LIBRARY)
    set(ROCTRACER_LIBRARIES ${ROCTRACER_LIBRARY} ${ROCTX_LIBRARY})
    set(ROCTRACER_FOUND TRUE)
    message(STATUS "Found ROCtracer: ${ROCTRACER_LIBRARIES}")
    mark_as_advanced(ROCTRACER_LIBRARIES ROCTRACER_LIBRARY ROCTX_LIBRARY ROCTRACER_INCLUDE_DIR)
endif()

find_package_handle_standard_args(
    ROCTracer
    DEFAULT_MSG
    ROCTRACER_LIBRARIES ROCTRACER_INCLUDE_DIR)
