# Once done, this will define
#
# CONDUITBLUEPRINTMPI_FOUND         - system has conduit blueprint
# CONDUITBLUEPRINTMPI_INCLUDE_DIR   - conduit blueprint include directory
# CONDUITBLUEPRINTMPI_LIBRARIES     - conduit blueprint library

include(FindPackageHandleStandardArgs)

if(TETON_BUILDING_WITH_PARMETIS)
   # Find Conduit using the Conduit installed package so we can get any library
   # dependencies for the conduit_blueprint_mpi library.
   unset(CONDUIT_FOUND)
   find_package(Conduit PATHS ${CONDUIT_ROOT})
   set(CONDUIT_FOUND TRUE)
   get_target_property(deps conduit_blueprint_mpi INTERFACE_LINK_LIBRARIES)
   set(CONDUITBLUEPRINTMPI_LIBRARIES conduit_blueprint_mpi;${deps})
   message(STATUS "CONDUITBLUEPRINTMPI_LIBRARIES=${CONDUITBLUEPRINTMPI_LIBRARIES}")

else()
   find_path(
     CONDUITBLUEPRINTMPI_INCLUDE_DIR
     NAMES conduit/conduit_blueprint_mpi.hpp
     PATHS ${CONDUIT_ROOT}
     PATH_SUFFIXES include
     NO_DEFAULT_PATH
     NO_CMAKE_ENVIRONMENT_PATH
     NO_CMAKE_PATH
     NO_SYSTEM_ENVIRONMENT_PATH
     NO_CMAKE_SYSTEM_PATH
   ) 

   find_library(
     CONDUITBLUEPRINTMPI_LIBRARIES
     NAMES conduit_blueprint_mpi
     PATHS ${CONDUIT_ROOT}
     PATH_SUFFIXES lib
     NO_DEFAULT_PATH
     NO_CMAKE_ENVIRONMENT_PATH
     NO_CMAKE_PATH
     NO_SYSTEM_ENVIRONMENT_PATH
     NO_CMAKE_SYSTEM_PATH
   )

   find_package_handle_standard_args(
       ConduitBlueprintMPI
       DEFAULT_MSG
       CONDUITBLUEPRINTMPI_LIBRARIES CONDUITBLUEPRINTMPI_INCLUDE_DIR)

   mark_as_advanced(CONDUITBLUEPRINTMPI_LIBRARIES CONDUITBLUEPRINTMPI_INCLUDE_DIR)

endif()
