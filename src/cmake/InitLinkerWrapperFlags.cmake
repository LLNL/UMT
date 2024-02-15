############################
# Setup linker wrapper flags, for compilers that are failing to do this.
#
# In these instances, we should submit a ticket to the compiler vendor, but until then this module provides a workaround.
# For more information see
# https://cmake.org/cmake/help/latest/command/target_link_options.html
# specifically at the section on 
# Handling Compiler Driver Differences
############################

get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)

if ("C" IN_LIST languages AND "${CMAKE_C_LINKER_WRAPPER_FLAG}" STREQUAL "")
  message(FATAL_ERROR "Toolchain file failed to initialize CMAKE_C_LINKER_WRAPPER_FLAG, please submit a ticket to vendor.")
endif()

if ("CXX" IN_LIST languages AND "${CMAKE_CXX_LINKER_WRAPPER_FLAG}" STREQUAL "")
  message(FATAL_ERROR "Toolchain file failed to initialize CMAKE_CXX_LINKER_WRAPPER_FLAG, please submit a ticket to vendor.")
endif()

if ("Fortran" IN_LIST languages AND "${CMAKE_Fortran_LINKER_WRAPPER_FLAG}" STREQUAL "")

  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
      message(WARNING "Cray Fortran toolchain failed to set CMAKE_Fortran_LINKER_WRAPPER_FLAGS.  Report this to vendor.  Setting flag to '-Wl,' manually as workaround.")
      set(CMAKE_Fortran_LINKER_WRAPPER_FLAG "-Wl,")
      set(CMAKE_Fortran_LINKER_WRAPPER_FLAG_SEP ",")
  elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang")
         message(WARNING "LLVM Flang toolchain failed to set CMAKE_Fortran_LINKER_WRAPPER_FLAGS.  Report this to vendor.  Setting flag to '-Wl,' manually as workaround.")
      set(CMAKE_Fortran_LINKER_WRAPPER_FLAG "-Wl,")
      set(CMAKE_Fortran_LINKER_WRAPPER_FLAG_SEP ",")
  else()
      message(FATAL_ERROR "Toolchain file failed to initialize CMAKE_Fortran_LINKER_WRAPPER_FLAG, please submit a ticket to vendor.")
  endif()
endif()
