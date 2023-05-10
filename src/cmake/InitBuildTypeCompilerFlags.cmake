############################
# Append our default compiler flags for CMake build types.
#
# NOTE: CMake will pre-populate the CMAKE_<lang>_FLAGS_<buildtype> flags from the CMAKE_<lang>_FLAGS_<buildtype>_INIT
# variable which is set by the vendor toolchain.
############################

if (NOT "${CMAKE_BUILD_TYPE}" STREQUAL "")
	message(STATUS "Detected CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")
   string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_UPPERCASE)

   #  --- Set C++ compiler flags ---
	if ("${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPERCASE}_INIT}" STREQUAL "")
   	message(WARNING "CMake toolchain failed to initialize CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPERCASE}_INIT, please submit a ticket to vendor.")
   endif()

	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

 	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "XLClang")
      set(CMAKE_CXX_FLAGS_RELEASE "-O3 -qstrict -qarch=auto -qtune=auto -qmaxmem=-1 -qsuppress=1500-036")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -qstrict -g -qarch=auto -qtune=auto -qmaxmem=-1 -qcheck -qflttrap=enable:nanq:invalid:zerodivide -qsuppress=1500-036")
      set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -qcheck -qfullpath -qflttrap=enable:nanq:invalid:zerodivide -qsuppress=1500-036")


	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
		set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
		set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG -diag-enable=remark -fp-trap-all=common -traceback")
		set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -diag-enable=remark -fp-trap-all=common -traceback")

	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")

	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")

	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")

   endif()

   # --- Set Fortran compiler flags ---
	if ("${CMAKE_Fortran_FLAGS_${BUILD_TYPE_UPPERCASE}_INIT}" STREQUAL "")
   	message(WARNING "CMake toolchain failed to initialize CMAKE_Fortran_FLAGS_${BUILD_TYPE_UPPERCASE}_INIT, please submit a ticket to vendor.")
   endif()

	if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
      set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -ffree-line-length-none")
      set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-fcheck=all -O3 -g -DNDEBUG -ffree-line-length-none")
      set(CMAKE_Fortran_FLAGS_DEBUG "-fcheck=all -O0 -g -ffree-line-length-none")

	elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Clang") # For Clang or AppleClang

 	elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
      # Enable F2003 support via the below -qxlf2003 flag list.  This behavior is the default if xlf2003 compiler is used, but not if xlf is used.
      set(CMAKE_Fortran_FLAGS_RELEASE "-qxlf2003=autorealloc,bozlitargs,nodynamicacval,nooldnaninf,polymorphic,signdzerointr,stopexcept,volatile -O3 -qstrict -qarch=auto -qtune=auto -qlargepage -qmaxmem=-1 -qsuppress=1500-036")
      set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-qxlf2003=autorealloc,bozlitargs,nodynamicacval,nooldnaninf,polymorphic,signdzerointr,stopexcept,volatile -O3 -qstrict -g -qarch=auto -qtune=auto -qlargepage -qmaxmem=-1 -qcheck -qflttrap=enable:nanq:invalid:zerodivide -qsuppress=1500-036")
      set(CMAKE_Fortran_FLAGS_DEBUG "-qxlf2003=autorealloc,bozlitargs,nodynamicacval,nooldnaninf,polymorphic,signdzerointr,stopexcept,volatile -O0 -g -qcheck -qfullpath -qflttrap=enable:nanq:invalid:zerodivide -qsuppress=1500-036")

# There's a new feature of Fortran 2018 called IMPLICIT NONE (EXTERNAL), which, if you specify it, requires that any procedure you call have the EXTERNAL attribute,
# which you typically get from an explicit interface.  Disable the warning on this, as we're not using F2018 across all our compilers.
#
# We're missing explicit interfaces on all procedures with a dummy argument that has the ALLOCATABLE, ASYNCHRONOUS, OPTIONAL, POINTER, TARGET, VALUE or VOLATILE attribute.
# This fails the Intel language standard check.  Disable checking the interfaces for now until resolved.  See
# https://rzlc.llnl.gov/gitlab/deterministic-transport/TRT/Teton/-/issues/296 
	elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
		set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
		set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -g -warn all,noexternal,nointerfaces -diag-enable=remark -fpe-all=0 -traceback")
		set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -warn all,noexternal,nointerfaces -diag-enable=remark -check all -fpe-all=0 -traceback")

	elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")

	elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")

   # Note : Cray Fortran completely fails to provide any initial set of flags for build types.  Ticket has been submitted to HPE. --black27
   # For now, don't append to existing flags ( since there are none ), just set the optimization and debug symbols flag ourselves.
	elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
      set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
      # G2 is the only level that doesn't disable OpenMP loop collapsing and still provides debug information.
      # A ticket has been submitted to ask HPE to update the -G# flag to be consistent with the "-g" flag in their C++ compiler.
      set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -G2 -DNDEBUG")
      set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -G2")
   endif()

	message(STATUS "Build type ${CMAKE_BUILD_TYPE} CXX flags: ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPERCASE}}")
	message(STATUS "Build type ${CMAKE_BUILD_TYPE} Fortran flags: ${CMAKE_Fortran_FLAGS_${BUILD_TYPE_UPPERCASE}}")

endif()
