# Determine some characteristics of the OpenMP target device (gpu).
# 
# For now, make some assumptions based on the compiler we are using.
# - XL - assume we are on Sierra
# - Cray - assume we are on El Cap
#
# Sets:
# This is the number of 'streaming multiprocessors' or 'compute units'.
# OMP_TARGET_NUM_PROCESSORS
# This is the number of threads supported per thread team/ block on the target device.
# OMP_TARGET_MAX_THREADS_PER_THREAD_TEAM

if (CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
  set(OMP_TARGET_NUM_PROCESSORS 80)
  set(OMP_TARGET_MAX_THREADS_PER_THREAD_TEAM 1024)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(OMP_TARGET_NUM_PROCESSORS 120)
  set(OMP_TARGET_MAX_THREADS_PER_THREAD_TEAM 1024)
else()
  set(OMP_TARGET_NUM_PROCESSORS 1)
  set(OMP_TARGET_MAX_THREADS_PER_THREAD_TEAM 1)
endif()

mark_as_advanced(OMP_TARGET_NUM_PROCESSORS OMP_TARGET_MAX_THREADS_PER_THREAD_TEAM)
