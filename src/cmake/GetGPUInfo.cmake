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

if (ENABLE_OPENMP_OFFLOAD)
  message(STATUS "Checking for GPU...")
  if (CMAKE_CUDA_ARCHITECTURES STREQUAL 70)
    message(STATUS "-- Detected NVIDIA Volta, setting device num processors to 70")
    set(OMP_DEVICE_NUM_PROCESSORS 80)
    set(GSET_MIN_SIZE 16)
    set(MAX_NUM_HYPER_DOMAINS 16)
    set(OMP_DEVICE_TEAM_THREAD_LIMIT 1024)

  elseif (CMAKE_HIP_ARCHITECTURES STREQUAL gfx90a)
# AMD MI250X - 110 CUs
    message(STATUS "-- Detected AMD MI250, setting device num processors to 110")
    set(OMP_DEVICE_NUM_PROCESSORS 110)
    set(GSET_MIN_SIZE 16)
    set(MAX_NUM_HYPER_DOMAINS 16)
    set(OMP_DEVICE_TEAM_THREAD_LIMIT 1024)

# AMD MI300 - 228 CUs
  elseif (CMAKE_HIP_ARCHITECTURES STREQUAL gfx942)
    message(STATUS "-- Detected AMD MI300, setting device num processors to 228")
    set(OMP_DEVICE_NUM_PROCESSORS 228)
    set(OPENMP_UNIFIED_MEMORY TRUE)
    set(GSET_MIN_SIZE 8)
    set(MAX_NUM_HYPER_DOMAINS 16)
    set(OMP_DEVICE_TEAM_THREAD_LIMIT 1024)

  else()
    message(ERROR "-- Unrecogized or unset value for CUDA or HIP architecture.")
  endif()

else()
  message(STATUS -- No GPU detected.)
  # These are only used if running the GPU kernels on the CPU for testing purposes.
  set(OMP_DEVICE_NUM_PROCESSORS 1)
  set(OMP_DEVICE_TEAM_THREAD_LIMIT 1)
  set(GSET_MIN_SIZE 1)
  set(MAX_NUM_HYPER_DOMAINS 1)
  set(OMP_DEVICE_TEAM_THREAD_LIMIT 1)
endif()

mark_as_advanced(OMP_DEVICE_NUM_PROCESSORS OMP_DEVICE_TEAM_THREAD_LIMIT)
