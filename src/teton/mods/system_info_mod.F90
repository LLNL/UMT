#include "macros.h"
!=======================================================================
! This module contains a variety of functions for querying system and
! hardware properties
!=======================================================================

module system_info_mod
   use, intrinsic :: iso_c_binding, only : C_SIZE_T, C_INT
   implicit none

contains

!-----------------------------------------------------------------------------------------
! Query active GPU for available memory and device id
!-----------------------------------------------------------------------------------------
subroutine getGPUMemInfo(free_bytes, total_bytes)

#if defined (TETON_ENABLE_CUDA)
   use cudafor, only : cudaMemGetInfo, cudaGetDevice
#endif

   integer(kind=C_SIZE_T), intent(inout) :: free_bytes, total_bytes
   integer :: status_code

#if defined (TETON_ENABLE_CUDA)
   status_code = cudaMemGetInfo(free_bytes, total_bytes)
#else
   free_bytes = 0
   total_bytes = 0
#endif

end subroutine getGPUMemInfo


subroutine printGPUMemInfo(rank)
   use mpi

   integer, intent(in) :: rank
   integer(kind=C_SIZE_T) :: gpu_total_bytes, gpu_free_bytes
   character(len=80) :: cuda_visible_devices
   character(len=80) :: host_name
   integer :: ierr, resultlen

   call getGPUMemInfo(gpu_free_bytes, gpu_total_bytes)
   call get_environment_variable("CUDA_VISIBLE_DEVICES", cuda_visible_devices)
   call mpi_get_processor_name(host_name, resultlen, ierr) 

   print *, "TETON GPU used mem rpt: rank ", rank, "host ", trim(host_name), ", CUDA_VISIBLE_DEVICES:", trim(cuda_visible_devices), &
     ", currently free: ", gpu_free_bytes/(2**20), "MB, currently used: ",  (gpu_total_bytes-gpu_free_bytes)/(2**20), "MB"

end subroutine printGPUMemInfo

subroutine printGPUMemRequired(rank)
   use mpi
   use MemoryAllocator_mod

#if defined(TETON_ENABLE_CUDA)
   use cuda_utils_mod
#endif

   integer, intent(in) :: rank

   integer(kind=C_SIZE_T) :: gpu_total_bytes, gpu_free_bytes, &
                             umpire_cpu_allocator_used_bytes, nlsolver_estimated_bytes
   character(len=80) :: cuda_visible_devices
   character(len=80) :: host_name
   integer :: ierr, resultlen

   nlsolver_estimated_bytes = 0
   umpire_cpu_allocator_used_bytes = 0

   ! Use UMPIRE pinned memory allocation size as an estimator for amount of device memory needed.
#if defined(TETON_ENABLE_UMPIRE)
   umpire_cpu_allocator_used_bytes = Allocator%umpire_host_allocator%get_current_size()
#endif

#if defined(TETON_ENABLE_CUDA)
#  if !defined(TETON_ENABLE_MINIAPP_BUILD)
   ! TODO - The NLsolver is not currently using UMPIRE, need to use specific call to get estimated memory.
   call getBCSolverMemEstimate(nlsolver_estimated_bytes)
#  endif
#endif

   call getGPUMemInfo(gpu_free_bytes, gpu_total_bytes)
   call get_environment_variable("CUDA_VISIBLE_DEVICES", cuda_visible_devices)
   call mpi_get_processor_name(host_name, resultlen, ierr) 

   print *, "TETON GPU mem estimate rpt: rank ", rank, "host ", trim(host_name), ", CUDA_VISIBLE_DEVICES:", trim(cuda_visible_devices), &
     ", Problem requires: Sweep ", (umpire_cpu_allocator_used_bytes)/(2**20), "MB + Scattering solver ", (nlsolver_estimated_bytes)/(2**20), &
     "MB = Total ", (umpire_cpu_allocator_used_bytes+nlsolver_estimated_bytes)/(2**20), "MB", &
     ", currently free: ", gpu_free_bytes/(2**20), "MB, currently used: ",  (gpu_total_bytes-gpu_free_bytes)/(2**20), "MB"

end subroutine printGPUMemRequired

end module system_info_mod
