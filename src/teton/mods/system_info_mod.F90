#include "macros.h"
!=======================================================================
! This module contains a variety of functions for querying system and
! hardware properties
!=======================================================================

module system_info_mod
   use iso_c_binding, only: c_int, c_size_t
   implicit none

contains

!-----------------------------------------------------------------------------------------
! Query active GPU for available memory and device id
!-----------------------------------------------------------------------------------------
subroutine getGPUMemInfo(free_bytes, total_bytes)

#if defined (TETON_ENABLE_CUDA)
   interface
      integer(kind=c_int) function gpuMemGetInfo(free_bytes, total_bytes) bind(c, name="cudaMemGetInfo")
         use iso_c_binding, only: c_int, c_size_t
         integer(kind=c_size_t), intent(inout) :: free_bytes, total_bytes
      end function gpuMemGetInfo
   end interface
#elif defined(TETON_ENABLE_HIP)
   interface
      integer(kind=c_int) function gpuMemGetInfo(free_bytes, total_bytes) bind(c, name="hipMemGetInfo")
         use iso_c_binding, only: c_int, c_size_t
         integer(kind=c_size_t), intent(inout) :: free_bytes, total_bytes
      end function gpuMemGetInfo
   end interface
#endif

   integer(kind=c_size_t), intent(inout) :: free_bytes, total_bytes
   integer(kind=c_int) :: status_code

#if defined (TETON_ENABLE_CUDA) || defined(TETON_ENABLE_HIP)
   status_code = gpuMemGetInfo(free_bytes, total_bytes)
#else
   free_bytes = 0
   total_bytes = 0
#endif

end subroutine getGPUMemInfo


subroutine printGPUMemInfo(rank)
   use mpi
   use MemoryAllocator_mod

   integer, intent(in) :: rank
   integer(kind=C_SIZE_T) :: gpu_total_bytes, gpu_free_bytes
   character(len=80) :: cuda_visible_devices
   character(len=80) :: host_name
   integer :: ierr, resultlen

   call getGPUMemInfo(gpu_free_bytes, gpu_total_bytes)
   call mpi_get_processor_name(host_name, resultlen, ierr) 

   print *, "TETON GPU used mem: rank ", rank, "host ", trim(host_name), ", currently free: ", gpu_free_bytes/(2**20), &
      "MB, currently used: ",  (gpu_total_bytes-gpu_free_bytes)/(2**20), "MB"

#if defined(TETON_ENABLE_UMPIRE)
   if (Allocator%umpire_host_allocator_id > -1) then
      print *, "TETON UMPIRE allocator size: ", Allocator%umpire_host_allocator%get_current_size() / (2**20), "MB"
   endif
   if (Allocator%umpire_device_allocator_id > -1) then
      print *, "TETON UMPIRE device allocator size: ", Allocator%umpire_device_allocator%get_current_size() / (2**20), "MB"
   endif
#endif

end subroutine printGPUMemInfo

subroutine printGPUMemRequired(rank)

   use mpi
   use MemoryAllocator_mod

#if defined(TETON_ENABLE_CUDA)
   use cuda_utils_mod
#endif

   integer, intent(in) :: rank

   integer(kind=C_SIZE_T) :: gpu_total_bytes, gpu_free_bytes, &
                             mem_estimate_bytes, cuda_bc_solver_bytes
   character(len=80) :: host_name
   integer :: ierr, resultlen

#if defined(TETON_ENABLE_UMPIRE)
   ! Use UMPIRE pinned memory allocation size as an estimator for amount of device memory needed.
   ! Skip estimating the memory usage if we are not using an umpire cpu allocator.
   if (Allocator%umpire_host_allocator_id > -1) then
      mem_estimate_bytes = Allocator%umpire_host_allocator%get_current_size()

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
#   if defined(TETON_ENABLE_CUDA)
      ! TODO - The NLsolver is not currently using UMPIRE, need to use specific call to get estimated memory.
      call getBCSolverMemEstimate(cuda_bc_solver_bytes)
      mem_estimate_bytes = mem_estimate_bytes + cuda_bc_solver_bytes
#   endif
#endif

      call getGPUMemInfo(gpu_free_bytes, gpu_total_bytes)
      call mpi_get_processor_name(host_name, resultlen, ierr) 

      print *, "TETON GPU mem estimate: rank ", rank, "host ", trim(host_name), ", Problem requires: ", &
         (mem_estimate_bytes)/(2**20), "MB.  Currently free: ", gpu_free_bytes/(2**20), "MB.  Currently used: ", &
         (gpu_total_bytes-gpu_free_bytes)/(2**20), "MB"

   endif
#endif
end subroutine printGPUMemRequired

end module system_info_mod
