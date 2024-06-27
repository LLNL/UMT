#include "macros.h"

module MemoryAllocator_mod

   use iso_c_binding, only : c_int, c_double, c_bool, C_SIZE_T, c_f_pointer, c_ptr, c_loc
   use Size_mod, only : Size
#if defined(TETON_ENABLE_UMPIRE)
   use umpire_mod
#endif
   implicit none

   private

   type, public :: AllocatorType

      private

#if defined(TETON_ENABLE_UMPIRE)
      ! Umpire allocators provided by host code.
      type(UmpireAllocator), public :: umpire_host_allocator
      type(UmpireAllocator), public :: umpire_device_allocator
#endif
      integer, public :: umpire_host_allocator_id = -1
      integer, public :: umpire_device_allocator_id = -1
      logical, public :: host_allocator_present = .FALSE.
      logical, public :: device_allocator_present = .FALSE.
      logical, public :: use_for_comm_data = .FALSE.

   contains

      procedure, private :: isGPUAwareMPIEnabled

      procedure, public :: construct
      procedure, public :: destruct

      procedure :: ar1 => allocate_host_real_c_double_1
      procedure :: ar2 => allocate_host_real_c_double_2
      procedure :: ar3 => allocate_host_real_c_double_3
      procedure :: ar4 => allocate_host_real_c_double_4
      procedure :: ai1 => allocate_host_integer_c_int_1
      procedure :: ai2 => allocate_host_integer_c_int_2
      procedure :: ai3 => allocate_host_integer_c_int_3
      procedure :: ai4 => allocate_host_integer_c_int_4
      procedure :: dr1 => deallocate_host_real_c_double_1
      procedure :: dr2 => deallocate_host_real_c_double_2
      procedure :: dr3 => deallocate_host_real_c_double_3
      procedure :: dr4 => deallocate_host_real_c_double_4
      procedure :: di1 => deallocate_host_integer_c_int_1
      procedure :: di2 => deallocate_host_integer_c_int_2
      procedure :: di3 => deallocate_host_integer_c_int_3
      procedure :: di4 => deallocate_host_integer_c_int_4

      generic, public :: allocate => ar1, ar2, ar3, ar4, ai1, ai2, ai3, ai4
      generic, public :: deallocate => dr1, dr2, dr3, dr4, di1, di2, di3, di4
      
   end type AllocatorType

   type(AllocatorType), pointer, public :: Allocator => null()

!-----------------------------------------------------------------------------------------
! construct
!-----------------------------------------------------------------------------------------

contains

   subroutine construct(self, umpire_host_allocator_id, umpire_device_allocator_id)
      use Options_mod, only : Options

      class(AllocatorType), intent(inout) :: self
      integer :: umpire_host_allocator_id
      integer :: umpire_device_allocator_id

#if defined(TETON_ENABLE_UMPIRE)
      type(UmpireResourceManager) :: umpire_resource_manager

      umpire_resource_manager = umpire_resource_manager%get_instance()

      if (umpire_host_allocator_id >= 0) then
         self%umpire_host_allocator = umpire_resource_manager%get_allocator_by_id(umpire_host_allocator_id)
!#if defined(TETON_ENABLE_OPENMP)
! Need a way to verify this is a thread safe umpire allocator...
!#endif
         self%umpire_host_allocator_id = umpire_host_allocator_id
         self%host_allocator_present = .TRUE.
      endif

      if (umpire_device_allocator_id >= 0) then
         self%umpire_device_allocator = umpire_resource_manager%get_allocator_by_id(umpire_device_allocator_id)
         self%umpire_device_allocator_id = umpire_device_allocator_id
         self%device_allocator_present = .TRUE.
      endif

      self%use_for_comm_data = self%isGPUAwareMPIEnabled()

      if (self%use_for_comm_data) then
         if ( Options%isRankVerbose() > 0 ) then
            print *, "Teton memory allocator: Detected env variable MPICH_GPU_SUPPORT_ENABLED=1, allocating comm data using umpire."
         endif
      endif
#endif

   end subroutine construct
   
!-----------------------------------------------------------------------------------------
! AllocatorType_destruct
!-----------------------------------------------------------------------------------------
   subroutine destruct(self)

      class(AllocatorType), intent(inout) :: self

   end subroutine destruct

!-----------------------------------------------------------------------------
! Some notes on useful Umpire allocator functions.  These all return
! integer(kind=C_SIZE_T)
!
! Get high memory watermark for this allocator.
! bytes = umpire_host/device_allocator%get_high_watermark()
!
! This is sum of the sizes of all the tracked allocations.
!  bytes = self%umpire_host/device_allocator%get_current_size()
!
! Return the actual size of this allocator.
! For non-pool allocators, this will be the same as getCurrentSize(). For pools
! this is the total amount of memory allocated for blocks managed by the pool.
! bytes = self%umpire_host/device_allocator%get_actual_size()

!-----------------------------------------------------------------------------------------
! Family of overloaded subroutines to handle real 1d, 2d, 3d array allocations
! See README.md for more information on the FTM macros.
!-----------------------------------------------------------------------------------------
!#define FTM_DEBUG

#define FTM_RANK 1
#include "MemoryAllocator_mod.F90.templates"

#define FTM_RANK 2
#include "MemoryAllocator_mod.F90.templates"

#define FTM_RANK 3
#include "MemoryAllocator_mod.F90.templates"

#define FTM_RANK 4 
#include "MemoryAllocator_mod.F90.templates"

!***********************************************************************
!    isGPUAwareMPIEnabled - Return true if gpu aware mpi support is enabled.
!
!    This functionality is for internal developer use only, for the El Cap EAS
!    nodes and should not be exposed to users.
!
!    The El Cap EAS nodes have an early MPI GPU-AWARE functionality
!    available for testing.  It currently requires setting the
!    MPICH_GPU_SUPPORT_ENABLED env var to enable it.
!
!    If Teton detects this env var, it needs to allocate its comm data
!    using umpire.  Otherwise, it should not allocate it via umpire.
!
!    NOTE: Teton can not default to allocating the comm data using umpire.
!    CPU based MPI functionality will crash if pointers allocated via hipMalloc
!    are passed to it, despite the unified memory.  If this is supported in
!    the future, this function should be removed and Teton should always use
!    the umpire allocator for comm data.
!***********************************************************************
  logical function isGPUAwareMPIEnabled(self) result(useUmpire)
    class(AllocatorType), intent(inout) :: self
    character(len=255) :: env_var_value
    integer :: env_status

    call get_environment_variable("MPICH_GPU_SUPPORT_ENABLED", env_var_value, status=env_status)
    if (env_status == 0 .AND. trim(env_var_value) == "1") then
       useUmpire = .TRUE.
    else
       useUmpire = .FALSE.
    endif

    return
  end function

end module MemoryAllocator_mod

!-----------------------------------------------------------------------------
! C interfaces 
!-----------------------------------------------------------------------------


   subroutine ConstructMemoryAllocator( umpire_host_allocator_id, &
                                        umpire_device_allocator_id ) BIND(C,NAME="teton_constructmemoryallocator")

   use iso_c_binding
   use MemoryAllocator_mod
   implicit none

   integer :: umpire_host_allocator_id
   integer :: umpire_device_allocator_id

!  Construct the Memory Allocator Module 
   allocate (Allocator)
   call Allocator%construct(umpire_host_allocator_id, umpire_device_allocator_id)

   return
   end subroutine ConstructMemoryAllocator


   subroutine DestructMemoryAllocator( ) BIND(C,NAME="teton_destructmemoryallocator")

   use iso_c_binding
   use MemoryAllocator_mod
   implicit none

!  Destruct the Memory Allocator Module
   call Allocator%destruct()
   deallocate (Allocator)

   return
   end subroutine DestructMemoryAllocator

