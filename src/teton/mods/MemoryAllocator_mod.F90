#include "macros.h"

module MemoryAllocator_mod

   use iso_c_binding, only : c_int, c_double, c_bool, C_SIZE_T, c_f_pointer, c_ptr, c_loc
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
      logical, public :: use_internal_host_allocator = .FALSE.

   contains

      procedure, public :: construct => AllocatorType_construct
      procedure, public :: destruct => AllocatorType_destruct

      procedure :: ar1 => allocate_host_real_c_double_1
      procedure :: ar2 => allocate_host_real_c_double_2
      procedure :: ar3 => allocate_host_real_c_double_3
      procedure :: ai1 => allocate_host_integer_c_int_1
      procedure :: ai2 => allocate_host_integer_c_int_2
      procedure :: ai3 => allocate_host_integer_c_int_3
      procedure :: dr1 => deallocate_host_real_c_double_1
      procedure :: dr2 => deallocate_host_real_c_double_2
      procedure :: dr3 => deallocate_host_real_c_double_3
      procedure :: di1 => deallocate_host_integer_c_int_1
      procedure :: di2 => deallocate_host_integer_c_int_2
      procedure :: di3 => deallocate_host_integer_c_int_3

      generic, public :: allocate => ar1, ar2, ar3, ai1, ai2, ai3
      generic, public :: deallocate => dr1, dr2, dr3, di1, di2, di3
      
   end type AllocatorType

   type(AllocatorType), pointer, public :: Allocator => null()

!-----------------------------------------------------------------------------------------
! AllocatorType_construct
!-----------------------------------------------------------------------------------------

contains

   subroutine AllocatorType_construct(self, umpire_host_allocator_id, umpire_device_allocator_id)
      use Size_mod, only: Size
      use Datastore_mod, only : theDatastore

      class(AllocatorType), intent(inout) :: self
      integer :: umpire_host_allocator_id
      integer :: umpire_device_allocator_id

#if defined(TETON_ENABLE_UMPIRE)
      type(UmpireResourceManager) :: umpire_resource_manager

      ! These are used if teton needs to create its own host pinned pool allocator.
      integer(kind=C_SIZE_T) :: initial_pool_size ! Initial size in bytes for a pool allocator
      integer(kind=C_SIZE_T) :: pool_growth_size ! Size of additional blocks to allocate, if pool needs to grow.

      initial_pool_size = 512*1024*1024 ! 512MB is default for Umpire
      pool_growth_size = 1*1024*1024 ! 1MB is default for Umpire

      umpire_resource_manager = umpire_resource_manager%get_instance()

      if (umpire_host_allocator_id >= 0) then
         self%umpire_host_allocator = umpire_resource_manager%get_allocator_by_id(umpire_host_allocator_id)
         self%umpire_host_allocator_id = umpire_host_allocator_id
      else
! If Teton is using OpenMP target offload kernel or CUDA kernels, it requires an
! umpire allocator to allocate memory in page-locked ( pinned ) memory.  If one
! was not provided then Teton will create its own.
!
! We should check the Size%useGPU variable to see if we are executing
! kernels on the GPU.  However, the Size module is not constructed until _after_
! the memory allocator module.  This needs to be refactored.  Until then, we'll
! default to constructing an allocator if Teton was compiled with gpu kernel
! support.
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
         self%use_internal_host_allocator = .TRUE.
         self%umpire_host_allocator = umpire_resource_manager%get_allocator_by_name("PINNED")
         self%umpire_host_allocator = umpire_resource_manager%make_allocator_pool("POOL", self%umpire_host_allocator, initial_pool_size, pool_growth_size)
         self%umpire_host_allocator_id = self%umpire_host_allocator%get_id()
#endif
      endif

      if (umpire_device_allocator_id >= 0) then
         self%umpire_device_allocator = umpire_resource_manager%get_allocator_by_id(umpire_device_allocator_id)
         self%umpire_device_allocator_id = umpire_device_allocator_id
      endif
#endif

      ! Make sure the datastore is initialized, as this allocator will
      ! automatically added entries to it for anything allocated.
      call theDatastore%initialize()

   end subroutine AllocatorType_construct
   
!-----------------------------------------------------------------------------------------
! AllocatorType_destruct
!-----------------------------------------------------------------------------------------
   subroutine AllocatorType_destruct(self)

      class(AllocatorType), intent(inout) :: self

#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
      !if (use_internal_host_allocator) then
      !  TODO Ask Umpire how to teardown allocators.
      !  Need to tear down the pinned pool allocator.
      !endif
#endif

   end subroutine AllocatorType_destruct

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

