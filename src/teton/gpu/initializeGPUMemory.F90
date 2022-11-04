#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  01/2018, PFN                    *
!                                                                      *
!   initializeGPUMemory  - Initializes GPU memory that resides for     *
!                          entire radiation step.                      *
!                                                                      *
!***********************************************************************   

   subroutine initializeGPUMemory

   use kind_mod
   use QuadratureList_mod
   use Size_mod
   use SetData_mod
   use OMPWrappers_mod
   use MemoryAllocator_mod

   implicit none

!  Local 

   integer :: setID
   integer :: nSets
   integer :: err_code

!  Constants

   nSets = getNumberOfSets(Quad)

   do setID=1,nSets

     if (Size% ndim == 2) then
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID) % PsiM)
     endif

     if (Allocator%umpire_device_allocator_id >= 0) then

       ! TODO - Move this documentation to an appropriate place.
       ! Maybe in the OpenMP wrapper module?  Or a README?
       !
       ! The target_alloc_and_pair_ptrs function allows OpenMP to make use
       ! of memory from an umpire device memory pool.  However, any data
       ! mapped this way will not work with map to/from pragmas.  They will
       ! essentially be no-ops.
       !
       ! This function only supports standalone pointers, on its own.
       ! To work with nested data, a follow up map is required (see below).
       err_code = target_alloc_and_pair_ptrs(Quad% SetDataPtr(setID) % Psi, "Psi")
       err_code = target_alloc_and_pair_ptrs(Quad% SetDataPtr(setID)% Psi1, "Psi1")
       err_code = target_alloc_and_pair_ptrs(Quad% SetDataPtr(setID)% PsiB, "PsiB")
       err_code = target_alloc_and_pair_ptrs(Quad% SetDataPtr(setID)% Q, "Q")
       err_code = target_alloc_and_pair_ptrs(Quad% SetDataPtr(setID)% S, "S")
       err_code = target_alloc_and_pair_ptrs(Quad% SetDataPtr(setID)% cyclePsi, "cyclePsi")

       ! The maps below are essentially no-ops ( they do not actually allocate
       ! data on the gpu ).  However, when combined with the 'always' keyword
       ! they force the runtime to perform the pointer attachment between
       ! the derived type and its member pointers.
       ! Accesses such as 'Set%Psi' on the gpu will then work.
       !
       ! See OpenMP spec "2.19.7.1 map Clause" for details on this behavior.
       TOMP(target enter data map(always,alloc:Quad% SetDataPtr(setID)% Psi))
       TOMP(target enter data map(always,alloc:Quad% SetDataPtr(setID)% Psi1))
       TOMP(target enter data map(always,alloc:Quad% SetDataPtr(setID)% PsiB))
       TOMP(target enter data map(always,alloc:Quad% SetDataPtr(setID)% Q))
       TOMP(target enter data map(always,alloc:Quad% SetDataPtr(setID)% S))
       TOMP(target enter data map(always,alloc:Quad% SetDataPtr(setID)% cyclePsi))

       TOMP(target update to(Quad% SetDataPtr(setID)% Psi))
       TOMP(target update to(Quad% SetDataPtr(setID)% Psi1))
       TOMP(target update to(Quad% SetDataPtr(setID)% PsiB))
       TOMP(target update to(Quad% SetDataPtr(setID)% Q))
       TOMP(target update to(Quad% SetDataPtr(setID)% S))
       TOMP(target update to(Quad% SetDataPtr(setID)% cyclePsi))

     else

       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% Psi)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% Psi1)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% PsiB)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% Q)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% S)
       TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% cyclePsi)

     endif

   enddo

   return
   end subroutine initializeGPUMemory

