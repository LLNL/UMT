#include "macros.h"
#include "omp_wrappers.h"
!***********************************************************************
!                        Last Update:  01/2018, PFN                    *
!                                                                      *
!   finalizeGPUMemory  - Releases GPU memory that at the conclusion    *
!                        of the  radiation step.                       *
!                                                                      *
!***********************************************************************

   subroutine finalizeGPUMemory(setID)

   use kind_mod
   use QuadratureList_mod
   use Size_mod
   use SetData_mod
   use OMPWrappers_mod
   use MemoryAllocator_mod

   implicit none

!  Arguments
   integer,       intent(in) :: setID

!  Locals
   integer :: err_code

!  Update Psi on the host before releasing it's memory on the device
   TOMP(target update from(Quad% SetDataPtr(setID)% Psi) )

!  Delete the arrays

   TOMP(target exit data map(release:Quad% SetDataPtr(setID)% AngleOrder))

   if (Size% ndim == 2) then

     TOMP(target exit data map(release:Quad% SetDataPtr(setID)% PsiM))

   endif

     if (Allocator%umpire_device_allocator_id >= 0) then
       err_code = target_free_and_unpair_ptrs(Quad% SetDataPtr(setID)% Psi)
       err_code = target_free_and_unpair_ptrs(Quad% SetDataPtr(setID)% Psi1)
       err_code = target_free_and_unpair_ptrs(Quad% SetDataPtr(setID)% PsiB)
       err_code = target_free_and_unpair_ptrs(Quad% SetDataPtr(setID)% Q)
       err_code = target_free_and_unpair_ptrs(Quad% SetDataPtr(setID)% S)
       err_code = target_free_and_unpair_ptrs(Quad% SetDataPtr(setID)% cyclePsi)

       ! I'm assuming these are necessary to force the runtime to perform the
       ! pointer unattachments between the derived types and their pointer
       ! members. --Aaron
       TOMP(target exit data map(always, release:Quad% SetDataPtr(setID)% Psi))
       TOMP(target exit data map(always, release:Quad% SetDataPtr(setID)% Psi1))
       TOMP(target exit data map(always, release:Quad% SetDataPtr(setID)% PsiB))
       TOMP(target exit data map(always, release:Quad% SetDataPtr(setID)% Q))
       TOMP(target exit data map(always, release:Quad% SetDataPtr(setID)% S))
       TOMP(target exit data map(always, release:Quad% SetDataPtr(setID)% cyclePsi))

     else

       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% Psi))
       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% Psi1))
       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% PsiB))
       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% Q))
       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% S))
       TOMP(target exit data map(release:Quad% SetDataPtr(setID)% cyclePsi))

     endif

   return
   end subroutine finalizeGPUMemory
