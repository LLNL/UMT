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

   if (Size% ndim == 2) then

     UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% PsiM)
     TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% PsiM))

   endif

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% Psi)
   TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% Psi))

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% Psi1)
   TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% Psi1))

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% PsiB)
   TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% PsiB))

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% Q)
   TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% Q))

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% S)
   TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% S))

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% cyclePsi)
   TOMP(target exit data map(always,release:Quad% SetDataPtr(setID)% cyclePsi))

   return
   end subroutine finalizeGPUMemory
