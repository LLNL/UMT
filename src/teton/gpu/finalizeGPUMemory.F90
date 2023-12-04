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

     TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% PsiM)

   endif

   TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% Psi)
   TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% Psi1)
   TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% PsiB)
   TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% Q)
   TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% S)
   TOMP_TARGET_EXIT_DATA_MAP_RELEASE(Quad% SetDataPtr(setID)% cyclePsi)

   return
   end subroutine finalizeGPUMemory
