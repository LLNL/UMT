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

     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% Psi)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% Psi1)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% PsiB)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% Q)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% S)
     TOMP_TARGET_ENTER_DATA_MAP_TO(Quad% SetDataPtr(setID)% cyclePsi)

   enddo

   return
   end subroutine initializeGPUMemory

