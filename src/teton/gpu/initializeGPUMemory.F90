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
       UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID) % PsiM)
       TOMP(target enter data map(always,to:Quad% SetDataPtr(setID) % PsiM))

     endif

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% Psi)
     TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% Psi))

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% Psi1)
     TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% Psi1))

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% PsiB)
     TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% PsiB))

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% Q)
     TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% Q))

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% S)
     TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% S))

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% cyclePsi)
     TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% cyclePsi))

   enddo

   return
   end subroutine initializeGPUMemory

