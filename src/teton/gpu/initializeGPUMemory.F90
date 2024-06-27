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
   use Options_mod
   use MemoryAllocator_mod

   implicit none

!  Local 

   integer  :: setID
   integer  :: dom
   integer  :: nSets
   integer  :: nHyperDomains 
   integer  :: sweepVersion

!  Constants

   nSets         = getNumberOfSets(Quad)
   nHyperDomains = getNumberOfHyperDomains(Quad,1)
   sweepVersion  = Options% getSweepVersion()

   do setID=1,nSets

     if (Size% ndim == 2) then
       UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID) % PsiM)
       TOMP_MAP(target enter data map(always,to:Quad% SetDataPtr(setID) % PsiM))
     endif

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% Psi)
     TOMP_MAP(target enter data map(always,to:Quad% SetDataPtr(setID)% Psi))

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% Psi1)
     TOMP_MAP(target enter data map(always,to:Quad% SetDataPtr(setID)% Psi1))

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% PsiB)
     TOMP_MAP(target enter data map(always,to:Quad% SetDataPtr(setID)% PsiB))

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% cyclePsi)
     TOMP_MAP(target enter data map(always,to:Quad% SetDataPtr(setID)% cyclePsi))

     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% PsiInt)
     TOMP_MAP(target enter data map(always,to:Quad% SetDataPtr(setID)% PsiInt))

!    Both Cray and XLF throw an error on the following line  PFN 02/07/2024
!     UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% SweepPtr)
     TOMP_MAP(target enter data map(always,to:Quad% SetDataPtr(setID)% SweepPtr))

     if ( sweepVersion == 0 ) then

       do dom=1,nHyperDomains
         UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% SweepPtr(dom)% Q)
         TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% SweepPtr(dom)% Q))

         UMPIRE_DEVICE_POOL_ALLOC(Quad% SetDataPtr(setID)% SweepPtr(dom)% S)
         TOMP(target enter data map(always,to:Quad% SetDataPtr(setID)% SweepPtr(dom)% S))
       enddo

     endif

   enddo

   return
   end subroutine initializeGPUMemory

