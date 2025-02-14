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
   use Options_mod
   use MemoryAllocator_mod

   implicit none

!  Arguments
   integer,       intent(in) :: setID

!  Local
   integer :: dom
   integer :: nHyperDomains 
   integer :: sweepVersion

   nHyperDomains = getNumberOfHyperDomains(Quad,1)
   sweepVersion  = Options% getSweepVersion()

!  Update Psi on the host before releasing it's memory on the device
   TOMP_UPDATE(target update from(Quad% SetDataPtr(setID)% Psi) )

!  Delete the arrays

   if (Size% ndim == 2) then

     UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% PsiM)
     TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% PsiM))

   endif

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% Psi)
   TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% Psi))

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% Psi1)
   TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% Psi1))

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% PsiB)
   TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% PsiB))

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% cyclePsi)
   TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% cyclePsi))

   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% PsiInt)
   TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% PsiInt))

! These are only used in the zone sweep.
   if ( sweepVersion == 1 ) then

     do dom=1,nHyperDomains
       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% SweepPtr(dom)% Q)
       TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% SweepPtr(dom)% Q))

       UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% SweepPtr(dom)% S)
       TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% SweepPtr(dom)% S))
     enddo

   endif

!   Both Cray and XLF throw an error on the following line  PFN 02/07/2024
!   UMPIRE_DEVICE_POOL_FREE(Quad% SetDataPtr(setID)% SweepPtr)
   TOMP_MAP(target exit data map(always,release:Quad% SetDataPtr(setID)% SweepPtr))


   return
   end subroutine finalizeGPUMemory
