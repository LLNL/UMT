#include "macros.h"
#include "omp_wrappers.h"
!******************************************************************************
!   OMPWrappers -  A set of subroutines that wrap OpenMP operations for       *
!                  mapping data structures to/from a target device, or        *
!                  updating data between the host and device.                 *
!******************************************************************************

module OMPWrappers_mod

#if defined(TETON_ENABLE_OPENMP)
   use omp_lib
#endif

!******************************************************************************
! If using the UMPIRE library with OpenMP offload kernels, a family of
! procedures is provided below to allow the code to make use of UMPIRE for
! device allocations.
!******************************************************************************
#if defined(TETON_ENABLE_OPENMP_OFFLOAD) && defined(TETON_ENABLE_UMPIRE)

   use iso_c_binding, only : c_double, c_int, c_bool, c_size_t, c_ptr, c_null_ptr, c_loc, c_associated
   use AngleSet_mod, only : AngleSet, HypPlane, BdyExit
   use Geometry_mod, only : Geometry
   use GreyAcceleration_mod, only : GreyAcceleration
   use GroupSet_mod, only : GroupSet
   use QuadratureList_mod, only : QuadratureList
   use SetData_mod, only : SetData
   use MemoryAllocator_mod

   implicit none

   interface target_alloc_and_pair_ptrs
      module procedure target_alloc_and_pair_ptrs_real_c_double_1, &
                       target_alloc_and_pair_ptrs_real_c_double_2, &
                       target_alloc_and_pair_ptrs_real_c_double_3, &
                       target_alloc_and_pair_ptrs_integer_c_int_1, &
                       target_alloc_and_pair_ptrs_integer_c_int_2, &
                       target_alloc_and_pair_ptrs_integer_c_int_3, &
                       target_alloc_and_pair_ptrs_logical_c_bool_1, &
                       target_alloc_and_pair_ptrs_type_AngleSet_1, &
                       target_alloc_and_pair_ptrs_type_BdyExit_1, &
                       target_alloc_and_pair_ptrs_type_GroupSet_1, &
                       target_alloc_and_pair_ptrs_type_HypPlane_1, &
                       target_alloc_and_pair_ptrs_type_SetData_1
   end interface

   interface target_free_and_unpair_ptrs
      module procedure target_free_and_unpair_ptrs_real_c_double_1, &
                       target_free_and_unpair_ptrs_real_c_double_2, &
                       target_free_and_unpair_ptrs_real_c_double_3, &
                       target_free_and_unpair_ptrs_integer_c_int_1, &
                       target_free_and_unpair_ptrs_integer_c_int_2, &
                       target_free_and_unpair_ptrs_integer_c_int_3, &
                       target_free_and_unpair_ptrs_logical_c_bool_1, &
                       target_free_and_unpair_ptrs_type_AngleSet_1, &
                       target_free_and_unpair_ptrs_type_BdyExit_1, &
                       target_free_and_unpair_ptrs_type_GroupSet_1, &
                       target_free_and_unpair_ptrs_type_HypPlane_1, &
                       target_free_and_unpair_ptrs_type_SetData_1
   end interface

#if !defined(TETON_OPENMP_HAS_FORTRAN_INTERFACE)

   public :: omp_target_alloc, omp_target_free, omp_target_associate_ptr, omp_target_disassociate_ptr, omp_target_is_present

   interface

      type(C_PTR) function omp_target_alloc( num_bytes, device_num ) bind (c)
        use iso_c_binding
        implicit none

        integer(C_SIZE_T), value :: num_bytes
        integer(C_INT), value :: device_num
      end function omp_target_alloc

      subroutine omp_target_free( h_ptr, device_num ) bind (c)
        use iso_c_binding
        implicit none

        type(C_PTR), value :: h_ptr
        integer(C_INT), value :: device_num
      end subroutine omp_target_free

      integer (C_INT) function omp_target_associate_ptr( h_ptr, d_ptr, num_bytes, offset, device_num) bind (c)
        use iso_c_binding
        implicit none

        type(C_PTR), value :: h_ptr, d_ptr
        integer(C_SIZE_T), value :: num_bytes, offset
        integer(C_INT), value :: device_num
      end function omp_target_associate_ptr

      integer (C_INT) function omp_target_disassociate_ptr( h_ptr, device_num) bind(c)
        use iso_c_binding
        implicit none

        type(C_PTR), value :: h_ptr
        integer(C_INT), value :: device_num
      end function omp_target_disassociate_ptr

      integer (C_INT) function omp_target_is_present( h_ptr, device_num) bind(c)
        use iso_c_binding
        implicit none

        type(C_PTR), value :: h_ptr
        integer(C_INT), value :: device_num
      end function omp_target_is_present

   end interface

#  endif

contains

   !--------------------------------------------------------------------------
   ! overloaded family of functions for array of supported types
   ! The #defines do not need to be unset here, that occurs inside the
   ! templates.
   !--------------------------------------------------------------------------
 
   ! ARRAYS OF DOUBLE PRECISION
#define FTM_TYPE real
#define FTM_KIND c_double
#define FTM_RANK 1
#include "OMPWrappers.F90.templates"

#define FTM_TYPE real
#define FTM_KIND c_double
#define FTM_RANK 2
#include "OMPWrappers.F90.templates"

#define FTM_TYPE real
#define FTM_KIND c_double
#define FTM_RANK 3
#include "OMPWrappers.F90.templates"

   ! ARRAYS OF INTEGERS
#define FTM_TYPE integer
#define FTM_KIND c_int
#define FTM_RANK 1
#include "OMPWrappers.F90.templates"

#define FTM_TYPE integer
#define FTM_KIND c_int
#define FTM_RANK 2
#include "OMPWrappers.F90.templates"

#define FTM_TYPE integer
#define FTM_KIND c_int
#define FTM_RANK 3
#include "OMPWrappers.F90.templates"

   ! ARRAYS OF LOGICALS
#define FTM_TYPE logical
#define FTM_KIND c_bool
#define FTM_RANK 1
#include "OMPWrappers.F90.templates"

   ! ARRAYS OF DERIVED TYPES
#define FTM_TYPE type
#define FTM_KIND AngleSet
#define FTM_RANK 1
#include "OMPWrappers.F90.templates"

#define FTM_TYPE type
#define FTM_KIND GroupSet
#define FTM_RANK 1
#include "OMPWrappers.F90.templates"

#define FTM_TYPE type
#define FTM_KIND SetData
#define FTM_RANK 1
#include "OMPWrappers.F90.templates"

#define FTM_TYPE type
#define FTM_KIND HypPlane
#define FTM_RANK 1
#include "OMPWrappers.F90.templates"

#define FTM_TYPE type
#define FTM_KIND BdyExit
#define FTM_RANK 1
#include "OMPWrappers.F90.templates"

#endif

end module OMPWrappers_mod
