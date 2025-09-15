#include "macros.h"

module OMPUtilities_mod
   use, intrinsic :: iso_c_binding
   implicit none

   contains

#if defined(TETON_ENABLE_OPENMP)
   subroutine print_thread_bindings()
      interface 
         subroutine teton_print_thread_bindings_c() bind(c)
         end subroutine teton_print_thread_bindings_c
      end interface

      call teton_print_thread_bindings_c()
   end subroutine print_thread_bindings
#endif

#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
   function get_gpu_processor_count() result(num_procs)
      integer(c_int) num_procs
      interface 
         function teton_get_gpu_processor_count_c() bind(c) result(res)
            use iso_c_binding, only : c_int
            integer(c_int) res
         end function
      end interface

      num_procs = teton_get_gpu_processor_count_c()

   end function get_gpu_processor_count
#endif


end module OMPUtilities_mod
