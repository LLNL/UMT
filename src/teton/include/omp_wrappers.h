!-----------------------------------------------------------------------------
! Defines utility macros for use in Teton
! Note: Do not indent the '#' symbol, it causes FPP to choke.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! These macros provide a mechanism for enabling/disabling OpenMP target pragma lines.
! Some compilers have trouble compiling newer OpenMP target offload pragmas.
!
! TETON_ENABLE_OPENMP_OFFLOAD - enables any lines annotated with the TOMP macros, should be set by the build system.
!
! TOMP -  Macro to enable/disable OMP target offload pragma line.
! TOMPC -  Macro to enable/disable continued OMP pragma line.
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
#  define TOMP(source) !$omp source
#  define TOMPC(source) !$omp& source
#  if defined(TETON_ENABLE_UMPIRE)
#    define UMPIRE_DEVICE_POOL_ALLOC(source) call target_alloc_and_pair_ptrs(source)
#    define UMPIRE_DEVICE_POOL_FREE(source) call target_free_and_unpair_ptrs(source)
#  else
#    define UMPIRE_DEVICE_POOL_ALLOC(source)
#    define UMPIRE_DEVICE_POOL_FREE(source)
#  endif

#else
#  define TOMP(source)
#  define TOMPC(source)
#  define UMPIRE_DEVICE_POOL_ALLOC(source)
#  define UMPIRE_DEVICE_POOL_FREE(source)
#endif
