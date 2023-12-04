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

#define TOMP_TARGET_ENTER_DATA_MAP_TO(source) \
call target_alloc_and_pair_ptrs(source, "source");\
!$omp target enter data map(always,alloc:source);\
!$omp target update to(source)

#define TOMP_TARGET_EXIT_DATA_MAP_RELEASE(source) \
call target_free_and_unpair_ptrs(source, "source"); \
!$omp target exit data map(always, release:source)

#  define TOMP_TARGET_ENTER_DATA_MAP_ALLOC(source) !$omp target enter data map(alloc:source)
#  define TOMP_TARGET_ENTER_DATA_MAP_FROM(source) !$omp target exit data map(from:source)
#else
#  define TOMP(source)
#  define TOMPC(source)
#  define TOMP_TARGET_ENTER_DATA_MAP_TO(source)
#  define TOMP_TARGET_ENTER_DATA_MAP_ALLOC(source)
#  define TOMP_TARGET_ENTER_DATA_MAP_FROM(source)
#  define TOMP_TARGET_EXIT_DATA_MAP_RELEASE(source)
#endif
