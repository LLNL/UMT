!-----------------------------------------------------------------------------
! Defines utility macros for use in Teton
! Note: Do not indent the '#' symbol, it causes FPP to choke.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! These macros provide a mechanism for enabling/disabling OpenMP target pragma lines
! and OpenMP target functionality, such as umpire host and device pool integration.
!
! If this macro logic gets any more complex, this might need to be reworked... -- black27 

! Enable basic openmp pragmas in code if openmp offloading enabled.
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
#  define TOMP(source) !$omp source
#  define TOMPC(source) !$omp& source
#else
#  define TOMP(source)
#  define TOMPC(source)
#endif

! Enable openmp data map and update pragmas, if openmp offloading enabled and not using unified cpu and gpu memory.
#if defined(TETON_ENABLE_OPENMP_OFFLOAD) && !defined(TETON_OPENMP_HAS_UNIFIED_MEMORY)
#  define TOMP_MAP(source) !$omp source
#  define TOMP_UPDATE(source) !$omp source
#else
#  define TOMP_MAP(source)
#  define TOMP_UPDATE(source)
#endif

! Enable Umpire integration for host and device pools, if openmp offloading enabled and not using unified cpu and gpu memory.
#if defined(TETON_ENABLE_UMPIRE) && defined(TETON_ENABLE_OPENMP_OFFLOAD) && !defined(TETON_OPENMP_HAS_UNIFIED_MEMORY)
#  define UMPIRE_DEVICE_POOL_ALLOC(source) call target_alloc_and_pair_ptrs(source)
#  define UMPIRE_DEVICE_POOL_FREE(source) call target_free_and_unpair_ptrs(source)
#else
#  define UMPIRE_DEVICE_POOL_ALLOC(source)
#  define UMPIRE_DEVICE_POOL_FREE(source)
#endif
