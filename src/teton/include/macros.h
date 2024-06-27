#if !defined(__TETON_MACROS_H__)
#define __TETON_MACROS_H__

!-----------------------------------------------------------------------------
! Defines utility macros for use in Teton
! Note: Do not indent the '#' symbol, it causes FPP to choke.
!-----------------------------------------------------------------------------
!#define TETON_CHECK_OUT_OF_BOUNDS_ARRAY_ACCESSES

!-----------------------------------------------------------------------------
! Suppress unused variable warnings
!-----------------------------------------------------------------------------
#define SUPPRESS_UNUSED_VAR_WARNING(x) if (.FALSE.) print *, x

!-----------------------------------------------------------------------------
! Code contract checks.
!-----------------------------------------------------------------------------
! - TETON_ASSERT - enabled in debug builds or non-performance critical builds when extra checking desired.
!                 Will conditionally emit a message and shut down code if provided logical check fails.
!                 Define TETON_COMPILE_ASSERTS to enable.
!
! - TETON_VERIFY - always enabled, should used to ensure correct problem input or critical areas of code.
!                 Will conditionally emit a message and shut down code if provided logical check fails.
!
! - TETON_FATAL  - always enabled, use when code state is in an unrecoverable state.
!                 Will unconditionally emit an error messaged and shut down code.
! - AB
#ifdef TETON_COMPILE_ASSERTS
#   define TETON_ASSERT(bool,s) call f90assert(bool,__FILE__,__LINE__,s)
#else
#   define TETON_ASSERT(bool,s)
#endif

#define TETON_VERIFY(bool,s) call f90verify(bool,__FILE__,__LINE__,s)
! The old f90fatal does only accepts a message, not file and line number info.
#define TETON_FATAL(s) call f90fatal2(__FILE__,__LINE__,s)

!-----------------------------------------------------------------------------
! Convenience macro to verify an array access won't be out-of-bounds
!-----------------------------------------------------------------------------
#if defined(TETON_CHECK_OUT_OF_BOUNDS_ARRAY_ACCESSES)
#    define TETON_CHECK_BOUNDS1(ptr, i) call f90verify(is_legal_access(ptr, i ), __FILE__, __LINE__, "Out of bounds array access.")
#    define TETON_CHECK_BOUNDS2(ptr, i, j) call f90verify(is_legal_access(ptr, i, j ), __FILE__, __LINE__, "Out of bounds array access.")
#    define TETON_CHECK_BOUNDS3(ptr, i, j, k) call f90verify(is_legal_access(ptr, i, j, k ), __FILE__, __LINE__, "Out of bounds array access.")
#else
#    define TETON_CHECK_BOUNDS1(ptr, i)
#    define TETON_CHECK_BOUNDS2(ptr, i, j)
#    define TETON_CHECK_BOUNDS3(ptr, i, j, k)
#endif

!-----------------------------------------------------------------------------
! Macro to start/stop a profiler range.
!
! Requires the Caliper library.  For NVIDIA GPU ranges, enable the NVPROF backend in Caliper.
!-----------------------------------------------------------------------------
#if defined(TETON_ENABLE_CALIPER)
#   define START_RANGE(name) call cali_begin_region(name)
#   define END_RANGE(name) call cali_end_region(name)
#else
#   define START_RANGE(name)
#   define END_RANGE(name)
#endif

! In implementations that support a preprocessor, the _OPENMP macro name is
! defined to have the decimal value yyyymm where yyyy and mm are the year and
! month designations of the version of the OpenMP API that the implementation supports.

! For quick reference, its:
! 200505 2.5
! 200805 3.0
! 201107 3.1
! 201307 4.0
! 201511 4.5
! 201811 5.0
! 202011 5.1
!
! This can be used to control what code is compiled across our compilers, but if
! we have a compiler that supports a newer feature, but is reporting an older
! implementation we have to rely on using a separate preprocessor define.

#if (_OPENMP > 201511)
# define TETON_OPENMP_HAS_FORTRAN_INTERFACE
# define TETON_OPENMP_HAS_USE_DEVICE_ADDR
#endif

#ifdef TETON_ENABLE_OPENACC
# define ATOMIC_UPDATE !$acc atomic update
# define ATOMIC_END !$acc end atomic
#elif defined(TETON_ENABLE_OPENMP_OFFLOAD)
# define ATOMIC_UPDATE !$omp atomic update
# define ATOMIC_END !$omp end atomic
#else
# define ATOMIC_UPDATE
# define ATOMIC_END
#endif

#endif
