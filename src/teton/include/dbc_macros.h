#if !defined(__TETON_DBC_MACROS_H__)
#define __TETON_DBC_MACROS_H__

#include "mpi.h"

//-----------------------------------------------------------------------------
// Code contract checks, similar to the Fortran versions.
//-----------------------------------------------------------------------------
// - TETON_ASSERT - enabled in debug builds or non-performance critical builds when extra checking desired.
//                 Will conditionally emit a message and shut down code if provided logical check fails.
//                 Define TETON_COMPILE_ASSERTS to enable.
//
// - TETON_VERIFY - always enabled, should used to ensure correct problem input or critical areas of code.
//                 Will conditionally emit a message and shut down code if provided logical check fails.
//
// - TETON_FATAL  - always enabled, use when code state is in an unrecoverable state.
//                 Will unconditionally emit an error messaged and shut down code.
//-----------------------------------------------------------------------------

#ifdef TETON_COMPILE_ASSERTS
#   define TETON_ASSERT_C(rank,condition,message) if (! (condition) ) { std::cerr << "Teton verification failed: " << message << " on rank: " << rank << ", in file: " << __FILE__ << ", at line: " << __LINE__ << std::endl; MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER); }
#else
#   define TETON_ASSERT_C(rank,condition,message)
#endif

// Per cppreference.com, any output sent to std::cerr stream objects is immediately flushed to the OS, so we
// shouldn't need a flush call on these.  Output sent to std::cerr will also cause std::cout to flush. -- black27
#define TETON_VERIFY_C(rank,condition,message) if (! (condition) ) { std::cerr << "Teton verification failed: " << message << ", on rank: " << rank << ", in file: " << __FILE__ << ", at line: " << __LINE__ << std::endl; MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER); }
#define TETON_FATAL_C(rank,message) std::cerr << "Teton fatal error: " << message << ", on rank: " << rank << ", in file: " << __FILE__ << ", at line: " << __LINE__ << std::endl; MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);

#endif
