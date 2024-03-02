System libraries
-----------------------
UMT Requires MPI to be installed on your system to provide a C++ and Fortran compiler.

To build UMT with LLVM Clang/Flang, you must also use them to build MPI so that
Fortran mod files are compatible.  For example, you might build OpenMPI as
follows:

```
$ CC=clang CXX=clang++ FC=flang-new \
  CFLAGS=-O3 CXXFLAGS=-O3 FCFLAGS=-O3 \
  ../configure --prefix=$PWD/../install-flang \
    --without-knem --without-ucx --enable-pretty-print-stacktrace \
    --enable-orterun-prefix-by-default --enable-mpi1-compatibility
$ make -j
$ make install
```

Reference: <https://github.com/llvm/llvm-project/issues/56348>

Third party libraries
-----------------------
The code depends on several libraries.

Required libraries are:
- Conduit, a io interchange library
  https://github.com/LLNL/conduit

Optional libraries:
- MFEM, a finite element methods library
  https://github.com/mfem/mfem

MFEM requires the additional libraries METIS and HYPRE.  See the MFEM github website for more information.
