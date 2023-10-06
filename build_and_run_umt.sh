#!/bin/sh -x
# If you have a UMT tarball, untar it.  Otherwise, git clone it from github.
# git clone https://github.com/LLNL/UMT.git

# This script expects the following variables to be set.
# CC = <path to C compiler>
# CXX = <path to C++ compiler>
# FC = <path to Fortran compiler>

# Default to GNU
CC=gcc
CXX=g++
FC=gfortran

FFLAGS=-fallow-argument-mismatch

# Intel example
# CC=icx
# CXX=icpx
# FC=ifx


# This script will compile a basic Release build of UMT.  Additional CMake options can be added to the command line args of this script, and they will be picked up and added to the UMT CMake command at the bottom of this script.
# For a list of supported CMake options, run 'ccmake /path/to/umt/src'.
# Do not copy this script out of the UMT repo directory, it assumes it is located next to the UMT source files in order to work.

# Get directory this script is located in.  This is assumed to be the UMT repo location.
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
UMT_REPO_PATH="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

# Create workspace for building umt and required libraries conduit, metis, hypre, mfem.
INSTALL_PATH=${PWD}/umt_workspace/install
mkdir -p ${INSTALL_PATH}
echo Libraries will be installed to: ${INSTALL_PATH}

cd umt_workspace

git clone --recurse-submodules  https://github.com/LLNL/conduit.git conduit
mkdir build_conduit
cd build_conduit
cmake ${PWD}/../conduit/src -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC} -DMPI_CXX_COMPILER=mpicxx -DMPI_Fortran_COMPILER=mpifort -DBUILD_SHARED_LIBS=OFF -DENABLE_TESTS=OFF -DENABLE_EXAMPLES=OFF -DENABLE_DOCS=OFF -DENABLE_FORTRAN=ON -DENABLE_MPI=ON -DENABLE_PYTHON=OFF
gmake -j install
cd ..

# Run CMake on UMT, compile, and install.
cmake ${UMT_REPO_PATH}/src -DCMAKE_Fortran_FLAGS=${FFLAGS} -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC} -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCONDUIT_ROOT=${INSTALL_PATH} $1
gmake -j install
cd ..

srun -n 8 ${INSTALL_PATH}/bin/test_driver -c 10 -B -d 8,8,0 --benchmark_problem 2
srun -n 8 ${INSTALL_PATH}/bin/test_driver -c 10 -B -d 4,4,4 --benchmark_problem 2

# Test UMT on SSP1 unstructured 3d mesh problem on two mpi ranks. Refine the mesh via -r and -R arguments.
