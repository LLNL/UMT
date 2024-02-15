#!/bin/bash -xe
# This script will compile a basic Release build of UMT.  Additional CMake options can be added to the command line args of this
# script, and they will be picked up and added to the UMT CMake command at the bottom of this script.
# For a list of supported CMake options, run 'ccmake /path/to/umt/src'.
# Do not copy this script out of the UMT repo directory, it assumes it is located next to the UMT source files in order to work.

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

# Set to 1 to optionally build UMT with UMPIRE support.
# For more information on UMPIRE CMake options, please see:
# https://umpire.readthedocs.io/en/develop/sphinx/advanced_configuration.html 
USE_UMPIRE=0

FFLAGS=-fallow-argument-mismatch
# Intel example
#CC=icx
#CXX=icpx
#FC=ifx

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

git clone --recurse-submodules  https://github.com/LLNL/conduit.git conduit -b v0.9.0
mkdir build_conduit
cd build_conduit
cmake ${PWD}/../conduit/src -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC} -DMPI_CXX_COMPILER=mpicxx -DMPI_Fortran_COMPILER=mpifort -DBUILD_SHARED_LIBS=OFF -DENABLE_TESTS=OFF -DENABLE_EXAMPLES=OFF -DENABLE_DOCS=OFF -DENABLE_FORTRAN=ON -DENABLE_MPI=ON -DENABLE_PYTHON=OFF
gmake -j install
cd ..

UMPIRE_CMAKE_ARGS=
UMPIRE_RUNLINE_ARGS=

if [ $USE_UMPIRE -eq 1 ]; then
  echo "Enabling UMPIRE support"
  # If building Umpire, enable it in the UMT CMake and provide the path to the installation.
  UMPIRE_CMAKE_ARGS="-DENABLE_UMPIRE=TRUE -DUMPIRE_ROOT=${INSTALL_PATH}"
  # If building Umpire, add the '-u 1' command line arg to the UMT test driver run line.
  # This will tell it to use an Umpire CPU memory pool, or if a GPU run, to use a CPU pinned memory pool.
  UMPIRE_RUNLINE_ARGS="-u 1"
  git clone --recurse-submodules https://github.com/LLNL/Umpire.git -b v2023.06.0

  mkdir build_umpire
  cd build_umpire
  cmake ${PWD}/../Umpire -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC} -DMPI_CXX_COMPILER=mpicxx -DMPI_Fortran_COMPILER=mpifort -DBUILD_SHARED_LIBS=OFF -DENABLE_TESTS=OFF -DENABLE_EXAMPLES=OFF -DENABLE_DOCS=OFF -DENABLE_FORTRAN=ON -DENABLE_MPI=ON
  gmake -j install
  cd ..
fi

# Run CMake on UMT, compile, and install.
cmake ${UMT_REPO_PATH}/src -DCMAKE_Fortran_FLAGS=${FFLAGS} -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC} -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCONDUIT_ROOT=${INSTALL_PATH} ${UMPIRE_CMAKE_ARGS} $1
gmake -j install
cd ..

# Run two smoke tests to verify executable, one on 2D 8x8 tiles mesh and one on 3D 4x4x4 tiles mesh.
srun -n 8 ${INSTALL_PATH}/bin/test_driver -c 10 -B local -d 8,8,0 --benchmark_problem 2 ${UMPIRE_RUNLINE_ARG}
srun -n 8 ${INSTALL_PATH}/bin/test_driver -c 10 -B local -d 4,4,4 --benchmark_problem 2 ${UMPIRE_RUNLINE_ARG}
