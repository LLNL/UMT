#!/bin/sh -x
# This script will compile a basic RelWithDebInfo build of UMT.  Additional CMake options can be added to the command line args of this script, and they will be picked up and added to the UMT CMake command at the bottom of this script.
# For a list of supported CMake options, run 'ccmake /path/to/umt/src'.
# Do not copy this script out of the UMT repo directory, it assumes it is located next to the UMT source files in order to work.

# If you have a UMT tarball, untar it.  Otherwise, git clone it from github.
# git clone https://github.com/LLNL/UMT.git

# This script expects the following variables to be set.
# CC = <path to C compiler>
# CXX = <path to C++ compiler>
# FC = <path to Fortran compiler>

# Alternatively, you can set the compilers to point to MPI compiler wrappers, if CMake has trouble finding your MPI installation location.

# Example using GNU
#CC=gcc
#CXX=g++
#FC=gfortran

# Example using Intel
# CC=icx
# CXX=icpx
# FC=ifx

# Example using MPI compiler wrappers
CC=mpicc
CXX=mpicxx
FC=mpif90

#If using CUDA, set this to point to your CUDA installation
#export CUDA_TOOLKIT_ROOT_DIR=/usr/tce/packages/cuda/cuda-11.8.0

# Desired installation path
INSTALL_PATH=${PWD}/umt_workspace/install

# Set common CMAKE args for Conduit and UMT.
CMAKE_COMMON_ARGS="-DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC}"

# Get directory this script is located in.  This is assumed to be the UMT repo location.
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
UMT_REPO_PATH="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

# Create workspace for building umt and required library conduit
mkdir -p ${INSTALL_PATH}
echo Libraries will be installed to: ${INSTALL_PATH}

cd umt_workspace

# Git latest conduit source from github
git clone --recurse-submodules  https://github.com/LLNL/conduit.git conduit
mkdir build_conduit
cd build_conduit
# Run CMake on Conduit, compile, and install.
cmake ${PWD}/../conduit/src ${CMAKE_COMMON_ARGS} -DBUILD_SHARED_LIBS=OFF -DENABLE_TESTS=OFF -DENABLE_EXAMPLES=OFF -DENABLE_DOCS=OFF -DENABLE_FORTRAN=ON -DENABLE_MPI=ON -DENABLE_PYTHON=OFF -DENABLE_UTILS=OFF -DENABLE_RELAY_WEBSERVER=OFF
gmake -j install
cd ..

mkdir build_umt
cd build_umt
# Run CMake on UMT, compile, and install.
cmake ${UMT_REPO_PATH}/src ${CMAKE_COMMON_ARGS} -DCONDUIT_ROOT=${INSTALL_PATH}
gmake -j install
cd ..

# Test on 2D and 3D problem
srun -n 8 ${INSTALL_PATH}/bin/test_driver -c 10 -B -d 8,8,0 --benchmark_problem 2
srun -n 8 ${INSTALL_PATH}/bin/test_driver -c 10 -B -d 4,4,4 --benchmark_problem 2
