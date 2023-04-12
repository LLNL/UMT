#!/bin/sh -x
# If you have a UMT tarball, untar it.  Otherwise, git clone it from github.
# git clone https://github.com/LLNL/UMT.git

# This script expects the following env variables to be set:
# CC = C compiler
# CXX = C++ compiler
# FC = Fortran compiler

# Examples of supported compilers
# Intel
# CC=icx
# CXX=icpx
# FC=ifx

# GNU
# CC=gcc
# CXX=g++
# FC=gfortran

# This script can either download tarballs for metis, conduit, hypre, and mfem, or it can copy them from a file system directory.  If you have the tarballs already present on the file system, then set this:
# TARBALLS=path/to/tarballs

# This script assumes that have a directory called 'UMT' under your current directory. If its in another location, modify this next line.
UMT_PATH=${PWD}/UMT

# Create workspace for building umt and required libraries conduit, metis, hypre, mfem.
INSTALL_PATH=${PWD}/umt_workspace/install
mkdir -p ${INSTALL_PATH}
echo Libraries will be installed to: ${INSTALL_PATH}

cd umt_workspace

# Download and build METIS
METIS=metis-5.1.0.tar.gz
if [ -f ${TARBALLS}/${METIS} ]; then
   cp ${TARBALLS}/${METIS} .
else
	wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/${METIS}
fi
tar xvzf ${METIS}
cd metis-5.1.0
cmake . -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DGKLIB_PATH=${PWD}/GKlib
gmake -j install
cd ..

# Download and build conduit
CONDUIT=conduit-v0.8.7-src-with-blt.tar.gz
if [ -f ${TARBALLS}/${CONDUIT} ]; then
   cp ${TARBALLS}/${CONDUIT} .
else
	wget https://github.com/LLNL/conduit/releases/download/v0.8.7/${CONDUIT}
fi
tar xvzf ${CONDUIT}
mkdir build_conduit
cd build_conduit
cmake ${PWD}/../conduit-v0.8.7/src -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC} -DMPI_CXX_COMPILER=mpicxx -DMPI_Fortran_COMPILER=mpifort -DBUILD_SHARED_LIBS=OFF -DENABLE_TESTS=OFF -DENABLE_EXAMPLES=OFF -DENABLE_DOCS=OFF -DENABLE_FORTRAN=ON -DENABLE_MPI=ON -DENABLE_PYTHON=OFF
gmake -j install
cd ..

# Download and build HYPRE
HYPRE=v2.24.0.tar.gz
if [ -f ${TARBALLS}/${HYPRE} ]; then
   cp ${TARBALLS}/${HYPRE} .
else
	wget https://github.com/hypre-space/hypre/archive/refs/tags/${HYPRE}
fi
tar xvzf ${HYPRE}
mkdir build_hypre
cd build_hypre
cmake ${PWD}/../hypre-2.24.0/src -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCMAKE_C_COMPILER=${CC}
gmake -j install
cd ..

# Download and build MFEM
MFEM=v4.4.zip
if [ -f ${TARBALLS}/${MFEM} ]; then
   cp ${TARBALLS}/${MFEM} .
else
	wget https://github.com/mfem/mfem/archive/refs/tags/${MFEM}
fi
unzip v4.4.zip
mkdir build_mfem
cd build_mfem
cmake ${PWD}/../mfem-4.4 -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DMPI_CXX_COMPILER=mpicxx -DMFEM_USE_MPI=TRUE -DMFEM_USE_CONDUIT=TRUE -DMFEM_USE_METIS_5=TRUE -DCMAKE_PREFIX_PATH=${INSTALL_PATH}}
gmake -j install
cd ..

# Build UMT
mkdir build_umt
cd build_umt

# If needed, additional cmake arguments for the UMT build can be provided on the script command line args.  This line will pick them up via the $1 on the next line.
cmake ${UMT_PATH}/src -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC} -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DMETIS_ROOT=${INSTALL_PATH} -DHYPRE_ROOT=${INSTALL_PATH} -DMFEM_ROOT=${INSTALL_PATH} -DCONDUIT_ROOT=${INSTALL_PATH} $1
gmake -j install
cd ..

# Test UMT on 10 cycles of an unstructured 3d mesh problem. Refine the mesh via -r and -R arguments.
srun -n1 ${INSTALL_PATH}/bin/makeUnstructuredBox
srun -n2 ${INSTALL_PATH}/bin/test_driver -i ./unstructBox3D.mesh -c 10 -r 1 -R 6
