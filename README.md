# Compiling and Running umt2016

## 0. Choose compilers. I use the following script to set up pgi and MPI paths
```
VERSION=17.4

export PGROUPD_LICENSE_FILE=/opt/pgi/license.pgi.$VERSION

MPI=/opt/ibm/spectrum_mpi
#MPI=/shared/comms/openmpi-2.0.1/gnu


export LD_LIBRARY_PATH=/opt/pgi/linuxpower/$VERSION/lib:/shared/lsf/10.1/linux3.10-glibc2.17-ppc64le/lib:/usr/local/cuda/lib64:/usr/local/lib

export LD_LIBRARY_PATH=$MPI/lib

export CPP="/opt/pgi/linuxpower/$VERSION/bin/pgcc -E"
export PGI=/opt/pgi
export PATH=/opt/pgi/linuxpower/$VERSION/bin:$MPI/bin:/shared/lsf/10.1/linux3.10-glibc2.17-ppc64le/etc:/shared/lsf/10.1/linux3.10-glibc2.17-ppc64le/bin::/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/cuda/bin:/opt/ibutils/bin:/usr/local/bin:/usr/local/sbin

export F90=/opt/pgi/linuxpower/$VERSION/bin/pgf90

export F77=/opt/pgi/linuxpower/$VERSION/bin/pgf77

export CXX=/opt/pgi/linuxpower/$VERSION/bin/pgc++

export FC=/opt/pgi/linuxpower/$VERSION/bin/pgfortran

export CC=/opt/pgi/linuxpower/$VERSION/bin/pgcc
```


## 1. from the top directory set the MPI compilers to use gcc:

```
source MPIdefs.sh
```

## 2. Then make in this top directory (make ibm timers first). If changes need to be made to compiler flags, they usually go in make.defs

```
cd ibmtimers
make
cd ..
make
```

## 3. Link by going to Teton directory and making SuOlsonTest target:

```
cd Teton
make SuOlsonTest
```

## 4. Run a test problem: 

### 4.1 See coral-runs for acceptence benchmark run scripts and output.

### 4.2 See runs directory for smaller, generic, and older examples.