# Compiling and Running umt2016

##0. Choose compilers. I use the following script to set up pgi and MPI paths
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


##1. from the top directory set the correct MPI compilers to use gcc:

```
source MPIdefs.sh
```

##2. Then make in this top directory. If changes need to be made to compiler flags, they usually go in make.defs

```
make
```

##3. Link by going to Teton directory and making SuOlsonTest target:

```
cd Teton
make SuOlsonTest
```

##4. Run a test problem:

###4.1 Create a run directory in the top level directory

```
cd ..
mkdir runtest
cd runtest
```

###4.1 you need a mesh descriptor file named 2x2x1_12.cmg Contents should be 

```
#processor block decomposition
sms(2,2,1)
#Always specify blocks in block base numbering
blk(on,0:1,0:1,0:0)

# tag boundary faces
tag("xMinFaces",face,(0:0,0:2,0:1))
tag("xMaxFaces",face,(2:2,0:2,0:1))
tag("yMinFaces",face,(0:2,0:0,0:1))
tag("yMaxFaces",face,(0:2,2:2,0:1))
tag("zMinFaces",face,(0:2,0:2,0:0))
tag("zMaxFaces",face,(0:2,0:2,1:1))

# define number of zones in each axis
numzones(12,12,12)

#Hex subdivisions
sub(10%,0:1, 0:1, 0:0,(7,0,0,0)) #7 hex
seed(10)
```

###4.2 Launch the job (uses 4 mpi ranks and 20 omp threads) on lsf with a script like this:

```
#!/bin/bash                                                                                                              
nodes=1
ppn=4
let nmpi=$nodes*$ppn
grid=2x2x1_12.cmg
order=16
groups=32
#groups=16                                                                                                               
type=2
polar=8
azim=4
#--------------------------------------                                                                                  
cat >batch.job <<EOF                                                                                                     
#BSUB -o %J.out                                                                                                          
#BSUB -e %J.err                                                                                                          
#BSUB -R "span[ptile=${ppn}]"                                                                                            
#BSUB -R "rusage[ngpus_shared=4]"                                                                                        
##BSUB -R "select[maxmem > 400000]"                                                                                       
#BSUB -n ${nmpi}                                                                                                         
#BSUB -x                                                                                                                 
#BSUB -q excl                                                                                                            
#BSUB -W 60                                                                                                             
#BSUB -env "all,LSB_START_JOB_MPS=N"                                                                                     
#---------------------------------------                                                                                 
export OMP_NUM_THREADS=20  
export OMP_WAIT_POLICY=active
export HPM_GROUP_LIST=10,21

mpirun --bind-to none -np $nmpi /home/walkup/bin/set_device_and_bind.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim
EOF                                                                                                                      
#---------------------------------------                                                                                 

bsub <batch.job
```