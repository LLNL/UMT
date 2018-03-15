#!/bin/bash
nodes=4
ppn=4
let nmpi=$nodes*$ppn
grid=4x2x2_20.cmg
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
#BSUB -n ${nmpi}
#BSUB -x
##BSUB -q excl_short
#BSUB -q excl
##BSUB -W 30
#BSUB -env "all,LSB_START_JOB_MPS=N"
#---------------------------------------
export OMP_NUM_THREADS=20
export OMP_WAIT_POLICY=active
export HPM_GROUP_LIST=10,21


mpirun --bind-to none -np $nmpi /home/walkup/bin/set_device_and_bind.sh nvprof -s -f -o 1node.%q{OMPI_COMM_WORLD_RANK}.nvprof ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim



#mpirun --mxm --bind-to none -np $nmpi ./bind.sh /home/walkup/bin/set_device.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim



#mpirun --mxm --bind-to none -np $nmpi ./bind.sh ~/bin/mps_helper.sh nvprof -s -f -o 1node.%q{OMPI_COMM_WORLD_RANK}.nvprof ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

#mpirun --mxm --bind-to none -np $nmpi ./bind.sh ~/bin/mps_helper.sh nvprof -s -f -o 1node.%q{OMPI_COMM_WORLD_RANK}.nvprof ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

#mpirun --bind-to socket -np $nmpi nvprof -s -o -f 1node.%q{OMPI_COMM_WORLD_RANK}.nvprof ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

#/gpfs/ess2fs0/walkup/openmpi-1.8.8/bin/mpirun -np $nmpi ./bind.sh cuda-memcheck --language fortran ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim
EOF
#---------------------------------------
#bsub -U apex_fri <batch.job
bsub <batch.job
