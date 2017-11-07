#!/bin/bash
nodes=1
ppn=1
let nmpi=$nodes*$ppn
grid=1x1x1_32.cmg
order=16
#groups=16
groups=32
type=2
polar=8
azim=4
#polar=8
#azim=4
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=${ppn}]"
#BSUB -R "rusage[ngpus_shared=4]"
##BSUB -R "select[maxmem > 400000]"
#BSUB -n ${nmpi}
#BSUB -x
##BSUB -q shared
#BSUB -q excl
#BSUB -W 20
#---------------------------------------
export OMP_NUM_THREADS=20
export OMP_WAIT_POLICY=active
export HPM_GROUP_LIST=10,21
#/gpfs/ess2fs0/walkup/openmpi-1.8.8/bin/mpirun --bind-to none -np $nmpi ./bind.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim
#/gpfs/ess2fs0/walkup/openmpi-1.8.8/bin/mpirun -np $nmpi ./bind.sh nvprof --print-gpu-trace ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim
#/gpfs/ess2fs0/walkup/openmpi-1.8.8/bin/mpirun -np $nmpi ./bind.sh nvprof -s ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim
#/gpfs/ess2fs0/walkup/openmpi-1.8.8/bin/mpirun --bind-to none -np $nmpi ./bind.sh nvprof -s ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

#mpirun --bind-to none -np $nmpi ./bind.sh nvprof -s ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

#mpirun --bind-to none -np $nmpi ./bind.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

mpirun --bind-to none -np $nmpi ./bind.sh  /home/walkup/bin/set_device.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

#mpirun --bind-to none -np $nmpi ./bind.sh nvprof --print-gpu-trace --metrics achieved_occupancy ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim
#mpirun --bind-to none -np $nmpi ./bind.sh nvprof --print-gpu-trace ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim


#mpirun --bind-to none -np $nmpi ./bind.sh cuda-memcheck --language fortran ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim
EOF
#---------------------------------------
#bsub -U apex_fri <batch.job
#bsub -U apex_tues <batch.job
bsub <batch.job
