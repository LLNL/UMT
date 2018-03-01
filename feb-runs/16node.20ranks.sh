#!/bin/bash
nodes=16
ppn=20
let nmpi=$nodes*$ppn
#grid=4x4x8_18.cmg
#grid=5x8x8_10.cmg
grid=5x8x8_20.cmg
order=16
groups=16
type=2
polar=8
azim=4
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -nnodes ${nodes}
#BSUB -q excl_6gpus
##BSUB -q excl_4gpus
#BSUB -W 30
#---------------------------------------

ulimit -s 10240
export BIND_THREADS=yes
export USE_GOMP=yes
export USE_MPS=yes
export RANKS_PER_NODE=${ppn}

export CUDA_LAUNCH_BLOCKING=0

echo 'starting jsrun'

# Use this one, even for 4 gpu runs (4/6 gpus used)

/opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun --rs_per_host 1 --tasks_per_rs ${ppn} --cpu_per_rs 44 --gpu_per_rs 6 --nrs ${nodes} -d plane:${ppn} ./helper_4gpu.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim


# Use this command with 4 GPU que and have to modify helper to have gpu 0,1,2,3 

#/opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun --rs_per_host 1 --tasks_per_rs ${ppn} --cpu_per_rs 44 --gpu_per_rs 4 --nrs ${nodes} -d plane:${ppn} ./helper_4gpu.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim


EOF
#---------------------------------------
#bsub -core_isolation y <batch.job
bsub <batch.job
