#!/bin/bash

nodes=16
gpn=4  #gpu per node
ppn=20 # processes per node
let tpr=$ppn/$gpn
#echo "tasks per rs = $tpr"
let nrs=$nodes*$gpn
let nmpi=$nodes*$ppn
let cores_per_rank=40/$ppn
let threads_per_rank=2*$cores_per_rank
let cpu_per_rs=$tpr*$cores_per_rank

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
##BSUB -alloc_flags gpumps
#BSUB -alloc_flags smt2
###BSUB -P VEN201
#BSUB -q excl_6gpus
#BSUB -core_isolation 1
#---------------------------------------
ulimit -s 10240

export OMP_NUM_THREADS=$threads_per_rank
export CUDA_LAUNCH_BLOCKING=0

export OMP_STACKSIZE=64M
export PAMI_ENABLE_STRIPING=1

export LD_PRELOAD=/gpfs/wscgpfs01/walkup/mpitrace/spectrum_mpi/libhpmprof.so
export SAVE_LIST="1"

echo 'starting jsrun'

jsrun -X 1 --stdio_mode prepended -D CUDA_VISIBLE_DEVICES -E OMP_NUM_THREADS=$threads_per_rank \\
--rs_per_host ${gpn} \\
--gpu_per_rs 1  \\
--tasks_per_rs ${tpr}  \\
--cpu_per_rs ${cpu_per_rs} \\
--nrs ${nrs}  \\
-b packed:${cores_per_rank} \\
-d plane:${tpr} \\
../../Teton/SuOlsonTest $grid $groups $type $order $polar $azim


EOF
#---------------------------------------
bsub <batch.job
