#!/bin/bash
# There are 21 (of 22) cores available to the application per socket
# Each core can go up to smt4 (smt2 gives best perf for umt).
nodes=4
gpus_per_socket=3 # number of gpus to use per socket
ranks_per_gpu=1 # ranks per gpu. If greater than 1, should use mps.
let ranks_per_socket=$gpus_per_socket*$ranks_per_gpu # needs to be evenly divisible by gpus_per_socket. 
let cores_per_rank=21/$ranks_per_socket # 21 avail cores divided into the ranks.
let nmpi=2*$ranks_per_socket*$nodes  # total number of mpi ranks
let cores_per_socket=$cores_per_rank*$ranks_per_socket # this is used cores per socket (not necessarily 21)
let num_sockets=$nodes*2 #nmpi/ranks_per_socket # total number of sockets
let threads_per_rank=2*$cores_per_rank

let res_sets=2*$ranks_per_socket*$nodes

echo "nodes = $nodes"
echo "gpus used per socket = $gpus_per_socket"
echo "ranks_per_socket = $ranks_per_socket"
echo "cores_per_rank = $cores_per_rank"
echo "used cores per socket = $cores_per_socket"
echo "threads per rank = $threads_per_rank"
#--------------------------------------
grid=4x2x3_38.cmg
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
###BSUB -alloc_flags gpumps
#BSUB -alloc_flags smt2
#BSUB -P VEN201
#BSUB -q batch
#BSUB -core_isolation 1
#BSUB -W 60
#---------------------------------------
ulimit -s 10240

export OMP_NUM_THREADS=$threads_per_rank
export CUDA_LAUNCH_BLOCKING=0

export OMP_STACKSIZE=64M
export PAMI_ENABLE_STRIPING=1

#export LD_PRELOAD=/gpfs/wscgpfs01/walkup/mpitrace/spectrum_mpi/libhpmprof.so
#export SAVE_LIST="1"

export PROFILE_RANK=-1
#export PROFILE_PATH="./nvp4_1.prof"

jsrun --stdio_mode=prepend -D CUDA_VISIBLE_DEVICES -E OMP_NUM_THREADS=$threads_per_rank --nrs ${res_sets}  --tasks_per_rs 1 --cpu_per_rs ${cores_per_rank}   --gpu_per_rs 1 --bind=proportional-packed:${cores_per_rank} -d plane:1 ./profile-helper.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim


EOF
#---------------------------------------
#bsub <batch.job
bsub  batch.job

