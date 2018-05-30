#!/bin/bash
# There are 21 (of 22) cores available to the application per socket
# Each core can go up to smt4
nodes=1
gpus_per_socket=3 # number of gpus to use per socket
ranks_per_gpu=1 # ranks per gpu. If greater than 1, should use mps.
let ranks_per_socket=$gpus_per_socket*$ranks_per_gpu # needs to be evenly divisible by gpus_per_socket. 
let cores_per_rank=21/$ranks_per_socket # 21 avail cores divided into the ranks.
let nmpi=2*$ranks_per_socket*$nodes  # total number of mpi ranks
let cores_per_socket=$cores_per_rank*$ranks_per_socket # this is used cores per socket (not necessarily 21)
let num_sockets=$nodes*2 #nmpi/ranks_per_socket # total number of sockets
let threads_per_rank=2*$cores_per_rank

echo "nodes = $nodes"
echo "gpus used per socket = $gpus_per_socket"
echo "ranks_per_socket = $ranks_per_socket"
echo "cores_per_rank = $cores_per_rank"
echo "used cores per socket = $cores_per_socket"
echo "threads per rank = $threads_per_rank"
#--------------------------------------
grid=3x2x1_12.cmg
#grid=1x2x3_38.cmg
#grid=1x2x3_40.cmg
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
###BSUB -P VEN201
#BSUB -q excl_6gpus
#BSUB -W 20
#---------------------------------------
ulimit -s 10240

export OMP_NUM_THREADS=$threads_per_rank
export CUDA_LAUNCH_BLOCKING=0

export OMP_STACKSIZE=64M
export PAMI_ENABLE_STRIPING=1
#export LD_PRELOAD=/ccs/home/walkup/logger/mpitrace/src/libmpitrace.so

echo 'starting jsrun with'
echo "nodes = $nodes"
echo "gpus used per socket = $gpus_per_socket"
echo "ranks_per_socket = $ranks_per_socket"
echo "cores_per_rank = $cores_per_rank"
echo "used cores per socket = $cores_per_socket"
echo "threads per rank = $threads_per_rank"

export RANKS_PER_SOCKET=$ranks_per_socket
export RANKS_PER_GPU=$ranks_per_gpu

# profiling stuff:
export PROFILE_RANK=-1  #rank where device-bind will run nvprof
#export PROFILE_PATH="/gpfs/alpinetds/scratch/dappelh/ven201/nvp216_3.prof"
echo "nvprof output at $PROFILE_PATH"


# -mxm

jsrun --smpiargs=" --mca btl_openib_warn_default_gid_prefix 0 --mca mpi_warn_on_fork 0" \
  --nrs ${num_sockets}  --tasks_per_rs ${ranks_per_socket} --cpu_per_rs ${cores_per_socket} \
  --gpu_per_rs ${gpus_per_socket} --bind=proportional-packed:${cores_per_rank} -d plane:${ranks_per_socket}  \
  ./device-bind.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

# jsrun --smpiargs="-mxm --mca btl_openib_warn_default_gid_prefix 0 --mca mpi_warn_on_fork 0" \
#   --nrs ${num_sockets}  --tasks_per_rs ${ranks_per_socket} --cpu_per_rs ${cores_per_socket} \
#   --gpu_per_rs ${gpus_per_socket} --bind=proportional-packed:${cores_per_rank} -d plane:${ranks_per_socket}  \
#   ../../sm60/Teton/SuOlsonTest $grid $groups $type $order $polar $azim


EOF
#---------------------------------------
bsub -core_isolation 1 <batch.job
#bsub  batch.job
