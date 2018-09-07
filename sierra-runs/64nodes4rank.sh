#!/bin/bash
# On Sierra there are 20 (of 22) cores available to the application per socket
# Each core can go up to smt4
nodes=64
gpus_per_socket=2 # number of gpus to use per socket
ranks_per_gpu=1 # ranks per gpu. If greater than 1, should use mps.
let ranks_per_socket=$gpus_per_socket*$ranks_per_gpu # needs to be evenly divisible by gpus_per_socket. 
let cores_per_rank=20/$ranks_per_socket # 21 avail cores divided into the ranks.
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
grid=8x8x4_32.cmg
order=16
groups=16
type=2
polar=8
azim=4
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -nnodes $nodes
##BSUB -csm y
##BSUB -R "1*{select[LN]} + 1344*{select[CN&&(hname!=c699c250)&&(type==any)]span[ptile=42]}"

##BSUB -alloc_flags gpumps
##BSUB -alloc_flags smt2
#BSUB -G guests
#BSUB -q pbatch
#BSUB -core_isolation 2
#BSUB -W 20
#---------------------------------------
ulimit -s 10240

export OMP_NUM_THREADS=$threads_per_rank
export CUDA_LAUNCH_BLOCKING=0

export OMP_STACKSIZE=64M
export PAMI_ENABLE_STRIPING=1

export OMPI_LD_PRELOAD_POSTPEND=${OPAL_PREFIX}/lib/libmpitrace.so

#export LD_PRELOAD=/ccs/home/walkup/logger/mpitrace/src/libmpitrace.so

echo 'starting jsrun with'
echo "nodes = $nodes"
echo "gpus used per socket = $gpus_per_socket"
echo "ranks_per_socket = $ranks_per_socket"
echo "cores_per_rank = $cores_per_rank"
echo "used cores per socket = $cores_per_socket"
echo "threads per rank = $threads_per_rank"

export RANKS_PER_SOCKET=2
export RANKS_PER_GPU=1

# profiling stuff:
export PROFILE_RANK=-1  #rank where device-bind will run nvprof
#export PROFILE_PATH="/gpfs/alpinetds/scratch/dappelh/ven201/nvp216_3.prof"
echo "nvprof output at "


# -mxm

tstamp=`date +%m_%d_%H_%M_%S`; jsrun -X 1 --stdio_mode prepended --progress ${tstamp}.progress -D CUDA_VISIBLE_DEVICES -E OMP_NUM_THREADS=$threads_per_rank --nrs ${res_sets}   --tasks_per_rs 1 --cpu_per_rs ${cores_per_rank}  --gpu_per_rs 1 --bind=proportional-packed:${cores_per_rank} -d plane:1  ../Teton/SuOlsonTest ${grid} 16 2 16 8 4

# jsrun --smpiargs="-mxm --mca btl_openib_warn_default_gid_prefix 0 --mca mpi_warn_on_fork 0" #   --nrs 2  --tasks_per_rs 2 --cpu_per_rs 20 #   --gpu_per_rs 2 --bind=proportional-packed:10 -d plane:2  #   ../../sm60/Teton/SuOlsonTest 2x2x1_34.cmg 16 2 16 8 4



EOF
#---------------------------------------
bsub < batch.job
#bsub  batch.job
