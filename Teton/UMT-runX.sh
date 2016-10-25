#!/bin/bash

##PBS -A CSC122
##PBS -l gres=widow2
#PBS -j oe
#PBS -k oe
#PBS -N UMT2013
#PBS -l walltime=01:00:00
#PBS -V


# set to number of GPUs per node you want to use 
export NUMGPUS=4
# set to number of ranks per node you want to use
export RPN=32
# enable if you want to capture an nvprof timeline
#export PROFILING=1

# dimensions for the grid
export Px=4
export Py=4
export Pz=2
export Xzones=12
export Yzones=12
export Zzones=12

Order=16
Groups=16
quadType=2
Polar=8
Azim=4
Ranks=$((${Px}*${Py}*${Pz}))

cd $PBS_O_WORKDIR

gridFileName=grid_${Ranks}MPI_${Xzones}x${Yzones}x${Zzones}.cmg
./gridfile.sh $gridFileName $Px $Py $Pz $Xzones $Yzones $Zzones

#
# Note: these mappings are specific to a particular machine
# configuration with two 8-core CPUs and two GPUs per CPU
#
if [[ $NUMGPUS -eq 1 ]]; then
    export MV2_CPU_MAPPING=0:1:2:3:4:5:6:7:16:17:18:19:20:21:22:23
elif [[ $NUMGPUS -eq 2 ]]; then
    export MV2_CPU_MAPPING=0:8:1:9:2:10:3:11:4:12:5:13:6:14:7:15:16:24:17:25:18:26:19:27:20:28:21:29:22:30:23:31
elif [[ $NUMGPUS -ge 3 ]]; then
    export MV2_CPU_MAPPING=0:1:8:9:2:3:10:11:4:5:12:13:6:7:14:15:16:17:24:25:18:19:26:27:20:21:28:29:22:23:30:31
fi

echo "Params: $NUMGPUS $RPN $Px $Py $Pz $Xzones $Yzones $Zzones"

mpirun -n $Ranks -hostfile ${PBS_NODEFILE} ./mps.sh ./SuOlsonTest $gridFileName $Groups $quadType $Order $Polar $Azim


