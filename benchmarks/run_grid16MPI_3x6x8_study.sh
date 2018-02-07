#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=4
#SBATCH --time=00:30:00
#SBATCH --job-name=umt-small

Px=4
Py=2
Pz=2

Groups=200

quadType=2
Order=16
Polar=9
Azim=10

Xzones=3
Yzones=6
Zzones=8
Ranks=$((${Px}*${Py}*${Pz}))

#
# Input: SuOlsonTest $gridFileName $Groups $quadType $Order $Polar $Azim
# Allowed values: 1 <= quadType <= 2; 1 <= Polar <= 18; 1 <= Azim <= 22 
#
gridFileName=grid16MPI_3x6x8.cmg
export KMP_AFFINITY=compact,granularity=core,1
for t in 1 2 4 6; do
  echo "---> Running with OMP_NUM_THREADS=$t"
  export OMP_NUM_THREADS=$t
  srun -n 16 -c 6 ../Teton/SuOlsonTest $gridFileName $Groups $quadType $Order $Polar $Azim
done

