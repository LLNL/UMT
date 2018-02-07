#!/bin/bash
#SBATCH -J UMT2015-ATS3
#SBATCH -N 1152
#SBATCH -p regular
#SBATCH --time=00:30:00

export OMP_NUM_THREADS=4

export Order=16
export Groups=200
export quadType=2
export Polar=9
export Azim=10

export gridFileName=grid6912MPI_3x6x8.cmg

echo "Running job in $SLURM_SUBMIT_DIR .. "
cd $SLURM_SUBMIT_DIR

echo "Running deck: $gridFileName Order=$Order Groups=$Groups quadType=$quadType Polar=$Polar Azim=$Aim Thr=$OMP_NUM_THREADS"

srun -N 1152 -n 6912 -c $OMP_NUM_THREADS ./SuOlsonTest $gridFileName $Groups $quadType $Order $Polar $Azim

# aprun -n 1 -d $OMP_NUM_THREADS -j 1 -S 1 -ss ./SuOlsonTest $gridFileName $Groups $quadType $Order $Polar $Azim
# srun -n 1 -N 1 ./SuOlsonTest $gridFileName $Groups $quadType $Order $Polar $Azim
#
# where:
#    gridFileName=grid1MPI_3x3x4.cmg   # this file is included.
#    Order=16
#    Groups=200
#    quadType=2
#    Polar=9
#    Azim=10

