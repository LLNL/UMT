#!/bin/bash
nodes=1
ppn=1
let nmpi=$nodes*$ppn
grid=1x1x1_12.cmg
order=16
#groups=192
groups=32
type=2
polar=4
azim=4
#--------------------------------------
export OMP_NUM_THREADS=20
export OMP_WAIT_POLICY=active
export HPM_GROUP_LIST=10,21


mpirun --bind-to none -np $nmpi /home/walkup/bin/set_device_and_bind.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim


#mpirun -np $nmpi ./bind.sh cuda-memcheck ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim 

#mpirun --bind-to none -np $nmpi ./bind.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim   

#mpirun --bind-to none -np $nmpi ./bind.sh nvprof -s -f -o 1node.%q{OMPI_COMM_WORLD_RANK}.nvprof ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

#nvidia-smi -q -l 1 -d memory

#mpirun --bind-to none -np $nmpi ./bind.sh nvprof -s ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim 

#mpirun --bind-to none -np $nmpi /home/walkup/bin/set_device_and_bind.sh gdb --args ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim    

#/gpfs/ess2fs0/walkup/openmpi-1.8.8/bin/mpirun -np $nmpi ./bind.sh cuda-gdb --args ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim

#---------------------------------------                                                                                                      
