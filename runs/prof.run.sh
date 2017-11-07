#!/bin/bash
nodes=1
ppn=1
let nmpi=$nodes*$ppn
grid=1x1x1profile.cmg
order=16
groups=32
type=2
polar=4
azim=4
#------------------------
export OMP_NUM_THREADS=20
export OMP_WAIT_POLICY=active
export HPM_GROUP_LIST=10,21
mpirun --bind-to none -np $nmpi ./bind.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim
