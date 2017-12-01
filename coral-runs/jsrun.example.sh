#!/bin/bash
nodes=1
ppn=12
let nmpi=$nodes*$ppn
grid=3x2x2_34.cmg
order=16
#groups=32
groups=16
type=2
polar=8
azim=4

export LIBRARY_PATH=/opt/pgi/linuxpower/17.9/lib:/shared/opt/xlf/16.1.0.0-171111a/xlf/16.1.0/lib:/shared/opt/xlf/16.1.0.0-171111a/lib:/shared/opt/xlC/14.1.0.0-171111a/xlC/14.1.0/lib:/shared/opt/xlC/14.1.0.0-171111a/lib:/opt/ibm/spectrum_mpi/lib

export LD_LIBRARY_PATH=/opt/ibm/spectrum_mpi/jsm_pmix/lib/:/opt/ibm/spectrum_mpi/jsm_pmix/../lib/:/opt/ibm/csm/lib/:/opt/ibm/spectrum_mpi/lib/pami_port:/opt/mellanox/hcoll/lib:/shared/lsf/10.1/linux3.10-glibc2.17-ppc64le/lib:/opt/pgi/linuxpower/17.9/lib:/shared/opt/xlf/16.1.0.0-171111a/xlf/16.1.0/lib:/shared/opt/xlf/16.1.0.0-171111a/lib:/shared/opt/xlC/14.1.0.0-171111a/xlC/14.1.0/lib:/shared/opt/xlC/14.1.0.0-171111a/lib:/opt/ibm/spectrum_mpi/lib:/usr/local/cuda/lib64

export OMP_NUM_THREADS=12
export USE_MPS=yes
export RANKS_PER_NODE=$ppn

jsrun -E LD_LIBRARY_PATH -E LIBRARY_PATH -E USE_MPS -E RANKS_PER_NODE --rs_per_host $ppn --np $nmpi ./umt_jsmhelper_6gpus.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim 2> jsm_6gpu.err 1> jsm_6gpu.out

# There seems to be a problem getting the correct environment passed through to jsm and jsrun. UMT starts but does not get far. May work better with correctly integrated CRM environment at the labs.
