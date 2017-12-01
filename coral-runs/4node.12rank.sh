#!/bin/bash
nodes=4
ppn=12
#ppn=6
let nmpi=$nodes*$ppn
grid=3x4x4_32.cmg
#grid=3x2x1_20.cmg
order=16
#groups=32
groups=16
type=2
polar=8
azim=4
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=${ppn}]"
#BSUB -R "rusage[ngpus_shared=6]"
#BSUB -n ${nmpi}
#BSUB -x
#BSUB -q excl_ws_dd21
##BSUB -R "select[fr=5]"
#BSUB -W 120
#BSUB -env "all,LSF_CPU_ISOLATION=on,LSF_IRQ_ISOLATION=on, LSF_START_JOBS_MPS=N"
#BSUB -R "select[hname!='c699c026']"
#BSUB -R "select[hname!='c699c072']"
#BSUB -R "select[hname!='c699c032']"
#BSUB -R "select[hname!='c699c007']"
#BSUB -R "select[hname!='c699c016']"
#BSUB -R "select[hname!='c699c017']"
#BSUB -R "select[hname!='c699c025']"
#BSUB -R "select[hname!='c699c013']"
#BSUB -R "select[hname!='c699c018']"

#---------------------------------------
#export OMP_NUM_THREADS=22
export OMP_WAIT_POLICY=active
#export HPM_GROUP_LIST=10,21

export PROFILE_RANK=-1
export USE_MPS=yes

export TMPDIR=/tmp

df -h /tmp
echo 'starting mpirun'

ls /dev/cpuset

mpirun --bind-to none -np $nmpi ./p9_helper_6gpu.sh ../Teton/SuOlsonTest $grid $groups $type $order $polar $azim


EOF
#---------------------------------------
bsub <batch.job
