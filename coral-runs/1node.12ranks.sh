#!/bin/bash
nodes=1
ppn=12
let nmpi=$nodes*$ppn
grid=3x2x2_32.cmg
order=16
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
#BSUB -W 80
#BSUB -env "all,LSF_CPU_ISOLATION=on,LSF_IRQ_ISOLATION=on, LSF_START_JOBS_MPS=N"
#BSUB -R "select[hname!='c699c017']"
#BSUB -R "select[hname!='c699c025']"
#BSUB -R "select[hname!='c699c026']"

# #BSUB -R "select[hname!='c699c019']"
# #BSUB -R "select[hname!='c699c032']"
# #BSUB -R "select[hname!='c699c037']"
# #BSUB -R "select[hname!='c699c040']"
# #BSUB -R "select[hname!='c699c071']"
# #BSUB -R "select[hname!='c699c082']"
# #BSUB -R "select[hname!='c699c083']"
# #BSUB -R "select[hname!='c699c084']"
# #BSUB -R "select[hname!='c699c085']"
# #BSUB -R "select[hname!='c699c086']"
# #BSUB -R "select[hname!='c699c092']"
# #BSUB -R "select[hname!='c699c094']"
# #BSUB -R "select[hname!='c699c095']"
# #BSUB -R "select[hname!='c699c098']"
# #BSUB -R "select[hname!='c699c110']"
# #BSUB -R "select[hname!='c699c117']"
# #BSUB -R "select[hname!='c699c128']"
# #BSUB -R "select[hname!='c699c130']"
# #BSUB -R "select[hname!='c699c134']"
# #BSUB -R "select[hname!='c699c145']"
# #BSUB -R "select[hname!='c699c172']"
#---------------------------------------

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
