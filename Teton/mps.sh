#! /bin/sh
#
# Wrapper script for UMT2013 that allows the use of
# CUDA Multi-Process Services (MPS) with multiple GPUs
# per node (not formally supported CUDA 6.0, but works)
#

#
# Determine local rank ID (method depends on MPI implementation)
#
MYID=
if   [ "$PMI_ID " != " " ] ; then
    MYID=$PMI_ID
elif [ "$OMPI_COMM_WORLD_LOCAL_RANK " != " " ] ; then
    MYID=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ "$MV2_COMM_WORLD_LOCAL_RANK " != " " ] ; then
    MYID=$MV2_COMM_WORLD_LOCAL_RANK
elif [ "$PBS_VNODENUM " != " " ] ; then
    MYID=$PBS_VNODENUM
fi

if [ "$MYID " == " " ]; then
   echo "Unable to determine local rank ID - exiting"
   exit 1
fi

#
# Query number of GPUs actually present; number you wish to
# use (specified with NUMGPUS) might be fewer.
#
NUMREALGPUS=`nvidia-smi -q | grep ^GPU | wc -l`

if [ "$NUMGPUS " == " " ]; then
   echo "Unknown number of GPUs - exiting"
   exit 1
fi

#
# Map ranks to GPUs
#
HOST=`hostname`
if [ `expr \( $NUMREALGPUS / $NUMGPUS \)` -lt 2 ]; then
   GPUID=`expr \( $MYID % $NUMGPUS \)`
else
   GPUID=`expr \( $MYID % $NUMGPUS \) \* \( $NUMREALGPUS / $NUMGPUS \)`
fi
export CUDA_MPS_PIPE_DIRECTORY=/tmp/$USER.mps$GPUID

#
# Configure the GPUs for persistence mode, max boost clocks, and
# compute-exclusive mode, and start up an instance of the CUDA MPS
# server for each GPU
#
if [ $MYID -lt $NUMGPUS ]; then
   echo "$HOST $MYID: enabling persistence mode"
   sudo nvidia-smi -pm 1 -g $GPUID
   echo "$HOST $MYID: setting clocks"
   sudo nvidia-smi -ac 3004,875 -g $GPUID
   echo "$HOST $MYID: setting exclusive mode"
   sudo nvidia-smi -c 3 -g $GPUID
   echo "$HOST $MYID: starting up CUDA MPS"
   export CUDA_VISIBLE_DEVICES=$GPUID
   mkdir -p $CUDA_MPS_PIPE_DIRECTORY
   nvidia-cuda-mps-control -d
   unset CUDA_VISIBLE_DEVICES

   if [ "$PROFILING " == "1 " ] && [ $MYID -eq 0 ]; then
      echo "$HOST $MYID: starting profiler"
      nvprof -o "UMT2013-$NUMGPUS-$RPN-%p.nvp" --profile-all-processes &
   fi
else
   echo "$HOST $MYID: waiting"
   sleep 5
fi

echo "$HOST: lrank $MYID, gpu $GPUID"

#
# Launch the application
#
$*

#
# Done -- stop profiling, shutdown MPS servers, restore state of GPUs
#
sleep 5
if [ "$PROFILING " == "1 " ] && [ $MYID -eq 0 ]; then
   echo "$HOST $MYID: stopping profiler"
   killall -HUP nvprof
fi

if [ $MYID -lt $NUMGPUS ]; then
   echo "$HOST $MYID: stopping CUDA MPS"
   echo quit | nvidia-cuda-mps-control
   rm -rf $CUDA_MPS_PIPE_DIRECTORY
   echo "$HOST $MYID: resetting default mode"
   sudo nvidia-smi -c 0 -g $MYID
   echo "$HOST $MYID: resetting clocks"
   sudo nvidia-smi -rac -g $MYID
fi
