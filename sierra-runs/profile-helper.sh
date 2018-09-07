#!/bin/bash

world_rank=$PMIX_RANK

export CUDA_CACHE_PATH=/dev/shm/$USER/nvcache_$PMIX_RANK

echo $CUDA_CACHE_PATH

executable=$1

shift

if [ $world_rank == $PROFILE_RANK ]; then
#    nvprof -f -o $PROFILE_PATH $executable "$@"
    cuda-memcheck $executable "$@"
    #nvprof --print-gpu-trace $executable "$@"
else
    $executable "$@"
fi

echo "Cleaning up nvcache path"

rm -rf /dev/shm/$USER/nvcache_$PMIX_RANK
