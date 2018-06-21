#!/bin/bash

world_rank=$PMIX_RANK

export CUDA_CACHE_PATH=/dev/shm/$USER/nvcache_$PMIX_RANK

echo $CUDA_CACHE_PATH

executable=$1

shift

if [ $world_rank == $PROFILE_RANK ]; then
    nvprof --profile-from-start off -f -o $PROFILE_PATH $executable "$@"
else
    $executable "$@"
fi

echo "Cleaning up nvcache path"

rm -rf /dev/shm/$USER/nvcache_$PMIX_RANK
