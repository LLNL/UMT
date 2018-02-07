# Run this in the top-level UMT2013 directory
#
# Make sure make.defs reflects the platform's compilers, compiler options, libraries, MPI wrappers etc. 
#
#use ic-14.0.211
#use ifort-14.0.211
use -q ic-13.1.163
use -q ifort-13.1.163

gmake clean
gmake

cd ./Teton
gmake SuOlsonTest 
