rzoz-login1{langer1}24: cd Teton/
/g/g17/langer1/UMT2013/Teton
rzoz-login1{langer1}25: gmake SuOlsonTest
mpiicpc -c   -I/usr/local/tools/mvapich-intel/include -I. -I../cmg2Kull/sources -I../CMG_CLEAN/src -I/g/g17/langer1/UMT2013/Teton/utilities -g -O3 -w  -no-prec-div -fPIC  -openmp   -DLINUX -DLinux -DUSE_MPI -DMPICH_SKIP_MPICXX  SuOlsonTest.cc -o SuOlsonTest.o
gmake -C geom
gmake[1]: Entering directory `/g/g17/langer1/UMT2013/Teton/geom'
gmake -C CMI
gmake[2]: Entering directory `/g/g17/langer1/UMT2013/Teton/geom/CMI'
gmake[2]: Nothing to be done for `all'.
gmake[2]: Leaving directory `/g/g17/langer1/UMT2013/Teton/geom/CMI'
gmake -C Field
gmake[2]: Entering directory `/g/g17/langer1/UMT2013/Teton/geom/Field'
gmake[2]: Nothing to be done for `all'.
gmake[2]: Leaving directory `/g/g17/langer1/UMT2013/Teton/geom/Field'
gmake -C Region
gmake[2]: Entering directory `/g/g17/langer1/UMT2013/Teton/geom/Region'
gmake[2]: Nothing to be done for `all'.
gmake[2]: Leaving directory `/g/g17/langer1/UMT2013/Teton/geom/Region'
gmake[1]: Leaving directory `/g/g17/langer1/UMT2013/Teton/geom'
gmake -C transport
gmake[1]: Entering directory `/g/g17/langer1/UMT2013/Teton/transport'
gmake -C Teton
gmake[2]: Entering directory `/g/g17/langer1/UMT2013/Teton/transport/Teton'
gmake[2]: Nothing to be done for `all'.
gmake[2]: Leaving directory `/g/g17/langer1/UMT2013/Teton/transport/Teton'
gmake -C TetonInterface
gmake[2]: Entering directory `/g/g17/langer1/UMT2013/Teton/transport/TetonInterface'
gmake[2]: Nothing to be done for `all'.
gmake[2]: Leaving directory `/g/g17/langer1/UMT2013/Teton/transport/TetonInterface'
gmake[1]: Leaving directory `/g/g17/langer1/UMT2013/Teton/transport'
gmake -C part
gmake[1]: Entering directory `/g/g17/langer1/UMT2013/Teton/part'
gmake[1]: Nothing to be done for `all'.
gmake[1]: Leaving directory `/g/g17/langer1/UMT2013/Teton/part'
gmake -C communication
gmake[1]: Entering directory `/g/g17/langer1/UMT2013/Teton/communication'
gmake[1]: Nothing to be done for `all'.
gmake[1]: Leaving directory `/g/g17/langer1/UMT2013/Teton/communication'
mpiicpc   -I/usr/local/tools/mvapich-intel/include -I. -I../cmg2Kull/sources -I../CMG_CLEAN/src -I/g/g17/langer1/UMT2013/Teton/utilities -g -O3 -w  -no-prec-div -fPIC  -openmp   SuOlsonTest.o -o SuOlsonTest   -openmp  -L../CMG_CLEAN/src -Wl,-rpath,/g/g17/langer1/UMT2013/Teton/../CMG_CLEAN/src -L../cmg2Kull/sources -Wl,-rpath,/g/g17/langer1/UMT2013/Teton/../cmg2Kull/sources -L../Teton -Wl,-rpath,/g/g17/langer1/UMT2013/Teton/../Teton -L/usr/local/tools/mvapich-intel/lib  -lfmpich -lmpichfarg -lmpich  -libverbs -lpthread -libumad  -ldl   -lInfrastructure -lTetonUtils -lc2k -lcmgp  -lifcore -lm -lc -lsvml  -lfmpich -lmpichfarg -lmpich  -libverbs -lpthread -libumad   
