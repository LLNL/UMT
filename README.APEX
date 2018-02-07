
===========================================================
UMT2015 Benchmark for APEX Procurement
===========================================================

The UMT2015 benchmark is an updated source code for the APEX procurement. The
code is developed and maintained by Lawrence Livermore National Laboratory
(LLNL), California.

===========================================================

== I.    Description: ==

UMT2015 is an LLNL ASC proxy application (mini-app) that performs
three-dimensional, non-linear, radiation transport calculations
using deterministic methods. The method of solution is as follws.
UMT performs the solution of time-dependent, energy-dependent,
discrete ordinates, and nonlinear radiation problems on three-dimensional
unstructured spatial grids on large distributed-memory multi-node
parallel computer systems with multiple cores per node.
To achieve extreme scalability, the application exploits both
spatial decomposition using message passing between nodes and a
threading algorithm in angle within the node.

The project addresses the need to obtain accurate transport solutions to
three-dimensional radiation transfer problems (i.e. the transfer of
thermal photons). The solution calls for the use of an unstructured grid.
This class of problems is characterized by tens of thousands of unknowns
per zone and upwards of millions of zones, thus requiring large, scalable,
parallel computing platforms. The package uses a combination of
message passing and threading, utilizes the large distributed memory
of the platform and provides unprecedented (weak) scaling
to very large core counts.

For further detail see:

P.F. Nowak and M.K. Nemanic, "Radiation Transport Calculations on Unstructured
Grids Using a Spatially Decomposed and Threaded Algorithm," Proc. Int. Conf.
Mathematics and Computation, Reactor Physics and Environmental Analysis in
Nuclear Applications, Madrid, Spain, September 27-30, 1999, Vol. 1, p. 379 (1999).

Paul Nowak, "Deterministic Methods for Radiation  Transport: Lessons Learned and 
Future Directions", ASC Workshop on Methods for Computational Physics and Modern
Software Practices, March 2004, LLNL document UCRL-CONF-202912.


== II.   Parallelism: ==

UMT2015 implements MPI-based multi-node parallelism and OpenMP threading within
each MPI rank. OpenMP parallel-for regions and reductions are used in the code
with critical regions where needed to maintain correctness.


== III. How to Build: ==

(1) Create a make.defs file in the root directory of the source code. Example
definitions are provided for Linux and BlueGene platforms. The file must be
named "make.defs"

(2) Run make clean to ensure dependency files are cleaned correctly

(3) Run make

(4) At this point check that the following libraries (either shared or static)
are available:

 ./CMG_CLEAN/src/libcmgp.[so, a]
 ./Teton/libTetonUtils.[so, a]
 ./Teton/libInfrastructure.[so, a]
 ./cmg2Kull/sources/libc2k.[so, a]

(5) cd to the ./Teton directory

(6) Run make SulOlsonTest

(7) Check the SuOlsonTest binary is produced from the make


== IV.  How to Run: ==

The UMT application reads the problem size and problem description from a file
during load, additional parameters must be supplied on the command line.

The file format is as follows:

sms(Nx-1,Ny-1,Nz-1)

blk(on,0:Nx,0:Ny,0:Nz)

tag("xMinFaces",face,(0:0,0:Ny,0:Nz))
tag("xMaxFaces",face,(Nx:Nx,0:Ny,0:Nz))
tag("yMinFaces",face,(0:Nx,0:0,0:Nz))
tag("yMaxFaces",face,(0:Nx,Ny:Ny,0:Nz))
tag("zMinFaces",face,(0:Nx,0:Ny,0:0))
tag("zMaxFaces",face,(0:Nx,0:Ny,Nz:Nz))

numzones(Zx,Zy,Zz)

sub(10%,0:Nx,0:Ny,0:Nz, (7,0,0,0))
seed(10)

Where Nx, Ny and Nz are MPI process ranks in each dimension. Note that the
sms line requires the numbers to be 1 less than the blk line. The numzones
variable sets the number of zones per domain in pure MPI runs. By setting
OpenMP to use >1 threads, the work can be strong scaled within a rank.

Once a suitable deck has been defined the UMT code can be run in the 
following way:

mpirun -n <nodes> ./SuOlsonTest \
	$gridFileName $Groups \
	$quadType $Order \
	$Polar $Azim

Where:
	Order=16
    	Groups=200
    	quadType=2
    	Polar=9
    	Azim=10

Vendors are free to modify the MPI decomposition as required for their
architecture, BUT IT IS REQUIRED TO KEEP THE PROBLEM SIZE FIXED to the
same total number of zones. As such, if the MPI decomposition is changed
it is also necessary to change numzones(Zx,Zy,Zz) in the .cmg file.

It is strongly desirable to maintain roughly even MPI decomposition in each 
dimension where possible.

Reference Decks (can be found in the benchmarks directory):

- Small: (MPI Only, Single Node) grid16MPI_3x6x8.cmg
  (2,304 total zones)
- Large (Reference): grid6912MPI_3x6x8.cmg
  (995,328 total zones)
- Grand Challenge (APEX Target, Single Job Problem): grid55296MPI_3x6x8.cmg
  (7,962,624 total zones)


== V.   Reporting Results ==

UMT will run to complete and then report:

<job output .. >
cumulativeIterationCount= 45 cumulativeWorkTime=253.32 s.
figure of merit = 3.3099e+11

The Figure of Merit value should be reported in the response along with
details of all output produced by the application, any modifications to
the Makefiles (or make.defs) and modifications to the source code. In 
addition, the vendor should supply the grid file used for the FOM reported
runs.

