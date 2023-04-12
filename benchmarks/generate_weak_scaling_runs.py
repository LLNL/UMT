# Generate runs for a weak scaling run.
#
# This python script can be used to assist in generating a series of run configurations for the UMT unstructured box 3d mesh.
# The series of runs are designed to weak scale the problem from 1 to n processors.
#
# node_memory - total memory on a node (gigabytes)
# node_memory_target - % of memory to use on node at max number of MPI ranks.
# memory_tolerance - tolerance for problem memory usage.  Will try to generate runs that use the same amount of memory per rank to within X%.

# Node has 128gb of memory
node_memory = 128
# Max # mpi ranks to use ( typically equal to # cores on node ).
max_ranks = 36
# Use 50% of node's memory, within +/- 5%.
node_memory_target = 0.50
memory_tolerance = 0.05

# Number of polar x azimuthal angles to use in the product quadrature set
a = 3
# Number of energy groups to use
g = 32

# upper and lower bound of memory usage to calculate series for, in bytes, per rank
mem_lower_bound = 128 * node_memory_target * (1 - memory_tolerance) * 1024 * 1024 * 1024 / max_ranks
mem_upper_bound = 128 * node_memory_target * (1 + memory_tolerance) * 1024 * 1024 * 1024 / max_ranks

print ("#---------------------------")
print ("#Angles: 3x3   Groups: 32")
print ("#---------------------------")
for r in range(1,1000):
   # 3d mesh
   zones = 12*(r**3)
   angles = a * a * 8
   totalDOFs = zones * 8 * g * angles 
   totalMem = totalDOFs * 8
   for p in range(1,max_ranks):
      mem_per_rank = totalMem / p
      if mem_per_rank >= mem_lower_bound and mem_per_rank <= mem_upper_bound:
         DOFsPerRank = int(totalDOFs / p)
         if 8*p < zones:
            print(f"srun -N1-1 -n {p} --exclusive test_driver -c 1 -G {g} -A {a} -P {a} -r 1 -R {r} -i unstructBox3D.mesh >& run.{a}x{a}_product_quad.{angles}_angles.{g}_groups.{p}_ranks.{DOFsPerRank}_dofs_per_rank.{zones}_zones.{totalDOFs}_total_dofs.log")
