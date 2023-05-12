# Generate runs for a strong scaling study.
#
# This python script can be used to assist in generating a series of run configurations for the UMT unstructured box 3d mesh.
# The series of runs are designed to strong scale the problem from 1 to n processors.
#
# node_memory - total memory on a node (gigabytes)
# node_memory_target - % of memory to use on node at max number of MPI ranks.
# memory_tolerance - tolerance for problem memory usage.  Will try to generate runs that use the same amount of memory per rank to within X%.

# Node has 256gb of memory
node_memory = 256
# Use 50% of node's memory, within +/- 5%.
node_memory_target = 0.50
memory_tolerance = 0.02
max_ranks = 112
output_dir = "."

# Number of polar by azimuthal angles to use in the product quadrature set
a = 3
# Number of energy groups to use
g = 128

# upper and lower bound of memory usage to calculate series for, in bytes.
mem_lower_bound = node_memory * node_memory_target * (1 - memory_tolerance) * 1024 * 1024 * 1024
mem_upper_bound = node_memory * node_memory_target * (1 + memory_tolerance) * 1024 * 1024 * 1024

print ("# Calculating size of mesh to use for your angle and group count...")
for r in range(1,1000):
   # 3d mesh
   zones = 12*(r**3)
   angles = a * a * 8
   totalDOFs = zones * 8 * g * angles 
   mem = totalDOFs * 8
   memGB = mem/1024/1024/1024
   if mem >= mem_lower_bound and mem <= mem_upper_bound:
      print(f"#You will need to refine the mesh {r} times and generate a mesh of {zones} zones to use {memGB} GB of memory for your problem configuration.")
      print(f"./makeUnstructuredBox")
      print(f"srun -n1 --exclusive ./test_driver -i unstructBox3D.mesh -c 0 -r 1 -R {r} -o {output_dir}")

for r in range(1,max_ranks+1):
   print(f"srun -n {r} --exclusive ./test_driver -c 1 -b 1 -i {output_dir}/refined_mesh.mesh -o {output_dir} >& run.spp1.{r}_ranks.{zones}_zones.{totalDOFs}_dofs.log")
