#include "TetonUtilities.hh"

#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_blueprint_mesh.hpp"
#include "conduit/conduit_relay_mpi.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib> // For std::getenv and std::system
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <sched.h> // For sched_getcpu()
#include <set>

#if defined(TETON_ENABLE_OPENMP)
#include <omp.h>
#endif

#if defined(TETON_ENABLE_CUDA)
#include <cuda_runtime.h>
#elif defined(TETON_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

extern "C"
{
void teton_print_thread_bindings_c()
{
   Teton::utilities::printThreadBindings();
}
int teton_get_gpu_processor_count_c()
{
   return Teton::utilities::getGPUProcessorCount();
}
}

namespace Teton
{

namespace utilities
{

void convert_int32(int /*rank*/, conduit::Node &root, const std::vector<std::string> &keys)
{
   for (const auto &path : keys)
   {
      if (root.has_path(path))
      {
         conduit::Node &n = root.fetch_existing(path);
         if (!n.dtype().is_int32())
         {
            // Convert the data to int32 and put it back in node n. I tried move/swap
            // to try and steal ifield's data rather than copying it again but that
            // caused n's name to be blank.
            conduit::Node ifield;
            n.to_int32_array(ifield);
            n.set(ifield);
         }
      }
   }
}

void find_dtype(const conduit::Node &n,
                const conduit::DataType &dtype,
                const std::string &path,
                std::vector<std::string> &paths)
{
   auto concat = [](const std::string &base, const std::string &name)
   {
      if (base.empty())
         return name;
      return base + "/" + name;
   };

   if (n.number_of_children() > 0)
   {
      // Make paths be relative to node n at the top level.
      for (conduit::index_t i = 0; i < n.number_of_children(); i++)
      {
         if (path == "")
            find_dtype(n[i], dtype, n[i].name(), paths);
         else
            find_dtype(n[i], dtype, concat(path, n[i].name()), paths);
      }
   }
   else
   {
      if (n.dtype().id() == dtype.id())
         paths.push_back(path);
   }
}

std::vector<std::string> find_int64(const conduit::Node &n)
{
   std::vector<std::string> paths;
   find_dtype(n, conduit::DataType::int64(), "", paths);
   return paths;
}

bool scan_field_values(int rank, const conduit::Node &n)
{
   bool retval = true;
   if (n.dtype().is_float64())
   {
      int count = 0;
      conduit::float64_array arr = n.value();
      for (conduit::index_t i = 0; i < arr.number_of_elements(); i++)
      {
         if (!std::isfinite(arr[i]))
         {
            std::cout << rank << ":" << n.path() << ": elem[" << i << "] is not a number. " << arr[i] << std::endl;
            retval = false;
            count++;
            if (count > 10)
               break;
         }
      }
   }
   return retval;
}

bool find_local_duplicate_points(int domainId,
                                 const conduit::Node &dom,
                                 const conduit::Node &coordset,
                                 conduit::Node &info)
{
   bool retval = false;
   // Make sure Conduit is new enough.
#if (CONDUIT_VERSION_MAJOR == 0 && CONDUIT_VERSION_MINOR >= 9) || (CONDUIT_VERSION_MAJOR > 0)
   using conduit::index_t;

   // See whether any of the points in the local domain are duplicated.
   // If a point's query result does not equal its query index then it
   // must have been defined once before.
   conduit::blueprint::mesh::utils::query::PointQuery localPQ(dom);
   const index_t npts = conduit::blueprint::mesh::coordset::length(coordset);
   for (index_t pi = 0; pi < npts; pi++)
   {
      auto pt = conduit::blueprint::mesh::utils::coordset::_explicit::coords(coordset, pi);
      double pt3[3];
      pt3[0] = pt[0];
      pt3[1] = (pt.size() > 1) ? pt[1] : 0.;
      pt3[2] = (pt.size() > 2) ? pt[2] : 0.;
      localPQ.add(domainId, pt3);
   }
   localPQ.execute(coordset.name());
   for (index_t pi = 0; pi < npts; pi++)
   {
      const auto &res = localPQ.results(static_cast<int>(domainId));
      if (res[pi] != pi)
      {
         const auto pts = localPQ.inputs(domainId);
         double pt3[3]{pts[3 * pi], pts[3 * pi + 1], pts[3 * pi + 2]};
         std::stringstream ss;
         ss << "Domain " << domainId << " duplicated point " << pi << " (" << pt3[0] << ", " << pt3[1] << ", " << pt3[2]
            << ") at " << res[pi] << ".";

         conduit::Node &vn = info.append();
         vn["message"].set(ss.str());
         vn["vertex"] = pi;
         vn["duplicate_vertex"] = res[pi];
         vn["coordinate"].set(pt3, 3);

         retval = true;
      }
   }
#endif
   return retval;
}

std::vector<std::string> globalizeStringVector(const std::vector<std::string> &vec, MPI_Comm comm)
{
   // Make a Conduit node from it.
   conduit::Node send_node;
   for (const auto &value : vec)
      send_node[value] = 1;

   // Send the data to all ranks.
   conduit::Node recv_node;
   conduit::relay::mpi::all_gather_using_schema(send_node, recv_node, comm);

   // Pick through the output and make the output string vector from the node names.
   std::set<std::string> unique;
   for (conduit::index_t i = 0; i < recv_node.number_of_children(); i++)
   {
      const conduit::Node &child = recv_node[i];
      for (conduit::index_t j = 0; j < child.number_of_children(); j++)
         unique.insert(child[j].name());
   }
   std::vector<std::string> retval;
   retval.insert(retval.begin(), unique.begin(), unique.end());

   return retval;
}

//------------------------------------------------------------------------------
int Banner::level = 0;

Banner::Banner(MPI_Comm c, const std::string &str) : comm(c), rank(0), name(str)
{
   MPI_Comm_rank(comm, &rank);
   MPI_Barrier(comm);
   if (rank == 0)
      printLine(name + " (start)");
   MPI_Barrier(comm);
   level++;
}

Banner::~Banner()
{
   level--;
   MPI_Barrier(comm);
   if (rank == 0)
      printLine(name + " (end)");
   MPI_Barrier(comm);
}

void Banner::printLine(const std::string s) const
{
   int n = std::max(2, (80 - 2 - static_cast<int>(s.size())) / 2);
   emit(' ', level * 2);
   emit('=', n);
   std::cout << " " << s << " ";
   emit('=', n);
   std::cout << std::endl;
}

void Banner::emit(char c, int n) const
{
   for (int i = 0; i < n; i++)
      std::cout << c;
}

// Utility function to scan a vector of integers and return
// the ranges in []'s as a string, for printing to output.
std::string getConsecutiveRanges(const std::vector<int> &cores)
{
   if (cores.empty())
      return "";

   std::string result;
   int start = cores[0];
   int end = cores[0];

   for (size_t i = 1; i < cores.size(); ++i)
   {
      if (cores[i] == end + 1)
      {
         // Extend the range
         end = cores[i];
      }
      else
      {
         // Append the current range to the result and start a new range
         if (start == end)
         {
            result += std::to_string(start) + ",";
         }
         else
         {
            result += std::to_string(start) + "-" + std::to_string(end) + ",";
         }
         start = cores[i];
         end = cores[i];
      }
   }

   // Append the final range
   if (start == end)
   {
      result += std::to_string(start);
   }
   else
   {
      result += std::to_string(start) + "-" + std::to_string(end);
   }

   return result;
}

// Function to check and print the thread bindings, as well as any visible GPUs
void printThreadBindings()
{
#if (TETON_ENABLE_OPENMP)
   // Retrieve environment variables
   const char *rocrVisibleDevices = std::getenv("ROCR_VISIBLE_DEVICES");
   const char *cudaVisibleDevices = std::getenv("CUDA_VISIBLE_DEVICES");

   std::string gpuVisibleDevices;

   if (rocrVisibleDevices)
   {
      gpuVisibleDevices = rocrVisibleDevices;
   }

   else if (cudaVisibleDevices)
   {
      gpuVisibleDevices = cudaVisibleDevices;
   }

   int num_threads = omp_get_max_threads();

   // Vector to store the CPU core each thread ran on
   std::vector<int> cpu_cores(num_threads, -1);

// Parallel region
#pragma omp parallel
   {
      // Get the thread ID
      int thread_id = omp_get_thread_num();

      // Get the CPU core the thread is running on
      int cpu_core = sched_getcpu();

      // Store the core information in the vector
      cpu_cores[thread_id] = cpu_core;
   }

   // Use a set to find all unique CPU cores
   std::set<int> unique_cpu_cores(cpu_cores.begin(), cpu_cores.end());

   // Convert the set to a sorted vector
   std::vector<int> sorted_cores(unique_cpu_cores.begin(), unique_cpu_cores.end());

   std::string core_ranges = getConsecutiveRanges(sorted_cores);

   // Print list of visible GPUs
   std::cout << "Threads bound to cores: " << core_ranges << ", visible GPU ids: " << gpuVisibleDevices << std::endl;
#endif
}

// Retrieves the number of processors on the GPU.  Currently checks the first GPU visible to a process.
int getGPUProcessorCount()
{
   int nProcs = 0;

#if defined(TETON_ENABLE_CUDA)
   cudaDeviceProp prop;
   cudaError_t err = cudaGetDeviceProperties(&prop, 0);
   if (err != cudaSuccess)
   {
      std::cerr << "Teton failed to query device # processors: CUDA error: " << cudaGetErrorString(err) << std::endl;
      abort();
   }
   nProcs = prop.multiProcessorCount;

#elif (TETON_ENABLE_HIP)
   hipDeviceProp_t prop;
   hipError_t err = hipGetDeviceProperties(&prop, 0);
   if (err != hipSuccess)
   {
      std::cerr << "Teton failed to query device # processors: HIP error: " << hipGetErrorString(err) << std::endl;
      abort();
   }
   nProcs = prop.multiProcessorCount;
#endif

   return nProcs;
}

// This performs the same operation as the energy check in the Fortran written by Ben Yee.
bool checkEnergyConservation(int rank, const conduit::Node &datastore)
{
   bool result = false;

   double energy_check_tolerance = datastore.fetch_existing("options/iteration/relativeTolerance").value();
   energy_check_tolerance *= 10.0;
   double energy_radiation = datastore.fetch_existing("rtedits/EnergyRadiation").value();
   double energy_check = datastore.fetch_existing("rtedits/EnergyCheck").value();
   double rel_energy_check_result = std::abs(energy_check / (energy_radiation + 1.0e-50));

   if (rel_energy_check_result <= energy_check_tolerance)
   {
      result = true;
   }
   else
   {
      result = false;

      std::cerr << "Teton: Failed energy conservation check on rank " << rank << ". Relative difference of "
                << std::scientific << rel_energy_check_result << " exceeds tolerance of " << std::scientific
                << energy_check_tolerance << std::endl;

      double power_incident = datastore.fetch_existing("rtedits/PowerIncident").value();
      double power_escape = datastore.fetch_existing("rtedits/PowerEscape").value();
      double power_absorbed = datastore.fetch_existing("rtedits/PowerAbsorbed").value();
      double power_emitted = datastore.fetch_existing("rtedits/PowerEmitted").value();

      // Energy deposited in material =   -2.4290858270E+07 ERad total =    7.5347327996E+08 Energy check =  -1.8626451492E-07
      // TODO - Check with Ben if I should make these prints better match the Fortran (example above in above comment ).
      std::cerr << "Teton:: Energy radiation: " << std::scientific << energy_radiation << std::endl;
      std::cerr << "Teton:: Power incident: " << std::scientific << power_incident << std::endl;
      std::cerr << "Teton:: Power escaped: " << std::scientific << power_escape << std::endl;
      std::cerr << "Teton:: Power absorbed: " << std::scientific << power_absorbed << std::endl;
      std::cerr << "Teton:: Power emitted: " << std::scientific << power_emitted << std::endl;
   }

   return result;
}

} // namespace utilities

} // namespace Teton
