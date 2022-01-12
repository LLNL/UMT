#include <sstream>
#include <getopt.h>
#include <unistd.h> //for getopt, access
#include "mpi.h"
#include <sys/stat.h> //for mkdir

// XLF signal handler function to emit a stack trace.
#if defined(__ibmxl__)
#include "signal.h"
extern "C" void xl__trce(int, siginfo_t *, void *);
#endif

#if defined(TETON_ENABLE_OPENMP)
#include "omp.h"
#endif

#if defined(TETON_ENABLE_MEMUSAGES)
#include "MemUsages.h"
#endif

#if defined(TETON_ENABLE_UMPIRE)
#include "umpire/Umpire.hpp"
#include "umpire/strategy/QuickPool.hpp"
#include "umpire/strategy/ThreadSafeAllocator.hpp"
#endif

#if defined(TETON_ENABLE_CUDA)
#include "cuda_runtime_api.h"
#endif

#include "TetonInterface.hh"
#include "TetonConduitInterface.hh"

#include "conduit/conduit.hpp"
#include "conduit/conduit_relay.hpp"
#include "conduit/conduit_blueprint.hpp"

#if defined(TETON_ENABLE_MFEM)
#include "TetonBlueprint.hh"
#include "mfem.hpp"
#endif

#if defined(TETON_ENABLE_CALIPER)
#include "adiak.hpp"
#include "caliper/cali.h"
#include "caliper/cali-mpi.h"
#include "caliper/cali-manager.h"
#endif

// Utility function, check if string ends with another string.
bool endsWith(std::string const &fullString, std::string const &ending)
{
   if (fullString.length() >= ending.length())
   {
      return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
   }
   else
   {
      return false;
   }
}

int print_bytes_as_gb(const char *label, size_t bytes)
{
   double gbmem = ((double) bytes / (1024.0 * 1024.0 * 1024.0));
   fprintf(stdout, "%s: %7.4f GB\n", label, gbmem);
   return (0);
}

int print_gpu_mem(const char *label)
{
#if defined(TETON_ENABLE_CUDA)
   size_t free_on_gpu = 0;
   size_t total_on_gpu = 0;
   double gbtotal, gbfree, gbused;

   if (cudaMemGetInfo(&free_on_gpu, &total_on_gpu) != cudaSuccess)
   {
      printf("cudeMemGetInfo failed for GPU 0");
      return (1);
   }

   gbtotal = ((double) total_on_gpu) / (1024.0 * 1024.0 * 1024.0);
   gbfree = ((double) free_on_gpu) / (1024.0 * 1024.0 * 1024.0);
   gbused = ((double) total_on_gpu - free_on_gpu) / (1024.0 * 1024.0 * 1024.0);
   fprintf(stdout, "%s: total %7.4f GB; free %7.3f GB; used %7.3f GB\n", label, gbtotal, gbfree, gbused);
   fflush(stdout);
#endif
   return (0);
}

//==========================================================
//  run syntax:
//
// ./test_driver -c <# cycles> -i <input path> -o <output path>
//==========================================================

int main(int argc, char *argv[])
{
   int myRank = 0;
   int opt;
   int cycles = 1;
   int numPhaseAngleSets = 0;
   int useUmpire = 1;
   int numOmpMaxThreads = -1; // Max number of CPU threads to use  If -1, use value from omp_get_max_threads()
#if defined(TETON_ENABLE_MFEM)
   int numSerialRefinementLevels = 0;
   int numParallelRefinementLevels = 0;
   int numPolar = -1;
   int numAzimuthal = -1;
   int numGroups = 1;
#endif

   // MPI
   MPI_Comm comm = MPI_COMM_WORLD;
   int provided = 0;
   int claimed = 0;
   int request = MPI_THREAD_SINGLE;
   int verbose = 11;

   bool useGPU = false;
   bool useCUDASweep = false;
   bool dumpCopyOfInput = false;
   bool freeOMPResources = false;
   bool useNewGTASolver = false;
   bool useDeviceAwareMPI = false; // Pass device addresses to MPI ( device-aware MPI ).

   double memory_mb_used = -1;

   std::string inputPath(".");
   std::string outputPath("");
   std::string label("");
   std::string caliper_config("runtime-report");

#if defined(TETON_ENABLE_CUDA)
   caliper_config += ",nvprof";
#endif

//==========================================================
// Initialize MPI and MPI thread support.
//==========================================================
#if defined(TETON_ENABLE_OPENMP)
   request = MPI_THREAD_MULTIPLE;
#endif

#if defined(TETON_ENABLE_CALIPER)
   cali_mpi_init();
   cali::ConfigManager mgr;
   mgr.set_default_parameter("aggregate_across_ranks", "true");
   mgr.set_default_parameter("calc.inclusive", "true");
   mgr.set_default_parameter("main_thread_only", "true");
   auto configs = mgr.available_config_specs();
#endif

   if (MPI_Init_thread(&argc, &argv, request, &provided) != MPI_SUCCESS)
   {
      exit(1);
   }

   MPI_Query_thread(&claimed);

   if (provided < request)
   {
      std::cerr << "Teton driver: MPI_Init_thread was only able to provided thread level support " << provided
                << ".  Teton requested level " << request;
      exit(1);
   }
   if (claimed < request)
   {
      std::cerr << "Teton driver: MPI_Init_thread was only able to provided thread level support " << claimed
                << ".  Teton requested level " << request;
      exit(1);
   }

#if defined(TETON_ENABLE_CALIPER)
   adiak::init((void *) &comm);
   adiak::user();
   adiak::launchdate();
   adiak::launchday();
   adiak::executable();
   adiak::clustername();
   adiak::jobsize();
   adiak::hostlist();
   adiak::walltime();
   adiak::systime();
   adiak::cputime();
   adiak::value("TetonVersion", teton_get_version(), adiak_general, "TetonBuildInfo");
   adiak::value("TetonSHA1", teton_get_git_sha1(), adiak_general, "TetonBuildInfo");
   adiak::value("TetonCxxCompiler", teton_get_cxx_compiler(), adiak_general, "TetonBuildInfo");
   adiak::value("TetonFortranCompiler", teton_get_fortran_compiler(), adiak_general, "TetonBuildInfo");

   if (!label.empty())
   {
      adiak::value("label", label, adiak_general, "RunInfo");
   }
#endif

   MPI_Comm_rank(comm, &myRank);

#if defined(TETON_ENABLE_MEMUSAGES)
   MemUsage_init(comm);
#endif

   //==========================================================
   // Get command line arguments
   //==========================================================

   while (1)
   {
      static struct option long_options[]
          = { {"apply_label", no_argument, 0, 'l'},
              {"caliper", required_argument, 0, 'p'},
              {"dump_input_copy", no_argument, 0, 'd'},
              {"free_omp_resources", no_argument, 0, 'f'},
              {"help", no_argument, 0, 'h'},
              {"input_path", required_argument, 0, 'i'},
              {"num_cycles", required_argument, 0, 'c'},
              {"num_phase_space_sets", required_argument, 0, 's'},
              {"num_threads", required_argument, 0, 't'},
              {"output_path", required_argument, 0, 'o'},
              {"umpire_mode", required_argument, 0, 'u'},
              {"use_device_aware_mpi", no_argument, 0, 'm'},
              {"use_cuda_sweep", no_argument, 0, 'e'},
              {"use_gpu_kernels", no_argument, 0, 'g'},
              {"use_new_gta", no_argument, 0, 'n'},
              {"verbose", required_argument, 0, 'v'},
#if defined(TETON_ENABLE_MFEM)
              {"serial_refinement_levels", required_argument, 0, 'r'},
              {"parallel_refinement_levels", required_argument, 0, 'z'},
              {"num_Polar", required_argument, 0, 'P'},
              {"num_Azimuth", required_argument, 0, 'A'},
              {"num_Groups", required_argument, 0, 'G'},
#endif
              {0, 0, 0, 0} };

      /* getopt_long stores the option index here. */
      int option_index = 0;

#if defined(TETON_ENABLE_MFEM)
      auto optString = "A:G:P:c:defghi:l:mno:p:r:s:t:u:v:z:";
#else
      auto optString = "c:defghi:l:mno:p:s:t:u:v:";
#endif

      opt = getopt_long(argc, argv, optString, long_options, &option_index);

      /* Detect the end of the options. */
      if (opt == EOF)
      {
         break;
      }

      switch (opt)
      {
         case 'c':
            cycles = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: cycles to execute: " << cycles << std::endl;
            break;
         case 'd':
            dumpCopyOfInput - true;
            break;
         case 'e':
            useCUDASweep = true;
            if (myRank == 0)
               std::cout << "Using experimental streaming CUDA sweep." << std::endl;
            break;
         case 'f':
            freeOMPResources = true;
            if (myRank == 0)
               std::cout
                   << "Enabling freeing of all OpenMP resources on device between cycles via omp_pause_resource_all call."
                   << std::endl;
            break;
         case 'g':
            useGPU = true;
            break;
         case 'h':
            if (myRank == 0)
            {
               std::cout << "Usage: " << argv[0] << "[OPTIONS]" << std::endl;
               std::cout << " -c -num_cycles <cycles>      Number of cycles to execute." << std::endl;
               std::cout << " -d -dump_input_copy          Dump a copy of the problem input to conduit json file."
                         << std::endl;
               std::cout
                   << " -e -use_cuda_sweep           Use experimental CUDA sweep.  Do not specify this option and -g at the same time."
                   << std::endl;
               std::cout
                   << " -f -free_omp_resources       Free all OpenMP runtime resources from device after each problem cycle."
                   << std::endl;
               std::cout
                   << " -g -use-gpu                  Run solvers on GPU and enable associated sub-options, where supported."
                   << std::endl;
               std::cout << " -h -help                     Print this help and exit." << std::endl;
               std::cout << " -i -input_path <path>        Path to input files." << std::endl;
               std::cout
                   << " -l -apply_label              Label this run.  This label will be used to identify this run in any caliper reports."
                   << std::endl;
               std::cout << " -m -use_device_aware_mpi     Use device-aware MPI for GPU runs." << std::endl;
               std::cout << " -n -use_new_gta              Use newer GTA solver." << std::endl;
               std::cout
                   << " -o -output_path <path>       Path to generate output files.  If not set, will disable output files."
                   << std::endl;
#if defined(TETON_ENABLE_CALIPER)
               std::cout
                   << " -p -caliper <string>         Caliper configuration profile.  Use '-p help' to get supported keywords.  'None' disabled caliper.  Default is 'runtime-report'."
                   << std::endl;
#endif
               std::cout << " -s -num_phase_space_sets <num_sets>  Number of phase-angle sets to construct."
                         << std::endl;
               std::cout << " -t -num_threads <threads>    Max number of threads for cpu OpenMP parallel regions."
                         << std::endl;
               std::cout
                   << " -u -umpire_mode <0,1,2>      0 - Disable umpire.  1 - Use Umpire for CPU allocations.  2 - Use Umpire for CPU and GPU allocations."
                   << std::endl;
               std::cout << " -v -verbose [0,1,2]    0 - quite  1 - informational  2 - really chatty (default)"
                         << std::endl;
#if defined(TETON_ENABLE_MFEM)
               std::cout
                   << " -r -serial_refinement_levels <int> Number of times to refine the MFEM mesh before doing parallel decomposition"
                   << std::endl;
               std::cout
                   << " -z -parallel_refinement_levels <int> Number of times to refine the MFEM mesh after doing parallel decomposition"
                   << std::endl;
               std::cout << " -A -num_Azimuthal <int>      Number azimuthal angles in an octant" << std::endl;
               std::cout << " -P -num_Polar <int>          Number polar angles in an octant" << std::endl;
               std::cout << " -G -num_Groups <int>         Number energy groups" << std::endl;
#endif
            }
            return (0);
         case 'i':
            inputPath = std::string(optarg);
            break;
         case 'l':
            label = std::string(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: this run will be identified as '" << label
                         << "' in any caliper spot reports." << std::endl;
            break;
         case 'm':
            useDeviceAwareMPI = true;
            break;
         case 'n':
            useNewGTASolver = true;
            if (myRank == 0)
               std::cout << "Teton driver: using new gta solver." << std::endl;
            break;
         case 'o':
            outputPath = std::string(optarg);
            break;
#if defined(TETON_ENABLE_CALIPER)
         case 'p':
            caliper_config = std::string(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: using caliper configuration: " << caliper_config << std::endl;
            if (caliper_config == "help")
            {
               std::cout << std::endl << "--- AVAILABLE CALIPER KEYWORDS ---" << std::endl;
               for (auto str : configs)
               {
                  std::cout << mgr.get_documentation_for_spec(str.c_str()) << std::endl;
               }
               std::cout << std::endl;
               std::cout << std::endl << "----------------------------------" << std::endl;
               return (0);
            }
#endif
         case 's':
            numPhaseAngleSets = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: number of phase-angle sets to create: " << numPhaseAngleSets << std::endl;
            break;
#if defined(TETON_ENABLE_MFEM)
         case 'r':
            numSerialRefinementLevels = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: number of serial refinement levels: " << numSerialRefinementLevels
                         << std::endl;
            break;
         case 'z':
            numParallelRefinementLevels = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: number of parallel refinement levels: " << numParallelRefinementLevels
                         << std::endl;
            break;
         case 'A':
            numAzimuthal = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: number of azimuthal angles: " << numAzimuthal << std::endl;
            break;
         case 'P':
            numPolar = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: number of polar angles: " << numPolar << std::endl;
            break;
         case 'G':
            numGroups = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: number of energy groups: " << numGroups << std::endl;
            break;
#endif
         case 't':
            numOmpMaxThreads = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: setting max # cpu threads to " << numOmpMaxThreads << std::endl;
            break;
         case 'u':
            useUmpire = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: setting useUmpire to " << useUmpire << std::endl;
            break;
         case 'v':
            verbose = atoi(optarg);
            if (myRank == 0)
               std::cout << "Teton driver: setting verbosity to " << verbose << std::endl;
            break;
         case '?':
            if (myRank == 0)
               std::cout << "Incorrect arguments, try -h to see help." << std::endl;
            break;
      }
   }

//==========================================================
// Start caliper
//==========================================================
#if defined(TETON_ENABLE_CALIPER)
   if (caliper_config != "none")
   {
      mgr.add(caliper_config.c_str());
      if (mgr.error())
      {
         std::cout << "Teton driver: Caliper config error: " << mgr.error_msg() << std::endl;
      }
      mgr.start();
   }
   if (!label.empty())
   {
      adiak::value("Label", label, adiak_general);
   }
#endif

#if defined(TETON_ENABLE_OPENMP)
   if (numOmpMaxThreads == -1)
   {
      numOmpMaxThreads = omp_get_max_threads();
   }
   if (myRank == 0)
      std::cout << "Teton driver: Threading enabled, max number of threads is " << numOmpMaxThreads << std::endl;
#endif

//==========================================================
// Initialize environment on GPU
//==========================================================
#if defined(TETON_ENABLE_OPENMP_OFFLOAD)

   if (myRank == 0)
      print_gpu_mem("Teton driver: Before hello world gpu kernel run.");

// It's necessary to run a small GPU kernel to initialize the GPU state so our timers get accurate benchmarks later.
#pragma omp target
   {
      printf("Teton driver: Hello World! GPU is now initialized.\n");
   }

   if (myRank == 0)
      print_gpu_mem("Teton driver: After hello world gpu kernel run.");

   if (freeOMPResources)
   {
      omp_pause_resource_all(omp_pause_hard);
      if (myRank == 0)
         print_gpu_mem("Teton driver: After omp_pause_resource_all.");
   }

#endif

//==========================================================
// Set up signal handler
// If compiling with IBM XL, use XLF's trce function to emit a code stack trace if a TRAP signal is caught.  This can be used to
// catch errors in any OpenMP kernels by setting'XLSMPOPTS=MSG_TRAP' in your environment.
//==========================================================
#if defined(__ibmxl__)
   struct sigaction sa;
   sa.sa_flags = SA_SIGINFO | SA_RESTART;
   sa.sa_sigaction = xl__trce;
   sigemptyset(&sa.sa_mask);
   sigaction(SIGTRAP, &sa, NULL);
   sigaction(SIGFPE, &sa, NULL);
#endif

   //==========================================================
   // Read in conduit nodes or mfem mesh with problem input
   //==========================================================
   ::Teton::Teton myTetonObject;
   conduit::Node &datastore = myTetonObject.getDatastore();
   conduit::Node &meshBlueprint = myTetonObject.getMeshBlueprint();
   conduit::Node &checkpoint = myTetonObject.getCheckpoint();

   //==========================================================
   // Check if inputPath is path to mfem mesh file.
   // This functionality is currently constrained to a hard-coded
   // crooked pipe problem.  It assumes a mfem mesh is provided
   // with two materials specific boundary conditions.
   //==========================================================
   if (endsWith(inputPath, ".mesh"))
   {
      if (access(inputPath.c_str(), F_OK) != -1)
      {
#if defined(TETON_ENABLE_MFEM)
         if (myRank == 0)
            std::cout << "Teton driver: converting mfem mesh..." << std::endl;

         mfem::Mesh *mesh = new mfem::Mesh(inputPath.c_str(), 1, 1);
         for (int l = 0; l < numSerialRefinementLevels; ++l)
         {
            mesh->UniformRefinement();
         }
         mfem::ParMesh *pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *mesh);
         for (int l = 0; l < numParallelRefinementLevels; ++l)
         {
            mesh->UniformRefinement();
         }

         pmesh->PrintCharacteristics();

         // Set up problem for crooked pipe mesh.  These are all hard-coded currently.
         // TODO: Allow settings these from an input file of some kind, where a value could be set per
         // material ( pmesh tag id ).
         int nelem = mesh->GetNE();

         // TODO: Allow number of group to be set by command line argument.
         int numGroups = 1;
         mfem::Vector gr_bounds(numGroups + 1);
         const double lowerBound = 1.0e-6;
         const double upperBound = 1.0e2;
         const double upperBoundLog = std::log(upperBound);
         const double lowerBoundLog = std::log(lowerBound);
         const double deltaLog = (upperBoundLog - lowerBoundLog) / numGroups;
         for (int g = 0; g < numGroups; ++g)
         {
            gr_bounds[g] = std::exp(lowerBoundLog + g * deltaLog);
         }
         gr_bounds[numGroups] = upperBound;

         mfem::Vector density(nelem);
         mfem::Vector heat_capacity(nelem);
         mfem::Vector rad_temp(nelem);
         mfem::Vector material_temp(nelem);
         mfem::Vector electron_density(nelem);
         mfem::Vector abs_opacity(numGroups * nelem);
         mfem::Vector scat_opacity(numGroups * nelem);

         for (int i = 0; i < nelem; ++i)
         {
            int attr_no = pmesh->GetAttribute(i);
            electron_density[i] = 4.16100608392217e+24;
            if (attr_no == 1) // crooked pipe, thin material
            {
               density[i] = .01;
               heat_capacity[i] = .001;
               rad_temp[i] = .0005;
               material_temp[i] = 5.0e-5;
               for( int g = 0; g <= numGroups; ++g)
               {
                  abs_opacity[i*numGroups + g] = 0.2;
               }
            }
            else if (attr_no == 2) // crooked pipe, thick material
            {
               density[i] = 10.0;
               heat_capacity[i] = 1.0;
               rad_temp[i] = .5;
               material_temp[i] = .05;
               for( int g = 0; g <= numGroups; ++g)
               {
                  abs_opacity[i*numGroups + g] = 200.0;
               }
            }
            else
            {
               std::cerr << "unknown attribute number on mesh element" << std::endl;
               exit(1);
            }
         }
         scat_opacity = 0.0;

         TetonBlueprint blueprint(myTetonObject);

         // Set boundary conditions. These attribute numbers
         // correspond to the tags in the TopHatMesh.py file and
         // and the attribute numbers in the crooked pipe Nurbs
         // mesh in the data directory.
         std::map<int, int> boundary_id_to_type;
         boundary_id_to_type[14] = 32; // bottom (reflecting Teton ID 32) for bottom wall (attribute no. 14)
         boundary_id_to_type[12] = 32; // top wall
         boundary_id_to_type[13] = 32; // right wall
         boundary_id_to_type[11] = 32; // left wall
         boundary_id_to_type[10] = 34; // source (source Teton ID 34) for left pipe wall (attribute no. 10)
         blueprint.SetBoundaryIDs(boundary_id_to_type);
         // set the temperature source on the source boundary
         double T_val = .3;
         blueprint.SetSourceBoundaryTemperature(T_val);

         // Form conduit nodes
         blueprint.OutputConduitMesh(
             pmesh, density, heat_capacity, rad_temp, material_temp, gr_bounds, abs_opacity, scat_opacity, electron_density);
         datastore.update(blueprint.GetConduitInputNode());
         meshBlueprint.update(blueprint.GetConduitMeshNode());

         datastore["sources/profile1/Values"] = 0.3;
         datastore["sources/profile1/NumTimes"] = 1;
         datastore["sources/profile1/NumValues"] = 1;
         datastore["sources/profile1/Multiplier"] = 1.0;

         if (numAzimuthal > 0)
         {
            datastore["quadrature/nazimu"] = numAzimuthal;
         }
         if (numPolar > 0)
         {
            datastore["quadrature/npolar"] = numPolar;
         }

         delete mesh;
         delete pmesh;
#else
         std::cerr << "Unable to open mfem mesh, test driver was not compiled with TETON_ENABLE_MFEM." << std::endl;
         exit(1);
#endif
      }
      else
      {
         std::cerr << "Couldn't find mfem mesh at " << inputPath << std::endl;
         exit(1);
      }
   }
   // Assume this is an input directory with a set of parameter and mesh files.
   else
   {
      try
      {
         std::string input_file_path_base = inputPath + "/parameters_rank" + std::to_string(myRank);
         std::string input_file_path_full = input_file_path_base + ".hdf5";

         // Check for parameters node file ending in .hdf5.
         if (access(input_file_path_full.c_str(), F_OK) != -1)
         {
            conduit::relay::io::load_merged(input_file_path_full, "hdf5", datastore);
         }
         else
         {
            // Check for parameters node file ending in .conduit_json.
            input_file_path_full = input_file_path_base + ".conduit_json";
            if (access(input_file_path_full.c_str(), F_OK) != -1)
            {
               conduit::relay::io::load_merged(input_file_path_full, "conduit_json", datastore);
            }
            else
            {
               std::cerr << "Couldn't find parameters node at " << inputPath << std::endl;
               exit(1);
            }
         }

         // Check for mesh node file ending in .hdf5.
         input_file_path_base = inputPath + "/mesh_rank" + std::to_string(myRank);
         input_file_path_full = input_file_path_base + ".hdf5";
         if (access(input_file_path_full.c_str(), F_OK) != -1)
         {
            conduit::relay::io::load_merged(input_file_path_full, "hdf5", meshBlueprint);
         }
         else
         {
            // Check for parameters node file ending in .conduit_json.
            input_file_path_full = input_file_path_base + ".conduit_json";
            if (access(input_file_path_full.c_str(), F_OK) != -1)
            {
               conduit::relay::io::load_merged(input_file_path_full, "conduit_json", meshBlueprint);
            }
            else
            {
               std::cerr << "Couldn't find mesh node at " << inputPath << std::endl;
               exit(1);
            }
         }
         if (dumpCopyOfInput)
         {
            conduit::relay::io::save(
                datastore, outputPath + "/parameters_rank" + std::to_string(myRank) + ".conduit_json", "conduit_json");
            conduit::relay::io::save(
                meshBlueprint, outputPath + "/mesh_rank" + std::to_string(myRank) + "hdf5", ".hdf5");
         }
      }
      catch (const std::exception &e)
      {
         std::cout << "Teton driver: Error loading conduit node files at " << inputPath << ": " << e.what();
      }
   }

   datastore["memory_allocator/umpire_host_allocator_id"] = -1;
   datastore["memory_allocator/umpire_device_allocator_id"] = -1;

   datastore["size/useCUDASweep"] = false;

   //==========================================================
   // Set problem options passed in via command line
   //==========================================================

   if (useDeviceAwareMPI == true)
   {
      datastore["mpi/useDeviceAddresses"] = 0;
   }

   bool useDeviceAddrForMPI = false; // Pass device addresses to MPI ( device-aware MPI ).

   if (useNewGTASolver == true)
   {
      datastore["size/useNewGTASolver"] = true;
   }

   if (useCUDASweep == true)
   {
      datastore["size/useCUDASweep"] = true;
      datastore["size/useGPU"] = false;

      if (useGPU == true)
      {
         std::cerr << "Select either the experimental CUDA sweep (useCUDASweep) or useGPU, not both." << std::endl;
         exit(1);
      }
   }
   else if (useGPU == true)
   {
      if (myRank == 0)
         std::cout << "Teton driver: -g arg detected, enabling gpu kernels and associated datastore." << std::endl;

      datastore["size/useGPU"] = useGPU;

      // If using the GPU, enable several sub-options.
      if (useUmpire > 0)
      {
#if defined(TETON_ENABLE_UMPIRE)
         if (myRank == 0)
            std::cout << "Teton driver: Enabling use of Umpire CPU and GPU memory pools..." << std::endl;

         auto &rm = umpire::ResourceManager::getInstance();

         // Create umpire allocators.
         auto host_pinned_pool
             = rm.makeAllocator<umpire::strategy::QuickPool>("HOST_PINNED_QUICK_POOL", rm.getAllocator("PINNED"));
         auto thread_safe_host_pinned_pool = rm.makeAllocator<umpire::strategy::ThreadSafeAllocator>(
             "THREAD_SAFE_PINNED_QUICK_POOL", host_pinned_pool);

         datastore["memory_allocator/umpire_host_allocator_id"] = thread_safe_host_pinned_pool.getId();
         if (useUmpire > 1)
         {
            auto device_pool
                = rm.makeAllocator<umpire::strategy::QuickPool>("DEVICE_QUICK_POOL", rm.getAllocator("DEVICE"));
            auto thread_safe_device_pool
                = rm.makeAllocator<umpire::strategy::ThreadSafeAllocator>("THREAD_SAFE_DEVICE_QUICK_POOL", device_pool);
            datastore["memory_allocator/umpire_device_allocator_id"] = thread_safe_device_pool.getId();
         }
         else
         {
            datastore["memory_allocator/umpire_device_allocator_id"] = -1;
         }
#else
         if (myRank == 0)
         {
            std::cerr
                << "Teton driver: Unable to enable Umpire CPU and GPU memory pools, code was not built with TETON_ENABLE_UMPIRE."
                << std::endl;
         }
         exit(1);
#endif
      }

      // Enable the GPU CUDA Boltzmann Compton solver ( only has an effect if using BC solver).
      datastore["size/useCUDASolver"] = true;
   }

   if (numPhaseAngleSets != 0)
   {
      datastore["quadrature/nSetsMaster"] = numPhaseAngleSets;
   }

   datastore["concurrency/nOmpMaxThreads"] = numOmpMaxThreads;

   //Set verbosity level ( default is 1 )
   datastore["verbose"] = verbose;

   myTetonObject.initialize();
   // Be sure to clear out the blueprint node, since initialize() is finished.
   // Only populate the mesh blueprint node after this point if you want to trigger mesh movement or
   // updates to materials or opacities, etc.
   meshBlueprint.reset();

   for (int cycle = 1; cycle <= cycles; cycle++)
   {
#if defined(TETON_ENABLE_MEMUSAGES)
      double memory_mb_used = MemUsage_get_curr_proc_mbytes();
      if (myRank == 0)
         std::cout << "Teton driver: cycle " << cycle << " CPU process memory (MB): " << memory_mb_used << std::endl;
#endif

      if (!outputPath.empty())
      {
         myTetonObject.dump(cycle, outputPath, "dump");
      }
      if (myRank == 0)
      {
         std::cout << "----------" << std::endl;
         std::cout << "CYCLE " << cycle << std::endl;
         std::cout << "----------" << std::endl;
      }
      double dtrad = datastore["iteration/dtrad"].value();
      myTetonObject.setTimeStep(cycle, dtrad);
      myTetonObject.step(cycle);

#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
      if (freeOMPResources)
      {
         omp_pause_resource_all(omp_pause_hard);
         if (myRank == 0)
            print_gpu_mem("Teton driver: After omp_pause_resource_all.");
      }
#endif
   }

   std::cout << "=================================================================\n";
   std::cout << "=================================================================\n";
   std::cout << "Teton finished cycling\n";
   std::cout << "=================================================================\n";
   std::cout << "=================================================================\n";

#if defined(TETON_ENABLE_UMPIRE)
   if (useGPU == 1 && useUmpire > 0)
   {
      if (myRank == 0)
         std::cout << "Teton driver: Deleting Umpire CPU and GPU memory pools..." << std::endl;

      auto &rm = umpire::ResourceManager::getInstance();

      // Release memory from umpire allocators
      auto thread_safe_host_pinned_pool = rm.getAllocator("THREAD_SAFE_PINNED_QUICK_POOL");
      if (myRank == 0)
         print_bytes_as_gb("Teton driver: Thread safe host pinned pool size: ",
                           thread_safe_host_pinned_pool.getActualSize());
      thread_safe_host_pinned_pool.release();
      auto host_pinned_pool = rm.getAllocator("HOST_PINNED_QUICK_POOL");
      if (myRank == 0)
         print_bytes_as_gb("Teton driver: Host pinned (parent) pool size: ", host_pinned_pool.getActualSize());
      host_pinned_pool.release();

      if (useUmpire > 1)
      {
         auto thread_safe_device_pool = rm.getAllocator("THREAD_SAFE_DEVICE_QUICK_POOL");
         if (myRank == 0)
            print_bytes_as_gb("Teton driver: Thread safe device pool size: ", thread_safe_device_pool.getActualSize());
         thread_safe_device_pool.release();
         auto device_pool = rm.getAllocator("DEVICE_QUICK_POOL");
         if (myRank == 0)
            print_bytes_as_gb("Teton driver: Device (parent) pool size: ", device_pool.getActualSize());
         device_pool.release();
      }
      if (myRank == 0)
         print_gpu_mem("Teton driver: After Umpire device pool release.");
   }
#endif

#if defined(TETON_ENABLE_OPENMP_OFFLOAD)
   omp_pause_resource_all(omp_pause_hard);
   if (myRank == 0)
      print_gpu_mem("Teton driver: After omp_pause_resource_all at end of run.");
#endif

#if defined(TETON_ENABLE_CALIPER)
   adiak::fini();
   mgr.flush();
#endif

   MPI_Finalize();
   return 0;
}
