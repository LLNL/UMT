#include "mpi.h"
#include <getopt.h>
#include <malloc.h>
#include <sstream>
#include <string>
#include <sys/stat.h> //for mkdir
#include <unistd.h>   //for getopt, access
#include <stdio.h> // for getCurrentRSS and getPeakRSS

// XLF signal handler function to emit a stack trace.
#if defined(__ibmxl__)
#include "signal.h"
extern "C" void xl__trce(int, siginfo_t *, void *);
#endif

#include <signal.h>
#if defined(__linux__)
#include <fenv.h>
#endif

#if defined(TETON_ENABLE_OPENMP)
#include "omp.h"
#endif

#if defined(TETON_ENABLE_UMPIRE)
#include "umpire/Umpire.hpp"
#include "umpire/strategy/QuickPool.hpp"
#include "umpire/strategy/ThreadSafeAllocator.hpp"
#endif

#if defined(TETON_ENABLE_CUDA)
#include "cuda_runtime_api.h"
#endif

#include "TetonConduitInterface.hh"
#include "TetonInterface.hh"

#include "conduit/conduit.hpp"
#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_blueprint_mesh_utils.hpp"
#include "conduit/conduit_relay.hpp"

#if defined(TETON_ENABLE_MFEM)
#include "mfem.hpp"
#endif

#if defined(TETON_ENABLE_CALIPER)
#include "adiak.hpp"
#include "caliper/cali-manager.h"
#include "caliper/cali-mpi.h"
#include "caliper/cali.h"
#else
#define CALI_MARK_BEGIN(label)
#define CALI_MARK_END(label)
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
#else
   (void) label;
#endif
   return (0);
}

// Returns the current resident set size (physical memory use) measured in kbytes, or zero if the value cannot be determined on this OS.
size_t getCurrentRSS( )
{
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);
}

int get_used_mem_in_bytes()
{
   struct mallinfo minfo = mallinfo();

   // uordblks - Total allocated space (bytes)
   // hblkhd   - Space allocated in mmapped regions
   return (minfo.uordblks + minfo.hblkhd);
}

//==========================================================
//  run syntax:
//
// ./test_driver -c <# cycles> -i <input path> -o <output path>
// Note: to use the new blueprint format, ./test_driver -b ...
//==========================================================

int main(int argc, char *argv[])
{
   int myRank = 0;
   int mySize = 0;
   int opt;
   unsigned int cycles = 0;
   int numPhaseAngleSets = 0;
   int useUmpire = 1;
   int numOmpMaxThreads = -1; // Max number of CPU threads to use  If -1, use value from omp_get_max_threads()
   double fixedDT = 0.0;
   bool dumpViz = false;

#if defined(TETON_ENABLE_MFEM)
   int numSerialRefinementFactor = 2;
   int numParallelRefinementFactor = 2;
   int numSerialRefinementLevels = 0;
   int numParallelRefinementLevels = 0;
   int numPolar = -1;
   int numAzimuthal = -1;
   int numGroups = -1;
   int benchmarkProblem = 0;
   mfem::ParMesh *pmesh = nullptr;
   mfem::ConduitDataCollection *conduit_data_collec = nullptr;
#endif

   // MPI
   MPI_Comm comm = MPI_COMM_WORLD;
   int provided = 0;
   int claimed = 0;
   int request = MPI_THREAD_SINGLE;
   int verbose = 1;

   bool useGPU = false;
   bool useCUDASweep = false;
   bool useNewGTASolver = false;
   bool useDeviceAwareMPI = false; // Pass device addresses to MPI ( device-aware MPI ).

   std::string inputPath(".");
   std::string outputPath(".");
   std::string label("");
   std::string colorFile("");
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

   CALI_MARK_BEGIN("Teton_Test_Driver");

   if (MPI_Init_thread(&argc, &argv, request, &provided) != MPI_SUCCESS)
   {
      std::cerr << "Error calling MPI_Init_thread. " << std::endl;
      exit(1);
   }

   MPI_Query_thread(&claimed);

   if (provided < request)
   {
      std::cerr << "Teton driver: MPI_Init_thread was only able to provided thread level support " << provided
                << ".  Teton requested level " << request << std::endl;
      exit(1);
   }
   if (claimed < request)
   {
      std::cerr << "Teton driver: MPI_Init_thread was only able to provided thread level support " << claimed
                << ".  Teton requested level " << request << std::endl;
      exit(1);
   }

#ifdef SIGSEGV
   signal(SIGSEGV, SIG_DFL);
#endif
#ifdef SIGILL
   signal(SIGILL, SIG_DFL);
#endif
#ifdef SIGFPE
   signal(SIGFPE, SIG_DFL);
#endif
#ifdef SIGINT
   signal(SIGINT, SIG_DFL);
#endif
#ifdef SIGABRT
   signal(SIGABRT, SIG_DFL);
#endif
#ifdef SIGTERM
   signal(SIGTERM, SIG_DFL);
#endif
#ifdef SIGQUIT
   signal(SIGQUIT, SIG_DFL);
#endif
#if defined(__linux__)
   // This is in here for supporting Linux's floating point exceptions.
   feenableexcept(FE_DIVBYZERO);
   feenableexcept(FE_INVALID);
   feenableexcept(FE_OVERFLOW);
#endif

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
   adiak::value("Version", teton_get_version(), adiak_general, "TetonBuildInfo");
   adiak::value("SHA1", teton_get_git_sha1(), adiak_general, "TetonBuildInfo");
   adiak::value("CxxCompiler", teton_get_cxx_compiler(), adiak_general, "TetonBuildInfo");
   adiak::value("FortranCompiler", teton_get_fortran_compiler(), adiak_general, "TetonBuildInfo");
#endif

   MPI_Comm_rank(comm, &myRank);
   MPI_Comm_size(comm, &mySize);

   if (myRank == 0)
   {
      std::cout << "Teton driver: number of MPI ranks: " << mySize << std::endl;
   }

   //==========================================================
   // Get command line arguments
   //==========================================================

   while (1)
   {
      static struct option long_options[] = {
         {"apply_label", no_argument, 0, 'l'},
         {"benchmark_problem", required_argument, 0, 'b'},
         {"caliper", required_argument, 0, 'p'},
         {"help", no_argument, 0, 'h'},
         {"input_path", required_argument, 0, 'i'},
         {"num_cycles", required_argument, 0, 'c'},
         {"dt", required_argument, 0, 'D'},
         {"num_phase_space_sets", required_argument, 0, 's'},
         {"num_threads", required_argument, 0, 't'},
         {"output_path", required_argument, 0, 'o'},
         {"umpire_mode", required_argument, 0, 'u'},
         {"use_device_aware_mpi", no_argument, 0, 'm'},
         {"use_cuda_sweep", no_argument, 0, 'e'},
         {"use_gpu_kernels", no_argument, 0, 'g'},
         {"use_new_gta", no_argument, 0, 'n'},
         {"verbose", required_argument, 0, 'v'},
         {"write_viz_file", no_argument, 0, 'V'},
#if defined(TETON_ENABLE_MFEM)
         {"serial_refinement_levels", required_argument, 0, 'r'},
         {"parallel_refinement_levels", required_argument, 0, 'z'},
         {"serial_refinement_factor", required_argument, 0, 'R'},
         {"parallel_refinement_factor", required_argument, 0, 'Z'},
         {"num_Polar", required_argument, 0, 'P'},
         {"num_Azimuthal", required_argument, 0, 'A'},
         {"num_Groups", required_argument, 0, 'G'},
         {"color_file", required_argument, 0, 'C'},
#endif
         {0, 0, 0, 0}
      };

      /* getopt_long stores the option index here. */
      int option_index = 0;

#if defined(TETON_ENABLE_MFEM)
      auto optString = "A:b:G:P:C:c:D:eghi:l:mno:p:r:R:s:t:u:Vv:z:Z:";
#else
      auto optString = "c:D:eghi:l:mno:p:s:t:u:Vv:";
#endif

      opt = getopt_long(argc, argv, optString, long_options, &option_index);

      /* Detect the end of the options. */
      if (opt == EOF)
      {
         break;
      }

      switch (opt)
      {
         case 'b':
            benchmarkProblem = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: Running predefined benchmark problem UMT SP#" << benchmarkProblem
                         << std::endl;
            }
            break;
         case 'c':
            cycles = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: cycles to execute: " << cycles << std::endl;
            }
            break;
         case 'D':
            fixedDT = atof(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: fixed dt selected: " << fixedDT << std::endl;
            }
            break;
         case 'e':
            useCUDASweep = true;
            if (myRank == 0)
            {
               std::cout << "Using experimental streaming CUDA sweep." << std::endl;
            }
            break;
         case 'g':
            useGPU = true;
            break;
         case 'h':
            if (myRank == 0)
            {
               std::cout << "Usage: " << argv[0] << "[OPTIONS]" << std::endl;
               std::cout
                  << " -b, --benchmark_problem <0,1,2>      Run predefined UMT benchmark problem.  0 = user specified # angles and # groups."
                  << std::endl;
               std::cout << " -c, --num_cycles <cycles>      Number of cycles to execute." << std::endl;
               std::cout
                  << " -e, --use_cuda_sweep           Use experimental CUDA sweep.  Do not specify this option and -g at the same time."
                  << std::endl;
               std::cout
                  << " -g, --use-gpu                  Run solvers on GPU and enable associated sub-options, where supported."
                  << std::endl;
               std::cout << " -h, --help                     Print this help and exit." << std::endl;
               std::cout << " -i, --input_path <path>        Path to input files." << std::endl;
               std::cout
                  << " -l, --apply_label              Label this run.  This label will be used to identify this run in any caliper reports."
                  << std::endl;
               std::cout << " -m --use_device_aware_mpi     Use device-aware MPI for GPU runs." << std::endl;
               std::cout << " -n --use_new_gta              Use newer GTA solver." << std::endl;
               std::cout
                  << " -o, --output_path <path>       Path to generate output files.  If not set, will disable output files."
                  << std::endl;
#if defined(TETON_ENABLE_CALIPER)
               std::cout << " -p, --caliper <string>         Caliper configuration profile.  Set <string> to 'help'"
                         << " to get supported keywords.  'None' disabled caliper.  Default is 'runtime-report'."
                         << std::endl;
#endif
               std::cout << " -s, --num_phase_space_sets <num_sets>  Number of phase-angle sets to construct."
                         << std::endl;
               std::cout << " -t, --num_threads <threads>    Max number of threads for cpu OpenMP parallel regions."
                         << std::endl;
               std::cout << " -u, --umpire_mode <0,1,2>      0 - Disable umpire.  1 - Use Umpire for CPU allocations."
                         << "  2 - Use Umpire for CPU and GPU allocations." << std::endl;
               std::cout
                  << " -v, --verbose [0,1,2]    0 - quite  1 - informational(default)  2 - really chatty and dump files"
                  << std::endl;
#if defined(TETON_ENABLE_MFEM)
               std::cout
                  << " -r, --serial_refinement_levels <int> Number of times to halve each edge the MFEM mesh before "
                  << "doing parallel decomposition.  Applied after the refinement_factor. (factor of 2^((r*dim) new zones)"
                  << std::endl;
               std::cout
                  << " -z, --parallel_refinement_levels <int> Number of times to halve each edge the MFEM mesh after "
                  << "doing parallel decomposition.  Applied after the refinement_factor. (factor of 2^(z*dim) new zones)"
                  << std::endl;
               std::cout
                  << " -R, --serial_refinement_factor <int> Number of subdivisions for each edge in the original MFEM "
                  << "mesh before doing parallel decomposition. (factor of (R+1)^dim new zones)" << std::endl;
               std::cout
                  << " -Z, --parallel_refinement_factor <int> Number of subdivisions for each edge in the original MFEM "
                  << "mesh after doing parallel decomposition. (factor of (Z+1)^dim new zones)" << std::endl;
               std::cout << " -A, --num_Azimuthal <int>      Number azimuthal angles in an octant" << std::endl;
               std::cout << " -P, --num_Polar <int>          Number polar angles in an octant" << std::endl;
               std::cout << " -G, --num_Groups <int>         Number energy groups" << std::endl;
               std::cout << " -C, --color_file <string>      color file for manual decomposition" << std::endl;
#endif
               std::cout << " -V, --write-viz-file     Output blueprint mesh vizualization file each cycle"
                         << std::endl;
            }
            return (0);
         case 'i':
            inputPath = std::string(optarg);
            break;
         case 'l':
            label = std::string(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: this run will be identified as '" << label
                         << "' in any caliper spot reports." << std::endl;
            }
            break;
         case 'm':
            useDeviceAwareMPI = true;
            break;
         case 'n':
            useNewGTASolver = true;
            if (myRank == 0)
            {
               std::cout << "Teton driver: using new gta solver." << std::endl;
            }
            break;
         case 'o':
            outputPath = std::string(optarg);
            break;
#if defined(TETON_ENABLE_CALIPER)
         case 'p':
            caliper_config = std::string(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: using caliper configuration: " << caliper_config << std::endl;
            }
            if (caliper_config == "help")
            {
               if (myRank == 0)
               {
                  std::cout << std::endl << "--- AVAILABLE CALIPER KEYWORDS ---" << std::endl;
                  for (auto str : configs)
                  {
                     std::cout << mgr.get_documentation_for_spec(str.c_str()) << std::endl;
                  }
                  std::cout << std::endl;
                  std::cout << std::endl << "----------------------------------" << std::endl;
               }
               return (0);
            }
#endif
         case 's':
            numPhaseAngleSets = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of phase-angle sets to create: " << numPhaseAngleSets << std::endl;
            }
            break;
#if defined(TETON_ENABLE_MFEM)
         case 'r':
            numSerialRefinementLevels = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of serial refinement levels: " << numSerialRefinementLevels
                         << std::endl;
            }
            break;
         case 'z':
            numParallelRefinementLevels = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of parallel refinement levels: " << numParallelRefinementLevels
                         << std::endl;
            }
            break;
         case 'R':
            numSerialRefinementFactor = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: serial refinement factor: " << numSerialRefinementFactor << std::endl;
            }
            break;
         case 'Z':
            numParallelRefinementFactor = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: parallel refinement factor: " << numParallelRefinementFactor << std::endl;
            }
            break;
         case 'A':
            numAzimuthal = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of azimuthal angles: " << numAzimuthal << std::endl;
            }
            break;
         case 'P':
            numPolar = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of polar angles: " << numPolar << std::endl;
            }
            break;
         case 'G':
            numGroups = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of energy groups: " << numGroups << std::endl;
            }
            break;
         case 'C':
            colorFile = std::string(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: Using color file for decomposition: " << colorFile << std::endl;
            }
            break;
#endif
         case 't':
            numOmpMaxThreads = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: setting max # cpu threads to " << numOmpMaxThreads << std::endl;
            }
            break;
         case 'u':
            useUmpire = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: setting useUmpire to " << useUmpire << std::endl;
            }
            break;
         case 'V':
            dumpViz = true;
            if (myRank == 0)
            {
               std::cout << "Teton driver: output mesh blueprint visualization file each cycle." << std::endl;
            }
            break;
         case 'v':
            verbose = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: setting verbosity to " << verbose << std::endl;
            }
            break;
         case '?':
            if (myRank == 0)
            {
               std::cout << "Incorrect arguments, try -h to see help." << std::endl;
            }
            break;
      }
   }

   //==========================================================
   // If benchmark problem specified, set the parameters.
   // UMT SP #1
   // 3 x 3 product quadrature
   // 128 groups
   //
   // UMT SP #2
   // 2 x 2 product quadrature
   // 16 groups
   //
   // This still supports overriding each value with individual
   // command line arguments for angles and # groups.
   //==========================================================
   if (benchmarkProblem > 0)
   {
      if (cycles < 1)
      {
         cycles = 10;
      }

      if (benchmarkProblem == 1)
      {
         numPolar = 3;
         numAzimuthal = 3;
         numGroups = 128;
      }
      else if (benchmarkProblem == 2)
      {
         numPolar = 2;
         numAzimuthal = 2;
         numGroups = 16;
      }
      else
      {
         std::cerr << "Teton driver: Unknown benchmark problem #" << benchmarkProblem << std::endl;
      }

      if (label.empty())
      {
         label = "UMTSP#" + std::to_string(benchmarkProblem);
         if (useGPU)
         {
            label += "_GPU";
         }
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
         if (myRank == 0)
         {
            std::cout << "Teton driver: Caliper config error: " << mgr.error_msg() << std::endl;
         }
      }
      mgr.start();
   }
   if (!label.empty())
   {
      adiak::value("ProblemName", label, adiak_general);
   }
#endif

#if defined(TETON_ENABLE_OPENMP)
   if (numOmpMaxThreads == -1)
   {
      numOmpMaxThreads = omp_get_max_threads();
   }
   if (myRank == 0)
   {
      std::cout << "Teton driver: Threading enabled, max number of threads is " << numOmpMaxThreads << std::endl;
   }
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
   conduit::Node &options = myTetonObject.getOptions();
   conduit::Node &meshBlueprint = myTetonObject.getMeshBlueprint();

   //==========================================================
   // Read in mesh from an mfem mesh file.  We set up a uniform
   // temperature problem with simple source boundary conditions.
   //
   // TODO - All this hard-coding can be moved into an input file that lives
   // alongside the mfem mesh file.
   // TODO - All this code for converting a mfem mesh and input to a blueprint
   // mesh should be moved to another function in another source file, so the
   // driver doesn't have all this here. -- black27
   //==========================================================
   //==========================================================

   if (endsWith(inputPath, ".mesh"))
   {
      CALI_MARK_BEGIN("Teton_Read_Mfem_Input");
      if (access(inputPath.c_str(), F_OK) != -1)
      {
#if defined(TETON_ENABLE_MFEM)
         if (myRank == 0)
         {
            std::cout << "Teton driver: reading mfem mesh: " << inputPath << std::endl;
         }

         mfem::Mesh mesh(inputPath.c_str(), 1, 1);

         CALI_MARK_BEGIN("Teton_Refine_Serial_Mesh");
         for (int l = 0; l < numSerialRefinementLevels; ++l)
         {
            if (myRank == 0)
            {
               std::cout << "Teton driver: Uniformly refining serial mesh, iteration " << l + 1
                         << ", factor = " << numSerialRefinementFactor << std::endl;
            }
            if (numSerialRefinementFactor == 2)
            {
               mesh.UniformRefinement();
            }
            else
            {
               int ref_type = mfem::BasisType::ClosedUniform;
               mesh = mfem::Mesh::MakeRefined(mesh, numSerialRefinementFactor, ref_type);
            }
         }
         CALI_MARK_END("Teton_Refine_Serial_Mesh");

         if (benchmarkProblem > -1 && (numSerialRefinementLevels > 0))
         {
            // Save refined mesh to file, if running benchmark problems ( for later use ).
            std::ofstream mesh_ofs("refined_mesh.mesh");
            mesh_ofs.precision(16);
            mesh.Print(mesh_ofs);
            mesh_ofs.close();
         }

         CALI_MARK_BEGIN("Teton_Create_Par_Mesh");
         // MFEM does not support parallel refinement on NURBs meshes.
         // Convert to high order mesh to enable basic refinement capability.
         // Also add a grid function if we didn't have one before.
         mesh.SetCurvature(1);

         if (myRank == 0)
         {
            std::cout << "Teton driver: decomposing serial mesh into parallel" << std::endl;
         }

         if (colorFile.size() > 0)
         {
            if (access(colorFile.c_str(), F_OK) != -1)
            {
               int nelem = mesh.GetNE();
               std::vector<int> colorData;
               colorData.reserve(nelem);

               int e = 0;
               std::ifstream colorFileStream(colorFile.c_str());
               while (!colorFileStream.eof() or e == nelem)
               {
                  int c = 0;
                  colorFileStream >> c;
                  colorData.push_back(c);
                  ++e;
               }
               if (e < nelem)
               {
                  std::cerr << "Not enough colors in " << colorFile << std::endl;
                  return (1);
               }
               colorFileStream.close();

               pmesh = new mfem::ParMesh(comm, mesh, colorData.data());
            }
            else
            {
               std::cerr << "could not open color file: " << colorFile << std::endl;
               return (1);
            }
         }
         else
         {
            // Only re-order without the color file, otherwise we won't know
            // order the elements are.

            // Sort the grid for better locality
            mfem::Array<int> ordering;
            mesh.GetHilbertElementOrdering(ordering);
            mesh.ReorderElements(ordering);

            // TODO Make optional.  This will use the space-filling curve for
            // partitioning.  It's both nicer and more horrific than Metis,
            // if that was even possible.
            // mesh.EnsureNCMesh();

            pmesh = new mfem::ParMesh(comm, mesh);
            if (int wrong = pmesh->CheckElementOrientation(true) > 0)
            {
               std::cout << "There were " << wrong
                         << " 3D mesh elements with the wrong orientation after reordering.\n";
            }
            if (int wrong = pmesh->CheckBdrElementOrientation(true) > 0)
            {
               std::cout << "There were " << wrong
                         << " 3D mesh boundary elements with the wrong orientation after reordering.\n";
            }
         }
         CALI_MARK_END("Teton_Create_Par_Mesh");

         CALI_MARK_BEGIN("Teton_Refine_Par_Mesh");
         for (int l = 0; l < numParallelRefinementLevels; ++l)
         {
            if (myRank == 0)
            {
               std::cout << "Teton driver: Uniformly refining parallel mesh, iteration " << l + 1
                         << ", factor = " << numParallelRefinementFactor << std::endl;
            }
            if (numParallelRefinementFactor == 2)
            {
               mesh.UniformRefinement();
            }
            else
            {
               int ref_type = mfem::BasisType::ClosedUniform;
               *pmesh = mfem::ParMesh::MakeRefined(*pmesh, numParallelRefinementFactor, ref_type);
            }
         }
         CALI_MARK_END("Teton_Refine_Par_Mesh");

         if (myRank == 0)
         {
            std::cout << "Teton driver: Final parallel mesh characteristics are:" << std::endl;
         }
         pmesh->PrintInfo();

         // This is local number of elements.
         int nelem = pmesh->GetNE();

         CALI_MARK_BEGIN("Teton_Create_BP_Mesh");
         // Create a blueprint node from the mfem mesh
         conduit_data_collec = new mfem::ConduitDataCollection("mfem_conduit_data_collection", pmesh);
         // Note - the mesh blueprint node contains pointers back into the mfem
         // mesh for some of the data.  For example, the coordinates.
         // ***************************
         // DO NOT DELETE the mfem mesh objects until this blueprint node is no
         // longer needed.
         // ***************************
         conduit_data_collec->MeshToBlueprintMesh(pmesh, meshBlueprint);

         CALI_MARK_END("Teton_Create_BP_Mesh");
         CALI_MARK_BEGIN("Teton_Init_BP_Fields");
         // Delete extra fields we don't need.  Some of these are not yet supported
         // by VisIt (https://wci.llnl.gov/simulation/computer-codes/visit)
         if (meshBlueprint.has_path("topologies/main/grid_function"))
         {
            meshBlueprint.remove("topologies/main/grid_function");
         }
         if (meshBlueprint.has_path("topologies/main/boundary_topology"))
         {
            meshBlueprint.remove("topologies/main/boundary_topology");
         }
         if (meshBlueprint.has_path("fields/mesh_nodes"))
         {
            meshBlueprint.remove("fields/mesh_nodes");
         }

         if (numGroups < 1)
         {
            std::cerr
               << "Teton driver: Must specify number of energy groups angles via '-G#' or by specifying a benchmark problem via '-b#'."
               << std::endl;
            exit(1);
         }
         std::vector<double> gr_bounds(numGroups + 1);
         const double lowerBound = 1.0e-6;
         const double upperBound = 1.0e2;
         const double upperBoundLog = std::log(upperBound);
         const double lowerBoundLog = std::log(lowerBound);
         const double deltaLog = (upperBoundLog - lowerBoundLog) / numGroups;
         for (int g = 0; g <= numGroups; ++g)
         {
            gr_bounds[g] = std::exp(lowerBoundLog + g * deltaLog);
         }
         gr_bounds[numGroups] = upperBound;

         //Energy groups and SN quadrature info
         int qtype = 2;
         int qorder = 10;
         if (numPolar < 1)
         {
            std::cerr
               << "Teton driver: Must specify number of polar angles via '-P#' or by specifying a benchmark problem via '-b#'."
               << std::endl;
            exit(1);
         }
         if (numAzimuthal < 1)
         {
            std::cerr
               << "Teton driver: Must specify number of azimuthal angles via '-A#' or by specifying a benchmark problem via '-b#'."
               << std::endl;
            exit(1);
         }
         int paxis = 1;
         options["quadrature/gnu"].set(gr_bounds.data(), gr_bounds.size());
         options["quadrature/qtype"] = qtype;
         options["quadrature/qorder"] = qorder;
         options["quadrature/npolar"] = numPolar;
         options["quadrature/nazimu"] = numAzimuthal;
         options["quadrature/paxis"] = paxis;
         options["quadrature/num_groups"] = numGroups;
         options["quadrature/gtaorder"] = 2;
         options["quadrature/nSetsMaster"] = -1;
         options["quadrature/nSets"] = 1;

         // TODO: Make MFEM Grid functions for ConduitDataCollection to handle
         // instead.

         //Material dependent fields
         std::vector<double> thermo_density(nelem);
         std::vector<double> electron_specific_heat(nelem);
         std::vector<double> radiation_temperature(nelem);
         std::vector<double> electron_temperature(nelem);
         std::vector<double> absorption_opacity(numGroups * nelem);
         std::vector<double> scattering_opacity(numGroups * nelem, 0.0);
         std::vector<double> electron_number_density(nelem, 4.16100608392217e+24);

         // Field value to initialize each material to.
         // Material # -> field -> value.
         std::map<int, std::map<std::string, double>> material_field_vals;

         //NOTE: some of these fields need to have the initial values updated.
         //These values are based on the old field meanings, before they were
         //renamed for the current mesh blueprint interface.  Need to consult
         //with haut3.  -- black27

         // We only support one material on an MFEM mesh
         material_field_vals[1]["thermo_density"] = 1.31;
         material_field_vals[1]["electron_specific_heat"] = 0.501;
         material_field_vals[1]["radiation_temperature"] = 0.05;
         material_field_vals[1]["electron_temperature"] = 0.5;
         material_field_vals[1]["absorption_opacity"] = 10.0;

         for (int i = 0; i < nelem; ++i)
         {
            // Ignore, we only support one.
            //int attr_no = pmesh->GetAttribute(i);
            int attr_no = 1;

            thermo_density[i] = material_field_vals[attr_no]["thermo_density"];
            electron_specific_heat[i] = material_field_vals[attr_no]["electron_specific_heat"];
            radiation_temperature[i] = material_field_vals[attr_no]["radiation_temperature"];
            electron_temperature[i] = material_field_vals[attr_no]["electron_temperature"];

            double abs_opacity = material_field_vals[attr_no]["absorption_opacity"];

            for (int g = 0; g < numGroups; ++g)
            {
               absorption_opacity[i * numGroups + g] = abs_opacity;
            }
         }

         // Store the various fields (density, material temperature, etc.)
         meshBlueprint["fields/thermo_density/association"] = "element";
         meshBlueprint["fields/thermo_density/topology"] = "main";
         meshBlueprint["fields/thermo_density/values"].set(thermo_density.data(), thermo_density.size());

         meshBlueprint["fields/electron_specific_heat/association"] = "element";
         meshBlueprint["fields/electron_specific_heat/topology"] = "main";
         meshBlueprint["fields/electron_specific_heat/values"].set(electron_specific_heat.data(),
                                                                   electron_specific_heat.size());

         meshBlueprint["fields/electron_temperature/association"] = "element";
         meshBlueprint["fields/electron_temperature/topology"] = "main";
         meshBlueprint["fields/electron_temperature/values"].set(electron_temperature.data(),
                                                                 electron_temperature.size());

         meshBlueprint["fields/radiation_temperature/association"] = "element";
         meshBlueprint["fields/radiation_temperature/topology"] = "main";
         meshBlueprint["fields/radiation_temperature/values"].set(radiation_temperature.data(),
                                                                  radiation_temperature.size());

         meshBlueprint["fields/absorption_opacity/association"] = "element";
         meshBlueprint["fields/absorption_opacity/topology"] = "main";
         meshBlueprint["fields/absorption_opacity/values"].set(absorption_opacity.data(), absorption_opacity.size());

         meshBlueprint["fields/scattering_opacity/association"] = "element";
         meshBlueprint["fields/scattering_opacity/topology"] = "main";
         meshBlueprint["fields/scattering_opacity/values"].set(scattering_opacity.data(), scattering_opacity.size());

         meshBlueprint["fields/electron_number_density/association"] = "element";
         meshBlueprint["fields/electron_number_density/topology"] = "main";
         meshBlueprint["fields/electron_number_density/values"].set(electron_number_density.data(),
                                                                    electron_number_density.size());

         // Make all boundaries vacuum
         std::map<int, int> boundary_id_to_type;
         for (int i = 0; i < pmesh->bdr_attributes.Size(); ++i)
         {
            int ba = pmesh->bdr_attributes[i];
            // Draco will tag untagged faces with 9999, even if they are
            // interior.  We need to skip those.
            if (ba != 9999)
            {
               int bctype = 35; // vacuum
               // Maybe if we'd run in RZ, we need to make the axis reflecting?
               //if( ba == 1 or ba == 4)
               //{
               //   int bctype = 32; // Reflecting
               //}
               boundary_id_to_type[pmesh->bdr_attributes[i]] = 35; // vacuum
            }
         }

         std::vector<int> keys, values;
         for (std::map<int, int>::iterator it = boundary_id_to_type.begin(); it != boundary_id_to_type.end(); ++it)
         {
            int k = it->first;
            int v = it->second;
            keys.push_back(k);
            // There's some memory error here, which is very, very odd.
            values.push_back(v);
         }
         options["boundary_conditions/id_to_type_map/ids"] = keys;
         options["boundary_conditions/id_to_type_map/types"] = values;

         // Needs to be re-done when blueprint interface for specifying profiles is updated.
         options["sources/profile1/Values"] = 0.3;
         options["sources/profile1/NumTimes"] = 1;
         options["sources/profile1/NumValues"] = 1;
         options["sources/profile1/Multiplier"] = 1.0;

         // Prune the boundary topology of all interior boundary elements.  MFEM creates a boundary topology over
         // both problem surface elements and interior shared boundaries between domains.
         // Teton only wants the surface elements.  Teton uses the adjacency lists to determine shared boundary elements.
         conduit::int_accessor bndry_vals = meshBlueprint.fetch_existing("fields/boundary_attribute/values").value();
         conduit::int_accessor element_points = meshBlueprint
                                                   .fetch_existing("topologies/boundary/elements/connectivity")
                                                   .value();
         std::string element_type = meshBlueprint.fetch_existing("topologies/boundary/elements/shape").as_string();

         size_t num_points_in_element;

         if (element_type == "point")
         {
            num_points_in_element = 1;
         }
         else if (element_type == "line")
         {
            num_points_in_element = 2;
         }
         else if (element_type == "quad")
         {
            num_points_in_element = 4;
         }
         else
         {
            std::cerr << "Unsupported element type of: " << element_type << std::endl;
            return (1);
         }

         std::vector<int> new_bnd_attribs;
         std::vector<int> new_connectivity;
         size_t current_element = 0;
         for (size_t i = 0; i < bndry_vals.number_of_elements(); i++)
         {
            if (std::find(keys.begin(), keys.end(), bndry_vals[i]) != keys.end())
            {
               new_bnd_attribs.push_back(bndry_vals[i]);
               for (size_t j = 0; j < num_points_in_element; j++)
               {
                  new_connectivity.push_back(element_points[current_element + j]);
               }
            }

            current_element += num_points_in_element;
         }

         conduit::Node &attrib_node = meshBlueprint.fetch_existing("fields/boundary_attribute/values");
         attrib_node.reset();
         attrib_node.set(new_bnd_attribs);

         conduit::Node &conn_node = meshBlueprint.fetch_existing("topologies/boundary/elements/connectivity");
         conn_node.reset();
         conn_node.set(new_connectivity);

         // Disable updating the mesh vertices each cycle.  This is unnecessary, as this test problem has fixed
         // vertex positions.
         options["mesh_motion"] = 0;
         CALI_MARK_END("Teton_Init_BP_Fields");

         mesh.Clear();
#else
         std::cerr << "Unable to open mfem mesh, test driver was not configured with CMake's '-DENABLE_MFEM=ON'."
                   << std::endl;
         exit(1);
#endif
      }
      else
      {
         std::cerr << "Couldn't find mfem mesh at " << inputPath << std::endl;
         exit(1);
      }
      CALI_MARK_END("Teton_Read_Mfem_Input");
   }
   // Assume this is an input directory with a set of conduit blueprint mesh files and problem parameter files.
   // Note: The parameter files currently have all the input duplicated for each rank.  Look into making a
   // single 'global' parameter file for global input.
   else
   {
      CALI_MARK_BEGIN("Teton_Read_Conduit_Input");
      try
      {
         std::string input_file_path_base = inputPath + "/parameters_input_" + std::to_string(myRank);
         std::string input_file_path_full = input_file_path_base + ".hdf5";

         // Check for parameters node file ending in .hdf5.
         if (access(input_file_path_full.c_str(), F_OK) != -1)
         {
            conduit::relay::io::load_merged(input_file_path_full, "hdf5", options);
         }
         else
         {
            // Check for parameters node file ending in .conduit_json.
            input_file_path_full = input_file_path_base + ".conduit_json";
            if (access(input_file_path_full.c_str(), F_OK) != -1)
            {
               conduit::relay::io::load_merged(input_file_path_full, "conduit_json", options);
            }
            else
            {
               std::cerr << "Couldn't find parameters node at " << inputPath << std::endl;
               exit(1);
            }
         }

         // Check for mesh node file ending in .hdf5.
         input_file_path_base = inputPath + "/mesh_input_" + std::to_string(myRank);
         input_file_path_full = input_file_path_base + ".hdf5";
         if (access(input_file_path_full.c_str(), F_OK) != -1)
         {
            conduit::relay::io::load_merged(input_file_path_full, "hdf5", meshBlueprint);
         }
         else
         {
            // Check for mesh node file ending in .conduit_json.
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
         // Create Teton's expected mesh format
         if (myRank == 0)
         {
            std::cout << "Teton driver: creating Teton mesh node from blueprint node "
                      << "\n";
         }
      }
      catch (const std::exception &e)
      {
         std::cout << "Teton driver: Error loading conduit node files at " << inputPath << ": " << e.what();
         throw;
      }
      CALI_MARK_END("Teton_Read_Conduit_Input");
   }

   options["memory_allocator/umpire_host_allocator_id"] = -1;
   options["memory_allocator/umpire_device_allocator_id"] = -1;

   options["size/useCUDASweep"] = false;

   //==========================================================
   // Set problem options passed in via command line
   //==========================================================

   if (useDeviceAwareMPI == true)
   {
      options["mpi/useDeviceAddresses"] = 0;
   }

   if (useNewGTASolver == true)
   {
      options["size/useNewGTASolver"] = true;
   }

   if (useCUDASweep == true)
   {
      options["size/useCUDASweep"] = true;
      options["size/useGPU"] = false;

      if (useGPU == true)
      {
         std::cerr << "Select either the experimental CUDA sweep (useCUDASweep) or useGPU, not both." << std::endl;
         exit(1);
      }
   }
   else if (useGPU == true)
   {
      if (myRank == 0)
      {
         std::cout << "Teton driver: -g arg detected, enabling gpu kernels and associated options." << std::endl;
      }

      options["size/useGPU"] = useGPU;

      // If using the GPU, enable several sub-options.
      if (useUmpire > 0)
      {
#if defined(TETON_ENABLE_UMPIRE)
         if (myRank == 0)
         {
            std::cout << "Teton driver: Enabling use of Umpire CPU and GPU memory pools..." << std::endl;
         }

         auto &rm = umpire::ResourceManager::getInstance();

         // Create umpire allocators.
         auto host_pinned_pool = rm.makeAllocator<umpire::strategy::QuickPool>("HOST_PINNED_QUICK_POOL",
                                                                               rm.getAllocator("PINNED"));
         auto thread_safe_host_pinned_pool = rm.makeAllocator<umpire::strategy::ThreadSafeAllocator>(
            "THREAD_SAFE_PINNED_QUICK_POOL",
            host_pinned_pool);

         options["memory_allocator/umpire_host_allocator_id"] = thread_safe_host_pinned_pool.getId();
         if (useUmpire > 1)
         {
            auto device_pool = rm.makeAllocator<umpire::strategy::QuickPool>("DEVICE_QUICK_POOL",
                                                                             rm.getAllocator("DEVICE"));
            auto thread_safe_device_pool = rm.makeAllocator<umpire::strategy::ThreadSafeAllocator>(
               "THREAD_SAFE_DEVICE_QUICK_POOL",
               device_pool);
            options["memory_allocator/umpire_device_allocator_id"] = thread_safe_device_pool.getId();
         }
         else
         {
            options["memory_allocator/umpire_device_allocator_id"] = -1;
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
      options["size/useCUDASolver"] = true;
   }

   if (numPhaseAngleSets != 0)
   {
      options["quadrature/nSetsMaster"] = numPhaseAngleSets;
   }

   options["concurrency/nOmpMaxThreads"] = numOmpMaxThreads;

   //Set verbosity level ( default is 1 )
   options["verbose"] = verbose;

   if (fixedDT > 0)
   {
      options["iteration/dtrad"] = fixedDT;
   }

   // Can remove this when MFEM is fixed to stop writing out empty adjacency sets.
   // -- black27
   if (meshBlueprint.has_path("adjsets/main_adjset/groups")
       && meshBlueprint.fetch_existing("adjsets/main_adjset/groups").number_of_children() == 0)
   {
      meshBlueprint.remove("adjsets");
   }

   // Verify the blueprint is valid.
   std::string protocol = "mesh";
   conduit::Node info;
   conduit::blueprint::verify(protocol, meshBlueprint, info);

   CALI_MARK_BEGIN("Teton_Initialize");
   myTetonObject.initialize(comm);
   CALI_MARK_END("Teton_Initialize");

   // If a dtrad wasn't provided in the input file, the Teton initialize()
   // call will populate it with a default value.
   double dtrad = options.fetch_existing("iteration/dtrad").value();
   double timerad = 0.0;
   meshBlueprint["state/cycle"] = 0;

   // Calculate size of PSI to provide the number of unknowns solved for benchmarking, throughput calculations, etc.
   // This is # corners * # angles * # energy group bins

   // Put code here that calculates the # unknowns being solved.
   // TODO:
   // Some of this code ( especially the code that calculates the # angles ) can go in a helper function later, as
   // opposed to bloating up the test driver code.
   // -- black27

   // Get total number of corners in problem.
   const conduit::Node &corner_topology = meshBlueprint.fetch_existing("topologies/main_corner");
   unsigned long local_num_corners = conduit::blueprint::mesh::utils::topology::length(corner_topology);
   unsigned long num_corners;
   unsigned long num_unknowns;
   unsigned long local_num_unknowns;
   // int MPI_Reduce(_In_ void *sendbuf, _Out_opt_ void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
   int error_code = MPI_Reduce(&local_num_corners, &num_corners, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, comm);
   if (error_code != MPI_SUCCESS)
   {
      //TODO - error out
   }

   if (myRank == 0)
   {
      unsigned int num_groups = options.fetch_existing("quadrature/num_groups").to_unsigned_int();
      unsigned int num_polar_angles = options.fetch_existing("quadrature/npolar").to_unsigned_int();
      unsigned int num_azimuthal_angles = options.fetch_existing("quadrature/nazimu").to_unsigned_int();
      unsigned int num_angles = 0;
      unsigned int ndims = conduit::blueprint::mesh::utils::topology::dims(corner_topology);
      int quadrature_type = options.fetch_existing("quadrature/qtype").value();
      if (quadrature_type == 1) //level-symmetric quadrature
      {
         int quadrature_order = options.fetch_existing("quadrature/qorder").value();
         // Assume RZ, as we don't support XY
         if (ndims == 2)
         {
            num_angles = quadrature_order * (quadrature_order + 6) / 2;
         }
         else if (ndims == 3)
         {
            num_angles = quadrature_order * (quadrature_order + 2);
         }
         else
         {
            //TODO - Error out.
         }
      }
      else if (quadrature_type == 2)
      {
         // Assume RZ, as we don't support XY
         if (ndims == 2)
         {
            // 2D has four quadrants, and RZ has additional starting/finishing angles, so add one to azimuthal angles.
            num_angles = num_polar_angles * (num_azimuthal_angles + 1) * 4;
         }
         else if (ndims == 3)
         {
            // 3D has eight quadrants.
            num_angles = num_polar_angles * num_azimuthal_angles * 8;
         }
         else
         {
            //TODO - Error out.
         }
      }

      num_unknowns = num_corners * num_angles * num_groups;
      local_num_unknowns = local_num_corners * num_angles * num_groups;

      std::cout << "=================================================================" << std::endl;
      std::cout << "=================================================================" << std::endl;
      std::cout << "Teton starting cycling\n";
      std::cout << "=================================================================" << std::endl;
      std::cout << "Solving for " << num_unknowns << " global unknowns." << std::endl;
      std::cout << "(" << num_corners << " spatial elements * " << num_angles << " directions (angles) * " << num_groups << " energy groups)" << std::endl;
      // TODO - could beef this up to be a global memory estimate and global memory used, if we are testing problems with unbalanced mesh partition sizes, but this is meant as a rough memory estimate.
      std::cout << "Estimated CPU memory needed (rank 0): " << local_num_unknowns * sizeof(double) / 1024.0 / 1024.0 << "MB" << std::endl;
      std::cout << "Current CPU memory use (rank 0): " << getCurrentRSS() / 1024.0 / 1024.0 << "MB" << std::endl;
      std::cout << "=================================================================" << std::endl;
      std::cout << std::endl;
   }
   CALI_MARK_BEGIN("Teton_Cycle_Loop");

   double start_time = MPI_Wtime();
   double inner_start_time, inner_end_time, inner_elapsed_time;
   for (int cycle = 1; cycle <= cycles; cycle++)
   {
      if (dumpViz)
      {
         if (myRank == 0)
         {
            std::cout << "Teton driver: Dumping mesh blueprint for viz purposes." << std::endl;
         }
         myTetonObject.dump(comm, outputPath);
      }
      if (myRank == 0)
      {
         std::cout << "----------" << std::endl;
         std::cout << "CYCLE " << cycle << std::endl;
         std::cout << "----------" << std::endl;
      }
      timerad = timerad + dtrad;

      inner_start_time = MPI_Wtime();
      dtrad = myTetonObject.step(cycle);
      inner_end_time = MPI_Wtime();
      inner_elapsed_time = inner_end_time - inner_start_time;
      if (myRank == 0)
      {
         std::cout << "Teton driver: CYCLE WALL TIME = " << inner_elapsed_time << " seconds." << std::endl;
         std::cout << "Teton driver: CYCLE UNKNOWNS/SECOND = " << num_unknowns / inner_elapsed_time << std::endl;
      }
      // Either setTimeStep(cycle, dtrad) can be called to update time step, or
      // these can be directly updated in the options.
      if (fixedDT > 0)
      {
         dtrad = fixedDT;
      }
      options["iteration/dtrad"] = dtrad;
      options["iteration/timerad"] = timerad;
   }
   double end_time = MPI_Wtime();
   CALI_MARK_END("Teton_Cycle_Loop");

   double elapsed_time = end_time - start_time;

   if (myRank == 0)
   {
      double avg_unknowns_per_second = num_unknowns * cycles / (end_time - start_time);

      std::cout << std::endl;
      std::cout << "=================================================================" << std::endl;
      std::cout << "=================================================================" << std::endl;
      std::cout << "Teton finished cycling\n";
      std::cout << "=================================================================" << std::endl;
      std::cout << "Average number of unknowns/second solved per cycle was " << avg_unknowns_per_second << std::endl;
      std::cout << "=================================================================" << std::endl;
      std::cout << std::endl;

      // Appends # ranks and # unknowns solved per second to a .csv file.  Useful for scaling runs.
      // Also appends some problem state.
      if (benchmarkProblem > -1)
      {
         conduit::Node &datastore = myTetonObject.getDatastore();


         double energy_radiation = datastore.fetch_existing("rtedits/EnergyRadiation").value();
         double max_electron_temp = datastore.fetch_existing("rtedits/TeMax").value();
         double max_radiation_temp = datastore.fetch_existing("rtedits/TrMax").value();

         double power_incident = datastore.fetch_existing("rtedits/PowerIncident").value();
         double power_escape = datastore.fetch_existing("rtedits/PowerEscape").value();

         double power_absorbed = datastore.fetch_existing("rtedits/PowerAbsorbed").value();
         double power_emitted = datastore.fetch_existing("rtedits/PowerEmitted").value();

         double estimatedMemory = local_num_unknowns * sizeof(double);

         std::ofstream outfile;
         outfile.precision(16);
         std::string filePath = outputPath + "/" + label + ".csv";
         std::string firstline;
         if (access(filePath.c_str(), F_OK) == -1)
         {
            firstline = "# mpi ranks, memory estimate (bytes), memory actual (bytes), # solver unknowns, solver throughput (# unknowns solved per second), energy in radiation field, maximum electron temperature, maximum radiation temperature, incident power, escaping power, power absorbed, power emitted";
         }
         outfile.open(filePath, std::ios_base::app);
         if (!firstline.empty())
         {
            outfile << firstline << "\n";
         }
         outfile << mySize << ", " << getCurrentRSS() << ", " << num_unknowns << ", " << avg_unknowns_per_second
                 << ", " << energy_radiation << ", " << max_electron_temp << ", " << max_radiation_temp << ", "
                 << power_incident << ", " << power_escape << ", " << power_absorbed << ", " << power_emitted << "\n";

         outfile.close();
      }
   }

   // Tolerances were determined by running a scaling study on a LLNL Toss Sapphire Rapids machines from 1 to 112 ranks.
   // Values at one rank and 112 ranks were:
   //           # unknowns   # unknowns/sec     radiation energy      incident power  escaped power
   // 1 rank    17414258688, 19658662.13023822, 1.01692587612477e-09, 0,              1.727415298342991e-06
   // 110 ranks 17414258688, 737860099.4103793, 1.02731706435098e-09, 0,              1.745002013950014e-06

#if defined(TETON_ENABLE_UMPIRE)
   if (useGPU == 1 && useUmpire > 0)
   {
      if (myRank == 0)
      {
         std::cout << "Teton driver: Deleting Umpire CPU and GPU memory pools..." << std::endl;
      }

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

#if defined(TETON_ENABLE_MFEM)
   if (pmesh != nullptr)
   {
      delete pmesh;
   }
   if (conduit_data_collec != nullptr)
   {
      delete conduit_data_collec;
   }
#endif

   CALI_MARK_END("Teton_Test_Driver");
   //   MPI_Finalize();
   return 0;
}
