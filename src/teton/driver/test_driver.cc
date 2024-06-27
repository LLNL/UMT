#include "mpi.h"
#include <cmath>
#include <cstring>
#include <exception>
#include <getopt.h>
#include <malloc.h>
#include <sstream>
#include <stdio.h> // for getCurrentRSS and getPeakRSS
#include <string>
#include <sys/stat.h> //for mkdir
#include <unistd.h>   //for getopt, access

// XLF signal handler function to emit a stack trace.
#if defined(__ibmxl__)
#include "signal.h"
extern "C" void xl__trce(int, siginfo_t *, void *);
#endif

#include <signal.h>
#if defined(__linux__)
#include <fenv.h>
#endif

#if defined(TETON_ENABLE_UMPIRE)
#include "umpire/Umpire.hpp"
#include "umpire/strategy/QuickPool.hpp"
#include "umpire/strategy/ThreadSafeAllocator.hpp"
#endif

#if defined(TETON_ENABLE_CUDA)
#include "cuda_runtime_api.h"
#endif

#if defined(TETON_ENABLE_HIP)
#include "hip/hip_runtime_api.h"
#endif

#include "TetonConduitInterface.hh"
#include "TetonInterface.hh"

#include "conduit/conduit.hpp"
#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_blueprint_mesh_utils.hpp"
#include "conduit/conduit_blueprint_mpi_mesh_utils.hpp"
#include "conduit/conduit_relay.hpp"
#include "conduit/conduit_relay_mpi_io_blueprint.hpp"

#if defined(TETON_ENABLE_CALIPER)
#include "adiak.hpp"
#include "caliper/cali-manager.h"
#include "caliper/cali-mpi.h"
#include "caliper/cali.h"
#else
#define CALI_MARK_BEGIN(label)
#define CALI_MARK_END(label)
#define CALI_CXX_MARK_SCOPE(label)
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

#if defined(TETON_ENABLE_CUDA)
#define GPU_MEMGETINFO(free, total) cudaMemGetInfo(free, total)
#define GPU_SUCCESS cudaSuccess
#elif defined(TETON_ENABLE_HIP)
#define GPU_MEMGETINFO(free, total) hipMemGetInfo(free, total)
#define GPU_SUCCESS hipSuccess
#else
#define GPU_MEMGETINFO(free, total) 1
#define GPU_SUCCESS 1
#endif

void print_gpu_mem(const char *label)
{
   size_t free = 0;
   size_t total = 0;
   double gbtotal, gbfree, gbused;

   if (GPU_MEMGETINFO(&free, &total) != GPU_SUCCESS)
   {
      printf("MemGetInfo failed for GPU 0");
   }

   gbtotal = ((double) total) / (1024.0 * 1024.0 * 1024.0);
   gbfree = ((double) free) / (1024.0 * 1024.0 * 1024.0);
   gbused = ((double) total - free) / (1024.0 * 1024.0 * 1024.0);
   fprintf(stdout, "%s: total %7.4f GB; free %7.3f GB; used %7.3f GB\n", label, gbtotal, gbfree, gbused);
   fflush(stdout);
}

// Returns the current resident set size (physical memory use) measured in kbytes, or zero if the value cannot be determined on this OS.
size_t getCurrentRSS()
{
   long rss = 0L;
   FILE *fp = NULL;
   if ((fp = fopen("/proc/self/statm", "r")) == NULL)
      return (size_t) 0L; /* Can't open? */
   if (fscanf(fp, "%*s%ld", &rss) != 1)
   {
      fclose(fp);
      return (size_t) 0L; /* Can't read? */
   }
   fclose(fp);
   return (size_t) rss * (size_t) sysconf(_SC_PAGESIZE);
}

//---------------------------------------------------------------------------
/**
 @brief This Conduit error handler is invoked when Conduit would otherwise
        throw an exception. The error handler blocks forever, on purpose.
        This makes it easier to see the call stack that lead to the Conduit
        exception so we can more easily track down the offending code in a
        debugger.
 */
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   std::cout << rank << ": s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
   // This is on purpose.
   while (1)
      ;
}

//---------------------------------------------------------------------------
/**
 @brief The TetonDriver class contains the driver state, sets up a problem, and
        runs Teton via the Conduit interface.
 */
class TetonDriver
{
  public:
   TetonDriver() = default;
   ~TetonDriver();

   void initialize();
   int processArguments(int argc, char *argv[]);
   int execute();
   void finalize();
   void print_umpire_usage();

  private:
   void printUsage(const std::string &argv0) const;
   void startCaliper(const std::string &label);
   void writeStartSummary(unsigned int ndims, unsigned long num_corners, unsigned long &num_unknowns) const;
   void initializeBlueprintFields(int nelem, int numPolar, int numAzimuthal, int numGroups);
   void readConduitInputs();
   void setOptions();
   void verifyMesh();
   void cycleLoop(double &dtrad, double &timerad);
   void buildBlueprintTiledMesh();
   void writeEndSummary(double end_time, double start_time, unsigned long num_unknowns);

   void release();

  private:
   int return_status{0};
   int myRank{0};
   int mySize{0};
   int cycles{0};
   int numPhaseAngleSets{0};
   int useUmpire{2};
   int numOmpMaxThreads{1}; // Max number of CPU threads to use.
   double fixedDT{0.0};
   bool dumpViz{false};
   bool partition{false};
   double energy_check_tolerance{1.0e-6};
   int input_sanitizer_level{1};

   unsigned int benchmarkProblem{0};
   int numPolarUser{-1};
   int numAzimuthalUser{-1};
   int numGroupsUser{0};

   // MPI
   MPI_Comm comm{MPI_COMM_WORLD};
   int verbose{1};

   bool useGPU{false};
   bool useCUDASweep{false};
   int gta_kernel{1};
   int sweep_kernel{-1};
   int sweep_numhyperdomains{-1};
   int gta_numhyperdomains{-1};
   std::string scattering_kernel{};

   ::Teton::Teton myTetonObject{};
   std::string inputPath{"."};
   std::string outputPath{"."};
   std::string label{"unnamed"};
   std::string colorFile{};
   std::string caliper_config;
   std::string meshOrdering{"kdtree"};
#if defined(TETON_ENABLE_CALIPER)
   cali::ConfigManager mgr{};
#endif
   int blueprintMesh{0};
   int dims[3]{10, 10, 10}; //!< Number of cells in blueprint mesh.
   unsigned int total_num_flux_iterations{0};
};

//---------------------------------------------------------------------------
TetonDriver::~TetonDriver()
{
   release();
}

//---------------------------------------------------------------------------
void TetonDriver::initialize()
{
#if defined(TETON_ENABLE_CALIPER)
   cali_mpi_init();
   mgr.set_default_parameter("aggregate_across_ranks", "true");
   mgr.set_default_parameter("calc.inclusive", "true");
   mgr.set_default_parameter("main_thread_only", "true");
   mgr.add("runtime-report");

#if defined(TETON_ENABLE_CUDA)
   mgr.add("nvtx");
#elif defined(TETON_ENABLE_HIP)
   //mgr.add("roxtx");
#endif

#endif

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
}

//---------------------------------------------------------------------------
int TetonDriver::processArguments(int argc, char *argv[])
{
   //==========================================================
   // Get command line arguments
   //==========================================================

   while (1)
   {
      static struct option long_options[] = {{"apply_label", required_argument, 0, 'l'},
                                             {"benchmark_problem", required_argument, 0, 'b'},
                                             {"blueprint", required_argument, 0, 'B'},
                                             {"dims", required_argument, 0, 'd'},
                                             {"caliper", required_argument, 0, 'p'},
                                             {"input_sanitizer_level", required_argument, 0, 'y'},
                                             {"handler", no_argument, 0, 'H'},
                                             {"help", no_argument, 0, 'h'},
                                             {"input_path", required_argument, 0, 'i'},
                                             {"num_cycles", required_argument, 0, 'c'},
                                             {"dt", required_argument, 0, 'D'},
                                             {"num_phase_space_sets", required_argument, 0, 's'},
                                             {"num_threads", required_argument, 0, 't'},
                                             {"output_path", required_argument, 0, 'o'},
                                             {"umpire_mode", required_argument, 0, 'u'},
                                             {"mesh_ordering", required_argument, 0, 'M'},
                                             {"use_device_aware_mpi", no_argument, 0, 'm'},
                                             {"use_cuda_sweep", no_argument, 0, 'e'},
                                             {"use_gpu_kernels", no_argument, 0, 'g'},
                                             {"gta_kernel", required_argument, 0, 'n'},
                                             {"scattering_kernel", required_argument, 0, 'k'},
                                             {"sweep_kernel", required_argument, 0, 'S'},
                                             {"sweep_numhyperdomains", required_argument, 0, 'Z'},
                                             {"gta_numhyperdomains", required_argument, 0, 'z'},
                                             {"verbose", required_argument, 0, 'v'},
                                             {"write_viz_file", no_argument, 0, 'V'},
                                             {"partition", no_argument, 0, 'x'},
                                             {"num_Polar", required_argument, 0, 'P'},
                                             {"num_Azimuthal", required_argument, 0, 'A'},
                                             {"num_Groups", required_argument, 0, 'G'},
                                             {0, 0, 0, 0}};

      /* getopt_long stores the option index here. */
      int option_index = 0;

      auto optString = "A:B:b:c:D:d:eG:gHhi:k:l:M:mn:o:P:p:s:S:t:u:Vv:xy:";

      int opt = getopt_long(argc, argv, optString, long_options, &option_index);

      /* Detect the end of the options. */
      if (opt == EOF)
      {
         break;
      }

      switch (opt)
      {
         case 'B':
            if (strcmp(optarg, "local") == 0)
               blueprintMesh = 1;
            else if (strcmp(optarg, "global") == 0)
               blueprintMesh = 2;
            else
               blueprintMesh = 0;
            break;
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
         case 'd':
         {
            int d[3] = {1, 1, 1};
            if (sscanf(optarg, "%d,%d,%d", &d[0], &d[1], &d[2]) == 3)
            {
               dims[0] = std::max(d[0], 1);
               dims[1] = std::max(d[1], 1);
               dims[2] = std::max(d[2], 0); // Allow 0 for 2D
            }
            else
            {
               throw std::runtime_error("Invalid dimensions " + std::string(optarg));
            }
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
         case 'H':
            // Install an alternate Conduit error handler for debugging.
            conduit::utils::set_error_handler(conduit_debug_err_handler);
            break;
         case 'h':
            if (myRank == 0)
            {
               printUsage(argv[0]);
            }
            throw std::runtime_error("");
         case 'i':
            inputPath = std::string(optarg);
            break;
         case 'k':
            scattering_kernel = std::string(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: using scattering kernel " << scattering_kernel << "." << std::endl;
            }
            break;
         case 'l':
            label = std::string(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: this run will be identified as '" << label
                         << "' in any caliper spot reports." << std::endl;
            }
            break;
         case 'M':
            meshOrdering = std::string(optarg);
            if (meshOrdering != "normal" && meshOrdering != "kdtree" && meshOrdering != "hilbert")
            {
               throw std::runtime_error("Unsupported mesh ordering " + meshOrdering);
            }
            break;
         case 'n':
            gta_kernel = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: using gta solver version " << gta_kernel << "." << std::endl;
            }
            break;
         case 'S':
            sweep_kernel = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: using sweep kernel version " << sweep_kernel
                         << ". (0=zone sweep, 1=corner sweep)" << std::endl;
            }
            break;
         case 'Z':
            sweep_numhyperdomains = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: setting sweep number of hyper-domains " << sweep_numhyperdomains
                         << ". (0=automatic setting, >=1 number of sweep hyper-domains)" << std::endl;
            }
            break;
         case 'z':
            gta_numhyperdomains = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: setting new GTA number of hyper-domains " << gta_numhyperdomains
                         << ". (0=automatic setting, >=1 number of new GTA hyper-domains)" << std::endl;
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
                  auto configs = mgr.available_config_specs();
                  for (auto str : configs)
                  {
                     std::cout << mgr.get_documentation_for_spec(str.c_str()) << std::endl;
                  }
                  std::cout << std::endl;
                  std::cout << std::endl << "----------------------------------" << std::endl;
               }
               throw std::runtime_error(""); // Early return
            }
            break;
#endif
         case 's':
            numPhaseAngleSets = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of phase-angle sets to create: " << numPhaseAngleSets << std::endl;
            }
            break;
         case 'A':
            numAzimuthalUser = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of azimuthal angles: " << numAzimuthalUser << std::endl;
            }
            break;
         case 'P':
            numPolarUser = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of polar angles: " << numPolarUser << std::endl;
            }
            break;
         case 'G':
            numGroupsUser = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: number of energy groups: " << numGroupsUser << std::endl;
            }
            break;
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
         case 'x':
            // We only partition if there are multiple ranks.
            partition = mySize > 1;
            if (myRank == 0)
            {
               std::cout << "Teton driver: partitioning enabled." << std::endl;
            }
            break;
         case 'y':
            input_sanitizer_level = atoi(optarg);
            if (myRank == 0)
            {
               std::cout << "Teton driver: setting input_sanitizer_level to " << input_sanitizer_level << std::endl;
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

#if defined(TETON_ENABLE_CALIPER)
   // Add output path for spot dump.  This should occur after the command line args are processed.
   std::string spot_line("spot");
   spot_line = spot_line + "(output=" + outputPath + "/" + label + ".cali)";
   mgr.add(spot_line.c_str());
#endif

   return 0;
}

//---------------------------------------------------------------------------
void TetonDriver::printUsage(const std::string &argv0) const
{
   std::cout << "Usage: " << argv0 << "[OPTIONS]" << std::endl;
   std::cout
      << " -b, --benchmark_problem <0,1,2>      Run predefined UMT benchmark problem.  0 = user specified # angles and # groups."
      << std::endl;
   std::cout << " -c, --num_cycles <cycles>      Number of cycles to execute." << std::endl;
   std::cout << " -D, --dt float      Set a fixed DT per cycle." << std::endl;
   std::cout
      << " -e, --use_cuda_sweep           Use experimental CUDA sweep.  Do not specify this option and -g at the same time."
      << std::endl;
   std::cout << " -g, --use_gpu_kernels          Run solvers on GPU and enable associated sub-options, where supported."
             << std::endl;
   std::cout << " -H, --handler                  Install an alternate Conduit error handler for debugging."
             << std::endl;
   std::cout << " -h, --help                     Print this help and exit." << std::endl;
   std::cout << " -i, --input_path <path>        Path to input files." << std::endl;
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   std::cout
      << " -k, --scattering_kernel <none, fp, bc>       Select compton scattering kernel. none=no scattering, fp=fokker planck, bc=boltzmann"
      << std::endl;
#endif
   std::cout
      << " -l, --apply_label              Label this run.  This label will be used to identify this run in any caliper reports."
      << std::endl;
   std::cout << " -m, --use_device_aware_mpi     Use device-aware MPI for GPU runs." << std::endl;
   std::cout << " -n, --gta_kernel <0,1,2>       Select GTA solver kernel version. 0=use default" << std::endl;
   std::cout
      << " -S, --sweep_kernel <0,1>       Select sweep kernel version. 0=zone, 1=corner.  Note: corner sweep only available as a GPU kernel, on 3D meshes."
      << std::endl;
   std::cout
      << " -Z, --sweep_numhyperdomains <0+>       Set number of sweep hyper-domains. 0=automatic, or >=1 number of sweep hyper-domains."
      << std::endl;
   std::cout
      << " -z, --gta_numhyperdomains <0+>       Set number of new GTA hyper-domains. 0=automatic, or >=1 number of new GTA hyper-domains."
      << std::endl;
   std::cout << " -o, --output_path <path>       Path to generate output files, including any caliper spot dumps."
             << std::endl;
#if defined(TETON_ENABLE_CALIPER)
   std::cout << " -p, --caliper <string>         Caliper configuration profile.  Set <string> to 'help'"
             << " to get supported keywords.  'None' disabled caliper.  Default is 'runtime-report,spot'." << std::endl;
#endif
   std::cout << " -s, --num_phase_space_sets <num_sets>  Number of phase-angle sets to construct." << std::endl;
   std::cout << " -t, --num_threads <threads>    Max number of threads for cpu OpenMP parallel regions." << std::endl;
   std::cout << " -u, --umpire_mode <0,1,2>      0 - Disable umpire.  1 - Use Umpire for CPU allocations."
             << "  2 - Use Umpire for CPU and GPU allocations." << std::endl;
   std::cout << " -V, --write_viz_file           Output blueprint mesh vizualization file each cycle" << std::endl;
   std::cout << " -v, --verbose [0,1,2]    0 - quite  1 - informational(default)  2 - really chatty and dump files"
             << std::endl;
   std::cout
      << " -y, --input_sanitizer_level  0 - don't check inputs\n 1 - print one message for each bad input category\n 2 - print one message for each bad value of each bad category"
      << std::endl;
   std::cout << " -A, --num_Azimuthal <int>      Number azimuthal angles in an octant" << std::endl;
   std::cout << " -P, --num_Polar <int>          Number polar angles in an octant" << std::endl;
   std::cout << " -G, --num_Groups <int>         Number energy groups" << std::endl;
   std::cout << " -B, --blueprint local|global   Generate Blueprint tiled mesh in memory using the specified scheme."
             << " The \"local\" scheme creates the same sized mesh on each MPI rank, allowing for weak scaling. The "
             << "\"global\" scheme creates the specified mesh size globally and decomposes that size over the available"
             << "MPI ranks, allowing for strong scaling." << std::endl;
   std::cout
      << " -d, --dims i,j,k               The size of the Blueprint mesh in tiles in i,j,k. k=0 builds a 2D mesh."
      << std::endl;
   std::cout << " -M, --mesh_ordering order      The name of the mesh ordering to use (normal or kdtree)." << std::endl;
}

//---------------------------------------------------------------------------
int TetonDriver::execute()
{
   {
      conduit::Node &options = myTetonObject.getOptions();
      conduit::Node &meshBlueprint = myTetonObject.getMeshBlueprint();

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
      // ! Iteration control defaults.
      //==========================================================
      if (benchmarkProblem > 0)
      {
// If running UMT, only the sweep kernel is active.  Tighten the tolerance and
// increase the number of allowed inner flux iterations to enable it to converge.
#if defined(TETON_ENABLE_MINIAPP_BUILD)
         if (myRank == 0)
         {
            std::cerr
               << "Detected UMT run, fixing temperature iterations to one and increasing max flux iterations to enable convergence."
               << std::endl;
         }
         options["iteration/incidentFluxMaxIt"] = 99;
         if (!options.has_path("iteration/relativeTolerance"))
         {
            options["iteration/relativeTolerance"] = energy_check_tolerance / 10.0;
         }
         energy_check_tolerance = 1e-9;
#endif
         if (fixedDT <= 0.0)
         {
            fixedDT = 1e-3;
         }
         if (cycles == 0)
         {
            cycles = 5;
         }

         if (benchmarkProblem == 1)
         {
            numPolarUser = 3;
            numAzimuthalUser = 3;
            numGroupsUser = 128;
         }
         else if (benchmarkProblem == 2)
         {
            numPolarUser = 2;
            numAzimuthalUser = 2;
            numGroupsUser = 16;
         }
         else
         {
            std::cerr << "Teton driver: Custom benchmark problem #" << benchmarkProblem << std::endl;
         }
      }

      // More initialization
      startCaliper(label);

      //==========================================================
      // Read in conduit nodes or generate a default blueprint mesh.
      //
      // TODO - All this hard-coding can be moved into an input file that lives
      // alongside the mfem mesh file.
      // TODO - All this code for converting a mfem mesh and input to a blueprint
      // mesh should be moved to another function in another source file, so the
      // driver doesn't have all this here. -- black27
      //==========================================================
      //==========================================================

      if (blueprintMesh > 0)
      {
         buildBlueprintTiledMesh();
      }
      // Assume this is an input directory with a set of conduit blueprint mesh files and problem parameter files.
      // Note: The parameter files currently have all the input duplicated for each rank.  Look into making a
      // single 'global' parameter file for global input.
      else
      {
         readConduitInputs();
      }

      //==========================================================
      // Set problem options passed in via command line
      //==========================================================
      setOptions();

      // Verify the mesh.
      verifyMesh();

      // Initialize Teton
      myTetonObject.initialize(comm);

      // Calculate size of PSI to provide the number of unknowns solved for benchmarking, throughput calculations, etc.
      // This is # corners * # angles * # energy group bins

      // Put code here that calculates the # unknowns being solved.
      // TODO:
      // Some of this code ( especially the code that calculates the # angles ) can go in a helper function later, as
      // opposed to bloating up the test driver code.
      // -- black27

      // Get total number of corners in problem. Use part mesh in case partitioning has occurred.
      unsigned long local_num_corners = 0;
      unsigned long num_corners = 0;
      unsigned int ndims = 1;
      const conduit::Node &part = myTetonObject.getMeshBlueprintPart();
      if (part.has_path("topologies/main_corner"))
      {
         const conduit::Node &corner_topology = part.fetch_existing("topologies/main_corner");
         local_num_corners = conduit::blueprint::mesh::utils::topology::length(corner_topology);
         ndims = conduit::blueprint::mesh::utils::topology::dims(corner_topology);
      }
      else
      {
         int nelem = 1;
         if (part.has_path("topologies/main"))
         {
            nelem = part.fetch_existing("topologies/main/elements/dims/i").value();
         }
         else // Take the first topology
         {
            conduit::NodeConstIterator topologies = part.fetch_existing("topologies").children();
            if (!topologies.has_next())
            {
               std::cout << "There must be at least one topology in your mesh!" << std::endl;
               exit(1);
            }
            nelem = topologies.next().fetch_existing("elements/dims/i").value();
         }
         local_num_corners = 2 * nelem;
      }

      int error_code = MPI_Reduce(&local_num_corners, &num_corners, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, comm);
      if (error_code != MPI_SUCCESS)
      {
         //TODO - error out
      }

      unsigned long num_unknowns = 0;
      writeStartSummary(ndims, num_corners, num_unknowns);

      // If a dtrad wasn't provided in the input file, the Teton initialize()
      // call will populate it with a default value.
      double dtrad = options.fetch_existing("iteration/dtrad").value();
      double timerad = 0.0;
      if (options.has_path("iteration/timerad"))
      {
         timerad = options.fetch_existing("iteration/timerad").value();
      }
      meshBlueprint["state/cycle"] = 0;

      double start_time = MPI_Wtime();
      cycleLoop(dtrad, timerad);
      double end_time = MPI_Wtime();

      myTetonObject.dumpTallyToJson();

      writeEndSummary(end_time, start_time, num_unknowns);
   }

   return return_status;
}

//---------------------------------------------------------------------------
void TetonDriver::startCaliper(const std::string &label2)
{
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
            std::cerr << "Teton driver: Caliper config error: " << mgr.error_msg() << std::endl;
            exit(1);
         }
      }
      mgr.start();
   }
   if (!label2.empty())
   {
      adiak::value("ProblemName", label2, adiak_general);
   }
#endif
}

//---------------------------------------------------------------------------
void TetonDriver::initializeBlueprintFields(int nelem, int numPolar, int numAzimuthal, int numGroups)
{
   conduit::Node &options = myTetonObject.getOptions();
   conduit::Node &meshBlueprint = myTetonObject.getMeshBlueprint();

   if (numGroups < 1)
   {
      throw std::runtime_error(
         "Teton driver: Must specify number of energy groups angles via '-G#' or by specifying a benchmark problem via '-b#'.");
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
   if (numPolar < 1)
   {
      throw std::runtime_error(
         "Teton driver: Must specify number of polar angles via '-P#' or by specifying a benchmark problem via '-b#'.");
   }
   if (numAzimuthal < 1)
   {
      throw std::runtime_error(
         "Teton driver: Must specify number of azimuthal angles via '-A#' or by specifying a benchmark problem via '-b#'.");
   }
   int qtype = 2;
   int qorder = 10;
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

   // TODO: Make  Grid functions for ConduitDataCollection to handle
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

   // We only support one material on an  mesh
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
   meshBlueprint["fields/electron_temperature/values"].set(electron_temperature.data(), electron_temperature.size());

   meshBlueprint["fields/radiation_temperature/association"] = "element";
   meshBlueprint["fields/radiation_temperature/topology"] = "main";
   meshBlueprint["fields/radiation_temperature/values"].set(radiation_temperature.data(), radiation_temperature.size());

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
}

//---------------------------------------------------------------------------
void TetonDriver::readConduitInputs()
{
   CALI_CXX_MARK_SCOPE("Teton_Read_Conduit_Input");

   conduit::Node &options = myTetonObject.getOptions();
   conduit::Node &meshBlueprint = myTetonObject.getMeshBlueprint();

   std::string input_file_path_base = inputPath + "/parameters_input_" + std::to_string(myRank);
   std::string input_file_path_full = input_file_path_base + ".hdf5";

   // Check for parameters node file ending in .hdf5.
   int has_data = 0;
   std::string conduitErr;
   try
   {
      if (access(input_file_path_full.c_str(), F_OK) != -1)
      {
         conduit::relay::io::load_merged(input_file_path_full, "hdf5", options);
         has_data = 1;
      }
      else
      {
         // Check for parameters node file ending in .conduit_json.
         input_file_path_full = input_file_path_base + ".conduit_json";
         if (access(input_file_path_full.c_str(), F_OK) != -1)
         {
            conduit::relay::io::load_merged(input_file_path_full, "conduit_json", options);
            has_data = 1;
         }
      }
   }
   catch (std::exception &e)
   {
      conduitErr = e.what();
   }
   // Make sure all rank had data.
   int ranks_with_data = 0;
   MPI_Allreduce(&has_data, &ranks_with_data, 1, MPI_INT, MPI_SUM, comm);
   if (ranks_with_data != mySize)
   {
      if (has_data)
         throw std::runtime_error("Other ranks couldn't find a parameters node.");
      else if (!conduitErr.empty())
         throw std::runtime_error("Couldn't find mesh parameters at " + input_file_path_full + ":" + conduitErr);
      else
         throw std::runtime_error("Couldn't find mesh parameters at " + input_file_path_full);
   }

   // Check for mesh node file ending in .hdf5.
   input_file_path_base = inputPath + "/mesh_input_" + std::to_string(myRank);
   input_file_path_full = input_file_path_base + ".hdf5";
   has_data = 0;
   try
   {
      if (access(input_file_path_full.c_str(), F_OK) != -1)
      {
         conduit::relay::io::load_merged(input_file_path_full, "hdf5", meshBlueprint);
         has_data = 1;
      }
      else
      {
         // Check for mesh node file ending in .conduit_json.
         input_file_path_full = input_file_path_base + ".conduit_json";
         if (access(input_file_path_full.c_str(), F_OK) != -1)
         {
            conduit::relay::io::load_merged(input_file_path_full, "conduit_json", meshBlueprint);
            has_data = 1;
         }
      }
   }
   catch (std::exception &e)
   {
      conduitErr = e.what();
   }
   // Make sure all rank had data.
   ranks_with_data = 0;
   MPI_Allreduce(&has_data, &ranks_with_data, 1, MPI_INT, MPI_SUM, comm);
   if (ranks_with_data != mySize)
   {
      if (has_data)
         throw std::runtime_error("Other ranks couldn't find a mesh node.");
      else if (!conduitErr.empty())
         throw std::runtime_error("Couldn't find mesh node at " + input_file_path_full + ":" + conduitErr);
      else
         throw std::runtime_error("Couldn't find mesh node at " + input_file_path_full);
   }

   // Create Teton's expected mesh format
   if (myRank == 0)
   {
      std::cout << "Teton driver: creating Teton mesh node from blueprint node "
                << "\n";
   }

   // If certain objects exist in the blueprint now, remove them since they will be regenerated.
   const std::vector<std::string> removals{"topologies/main_corner",
                                           "topologies/main_face",
                                           "fields/face_attribute",
                                           "adjsets/main_corner",
                                           "adjsets/main_face"};
   for (const auto &path : removals)
   {
      if (meshBlueprint.has_path(path))
         meshBlueprint.remove(path);
   }
}

//---------------------------------------------------------------------------
void TetonDriver::setOptions()
{
   conduit::Node &options = myTetonObject.getOptions();
   options["memory_allocator/umpire_host_allocator_id"] = -1;
   options["memory_allocator/umpire_device_allocator_id"] = -1;

   options["size/useCUDASweep"] = false;

   if (!scattering_kernel.empty())
   {
      int tetonComptonFlag;
      if (scattering_kernel == "none")
      {
         tetonComptonFlag = 61;
      }
      else if (scattering_kernel == "fp")
      {
         tetonComptonFlag = 62;
      }
      else if (scattering_kernel == "bc")
      {
         tetonComptonFlag = 63;
      }
      else
      {
         std::cerr << "Unsupported option for scattering kernel of '" << scattering_kernel << "'." << std::endl;
         exit(1);
      }
      options["compton/internalComptonFlag"] = tetonComptonFlag;
   }

   if (gta_kernel > 0)
   {
      if (gta_kernel == 1)
      {
         options["size/useNewGTASolver"] = false;
         if (myRank == 0)
            std::cout << "Teton driver: Using older GTA kernel, version 1." << std::endl;
      }
      else if (gta_kernel == 2)
      {
         options["size/useNewGTASolver"] = true;
         if (myRank == 0)
            std::cout << "Teton driver: Using newer GTA kernel, version 2." << std::endl;
      }
      else
      {
         if (myRank == 0)
            std::cerr << "Invalid value for gta solver kernel version of " << gta_kernel << "." << std::endl;
         exit(1);
      }
   }

   if (sweep_kernel > -1)
   {
      options["sweep/kernel/version"] = sweep_kernel;
   }

   if (sweep_numhyperdomains > -1)
   {
      options["sweep/sn/numhyperdomains"] = sweep_numhyperdomains;
   }

   if (gta_numhyperdomains > -1)
   {
      options["sweep/gta/numhyperdomains"] = gta_numhyperdomains;
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
            if (useUmpire == 1)
            {
               std::cout
                  << "Teton driver: Enabling use of Umpire for single memory pool backed by CPU native allocator."
                  << std::endl;
            }
            else if (useUmpire == 2)
            {
               std::cout
                  << "Teton driver: Enabling use of Umpire for separate memory pools for host and accelerator (cuda or hip)."
                  << std::endl;
            }
            else if (useUmpire == 3)
            {
               std::cout
                  << "Teton driver: Enabling use of Umpire for single memory pool backed by accelerator (cuda or hip) allocator."
                  << std::endl;
            }
         }

         auto &rm = umpire::ResourceManager::getInstance();
         if (useUmpire == 1 || useUmpire == 2)
         {
            // Create umpire allocators.
            auto host_pinned_pool = rm.makeAllocator<umpire::strategy::QuickPool>("HOST_PINNED_QUICK_POOL",
                                                                                  rm.getAllocator("PINNED"));
            auto thread_safe_host_pinned_pool = rm.makeAllocator<umpire::strategy::ThreadSafeAllocator>(
               "THREAD_SAFE_PINNED_QUICK_POOL",
               host_pinned_pool);

            options["memory_allocator/umpire_host_allocator_id"] = thread_safe_host_pinned_pool.getId();
         }
         if (useUmpire == 2 || 3)
         {
            auto device_pool = rm.makeAllocator<umpire::strategy::QuickPool>("DEVICE_QUICK_POOL",
                                                                             rm.getAllocator("DEVICE"));
            auto thread_safe_device_pool = rm.makeAllocator<umpire::strategy::ThreadSafeAllocator>(
               "THREAD_SAFE_DEVICE_QUICK_POOL",
               device_pool);

            if (useUmpire == 2)
            {
               options["memory_allocator/umpire_device_allocator_id"] = thread_safe_device_pool.getId();
            }
            else
            {
               options["memory_allocator/umpire_host_allocator_id"] = thread_safe_device_pool.getId();
            }
         }
         print_umpire_usage();

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

   options["concurrency/omp_cpu_max_threads"] = numOmpMaxThreads;

   //Set verbosity level ( default is 1 )
   options["verbose"] = verbose;

   // Set sanitizer options:
   options["iteration/sanitizer/level"] = input_sanitizer_level;
   options["iteration/sanitizer/kill_if_bad"] = 1; // Have TetonConduitInterface kill the code

   if (fixedDT > 0)
   {
      options["iteration/dtrad"] = fixedDT;
   }

   options["partitioning"] = partition;
}

//---------------------------------------------------------------------------
void TetonDriver::verifyMesh()
{
   conduit::Node &meshBlueprint = myTetonObject.getMeshBlueprint();

   // Verify the blueprint is valid.
   std::string protocol = "mesh";
   conduit::Node info;
   conduit::blueprint::verify(protocol, meshBlueprint, info);
}

//---------------------------------------------------------------------------
void TetonDriver::cycleLoop(double &dtrad, double &timerad)
{
   CALI_CXX_MARK_SCOPE("Teton_Cycle_Loop");

   conduit::Node &options = myTetonObject.getOptions();
   const conduit::Node &datastore = myTetonObject.getDatastore();
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
      timerad = timerad + dtrad;
      options["iteration/timerad"] = timerad;

      double inner_start_time = MPI_Wtime();
      dtrad = myTetonObject.step(cycle);
      double inner_end_time = MPI_Wtime();
      double inner_elapsed_time = inner_end_time - inner_start_time;
      if (myRank == 0)
      {
         std::cout << "Teton driver: CPU MEM USE (rank 0): " << getCurrentRSS() / 1024.0 / 1024.0 << "MB" << std::endl;
         print_umpire_usage();
      }
      // Either setTimeStep(cycle, dtrad) can be called to update time step, or
      // these can be directly updated in the options.
      if (fixedDT > 0)
      {
         dtrad = fixedDT;
      }
      options["iteration/dtrad"] = dtrad;

      total_num_flux_iterations += datastore.fetch_existing("rtedits/ninrt").as_int();
   }
}

//---------------------------------------------------------------------------
/**
 @brief Factor an integer and return a vector of its factors. If the number
        is not prime, we omit 1, num from the vector.

 @param num The number to factor.

 @return A vector containing the factors of the number. The product of these
         values will equal the original number.
 */
std::vector<int> factor(int num)
{
   std::vector<int> factors;
   int prod = 1;
   int n = num;
   for (int f = 2; f <= num && prod < num; f++)
   {
      while (n % f == 0)
      {
         factors.push_back(f);
         n /= f;
         prod *= f;
      }
   }
   if (factors.empty())
   {
      factors.push_back(1);
      if (num > 1)
         factors.push_back(num);
   }
   return factors;
}

//---------------------------------------------------------------------------
/**
 @brief Decompose size into a suitable I,J,K decomposition (domains) and
        compute domainid, the current rank's domain, within the domains.

 @param rank The current MPI rank.
 @param size The number of MPI ranks.
 @param ndims The number of dimensions we're allowed to decompose.
 @param[out] domainid The I,J,K value for the rank within the domains decomposition.
 @paran[out] domains The I,J,K decomposition of size.
 */
void decompose(int rank, int size, int ndims, int domainid[3], int domains[3])
{
   domains[2] = domains[1] = domains[0] = 1;
   auto factors = factor(size);
   for (size_t i = 0; i < factors.size(); i++)
   {
      int dim = (ndims - 1) - (i % ndims);
      domains[dim] *= factors[factors.size() - 1 - i];
   }

   // Convert rank to an I,J,K domainid.
   domainid[0] = rank % domains[0];
   domainid[1] = (rank % (domains[0] * domains[1])) / domains[0];
   domainid[2] = rank / (domains[0] * domains[1]);
}

//---------------------------------------------------------------------------
void TetonDriver::buildBlueprintTiledMesh()
{
   CALI_CXX_MARK_SCOPE("Teton_Build_Tiled_Mesh");
   conduit::Node &options = myTetonObject.getOptions();
   conduit::Node &meshBlueprint = myTetonObject.getMeshBlueprint();

   // Entire problem domain is 1x1x1
   const double extents[] = {0., 1., 0., 1., 0., 1.};

   // Common options.
   conduit::Node bopts;
   bopts["meshname"] = "main";
   bopts["datatype"] = "int32";
   // If the inputPath was set to a file path and the file exists, try reading it
   // as a new tile definition.
   if ((endsWith(inputPath, ".yaml") || endsWith(inputPath, ".json")) && access(inputPath.c_str(), F_OK) == 0)
   {
      conduit::relay::io::load(inputPath, bopts["tile"]);
   }

   if (mySize > 1 && blueprintMesh == 2) // global - strong scaling mesh
   {
      // Make options
      bopts["numDomains"] = mySize;
      bopts["extents"].set(extents, 6);
      bopts["selectedDomains"].set(std::vector<int>{myRank}); // Just build the local domain
      bopts["curveSplitting"] = 0;                            // Disabled for now.

      // Build mesh
      conduit::blueprint::mesh::examples::tiled(dims[0], dims[1], dims[2], meshBlueprint, bopts);
   }
   else // local - weak scaling mesh
   {
      // Figure out a suitable domain decomposition for the number of ranks and
      // where the current rank exists within it.
      int ndims = (dims[2] > 0) ? 3 : 2;
      int domain[] = {0, 0, 0};
      int domains[] = {1, 1, myRank};
      decompose(myRank, mySize, ndims, domain, domains);

      // Figure out this domain's extents.
      double sideX = (extents[1] - extents[0]) / static_cast<double>(domains[0]);
      double sideY = (extents[3] - extents[2]) / static_cast<double>(domains[1]);
      double sideZ = (extents[5] - extents[4]) / static_cast<double>(domains[2]);
      double domainExt[] = {extents[0] + domain[0] * sideX,
                            extents[0] + (domain[0] + 1) * sideX,
                            extents[2] + domain[1] * sideY,
                            extents[2] + (domain[1] + 1) * sideY,
                            extents[4] + domain[2] * sideZ,
                            extents[4] + (domain[2] + 1) * sideZ};

      // Make options.
      bopts["domain"].set(domain, 3);
      bopts["domains"].set(domains, 3);
      bopts["extents"].set(domainExt, 6);
      bopts["reorder"] = meshOrdering;

      // Build the mesh
      conduit::blueprint::mesh::examples::tiled(dims[0], dims[1], dims[2], meshBlueprint, bopts);
   }

   const conduit::Node &main_topo = meshBlueprint.fetch_existing("topologies/main");
   auto nelem = conduit::blueprint::mesh::topology::length(main_topo);
   // Try to use the values provided by the user, if any.
   const int numPolar = (numPolarUser > -1) ? numPolarUser : 2;
   const int numAzimuthal = (numAzimuthalUser > -1) ? numAzimuthalUser : 2;
   const int numGroups = (numGroupsUser > 0) ? numGroupsUser : 2;
   initializeBlueprintFields(nelem, numPolar, numAzimuthal, numGroups);

   // Set up some boundary condition map information so Teton can interpret
   // the boundary_attribute field created by the tiled() function.
   std::vector<int> keys{1, 2, 3, 4, 5, 6};         // left+1, right+1, ... used by tiler.
   std::vector<int> values{35, 35, 35, 35, 35, 35}; // All vacuum
   options["boundary_conditions/id_to_type_map/ids"] = keys;
   options["boundary_conditions/id_to_type_map/types"] = values;

   // Fill in some options.
   options["sources/profile1/Values"] = 0.3;
   options["sources/profile1/NumTimes"] = 1;
   options["sources/profile1/NumValues"] = 1;
   options["sources/profile1/Multiplier"] = 1.0;

   // Disable updating the mesh vertices each cycle.
   options["mesh_motion"] = 0;
}

//---------------------------------------------------------------------------
void TetonDriver::writeStartSummary(unsigned int ndims, unsigned long num_corners, unsigned long &num_unknowns) const
{
   if (myRank == 0)
   {
      const conduit::Node &options = myTetonObject.getOptions();
      unsigned int num_angles = 0;
      unsigned int num_groups = options.fetch_existing("quadrature/num_groups").to_unsigned_int();
      if (ndims > 1)
      {
         unsigned int num_polar_angles = options.fetch_existing("quadrature/npolar").to_unsigned_int();
         unsigned int num_azimuthal_angles = options.fetch_existing("quadrature/nazimu").to_unsigned_int();
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
      }
      else //ndim == 1
      {
         int quadrature_order = options.fetch_existing("quadrature/qorder").value();
         num_angles = quadrature_order + 2; // assuming spherical
      }

      num_unknowns = num_corners * num_angles * num_groups;

      std::cout << "=================================================================" << std::endl;
      std::cout << "=================================================================" << std::endl;
      std::cout << "Test driver starting time steps\n";
      std::cout << "=================================================================" << std::endl;
      std::cout << "Solving for " << num_unknowns << " global unknowns." << std::endl;
      std::cout << "(" << num_corners << " spatial elements * " << num_angles << " directions (angles) * " << num_groups
                << " energy groups)" << std::endl;
      // TODO - could beef this up to be a global memory estimate and global memory used, if we are testing problems with unbalanced mesh partition sizes, but this is meant as a rough memory estimate.
      std::cout << "CPU memory needed per rank (average) for radiation intensity (PSI): "
                << num_unknowns / mySize * sizeof(double) / 1024.0 / 1024.0 << "MB" << std::endl;
      std::cout << "Current CPU memory use (rank 0): " << getCurrentRSS() / 1024.0 / 1024.0 << "MB" << std::endl;
      if (options.has_path("iteration/relativeTolerance"))
      {
         double relative_tol = options.fetch_existing("iteration/relativeTolerance").value();
         std::cout << "Iteration control: relative tolerance set to " << relative_tol << "." << std::endl;
      }
      std::cout << "=================================================================" << std::endl;
      std::cout << std::endl;
   }
}

//---------------------------------------------------------------------------
void TetonDriver::writeEndSummary(double end_time, double start_time, unsigned long num_unknowns)
{
   if (myRank == 0)
   {
      const conduit::Node &datastore = myTetonObject.getDatastore();
      double avg_unknowns_per_second = num_unknowns * cycles / (end_time - start_time);
      double avg_unknowns_per_second_per_iteration = num_unknowns * total_num_flux_iterations / (end_time - start_time);

      std::cout << std::endl;
      std::cout << "=================================================================" << std::endl;
      std::cout << "=================================================================" << std::endl;
      std::cout << "Test driver finished time steps\n";
      std::cout << "=================================================================" << std::endl;
      std::cout << "Average throughput of single iteration of iterative solver was "
                << avg_unknowns_per_second_per_iteration << " unknowns calculated per second." << std::endl;
      std::cout << "Throughput of iterative solver was " << avg_unknowns_per_second
                << " unknowns calculated per second." << std::endl;
      std::cout << "(average throughput of single iteration * # iterations for solver to produce answer" << std::endl;
      std::cout << std::endl;
      std::cout << "Total number of flux solver iterations for run: " << total_num_flux_iterations << std::endl;
      std::cout << "Total wall time for run: " << end_time - start_time << " seconds." << std::endl;
      std::cout << "=================================================================" << std::endl;
      std::cout << std::endl;

      // Appends # ranks and # unknowns solved per second to a .csv file.  Useful for scaling runs.
      // Also appends some problem state.
      if (benchmarkProblem > 0)
      {
         double energy_radiation = datastore.fetch_existing("rtedits/EnergyRadiation").value();
         double max_electron_temp = datastore.fetch_existing("rtedits/TeMax").value();
         double max_radiation_temp = datastore.fetch_existing("rtedits/TrMax").value();

         double power_incident = datastore.fetch_existing("rtedits/PowerIncident").value();
         double power_escape = datastore.fetch_existing("rtedits/PowerEscape").value();

         double power_absorbed = datastore.fetch_existing("rtedits/PowerAbsorbed").value();
         double power_emitted = datastore.fetch_existing("rtedits/PowerEmitted").value();

         double energy_check = datastore.fetch_existing("rtedits/EnergyCheck").value();

         double memForPSI = num_unknowns / mySize * sizeof(double);

         std::ofstream outfile;
         outfile.precision(16);
         std::string filePath = outputPath + "/" + label + ".csv";

         // Check if file exists.  If not, need to add header line for .csv, otherwise just append new rows.
         bool isNewFile = (access(filePath.c_str(), F_OK) == -1);

         const std::string header = "# mpi ranks, "
                                    "Mem for PSI (kb), "
                                    "process rss mem (kb), "
                                    "# solver unknowns (extents of PSI), "
                                    "total # flux iterations, "
                                    "# time steps, "
                                    "walltime(seconds),"
                                    "energy check, "
                                    "energy in radiation field, "
                                    "maximum electron temperature, "
                                    "maximum radiation temperature, "
                                    "incident power, "
                                    "escaping power, "
                                    "power absorbed, "
                                    "power emitted"
                                    "\n";

         outfile.open(filePath, std::ios_base::app);
         if (isNewFile)
         {
            outfile << header;
         }
         outfile << mySize << ", " << memForPSI / 1024.0 << ", " << getCurrentRSS() / 1024.0 << ", " << num_unknowns
                 << ", " << total_num_flux_iterations << ", " << cycles << ", " << end_time - start_time << ", "
                 << energy_check << ", " << energy_radiation << ", " << max_electron_temp << ", " << max_radiation_temp
                 << ", " << power_incident << ", " << power_escape << ", " << power_absorbed << ", " << power_emitted
                 << "\n";

         outfile.close();

         double rel_energy_check_result = std::abs(energy_check / (energy_radiation + 1.0e-50));
         if (rel_energy_check_result <= energy_check_tolerance)
         {
            std::cout << "RESULT CHECK PASSED: Energy check (this is relative to total energy) "
                      << rel_energy_check_result << " within tolerance of +/- " << energy_check_tolerance << "; check '"
                      << filePath << "' for tally details\n"
                      << std::endl;
         }
         else
         {
            std::cerr << "RESULT CHECK FAILED: Energy check (this is relative to total energy) "
                      << rel_energy_check_result << " exceeded tolerance of +/- " << energy_check_tolerance
                      << "; check '" << filePath << "' for tally details\n"
                      << std::endl;

            std::cerr << "Energy radiation: " << energy_radiation << std::endl;
            std::cerr << "Power incident: " << power_incident << std::endl;
            std::cerr << "Power escaped: " << power_escape << std::endl;
            std::cerr << "Power absorbed: " << power_absorbed << std::endl;
            std::cerr << "Power emitted: " << power_emitted << std::endl << std::endl;

            return_status = 1;
         }
      }
   }
}

//---------------------------------------------------------------------------
void TetonDriver::print_umpire_usage()
{
#if defined(TETON_ENABLE_UMPIRE)
   if (useUmpire > 0)
   {
      conduit::Node &options = myTetonObject.getOptions();

      auto &rm = umpire::ResourceManager::getInstance();

      int host_allocator_id = options.fetch_existing("memory_allocator/umpire_host_allocator_id").value();
      int device_allocator_id = options.fetch_existing("memory_allocator/umpire_device_allocator_id").value();

      if (host_allocator_id != -1)
      {
         auto pool = rm.getAllocator(host_allocator_id);
         if (myRank == 0)
         {
            print_bytes_as_gb("Teton driver: Umpire pool size: ", pool.getActualSize());
         }
      }

      if (device_allocator_id != -1)
      {
         auto pool = rm.getAllocator(device_allocator_id);
         if (myRank == 0)
         {
            print_bytes_as_gb("Teton driver: Umpire device pool size: ", pool.getActualSize());
         }
      }
   }
#endif
}

//---------------------------------------------------------------------------
void TetonDriver::finalize()
{
#if defined(TETON_ENABLE_UMPIRE)
   if (useUmpire > 0)
   {
      auto &rm = umpire::ResourceManager::getInstance();
      conduit::Node &options = myTetonObject.getOptions();

      int host_allocator_id = options.fetch_existing("memory_allocator/umpire_host_allocator_id").value();
      int device_allocator_id = options.fetch_existing("memory_allocator/umpire_device_allocator_id").value();

      if (host_allocator_id != -1)
      {
         rm.getAllocator(host_allocator_id).release();
      }

      if (device_allocator_id != -1)
      {
         rm.getAllocator(host_allocator_id).release();
      }
   }
#endif

   release();

#if defined(TETON_ENABLE_CALIPER)
   adiak::fini();
   mgr.flush();
#endif
}

//---------------------------------------------------------------------------
void TetonDriver::release()
{
}

//==========================================================
//  run syntax:
//
// ./test_driver -c <# cycles> -i <input path> -o <output path>
// Note: to use the new blueprint format, ./test_driver -b ...
//==========================================================

int main(int argc, char *argv[])
{
   int retval = 0;

   //==========================================================
   // Initialize MPI and MPI thread support.
   //==========================================================
   int claimed = 0, provided = 0;
#if defined(TETON_ENABLE_OPENMP)
   int request = MPI_THREAD_MULTIPLE;
#else
   int request = MPI_THREAD_SINGLE;
#endif

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

   try
   {
      CALI_CXX_MARK_SCOPE("Teton_Test_Driver");
      TetonDriver driver;
      driver.initialize();
      driver.processArguments(argc, argv);
      retval = driver.execute();
      driver.finalize();
   }
   catch (const std::exception &e)
   {
      std::string w(e.what());
      if (!w.empty())
         std::cerr << "Teton driver: " << w << std::endl;
      retval = 1;
   }

   // NOTE: The teton object must be destroyed first, before MPI is finalized.  Teton uses MPI during its destructor.
   MPI_Finalize();

   return retval;
}
