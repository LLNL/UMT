#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_blueprint_mesh.hpp"
#include "conduit/conduit_blueprint_mesh_utils.hpp"
#include "conduit/conduit_blueprint_mpi_mesh.hpp"
#include "conduit/conduit_config.h"
#include "conduit/conduit_relay.hpp"
#include "conduit/conduit_relay_config.h"
#include "conduit/conduit_relay_mpi.hpp"
#include "conduit/conduit_relay_mpi_io_blueprint.hpp"
#if defined(CONDUIT_USE_PARMETIS)
#include "conduit/conduit_blueprint_mesh_topology_metadata.hpp"
// We can only enable partitioning right now if Conduit includes Parmetis.
#include "conduit/conduit_blueprint_mpi_mesh_parmetis.hpp"
#pragma message "Teton built with partitioning support."
#define TETON_PARTITIONING
#endif

#include "TetonBlueprint.hh"
#include "TetonConduitInterface.hh"
#include "TetonInterface.hh"
#include "TetonNDAccessor.hh"
#include "TetonSurfaceTallies.hh"
#include "TetonTesting.hh"
#include "TetonUtilities.hh"
#include "dbc_macros.h"

#if defined(TETON_USE_CUDA)
#include "cuda.h"
#include "cuda_runtime_api.h"
#endif

#if defined(TETON_ENABLE_CALIPER)
#pragma message "Teton built with Caliper support."
#include "caliper/cali.h"
#else
#define CALI_MARK_BEGIN(label)
#define CALI_MARK_END(label)
#define CALI_CXX_MARK_SCOPE(name)
#define CALI_CXX_MARK_FUNCTION
#endif

// Uncomment to enable partition debugging console output.
// #define PARTITION_DEBUG

extern "C"
{
extern conduit::Node *teton_get_datastore_cptr();
extern conduit::Node *teton_conduitcheckpoint_get_cptr();

extern void teton_conduitcheckpoint_prep_for_load();
extern void teton_conduitcheckpoint_data_loaded();
extern void teton_conduitcheckpoint_external_data_loaded();
extern void teton_conduitcheckpoint_prep_for_save();
extern void teton_conduitcheckpoint_teardown();
}

namespace Teton
{

std::string field_path(const std::string &fieldName)
{
   return "fields/" + fieldName;
}

std::string field_values(const std::string &fieldName)
{
   return field_path(fieldName) + "/values";
}

//---------------------------------------------------------------------------
// Teton
const std::string Teton::PREFIX("__teton__");
const std::string Teton::MCARRAY_PREFIX(Teton::PREFIX + "mcarray_");
const std::string Teton::PARTITION_FIELD(Teton::PREFIX + "parmetis_result");
const std::string Teton::PARTITION_FIELD_BOUNDARY(Teton::PREFIX + "parmetis_result_boundary");

const std::string Teton::FIELD_ELECTRON_ENERGY_DEPOSITED("electron_energy_deposited");
const std::string Teton::FIELD_RADIATION_ENERGY_DENSITY("radiation_energy_density");
// We gather radiation_temperature results into this field so values may be queried by
// getRadiationTemperature(). Prepend the prefix so we don't disturb the radiation_temperature
// field if the host code happens to provide one. We take this approach for some of the
// other fields below as well.
const std::string Teton::FIELD_RADIATION_TEMPERATURE(Teton::PREFIX + "radiation_temperature");

const std::string Teton::FIELD_RADIATION_FORCE_X("radiation_force_x");
const std::string Teton::FIELD_RADIATION_FORCE_Y("radiation_force_y");
const std::string Teton::FIELD_RADIATION_FORCE_Z("radiation_force_z");
const std::string Teton::FIELD_RADIATION_FORCE_R("radiation_force_r");
const std::string Teton::FIELD_CORNER_VOLUME_SUMS(Teton::PREFIX + "cornerVolumeSums");

const std::string Teton::FIELD_RADIATION_FLUX_X(Teton::PREFIX + "radiation_flux_x");
const std::string Teton::FIELD_RADIATION_FLUX_Y(Teton::PREFIX + "radiation_flux_y");
const std::string Teton::FIELD_RADIATION_FLUX_Z(Teton::PREFIX + "radiation_flux_z");
const std::string Teton::FIELD_RADIATION_FLUX_R(Teton::PREFIX + "radiation_flux_r");

// This field is handled similiar to FIELD_RADIATION_TEMPERATURE.
const std::string Teton::FIELD_MATERIAL_TEMPERATURE(Teton::PREFIX + "material_temperature");

const std::string Teton::TOPO_MAIN("main");
const std::string Teton::TOPO_BOUNDARY("boundary");

Teton::Teton()
   : mDTrad(0.),
     areSourceProfilesSet(false),
     mIsInitialized(false),
     mGTAorder(2),
     mInternalComptonFlag((int) tetonComptonFlag::none),
     mCommunicator(MPI_COMM_WORLD),
     mRank(0),
     mCornerToVertex(),
     mZoneToNCorners(),
     mZoneToCorners(),
     mCornerToZone(),
     mMapBackFields(),
     mMCArrays(),
     mRadiationForceDensityFields(),
     mRadiationFluxFields()
{
}

Teton::~Teton()
{
   if (mIsInitialized)
   {
      bool enableNLTE = false;
      if (getDatastore().has_path("options"))
      {
         conduit::Node &options = getOptions();
         if (options.has_path("size/enableNLTE"))
         {
            enableNLTE = options.fetch_existing("size/enableNLTE").as_int();
         }
      }
      teton_destructmeshdata(&enableNLTE);

      teton_destructmemoryallocator();
   }
}

int Teton::getVerbose() const
{
   const conduit::Node &options = getOptions();
   int verbose = 0;
   if (options.has_path("verbose"))
   {
      verbose = options.fetch_existing("verbose").value();
      // Initialize from the environment if TETON_VERBOSE is set.
      if (getenv("TETON_VERBOSE") != nullptr)
         verbose = atoi(getenv("TETON_VERBOSE"));
   }
   else
   {
      // Initialize from the environment if TETON_VERBOSE is set.
      if (getenv("TETON_VERBOSE") != nullptr)
      {
         verbose = atoi(getenv("TETON_VERBOSE"));
      }
   }
   return verbose;
}

void Teton::initialize(MPI_Comm communicator, bool fromRestart)
{
   CALI_MARK_BEGIN("Teton_Initialize");

   mCommunicator = communicator;
   MPI_Fint fcomm = MPI_Comm_c2f(communicator);

   MPI_Comm_rank(communicator, &mRank);
   mSourceManager.SetRank(mRank);

   conduit::Node &datastore = getDatastore();
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();
   conduit::Node &part = getMeshBlueprintPart();

   int verbose = 0;
   options["verbose"] = verbose = getVerbose();
   if (verbose && mRank == 0)
      std::cout << "Teton: setting verbosity to " << verbose << std::endl;

   if (verbose >= 2)
   {
      // Save parameters.
      if (mRank == 0)
      {
         std::cerr << "Teton: Dump copy of input..." << std::endl;
      }
      conduit::relay::io::save(options, "parameters_input_" + std::to_string(mRank) + ".conduit_json", "conduit_json");
      conduit::relay::io::save(options, "parameters_input_" + std::to_string(mRank) + ".json", "json");

      // Save mesh.
      if (verbose >= 3)
         conduit::blueprint::mesh::paint_adjset("main_adjset", "main_adjset", blueprint);
      conduit::relay::io::save(blueprint, "mesh_input_" + std::to_string(mRank) + ".conduit_json", "conduit_json");
      conduit::relay::io::save(blueprint, "mesh_input_" + std::to_string(mRank) + ".json", "json");
#if defined(PARTITION_DEBUG) && defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
#pragma message "Saving blueprint_initialize."
      if (mRank == 0)
      {
         std::cerr << "Teton: Save mesh..." << std::endl;
      }
      MPI_Barrier(communicator);
      // Save to a plottable file.
      conduit::relay::mpi::io::blueprint::save_mesh(blueprint, "blueprint_initialize", "hdf5", communicator);
#endif
   }

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   // Make sure that the blueprint mesh has temperature fields on it because
   // we need to map these values back.
   createRadiationTemperature();
   createMaterialTemperature();
#endif
   initializeRadiationFluxFieldNames();
   initializeRadiationForceDensityFieldNames();

   // Partition the mesh, if necessary. This migrates all fields to the partition mesh.
   partition(fromRestart);
   // The "part" node now contains partitioned data.

   // Create secondary (corner) mesh topology and connectivity arrays, using the part mesh.
   CALI_MARK_BEGIN("Teton_Construct_Corner_Mesh");
   TetonBlueprint blueprintHelper(part, options);
   blueprintHelper.OutputTetonMesh(mRank, mCommunicator);
   CALI_MARK_END("Teton_Construct_Corner_Mesh");

   if (verbose >= 2)
   {
      if (verbose >= 3)
         conduit::blueprint::mesh::paint_adjset("main_corner", "corner_adjset", part);
      if (mRank == 0)
      {
         std::cerr << "Teton: Dump blueprint with generated topologies..." << std::endl;
      }
      dump(communicator, ".");
   }

   // Set the Dt controls early, as tfloor is set here and needed by constructSize in the line below.
   constructDtControls();
   constructSize();
   constructMemoryAllocator();
   constructQuadrature();
   constructBoundaries();
   setSourceProfiles();
   teton_constructgeometry();
   setMeshConnectivity();

   MPI_Barrier(communicator);

   int ndim = options.fetch_existing("size/ndim").value();
   if (ndim > 1)
   {
      teton_setoppositeface(); //Prerequisite for calling setMeshSizeAndPositions() in 2D/3D
      setCommunication();      //Prerequisite for calling setMeshSizeAndPositions()
   }
   setMeshSizeAndPositions();

   // This is an awful hack to make sure Z%VolumeOld is initialized
   teton_getvolume();
   teton_getvolume();
   teton_setcommunicationgroup(&fcomm);
   teton_checksharedboundary();

   bool enableNLTE = options.fetch_existing("size/enableNLTE").as_int();
   teton_constructmaterial(&enableNLTE);

   teton_constructphasespacesets(&fromRestart);

   // This initializes the various zone fields, including the zonal electron
   // temperatures.
   // // NOTE: Teton's corner electron temperature is NOT set yet!
   // Do not call anything that needs this, until teton_initteton is called
   // below.
   setMaterials();

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   constructComptonControl();
#endif

   teton_constructradintensity();

   double EnergyRadiation = 0.0;

   // The corner temperature field that is passed in here will get set to the zonal
   // temperature field in Teton.  We're going to pass in Teton's own corner field
   // so we can get it initialized properly to the zonal field!
   double *tec = datastore.fetch_existing("material/Tec").value();
   teton_initteton(&EnergyRadiation, tec);
   datastore["rtedits/EnergyRadiation"] = EnergyRadiation;

   constructEdits();
   constructIterationControls();

   // Initialize default sweep kernel selection, if not provided.
   // Valid values are:
   // 0 - a sweep implementation that loops over the zones in each hyperplane, and the corners within those zones.
   // 1 - a sweep implementation that loops directly over the corners of each hyperplane.
   if (!options.has_path("sweep/kernel/version"))
   {
      // Default to the original 'loop over zones' version.
      options["sweep/kernel/version"] = 0;
   }

   // Initialize default numbers of hyper domains, if not provided.
   // Valid values are:
   // 0 - Indicates user did not set number of hyper-domains. Teton will automatically determine a number. (default)
   // >= 1 - Overrides Teton's automatic setting with this number.
   // Number of sweep hyper-domains:
   if (!options.has_path("sweep/sn/numhyperdomains"))
   {
      // Default to 0 to allow Teton's automatic algorithm to determine a number.
      options["sweep/sn/numhyperdomains"] = 0;
   }
   else if (options["sweep/sn/numhyperdomains"].as_int() < 0)
   {
      std::cerr
         << "Teton: Invalid number of sweep hyper-domains, please set options/sweep/sn/numhyperdomains to 0 for automatic setting, or a value >= 1 ..."
         << std::endl;
   }

   // Number of new GTA hyper-domains:
   if (!options.has_path("sweep/gta/numhyperdomains"))
   {
      // Default to 0 to allow Teton's automatic algorithm to determine a number.
      options["sweep/gta/numhyperdomains"] = 0;
   }
   else if (options["sweep/gta/numhyperdomains"].as_int() < 0)
   {
      std::cerr
         << "Teton: Invalid number of new GTA hyper-domains, please set options/sweep/gta/numhyperdomains to 0 for automatic setting, or a value >= 1 ..."
         << std::endl;
   }

   if (blueprint.has_path("fields/absorption_opacity/values"))
   {
      updateOpacity();
   }

   // Store mesh data needed for corner forces and update
   // of the mesh coordinates. In particular, the array corner_to_vertex
   // is stored, which associates the Teton corner id with the mesh vertex id
   storeMeshData();

   mIsInitialized = true;

   CALI_MARK_END("Teton_Initialize");
}

void Teton::storeMeshData()
{
   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();

   // To compute the radiation forces, Teton needs to hang on to
   // this connectivity array
   if (part.has_path("arrays/corner_to_vertex"))
   {
      int ncornr = options.fetch_existing("size/ncornr").to_int();
      mCornerToVertex.resize(ncornr);
      int *corner_to_vert_ptr = part.fetch_existing("arrays/corner_to_vertex").value();
      for (int c = 0; c < ncornr; ++c)
      {
         // Store the vertex ID corresponding to this corner ID.
         // !!!! NOTE: THIS WILL NOT WORK FOR AMR MESHES !!!!
         mCornerToVertex[c] = corner_to_vert_ptr[c];
      }
   }
   // To compute the radiation forces, Teton also needs to hang on to
   // this connectivity array
   if (part.has_path("relations/corner_to_zone"))
   {
      int ncornr = options.fetch_existing("size/ncornr").to_int();
      int *corner_to_zone_ptr = part.fetch_existing("relations/corner_to_zone").value();
      mCornerToZone.resize(ncornr);
      for (int c = 0; c < ncornr; ++c)
      {
         mCornerToZone[c] = corner_to_zone_ptr[c];
      }
   }

   if (part.has_path("arrays/zone_to_ncorners"))
   {
      const conduit::Node &zones_to_ncorners = part.fetch_existing("arrays/zone_to_ncorners");
      const int *zone_to_ncorner_ptr = zones_to_ncorners.value();
      auto n = static_cast<size_t>(zones_to_ncorners.dtype().number_of_elements());
      mZoneToNCorners.resize(n);
      memcpy(&mZoneToNCorners[0], zone_to_ncorner_ptr, n * sizeof(int));
   }
   if (part.has_path("arrays/zone_to_corners"))
   {
      const conduit::Node &zones_to_corners = part.fetch_existing("arrays/zone_to_corners");
      const int *zone_to_corner_ptr = zones_to_corners.value();
      auto n = static_cast<size_t>(zones_to_corners.dtype().number_of_elements());
      mZoneToCorners.resize(n);
      memcpy(&mZoneToCorners[0], zone_to_corner_ptr, n * sizeof(int));
   }
}

int Teton::checkInputSanity(const conduit::Node &sanitizer_node) const
{
   // level = 0 --> Don't run it
   // level = 1 --> print one complaint per problematic category
   // level = 2 --> print one complaint per problematic zone/corner
   const int sanitizer_level = sanitizer_node.fetch_existing("level").to_int();
   if (sanitizer_level < 1)
      return 0;

   // We'll kill the code in the C++ rather than the Fortran
   bool kill_if_bad = false;

   //Check all categories except scattering opacity:
   int num_cats = 6;
   const int default_cat_list[] = {2, 3, 4, 5, 6, 7};
   const int *cat_list_ptr = default_cat_list;
   if (sanitizer_node.has_path("cat_list"))
   {
      num_cats = sanitizer_node.fetch_existing("cat_list").dtype().number_of_elements();
      cat_list_ptr = sanitizer_node.fetch_existing("cat_list").as_int_ptr();
   }

   int num_bad_cats = 0;
   teton_checkinputsanity(&kill_if_bad, &sanitizer_level, &num_cats, cat_list_ptr, &num_bad_cats);

   if (sanitizer_node.has_path("kill_if_bad") && sanitizer_node.fetch_existing("kill_if_bad").to_int())
   {
      TETON_VERIFY_C(mRank, num_bad_cats == 0, "Bad inputs detected!");
   }

   return num_bad_cats;
}

void Teton::constructBoundaries()
{
   conduit::Node &options = getOptions();

   int nrefl = options.fetch_existing("boundary_conditions/num_reflecting").value();
   int nvac = options.fetch_existing("boundary_conditions/num_vacuum").value();
   int nsrc = options.fetch_existing("boundary_conditions/num_source").value();
   int num_comm = options.fetch_existing("boundary_conditions/num_comm").value();

   teton_constructboundary(&nrefl, &nvac, &nsrc, &num_comm);

   int ndim = options.fetch_existing("size/ndim").value();
   if (ndim > 1)
   {
      int numBCTotal = options.fetch_existing("boundary_conditions/num_total").value();
      int *BCTypeInt = options.fetch_existing("boundary_conditions/type").value();
      int *BCCornerFaces = options.fetch_existing("boundary_conditions/corner_face_ids").value();
      int *BCNeighborID = options.fetch_existing("boundary_conditions/neighbor_ids").value();

      TETON_VERIFY_C(mRank, (numBCTotal > 0), "No boundary conditions defined.");

      teton_addboundary(&numBCTotal, &BCTypeInt[0], &BCCornerFaces[0], &BCNeighborID[0]);
   }
   else
   {
      int *BCTypeInt = options.fetch_existing("boundary_conditions/type").value();
      int *BCNeighborID = options.fetch_existing("boundary_conditions/neighbor_ids").value();
      int *BCCornerFaces = options.fetch_existing("boundary_conditions/bc_ncorner_faces").value();
      int numBCTotal = 2;

      teton_addboundary(&numBCTotal, &BCTypeInt[0], &BCCornerFaces[0], &BCNeighborID[0]);
   }
}

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
void Teton::constructComptonControl()
{
   conduit::Node &options = getOptions();

   if (options.has_path("compton"))
   {
      mInternalComptonFlag = options.fetch_existing("compton/internalComptonFlag").value();
   }
   else
   {
      mInternalComptonFlag = (int) tetonComptonFlag::none;
   }

   teton_constructcomptoncontrol(&mInternalComptonFlag);
}
#endif

void Teton::constructSize()
{
   conduit::Node &options = getOptions();
   conduit::Node &node = options.fetch_existing("size");

   // These are required in the node.
   int nzones = node.fetch_existing("nzones").value();
   int ncornr = node.fetch_existing("ncornr").value();
   int nsides = node.fetch_existing("nsides").value();
   int nbelem = node.fetch_existing("nbelem").value();
   int maxcf = node.fetch_existing("maxcf").value();
   int maxCorner = node.fetch_existing("maxCorner").value();
   int ncomm = node.fetch_existing("ncomm").value();
   int ndim = node.fetch_existing("ndim").value();
   int ngr = options.fetch_existing("quadrature/num_groups").to_int(); // coerce from size_t or unsigned long
   int igeomToFortran = node.fetch_existing("geomType").value();

   // These are optional, and will be set to a default value.
   int enableNLTE = 0;
   int functionRNLTE = 2;
   int zoneBatchSize = 500;    // specific to CUDA BC solver.  Deprecate when OMP version of solver available.
   int nConcurrentBatches = 3; // specific to CUDA BC solver.  Deprecate when OMP version of solver available.
   double betaNLTE = 0.2;
   double gammaNLTE = 4.0;
   double radForceMultiplier = 1.0;
   bool dopplerShiftOn = true;
   bool usePWLD = false;
   bool useSurfaceMassLumping = false;
   bool useNewNonLinearSolver = false;
   bool useNewGTASolver = false;
   bool useGPU = false;
   bool useCUDASolver = false;
   bool useCUDASweep = false;

   if (!node.has_path("enableNLTE"))
   {
      node["enableNLTE"] = enableNLTE;
   }

   if (node.has_path("functionRNLTE"))
   {
      functionRNLTE = node.fetch_existing("functionRNLTE").value();
   }

   if (node.has_path("betaNLTE"))
   {
      betaNLTE = node.fetch_existing("betaNLTE").value();
   }

   if (node.has_path("gammaNLTE"))
   {
      betaNLTE = node.fetch_existing("gammaNLTE").value();
   }

   if (node.has_path("radForceMultiplier"))
   {
      radForceMultiplier = node.fetch_existing("radForceMultiplier").value();
   }

   if (node.has_path("DopplerShiftOn"))
   {
      dopplerShiftOn = node.fetch_existing("DopplerShiftOn").to_int();
   }

   if (node.has_path("usePWLD"))
   {
      usePWLD = node.fetch_existing("usePWLD").to_int();
   }

   if (node.has_path("useSurfaceMassLumping"))
   {
      useSurfaceMassLumping = node.fetch_existing("useSurfaceMassLumping").to_int();
   }

   if (node.has_path("useNewNonLinearSolver"))
   {
      useNewNonLinearSolver = node.fetch_existing("useNewNonLinearSolver").to_int();
   }

   if (node.has_path("useNewNonLinearSolver"))
   {
      useNewNonLinearSolver = node.fetch_existing("useNewNonLinearSolver").to_int();
   }

   if (node.has_path("useNewGTASolver"))
   {
      useNewGTASolver = node.fetch_existing("useNewGTASolver").to_int();
   }

   if (node.has_path("useGPU"))
   {
      useGPU = node.fetch_existing("useGPU").to_int();
   }

   if (node.has_path("useCUDASolver"))
   {
      useCUDASolver = node.fetch_existing("useCUDASolver").to_int();
   }

   if (node.has_path("useCUDASweep"))
   {
      useCUDASweep = node.fetch_existing("useCUDASweep").to_int();
   }

   //Temperature floor is already set, from constructDtControls
   double tfloor = options.fetch_existing("iteration/tfloor").value();

   teton_constructsize(&mRank,
                       &nzones,
                       &ncornr,
                       &nsides,
                       &nbelem,
                       &maxcf,
                       &maxCorner,
                       &ncomm,
                       &ndim,
                       &ngr,
                       &functionRNLTE,
                       &tfloor,
                       &radForceMultiplier,
                       &betaNLTE,
                       &gammaNLTE,
                       &dopplerShiftOn,
                       &useNewNonLinearSolver,
                       &useNewGTASolver,
                       &usePWLD,
                       &useSurfaceMassLumping,
                       &useGPU,
                       &useCUDASweep,
                       &useCUDASolver,
                       &zoneBatchSize,
                       &nConcurrentBatches,
                       &igeomToFortran);
}

void Teton::constructMemoryAllocator()
{
   conduit::Node &options = getOptions();

   int umpire_host_pinned_pool_allocator_id = -1;
   int umpire_device_pool_allocator_id = -1;

   if (options.has_path("memory_allocator/umpire_host_allocator_id"))
   {
      umpire_host_pinned_pool_allocator_id = options.fetch_existing("memory_allocator/umpire_host_allocator_id")
                                                .value();
   }
   if (options.has_path("memory_allocator/umpire_device_allocator_id"))
   {
      umpire_device_pool_allocator_id = options.fetch_existing("memory_allocator/umpire_device_allocator_id").value();
   }

   teton_constructmemoryallocator(&umpire_host_pinned_pool_allocator_id, &umpire_device_pool_allocator_id);
}

void Teton::dump(MPI_Comm communicator, std::string path)
{
// This is defined in conduit_relay_config.h
#if defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
   // NOTE: this routine saves the partitioned mesh given to Teton.
   conduit::Node &part = getMeshBlueprintPart();
   std::string file_protocol = "hdf5";
   conduit::relay::mpi::io::blueprint::save_mesh(part, path + "/blueprint_mesh", file_protocol, communicator);
#else
   std::cerr << " Teton: Unable to dump mesh blueprint viz file.  Conduit was not built with HDF5 support."
             << std::endl;
#endif
}
// ------------------------------------------------------------
//   step() - advance one cycle
// ------------------------------------------------------------

double Teton::step(int cycle)
{
#if defined(PARTITION_DEBUG)
   MPI_Barrier(mCommunicator);
   std::stringstream cs;
   cs << "Teton::step " << cycle;
   utilities::Banner b(mCommunicator, cs.str());
#endif
   conduit::Node &datastore = getDatastore();
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();
   conduit::Node &part = getMeshBlueprintPart();

   // TODO - These should be moved and made defaults in conduit node.
   int maxOSComptonChangeCorner = 1;
   double maxOSComptonChange = 0.0;

   if (options.has_path("iteration/outerMaxIt"))
   {
      int i = options.fetch_existing("iteration/outerMaxIt").value();
      teton_adjust_temperature_maxits(&i);
   }
   if (options.has_path("iteration/greyMaxIt"))
   {
      int i = options.fetch_existing("iteration/greyMaxIt").value();
      teton_adjust_grey_maxits(&i);
   }
   if (options.has_path("iteration/incidentFluxMaxIt"))
   {
      int i = options.fetch_existing("iteration/incidentFluxMaxIt").value();
      teton_adjust_fluxexchange_maxits(&i);
   }
   if (options.has_path("iteration/innerNLMaxIt"))
   {
      int i = options.fetch_existing("iteration/innerNLMaxIt").value();
      teton_adjust_nonlinear_maxits(&i);
   }

   // Global tolerance that gets set first, so others can override it.
   // This is the prefered option to set, and will set the others
   // consistently.
   if (options.has_path("iteration/relativeTolerance"))
   {
      double x = options.fetch_existing("iteration/relativeTolerance").value();
      teton_adjust_relative_tolerance(&x);
   }

   if (options.has_path("iteration/outerTempRelTol"))
   {
      double x = options.fetch_existing("iteration/outerTempRelTol").value();
      teton_adjust_temperature_reltol(&x);
   }
   if (options.has_path("iteration/outerPhiRelTol"))
   {
      double x = options.fetch_existing("iteration/outerPhiRelTol").value();
      teton_adjust_radenergydensity_reltol(&x);
   }
   if (options.has_path("iteration/incidentFluxRelTol"))
   {
      double x = options.fetch_existing("iteration/incidentFluxRelTol").value();
      teton_adjust_fluxexchange_reltol(&x);
   }
   if (options.has_path("iteration/innerNLRelTol"))
   {
      double x = options.fetch_existing("iteration/innerNLRelTol").value();
      teton_adjust_fluxexchange_reltol(&x);
   }
   if (options.has_path("iteration/greyRelTol"))
   {
      double x = options.fetch_existing("iteration/greyRelTol").value();
      teton_adjust_fluxexchange_reltol(&x);
   }

   // ------------------------------------------------------------
   // Update the mesh positions, material info, and opacity info
   // ------------------------------------------------------------

   // Updating the positions is not working when using the test_driver
   // on an mfem mesh.
   // At cycle two, the positions are wrong, and look like they are
   // either getting overwritten or the memory freed in the blueprint
   // node coordinate arrays.  This functionality works fine on
   // conduit mesh files.
   // For now, have the test_driver set the mesh_motion = 0 to skip
   // updating the positions, as we're only running static meshes
   // in our standalone meshes at present anyways.
   // -- black27
   //
   // Update zone vertex coordinates this cycle.  This should be done
   // after any mesh motion.
   //
   // Defaults to true, but can be disabled by setting
   // 'mesh_motion' = 0 in the options node.
   //
   // Calling updateMeshPositions will cause the volume difference from
   // last cycle to current cycle, that Teton tracks, to be updated.
   std::vector<std::string> updateFields;
   int mesh_motion = 1;
   if (options.has_path("mesh_motion"))
   {
      mesh_motion = options.fetch_existing("mesh_motion").value();
   }
   // Determine whether materials are present and required fields.
   int materials = 0;
   if (blueprint.has_path("fields/thermo_density/values"))
   {
      materials = 1;

      // Fields required in setMaterial.
      updateFields.push_back("thermo_density");
      updateFields.push_back("electron_specific_heat");
      updateFields.push_back("electron_temperature");
      updateFields.push_back("radiation_temperature");
      updateFields.push_back("electron_number_density");
      if (blueprint.has_path("fields/specific_energy_source"))
      {
         updateFields.push_back("specific_energy_source");
      }
   }
   int opacity = 0;
   if (blueprint.has_path("fields/absorption_opacity/values"))
   {
      opacity = 1;

      // Fields required in updateOpacities
      updateFields.push_back("absorption_opacity");
      updateFields.push_back("scattering_opacity");
   }

   // Update some fields, sending them through the partitioner from the blueprint
   // mesh to the partitioned mesh.
   std::string mainTopologyName(getMainTopology(part).name());
   sendFieldsOrig2Part(mainTopologyName, updateFields, mesh_motion == 1);

   // Now, do the updates using the partitiond data.
   if (mesh_motion)
   {
      // Partitioning is not necessary since if mesh coordinates were
      // updated, it happened above in sendFieldsOrig2Part.
      const bool dopartition = false;
      updateMeshPositions(dopartition);
   }
   //setMeshVelocity();
   // This updates the material properties (other than the opacities)

   // TODO Add something better than this to check for whether or not
   // new field values have been provided.
   if (materials)
   {
      setMaterials();
   }

   if (opacity)
   {
      updateOpacity();
   }

   // Set the time step information
   // A host code can either set these values in conduit, or can
   // setTimeStep() and that function will add these entries.
   double dtrad = options.fetch_existing("iteration/dtrad").value();
   double timerad = options.fetch_existing("iteration/timerad").value();
   double tfloor = options.fetch_existing("iteration/tfloor").value();

   // Update volumetric source profiles:
   mSourceManager.UpdateSources(timerad, dtrad);
   mSourceManager.UpdatePsiWithSources();

   // ------------------------------------------------------------
   // Run the step
   // ------------------------------------------------------------

   teton_settimestep(&cycle, &dtrad, &timerad, &tfloor);

   // Update cycle number in mesh blueprint.  This is used by conduit or Visit if we dump this mesh.
   blueprint["state/cycle"] = cycle;
   part["state/cycle"] = cycle;

   // If sanitizer level is provided:
   if (options.has_path("iteration/sanitizer/level"))
   {
      int num_bad_cats = checkInputSanity(options.fetch_existing("iteration/sanitizer"));
      options["iteration/sanitizer/num_bad_cats"] = num_bad_cats;

      if (num_bad_cats > 0)
      {
         // If we reached this point, we're letting the host code decide whether or not to kill the code.
         // We don't want to proceed with this time step because it'll most likely crash or hang, so we'll do nothing and return.
         // Host code should check the value of options["iteration/sanitizer/num_bad_cats"] after each step call.
         if (mRank == 0)
            std::cout << "Teton: Bad inputs found! Skipping time step." << std::endl;
         return mDTrad;
      }
   }

   // Main function in Teton to take a radiation step
   teton_radtr();

   // Update the radiation force (if the field is present)
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   // Note that, if the radiation force is present, there will always be
   // a z component in 2D or 3D ((r,z) or (x,y,z) coordinates). Also note
   // that we make this check on the original blueprint mesh since that is
   // where the host will have declared its field association.
   // TODO: fix when 1D is added
   mMapBackFields.clear();

   bool has_rad_force = blueprint.has_path("fields/radiation_force_z");
   std::string rad_force_type;
   bool elementAssociation = false;
   if (has_rad_force)
   {
      // The host provided radiation_force_z so we can check association.
      rad_force_type = blueprint.fetch_existing("fields/radiation_force_z/association").as_string();
      elementAssociation = rad_force_type == "element";
   }

   // The radiation_force_* fields are always computed so the getRadiationForceDensity()
   // method can work when partitioning is enabled.
   //
   // Create radiation_force_* fields on part mesh if they do not exist.
   // This does nothing if they already exist.
   createRadiationForceDensity(part, elementAssociation);
   if (elementAssociation)
   {
      updateZonalRadiationForce();
   }
   else
   {
      updateRadiationForce();
   }

   // Count number of zones in part topo.
   const conduit::Node &part_topo = getMainTopology(part);
   const int npart_zones = static_cast<int>(conduit::blueprint::mesh::utils::topology::length(part_topo));

   // Update the radiation energy deposited to the material.
   // Always compute this field so the getRadiationDeposited() method can work when
   // partitioning is enabled.
   //
   // Create field on the part mesh if it does not exist.
   createZonalField(part, mainTopologyName, FIELD_ELECTRON_ENERGY_DEPOSITED, npart_zones);

   std::string path(field_values(FIELD_ELECTRON_ENERGY_DEPOSITED));
   double *electron_energy_deposited = part.fetch_existing(path).value();
   getRadEnergyDeposited(electron_energy_deposited, npart_zones);
   mMapBackFields.push_back(FIELD_ELECTRON_ENERGY_DEPOSITED);

   // Update the radiation energy density
   // Check the blueprint mesh in case the host added this field.
   if (blueprint.has_path(field_values(FIELD_RADIATION_ENERGY_DENSITY)))
   {
      // During partitioning, FIELD_RADIATION_ENERGY_DENSITY would have been wrapped
      // as an mcarray due to it being a "multigroup" field. However, the partitioner
      // output would only be nzones in length. Since we're gathering data for mapback,
      // and Teton expects a large contiguous buffer, make sure it is large enough. This
      // should be ok since it is the partitioned mesh and we're immediately filling its
      // values from Teton.
      conduit::Node &red = part[field_path(FIELD_RADIATION_ENERGY_DENSITY)];
      if (doPartitioning())
      {
         const int ngr = options.fetch_existing("quadrature/num_groups").to_int();
         const auto expected_elements = static_cast<conduit::index_t>(npart_zones * ngr);
         if (red["values"].dtype().number_of_elements() < expected_elements)
         {
            red["values"].set(conduit::DataType::float64(expected_elements));
         }
         // If this is the first time filling out "red" then we may need to also set
         // some additional fields.
         red["association"] = "element";
         red["topology"] = mainTopologyName;
      }
      double *radiation_energy_density = red["values"].value();
      teton_getradiationenergydensity(radiation_energy_density);
      mMapBackFields.push_back(FIELD_RADIATION_ENERGY_DENSITY);
   }
#endif

   // Updates the counters and statistics needed by teton_printedits()
   // Also provides a copy of the corner temp field, but this use case has
   // been deprecated for a while.  For now, we pass Teton its own corner temp
   // field so we can still call this function.
   double *tec = datastore.fetch_existing("material/Tec").value();
   teton_rtedit(tec);

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   computeGenericSurfaceFluxTally();
#endif

   // Compute the recommended time step
   teton_dtnew(&maxOSComptonChangeCorner, &maxOSComptonChange);

   // put Teton's various edits in to its internal conduit node
   // This also puts the recommended timestep for the next iteration in mDTrad
   teton_publishedits(&mDTrad);

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   // Update the temperature result fields. NOTE: this had to come after teton_publishedits.

   // Create field on the part mesh if it does not exist.
   createZonalField(part, mainTopologyName, FIELD_RADIATION_TEMPERATURE, npart_zones);
   createZonalField(part, mainTopologyName, FIELD_MATERIAL_TEMPERATURE, npart_zones);

   double *radiation_temperature = part.fetch_existing(field_values(FIELD_RADIATION_TEMPERATURE)).value();
   getRadiationTemperature(radiation_temperature, npart_zones);
   mMapBackFields.push_back(FIELD_RADIATION_TEMPERATURE);

   double *material_temperature = part.fetch_existing(field_values(FIELD_MATERIAL_TEMPERATURE)).value();
   getMaterialTemperature(material_temperature, npart_zones);
   mMapBackFields.push_back(FIELD_MATERIAL_TEMPERATURE);
#endif

   // Migrate partition results to original mesh.
   sendFieldsPart2Orig(mainTopologyName, mMapBackFields);
   mMapBackFields.clear();

   double MatCoupTimeTotal, SweepTimeTotal, GPUSweepTimeTotal, GTATimeTotal;
   double RadtrTimeTotal, InitTimeTotal, FinalTimeTotal, timeNonRad = 0.0, timeOther = 0.0;

   teton_getrunstats(&MatCoupTimeTotal,
                     &SweepTimeTotal,
                     &GPUSweepTimeTotal,
                     &GTATimeTotal,
                     &RadtrTimeTotal,
                     &InitTimeTotal,
                     &FinalTimeTotal,
                     &timeNonRad,
                     &timeOther);

   options["iteration/dtcontrol/dtrec"] = mDTrad;
   int dtControlReason, dtControlProcess, dtControlZone;
   teton_getdtcontrolinfo(&dtControlReason, &dtControlProcess, &dtControlZone);
   options["iteration/dtcontrol/flag"] = dtControlReason;
   options["iteration/dtcontrol/process"] = dtControlProcess;
   options["iteration/dtcontrol/zone"] = dtControlZone;
   char *dtmessage_ptr;
   teton_getdtmessage(&dtmessage_ptr);
   std::string dtmsg = dtmessage_ptr;
   options["iteration/dtcontrol/message"] = dtmsg;

#if defined(PARTITION_DEBUG)
   MPI_Barrier(mCommunicator);
   bool testing = false;
   if (getenv("TETON_TESTING") != nullptr)
   {
      testing = atoi(getenv("TETON_TESTING")) > 0;
      if (testing)
      {
         // Test the results that have been computed in the cycle, store them.
         bool makeBaselines = getenv("TETON_TESTING_MAKE_BASELINES") != nullptr;
         const int flags = Test_RadiationForceDensity | Test_RadiationTemperature | Test_ReconstructPsi;
         conduit::Node n;
         const std::string fileBase = makeTestNode(n, getDatastore(), getMeshBlueprint(), getOptions(), flags);
         testing::test(n, fileBase, cycle, makeBaselines, mCommunicator);

         // Save the blueprint in a form we can look at in VisIt so we can compare baseline vs current.
         int verbose = getVerbose();
         if (verbose >= 2)
         {
            std::string name = makeBaselines ? "baseline" : "current";
            add_mcarray_fields(blueprint);
            conduit::relay::mpi::io::blueprint::save_mesh(blueprint, name, "hdf5", mCommunicator);
            remove_mcarray_fields(blueprint);

            if (doPartitioning())
            {
               std::string namep = makeBaselines ? "baseline_part" : "current_part";
               conduit::relay::mpi::io::blueprint::save_mesh(part, namep, "hdf5", mCommunicator);
            }
         }
      }
   }
#endif

   return mDTrad;
}

void Teton::constructEdits()
{
   conduit::Node &datastore = getDatastore();
   conduit::Node &options = getOptions();

   // Make a bunch of things zero before we get moving.
   datastore["rtedits/noutrt"] = 0;
   datastore["rtedits/ninrt"] = 0;
   datastore["rtedits/ngdart"] = 0;
   datastore["rtedits/nNLIters"] = 0;
   datastore["rtedits/maxNLIters"] = 0;
   datastore["rtedits/TrMaxZone"] = 0;
   datastore["rtedits/TeMaxZone"] = 0;
   datastore["rtedits/TrMaxProcess"] = 0;
   datastore["rtedits/TeMaxProcess"] = 0;
   datastore["rtedits/TeMax"] = 0.0;
   datastore["rtedits/TrMax"] = 0.0;
   // Total rad energy was already set when psi was initialized.
   datastore["rtedits/PowerIncident"] = 0.0;
   datastore["rtedits/PowerEscape"] = 0.0;
   datastore["rtedits/PowerAbsorbed"] = 0.0;
   datastore["rtedits/PowerEmitted"] = 0.0;
   datastore["rtedits/PowerExtSources"] = 0.0;
   datastore["rtedits/PowerCompton"] = 0.0;

   int ngr = options.fetch_existing("quadrature/num_groups").to_int(); //coerce from unsigned int or size_t

   // If boundary edits not provided, then create some defaults.
   // The Fortran code expects these arrays to be allocated externally, so hold them in the conduit tree.
   if (!options.has_path("boundary_edits"))
   {
      int numSpectrumAngleBins = 1;
      std::vector<double> spectrumAngleBinBoundaries{-1.0, 1.0};
      std::vector<double> RadPowerEscape(ngr, 0.0);
      std::vector<double> RadPowerIncident(ngr, 0.0);
      std::vector<double> PolarSectorPowerEscape(numSpectrumAngleBins * ngr, 0.0);

      options["boundary_edits/numSpectrumAngleBins"] = 1;
      options["boundary_edits/spectrumAngleBinBoundaryList"].set(spectrumAngleBinBoundaries);
      options["boundary_edits/RadPowerEscape"].set(RadPowerEscape);
      options["boundary_edits/RadPowerIncident"].set(RadPowerIncident);
      options["boundary_edits/PolarSectorPowerEscape"].set(PolarSectorPowerEscape);
   }

   int numSpectrumAngleBins = options.fetch_existing("boundary_edits/numSpectrumAngleBins").value();
   double *spectrumAngleBinBoundaryList_ptr = options.fetch_existing("boundary_edits/spectrumAngleBinBoundaryList")
                                                 .value();
   double *RadPowerEscape_ptr = options.fetch_existing("boundary_edits/RadPowerEscape").value();
   double *RadPowerIncident_ptr = options.fetch_existing("boundary_edits/RadPowerIncident").value();
   double *PolarSectorPowerEscape_ptr = options.fetch_existing("boundary_edits/PolarSectorPowerEscape").value();

   teton_constructeditor(&ngr,
                         &numSpectrumAngleBins,
                         &spectrumAngleBinBoundaryList_ptr[0],
                         &RadPowerEscape_ptr[0],
                         &RadPowerIncident_ptr[0],
                         &PolarSectorPowerEscape_ptr[0]);
}

// ------------------------------------------------------------
// computeGenericSurfaceFluxTally
//
// Constructs surface flux tallies
// ------------------------------------------------------------
void Teton::computeGenericSurfaceFluxTally()
{
   // The tally definition is split into two places.
   // The SURFACE information is in blueprint.
   // The other details of the tally (shape, groups, frame, etc.) live in options.
   conduit::Node &blueprint = getMeshBlueprint();
   conduit::Node &part = getMeshBlueprintPart();
   conduit::Node &options = getOptions();

   // NOTE: teton/surface_edits is created in TetonBlueprint::ProcessSurfaceEdits and
   //       all ranks will contain teton/surface_edits in the partitioned mesh, though
   //       some of the arrays may be empty if the rank had no faces.
   if (part.has_path("teton/surface_edits"))
   {
      conduit::Node &surface_edit_options_all = options.fetch_existing("surface_edits");
      conduit::Node &surface_edit_blueprint_all = part.fetch_existing("teton/surface_edits");
      conduit::NodeConstIterator surface_edit_blueprint_it = surface_edit_blueprint_all.children();
      while (surface_edit_blueprint_it.has_next())
      {
         const conduit::Node &surface_info = surface_edit_blueprint_it.next(); // surface info
         std::string surface_edit_name = surface_info.name();
         const conduit::Node &surface_edit_option = surface_edit_options_all[surface_edit_name]; // options for tallying

         int tmp = 0;
         const conduit::Node &corners = surface_info.fetch_existing("corners");
         const int num_corner_faces = corners.dtype().number_of_elements();
         const int *corners_ptr = (num_corner_faces > 0) ? corners.as_int_ptr() : &tmp;

         const conduit::Node &local_zone_faces = surface_info.fetch_existing("local_zone_faces");
         const int *local_zone_faces_ptr = (local_zone_faces.dtype().number_of_elements() > 0)
                                              ? local_zone_faces.as_int_ptr()
                                              : &tmp;

         const bool transform_to_lab_frame = surface_edit_option["transform_to_lab_frame"].as_int();
         const bool apply_time_shift = surface_edit_option["apply_time_shift"].as_int();

         // if not integrating over angles, then this is equal to the number of polar levels
         const bool integrate_over_angles = surface_edit_option["integrate_over_angles"].as_int();
         int num_angle_bins = 1;
         if (!integrate_over_angles)
            teton_getnumanglebins(&num_angle_bins);

         // if not integrating over all groups, then this is equal to the number of groups
         const bool integrate_over_all_groups = surface_edit_option["integrate_over_all_groups"].as_int();
         int num_groups = 1;
         if (!integrate_over_all_groups)
            num_groups = options["quadrature/num_groups"].as_int();

         const double *time_bin_boundaries = surface_edit_option["time_bin_boundaries"].as_double_ptr();
         const int num_time_bins = surface_edit_option["time_bin_boundaries"].dtype().number_of_elements() - 1;
         // TODO: make sure don't have to swap coordinate values in rz (since host code might
         // be in (z,r) coordinates
         const double *center_point = surface_edit_option.fetch_existing("center_point").as_double_ptr();
         const bool calculate_incident = surface_edit_option.fetch_existing("calculate_incident").as_int();
         const double scale_tally = surface_edit_option.fetch_existing("scale_tally").as_double();
         const bool calculate_error_metrics = surface_edit_option.fetch_existing("calculate_error_metrics").as_int();

         // The tally result arrays were passed in from the host code as fields on the
         // surface mesh, though it is not really a field for the surface mesh. Each
         // field is a non-spatial results array that is the same size on all ranks.
         double *tally = blueprint["fields/" + surface_edit_name + "_tallies/values"].as_double_ptr();
         double *tally_incident = nullptr;
         double *error_est_shift = nullptr;
         double *error_est_src_size = nullptr;
         if (calculate_incident)
         {
            tally_incident = blueprint["fields/" + surface_edit_name + "_tallies_incident/values"].as_double_ptr();
         }
         if (apply_time_shift && calculate_error_metrics)
         {
            error_est_shift = blueprint["fields/" + surface_edit_name + "_error_est_shift/values"].as_double_ptr();
            error_est_src_size = blueprint["fields/" + surface_edit_name + "_error_est_src_size/values"]
                                    .as_double_ptr();
         }

         // NOTE: This function performs global reductions on the output fields.
         teton_surfaceedit(&num_corner_faces,
                           &transform_to_lab_frame,
                           corners_ptr,
                           local_zone_faces_ptr,
                           &apply_time_shift,
                           center_point,
                           &num_angle_bins,
                           &num_groups,
                           &num_time_bins,
                           time_bin_boundaries,
                           &calculate_incident,
                           &scale_tally,
                           &calculate_error_metrics,
                           tally,
                           tally_incident,
                           error_est_shift,
                           error_est_src_size);
      }
   }
}

void Teton::dumpTallyToJson() const
{
   // The tally definition is split into two places.
   // The surface information associated with the tally is in blueprint.
   // The other details of the tally (shape, groups, frame, etc.) live in options.

   const conduit::Node &blueprint = getMeshBlueprint();
   const conduit::Node &options = getOptions();

   if (mRank == 0 && blueprint.has_path("teton/surface_edits"))
   {
      TetonSurfaceTallies::dumpTallyToJson(blueprint, options, mRank);
   }

   MPI_Barrier(mCommunicator);
} // end SnRad_dumpTally()

// ------------------------------------------------------------
// constructQuadrature
//
// Associate frequency group info with Teton (call once)
// ------------------------------------------------------------

void Teton::constructQuadrature()
{
   conduit::Node &options = getOptions();

   int nSets = options.fetch_existing("quadrature/nSets").value();
   int nSetsMaster = options.fetch_existing("quadrature/nSetsMaster").value();
   double *gnu_vals = options.fetch_existing("quadrature/gnu").value();
   int qtype = options.fetch_existing("quadrature/qtype").value();
   int qorder = options.fetch_existing("quadrature/qorder").value();
   int npolar = options.fetch_existing("quadrature/npolar").value();
   int nazimu = options.fetch_existing("quadrature/nazimu").value();
   int paxis = options.fetch_existing("quadrature/paxis").value();
   int ngr = options.fetch_existing("quadrature/num_groups").to_int(); //coerce from unsigned int or size_t
   int gtaOrder = options.fetch_existing("quadrature/gtaorder").value();
   int group = 0;
   std::vector<double> gnu;
   std::vector<int> quaddef;

   // Allocate memory for energy group bounds and quadrature definitions
   gnu.resize(ngr + 1);
   quaddef.resize(12);

   // Set energy group bounds
   for (group = 0; group < ngr; ++group)
   {
      gnu[group] = gnu_vals[group];
      gnu[group + 1] = gnu_vals[group + 1];
   }

   // Configure the quadrature for high-order Sn sweeps (all groups have the same quadrature)
   quaddef[0] = qtype;
   quaddef[1] = qorder;
   quaddef[2] = npolar;
   quaddef[3] = nazimu;
   quaddef[4] = paxis;
   quaddef[5] = -1; // Number of total angles (output). The Fortran populates this value.

   // Configure the quadrature for GTA sweeps ( reduced # angles and groups ).
   quaddef[6] = 1;
   quaddef[7] = gtaOrder;
   quaddef[8] = 1;
   quaddef[9] = 1;
   quaddef[10] = 1;
   quaddef[11] = -1; // Number of total angles (output).  The Fortran populates this.

   teton_constructquadrature_new(&nSetsMaster, &nSets, &quaddef[0], &gnu[0]);

   // Retrieve the # angles that the Fortran populated and add it to conduit so its accessible.
   // Just get the first group's value.  The Fortran performs verify checks to ensure the quadrature
   // definition is identical across all the groups, including the number of angles.
   // (excluding the grey acceleration set)
   int num_angles = quaddef[5];
   options["quadrature/num_angles"] = num_angles;
}

void Teton::resetSourceProfiles()
{
   conduit::Node &options = getOptions();

   if (!areSourceProfilesSet)
   {
      std::cerr << "setSourceProfiles must be called before resetSourceProfiles!" << std::endl;
      exit(1);
   }

   int nsrc = options.fetch_existing("boundary_conditions/num_source").value();

   for (int j = 0; j < nsrc; ++j)
   {
      std::string top = "sources/profile" + std::to_string(j + 1) + "/";
      int TetonProfileID = options.fetch_existing(top + "TetonProfileID").value();
      int NumTimes = options.fetch_existing(top + "NumTimes").value();
      int NumValues = options.fetch_existing(top + "NumValues").value();

      double *values_ptr = options.fetch_existing(top + "Values").value();
      std::vector<double> Values(NumValues);
      for (int k = 0; k < NumValues; ++k)
      {
         Values[k] = values_ptr[k];
      }

      if (NumTimes == 1 && NumValues == 1)
      {
         teton_resetprofile(&TetonProfileID, Values[0]);
      }
      else if (NumTimes == 1)
      {
         double Multiplier = options.fetch_existing(top + "Multiplier").value();
         teton_resetprofile(&TetonProfileID, &Multiplier, Values);
      }
      else
      {
         double Multiplier = options.fetch_existing(top + "Multiplier").value();

         double *times_ptr = options.fetch_existing(top + "Times").value();
         std::vector<double> Times(NumTimes);
         for (int k = 0; k < NumTimes; ++k)
         {
            Times[k] = times_ptr[k];
         }

         teton_resetprofile(&TetonProfileID, &NumTimes, &NumValues, &Multiplier, &Times[0], &Values[0]);
      }
   }
}

void Teton::setSourceProfiles()
{
   conduit::Node &options = getOptions();

   int nsrc_bdry = options.fetch_existing("boundary_conditions/num_source").value();

   for (int j = 0; j < nsrc_bdry; ++j)
   {
      std::string top = "sources/profile" + std::to_string(j + 1) + "/";
      int NumTimes = options.fetch_existing(top + "NumTimes").value();
      int NumValues = options.fetch_existing(top + "NumValues").value();

      double *values_ptr = options.fetch_existing(top + "Values").value();
      std::vector<double> Values(NumValues);
      for (int k = 0; k < NumValues; ++k)
      {
         Values[k] = values_ptr[k];
      }

      int TetonProfileID = -1;

      if (NumTimes == 1 && NumValues == 1)
      {
         teton_addprofile(Values[0], &TetonProfileID);
      }
      else if (NumTimes == 1)
      {
         double Multiplier = options.fetch_existing(top + "Multiplier").value();
         teton_addprofile(&Multiplier, Values, &TetonProfileID);
      }
      else
      {
         double *times_ptr = options.fetch_existing(top + "Times").value();
         std::vector<double> Times(NumTimes);
         for (int k = 0; k < NumTimes; ++k)
         {
            Times[k] = times_ptr[k];
         }

         double Multiplier = options.fetch_existing(top + "Multiplier").value();
         bool blackBody = options.fetch_existing(top + "blackBody")
                             .to_int(); // Conduit doesn't support a 'bool' data type.
         bool isotropic = options.fetch_existing(top + "isotropic")
                             .to_int(); // Conduit doesn't support a 'bool' data type.

         teton_addprofile(&NumTimes,
                          &NumValues,
                          &Multiplier,
                          &blackBody,
                          &isotropic,
                          &Times[0],
                          &Values[0],
                          &TetonProfileID);
      }

      // Save the TetonProfileID for later use:
      options[top + "TetonProfileID"] = TetonProfileID;
   }

   areSourceProfilesSet = true;

   if (!options.has_path("sources/interior_sources"))
      return;

   // TODO move this into TetonSources.cc?
   // all interior (volumetric?) sources:
   const conduit::Node &profiles_node = options["sources/profiles"];
   conduit::NodeConstIterator interior_sources_it = options["sources/interior_sources"].children();
   while (interior_sources_it.has_next())
   {
      const conduit::Node &src_node = interior_sources_it.next();
      std::string spatial_shape = src_node["spatial_shape"].as_string();
      std::string profile_str = src_node["profile"].as_string();
      if (spatial_shape == "point") // point source from tally
      {
         const conduit::Node &profile_node = profiles_node[profile_str];
         std::string profile_type = profile_node.fetch_existing("type").as_string();

         // const double* location = src_node["location"].as_double_ptr(); // TODO convert coordinate to zone index
         int source_rank = src_node.fetch_existing("rank").to_int(); // Rank that contains the point source
         int teton_zone_index = src_node.fetch_existing("zone_index").to_int();
         int teton_part_zone_index = -1;

         if (doPartitioning())
         {
            // The original rank and zone index were passed in for the source. The
            // zone index identifies the zone that contains the source point. If we
            // repartitioned, this could be a different rank and zone. We need to
            // let the new owner rank return the new zone index and have everyone
            // else return -1 for the zone index since they do not own it.
            //
            // NOTE: Any new sources that get implemented here would need to map
            //       their zone ids too to be compatible with partitioning.
            int origDomZone[2] = {source_rank, teton_zone_index};
            int partDomZone[2] = {-1, -1};
            zoneLookupOrig2Part(origDomZone, partDomZone);
            teton_part_zone_index = partDomZone[1];
         }
         else
         {
            teton_part_zone_index = (mRank == source_rank) ? teton_zone_index : -1;
         }

         double multiplier = 1.0;
         if (src_node.has_path("multiplier"))
            multiplier = src_node.fetch_existing("multiplier").to_double();

         if (profile_type == "tallyfile")
         {
            std::string tally_file_name = profile_node["filename"].as_string();
            std::string tally_name = profile_node["tallyname"].as_string();
            mSourceManager.AddPointSourceFromTally(teton_part_zone_index, tally_file_name, tally_name, multiplier);
         }
         // else if (profile_type == "isotropic")
         // { // TODO generic isotropic group-dependent
         // }
         else
         {
            std::cerr << "Unsupported source profile type " << profile_type << std::endl;
            exit(1);
         }
      }
      else
      {
         std::cerr << "Unsupported source spatial shape " << spatial_shape << std::endl;
         exit(1);
      }
   }
}

void Teton::zoneLookupOrig2Part(int originalDomZone[2], int partDomZone[2]) const
{
   if (doPartitioning())
   {
      const conduit::Node &part = getMeshBlueprintPart();
      const std::string mainTopologyName(getMainTopology(part).name());
      std::string vkey = "fields/" + mainTopologyName + "_original_element_ids/values";
      const conduit::Node &vnode = part.fetch_existing(vkey);
      const auto orig_domains = vnode.fetch_existing("domains").as_int_accessor();
      const auto orig_zones = vnode.fetch_existing("ids").as_int_accessor();
      const conduit::index_t n = orig_domains.number_of_elements();

      // Indicate not found.
      partDomZone[0] = -1;
      partDomZone[1] = -1;

      const int originalZone0 = originalDomZone[1] - 1;
      for (conduit::index_t i = 0; i < n; i++)
      {
         // Compare orig_zones against zero-origin zone number since the array
         // stores zero-origin zone ids.
         if (orig_domains[i] == originalDomZone[0] && orig_zones[i] == originalZone0)
         {
            // This rank contains the zone we're looking for. Return the new rank, zone index.
            partDomZone[0] = mRank;
            partDomZone[1] = static_cast<int>(i) + 1; // 1-origin zone
            break;
         }
      }
   }
   else
   {
      // Return the inputs as partitioning did not occur.
      partDomZone[0] = originalDomZone[0];
      partDomZone[1] = originalDomZone[1];
   }
}

void Teton::setMeshSizeAndPositions()
{
   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();

   int ndim = options.fetch_existing("size/ndim").value();
   int nzones = options.fetch_existing("size/nzones").value();

   if (ndim > 1)
   {
      double *zone_verts_ptr = part.fetch_existing("arrays/zone_verts").value();
      int *ncorners_ptr = part.fetch_existing("arrays/zone_to_ncorners").value();
      int ndim = options.fetch_existing("size/ndim").value();
      int maxCorner = options.fetch_existing("size/maxCorner").value();

      int off_set = 0;
      std::vector<double> zoneCoordinates(ndim * maxCorner);

      for (int zone = 0; zone < nzones; ++zone)
      {
         int zoneID = zone + 1;
         int ncorners = ncorners_ptr[zone];
         for (int c = 0; c < ncorners; ++c)
         {
            for (int i = 0; i < ndim; i++)
            {
               zoneCoordinates[ndim * c + i] = zone_verts_ptr[off_set];
               off_set += 1;
            }
         }
         teton_setnodeposition(&zoneID, &zoneCoordinates[0]);
      }
   }
   else
   {
      int nvertices = part["coordsets/coords/values/x"].dtype().number_of_elements();
      double *vertex_coords = part["coordsets/coords/values/x"].value();
      int nzones = nvertices - 1;
      std::vector<double> zoneCoordinates(2);
      for (int zone = 0; zone < nzones; ++zone)
      {
         int zoneID = zone + 1;
         zoneCoordinates[0] = vertex_coords[zone];
         zoneCoordinates[1] = vertex_coords[zone + 1];
         // TODO: add asset that zoneCoordinates[0] < zoneCoordinates[1]
         teton_setnodeposition(&zoneID, &zoneCoordinates[0]);
      }
   }
}

void Teton::setMeshVelocity()
{
   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();

   int nzones = options.fetch_existing("size/nzones").value();
   // TODO: change this to conform to blueprint standard
   double *velocities_ptr = part.fetch_existing("fields/velocity_at_corners").as_double_ptr();
   int *ncorners_ptr = part.fetch_existing("arrays/zone_to_ncorners").value();
   int ndim = options.fetch_existing("size/ndim").value();
   int maxCorner = options.fetch_existing("size/maxCorner").value();

   int off_set = 0;
   std::vector<double> velocity(ndim * maxCorner);

   for (int zone = 0; zone < nzones; ++zone)
   {
      int zoneID = zone + 1;
      int ncorners = ncorners_ptr[zone];
      for (int c = 0; c < ncorners; ++c)
      {
         for (int i = 0; i < ndim; i++)
         {
            velocity[ndim * c + i] = velocities_ptr[off_set];
            off_set += 1;
         }
      }
      teton_setnodevelocity(&zoneID, &velocity[0]);
   }

   return;
}

void Teton::setCommunication()
{
   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();

   int nsfaces;
   int *shared_faces_ptr = nullptr;

   if (part.has_path("shared_boundaries/nsfaces"))
   {
      nsfaces = part.fetch_existing("shared_boundaries/nsfaces").to_int();
      if (nsfaces > 0)
      {
         shared_faces_ptr = part.fetch_existing("shared_boundaries/shared_faces").value();
      }
   }
   else // if (options.has_path("shared_boundaries/nsfaces"))
   {
      // For backward compatbility
      nsfaces = options.fetch_existing("shared_boundaries/nsfaces").to_int();
      if (nsfaces > 0)
      {
         shared_faces_ptr = options.fetch_existing("shared_boundaries/shared_faces").value();
      }
   }
   int ndim = options.fetch_existing("size/ndim").value();

   int sface_offset = 0;
   for (int j = 0; j < nsfaces; ++j)
   {
      int bcID = shared_faces_ptr[sface_offset];
      int zoneID = shared_faces_ptr[sface_offset + 1];
      int faceLIDTeton = shared_faces_ptr[sface_offset + 2];
      int cid = shared_faces_ptr[sface_offset + 3];
      int cid2;
      if (ndim == 2)
      {
         cid2 = shared_faces_ptr[sface_offset + 4];
         sface_offset += 5;
      }
      else
      {
         sface_offset += 4;
      }

      teton_setsharedface(&bcID, &zoneID, &faceLIDTeton, &cid);
      if (ndim == 2)
      {
         teton_setsharedface(&bcID, &zoneID, &faceLIDTeton, &cid2);
      }
   }
}

void Teton::setMeshConnectivity()
{
   CALI_CXX_MARK_FUNCTION;

   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();

   int nzones = options.fetch_existing("size/nzones").value();

   std::string coord_type = part.fetch_existing("coordsets/coords/type").as_string();
   if (coord_type == "rectilinear")
   {
      // Teton expects two boundary conditions, one for zone1D == 1 and one for zoneID == nzones
      // Here BCZoneID[0] == 1, BCZoneID[1] == nzones or BCZoneID[0] == nzones and BCZoneID[1] == 1
      constexpr int numBCTotal = 2;
      int *BCZoneID = options.fetch_existing("boundary_conditions/zone_ids").value();
      for (int zone = 0; zone < nzones; ++zone)
      {
         int zoneID = zone + 1;
         teton_setzone1d(&zoneID, &numBCTotal, &BCZoneID[0]);
      }
   }
   else
   {
      int connect_off_set = 0;
      int *connectivity_ptr = part.fetch_existing("teton/arrays/corner_connectivity").value();
      for (int zone = 0; zone < nzones; ++zone)
      {
         int zoneID = connectivity_ptr[connect_off_set];
         int corner0 = connectivity_ptr[connect_off_set + 1];
         int zoneFaces = connectivity_ptr[connect_off_set + 2];
         int cornerFaces = connectivity_ptr[connect_off_set + 3];
         int zoneNCorner = connectivity_ptr[connect_off_set + 4];
         connect_off_set += 5;

         std::vector<int> zoneOpp(zoneFaces);
         std::vector<int> CornerID(cornerFaces);
         std::vector<int> CornerOpp(cornerFaces);
         std::vector<int> nCPerFace(zoneFaces);
         std::vector<int> FaceToBCList(zoneFaces);

         for (int j = 0; j < zoneFaces; ++j)
         {
            zoneOpp[j] = connectivity_ptr[connect_off_set + j];
         }
         connect_off_set += zoneFaces;

         for (int j = 0; j < cornerFaces; ++j)
         {
            CornerID[j] = connectivity_ptr[connect_off_set + j];
         }
         connect_off_set += cornerFaces;

         for (int j = 0; j < cornerFaces; ++j)
         {
            CornerOpp[j] = connectivity_ptr[connect_off_set + j];
         }
         connect_off_set += cornerFaces;

         for (int j = 0; j < zoneFaces; ++j)
         {
            nCPerFace[j] = connectivity_ptr[connect_off_set + j];
         }
         connect_off_set += zoneFaces;

         for (int j = 0; j < zoneFaces; ++j)
         {
            FaceToBCList[j] = connectivity_ptr[connect_off_set + j];
         }
         connect_off_set += zoneFaces;

         teton_setzone(&zoneID,
                       &corner0,
                       &zoneFaces,
                       &cornerFaces,
                       &zoneNCorner,
                       &zoneOpp[0],
                       &CornerID[0],
                       &CornerOpp[0],
                       &nCPerFace[0],
                       &FaceToBCList[0]);
      }

      part.remove("teton/arrays/corner_connectivity");
   }
}

void Teton::setMaterials()
{
   conduit::Node &datastore = getDatastore();
   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();

   int nzones = options.fetch_existing("size/nzones").value();

   double *density_ptr = part.fetch_existing("fields/thermo_density/values").value();
   double *cv_ptr = part.fetch_existing("fields/electron_specific_heat/values").value();
   double *tez_ptr = part.fetch_existing("fields/electron_temperature/values").value();
   double *trz_ptr = part.fetch_existing("fields/radiation_temperature/values").value();
   double *nez_ptr = part.fetch_existing("fields/electron_number_density/values").value();

   // Really the effective electron specific energy source.
   if (part.has_path("fields/specific_energy_source"))
   {
      double *matSource = part.fetch_existing("fields/specific_energy_source/values").value();
      teton_setmaterialsource(matSource);
   }

   // Initialize arrays to handle multi-material zones
   double *tec = datastore.fetch_existing("material/Tec").value();
   teton_initmaterial(tec);

   double scm = 1.;
   if (options.has_path("compton/stim_compton_mult"))
   {
      scm = options.fetch_existing("compton/stim_compton_mult").value();
   }

   for (int zone = 0; zone < nzones; ++zone)
   {
      int zoneID = zone + 1;
      double rho = density_ptr[zone];
      double cv = cv_ptr[zone];
      double tez = tez_ptr[zone];
      double trz = trz_ptr[zone];
      double nez = nez_ptr[zone];
      teton_setmaterial(&zoneID, &cv, &rho, &tez, &trz, &nez, &scm);
   }
}

void Teton::updateOpacity()
{
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();

   // NOTE: nzones in this case is on the part mesh.
   conduit::index_t nzones = options.fetch_existing("size/nzones").to_index_t();
   conduit::index_t ngroups = options.fetch_existing("quadrature/num_groups").to_index_t();

   bool useInternalSigmaS = false;
   if (options.has_path("compton/use_internal_sigma_s"))
   {
      useInternalSigmaS = options.fetch_existing("compton/use_internal_sigma_s").as_int();
   }
   bool useTableSigmaS = (not useInternalSigmaS);

   // zero out opacities
   teton_initopacity();

   // These fields would have been turned into mcarrays if partitioning. We
   // get them as though they might be mcarrays and then use NDAccessor to
   // get the (zone,ig) element data. Note though that we are getting them
   // from the part mesh and even if the mcarrays do not exist on the blueprint
   // mesh at present, they will be on the part mesh.
   auto n_absorption_opacity = const_cast<conduit::Node &>(fetch_mcarray(part, "absorption_opacity"));
   std::vector<utilities::NDDimension> dims{{"zone", nzones}, {"group", ngroups}};
   utilities::NDAccessor absorption_opacity(n_absorption_opacity["values"], dims, doInterleave("absorption_opacity"));
   std::vector<double> siga_loc(ngroups, 0), sigs_loc(ngroups, 0);
   if (useTableSigmaS)
   {
      auto n_scattering_opacity = const_cast<conduit::Node &>(fetch_mcarray(part, "scattering_opacity"));
      utilities::NDAccessor scattering_opacity(n_scattering_opacity["values"],
                                               dims,
                                               doInterleave("scattering_opacity"));
      for (conduit::index_t zone = 0; zone < nzones; zone++)
      {
         for (conduit::index_t ig = 0; ig < ngroups; ig++)
         {
            std::vector<conduit::index_t> idx{zone, ig};
            siga_loc[ig] = absorption_opacity(idx);
            sigs_loc[ig] = scattering_opacity(idx);
         }
         int zoneID = zone + 1;
         teton_setopacity(&zoneID, &siga_loc[0], &sigs_loc[0], &useTableSigmaS);
      }
   }
   else
   {
      for (conduit::index_t zone = 0; zone < nzones; zone++)
      {
         for (conduit::index_t ig = 0; ig < ngroups; ig++)
         {
            siga_loc[ig] = absorption_opacity(std::vector<conduit::index_t>{zone, ig});
            sigs_loc[ig] = 0.;
         }
         int zoneID = zone + 1;
         teton_setopacity(&zoneID, &siga_loc[0], &sigs_loc[0], &useTableSigmaS);
      }
   }

   if (not useTableSigmaS)
   {
      teton_setscatteringopacity(&mInternalComptonFlag);
   }
#endif
}

void Teton::constructIterationControls()
{
   // These are constructed with default values.
   // TODO - investigate pulling all the default values out of the Fortran and up into the C++, then
   // passing them into the older Fortran API.  Like constructDtControls.
   //
   // That would require codes using the older API to switch from
   // teton_constructitercontrols
   // to
   // TetonConduitInterface::constructIterationControls
   // -- Aaron
   teton_constructitercontrols();
}

void Teton::constructDtControls()
{
   conduit::Node &options = getOptions();

   // Default values for dt controls.
   int cycle = 0;
   double dtrad = 0.001;
   double dtrmn = 1e-40;
   double dtrmx = 0.1;
   double delte = 0.4;
   double deltr = 0.4;
   double tfloor = 1.0e-5;
   double timerad = 0.0;

   if (options.has_path("iteration/dtrad"))
   {
      dtrad = options.fetch_existing("iteration/dtrad").value();
   }
   else
   {
      options["iteration/dtrad"] = dtrad;
   }

   if (options.has_path("iteration/dtrmn"))
   {
      dtrmn = options.fetch_existing("iteration/dtrmn").value();
   }
   else
   {
      options["iteration/dtrmn"] = dtrmn;
   }

   if (options.has_path("iteration/dtrmx"))
   {
      dtrmx = options.fetch_existing("iteration/dtrmx").value();
   }
   else
   {
      options["iteration/dtrmx"] = dtrmx;
   }

   if (options.has_path("iteration/delte"))
   {
      delte = options.fetch_existing("iteration/delte").value();
   }
   else
   {
      options["iteration/delte"] = delte;
   }

   if (options.has_path("iteration/deltr"))
   {
      deltr = options.fetch_existing("iteration/deltr").value();
   }
   else
   {
      options["iteration/deltr"] = deltr;
   }

   teton_constructdtcontrols(&dtrad, &dtrmn, &dtrmx, &delte, &deltr);

   // These used later in the step() function, along with dtrad.
   if (!options.has_path("iteration/cycle"))
   {
      options["iteration/cycle"] = cycle;
   }

   if (!options.has_path("iteration/timerad"))
   {
      options["iteration/timerad"] = timerad;
   }

   if (!options.has_path("iteration/tfloor"))
   {
      options["iteration/tfloor"] = tfloor;
   }
}

// ---------------------------------------------------------------------------
// Function pertaining to checkpoints/restarts
// ---------------------------------------------------------------------------
conduit::Node &Teton::getCheckpoint()
{
   return *teton_conduitcheckpoint_get_cptr();
}

void Teton::checkpointPrepareForLoad()
{
   teton_conduitcheckpoint_prep_for_load();
}

void Teton::checkpointPrepareForSave()
{
   teton_conduitcheckpoint_prep_for_save();
}

void Teton::checkpointDataLoaded()
{
   teton_conduitcheckpoint_data_loaded();
}

void Teton::checkpointExternalDataLoaded()
{
   teton_conduitcheckpoint_external_data_loaded();
   // This may not be necessary, but we'll put it in here for now.
   // It updates the zonal electron and rad temperatures from the electron corner temps.
   conduit::Node &datastore = getDatastore();
   double *tec = datastore.fetch_existing("material/Tec").value();
   teton_rtedit(tec);
}

void Teton::checkpointFinished()
{
   // not implemented
   teton_conduitcheckpoint_teardown();

   // clear the node tree here instead of in fortran.
   // TODO: wrap this function so fortran can call it instead in 'teardown'
   conduit::Node &node = getCheckpoint();
   node.reset();
}

conduit::Node &Teton::getDatastore()
{
   return *teton_get_datastore_cptr();
}

const conduit::Node &Teton::getDatastore() const
{
   return *teton_get_datastore_cptr();
}

double Teton::getRadiationTemperature(int zone) const
{
   // This used to call teton_getradiationtemperature directly but we get the
   // results from the field to support partitioned meshes.

   const conduit::Node &blueprint = getMeshBlueprint();
   auto acc = blueprint.fetch_existing(field_values(FIELD_RADIATION_TEMPERATURE)).as_double_accessor();
   int zone0 = zone - 1;
   return acc[zone0];
}

double Teton::getMaterialTemperature(int zone) const
{
   // This used to call teton_getmaterialtemperature directly but we get the
   // results from the field to support partitioned meshes.

   const conduit::Node &blueprint = getMeshBlueprint();
   auto acc = blueprint.fetch_existing(field_values(FIELD_MATERIAL_TEMPERATURE)).as_double_accessor();
   int zone0 = zone - 1;
   return acc[zone0];
}

double Teton::getRadiationDeposited(int zone) const
{
   // This used to call teton_getradiationdeposited directly but we get the
   // results from the field to support partitioned meshes.

   const conduit::Node &blueprint = getMeshBlueprint();
   auto acc = blueprint.fetch_existing(field_values(FIELD_ELECTRON_ENERGY_DEPOSITED)).as_double_accessor();
   int zone0 = zone - 1;
   return acc[zone0];
}

void Teton::setTimeStep(int cycle, double dtrad, double timerad)
{
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();
   conduit::Node &part = getMeshBlueprintPart();

   options["iteration/cycle"] = cycle;
   options["iteration/dtrad"] = dtrad;
   options["iteration/timerad"] = timerad;

   // Used by conduit or Visit if mesh is dumped for viz purposes.
   blueprint["state/cycle"] = cycle;
   part["state/cycle"] = cycle;
}

void Teton::updateMeshPositions()
{
   // This method is public it gets called by client codes. We need to ensure
   // that the part mesh gets its coordinates updated from the blueprint mesh.
   const bool doPartition = true;
   updateMeshPositions(doPartition);
}

void Teton::updateMeshPositions(bool doPartition)
{
   CALI_CXX_MARK_FUNCTION;

   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();

   // If we're partitioning, update the coordinates.
   if (doPartition)
   {
      std::string mainTopologyName(getMainTopology(part).name());
      sendFieldsOrig2Part(mainTopologyName, std::vector<std::string>{}, true);
   }

   int nzones = mZoneToNCorners.size();
   int ndim = options.fetch_existing("size/ndim").value();

   // Best practice for updating the mesh zone vertices is for codes to update the blueprint coordinate arrays.
   //
   // MARBL is providing a zone_verts array directly, at the moment.  This array is a listing of the zone vertices,
   // in the same order as the corners in each zone.  Use that if it is present, otherwise generate it from the
   // blueprint coords.
   //
   // NOTE: If doing partitioning then there will be no arrays/zone_verts since
   //       it would have been supplied in the "blueprint" node instead of the
   //       "blueprint_partition" node. Thus, we'll make the node using data in
   //       the part mesh.
   if (!part.has_path("arrays/zone_verts"))
   {
      int corner_counter = 0;
      int zoneVertsSize = 0;

      const double *m_x = nullptr;
      const double *m_y = nullptr;
      const double *m_r = nullptr;
      const double *m_z = nullptr;

      if (ndim == 1)
      {
         m_x = part.fetch_existing("coordsets/coords/values/x").value();
      }
      else if (ndim == 2)
      {
         if (part.has_path("coordsets/coords/values/r"))
         {
            m_r = part.fetch_existing("coordsets/coords/values/r").value();
         }
         else
         { // assuming zr ordering for marbl/ares as fallback.  EVERYONE should just specify r and z directly.
            m_r = part.fetch_existing("coordsets/coords/values/y").value();
         }
         if (part.has_path("coordsets/coords/values/z"))
         {
            m_z = part.fetch_existing("coordsets/coords/values/z").value();
         }
         else
         { // assuming zr ordering for marbl/ares as fallback.  EVERYONE should just specify r and z directly.  For now, issue a warning.
            m_z = part.fetch_existing("coordsets/coords/values/x").value();
         }
      }
      else if (ndim == 3)
      {
         m_x = part.fetch_existing("coordsets/coords/values/x").value();
         m_y = part.fetch_existing("coordsets/coords/values/y").value();
         m_z = part.fetch_existing("coordsets/coords/values/z").value();
      }
      else
      {
         std::cerr << "Invalid number of dimensions." << std::endl;
         exit(1);
      }

      // loop over mesh elements
      for (int zone = 0; zone < nzones; ++zone)
      {
         int ncorners = mZoneToNCorners[zone];
         zoneVertsSize += ncorners * ndim;
      }

      std::vector<double> zoneVerts;
      zoneVerts.reserve(zoneVertsSize);
      for (int zone = 0; zone < nzones; ++zone)
      {
         int ncorners = mZoneToNCorners[zone];
         // loop over Teton corners in the element
         for (int c = 0; c < ncorners; ++c)
         {
            // get the index of the vertex in the coord array, corresponding to this corner
            int corner = mZoneToCorners[corner_counter];
            int v = mCornerToVertex[corner];

            // store the vertex coordinates
            if (ndim == 1)
            {
               zoneVerts.push_back(m_x[v]);
            }
            else if (ndim == 2)
            {
               zoneVerts.push_back(m_r[v]);
               zoneVerts.push_back(m_z[v]);
            }
            else
            {
               zoneVerts.push_back(m_x[v]);
               zoneVerts.push_back(m_y[v]);
               zoneVerts.push_back(m_z[v]);
            }
            corner_counter += 1;
         }
      }
      part["arrays/zone_verts"].set(zoneVerts.data(), zoneVerts.size());
   }

   setMeshSizeAndPositions();

   // Update Teton geometry
   teton_getvolume();

   // We're done updating the node positions, we shouldn't need zone_verts anymore.
   part.remove("arrays/zone_verts");

   return;
}

const std::vector<std::string> &Teton::radiationForceDensityFields() const
{
   return mRadiationForceDensityFields;
}

void Teton::initializeRadiationForceDensityFieldNames()
{
   // Get the number of dimensions from the blueprint mesh since it might not be
   // in the options yet.
   const conduit::Node &blueprint = getMeshBlueprint();
   std::string csname(getMainTopology(blueprint).fetch_existing("coordset").as_string());
   const conduit::Node &coordset = blueprint.fetch_existing("coordsets/" + csname);
   const int ndim = static_cast<int>(conduit::blueprint::mesh::coordset::dims(coordset));

   mRadiationForceDensityFields.clear();
   mRadiationForceDensityFields.reserve(ndim);
   if (ndim == 1)
   {
      mRadiationForceDensityFields.emplace_back(FIELD_RADIATION_FORCE_X);
   }
   else if (ndim == 2)
   {
      mRadiationForceDensityFields.emplace_back(FIELD_RADIATION_FORCE_Z);
      mRadiationForceDensityFields.emplace_back(FIELD_RADIATION_FORCE_R);
   }
   else if (ndim == 3)
   {
      mRadiationForceDensityFields.emplace_back(FIELD_RADIATION_FORCE_X);
      mRadiationForceDensityFields.emplace_back(FIELD_RADIATION_FORCE_Y);
      mRadiationForceDensityFields.emplace_back(FIELD_RADIATION_FORCE_Z);
   }
}

std::vector<double *> Teton::radiationForceDensity(conduit::Node &root) const
{
   std::vector<double *> ptrs;
   const auto names = radiationForceDensityFields();
   for (const auto &name : names)
   {
      double *d = root.fetch_existing(field_values(name)).value();
      ptrs.push_back(d);
   }
   return ptrs;
}

void Teton::createRadiationForceDensity(conduit::Node &root, bool elementAssociation)
{
   // NOTE: If we actually create these force fields then it could cause
   //       updateRadiationForce() to be called when it otherwise might
   //       not have been, as in the case where there were no force fields.

   // Make sure the radiation_force_ paths exist on the blueprint mesh since
   // they will be queried in getRadiationForceDensity().
   std::string mainTopologyName(getMainTopology(root).name());
   const auto names = radiationForceDensityFields();
   conduit::index_t nnodes = 0, nzones = 0;
   for (const auto &name : names)
   {
      // Determine number of nodes.
      if (nnodes == 0)
      {
         const conduit::Node &topo = root.fetch_existing("topologies/" + mainTopologyName);
         std::string csname(topo.fetch_existing("coordset").as_string());
         const conduit::Node &coordset = root.fetch_existing("coordsets/" + csname);
         nnodes = conduit::blueprint::mesh::coordset::length(coordset);
         nzones = conduit::blueprint::mesh::topology::length(topo);
      }

      const auto path = field_path(name);
      if (!root.has_path(path) && nnodes > 0 && nzones > 0)
      {
         conduit::index_t nvalues = elementAssociation ? nzones : nnodes;
         conduit::Node &f = root[path];
         f["association"] = elementAssociation ? "element" : "vertex";
         f["topology"] = mainTopologyName;
         f["values"].set(conduit::DataType::float64(nvalues));
         memset(f["values"].as_float64_ptr(), 0, nvalues * sizeof(conduit::float64));
      }
   }

   // Create a vertex field for the corner volume sums. These get used in
   // getRadiationForceDensity.
   if (!elementAssociation && !root.has_path(field_path(FIELD_CORNER_VOLUME_SUMS)) && nnodes > 0)
   {
      conduit::Node &f = root[field_path(FIELD_CORNER_VOLUME_SUMS)];
      f["association"] = "vertex";
      f["topology"] = mainTopologyName;
      f["values"].set(conduit::DataType::float64(nnodes));
      memset(f["values"].as_float64_ptr(), 0, nnodes * sizeof(conduit::float64));
   }
}

void Teton::SumSharedNodalValues(conduit::Node &root, double *nodal_field)
{
   const conduit::Node &options = getOptions();

   if (root.has_path("adjsets"))
   {
      int ndim = options.fetch_existing("size/ndim").value();
      std::string adjset_name = ndim > 1 ? "adjsets/main_adjset" : "adjsets/mesh";
      const conduit::Node &vertex_adjset = root[adjset_name];
      conduit::NodeConstIterator groups_it = vertex_adjset["groups"].children();
      const int num_vertex_groups = vertex_adjset["groups"].number_of_children();
      const int num_vertices = root.fetch_existing("coordsets/coords/values/x").dtype().number_of_elements();

      while (groups_it.has_next())
      {
         const conduit::Node &vertex_group = groups_it.next();
         const auto group_neighbors = vertex_group.fetch_existing("neighbors").as_int_accessor();
         const auto group_vertices = vertex_group.fetch_existing("values").as_int_accessor();
         const int num_neighbors = static_cast<int>(group_neighbors.number_of_elements());
         const int num_vertices = static_cast<int>(group_vertices.number_of_elements());

         std::vector<MPI_Request> requests_vec(2 * num_neighbors);
         MPI_Request *send_requests = requests_vec.data();
         MPI_Request *recv_requests = requests_vec.data() + num_neighbors;
         std::vector<MPI_Status> statuses_vec(num_neighbors);
         MPI_Status *statuses = statuses_vec.data();

         std::vector<std::vector<double>> fields_to_send(num_neighbors);
         std::vector<std::vector<double>> fields_to_recv(num_neighbors);
         for (int vn = 0; vn < num_neighbors; ++vn)
         {
            fields_to_send[vn].resize(num_vertices);
            fields_to_recv[vn].resize(num_vertices);
         }

         for (int vn = 0; vn < num_neighbors; ++vn)
         {
            const int nbr_rank = group_neighbors[vn];

            for (int j = 0; j < num_vertices; ++j)
            {
               const int vid = group_vertices[j];
               fields_to_send[vn][j] = nodal_field[vid];
            }

            int tag = 0;
            MPI_Isend(&fields_to_send[vn][0],
                      num_vertices,
                      MPI_DOUBLE,
                      nbr_rank,
                      tag,
                      mCommunicator,
                      &send_requests[vn]);
            MPI_Irecv(&fields_to_recv[vn][0],
                      num_vertices,
                      MPI_DOUBLE,
                      nbr_rank,
                      tag,
                      mCommunicator,
                      &recv_requests[vn]);
         }
         MPI_Waitall(num_neighbors, send_requests, statuses);
         MPI_Waitall(num_neighbors, recv_requests, statuses);

         // Add neighboring contributions to nodal field
         for (int vn = 0; vn < num_neighbors; ++vn)
         {
            for (int j = 0; j < num_vertices; ++j)
            {
               const int vid = group_vertices[j];
               nodal_field[vid] += fields_to_recv[vn][j];
            }
         }
      }
   }
}

// NOTE: the Vectors RadiationForceXTotal, ..., must
//       already be sized to the number of mesh vertices
void Teton::getRadiationForceDensity1D(double *RadiationForceDensityX)
{
   getRadiationForceDensity(RadiationForceDensityX, nullptr, nullptr);
}

// NOTE: the Vectors RadiationForceXTotal, ..., must
//       already be sized to the number of mesh vertices
void Teton::getRadiationForceDensity(double *RadiationForceDensityX,
                                     double *RadiationForceDensityY,
                                     double *RadiationForceDensityZ)
{
   // The data arrays we're copying from in the blueprint fields were updated
   // in updateRadiationForce() during step().
   const conduit::Node &blueprint = getMeshBlueprint();
   const auto fieldNames = radiationForceDensityFields();
   double *dest[] = {RadiationForceDensityX, RadiationForceDensityY, RadiationForceDensityZ};
   conduit::index_t n{};
   for (size_t c = 0; c < 3; c++)
   {
      if (c < fieldNames.size())
      {
         const conduit::Node &n_cvs = blueprint.fetch_existing(field_values(FIELD_CORNER_VOLUME_SUMS));
         const conduit::Node &n_comp = blueprint.fetch_existing(field_values(fieldNames[c]));
         const auto cvs = n_cvs.as_double_accessor();
         const auto acc = n_comp.as_double_accessor();
         n = acc.number_of_elements();
         for (conduit::index_t i = 0; i < n; i++)
            dest[c][i] = acc[i] / cvs[i];
      }
      else if (dest[c] != nullptr)
      {
         // Zero out this component. Relies on previous iteration setting n.
         memset(dest[c], 0, n * sizeof(double));
      }
   }
}

void Teton::updateRadiationForce()
{
   CALI_CXX_MARK_FUNCTION;

   // Compute the radiation force internally in Teton
   // for each zone and corner
   teton_setradiationforce();

   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();
   int ndim = options.fetch_existing("size/ndim").value();
   int maxCorner = options.fetch_existing("size/maxCorner").value();
   maxCorner = std::max(maxCorner, 2);
   int nzones = options.fetch_existing("size/nzones").value();
   int nverts = options.fetch_existing("size/nverts").value();
   std::vector<double> RadiationForce(ndim * maxCorner, 0.);
   std::vector<double> CornerVolumes(maxCorner, 0.);
   int corner_counter = 0;

   const auto fieldNames = radiationForceDensityFields();
   auto radiationForce = radiationForceDensity(part);
   conduit::Node &n_cvs = part.fetch_existing(field_values(FIELD_CORNER_VOLUME_SUMS));
   const conduit::index_t nvalues = n_cvs.dtype().number_of_elements();

   // We need to map these fields back to the blueprint mesh because they are
   // queried as results.
   for (const auto &f : fieldNames)
      mMapBackFields.push_back(f);
   mMapBackFields.push_back(FIELD_CORNER_VOLUME_SUMS);

   // Zero out the radiation force components.
   for (auto &ptr : radiationForce)
      memset(ptr, 0, nvalues * sizeof(double));
   auto nc = static_cast<int>(radiationForce.size());

   // Zero out the corner volume sums. These are needed for getRadiationForceDensity
   // to return the right values d=m/v.
   double *cornerVolumeSums = n_cvs.value();
   memset(cornerVolumeSums, 0, nvalues * sizeof(double));

   if (ndim == 1)
   {
      // Note: 1D does not involve certain mapping arrays.
      for (int zone = 0; zone < nzones; ++zone)
      {
         // Get the radiation force and volume on each corner of each zone
         int zoneID = zone + 1;
         teton_getradiationforce(&zoneID, &RadiationForce[0]);
         teton_getcornervolumes(&zoneID, &CornerVolumes[0]);

         int v1 = zone;
         int v2 = zone + 1;
         radiationForce[0][v1] += RadiationForce[0];
         radiationForce[0][v2] += RadiationForce[1];
         cornerVolumeSums[v1] += CornerVolumes[0];
         cornerVolumeSums[v2] += CornerVolumes[1];
      }
   }
   else
   {
      for (int zone = 0; zone < nzones; ++zone)
      {
         // Get the radiation force and volume on each corner of each zone
         int zoneID = zone + 1;
         teton_getradiationforce(&zoneID, &RadiationForce[0]);
         teton_getcornervolumes(&zoneID, &CornerVolumes[0]);

         // Average the radiation force around vertices
         int ncorners = mZoneToNCorners[zone];
         for (int c = 0; c < ncorners; ++c)
         {
            int cornerID = mZoneToCorners[corner_counter];
            int vertexID = mCornerToVertex[cornerID];
            corner_counter += 1;

            for (int comp = 0; comp < nc; comp++)
               radiationForce[comp][vertexID] += RadiationForce[c * ndim + comp];

            cornerVolumeSums[vertexID] += CornerVolumes[c];
         }
      }
   }

   // Sum shared vertex values across processor domains.
   SumSharedNodalValues(part, cornerVolumeSums);
   for (double *forceComponent : radiationForce)
      SumSharedNodalValues(part, forceComponent);
}

void Teton::updateZonalRadiationForce()
{
   CALI_CXX_MARK_FUNCTION;

   // Compute the radiation force internally in Teton
   // for each zone and corner
   teton_setradiationforce();

   conduit::Node &options = getOptions();
   conduit::Node &part = getMeshBlueprintPart();
   int ndim = options.fetch_existing("size/ndim").value();
   int maxCorner = options.fetch_existing("size/maxCorner").value();
   int nzones = options.fetch_existing("size/nzones").value();
   std::vector<double> RadiationForce(ndim * maxCorner);

   const auto fieldNames = radiationForceDensityFields();
   auto radiationForce = radiationForceDensity(part);

   // We need to map these fields back to the blueprint mesh because they are
   // queried as results.
   for (const auto &f : fieldNames)
      mMapBackFields.push_back(f);

   // Zero out the radiation force components.
   for (auto &ptr : radiationForce)
      memset(ptr, 0, nzones * sizeof(double));
   auto nc = static_cast<int>(radiationForce.size());

   for (int zone = 0; zone < nzones; ++zone)
   {
      // Get the radiation force and volume on each corner of each zone
      int zoneID = zone + 1;
      teton_getradiationforce(&zoneID, &RadiationForce[0]);
      // Sum the corner radiation force
      int ncorners = mZoneToNCorners[zone];
      for (int c = 0; c < ncorners; ++c)
      {
         for (int comp = 0; comp < nc; comp++)
            radiationForce[comp][zone] += RadiationForce[c * ndim + comp];
      }
   }
}

void Teton::getRadEnergyDeposited(double *RadEnergyDeposited, int nzones) const
{
   for (int zone = 0; zone < nzones; ++zone)
   {
      // Get the radiation energy deposited
      double rad_temp;
      int zoneID = zone + 1;
      teton_getradiationdeposited(&zoneID, &RadEnergyDeposited[zone], &rad_temp);
   }
}

conduit::Node &Teton::getMainTopology(conduit::Node &root)
{
   conduit::Node &topologies = root.fetch_existing("topologies");
   return topologies.child(0);
}

const conduit::Node &Teton::getMainTopology(const conduit::Node &root) const
{
   const conduit::Node &topologies = root.fetch_existing("topologies");
   return topologies.child(0);
}

void Teton::createZonalField(conduit::Node &root, const std::string &topoName, const std::string &fieldName, int nzones)
{
   std::string path(field_path(fieldName));
   if (!root.has_path(path))
   {
      conduit::Node &f = root[path];
      f["topology"] = topoName;
      f["association"] = "element";
      f["values"].set(conduit::DataType::float64(nzones));
      memset(f["values"].data_ptr(), 0, sizeof(conduit::float64) * nzones);
   }
}

void Teton::getRadiationTemperature(double *RadTemp, int nzones) const
{
   for (int zone = 0; zone < nzones; ++zone)
   {
      // Get the radiation temperature
      int zoneID = zone + 1;
      RadTemp[zone] = 0.;
      teton_getradiationtemperature(&zoneID, &RadTemp[zone]);
   }
}

void Teton::getMaterialTemperature(double *MatTemp, int nzones) const
{
   for (int zone = 0; zone < nzones; ++zone)
   {
      // Get the material temperature
      int zoneID = zone + 1;
      MatTemp[zone] = 0.;
      teton_getmaterialtemperature(&zoneID, &MatTemp[zone]);
   }
}

void Teton::reconstructPsi(double *rad_energy, const double *rad_energy_density)
{
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   CALI_CXX_MARK_FUNCTION;

   // Determine nzones, ngroups
   conduit::Node &blueprint = getMeshBlueprint();
   conduit::Node &options = getOptions();

   const conduit::Node &main_topo = getMainTopology(blueprint);
   conduit::index_t nzones = conduit::blueprint::mesh::utils::topology::length(main_topo);
   conduit::index_t ngroups = options.fetch_existing("quadrature/num_groups").to_index_t();

   if (doPartitioning())
   {
      conduit::Node &part = getMeshBlueprintPart();
      std::string mainTopologyName(main_topo.name());

      // Add rad_energy_density as an mcarray on the blueprint mesh.
      std::string fieldName(MCARRAY_PREFIX + "rad_energy_density");
      conduit::Node &fields = blueprint["fields"];
      conduit::Node &n_f = fields[fieldName];
      n_f["topology"] = mainTopologyName;
      n_f["association"] = "element";
      conduit::Node &values = n_f["values"];
      // radiation_energy_density is shaped double[ngroups][nzones] so we do not need to interleave.
      bool interleave = false;
      utilities::NDAccessor acc(values, {{"zone", nzones}, {"group", ngroups}}, interleave);
      acc.set_external(rad_energy_density);

      // Partition to get rad_energy_density on part mesh.
      const std::vector<std::string> fieldNames{fieldName};
      sendFieldsOrig2Part(mainTopologyName, fieldNames, false);

      // On part mesh, get various mcarray components, put back together into a
      // contiguous block [ngroups][nzones] so we can call Teton.
      conduit::Node &partFields = part["fields"];
      conduit::Node &partValues = partFields.fetch_existing(fieldName + "/values");
      const conduit::Node &part_main_topo = getMainTopology(part);
      conduit::index_t npartZones = conduit::blueprint::mesh::utils::topology::length(part_main_topo);
      utilities::NDAccessor accp(partValues, {{"zone", npartZones}, {"group", ngroups}}, interleave);
      std::vector<double> part_rad_energy_density(ngroups * npartZones);
      accp.to_contiguous(&part_rad_energy_density[0]);

      // Call Teton
      teton_reconstructpsi(rad_energy, &part_rad_energy_density[0]);

      // Clean up
      fields.remove_child(fieldName);
      partFields.remove_child(fieldName);
   }
   else
   {
      teton_reconstructpsi(rad_energy, const_cast<double *>(rad_energy_density));
   }
#endif
}

void Teton::reconstructPsiFromdV()
{
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   teton_reconstructpsifromdv();
#endif
}

void Teton::getZonalPsi(int numAngles, double *psi)
{
#if defined(TETON_PARTITIONING) && !defined(TETON_ENABLE_MINIAPP_BUILD)
   CALI_CXX_MARK_FUNCTION;

   if (doPartitioning())
   {
      conduit::Node &blueprint = getMeshBlueprint();
      conduit::Node &options = getOptions();
      conduit::Node &part = getMeshBlueprintPart();

      const conduit::Node &main_topo = getMainTopology(part);
      std::string mainTopologyName(main_topo.name());
      conduit::index_t nzones = conduit::blueprint::mesh::utils::topology::length(main_topo);
      conduit::index_t ngroups = options.fetch_existing("quadrature/num_groups").to_index_t();
      //conduit::index_t nangles = options.fetch_existing("quadrature/nSetsMaster").to_index_t();
      conduit::index_t nangles = numAngles;

      // Get the zonalpsi for the part mesh.
      double *part_psi = new double[nangles * ngroups * nzones];
      teton_getzonalpsi(&numAngles, part_psi);

      // Add zonalpsi as an mcarray on the part mesh (external).
      const std::string fieldName(MCARRAY_PREFIX + "zonalpsi");
      const std::vector<std::string> fieldNames{fieldName};
      conduit::Node &partfields = part.fetch_existing("fields");
      conduit::Node &pzpsi = partfields[fieldName];
      pzpsi["association"] = "element";
      pzpsi["topology"] = mainTopologyName;
      bool interleave = false; // since zones is the fastest changing in how the data are stored.
      utilities::NDAccessor accp(pzpsi["values"],
                                 {{"zone", nzones}, {"group", ngroups}, {"angle", nangles}},
                                 interleave);
      accp.set_external(part_psi);

      // Send zonalpsi mcarray back to the blueprint mesh.
      sendFieldsPart2Orig(mainTopologyName, fieldNames);

      // Cleanup part fields.
      partfields.remove_child(fieldName);
      delete[] part_psi;

      // Construct an accessor on the blueprint mesh's zonalpsi field and then iterate
      // over it, copying the data into the psi array.
      const conduit::Node &bp_main_topo = getMainTopology(blueprint);
      conduit::index_t nbpzones = conduit::blueprint::mesh::utils::topology::length(bp_main_topo);
      conduit::Node &fields = blueprint.fetch_existing("fields");
      conduit::Node &zpsi = fields[fieldName];
      utilities::NDAccessor acc(zpsi["values"],
                                {{"zone", nbpzones}, {"group", ngroups}, {"angle", nangles}},
                                interleave);
#if 1
      acc.to_contiguous(psi);
#else
      // For reference
      double *dptr = psi;
      for (conduit::index_t a = 0; a < nangles; a++)
         for (conduit::index_t g = 0; g < ngroups; g++)
            for (conduit::index_t z = 0; z < nzones; z++)
            {
               *dptr++ = acc(std::vector<conduit::index_t>{z, g, a});
            }
#endif
      // Cleanup blueprint fields.
      fields.remove_child(fieldName);
   }
   else
   {
      teton_getzonalpsi(&numAngles, psi);
   }
#else
   teton_getzonalpsi(&numAngles, psi);
#endif
}

//---------------------------------------------------------------------------

bool Teton::doPartitioning() const
{
   bool p = false;
#if defined(TETON_PARTITIONING)
   if (getOptions().has_path("partitioning"))
   {
      int value = getOptions().fetch_existing("partitioning").to_int();
      p = value != 0;
   }
   if (getenv("TETON_PARTITION") != nullptr)
   {
      int value = atoi(getenv("TETON_PARTITION"));
      p = value != 0;
   }
#endif
   return p;
}

conduit::Node &Teton::getMeshBlueprintPart()
{
#if defined(TETON_PARTITIONING)
   std::string meshKey(doPartitioning() ? "blueprint_partitioned" : "blueprint");
   return getDatastore()[meshKey];
#else
   return getDatastore()["blueprint"];
#endif
}

const conduit::Node &Teton::getMeshBlueprintPart() const
{
#if defined(TETON_PARTITIONING)
   std::string meshKey(doPartitioning() ? "blueprint_partitioned" : "blueprint");
   return getDatastore()[meshKey];
#else
   return getDatastore()["blueprint"];
#endif
}

std::vector<std::string> Teton::createPartitionFields(conduit::Node &mesh, const std::vector<std::string> &topoNames)
{
   std::vector<std::string> fieldNames;
#if defined(TETON_PARTITIONING)
   CALI_CXX_MARK_FUNCTION;

   // Make a list of field names
   for (const auto &tname : topoNames)
      fieldNames.push_back(Teton::PREFIX + "parmetis_result_" + tname);

   auto doms = conduit::blueprint::mesh::domains(mesh);
   for (const auto &domptr : doms)
   {
      conduit::Node &dom = *domptr;

      // Look in the domain to see if any of the fields need to be created.
      bool buildFields = true;
      if (dom.has_path("fields"))
      {
         const conduit::Node &fields = dom.fetch_existing("fields");
         int missingCount = 0;
         for (const auto &fname : fieldNames)
            missingCount = fields.has_path(fname) ? 0 : 1;
         buildFields = missingCount > 0;
      }

      // If any fields need to be created, do it.
      if (buildFields)
      {
         const conduit::Node &topo = getMainTopology(dom);
         const conduit::Node &coordset = dom.fetch_existing("coordsets/" + topo["coordset"].as_string());
         // Assume all of the input topologies will have the same topological dimension and
         // that they will match the "boundary" topology. It is possible the boundary topology
         // does not exist. Let's assume that if the boundary topology does not exist then
         // neither will the others. Any other topologies we're dealing with will have the
         // same topological dimension as the boundary topology.
         if (dom.has_path("topologies/" + TOPO_BOUNDARY))
         {
            const conduit::Node &btopo = dom.fetch_existing("topologies/" + TOPO_BOUNDARY);
            auto d = static_cast<size_t>(conduit::blueprint::mesh::utils::topology::dims(btopo));
#if defined(PARTITION_DEBUG)
            if (mRank == 0)
               std::cout << "Teton: partition - create partition field - build hash" << std::endl;
#endif
            // Produce the external "faces" of the domain and for each "face", hash
            // its node ids and associate that hash with the parent zone for the face.
            std::vector<std::pair<size_t, size_t>> desired_maps{{d, d + 1}};
            conduit::blueprint::mesh::utils::TopologyMetadata md(topo, coordset, d, desired_maps);
            const conduit::Node &dtopo = md.get_topology(d);
            auto nent = md.get_topology_length(d);
            std::map<conduit::uint64, int> hashToZone;
            for (conduit::index_t ei = 0; ei < nent; ei++)
            {
               auto vv = md.get_global_association(ei, d, d + 1);
               if (vv.size() == 1)
               {
                  // Get the ids that make up the entity and hash them.
                  auto ids = conduit::blueprint::mesh::utils::topology::unstructured::points(dtopo, ei);
                  std::sort(ids.begin(), ids.end());
                  conduit::uint64 h = conduit::utils::hash(&ids[0], static_cast<unsigned int>(ids.size()));

                  // Save hash to parent zone.
                  hashToZone[h] = vv[0];
               }
            }

            // Get the partition field for the main topology.
            const conduit::Node &pf = dom.fetch_existing("fields/" + PARTITION_FIELD + "/values");
            auto f = pf.as_int32_accessor();

            // Now, iterate through the secondary topologies, hash each entity's ids
            // and try to look up the parent zone. The hashToZone map should contain
            // all possible external faces for the domain so the boundary should be
            // a subset of that.
            for (size_t ti = 0; ti < topoNames.size(); ti++)
            {
               std::string fieldKey("fields/" + fieldNames[ti]);
               if (!dom.has_path(fieldKey))
               {
                  // Only make the field if the topology exists on this rank.
                  std::string topoKey("topologies/" + topoNames[ti]);
                  if (dom.has_path(topoKey))
                  {
                     const conduit::Node &topo = dom.fetch_existing(topoKey);
                     auto blen = conduit::blueprint::mesh::topology::length(topo);
#if defined(PARTITION_DEBUG)
                     if (mRank == 0)
                        std::cout << "Teton: partition - create partition field - " << fieldNames[ti] << std::endl;
#endif
                     conduit::Node &newfield = dom[fieldKey];
                     newfield["association"] = "element";
                     newfield["topology"] = topoNames[ti];
                     newfield["values"].set(conduit::DataType::int32(blen));
                     auto topoPartition = newfield["values"].as_int32_ptr();
                     for (conduit::index_t ei = 0; ei < blen; ei++)
                     {
                        // Get the ids that make up the entity and hash them.
                        auto ids = conduit::blueprint::mesh::utils::topology::unstructured::points(topo, ei);
                        std::sort(ids.begin(), ids.end());
                        conduit::uint64 h = conduit::utils::hash(&ids[0], static_cast<unsigned int>(ids.size()));

                        // Look up the zone id and map it through the partition field.
                        const auto it = hashToZone.find(h);
                        topoPartition[ei] = (it != hashToZone.end()) ? f[it->second] : 0;
                     }
                  }
               }
            }
         }
      }
   }
#endif
   return fieldNames;
}

void Teton::assimilateTopology(conduit::Node &partmesh,
                               const std::string &topoName,
                               conduit::Node &secondPartmesh,
                               const std::string &secondTopoName)
{
#if defined(TETON_PARTITIONING)
   CALI_CXX_MARK_FUNCTION;

   if (partmesh.dtype().is_empty())
   {
      return;
   }
   if (secondPartmesh.dtype().is_empty())
   {
      return;
   }

   // Get the coordset and topo for the volume mesh
   auto domains = conduit::blueprint::mesh::domains(partmesh);
   if (domains.size() < 1)
   {
      std::cout << "assimilateTopology: Must have at least one domain: " << domains.size() << std::endl;
      return;
   }

   // Get the coordset and topo for the boundary mesh
   auto bdomains = conduit::blueprint::mesh::domains(secondPartmesh);
   if (bdomains.size() < 1)
   {
      std::cout << "assimilateTopology: Must have at least one domain: " << bdomains.size() << std::endl;
      return;
   }

   if (domains.size() != bdomains.size())
   {
      std::cout << "assimilateTopology: Incompatible numbers of domains " << domains.size() << ", " << bdomains.size()
                << std::endl;
      return;
   }

   for (size_t domid = 0; domid < domains.size(); domid++)
   {
      conduit::Node &mesh = *domains[domid];
      const conduit::Node &topo = mesh["topologies/" + topoName];
      const conduit::Node &coordset = conduit::blueprint::mesh::utils::topology::coordset(topo);

      conduit::Node &bmesh = *bdomains[domid];
      const conduit::Node &btopo = bmesh["topologies/" + secondTopoName];
      const conduit::Node &bcoordset = conduit::blueprint::mesh::utils::topology::coordset(btopo);

      // Iterate over the boundary mesh coordinates and look them up in the
      // volume mesh's coordset.
      conduit::blueprint::mesh::utils::query::PointQuery Q(mesh);
      const auto axes = conduit::blueprint::mesh::utils::coordset::axes(bcoordset);
      const auto ndims = axes.size();
      const conduit::Node &bcvalues = bcoordset.fetch_existing("values");
      const int domain_id = static_cast<int>(domid);
      const auto bx = bcvalues[axes[0]].as_double_accessor();
      const auto by = bcvalues[axes[ndims > 1 ? 1 : 0]].as_double_accessor();
      const auto bz = bcvalues[axes[ndims > 2 ? 2 : 0]].as_double_accessor();
      conduit::index_t nSearchPoints = bx.number_of_elements();
      for (conduit::index_t i = 0; i < nSearchPoints; i++)
      {
         double pt[3];
         pt[0] = bx[i];
         pt[1] = ndims > 1 ? by[i] : 0.;
         pt[2] = ndims > 2 ? bz[i] : 0.;
         Q.add(domain_id, pt);
      }
      Q.execute(coordset.name());

      // Make a new the topology that uses the volume mesh coordset.
      // We remap the connectivity.
      const auto &res = Q.results(domain_id);
      auto bconnSrc = btopo["elements/connectivity"].as_int32_accessor();
      conduit::index_t nbconn = bconnSrc.number_of_elements();
      conduit::Node &newtopo = mesh["topologies/" + secondTopoName];
      newtopo["type"] = btopo["type"];
      newtopo["coordset"] = coordset.name();
      newtopo["elements/shape"] = btopo["elements/shape"];
      newtopo["elements/connectivity"].set(conduit::DataType::int32(nbconn));
      auto bconnNew = newtopo["elements/connectivity"].as_int32_array();
      for (conduit::index_t i = 0; i < nbconn; i++)
      {
         bconnNew[i] = res[bconnSrc[i]];
      }
      if (btopo.has_path("elements/sizes"))
         btopo["elements/sizes"].to_data_type(conduit::DataType::int32().id(), newtopo["elements/sizes"]);
      if (btopo.has_path("elements/offsets"))
         btopo["elements/offsets"].to_data_type(conduit::DataType::int32().id(), newtopo["elements/offsets"]);

      // Iterate over the boundary mesh's fields and steal them for the partmesh.
      if (bmesh.has_child("fields"))
      {
         conduit::Node &srcFields = bmesh["fields"];
         conduit::Node &destFields = mesh["fields"];
         for (conduit::index_t i = 0; i < srcFields.number_of_children(); i++)
         {
            conduit::Node &f = srcFields[i];
            destFields[f.name()].set(f);
         }
      }
   }
#endif
}

void Teton::createRadiationTemperature()
{
   // Checking pre-partition so we use blueprint instead of part.
   conduit::Node &blueprint = getMeshBlueprint();

   // This field is used to return radiation temperature values from Teton as a
   // field on the mesh that will be queried via getRadiationTemperature().
   if (doPartitioning())
   {
      const conduit::Node &main_topo = getMainTopology(blueprint);
      std::string mainTopologyName(main_topo.name());
      const conduit::index_t nzones = conduit::blueprint::mesh::utils::topology::length(main_topo);
      createZonalField(blueprint, mainTopologyName, FIELD_RADIATION_TEMPERATURE, nzones);
   }
}

void Teton::createMaterialTemperature()
{
   // Checking pre-partition so we use blueprint instead of part.
   conduit::Node &blueprint = getMeshBlueprint();

   // This field is used to return radiation temperature values from Teton as a
   // field on the mesh that will be queried via getMaterialTemperature().
   if (doPartitioning())
   {
      const conduit::Node &main_topo = getMainTopology(blueprint);
      std::string mainTopologyName(main_topo.name());
      const conduit::index_t nzones = conduit::blueprint::mesh::utils::topology::length(main_topo);
      createZonalField(blueprint, mainTopologyName, FIELD_MATERIAL_TEMPERATURE, nzones);
   }
}

bool Teton::doInterleave(const std::string &fieldName) const
{
   bool retval = fieldName != FIELD_RADIATION_ENERGY_DENSITY;
   return retval;
}

void Teton::add_mcarray_fields(conduit::Node &root)
{
#if defined(TETON_PARTITIONING)
   conduit::Node &options = getOptions();
   const conduit::Node &main_topo = getMainTopology(root);
   conduit::Node &fields = root.fetch_existing("fields");

   conduit::index_t nzones = conduit::blueprint::mesh::utils::topology::length(main_topo);
   conduit::index_t ngroups = options.fetch_existing("quadrature/num_groups").to_index_t();

   // Some fields have been provided as a single oversized buffer that has
   // multiple components but is not actually an MCArray. We'll make an MCArray
   // for that field, giving it a new name. This is done so we can have
   // the field pass through the partitioner while preserving all of the data.

   // Make new mcarray fields for the fields that look like mcarrays.
   const std::vector<std::string> copy_keys{"association", "topology"};
   utilities::iterate_mcarray_candidates(root,
                                         main_topo.name(),
                                         options,
                                         std::vector<std::string>{},
                                         [&](const conduit::Node &srcField)
                                         {
      std::string fieldName(srcField.name());
      std::string newFieldName(MCARRAY_PREFIX + srcField.name());

      // Copy basic attributes
      conduit::Node &newField = fields[newFieldName];
      for (const auto &k : copy_keys)
      {
         if (srcField.has_child(k))
            newField[k].set(srcField.fetch_existing(k));
      }

      // Make the mcarray so it points to the original field's data. These fields are interleaved.
      const conduit::Node &srcValues = srcField.fetch_existing("values");
      const double *srcData = srcValues.as_float64_ptr();
      utilities::NDAccessor acc(newField["values"], {{"zone", nzones}, {"group", ngroups}}, doInterleave(fieldName));
      acc.set_external(srcData);

      // Record that we made a new mcarray.
      mMCArrays[fieldName] = newFieldName;
   });
#endif
}

void Teton::remove_mcarray_fields(conduit::Node &root)
{
#if defined(TETON_PARTITIONING)
   if (doPartitioning())
   {
      conduit::Node &fields = root["fields"];
      for (auto it = mMCArrays.begin(); it != mMCArrays.end(); it++)
      {
         if (fields.has_path(it->second))
         {
            fields.remove(it->second);
         }
      }
   }
#endif
}

const conduit::Node &Teton::fetch_mcarray(const conduit::Node &root, const std::string &fieldName) const
{
   const conduit::Node &fields = root.fetch_existing("fields");

   // See if the requested field is an mcarray we created.
   auto it = mMCArrays.find(fieldName);
   if (it != mMCArrays.end())
   {
      // Return the mcarray for the requested old field name.
      return fields.fetch_existing(it->second);
   }
   // We did not find it. Assume that fieldName is a regular field.
   return fields.fetch_existing(fieldName);
}

void Teton::partition(bool fromRestart)
{
#if defined(TETON_PARTITIONING)
   // Make a lambda that helps display some dtype information for the mesh.
#if defined(PARTITION_DEBUG)
   auto check_widest_dtype = [](const conduit::Node &mesh, int rank, const std::string &caption)
   {
      auto dtype = conduit::blueprint::mesh::utils::find_widest_dtype(mesh, conduit::DataType::int32());
      std::cout << "  Widest " << caption << " int type: " << dtype.name() << std::endl;
      if (!dtype.is_int32())
      {
         const auto paths = utilities::find_int64(mesh);
         if (!paths.empty())
         {
            if (rank == 0)
            {
               std::cout << "int64 paths: ";
               for (const auto &p : paths)
               {
                  std::cout << p << ", ";
               }
               std::cout << std::endl;
            }
         }
      }
   };

   utilities::Banner b(mCommunicator, "Teton::partition");
#endif
   conduit::Node &blueprint = getMeshBlueprint();
   conduit::Node &part = getMeshBlueprintPart();
   bool alreadyPartitioned = part.has_child("partition_options_main");
   if (doPartitioning() && (!alreadyPartitioned || fromRestart))
   {
      CALI_CXX_MARK_SCOPE("Teton_Partition");

      int rank = 0, size = 1;
      MPI_Comm_rank(mCommunicator, &rank);
      MPI_Comm_size(mCommunicator, &size);

      std::string mainTopoName(getMainTopology(blueprint).name());

      // Create a new field on the blueprint mesh that we're partitioning.
      conduit::Node opts;
      opts["topology"] = mainTopoName;
      opts["field_prefix"] = PREFIX;
      if (blueprint.has_path("adjsets/main_adjset"))
         opts["adjset"] = "main_adjset"; // plays a role in global node id generation
      if (rank == 0)
         std::cout << "Teton: partition - make partition field." << std::endl;
      conduit::blueprint::mpi::mesh::generate_partition_field(blueprint, opts, mCommunicator);

      // The partition field will be index_t, which can make Conduit start
      // generating index_t for other things like topology maps. Force int32
      // because other int types are catastrophic in the interface at this time.
      const std::vector<std::string> replacements{"fields/" + PARTITION_FIELD + "/values",
                                                  "fields/" + PREFIX + "global_element_ids/values",
                                                  "fields/" + PREFIX + "global_vertex_ids/values"};
      utilities::convert_int32(rank, blueprint, replacements);

      // There are no int64/index_t in blueprint mesh now.
#if defined(PARTITION_DEBUG)
      auto dtype = conduit::blueprint::mesh::utils::find_widest_dtype(blueprint, conduit::DataType::int32());
      if (rank == 0)
         std::cout << "  Widest blueprint int type: " << dtype.name() << std::endl;
      MPI_Barrier(mCommunicator);
#endif
      // Create a field selection for the main topology.
      conduit::Node partopts;
      partopts["mapping"] = 1;
      partopts["original_element_ids"] = "main_original_element_ids";
      partopts["original_vertex_ids"] = "main_vertex_element_ids";
      conduit::Node &selections = partopts["selections"];
      conduit::Node &sel1 = partopts["selections"].append();
      sel1["type"] = "field";
      sel1["domain_id"] = "any";
      sel1["field"] = PARTITION_FIELD;
      sel1["topology"] = mainTopoName;
      sel1["destination_ranks"].set(conduit::DataType::int32(size));
      auto ranks = sel1["destination_ranks"].as_int32_ptr();
      for (int i = 0; i < size; i++)
         ranks[i] = i;

#if defined(PARTITION_DEBUG)
      // Print the partitioning options.
      if (rank == 0)
      {
         std::cout << "Teton: partition enabled." << std::endl;
         std::cout << "Teton: part.path=" << part.path() << std::endl;
         std::cout << "Teton: partops:" << std::endl;
         partopts.print();
      }

      // Save out the partition mesh and parameters to YAML
      const std::string protocol("yaml");
      std::stringstream ss, ss2;
      ss << "partition_mesh." << rank << "." << protocol;
      std::string meshFilename(ss.str());
      ss2 << "partition_options_main." << rank << "." << protocol;
      std::string optsFilename(ss2.str());
#pragma message "Partition mesh input will be saved."
      conduit::relay::io::save(blueprint, meshFilename, protocol);
      conduit::relay::io::save(partopts, optsFilename, protocol);
#endif

      // Partition the blueprint mesh and store the results in part.
      if (rank == 0)
         std::cout << "Teton: partition - partition main" << std::endl;
      add_mcarray_fields(blueprint);
      conduit::blueprint::mpi::mesh::partition(blueprint, partopts, part, mCommunicator);
      remove_mcarray_fields(blueprint);

      // If we sent fields to the part mesh that were mcarray, keep only the
      // mcarray version of the field. The non-mcarray version is probably not
      // the right size.
      for (auto it = mMCArrays.begin(); it != mMCArrays.end(); it++)
      {
         std::string origFieldKey = "fields/" + it->first;
         if (part.has_path(origFieldKey))
         {
            part.remove(origFieldKey);
         }
      }

      // The partitioner may produce index_t data, even if the inputs were not
      // index_t. We must convert it to int32 or other parts of Conduit may
      // start inserting index_t data because they find the widest dtype. This
      // causes problems for Teton down the line since it requires int32.
      auto repkeys = utilities::find_int64(part);
#if defined(PARTITION_DEBUG)
      if (rank == 0)
         check_widest_dtype(part, rank, "part");
      MPI_Barrier(mCommunicator);
#endif

      // Do the conversion
      utilities::convert_int32(rank, part, repkeys);

#if defined(PARTITION_DEBUG)
      // Double-check that the dtype is int32.
      if (rank == 0)
         check_widest_dtype(part, rank, "part");
      MPI_Barrier(mCommunicator);
#endif
#if defined(PARTITION_DEBUG) && defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
      // Save the partitioned mesh to a file that can be visualized.
      const std::string file_protocol = "hdf5";
      if (rank == 0)
         std::cout << "Teton: partition - save part node to " << file_protocol << std::endl;
      conduit::relay::mpi::io::blueprint::save_mesh(part, "part", file_protocol, mCommunicator);
#endif
      // Get the names of the topos that we need to partition in addition to main,
      // the "secondary" topologies. These include the boundary topology.
      auto partitionTopos = getPartitionTopologies(blueprint);

      // Create partition field for secondary topologies.
      if (rank == 0)
         std::cout << "Teton: partition - create partition fields" << std::endl;
      auto partitionFields = createPartitionFields(blueprint, partitionTopos);

      // Process each secondary topology.
      for (size_t ti = 0; ti < partitionTopos.size(); ti++)
      {
         // Now repartition the secondary topology.
         conduit::Node tpartopts;
         conduit::Node &tsel = tpartopts["selections"].append();
         tsel["type"] = "field";
         tsel["domain_id"] = "any";
         tsel["field"] = partitionFields[ti];
         tsel["topology"] = partitionTopos[ti];
         tsel["destination_ranks"].set(sel1["destination_ranks"]); // same rank map
         tpartopts["mapping"] = 0;                                 // If we ever need to map fields back, set this to 1.
         tpartopts["original_element_ids"] = partitionTopos[ti] + "_original_element_ids";
         tpartopts["original_vertex_ids"] = partitionTopos[ti] + "_vertex_element_ids";
         tpartopts["build_adjsets"] = 0;
         if (partitionTopos[ti] != TOPO_BOUNDARY)
         {
            // Restrict fields mapped to one that does not exist. Map no fields.
            // This is done because surface flux topologies fields are not really
            // fields on that mesh.
            conduit::Node &fields = tpartopts["fields"];
            fields[PREFIX + "_impossible_to_find_123456789"] = 1;
         }

         if (rank == 0)
         {
            std::cout << "Teton: partition - partition " << partitionTopos[ti] << std::endl;
#if defined(PARTITION_DEBUG)
            tpartopts.print();
#endif
         }
         conduit::Node newpart;
         conduit::blueprint::mpi::mesh::partition(blueprint, tpartopts, newpart, mCommunicator);

         // Merge the secondary topology from newpart into part.
         assimilateTopology(part, mainTopoName, newpart, partitionTopos[ti]);

         // Save partitioning options for later.
         part["partition_options_" + partitionTopos[ti]].move(tpartopts);
      }

      // Do the conversion
      repkeys = utilities::find_int64(part);
      utilities::convert_int32(rank, part, repkeys);

      // Save partitioning options for later. We have to do this after calls
      // to mesh::partition() because that method resets the input node.
      part["partition_options_main"].move(partopts);

#if defined(PARTITION_DEBUG)
      check_widest_dtype(part, rank, "combined part");
      MPI_Barrier(mCommunicator);

#if defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
      // Save the partitioned mesh, plus boundary to a file that can be visualized.
      if (rank == 0)
         std::cout << "Teton: partition - save part_with_boundary to " << file_protocol << std::endl;
      conduit::relay::mpi::io::blueprint::save_mesh(part, "part_with_boundary", file_protocol, mCommunicator);
#endif

      // Check whether there are duplicated local points. We hope not.
      conduit::Node info;
      bool dups = utilities::find_local_duplicate_points(rank, part, part["coordsets/coords"], info);
      if (dups)
      {
         info.print();
      }
      MPI_Barrier(mCommunicator);
#endif
   }
#endif
}

std::vector<std::string> Teton::getPartitionTopologies(const conduit::Node &root) const
{
   std::vector<std::string> topoNames;

   // We'll exclude any of these names. The "main" topology is handled explicitly
   // and we do not want any of the derived topologies to go through partitioning.
   const std::string mainTopologyName(getMainTopology(root).name());
   const std::vector<std::string> exclusions{mainTopologyName, "main_corner", "main_face"};

   const conduit::Node &topologies = root.fetch_existing("topologies");
   const std::string mainCoordset(getMainTopology(root).fetch_existing("coordset").as_string());
   // Make a vector of all of the topology names that share a coordset with main.
   for (conduit::index_t i = 0; i < topologies.number_of_children(); i++)
   {
      const conduit::Node &topo = topologies[i];

      if (std::find(exclusions.begin(), exclusions.end(), topo.name()) == exclusions.end())
      {
         const std::string coordset = topo.fetch_existing("coordset").as_string();
         if (coordset == mainCoordset)
         {
            topoNames.push_back(topo.name());
         }
      }
   }
   return utilities::globalizeStringVector(topoNames, mCommunicator);
}

void Teton::sendFieldsOrig2Part(const std::string &topoName,
                                const std::vector<std::string> &fieldNames,
                                bool updateCoords)
{
#if defined(TETON_PARTITIONING)
   CALI_CXX_MARK_FUNCTION;

   if (doPartitioning() && (!fieldNames.empty() || updateCoords))
   {
#if defined(PARTITION_DEBUG)
      utilities::Banner b(mCommunicator, "sendFieldsOrig2Part");
      if (mRank == 0)
      {
         std::cout << "topoName: " << topoName << "\n";
         std::cout << "updateCoords: " << updateCoords << "\n";
         std::cout << "fieldNames:\n";
         for (const auto &name : fieldNames)
            std::cout << "  - \"" << name << "\"" << std::endl;
      }
      MPI_Barrier(mCommunicator);
#endif
      conduit::Node &blueprint = getMeshBlueprint();
      conduit::Node &part = getMeshBlueprintPart();
      const conduit::Node &options = getOptions();
      conduit::Node &partopts = part["partition_options_" + topoName];

      // Make sure mcarray fields are up to date.
      add_mcarray_fields(blueprint);

      conduit::Node newpartmesh, updateopts;
      // Copy the partition options for the topology and restrict the fields
      // (add fields only if they exist in the blueprint mesh).
      updateopts.set(partopts);
      updateopts["mapping"] = 0;
      updateopts["build_adjsets"] = 0;
      conduit::Node &fields = updateopts["fields"];
      for (const auto &f : fieldNames)
      {
         // If the field happens to be the name of an mcarray then we want to send
         // the mcarray instead of the original field name.
         const conduit::Node &mcf = fetch_mcarray(blueprint, f);
         std::string actualField(mcf.name());
         if (blueprint.has_path("fields/" + actualField))
            fields[actualField] = 1;
      }
      // Without a "fields" node, the partitioner will attempt to map all fields to
      // the partitioned mesh. We really do not want to do any fields here since the
      // host code could have deallocated the field memory that we know about in the
      // blueprint node. This actually happened! Since we need a fields node, make up
      // a field name that would never exist so the partitioner will not find any of
      // the field names in this list.
      if (fieldNames.empty())
         fields[PREFIX + "_impossible_to_find_123456789"] = 1;

      // Partition the mesh again.
      conduit::blueprint::mpi::mesh::partition(blueprint, updateopts, newpartmesh, mCommunicator);

      // Iterate through the fields in each domain and move them over to the part mesh.
      // This assumes that the newpartmesh was partitioned the same way as the
      // original part mesh, which should be true.
      auto destDoms = conduit::blueprint::mesh::domains(part);
      auto srcDoms = conduit::blueprint::mesh::domains(newpartmesh);
      assert(destDoms.size() == srcDoms.size());
      for (size_t i = 0; i < srcDoms.size(); i++)
      {
         conduit::Node &srcFields = srcDoms[i]->fetch_existing("fields");
         conduit::Node &destFields = destDoms[i]->fetch_existing("fields");
         for (conduit::index_t fi = 0; fi < srcFields.number_of_children(); fi++)
         {
            conduit::Node &src = srcFields[fi];
            std::string fname(src.name());
            // Move the field over if it was one we wanted in the options.
            if (fields.has_child(fname))
            {
               // BJW: There is a chance that the partitioned domain was passed straight
               //      through the partitioner without modification. For those cases, we
               //      do not want to disturb the destFields/srcFields.
               conduit::Node &dest = destFields[fname];

               // If the src and dest field are not the same then we need to copy.
               bool doCopy = true;
               if (dest.has_path("values"))
               {
                  conduit::Node &srcValues = src["values"];
                  conduit::Node &destValues = dest["values"];
                  if (srcValues.number_of_children() > 0 && destValues.number_of_children() > 0)
                     doCopy = destValues[0].data_ptr() != srcValues[0].data_ptr();
                  else
                     doCopy = destValues.data_ptr() != srcValues.data_ptr();
               }
               if (doCopy)
               {
                  // Copy the field to destFields.
                  dest.reset();
                  dest["association"] = src["association"];
                  dest["topology"] = src["topology"];
                  dest["values"] = src["values"];
               }
            }
         }

         // Move coordset from src to dest.
         if (updateCoords)
         {
            conduit::Node &srcCoords = srcDoms[i]->fetch_existing("coordsets/coords");
            conduit::Node &destCoords = destDoms[i]->fetch_existing("coordsets/coords");
            destCoords.set(srcCoords);
         }
      }

      // Remove mcarray fields
      remove_mcarray_fields(blueprint);
   }
#endif
}

void Teton::sendFieldsPart2Orig(const std::string &topoName, const std::vector<std::string> &fieldNames)
{
#if defined(TETON_PARTITIONING)
   CALI_CXX_MARK_FUNCTION;

   // The plan here is to send fields from the part mesh back to the original mesh.
   // For normal fields, we can do this no problem. For oversize "mcarray" fields,
   // we send back the mcarray field names and then recopy their data into the
   // field they represent. On the partmesh, the mcarray components are often separate
   // fields since they likely came from the original mesh in the first place and
   // were not reaggregated into contiguous memory.

   if (doPartitioning())
   {
      conduit::Node &blueprint = getMeshBlueprint();
      conduit::Node &part = getMeshBlueprintPart();
      const conduit::Node &options = getOptions();
      const std::string mainTopologyName(getMainTopology(blueprint).name());

#if defined(PARTITION_DEBUG)
      utilities::Banner b(mCommunicator, "sendFieldsPart2Orig");
      if (mRank == 0)
      {
         std::cout << "topoName: " << topoName << std::endl;
         std::cout << "fieldNames:\n";
         for (const auto &name : fieldNames)
            std::cout << "  - \"" << name << "\"" << std::endl;
      }
      MPI_Barrier(mCommunicator);
#endif

      // Make sure the part mesh has its mcarray fields wrapped to send back.
      add_mcarray_fields(part);

      // Get the list of fields that we think need to be mcarrays. If we're sending
      // back one of these, send back the mcarray instead since the mcarrays are often
      // set_external'd from the host. Also, much of the time, we send the mcarray
      // data
      std::map<std::string, std::string> normal2mcarray;
      utilities::iterate_mcarray_candidates(blueprint,
                                            mainTopologyName,
                                            options,
                                            fieldNames,
                                            [&](const conduit::Node &f)
                                            { normal2mcarray[f.name()] = MCARRAY_PREFIX + f.name(); });

      // Build up the mapback options.
      conduit::Node mbopts;
      mbopts["field_prefix"] = PREFIX;
      mbopts["original_element_ids"] = topoName + "_original_element_ids";
      mbopts["original_vertex_ids"] = topoName + "_original_vertex_ids";
      std::vector<std::string> sname;
      for (const auto &f : fieldNames)
      {
         // If the field we want to send back seems like an mcarray, we will send back
         // the mcarray instead.
         std::string sendName(f);
         const auto it = normal2mcarray.find(f);
         if (it != normal2mcarray.end())
            sendName = it->second;

         mbopts["fields"].append().set(sendName);

         // Make sure the part mesh contains the field name.
         if (part.has_path("fields/" + sendName))
         {
            sname.push_back(sendName);
         }
      }

      if (mbopts.has_path("fields") && mbopts["fields"].number_of_children() > 0)
      {
         // Move selected fields from the partitioned mesh back to the original mesh.
         conduit::blueprint::mpi::mesh::partition_map_back(part, mbopts, blueprint, mCommunicator);

         if (topoName == mainTopologyName)
         {
            // At this point, fields have been mapped from part fields back onto blueprint
            // fields. Copy the mcarray fields back to their single-buffer original variable.
            const conduit::Node &main_topo = getMainTopology(blueprint);
            conduit::Node &fields = blueprint.fetch_existing("fields");
            conduit::index_t nzones = conduit::blueprint::mesh::utils::topology::length(main_topo);
            conduit::index_t ngroups = options.fetch_existing("quadrature/num_groups").to_index_t();

            for (auto it = normal2mcarray.begin(); it != normal2mcarray.end(); it++)
            {
               conduit::Node &n_origField = fields[it->first];
               conduit::Node &n_origValues = n_origField["values"];
               double *origValues = n_origValues.as_float64_ptr();

               // Get the mcarray field we're copying into the original field.
               const conduit::Node &n_srcField = fields[it->second];
               auto values = const_cast<conduit::Node &>(n_srcField.fetch_existing("values"));

               // Copy mcarray components back into original contiguous field. This is
               // compatible with original fields that are set_external.
               utilities::NDAccessor src(values, {{"zone", nzones}, {"group", ngroups}}, doInterleave(it->first));
               src.to_contiguous(origValues);

               // We don't need the mcarray data anymore.
               fields.remove(it->second);
            }
         }
      }
   }
#endif
}

void Teton::initializeRadiationFluxFieldNames()
{
   // Get the number of dimensions from the blueprint mesh since it might not be
   // in the options yet.
   const conduit::Node &blueprint = getMeshBlueprint();
   std::string csname(getMainTopology(blueprint).fetch_existing("coordset").as_string());
   const conduit::Node &coordset = blueprint.fetch_existing("coordsets/" + csname);
   const int ndim = static_cast<int>(conduit::blueprint::mesh::coordset::dims(coordset));

   mRadiationFluxFields.clear();
   mRadiationFluxFields.reserve(ndim);
   if (ndim == 1)
   {
      mRadiationFluxFields.emplace_back(FIELD_RADIATION_FLUX_X);
   }
   else if (ndim == 2)
   {
      mRadiationFluxFields.emplace_back(FIELD_RADIATION_FLUX_Z);
      mRadiationFluxFields.emplace_back(FIELD_RADIATION_FLUX_R);
   }
   else if (ndim == 3)
   {
      mRadiationFluxFields.emplace_back(FIELD_RADIATION_FLUX_X);
      mRadiationFluxFields.emplace_back(FIELD_RADIATION_FLUX_Y);
      mRadiationFluxFields.emplace_back(FIELD_RADIATION_FLUX_Z);
   }
}

const std::vector<std::string> &Teton::getRadiationFluxFields() const
{
   return mRadiationFluxFields;
}

void Teton::setRadiationFlux()
{
   CALI_CXX_MARK_FUNCTION;

   // Instruct Teton to compute the radiation flux prior to retrieval
   teton_setradiationflux();

#if defined(TETON_PARTITIONING)
   if (doPartitioning())
   {
      conduit::Node &options = getOptions();
      conduit::Node &part = getMeshBlueprintPart();

      // Step 1, make some fields on the part mesh and store radiation flux into them.
      const conduit::Node &main_topo = getMainTopology(part);
      const std::string mainTopologyName(main_topo.name());
      conduit::Node &fields = part.fetch_existing("fields");
      conduit::index_t nzones = conduit::blueprint::mesh::utils::topology::length(main_topo);
      conduit::index_t ngroups = options.fetch_existing("quadrature/num_groups").to_index_t();

      // We arrange the data this way so we can more easily send it through the partitioner.
      std::vector<conduit::float64 *> dimGroups[3];
      const auto fieldNames = getRadiationFluxFields();
      int ndims = static_cast<int>(fieldNames.size());
      for (int dim = 0; dim < ndims; dim++)
      {
         conduit::Node &n = fields[fieldNames[dim]];
         n["topology"] = mainTopologyName;
         n["association"] = "element";
         conduit::Node &values = n["values"];
         for (int g = 0; g < ngroups; g++)
         {
            std::stringstream gs;
            gs << "group" << g;
            std::string gname(gs.str());

            values[gname].set(conduit::DataType::float64(nzones));
            dimGroups[dim].push_back(values[gname].as_float64_ptr());
         }
      }

      // step 2, get data from Teton into Conduit fields
      std::vector<double> zflux(ngroups * ndims);
      for (int zone = 0; zone < nzones; zone++)
      {
         int zone1 = zone + 1;
         teton_getradiationflux(&zone1, &zflux[0]);

         int idx = 0;
         for (int g = 0; g < ngroups; g++)
         {
            for (int dim = 0; dim < ndims; dim++)
            {
               dimGroups[dim][g][zone] = zflux[idx++];
            }
         }
      }

      // step 3, map back fields to main mesh.
      sendFieldsPart2Orig(mainTopologyName, fieldNames);
   }
#endif
}

void Teton::getRadiationFlux(int zone, double *zflux) const
{
#if defined(TETON_PARTITIONING)
   if (doPartitioning())
   {
      const conduit::Node &options = getOptions();
      const conduit::Node &blueprint = getMeshBlueprint();
      const conduit::Node &fields = blueprint.fetch_existing("fields");
      const conduit::index_t ngroups = options.fetch_existing("quadrature/num_groups").to_index_t();
      int zone0 = zone - 1;

      // Pull the data out from the Conduit fields and return in the order that
      // Teton would have returned it.
      const auto fieldNames = getRadiationFluxFields();
      int ndims = static_cast<int>(fieldNames.size());
      for (int dim = 0; dim < ndims; dim++)
      {
         const conduit::Node &f = fields.fetch_existing(fieldNames[dim]);
         const conduit::Node &values = f["values"];
         for (conduit::index_t g = 0; g < values.number_of_children(); g++)
         {
            const auto zonal_array = values[g].as_float64_ptr();
            zflux[g * ndims + dim] = zonal_array[zone0];
         }
      }
   }
   else
   {
      teton_getradiationflux(&zone, zflux);
   }
#endif
}

void Teton::getEdits(int &noutrt,
                     int &ninrt,
                     int &ngdart,
                     int &nNLIters,
                     int &maxNLIters,
                     int &TrMaxZone,
                     int &TeMaxZone,
                     int &TrMaxProcess,
                     int &TeMaxProcess,
                     double &dtused,
                     double &dtrad,
                     double &TrMax,
                     double &TeMax,
                     double &EnergyRadiation,
                     double &PowerIncident,
                     double &PowerEscape,
                     double &PowerAbsorbed,
                     double &PowerEmitted,
                     double &PowerExtSources,
                     double &PowerCompton,
                     double &EnergyCheck) const
{
   teton_getedits(&noutrt,
                  &ninrt,
                  &ngdart,
                  &nNLIters,
                  &maxNLIters,
                  &TrMaxZone,
                  &TeMaxZone,
                  &TrMaxProcess,
                  &TeMaxProcess,
                  &dtused,
                  &dtrad,
                  &TrMax,
                  &TeMax,
                  &EnergyRadiation,
                  &PowerIncident,
                  &PowerEscape,
                  &PowerAbsorbed,
                  &PowerEmitted,
                  &PowerExtSources,
                  &PowerCompton,
                  &EnergyCheck);

   // If partitioning occurred then we need to fix up some max zone/process information.
   if (doPartitioning())
   {
      const conduit::Node &part = getMeshBlueprintPart();
      const std::string mainTopologyName(getMainTopology(part).name());
      std::string vkey = "fields/" + mainTopologyName + "_original_element_ids/values";
      const conduit::Node &vnode = part.fetch_existing(vkey);
      const auto orig_domains = vnode.fetch_existing("domains").as_int_accessor();
      const auto orig_zones = vnode.fetch_existing("ids").as_int_accessor();

      int maxvals[4] = {0, 0, 0, 0}, finalmaxvals[4] = {0, 0, 0, 0};
      if (mRank == TrMaxProcess)
      {
         // This rank owns the partitioned zone so it knows where it came from.
         maxvals[0] = orig_zones[TrMaxZone - 1] + 1;
         maxvals[1] = orig_domains[TrMaxZone - 1];
      }
      if (mRank == TeMaxProcess)
      {
         // This rank owns the partitioned zone so it knows where it came from.
         maxvals[2] = orig_zones[TeMaxZone - 1] + 1;
         maxvals[3] = orig_domains[TeMaxZone - 1];
      }
      MPI_Allreduce(maxvals, finalmaxvals, 4, MPI_INT, MPI_MAX, mCommunicator);

      TrMaxZone = finalmaxvals[0];
      TrMaxProcess = finalmaxvals[1];
      TeMaxZone = finalmaxvals[2];
      TeMaxProcess = finalmaxvals[3];
   }
}

void Teton::getDtControls(int &flag, int &process, int &zone, std::string &message) const
{
   const conduit::Node &options = getOptions();
   flag = options.fetch_existing("iteration/dtcontrol/flag").value();
   process = options.fetch_existing("iteration/dtcontrol/process").value();
   zone = options.fetch_existing("iteration/dtcontrol/zone").value();
   message = options.fetch_existing("iteration/dtcontrol/message").as_string();

   // If partitioning occurred then we need to fix up some max zone/process information.
   if (doPartitioning())
   {
      const conduit::Node &part = getMeshBlueprintPart();
      const std::string mainTopologyName(getMainTopology(part).name());
      // Check this in case partitioning has not actually happened yet.
      std::string vkey = "fields/" + mainTopologyName + "_original_element_ids/values";
      if (part.has_path(vkey))
      {
         const conduit::Node &vnode = part.fetch_existing(vkey);
         const auto orig_domains = vnode.fetch_existing("domains").as_int_accessor();
         const auto orig_zones = vnode.fetch_existing("ids").as_int_accessor();

         int maxvals[2] = {0, 0}, finalmaxvals[2] = {0, 0};
         if (mRank == process)
         {
            // This rank owns the partitioned zone so it knows where it came from.
            maxvals[0] = orig_zones[zone - 1] + 1;
            maxvals[1] = orig_domains[zone - 1];
         }
         MPI_Allreduce(maxvals, finalmaxvals, 2, MPI_INT, MPI_MAX, mCommunicator);

         zone = finalmaxvals[0];
         process = finalmaxvals[1];

         message += " Note - Teton repartitioned the mesh so process/zone may differ.";
      }
   }
}

void Teton::partitionCleanup()
{
#if defined(TETON_PARTITIONING)
   CALI_CXX_MARK_FUNCTION;

   if (doPartitioning())
   {
      conduit::Node &blueprint = getMeshBlueprint();
      conduit::Node &part = getMeshBlueprintPart();
#if defined(CLEANUP_PARTITION_TOPOLOGY)
      // Totally clear out the partitioned mesh.
      part.reset();
#endif
      // Remove some fields that we added to the original mesh.
      conduit::Node &fields = blueprint["fields"];
      std::vector<std::string> removals{PARTITION_FIELD, PARTITION_FIELD_BOUNDARY};
      for (conduit::index_t i = 0; i < fields.number_of_children(); i++)
      {
         if (fields[i].name().find("original_element_ids") != std::string::npos)
            removals.push_back(fields[i].name());
         if (fields[i].name().find("original_vertex_ids") != std::string::npos)
            removals.push_back(fields[i].name());
         // Remove any field that begins with PREFIX.
         if (fields[i].name().find(PREFIX) == 0)
            removals.push_back(fields[i].name());
      }
      for (const auto &name : removals)
      {
         if (fields.has_child(name))
         {
            fields.remove(name);
         }
      }
   }
#endif
}

std::string Teton::makeTestNode(conduit::Node &n,
                                conduit::Node &datastore,
                                conduit::Node &bp,
                                conduit::Node &options,
                                int flags)
{
#if defined(PARTITION_DEBUG)
   utilities::Banner b(mCommunicator, "test");
#endif
   const conduit::Node &main_topo = getMainTopology(bp);
   std::string mainTopologyName(main_topo.name());
   conduit::index_t nzones = conduit::blueprint::mesh::utils::topology::length(main_topo);
   conduit::index_t ngroups = options.fetch_existing("quadrature/num_groups").to_index_t();
   conduit::index_t nangles = options.fetch_existing("quadrature/num_angles").to_index_t();

   // Make a node that we'll check for validity. We make the node ourselves so
   // it only has the things we want in it.
   const std::vector<std::string> names{
      // inputs
      "fields/thermo_density",
      "fields/electron_temperature",
      "fields/radiation_temperature",
      "fields/electron_number_density",
      "fields/electron_specific_heat",
      "fields/absorption_opacity",
      "fields/scattering_opacity",
      // Some outputs that may have been registered as fields
      "fields/radiation_energy_density",
      "fields/electron_energy_deposited",
      "fields/radiation_force_x",
      "fields/radiation_force_y",
      "fields/radiation_force_z",
      "fields/radiation_force_r",
   };

   for (const auto &name : names)
   {
      if (bp.has_path(name))
         n[name].set_external(bp.fetch_existing(name));
   }

   // Get the radiation flux outputs.
   setRadiationFlux();
   for (const auto &name : getRadiationFluxFields())
   {
      if (!n.has_path(name) && bp.has_path(name))
         n[name].set_external(bp.fetch_existing(name));
   }

   // Get some radiation force outputs.
   if (((flags & Test_RadiationForceDensity) > 0) && bp.has_path(field_values(FIELD_CORNER_VOLUME_SUMS)))
   {
      conduit::index_t nnodes = conduit::blueprint::mesh::coordset::length(bp.fetch_existing("coordsets/coords"));
      n["fields/__result__radiation_force_x/association"] = "element";
      n["fields/__result__radiation_force_x/topology"] = mainTopologyName;
      n["fields/__result__radiation_force_x/values"].set(conduit::DataType::float64(nnodes));
      double *fx = n["fields/__result__radiation_force_x/values"].value();
      n["fields/__result__radiation_force_y/association"] = "element";
      n["fields/__result__radiation_force_y/topology"] = mainTopologyName;
      n["fields/__result__radiation_force_y/values"].set(conduit::DataType::float64(nnodes));
      double *fy = n["fields/__result__radiation_force_y/values"].value();
      n["fields/__result__radiation_force_z/association"] = "element";
      n["fields/__result__radiation_force_z/topology"] = mainTopologyName;
      n["fields/__result__radiation_force_z/values"].set(conduit::DataType::float64(nnodes));
      double *fz = n["fields/__result__radiation_force_z/values"].value();
      getRadiationForceDensity(fx, fy, fz);
   }

   if ((flags & Test_ZonalPsi) > 0)
   {
      // Call getZonalPsi and get the values out.
      auto nvalues = nzones * ngroups * nangles;
      n["fields/__result__zonal_psi/association"] = "element";
      n["fields/__result__zonal_psi/topology"] = mainTopologyName;
      n["fields/__result__zonal_psi/values"].set(conduit::DataType::float64(nvalues));
      double *values_ptr = n["fields/__result__zonal_psi/values"].value();
      getZonalPsi(nangles, values_ptr);
   }

   // NOTE: getMaterialTemperature is not backed by a Conduit field with the right length.
   // So, with partitioning, we can't ask for nzones values without going out of bounds
   if ((flags & Test_MaterialTemperature) > 0)
   {
      // Call getMaterialTemperature and get the values out.
      n["fields/__result__material_temperature/association"] = "element";
      n["fields/__result__material_temperature/topology"] = mainTopologyName;
      n["fields/__result__material_temperature/values"].set(conduit::DataType::float64(nzones));
      double *mt = n["fields/__result__material_temperature/values"].value();
      for (int zid = 0; zid < nzones; zid++)
         mt[zid] = getMaterialTemperature(zid + 1);
   }

   if ((flags & Test_RadiationTemperature) > 0)
   {
      // Call getRadiationTemperature and get the values out.
      n["fields/__result__radiation_temperature/association"] = "element";
      n["fields/__result__radiation_temperature/topology"] = mainTopologyName;
      n["fields/__result__radiation_temperature/values"].set(conduit::DataType::float64(nzones));
      double *rt = n["fields/__result__radiation_temperature/values"].value();
      for (int zid = 0; zid < nzones; zid++)
         rt[zid] = getRadiationTemperature(zid + 1);
   }

   if ((flags & Test_RadiationDeposited) > 0)
   {
      n["fields/__result__rad_energy_deposited/association"] = "element";
      n["fields/__result__rad_energy_deposited/topology"] = mainTopologyName;
      n["fields/__result__rad_energy_deposited/values"].set(conduit::DataType::float64(nzones));
      double *red = n["fields/__result__rad_energy_deposited/values"].value();
      for (int z = 0; z < nzones; z++)
         red[z] = getRadiationDeposited(z + 1);
   }

   if ((flags & Test_ReconstructPsi) > 0)
   {
      double erad = 0.;
      n["fields/__result__reconstructpsi_radEnergyDensity/association"] = "element";
      n["fields/__result__reconstructpsi_radEnergyDensity/topology"] = mainTopologyName;
      n["fields/__result__reconstructpsi_radEnergyDensity/values"].set(conduit::DataType::float64(nzones * ngroups));
      double *psired = n["fields/__result__reconstructpsi_radEnergyDensity/values"].value();
      reconstructPsi(&erad, psired);
      n["__result__reconstructpsi_erad"] = erad;
   }

   // Add in the options and edits so we can check their values too.
   n["options"].set_external(options);
   n["rtedits"].set_external(datastore.fetch_existing("rtedits"));

   int cycle = bp["state/cycle"].to_int();
   std::stringstream ss;
   ss << "_cycle=" << cycle << "_rank=" << mRank << "_a=" << nangles << "_g=" << ngroups << "_z=" << nzones;

   return ss.str();
}

} // namespace Teton
