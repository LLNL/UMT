#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <fstream>

#include "conduit/conduit_relay.hpp"
#include "conduit/conduit_relay_mpi_io_blueprint.hpp"
#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_relay_config.h"

#include "TetonInterface.hh"
#include "TetonConduitInterface.hh"
#include "TetonBlueprint.hh"
#include "dbc_macros.h"

#if defined(TETON_USE_CUDA)
#include "cuda.h"
#include "cuda_runtime_api.h"
#endif

#if defined(TETON_ENABLE_CALIPER)
#include "caliper/cali.h"
#else
#define CALI_MARK_BEGIN(label)
#define CALI_MARK_END(label)
#endif

extern "C" {
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

Teton::~Teton()
{
   conduit::Node &options = getOptions();
   bool enableNLTE = options.fetch_existing("size/enableNLTE").as_int();

   teton_destructmeshdata(&enableNLTE);
}

void Teton::initialize(MPI_Comm communicator, bool fromRestart)
{
   int myRank;
   CALI_MARK_BEGIN("Teton_Initialize");

   MPI_Fint fcomm = MPI_Comm_c2f(communicator);

   MPI_Comm_rank(communicator, &myRank);

   conduit::Node &datastore = getDatastore();
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();

   int verbose = 0;
   if (options.has_path("verbose"))
   {
      verbose = options.fetch_existing("verbose").value();
      if (verbose && myRank == 0)
         std::cout << "Teton: setting verbosity to " << verbose << std::endl;
   }
   else
   {
      options["verbose"] = verbose;
   }

   if (verbose >= 2)
   {
      std::cerr << "Teton: Dump copy of input..." << std::endl;
      conduit::relay::io::save(options, "parameters_input_" + std::to_string(myRank) + ".conduit_json", "conduit_json");
      conduit::relay::io::save(options, "parameters_input_" + std::to_string(myRank) + ".json", "json");
// The conduit install in /usr/workspace/teton/libraries/mapp has this added. Remove when new conduit version comes out.
#if defined(CONDUIT_CUSTOMINSTALL_SUPPORT_PAINT_ADJSET)
      conduit::blueprint::mesh::paint_adjset("main_adjset", "main_adjset", blueprint);
#endif
      conduit::relay::io::save(blueprint, "mesh_input_" + std::to_string(myRank) + ".conduit_json", "conduit_json");
      conduit::relay::io::save(blueprint, "mesh_input_" + std::to_string(myRank) + ".json", "json");
   }

   // Create secondary (corner) mesh topology and connectivity arrays.
   CALI_MARK_BEGIN("Teton_Construct_Corner_Mesh");
   TetonBlueprint blueprintHelper(blueprint, options);
   blueprintHelper.OutputTetonMesh(myRank, communicator);
   CALI_MARK_END("Teton_Construct_Corner_Mesh");

   if (verbose >= 2)
   {
// The conduit install in /usr/workspace/teton/libraries/mapp has this added. Remove when new conduit version comes out.
#if defined(CONDUIT_CUSTOMINSTALL_SUPPORT_PAINT_ADJSET)
      conduit::blueprint::mesh::paint_adjset("main_corner", "corner_adjset", blueprint);
#endif
      std::cerr << "Teton: Dump blueprint with generated topologies..." << std::endl;
#if defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
      std::string file_protocol = "hdf5";
      // TODO: fix
      // REMOVE uncomment below //
      // REMOVE uncomment below //
      //conduit::relay::mpi::io::blueprint::save_mesh(blueprint,
      //                              "./blueprint_mesh",
      //                              file_protocol,
      //                              communicator);
#endif

      conduit::relay::io::save(
          options, "parameters_post_processed_" + std::to_string(myRank) + ".conduit_json", "conduit_json");
      conduit::relay::io::save(blueprint, "mesh_input_post_processed_" + std::to_string(myRank) + ".json", "json");
   }

   // Set the Dt controls early, as tfloor is set here and needed by constructSize in the line below.
   constructDtControls();

   constructSize(myRank);
   constructMemoryAllocator();
   constructQuadrature();
   constructBoundaries(myRank);
   setSourceProfiles();
   teton_constructgeometry();
   setMeshConnectivity();

   MPI_Barrier(communicator);

   teton_setoppositeface(); //Prerequisite for calling setMeshSizeAndPositions()
   setCommunication();      //Prerequisite for calling setMeshSizeAndPositions()

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

   if (blueprint.has_path("fields/absorption_opacity/values"))
   {
      updateOpacity();
   }

   // Store mesh data needed for corner forces and update
   // of the mesh coordinates. In particular, the array corner_to_vertex
   // is stored, which associates the Teton corner id with the mesh vertex id
   storeMeshData();

   CALI_MARK_END("Teton_Initialize");
}

void Teton::storeMeshData()
{
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();

   // To compute the radiation forces, Teton needs to hang on to
   // this connectivity array
   if (blueprint.has_path("arrays/corner_to_vertex"))
   {
      int ncornr = options.fetch_existing("size/ncornr").to_int();
      mCornerToVertex.resize(ncornr);
      int *corner_to_vert_ptr = blueprint.fetch_existing("arrays/corner_to_vertex").value();
      for (int c = 0; c < ncornr; ++c)
      {
         // Store the vertex ID corresponding to this corner ID.
         // !!!! NOTE: THIS WILL NOT WORK FOR AMR MESHES !!!!
         mCornerToVertex[c] = corner_to_vert_ptr[c];
      }
   }
   // To compute the radiation forces, Teton also needs to hang on to
   // this connectivity array
   if (blueprint.has_path("relations/corner_to_zone"))
   {
      int ncornr = options.fetch_existing("size/ncornr").to_int();
      int *corner_to_zone_ptr = blueprint.fetch_existing("relations/corner_to_zone").value();
      mCornerToZone.resize(ncornr);
      for (int c = 0; c < ncornr; ++c)
      {
         mCornerToZone[c] = corner_to_zone_ptr[c];
      }
   }
   if (blueprint.has_path("arrays/zone_to_ncorners"))
   {
      int *zone_to_ncorner_ptr = blueprint.fetch_existing("arrays/zone_to_ncorners").value();
      int nzones = options.fetch_existing("size/nzones").to_int();
      mZoneToNCorners.resize(nzones);
      for (int zone = 0; zone < nzones; ++zone)
      {
         int ncorners = zone_to_ncorner_ptr[zone];
         mZoneToNCorners[zone] = ncorners;
      }
   }
   if (blueprint.has_path("arrays/zone_to_corners"))
   {
      int *zone_to_corner_ptr = blueprint.fetch_existing("arrays/zone_to_corners").value();
      int corner_counter = 0;
      int nzones = options.fetch_existing("size/nzones").to_int();
      for (int zone = 0; zone < nzones; ++zone)
      {
         int ncorners = mZoneToNCorners[zone];
         for (int c = 0; c < ncorners; ++c)
         {
            int corner = zone_to_corner_ptr[corner_counter];
            mZoneToCorners.push_back(corner);
            corner_counter += 1;
         }
      }
   }
}

void Teton::constructBoundaries(int rank)
{
   conduit::Node &options = getOptions();

   int nrefl = options.fetch_existing("boundary_conditions/num_reflecting").value();
   int nvac = options.fetch_existing("boundary_conditions/num_vacuum").value();
   int nsrc = options.fetch_existing("boundary_conditions/num_source").value();
   int num_comm = options.fetch_existing("boundary_conditions/num_comm").value();

   teton_constructboundary(&nrefl, &nvac, &nsrc, &num_comm);

   int numBCTotal = options.fetch_existing("boundary_conditions/num_total").value();
   int *BCTypeInt64 = options.fetch_existing("boundary_conditions/type").value();
   int *BCCornerFaces64 = options.fetch_existing("boundary_conditions/corner_face_ids").value();
   int *BCNeighborID64 = options.fetch_existing("boundary_conditions/neighbor_ids").value();

   // Ask conduit team how to pull out these items as integer pointers, instead of needing to do this.
   std::vector<int> BCTypeInt(BCTypeInt64, BCTypeInt64 + numBCTotal);
   std::vector<int> BCCornerFaces(BCCornerFaces64, BCCornerFaces64 + numBCTotal);
   std::vector<int> BCNeighborID(BCNeighborID64, BCNeighborID64 + numBCTotal);

   TETON_VERIFY_C(rank, (numBCTotal > 0), "No boundary conditions defined.");

   teton_addboundary(&numBCTotal, &BCTypeInt[0], &BCCornerFaces[0], &BCNeighborID[0]);
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

void Teton::constructSize(int rank)
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

   teton_constructsize(&rank,
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
      umpire_host_pinned_pool_allocator_id
          = options.fetch_existing("memory_allocator/umpire_host_allocator_id").value();
   }
   if (options.has_path("memory_allocator/umpire_device_allocator_id"))
   {
      umpire_device_pool_allocator_id = options.fetch_existing("memory_allocator/umpire_device_allocator_id").value();
   }

   teton_constructmemoryallocator(&umpire_host_pinned_pool_allocator_id, &umpire_device_pool_allocator_id);
}

void Teton::dump(int cycle, std::string path)
{
   conduit::Node &blueprint = getMeshBlueprint();

   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// This is defined in conduit_relay_config.h
#if defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
   std::string file_protocol = "hdf5";
#else
   std::string file_protocol = "conduit_json";
#endif

   conduit::relay::io::save(blueprint,
                            path + "/blueprint_mesh.rank" + std::to_string(rank) + ".cycle" + std::to_string(cycle) + "." + file_protocol,
                            file_protocol);
}
// ------------------------------------------------------------
//   step() - advance one cycle
// ------------------------------------------------------------

double Teton::step(int cycle)
{
   conduit::Node &datastore = getDatastore();
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();

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

   // TODO: need this for mesh, material, and opacity changes
   // Update zone vertex coordinates after hydro
   // Calling updateMeshPositions will cause the volume difference from
   // last cycle to current cycle, that Teton tracks, to be updated.
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   // Having trouble with UMT with this routine, with mfem meshes.
   // Need to troubleshoot later, but this is lower priority as UMT
   // problems don't have mesh motion.
   updateMeshPositions();
#endif
   //setMeshVelocity();
   // This updates the material properties (other than the opacities)

   // TODO Add something better than this to check for whether or not
   // new field values have been provided.
   if (blueprint.has_path("fields/thermo_density/values"))
   {
      setMaterials();
   }

   if (blueprint.has_path("fields/absorption_opacity/values"))
   {
      updateOpacity();
   }

   // ------------------------------------------------------------
   // Run the step
   // ------------------------------------------------------------

   // Set the time step information
   // A host code can either set these values in conduit, or can
   // setTimeStep() and that function will add these entries.
   // cycle = options.fetch_existing("iteration/cycle").value();
   double dtrad = options.fetch_existing("iteration/dtrad").value();
   double timerad = options.fetch_existing("iteration/timerad").value();
   double tfloor = options.fetch_existing("iteration/tfloor").value();

   teton_settimestep(&cycle, &dtrad, &timerad, &tfloor);

   // Main function in Teton to take a radiation step
   teton_radtr();

   // Update the radiation force (if the field is present)
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   // Note that, if the radiation force is present, there will always be
   // a z component in 2D or 3D ((r,z) or (x,y,z) coordinates).
   // TODO: fix when 1D is added
   bool has_rad_force = blueprint.has_path("fields/radiation_force_z");
   std::string rad_force_type;
   if (has_rad_force) 
   {
      rad_force_type = blueprint["fields/radiation_force_z/association"].as_string();
   }
   if (has_rad_force && rad_force_type == "element") 
   {
      updateZonalRadiationForce();
   }
   if (has_rad_force && rad_force_type != "element") 
   {
      updateRadiationForce();
   }

   // Update the radiation energy deposited to the material
   if (blueprint.has_path("fields/electron_energy_deposited/values")) 
   {
      double *electron_energy_deposited = blueprint.fetch_existing("fields/electron_energy_deposited/values").value();
      getRadEnergyDeposited(electron_energy_deposited);
   }

   // Update the radiation energy density
   if (blueprint.has_path("fields/radiation_energy_density/values")) 
   {
      double *radiation_energy_density = blueprint.fetch_existing("fields/radiation_energy_density/values").value();
      teton_getradiationenergydensity(radiation_energy_density);
   }
#endif

   // Updates the counters and statistics needed by teton_printedits()
   // Also provides a copy of the corner temp field, but this use case has
   // been deprecated for a while.  For now, we pass Teton its own corner temp
   // field so we can still call this function.
   double *tec = datastore.fetch_existing("material/Tec").value();
   teton_rtedit(tec);

   // Compute the recommended time step
   teton_dtnew(&maxOSComptonChangeCorner, &maxOSComptonChange);

   // put Teton's various edits in to its internal conduit node
   teton_publishedits(&mDTrad);

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
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
   double *spectrumAngleBinBoundaryList_ptr
       = options.fetch_existing("boundary_edits/spectrumAngleBinBoundaryList").value();
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
   int group = 0, offset = 0;
   int D_ngr6;
   std::vector<double> gnu;
   std::vector<int> quaddef;


   D_ngr6 = std::max(6 * (ngr + 1), 0);

   // Make sure there is room for this group
   gnu.resize(ngr + 1);
   quaddef.resize(D_ngr6);

// The Fortran expects the quadrature definition to be provided per-group, but should all have the
// same definition.
   for (group = 0; group < ngr; ++group)
   {
      gnu[group] = gnu_vals[group];
      gnu[group + 1] = gnu_vals[group + 1];
      quaddef[offset] = qtype;
      quaddef[offset + 1] = qorder;
      quaddef[offset + 2] = npolar;
      quaddef[offset + 3] = nazimu;
      quaddef[offset + 4] = paxis;
      quaddef[offset + 5] = -1; //Number of total angles (output).  The Fortran populates this value.
      offset = offset + 6;
   }

// Configure the quadrature for the grey acceleration ( reduced # angles and groups ).
   quaddef[offset] = 1;
   quaddef[offset + 1] = gtaOrder;
   quaddef[offset + 2] = 1;
   quaddef[offset + 3] = 1;
   quaddef[offset + 4] = 1;
   quaddef[offset + 5] = -1; // Number of total angles (output).  The Fortran populates this.

   teton_constructquadrature(&nSetsMaster, &nSets, &quaddef[0], &gnu[0]);

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

   int nsrc = options.fetch_existing("boundary_conditions/num_source").value();

   for (int j = 0; j < nsrc; ++j)
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
         bool blackBody
             = options.fetch_existing(top + "blackBody").to_int(); // Conduit doesn't support a 'bool' data type.
         bool isotropic
             = options.fetch_existing(top + "isotropic").to_int(); // Conduit doesn't support a 'bool' data type.

         teton_addprofile(
             &NumTimes, &NumValues, &Multiplier, &blackBody, &isotropic, &Times[0], &Values[0], &TetonProfileID);
      }

      // Save the TetonProfileID for later use:
      options[top + "TetonProfileID"] = TetonProfileID;
   }

   areSourceProfilesSet = true;
}

void Teton::setMeshSizeAndPositions()
{
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();

   int nzones = options.fetch_existing("size/nzones").value();
   double *zone_verts_ptr = blueprint.fetch_existing("arrays/zone_verts").value();
   //int *ncorners_ptr = blueprint.fetch_existing("fields/ncorners").value();
   int *ncorners_ptr = blueprint.fetch_existing("arrays/zone_to_ncorners").value();
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

   return;
}

void Teton::setMeshVelocity()
{
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();

   int nzones = options.fetch_existing("size/nzones").value();
   // TODO: change this to conform to blueprint standard
   double *velocities_ptr = blueprint.fetch_existing("fields/velocity_at_corners").as_double_ptr();
   int *ncorners_ptr = blueprint.fetch_existing("arrays/zone_to_ncorners").value();
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
   conduit::Node &blueprint = getMeshBlueprint();

   int nsfaces;
   int *shared_faces_ptr = nullptr;

   if (blueprint.has_path("shared_boundaries/nsfaces"))
   {
      nsfaces = blueprint.fetch_existing("shared_boundaries/nsfaces").to_int();
      if (nsfaces > 0)
      {
         shared_faces_ptr = blueprint.fetch_existing("shared_boundaries/shared_faces").value();
      }
   }
   else // if (options.has_path("shared_boundaries/nsfaces"))
   {
      // For backward compatbility
      nsfaces = options.fetch_existing("shared_boundaries/nsfaces").value();
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
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();

   int *connectivity_ptr = blueprint.fetch_existing("teton/arrays/corner_connectivity").value();

   int nzones = options.fetch_existing("size/nzones").value();

   // ADDED //
   int connect_off_set = 0;
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

   blueprint.remove("teton/arrays/corner_connectivity");
}

void Teton::setMaterials()
{
   conduit::Node &datastore = getDatastore();
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();

   int nzones = options.fetch_existing("size/nzones").value();

   double *density_ptr = blueprint.fetch_existing("fields/thermo_density/values").value();
   double *cv_ptr = blueprint.fetch_existing("fields/electron_specific_heat/values").value();
   double *tez_ptr = blueprint.fetch_existing("fields/electron_temperature/values").value();
   double *trz_ptr = blueprint.fetch_existing("fields/radiation_temperature/values").value();
   double *nez_ptr = blueprint.fetch_existing("fields/electron_number_density/values").value();

   // Really the effective electron specific energy source.
   if (blueprint.has_path("fields/specific_energy_source"))
   {
      double *matSource = blueprint.fetch_existing("fields/specific_energy_source/values").value();
      teton_setmaterialsource(matSource);
   }

   // Initialize arrays to handle multi-material zones
   double *tec = datastore.fetch_existing("material/Tec").value();
   teton_initmaterial(tec);

   for (int zone = 0; zone < nzones; ++zone)
   {
      double scm = 0.;
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
   conduit::Node &blueprint = getMeshBlueprint();

   int ngr = options.fetch_existing("quadrature/num_groups").to_int(); //coerce from unsigned int or size_t
   int ig;
   std::vector<double> siga_loc;
   std::vector<double> sigs_loc;
   siga_loc.resize(ngr);
   sigs_loc.resize(ngr);

   bool useInternalSigmaS = false;
   if (options.has_path("compton/use_internal_sigma_s"))
   {
      useInternalSigmaS = options.fetch_existing("compton/use_internal_sigma_s").as_int();
   }
   bool useTableSigmaS = (not useInternalSigmaS);

   // zero out opacities
   teton_initopacity();

   int nzones = options.fetch_existing("size/nzones").value();
   double *absorption_opacity_ptr = blueprint.fetch_existing("fields/absorption_opacity/values").value();
   double *scattering_opacity_ptr = nullptr;
   if (useTableSigmaS)
   {
      scattering_opacity_ptr = blueprint.fetch_existing("fields/scattering_opacity/values").value();
   }

   // Initialize opacities to handle multi-material zones
   int offset = 0;
   for (int zone = 0; zone < nzones; zone++)
   {
      for (ig = 0; ig < ngr; ig++)
      {
         siga_loc[ig] = absorption_opacity_ptr[offset];
         sigs_loc[ig] = useTableSigmaS ? scattering_opacity_ptr[offset] : 0.;
         offset += 1;
      }
      int zoneID = zone + 1;
      teton_setopacity(&zoneID, &siga_loc[0], &sigs_loc[0], &useTableSigmaS);
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

double Teton::getMaterialTemperature(int zone)
{
   double matTemp;
   teton_getmaterialtemperature(&zone, &matTemp);
   return matTemp;
}

conduit::Node &Teton::getDatastore()
{
   return *teton_get_datastore_cptr();
}

double Teton::getRadiationTemperature(int zone)
{
   double radTemp;
   teton_getradiationtemperature(&zone, &radTemp);
   return radTemp;
}

double Teton::getRadiationDeposited(int zone)
{
   double eDep, tRad;
   teton_getradiationdeposited(&zone, &eDep, &tRad);
   return eDep;
}

// TODO - Can we deprecate this function version?  Why isn't timerad being provided?
void Teton::setTimeStep(int cycle, double dtrad)
{
   double timerad = 0.0;

   setTimeStep(cycle, dtrad, timerad);
}

void Teton::setTimeStep(int cycle, double dtrad, double timerad)
{
   conduit::Node &options = getOptions();

   options["iteration/cycle"] = cycle;
   options["iteration/dtrad"] = dtrad;
   options["iteration/timerad"] = timerad;
}

void Teton::updateMeshPositions()
{
   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();

   int nzones = mZoneToNCorners.size();
   int ndim = options.fetch_existing("size/ndim").value();

   // Best practice for updating the mesh zone vertices is for codes to update the blueprint coordinate arrays.
   //
   // MARBL is providing a zone_verts array directly, at the moment.  This array is a listing of the zone vertices,
   // in the same order as the corners in each zone.  Use that if it is present, otherwise generate it from the
   // blueprint coords.
   if (!blueprint.has_path("arrays/zone_verts"))
   {
      int corner_counter  = 0;
      int zoneVertsSize   = 0;

      const double *m_x = nullptr;
      const double *m_y = nullptr;
      const double *m_r = nullptr;
      const double *m_z = nullptr;
      if (ndim == 1)
      {
         m_x = blueprint.fetch_existing("coordsets/coords/values/x").value();
      }
      else if (ndim == 2)
      {
         if (blueprint.has_path("coordsets/coords/values/r"))
         {
            m_r = blueprint.fetch_existing("coordsets/coords/values/r").value();
         }
         else
         { // assuming zr ordering for marbl/ares as fallback.  EVERYONE should just specify r and z directly.
            m_r = blueprint.fetch_existing("coordsets/coords/values/y").value();
         }
         if (blueprint.has_path("coordsets/coords/values/z"))
         {
            m_z = blueprint.fetch_existing("coordsets/coords/values/z").value();
         }
         else
         { // assuming zr ordering for marbl/ares as fallback.  EVERYONE should just specify r and z directly.  For now, issue a warning.
            m_z = blueprint.fetch_existing("coordsets/coords/values/x").value();
         }
      }
      else if (ndim == 3)
      {
         m_x = blueprint.fetch_existing("coordsets/coords/values/x").value();
         m_y = blueprint.fetch_existing("coordsets/coords/values/y").value();
         m_z = blueprint.fetch_existing("coordsets/coords/values/z").value();
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
            int v      = mCornerToVertex[corner];

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
      blueprint["arrays/zone_verts"].set(&zoneVerts[0], zoneVerts.size());
   }

   setMeshSizeAndPositions();

   // Update Teton geometry
   teton_getvolume();

   // We're done updating the node positions, we shouldn't need zone_verts anymore.
   blueprint.remove("arrays/zone_verts");

   return;
}

// NOTE: the Vectors RadiationForceXTotal, ..., must
//       already be sized to the number of mesh vertices
void Teton::getRadiationForceDensity(double *RadiationForceDensityX,
                                     double *RadiationForceDensityY,
                                     double *RadiationForceDensityZ)
{
   // Compute the radiation force internally in Teton
   // for each zone and corner
   teton_setradiationforce();

   conduit::Node &options = getOptions();
   int ndim = options.fetch_existing("size/ndim").value();
   int maxCorner = options.fetch_existing("size/maxCorner").value();
   int nzones = options.fetch_existing("size/nzones").value();
   int nverts = options.fetch_existing("size/nverts").value();
   std::vector<double> RadiationForce(ndim * maxCorner);
   std::vector<double> CornerVolumes(maxCorner);
   std::vector<double> CornerVolumeSumsAtVertex(nverts);
   int corner_counter = 0;

   // TODO: finish this case
   if (ndim == 1)
   {
      std::cerr << "1D getRadiationForceDensity not yet implemented! Teton is exiting . . ." << std::endl;
      exit(1);
   }

   for (int v = 0; v < nverts; ++v)
   {
      CornerVolumeSumsAtVertex[v] = 0.0;
      RadiationForceDensityX[v] = 0.0;
      RadiationForceDensityY[v] = 0.0;
      if (ndim == 3)
         RadiationForceDensityZ[v] = 0.0;
   }

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
         RadiationForceDensityX[vertexID] += RadiationForce[c * ndim + 0];
         RadiationForceDensityY[vertexID] += RadiationForce[c * ndim + 1];
         if (ndim == 3)
            RadiationForceDensityZ[vertexID] += RadiationForce[c * ndim + 2];
         CornerVolumeSumsAtVertex[vertexID] += CornerVolumes[c];
      }
   }

   for (int v = 0; v < nverts; ++v)
   {
      RadiationForceDensityX[v] /= CornerVolumeSumsAtVertex[v];
      RadiationForceDensityY[v] /= CornerVolumeSumsAtVertex[v];
      if (ndim == 3)
         RadiationForceDensityZ[v] /= CornerVolumeSumsAtVertex[v];
   }
}

void Teton::updateRadiationForce()
{
   // Compute the radiation force internally in Teton
   // for each zone and corner
   teton_setradiationforce();

   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();
   int ndim = options.fetch_existing("size/ndim").value();
   int maxCorner = options.fetch_existing("size/maxCorner").value();
   int nzones = options.fetch_existing("size/nzones").value();
   int nverts = options.fetch_existing("size/nverts").value();
   std::vector<double> RadiationForce(ndim * maxCorner);
   int corner_counter = 0;

   double *radiation_force_x;
   double *radiation_force_y;
   double *radiation_force_z;
   if (ndim == 2)
   {
      radiation_force_x = blueprint.fetch_existing("fields/radiation_force_r/values").value();
      radiation_force_y = blueprint.fetch_existing("fields/radiation_force_z/values").value();
   }
   else if (ndim == 3) 
   {
      radiation_force_x = blueprint.fetch_existing("fields/radiation_force_x/values").value();
      radiation_force_y = blueprint.fetch_existing("fields/radiation_force_y/values").value();
      radiation_force_z = blueprint.fetch_existing("fields/radiation_force_z/values").value();
   }
   // TODO: finish this case
   else 
   {
      std::cerr << "1D updateRadiationForce not yet implemented! Teton is exiting . . ." << std::endl;
      exit(1);
   }

   for (int v = 0; v < nverts; ++v)
   {
      radiation_force_x[v] = 0.0;
      radiation_force_y[v] = 0.0;
      if (ndim == 3)
         radiation_force_z[v] = 0.0;
   }

   for (int zone = 0; zone < nzones; ++zone)
   {
      // Get the radiation force and volume on each corner of each zone
      int zoneID = zone + 1;
      teton_getradiationforce(&zoneID, &RadiationForce[0]);

      // Average the radiation force around vertices
      int ncorners = mZoneToNCorners[zone];
      for (int c = 0; c < ncorners; ++c)
      {
         int cornerID = mZoneToCorners[corner_counter];
         int vertexID = mCornerToVertex[cornerID];
         corner_counter += 1;
         radiation_force_x[vertexID] += RadiationForce[c * ndim + 0];
         radiation_force_y[vertexID] += RadiationForce[c * ndim + 1];
         if (ndim == 3)
            radiation_force_z[vertexID] += RadiationForce[c * ndim + 2];
      }
   }
}

void Teton::updateZonalRadiationForce()
{
   // Compute the radiation force internally in Teton
   // for each zone and corner
   teton_setradiationforce();

   conduit::Node &options = getOptions();
   conduit::Node &blueprint = getMeshBlueprint();
   int ndim = options.fetch_existing("size/ndim").value();
   int maxCorner = options.fetch_existing("size/maxCorner").value();
   int nzones = options.fetch_existing("size/nzones").value();
   std::vector<double> RadiationForce(ndim * maxCorner);
   int corner_counter = 0;

   double *radiation_force_x;
   double *radiation_force_y;
   double *radiation_force_z;
   if (ndim == 2)
   {
      radiation_force_x = blueprint.fetch_existing("fields/radiation_force_r/values").value();
      radiation_force_y = blueprint.fetch_existing("fields/radiation_force_z/values").value();
   }
   else if (ndim == 3) 
   {
      radiation_force_x = blueprint.fetch_existing("fields/radiation_force_x/values").value();
      radiation_force_y = blueprint.fetch_existing("fields/radiation_force_y/values").value();
      radiation_force_z = blueprint.fetch_existing("fields/radiation_force_z/values").value();
   }
   // TODO: finish this case
   else 
   {
      std::cerr << "1D updateZonalRadiationForce not yet implemented! Teton is exiting . . ." << std::endl;
      exit(1);
   }

   for (int zone = 0; zone < nzones; ++zone)
   {
      radiation_force_x[zone] = 0.0; 
      radiation_force_y[zone] = 0.0; 
   }
   if (ndim == 3)
   {
      for (int zone = 0; zone < nzones; ++zone)
      {
         radiation_force_z[zone] = 0.0; 
      }
   }

   for (int zone = 0; zone < nzones; ++zone)
   {
      // Get the radiation force and volume on each corner of each zone
      int zoneID = zone + 1;
      teton_getradiationforce(&zoneID, &RadiationForce[0]);
      // Sum the corner radiation force
      int ncorners = mZoneToNCorners[zone];
      for (int c = 0; c < ncorners; ++c)
      {
         radiation_force_x[zone] += RadiationForce[c * ndim + 0];
         radiation_force_y[zone] += RadiationForce[c * ndim + 1];
         if (ndim == 3)
            radiation_force_z[zone] += RadiationForce[c * ndim + 2];
      }
   }
}

void Teton::getRadEnergyDeposited(double *RadEnergyDeposited)
{
   conduit::Node &options = getOptions();
   int nzones = options.fetch_existing("size/nzones").value();
   for (int zone = 0; zone < nzones; ++zone)
   {
      RadEnergyDeposited[zone] = 0.0; 
   }
   for (int zone = 0; zone < nzones; ++zone)
   {
      // Get the radiation energy deposited 
      double rad_temp;
      int zoneID = zone + 1;
      teton_getradiationdeposited(&zoneID, &RadEnergyDeposited[zone], &rad_temp);
   }
}

void Teton::reconstructPsi(double *rad_energy, double *rad_energy_density)
{
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   teton_reconstructpsi(rad_energy, rad_energy_density);
#endif
}

void Teton::reconstructPsiFromdV()
{
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   teton_reconstructpsifromdv();
#endif
}

} // namespace Teton
