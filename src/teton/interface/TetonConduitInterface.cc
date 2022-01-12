#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <stdexcept>
#include <iostream>
#include <fstream>

#include "conduit/conduit_relay.hpp"
#include "conduit/conduit_blueprint.hpp"

#include "TetonInterface.hh"
#include "TetonConduitInterface.hh"

#if defined(TETON_USE_CUDA)
#include "cuda.h"
#include "cuda_runtime_api.h"
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
void Teton::initialize()
{
   int myRank;
   // TODO: remove MPI_COMM_WORLD. See
   // https://rzlc.llnl.gov/gitlab/deterministic-transport/TRT/Teton/-/issues/82
   MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);

   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   conduit::Node &datastore = getDatastore();

   int verbose = 0;
   if (datastore.has_path("verbose"))
   {
      verbose = datastore["verbose"].value();
      if (myRank == 0)
         std::cout << "Teton driver: setting verbosity to " << verbose << std::endl;
   }
   else
   {
      datastore["verbose"] = verbose;
   }

   double *cornerTe_ptr = mMeshBlueprint["fields/corner_temperature"].value();

   int ncornr = datastore["size/ncornr"].value();
   mCornerTemperature.resize(ncornr);

   int nzones = datastore["size/nzones"].value();
   int ngr = datastore["quadrature/num_groups"].to_int(); //coerse from unsigned int or size_t
   mRadiationEnergyDensity.resize(nzones * ngr);
   if (myRank == 0)
      std::cout << "Teton driver: constructing Size..." << std::endl;
   constructSize();
   if (myRank == 0)
      std::cout << "Teton driver: constructing memory allocator..." << std::endl;
   constructMemoryAllocator();
   if (myRank == 0)
      std::cout << "Teton driver: constructing quadrature..." << std::endl;
   constructQuadrature();
   if (myRank == 0)
      std::cout << "Teton driver: constructing boundaries..." << std::endl;
   constructBoundaries();
   if (myRank == 0)
      std::cout << "Teton driver: setting source profiles..." << std::endl;
   setSourceProfiles();
   if (myRank == 0)
      std::cout << "Teton driver: constructing geometry..." << std::endl;
   teton_constructgeometry();
   if (myRank == 0)
      std::cout << "Teton driver: setting mesh connectivity..." << std::endl;
   setMeshConnectivity();

   MPI_Barrier(MPI_COMM_WORLD);

   if (myRank == 0)
      std::cout << "Teton driver: setting opposite faces..." << std::endl;
   teton_setoppositeface(); //Prerequisite for calling setMeshSizeAndPositions()
   if (myRank == 0)
      std::cout << "Teton driver: setting comms..." << std::endl;
   setCommunication(); //Prerequisite for calling setMeshSizeAndPositions()

   if (myRank == 0)
      std::cout << "Teton driver: setting mesh sizes and positions..." << std::endl;
   setMeshSizeAndPositions();

   // This is an awful hack to make sure Z%VolumeOld is initialized
   teton_getvolume();
   teton_getvolume();

   if (myRank == 0)
      std::cout << "Teton driver: setting communication groups..." << std::endl;
   teton_setcommunicationgroup(&fcomm);

   if (myRank == 0)
      std::cout << "Teton driver: checking shared boundaries..." << std::endl;
   teton_checksharedboundary();

   bool enableNLTE = datastore["size/enableNLTE"].as_int();
   if (myRank == 0)
      std::cout << "Teton driver: constructing materials..." << std::endl;
   teton_constructmaterial(&enableNLTE);

   bool fromRestart = false;
   if (myRank == 0)
      std::cout << "Teton driver: constructing phase space sets..." << std::endl;
   teton_constructphasespacesets(&fromRestart);

   if (myRank == 0)
      std::cout << "Teton driver: initializing material fields..." << std::endl;
   setMaterials(cornerTe_ptr);

   mInternalComptonFlag = (int) tetonComptonFlag::invalid;
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   if (myRank == 0)
      std::cout << "Teton driver: constructing compton controller..." << std::endl;
   constructComptonControl();
#endif

   if (myRank == 0)
      std::cout << "Teton driver: constructing radiation intensity..." << std::endl;
   teton_constructradintensity();

   double EnergyRadiation = datastore["radiation_energy"].value();
   if (myRank == 0)
      std::cout << "Teton driver: init'ing fields..." << std::endl;
   teton_initteton(&EnergyRadiation, &cornerTe_ptr[0]);
   if (myRank == 0)
      std::cout << "Teton driver: constructing edits..." << std::endl;
   constructEdits();
   if (myRank == 0)
      std::cout << "Teton driver: constructing iteration controls..." << std::endl;
   constructIterationControls();
   if (myRank == 0)
      std::cout << "Teton driver: constructing dt controls..." << std::endl;
   constructDtControls();

   if (myRank == 0)
      std::cout << "Teton driver: updating opacities..." << std::endl;
   updateOpacity();

   // Store mesh data needed for corner forces and update
   // of the mesh coordinates. In particular, the array corner_to_vertex
   // is stored, which associates the Teton corner id with the mesh vertex id
   storeMeshData();
}

void Teton::storeMeshData()
{
   conduit::Node &datastore = getDatastore();

   // To compute the radiation forces, Teton needs to hang on to
   // this connectivity array
   if (mMeshBlueprint.has_path("topologies/mesh/corners/corner_to_vertex"))
   {
      int ncornr = datastore["size/ncornr"].to_int();
      mCornerToVertex.resize(ncornr);
      int *corner_to_vert_ptr = mMeshBlueprint["topologies/mesh/corners/corner_to_vertex"].value();
      for (int c = 0; c < ncornr; ++c)
      {
         // Store the vertex ID corresponding to this corner ID.
         // !!!! NOTE: THIS WILL NOT WORK FOR AMR MESHES !!!!
         mCornerToVertex[c] = corner_to_vert_ptr[c];
      }
   }
   if (mMeshBlueprint.has_path("topologies/mesh/elements/zone_to_corners"))
   {
      int *zone_to_corner_ptr = mMeshBlueprint["topologies/mesh/elements/zone_to_corners"].value();
      int offset = 0;
      int nzones = datastore["size/nzones"].to_int();
      mZoneToNCorners.resize(nzones);
      for (int zone = 0; zone < nzones; ++zone)
      {
         int ncorners = zone_to_corner_ptr[offset];
         offset += 1;
         mZoneToNCorners[zone] = ncorners;
         offset += ncorners;
      }
   }
}

void Teton::constructBoundaries()
{
   conduit::Node &datastore = getDatastore();

   int nrefl = datastore["boundary_conditions/num_reflecting"].value();
   int nvac = datastore["boundary_conditions/num_vacuum"].value();
   int nsrc = datastore["boundary_conditions/num_source"].value();
   int num_comm = datastore["boundary_conditions/num_comm"].value();

   teton_constructboundary(&nrefl, &nvac, &nsrc, &num_comm);

   int numBCTotal = datastore["boundary_conditions/num_total"].value();
   int *BCTypeInt64 = datastore["boundary_conditions/type"].value();
   int *BCCornerFaces64 = datastore["boundary_conditions/corner_face_ids"].value();
   int *BCNeighborID64 = datastore["boundary_conditions/neighbor_ids"].value();

   // Ask conduit team how to pull out these items as integer pointers, instead of needing to do this.
   std::vector<int> BCTypeInt(BCTypeInt64, BCTypeInt64 + numBCTotal);
   std::vector<int> BCCornerFaces(BCCornerFaces64, BCCornerFaces64 + numBCTotal);
   std::vector<int> BCNeighborID(BCNeighborID64, BCNeighborID64 + numBCTotal);

   teton_addboundary(&numBCTotal, &BCTypeInt[0], &BCCornerFaces[0], &BCNeighborID[0]);
}

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
void Teton::constructComptonControl()
{
   conduit::Node &datastore = getDatastore();

   if (datastore.has_path("compton/internalComptonFlag"))
   {
      mInternalComptonFlag = datastore["compton/internalComptonFlag"].value();
   }
   else
   {
      bool isoThomsonCompton = datastore["compton/use_thomas_compton"].as_int();
      bool internalCompton = datastore["compton/use_internal_compton"].as_int();

      if (isoThomsonCompton)
      {
         mInternalComptonFlag = (int) tetonComptonFlag::thomson;
      }
      else if (internalCompton)
      {
         bool isoBoltzCompton = datastore["compton/use_boltzmann_compton"].as_int();
         mInternalComptonFlag = isoBoltzCompton ? (int) tetonComptonFlag::boltzmann : (int) tetonComptonFlag::fp;
      }
      else
      {
         mInternalComptonFlag = (int) tetonComptonFlag::none;
      }
   }

   teton_constructcomptoncontrol(&mInternalComptonFlag);
}
#endif

void Teton::constructSize()
{
   conduit::Node &datastore = getDatastore();

   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   int nzones = datastore["size/nzones"].value();
   int ncornr = datastore["size/ncornr"].value();
   int nsides = datastore["size/nsides"].value();
   int nbelem = datastore["size/nbelem"].value();
   int maxcf = datastore["size/maxcf"].value();
   int maxCorner = datastore["size/maxCorner"].value();
   int ncomm = datastore["size/ncomm"].value();
   int ndim = datastore["size/ndim"].value();
   int ngr = datastore["quadrature/num_groups"].to_int(); // coerce from size_t or unsigned long
   int functionRNLTE = datastore["size/functionRNLTE"].value();
   double tfloor = datastore["size/tfloor"].value();
   double radForceMultiplier = datastore["size/radForceMultiplier"].value();
   double betaNLTE = datastore["size/betaNLTE"].value();
   double gammaNLTE = datastore["size/gammaNLTE"].value();
   bool dopplerShiftOn = datastore["size/DopplerShiftOn"].to_int();
   bool useNewNonLinearSolver = datastore["size/useNewNonLinearSolver"].as_int();
   bool useNewGTASolver = datastore["size/useNewGTASolver"].as_int();
   bool usePWLD = datastore["size/usePWLD"].as_int();
   bool useSurfaceMassLumping = datastore["size/useSurfaceMassLumping"].as_int();
   bool useGPU = datastore["size/useGPU"].as_int();
   bool useCUDASweep = datastore["size/useCUDASweep"].as_int();
   bool useCUDASolver = datastore["size/useCUDASolver"].as_int();
   int zoneBatchSize = datastore["size/zoneBatchSize"].value();
   int nConcurrentBatches = datastore["size/nConcurrentBatches"].value();
   int igeomToFortran = datastore["size/geomType"].value();

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
   conduit::Node &datastore = getDatastore();

   int umpire_host_pinned_pool_allocator_id = datastore["memory_allocator/umpire_host_allocator_id"].value();
   int umpire_device_pool_allocator_id = datastore["memory_allocator/umpire_device_allocator_id"].value();

   teton_constructmemoryallocator(&umpire_host_pinned_pool_allocator_id, &umpire_device_pool_allocator_id);
}

void Teton::dump(int cycle, std::string path, std::string name, bool allFields)
{
   conduit::Node &datastore = getDatastore();

   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (allFields)
   {
      conduit::relay::io::save(datastore,
                               path + "/" + name + ".rank" + std::to_string(rank) + ".cycle" + std::to_string(cycle)
                                   + ".hdf5",
                               "hdf5");
   }
   // Only dump a subset of fields commonly used for diagnostics
   else
   {
      conduit::Node subset;
      subset["geometry/RadEnergyDensity"] = datastore["geometry/RadEnergyDensity"];
      subset["material/Tec"] = datastore["material/Tec"];
      conduit::relay::io::save(subset,
                               path + "/" + name + ".rank" + std::to_string(rank) + ".cycle" + std::to_string(cycle)
                                   + ".hdf5",
                               "hdf5");
   }
}
// ------------------------------------------------------------
//   step() - advance one cycle
// ------------------------------------------------------------

double Teton::step(int cycle)
{
   conduit::Node &datastore = getDatastore();

   int verbose, noutrt, ninrt, ngdart, nNLIters, maxNLIters;
   int TrMaxZone, TeMaxZone, TrMaxProcess, TeMaxProcess;
   int maxOSComptonChangeCorner = 1;

   double TrMax, TeMax, EnergyRadiation;
   double PowerIncident, PowerEscape, PowerAbsorbed;
   double PowerEmitted, PowerExtSources, PowerCompton;
   double MatCoupTimeTotal, SweepTimeTotal, GPUSweepTimeTotal, GTATimeTotal;
   double RadtrTimeTotal, InitTimeTotal, FinalTimeTotal, timeNonRad = 0.0, timeOther = 0.0;
   double maxOSComptonChange = 0.0;
   double timerad = 0.0;

   double dtused = 0.0; // make this a class member, get first value from JSON input and update each step.
   verbose = 0;         // Deprecated, use teton_setverbose(int) instead.

   if (datastore.has_path("iteration/outerMaxIt"))
   {
      int i = datastore["iteration/outerMaxIt"].value();
      teton_adjust_temperature_maxits(&i);
   }
   if (datastore.has_path("iteration/greyMaxIt"))
   {
      int i = datastore["iteration/greyMaxIt"].value();
      teton_adjust_grey_maxits(&i);
   }
   if (datastore.has_path("iteration/incidentFluxMaxIt"))
   {
      int i = datastore["iteration/incidentFluxMaxIt"].value();
      teton_adjust_fluxexchange_maxits(&i);
   }
   if (datastore.has_path("iteration/innerNLMaxIt"))
   {
      int i = datastore["iteration/innerNLMaxIt"].value();
      teton_adjust_nonlinear_maxits(&i);
   }

   // Global tolerance that gets set first, so others can override it.
   // This is the prefered option to set, and will set the others
   // consistently.
   if (datastore.has_path("iteration/relativeTolerance"))
   {
      double x = datastore["iteration/relativeTolerance"].value();
      teton_adjust_relative_tolerance(&x);
   }

   if (datastore.has_path("iteration/outerTempRelTol"))
   {
      double x = datastore["iteration/outerTempRelTol"].value();
      teton_adjust_temperature_reltol(&x);
   }
   if (datastore.has_path("iteration/outerPhiRelTol"))
   {
      double x = datastore["iteration/outerPhiRelTol"].value();
      teton_adjust_radenergydensity_reltol(&x);
   }
   if (datastore.has_path("iteration/incidentFluxRelTol"))
   {
      double x = datastore["iteration/incidentFluxRelTol"].value();
      teton_adjust_fluxexchange_reltol(&x);
   }
   if (datastore.has_path("iteration/innerNLRelTol"))
   {
      double x = datastore["iteration/innerNLRelTol"].value();
      teton_adjust_fluxexchange_reltol(&x);
   }
   if (datastore.has_path("iteration/greyRelTol"))
   {
      double x = datastore["iteration/greyRelTol"].value();
      teton_adjust_fluxexchange_reltol(&x);
   }

   // ------------------------------------------------------------
   // Update the mesh positions, material info, and opacity info
   // ------------------------------------------------------------

   // TODO: need this for mesh, material, and opacity changes
   if (cycle > 0)
   {
      // Update zone vertex coordinates after hydro
      updateMeshPositions();
      //setMeshVelocity();
      // TODO: need to expose this
      // NOTE: mCornerTemperature gets updated via the call to teton_rtedit
      // after the call to teton_radtr
      // This updates the material properties (other than the opacities)
      setMaterials(mCornerTemperature.data());
      updateOpacity();
   }

   // ------------------------------------------------------------
   // Run the step
   // ------------------------------------------------------------

   // Let Teton set the timestep
   // TODO - add support using timesteps file if one present
   // double dtrad = datastore["iteration/dtrad"].value();
   // teton_settimestep(&cycle, &dtrad, &timerad, &tfloor);

   // Main function in Teton to take a radiation step
   teton_radtr();

   // Update corner temperature field in conduit node.
   // Also updates the counters and statistics needed by teton_printedits()
   teton_rtedit(mCornerTemperature.data());
   //   mMeshBlueprint["fields/corner_temperature"].set(mCornerTemperature.data(), mCornerTemperature.size());

   // Compute the recommended time step
   teton_dtnew(&maxOSComptonChangeCorner, &maxOSComptonChange);

   // put Teton's various edits in to its internal conduit node
   teton_publishedits(&dtrad);

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
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
   return dtrad;
}

void Teton::constructEdits()
{
   conduit::Node &datastore = getDatastore();

   int ngr = datastore["quadrature/num_groups"].to_int(); //coerce from unsigned int or size_t
   int numSpectrumAngleBins = datastore["boundary_edits/numSpectrumAngleBins"].value();

   double *spectrumAngleBinBoundaryList_ptr = datastore["boundary_edits/spectrumAngleBinBoundaryList"].value();
   double *RadPowerEscape_ptr = datastore["boundary_edits/RadPowerEscape"].value();
   double *RadPowerIncident_ptr = datastore["boundary_edits/RadPowerIncident"].value();
   double *PolarSectorPowerEscape_ptr = datastore["boundary_edits/PolarSectorPowerEscape"].value();

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
   conduit::Node &datastore = getDatastore();

   int nSets = datastore["quadrature/nSets"].value();
   int nSetsMaster = datastore["quadrature/nSetsMaster"].value();
   double *gnu_vals = datastore["quadrature/gnu"].value();
   int qtype = datastore["quadrature/qtype"].value();
   int qorder = datastore["quadrature/qorder"].value();
   int npolar = datastore["quadrature/npolar"].value();
   int nazimu = datastore["quadrature/nazimu"].value();
   int paxis = datastore["quadrature/paxis"].value();
   int ngr = datastore["quadrature/num_groups"].to_int(); //coerce from unsigned int or size_t
   int gtaOrder = datastore["quadrature/gtaorder"].value();
   int group = 0, offset = 0;
   int D_ngr6;
   std::vector<double> gnu;
   std::vector<int> quaddef;

   D_ngr6 = std::max(6 * (ngr + 1), 0);

   // Make sure there is room for this group
   gnu.resize(ngr + 1);
   quaddef.resize(D_ngr6);

   for (group = 0; group < ngr; ++group)
   {
      gnu[group] = gnu_vals[group];
      gnu[group + 1] = gnu_vals[group + 1];
      quaddef[offset] = qtype;
      quaddef[offset + 1] = qorder;
      quaddef[offset + 2] = npolar;
      quaddef[offset + 3] = nazimu;
      quaddef[offset + 4] = paxis;
      quaddef[offset + 5] = -1; //unused
      offset = offset + 6;
   }

   // Set quadrature definition for acceleration
   quaddef[offset] = 1;
   quaddef[offset + 1] = gtaOrder;
   quaddef[offset + 2] = 1;
   quaddef[offset + 3] = 1;
   quaddef[offset + 4] = 1;
   quaddef[offset + 5] = -1; // unused

   teton_constructquadrature(&nSetsMaster, &nSets, &quaddef[0], &gnu[0]);
}

void Teton::resetSourceProfiles()
{
   conduit::Node &datastore = getDatastore();

   if (!areSourceProfilesSet)
   {
      std::cerr << "setSourceProfiles must be called before resetSourceProfiles!" << std::endl;
      exit(1);
   }

   int nsrc = datastore["boundary_conditions/num_source"].value();

   for (int j = 0; j < nsrc; ++j)
   {
      std::string top = "sources/profile" + std::to_string(j + 1) + "/";
      int TetonProfileID = datastore[top + "TetonProfileID"].value();
      int NumTimes = datastore[top + "NumTimes"].value();
      int NumValues = datastore[top + "NumValues"].value();

      double *values_ptr = datastore[top + "Values"].value();
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
         double Multiplier = datastore[top + "Multiplier"].value();
         teton_resetprofile(&TetonProfileID, &Multiplier, Values);
      }
      else
      {
         double Multiplier = datastore[top + "Multiplier"].value();

         double *times_ptr = datastore[top + "Times"].value();
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
   conduit::Node &datastore = getDatastore();

   int nsrc = datastore["boundary_conditions/num_source"].value();

   for (int j = 0; j < nsrc; ++j)
   {
      std::string top = "sources/profile" + std::to_string(j + 1) + "/";
      int NumTimes = datastore[top + "NumTimes"].value();
      int NumValues = datastore[top + "NumValues"].value();

      double *values_ptr = datastore[top + "Values"].value();
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
         double Multiplier = datastore[top + "Multiplier"].value();
         teton_addprofile(&Multiplier, Values, &TetonProfileID);
      }
      else
      {
         double *times_ptr = datastore[top + "Times"].value();
         std::vector<double> Times(NumTimes);
         for (int k = 0; k < NumTimes; ++k)
         {
            Times[k] = times_ptr[k];
         }

         double Multiplier = datastore[top + "Multiplier"].value();
         bool blackBody = datastore[top + "blackBody"].to_int(); // Conduit doesn't support a 'bool' data type.
         bool isotropic = datastore[top + "isotropic"].to_int(); // Conduit doesn't support a 'bool' data type.

         teton_addprofile(
             &NumTimes, &NumValues, &Multiplier, &blackBody, &isotropic, &Times[0], &Values[0], &TetonProfileID);
      }

      // Save the TetonProfileID for later use:
      datastore[top + "TetonProfileID"] = TetonProfileID;
   }

   areSourceProfilesSet = true;
}

void Teton::setMeshSizeAndPositions()
{
   conduit::Node &datastore = getDatastore();

   int nzones = datastore["size/nzones"].value();
   double *zone_verts_ptr = mMeshBlueprint["fields/zone_verts"].value();
   int *ncorners_ptr = mMeshBlueprint["fields/ncorners"].value();
   int ndim = datastore["size/ndim"].value();
   int maxCorner = datastore["size/maxCorner"].value();

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
   conduit::Node &datastore = getDatastore();

   int nzones = datastore["size/nzones"].value();
   double *velocities_ptr = mMeshBlueprint["fields/velocity_at_corners"].as_double_ptr();
   int *ncorners_ptr = mMeshBlueprint["fields/ncorners"].value();
   int ndim = datastore["size/ndim"].value();
   int maxCorner = datastore["size/maxCorner"].value();

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
   conduit::Node &datastore = getDatastore();

   int nsfaces;
   int *shared_faces_ptr = nullptr;

   if (mMeshBlueprint.has_path("shared_boundaries/nsfaces"))
   {
      nsfaces = mMeshBlueprint["shared_boundaries/nsfaces"].to_int();
      if (nsfaces > 0)
      {
         shared_faces_ptr = mMeshBlueprint["shared_boundaries/shared_faces"].value();
      }
   }
   else // if (datastore.has_path("shared_boundaries/nsfaces"))
   {
      // For backward compatbility
      nsfaces = datastore["shared_boundaries/nsfaces"].value();
      if (nsfaces > 0)
      {
         shared_faces_ptr = datastore["shared_boundaries/shared_faces"].value();
      }
   }
   int ndim = datastore["size/ndim"].value();

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
   conduit::Node &datastore = getDatastore();

   int *connectivity_ptr = mMeshBlueprint["fields/corner_connectivity"].value();
   int *boundaries_types_ptr = datastore["boundary_conditions/type"].value();

   int nzones = datastore["size/nzones"].value();
   int maxCorner = datastore["size/maxCorner"].value();
   int maxcf = datastore["size/maxcf"].value();

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
}

void Teton::setMaterials(double *cornerTe_ptr)
{
   conduit::Node &datastore = getDatastore();

   // TODO add something better than this to check for whether or not we need to move the mesh?
   if (!(mMeshBlueprint.has_path("fields/density")))
   {
      return;
   }

   int nzones = datastore["size/nzones"].value();

   double *density_ptr = mMeshBlueprint["fields/density"].value();
   double *cv_ptr = mMeshBlueprint["fields/heat_capacity"].value();
   double *tez_ptr = mMeshBlueprint["fields/electron_temp"].value();
   double *trz_ptr = mMeshBlueprint["fields/rad_temp"].value();
   double *nez_ptr = mMeshBlueprint["fields/electron_density"].value();

   // Really the effective electron specific energy source.
   if( mMeshBlueprint.has_path("fields/specific_energy_source"))
   {
      double *matSource = mMeshBlueprint["fields/specific_energy_source"].value();
      teton_setmaterialsource(matSource);
   }

   // Initialize arrays to handle multi-material zones
   teton_initmaterial(&cornerTe_ptr[0]);

   for (int zone = 0; zone < nzones; ++zone)
   {
      // TODO: change
      double scm = 0.;
      int zoneID = zone + 1;
      double rho = density_ptr[zone];
      double cv = cv_ptr[zone];
      double tez = tez_ptr[zone];
      double trz = trz_ptr[zone];
      double nez = nez_ptr[zone];
      teton_setmaterial(&zoneID, &cv, &rho, &tez, &trz, &nez, &scm);
   }
   teton_normalizematerial();
}

void Teton::updateOpacity()
{
   conduit::Node &datastore = getDatastore();

   // TODO add something better than this to check for whether or not we need to move the mesh?
   if (!(mMeshBlueprint.has_path("fields/absorption_opacity")))
   {
      return;
   }

   int ngr = datastore["quadrature/num_groups"].to_int(); //coerce from unsigned int or size_t
   int ig;
   std::vector<double> siga_loc;
   std::vector<double> sigs_loc;
   siga_loc.resize(ngr);
   sigs_loc.resize(ngr);

   bool useInternalSigmaS = false;
   if (datastore.has_path("compton/use_internal_sigma_s"))
   {
      useInternalSigmaS = datastore["compton/use_internal_sigma_s"].as_int();
   }
   bool useTableSigmaS = (not useInternalSigmaS);

   // zero out opacities
   teton_initopacity();

   int nzones = datastore["size/nzones"].value();
   double *absorption_opacity_ptr = mMeshBlueprint["fields/absorption_opacity"].value();
   double *scattering_opacity_ptr = nullptr;
   if (useTableSigmaS)
   {
      scattering_opacity_ptr = mMeshBlueprint["fields/scattering_opacity"].value();
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

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
   if (not useTableSigmaS)
   {
      teton_setscatteringopacity(&mInternalComptonFlag);
   }
#endif
}

void Teton::constructIterationControls()
{
   conduit::Node &datastore = getDatastore();

   // These are constructed with default values.
   teton_constructitercontrols();

   // TODO: remove MPI_COMM_WORLD. See
   // https://rzlc.llnl.gov/gitlab/deterministic-transport/TRT/Teton/-/issues/82
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // These are backward-compatible old ways to set the data.
   // TODO: Remove soon (Brunner 10/2021)
   if (datastore.has_path("iteration/noutmx"))
   {
      if (rank == 0)
         std::cout
             << "WARNING: TETON Setting 'iteration/noutmx' in Conduit controls is deprecated, use 'iteration/outerMaxIt' instead.\n";
      int i = datastore["iteration/noutmx"].value();
      teton_adjust_temperature_maxits(&i);
   }
   if (datastore.has_path("iteration/ngdamx"))
   {
      if (rank == 0)
         std::cout
             << "WARNING: TETON Setting 'iteration/ngdamx' in Conduit controls is deprecated, use 'iteration/greyMaxIt' instead.\n";
      int i = datastore["iteration/ngdamx"].value();
      // The old interface that this supplies took sweeps.
      // The new interface takes iterations, so we convert here
      // so that we don't have to explain it in user manuals.
      int greyIters = (i - 1) / 2;
      teton_adjust_grey_maxits(&greyIters);
   }
   if (datastore.has_path("iteration/ndynmx"))
   {
      if (rank == 0)
         std::cout
             << "WARNING: TETON Setting 'iteration/ndynmx' in Conduit controls is deprecated, use 'iteration/incidentFluxMaxIt' instead.\n";
      int i = datastore["iteration/ndynmx"].value();
      teton_adjust_fluxexchange_maxits(&i);
   }

   if (datastore.has_path("iteration/epstmp"))
   {
      if (rank == 0)
         std::cout
             << "WARNING: TETON Setting 'iteration/epstmp' in Conduit controls is deprecated, use 'iteration/outerTempRelTol' instead.\n";
      double x = datastore["iteration/epstmp"].value();
      teton_adjust_temperature_reltol(&x);
   }
   if (datastore.has_path("iteration/epsgda"))
   {
      if (rank == 0)
         std::cout
             << "WARNING: TETON Setting 'iteration/epsgda' in Conduit controls is deprecated, use 'iteration/greyRelTol' instead.\n";
      double x = datastore["iteration/epsgda"].value();
      teton_adjust_grey_reltol(&x);
   }
   if (datastore.has_path("iteration/epsinr"))
   {
      if (rank == 0)
         std::cout
             << "WARNING: TETON Setting 'iteration/epsinr' in Conduit controls is deprecated, use 'iteration/outerPhiRelTol' instead.\n";
      double x = datastore["iteration/epsinr"].value();
      teton_adjust_radenergydensity_reltol(&x);
   }
   if (datastore.has_path("iteration/epsinc"))
   {
      if (rank == 0)
         std::cout
             << "WARNING: TETON Setting 'iteration/epsinc' in Conduit controls is deprecated, use 'iteration/incidentFluxRelTol' instead.\n";
      double x = datastore["iteration/epsinc"].value();
      teton_adjust_fluxexchange_reltol(&x);
   }
}

void Teton::constructDtControls()
{
   conduit::Node &datastore = getDatastore();

   double dtrad = datastore["iteration/dtrad"].value();
   double dtrmn = datastore["iteration/dtrmn"].value();
   double dtrmx = datastore["iteration/dtrmx"].value();
   double delte = datastore["iteration/delte"].value();
   double deltr = datastore["iteration/deltr"].value();

   teton_constructdtcontrols(&dtrad, &dtrmn, &dtrmx, &delte, &deltr);
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
   // Update the C++ copy of the corner temperatures.
   // TODO - Evaluate getting rid of the C++ copy of the corner temperatures and just accessing the field directly
   // from Teton.

   teton_rtedit(mCornerTemperature.data());
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

double *Teton::getCornerTemperature()
{
   return mCornerTemperature.data();
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

void Teton::setTimeStep(int cycle, double dtrad)
{
   double timerad = 0.0;
   double tfloor = 1.0e-5;
   teton_settimestep(&cycle, &dtrad, &timerad, &tfloor);
}

void Teton::setTimeStep(int cycle, double dtrad, double timerad)
{
   double tfloor = 1.0e-5;
   teton_settimestep(&cycle, &dtrad, &timerad, &tfloor);
}

void Teton::updateMeshPositions()
{
   conduit::Node &datastore = getDatastore();

   // TODO add something better than this to check for whether or not we need to move the mesh?
   if (mMeshBlueprint.has_path("fields/zone_verts"))
   {
      int nzones = mZoneToNCorners.size();
      double *zone_verts_ptr = mMeshBlueprint["fields/zone_verts"].value();
      int ndim = datastore["size/ndim"].value();
      int maxCorner = datastore["size/maxCorner"].to_int();

      int off_set = 0;
      std::vector<double> zoneCoordinates(ndim * maxCorner);

      for (int zone = 0; zone < nzones; ++zone)
      {
         int zoneID = zone + 1;
         int ncorners = mZoneToNCorners[zone];
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

      // Update Teton geometry
      teton_getvolume();
   }

   return;
}

// TODO: !!!!! NEED TO CHANGE THIS TO WORK WITH JOE'S NEW BLUEPRINT
//             GENERATE_* FUNCTIONS !!!!!!
// NOTE: the Vectors RadiationForceXTotal, ..., must
//       already be sized to the number of mesh vertices
void Teton::getRadiationForceOnVerts(double *RadiationForceXTotal,
                                     double *RadiationForceYTotal,
                                     double *RadiationForceZTotal)
{
   conduit::Node &datastore = getDatastore();

   std::vector<double> RadiationForce;

   // ??? ASK Tom: why is this here ????
   //updateOpacity();
   //setMeshPositions();
   // ??? ASK Tom: what is this?
   // double volumeRatio = M.getOptions().getGeometryFraction();

   // Compute the radiation force internally in Teton
   // for each zone and corner
   teton_setradiationforce();

   int ndim = datastore["size/ndim"].value();
   int maxCorner = datastore["size/maxCorner"].value();
   RadiationForce.resize(ndim * maxCorner);
   int corner_counter = 0;
   // !!!! NOTE: We assume that the corner ordering is encoded
   //            in corner_to_vertex and that the corner IDs are ordered
   //            in a way that corresponds to looping over zones and
   //            then vertices within the zone
   int nzones = datastore["size/nzones"].value();
   int nverts = datastore["size/nverts"].value();

   for (int v = 0; v < nverts; ++v)
   {
      if (ndim == 3)
      {
         RadiationForceXTotal[v] = 0.0;
         RadiationForceYTotal[v] = 0.0;
         RadiationForceZTotal[v] = 0.0;
      }
      else if (ndim == 2)
      {
         RadiationForceXTotal[v] = 0.0;
         RadiationForceYTotal[v] = 0.0;
      }
   }

   for (int zone = 0; zone < nzones; ++zone)
   {
      int zoneID = zone + 1;
      // Get the radiation force on each corner of each zone
      teton_getradiationforce(&zoneID, &RadiationForce[0]);
      int ncorners = mZoneToNCorners[zone];
      for (int c = 0; c < ncorners; ++c)
      {
         int vertexID = mCornerToVertex[corner_counter];
         corner_counter += 1;
         if (ndim == 3)
         {
            RadiationForceXTotal[vertexID] += RadiationForce[c * 3];
            RadiationForceYTotal[vertexID] += RadiationForce[c * 3 + 1];
            RadiationForceZTotal[vertexID] += RadiationForce[c * 3 + 2];
         }
         else if (ndim == 2)
         {
            RadiationForceXTotal[vertexID] += RadiationForce[c * 2];
            RadiationForceYTotal[vertexID] += RadiationForce[c * 2 + 1];
         }
         // TODO: finish this case
         else if (ndim == 1)
         {
            std::cerr << "1D GetRadiationForceOnVerts not yet implemented! Teton is exiting . . ." << std::endl;
            exit(1);
            // ASK Tom what does mRadialVectorPtr do ????
            //   CCVF &radialVector = *mRadialVectorPtr;

            //   if ( cIDLocal != Geometry::LOCAL_ID_NONE() ) {

            //      if ( ndimReal == 3 )
            //      {
            //         RadiationForce1D = 0.25*volumeRatio*RadiationForce[cIDLocal];

            //         cornerRadForce[*cidP].SetX( RadiationForce1D*radialVector[*cidP].GetX() );
            //         cornerRadForce[*cidP].SetY( RadiationForce1D*radialVector[*cidP].GetY() );
            //         cornerRadForce[*cidP].SetZ( RadiationForce1D*radialVector[*cidP].GetZ() );
            //      }
            //      else
            //      {
            //         RadiationForce1D = 0.5*volumeRatio*RadiationForce[cIDLocal];

            //         cornerRadForce[*cidP].SetX( RadiationForce1D*radialVector[*cidP].GetX() );
            //         cornerRadForce[*cidP].SetZ( RadiationForce1D*radialVector[*cidP].GetZ() );
            //      }
            //   }
            //   else
            //   {
            //      RadiationForce1D = 0.0;
            //      if ( ndimReal == 3 ) {
            //         cornerRadForce[*cidP].SetX( RadiationForce1D );
            //         cornerRadForce[*cidP].SetY( RadiationForce1D );
            //         cornerRadForce[*cidP].SetZ( RadiationForce1D );
            //      } else {
            //         cornerRadForce[*cidP].SetX( RadiationForce1D );
            //         cornerRadForce[*cidP].SetZ( RadiationForce1D );
            //      }

            //   }
         }
      }
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
