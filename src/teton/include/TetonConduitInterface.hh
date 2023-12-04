//--------------------------------------------------------------------------//
// TetonConduitInterface.hh
//
// This class provides a C++ object oriented interface for setting up Teton
// that relies on input passed in via a conduit node.
//
// This interface is still under development.  If you are interested in
// setting up  Teton via a conduit node contact the Teton team for more info.
//
// Using this interface requires compiling Teton with Conduit support.
//--------------------------------------------------------------------------//

#ifndef __TETON_CONDUIT_INTERFACE_HH__
#define __TETON_CONDUIT_INTERFACE_HH__

#include "TetonSources.hh"
#include "conduit/conduit.hpp"
#include <mpi.h>
#include <string>

namespace Teton
{
class Teton
{
  public:
   Teton() : areSourceProfilesSet(false), mIsInitialized(false)
   {
   }

   // Teton destructor.  Important note:  This function requires MPI to still be
   // active ( not finalized ), as it uses MPI when tearing down Teton.
   ~Teton();

   void initialize(MPI_Comm communicator, bool fromRestart = false);

   // This stores the needed mesh data needed to compute
   // forces on the vertices
   void storeMeshData();

   // sanitizer_node must have `level`:
   //   0 - no sanitizer (does nothing and returns
   //   1 - quieter sanitizer
   //   2 - noisy sanitizer
   // Optional entries:
   //   -`cat_list` (defaults to all inputs except \sigma_s)
   //   -`kill_if_bad` (defaults to false)
   // Returns the number of bad input categories
   int checkInputSanity(const conduit::Node &sanitizer_node) const;

   void constructBoundaries();

   void constructComptonControl();

   void constructEdits();

   void computeGenericSurfaceFluxTally();
   void dumpTallyToJson() const;

   void constructSize();

   void constructMemoryAllocator();

   // Advance a radiation step, returns dt recommended by Teton for the next time step
   double step(int cycle);

   void constructQuadrature();

   // set Teton node positions
   void setMeshSizeAndPositions();

   // set corner velocities
   void setMeshVelocity();

   // communication
   void setCommunication();

   void setSourceProfiles();
   // Do not call this unless setSourceProfiles has been set:
   void resetSourceProfiles();

   // set up the mesh connectivity arrays used in Teton.
   void setMeshConnectivity();

   // Initializes material properties used by Teton at the beginning
   // of a run or after a restart.
   void setMaterials();

   // create IterControl/DtControl object used in Teton for
   // iteration/time-step control etc.
   void constructIterationControls();
   void constructDtControls();

   // a member function to update opacity
   void updateOpacity();

   // updates the radiation force if the fields
   // "radiation_force_r" (dim == 2) or "radiation_force_x"
   // (dim == 3) fields are present in the conduit blueprint node
   void updateRadiationForce();
   void updateZonalRadiationForce();

   // Requires the mesh blueprint node to have 'state/cycle' populated.
   void dump(MPI_Comm communicator, std::string path = ".");

   double *getCornerTemperature();

   double getMaterialTemperature(int zone);

   double getRadiationTemperature(int zone);

   double getRadiationDeposited(int zone);

   void setTimeStep(int cycle, double dtrad, double timerad);
   // This updates Teton's zone vertex coordinates based on
   // changes to the mesh nodes (from hydro)
   void updateMeshPositions();

   // TODO: remove these once all host codes swich to getting the
   // force density fields from the conduit node
   void getRadiationForceDensity1D(double *RadiationForceDensityX);
   void getRadiationForceDensity(double *RadiationForceDensityX,
                                 double *RadiationForceDensityY,
                                 double *RadiationForceDensityZ);

   // Updates the field "fields/rad_energy_deposited/values" in the
   // blueprint conduit node if it exists
   void getRadEnergyDeposited(double *RadEnergyDeposited);
   void updateRadEnergyDeposited();

   int *getCornerToVertexArray()
   {
      return &mCornerToVertex[0];
   }
   int *getZoneToNCornersArray()
   {
      return &mZoneToNCorners[0];
   }
   int *getZoneToCornersArray()
   {
      return &mZoneToCorners[0];
   }
   // This is used for the post-ALE step of rescaling psi
   // based on the remapped radiation energy density.
   // Here the array rad_energy_density needs to be sized
   // to (ngrousps * nzones) before being passed
   void reconstructPsi(double *rad_energy, double *rad_energy_density);
   // This is used to update the angular intensity
   // to be consistent with changes in the corner volumes
   // of the mesh from the Lagrange motion. That is,
   // psi is rescale so that the total radiation energy in
   // the zone remains constant.
   void reconstructPsiFromdV();

   conduit::Node &getMeshBlueprint()
   {
      return getDatastore()["blueprint"];
   }
   conduit::Node &getDatastore();
   conduit::Node &getOptions()
   {
      return getDatastore()["options"];
   }
   // Const version of above functions:
   const conduit::Node &getMeshBlueprint() const
   {
      return getDatastore()["blueprint"];
   }
   const conduit::Node &getDatastore() const;
   const conduit::Node &getOptions() const
   {
      return getDatastore()["options"];
   }

   // ---------------------------------------------------------------------------
   // Functions pertaining to checkpoints/restarts
   // ---------------------------------------------------------------------------
   conduit::Node &getCheckpoint();
   void checkpointPrepareForLoad();
   void checkpointPrepareForSave();
   void checkpointDataLoaded();
   void checkpointExternalDataLoaded();
   void checkpointFinished();

   double mDTrad;

  private:
   bool areSourceProfilesSet; // Whether or not setSourceProfiles has been called
   bool mIsInitialized;       // Whether or not Teton::initialize has been called

   int mGTAorder; // quadrature order used for grey transport acceleration (def=2 for s2 acc)
   int mInternalComptonFlag;

   // Cached MPI communicator details:
   MPI_Comm mCommunicator;
   int mRank;

   // list of sources to be appended to right hand side before each time step
   // these could be point sources, MMS, etc.
   TetonSourceManager mSourceManager;

   // To compute radiation forces on the vertices, Teton
   // needs to hang on to this connectivity array
   // !!!! NOTE: THIS WILL NOT WORK FOR AMR MESHES, WHERE A
   //      CORNER CAN CORRESPOND TO A HANGING NODE !!!!
   // !!!! TODO: Make this work for AMR meshes
   std::vector<int> mCornerToVertex;
   std::vector<int> mZoneToNCorners;
   std::vector<int> mZoneToCorners;
   std::vector<int> mCornerToZone;
};
} //end namespace Teton

#endif // __TETON_CONDUIT_INTERFACE_HH__
