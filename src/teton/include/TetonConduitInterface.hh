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
#include <map>
#include <mpi.h>
#include <string>

namespace Teton
{
class Teton
{
   // ---------------------------------------------------------------------------
   // Some internal constants.
   // ---------------------------------------------------------------------------
   static const std::string PREFIX;
   static const std::string MCARRAY_PREFIX;
   static const std::string PARTITION_FIELD;
   static const std::string PARTITION_FIELD_BOUNDARY;

  public:
   Teton();

   ~Teton();

   conduit::Node &getMeshBlueprint()
   {
      return getDatastore()["blueprint"];
   }
   const conduit::Node &getMeshBlueprint() const
   {
      return getDatastore()["blueprint"];
   }

   conduit::Node &getDatastore();
   const conduit::Node &getDatastore() const;

   conduit::Node &getOptions()
   {
      return getDatastore()["options"];
   }
   const conduit::Node &getOptions() const
   {
      return getDatastore()["options"];
   }

   /*!
    * \brief Returns a node that contains the partitioned mesh. If partitioning
    *        is not being done then the blueprint node is returned.
    *
    * \return A Conduit node that contains the partitioned mesh, or the blueprint
    *         mesh if no partitioning is being done.
    */
   conduit::Node &getMeshBlueprintPart();
   const conduit::Node &getMeshBlueprintPart() const;

   /*!
    * \brief Get whether verbose output is selected.
    *
    * \return An integer that indicates the verbosity level.
    *
    * \note The value is retrieved from the options, if the value exists or from
    *       the environment. The environment can override the options.
    */
   int getVerbose() const;

   void initialize(MPI_Comm communicator, bool fromRestart = false);

   // Advance a radiation step, returns dt recommended by Teton for the next time step
   double step(int cycle);

   // Requires the mesh blueprint node to have 'state/cycle' populated.
   void dump(MPI_Comm communicator, std::string path = ".");

   void setTimeStep(int cycle, double dtrad, double timerad);

   /*!
    * \brief Update the mesh positions for Teton, possibly repartitioning.
    */
   void updateMeshPositions();

   //------------------------------------------------------------------------
   // Result-getting functions
   //------------------------------------------------------------------------

   /*!
    \brief Get the radiation temperature for a zone. 

    \param zone A 1-origin zone id valid for the blueprint mesh.

    \note This value is obtained from the blueprint mesh's radiation_energy_density
          field, which must exist in order to return a valid value.

    \return The radiation temperature.
    */
   double getRadiationTemperature(int zone) const;

   /*!
    \brief Get the radiation deposited for a zone. 

    \param zone A 1-origin zone id valid for the blueprint mesh.

    \note This value is obtained from the blueprint mesh's electron_energy_deposited
          field, which must exist in order to return a valid value.

    \return The radiation deposited.
    */
   double getRadiationDeposited(int zone) const;

   /*!
    \brief Get the material temperature for a zone.

    \param zone A 1-origin zone id valid for the blueprint mesh.

    \return The material temperature.
    */
   double getMaterialTemperature(int zone) const;

   // TODO: remove these once all host codes swich to getting the
   // force density fields from the conduit node
   void getRadiationForceDensity1D(double *RadiationForceDensityX);

   /*!
    \brief This method copies the radiation_force_{x,y,x} or {z,r} fields into
           the supplied data arrays. The arrays should be sized the same as
           the radiation_force_z field.

    \param[out] RadiationForceDensityX The array to hold x (or z data in 2D).
    \param[out] RadiationForceDensityY The array to hold y (or r data in 2D).
    \param[out] RadiationForceDensityZ The array to hold z data.

    \note This method simply returns the radiation force components in the
          supplied output arrays. This same thing can be achieved using the
          radiation_force_* fields from blueprint directly.

    \note remove once host codes swich to getting the data from Conduit node.
    */
   void getRadiationForceDensity(double *RadiationForceDensityX,
                                 double *RadiationForceDensityY,
                                 double *RadiationForceDensityZ);

   /*!
     \brief This is used for the post-ALE step of rescaling psi
            based on the remapped radiation energy density.
            Here the array rad_energy_density needs to be sized
            to (ngroups * nzones) before being passed.

     \note The underlying teton_reconstructpsi function does bookkeeping
           that involves sets and volume of the geometry it knows about.
           It also sets Mat%trz for each zone. This says that this function
           should be called on the partitioned mesh.

           However, the host code will be passing a rad_energy_density that
           is valid for the blueprint mesh and not the partitiond mesh. This
           suggests another partition step internally to send the field to
           the partition mesh.

           The returned rad_energy is valid for the partitioned mesh. That
           computation can't be done on the blueprint mesh without knowing
           the zone volume. It might not matter though since the host code
           appears to reduce(sum) over all ranks.

     \param[out] rad_energy        Total rad energy.
     \param[in] rad_energy_density Array sized double[ngroups][nzones] that contains
                                   the radiation energy density (nzones varies fastest).
    */
   void reconstructPsi(double *rad_energy, const double *rad_energy_density);

   /*!
    \brief This is used to update the angular intensity to be consistent with
           changes in the corner volumes of the mesh from the Lagrange motion.
           That is, psi is rescale so that the total radiation energy in the
           zone remains constant.
    */
   void reconstructPsiFromdV();

   /*!
    \brief Copies the zonal psi values into the \a psi array. The psi array
           should be sized double[nAngles][ngroups][nzones] (nzones varying fastest).
           If the Teton operated on a partitiond mesh, this involves mapping the
           partitioned psi back onto the blueprint mesh.

    \param[in] numAngles The number of angles. This must match what Teton was given previously.
    \param[out] psi The zonal psi array that will contain the values.

    \note Ideally we would not pass numAngles as Teton should already know it.
    */
   void getZonalPsi(int numAngles, double *psi);

   /*!
    * \brief Tell Teton to compute the radiation flux and make it available to retrieve
    *        via getRadiationFlux.
    *
    * \note Makes a bunch of Conduit fields to store the different groups of the radiation
    *       flux and migrates these fields back to the original mesh where they are left
    *       as separate scalar fields, rather than being reassembled.
    */
   void setRadiationFlux();

   /*!
    * \brief Reads values from Teton into the supplied array. The values are obtained
    *       from teton_getradiationflux.
    *
    * \param zone A 1-origin zone index.
    * \param[out] zflux Array sized double[ngroups][ndims]
    *
    * \note The Conduit scalars that were produced in setRadiationFlux are sampled for
    *       the specified zone and data are stored into zflux.
    */
   void getRadiationFlux(int zone, double *zflux) const;

   /*!
    * \brief Get information about the Teton run.
    *
    * \param[out] noutrt,          number of thermal [outer] iterations in this cycle
    * \param[out] ninrt            number of inner [transport] iterations in this cycle
    * \param[out] ngdart           number of grey diffusion iterations
    * \param[out] nNLIters         number of nonlinear iterations
    * \param[out] maxNLIters       maximum nonlinear iterations
    * \param[out] TrMaxZone        zone id with maximum T rad
    * \param[out] TeMaxZone        zone id with maximum electron temperature
    * \param[out] TrMaxProcess     MPI process with maximum T rad
    * \param[out] TeMaxProcess     MPI process with maximum electron temperature
    * \param[out] dtused           Teton returns the time step used in cycle just completed
    * \param[out] dtrad            radiation vote for next time step
    * \param[out] TrMax            T rad maximum value
    * \param[out] TeMax            T electron maximum value
    * \param[out] EnergyRadiation  energy contained in the radiation field
    * \param[out] PowerIncident    power of photons incident
    * \param[out] PowerEscape      power of photons escaping
    * \param[out] PowerAbsorbed    power of energy absorbed
    * \param[out] PowerEmitted     power of photons emitted
    * \param[out] PowerExtSources  power of photons from fixed volumetric sources
    * \param[out] PowerCompton     power of Compton scattering photons
    * \param[out] EnergyCheck      energy not accounted for this cycle.
    */
   void getEdits(int &noutrt,
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
                 double &EnergyCheck) const;

   /*!
    * \brief Get iteration/dtcontrol values, taking into account possible partitioning.
    *
    * \param[out] flag A reason indicating a change in time step.
    * \param[out] process The MPI rank of the process that is causing time step change.
    * \param[out] zone The 1-origin zone number that is causing time step change.
    * \param[out] message An informative message for the time step change.
    */
   void getDtControls(int &flag, int &process, int &zone, std::string &message) const;

   void dumpTallyToJson() const;

   void setSourceProfiles();

   // Do not call this unless setSourceProfiles has been set:
   void resetSourceProfiles();

   // NOTE: These functions return information for the partitioned corner mesh.
   //       Codes that would be calling these expect to deal with the blueprint
   //       mesh. Can we eliminate the reason for why a code would need to call
   //       these? These methods are not compatible with partitioning.
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

   // ---------------------------------------------------------------------------
   // Functions pertaining to checkpoints/restarts
   // ---------------------------------------------------------------------------
   conduit::Node &getCheckpoint();
   void checkpointPrepareForLoad();
   void checkpointPrepareForSave();
   void checkpointDataLoaded();
   void checkpointExternalDataLoaded();
   void checkpointFinished();

   // ---------------------------------------------------------------------------
   // Some relevant field and topo names.
   // ---------------------------------------------------------------------------
   static const std::string FIELD_ELECTRON_ENERGY_DEPOSITED;
   static const std::string FIELD_RADIATION_ENERGY_DENSITY;
   static const std::string FIELD_RADIATION_TEMPERATURE;
   static const std::string FIELD_RADIATION_FORCE_X;
   static const std::string FIELD_RADIATION_FORCE_Y;
   static const std::string FIELD_RADIATION_FORCE_Z;
   static const std::string FIELD_RADIATION_FORCE_R;
   static const std::string FIELD_CORNER_VOLUME_SUMS;
   static const std::string FIELD_RADIATION_FLUX_X;
   static const std::string FIELD_RADIATION_FLUX_Y;
   static const std::string FIELD_RADIATION_FLUX_Z;
   static const std::string FIELD_RADIATION_FLUX_R;
   static const std::string FIELD_MATERIAL_TEMPERATURE;
   static const std::string TOPO_MAIN;
   static const std::string TOPO_BOUNDARY;

  private:
   // ---------------------------------------------------------------------------
   // Internal helpers
   // ---------------------------------------------------------------------------

   /*!
    * \brief Get the main topology.

    * \param root The node through which we'll get the main topology.

    * \return A reference to the main topology.
    */
   conduit::Node &getMainTopology(conduit::Node &root);
   const conduit::Node &getMainTopology(const conduit::Node &root) const;

   /*!
    \brief Creates a new zonal field it does not exist.

    \param root The root node into which we'll create the field.
    \param topoName The name of the topology.
    \param fieldName The name of the field to make.
    \param nzones The number of zones.
    */
   void createZonalField(conduit::Node &root, const std::string &topoName, const std::string &fieldName, int nzones);

   /*!
    \brief Reads values from Teton into the supplied array. The values are obtained
           from teton_getradiationdeposited.

    \param[out] RadEnergyDeposited The destination array that holds nzones values.
    \param nzones The number of zones.
    */
   void getRadEnergyDeposited(double *RadEnergyDeposited, int nzones) const;

   /*!
    \brief Reads values from Teton into the supplied array. The values are obtained
           from teton_getradiationtemperature.

    \param[out] RadTemp The destination array that holds nzones values.
    \param nzones The number of zones.
    */
   void getRadiationTemperature(double *RadTemp, int nzones) const;

   /*!
    \brief Reads values from Teton into the supplied array. The values are obtained
           from teton_getmaterialtemperature.

    \param[out] MatTemp The destination array that holds nzones values.
    \param nzones The number of zones.
    */
   void getMaterialTemperature(double *MatTemp, int nzones) const;

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

   /*!
    * \brief Constructs surface flux tallies. Results are output into fields that
    *        are present on the tally surfaces in the blueprint node (even when 
    *        partitioning is enabled).
    */
   void computeGenericSurfaceFluxTally();

   void constructSize();

   void constructMemoryAllocator();

   void constructQuadrature();

   // set Teton node positions
   void setMeshSizeAndPositions();

   // set corner velocities
   void setMeshVelocity();

   // communication
   void setCommunication();

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

   void SumSharedNodalValues(conduit::Node &root, double *nodal_field);

   // updates the radiation force if the fields
   // "radiation_force_r" (dim == 2) or "radiation_force_x"
   // (dim == 3) fields are present in the conduit blueprint node
   void updateRadiationForce();
   void updateZonalRadiationForce();

   double *getCornerTemperature();

   void updateRadEnergyDeposited();

   // ---------------------------------------------------------------------------
   // Internal functions pertaining to partitioning
   // ---------------------------------------------------------------------------

   /*!
    * \brief Consult the options and return whether partitioning is being done.
    *
    * \return True of partitioning is enabled; False otherwise.
    */
   bool doPartitioning() const;

   /*!
    * \brief Creates a mapping between boundary elements and main topology elements
    *        and then uses the partition field on the main topology to generate new
    *        partition fields for the supplied topology names.
    *
    * \param mesh The input mesh.
    * \param topoNames A vector of topology names.
    *
    * \return A vector containing new partition field names.
    */
   std::vector<std::string> createPartitionFields(conduit::Node &mesh, const std::vector<std::string> &topoNames);

   /*!
    * \brief Take a partitioned secondary mesh and combine it into the partitioned
    *        mesh node, rewriting its connectivity so it shares coordinates with
    *        the partmesh instead of having its own coordset.
    *
    * \param partmesh The node that contains the partitioned volume mesh.
    * \param topoName The name of the volume topology.
    * \param secondPartmesh The partitioned secondary mesh.
    * \param secondTopoName The name of the secondary topology.
    *
    * \note We assume that the number of domains in partmesh and bpartmesh is
    *       the same and the domains ids are in the same order.
    */
   void assimilateTopology(conduit::Node &partmesh,
                           const std::string &topoName,
                           conduit::Node &secondPartmesh,
                           const std::string &secondTopoName);

   /*!
    * \brief Returns whether the field (assumed to be size ngroup*nzones) needs
    *        interleaving or not to be represented as a Blueprint mcarray. An interleaved
    *        field would be sized array[nzones][ngroup] whereas an array sized
    *        array[ngroup][nzones] would not need interleaving.
    *
    * \param fieldName The name of the field we're checking.
    *
    * \return True if the field needs interleaving, false otherwise.
    */
   bool doInterleave(const std::string &fieldName) const;

   /*!
    * \brief Scan through fields in the source mesh and identify those
    *        which have too many elements (multiple of numGroups). Create
    *        alternate mcarray zero-copy representations for those fields so
    *        they can pass through the partitioner.
    *
    * \param root The root, which is usually the blueprint mesh. When mapping back,
    *             we can pass the part mesh.
    */
   void add_mcarray_fields(conduit::Node &root);

   /*!
    * \brief Remove any mcarray fields from the root node.
    *
    * \param root The root, which is usually the blueprint mesh. When mapping back,
    *             we can pass the part mesh.
    */
   void remove_mcarray_fields(conduit::Node &root);

   /*!
    * \brief Fetch a conduit mcarray field node given the original name of the field.
    *        If \a fieldName is not an mcarray created by the Teton interface, return
    *        the regular field node. If the node is not located, a Conduit exception
    *        is thrown.
    *
    * \param root A top-level node like the blueprint or partition node that contains
    *             a "fields" node.
    * \param fieldName The name of the field that may have been converted to mcarray.

    * \return The field node that contains the mcarray data.
    */
   const conduit::Node &fetch_mcarray(const conduit::Node &root, const std::string &fieldName) const;

   /*!
    * \brief Takes the mesh from the blueprint node and partitions it and stores
    *        the partitiond mesh in the blueprint_partitioned node.
    *
    * In our scheme, the host code will use the getMeshBlueprint() method to store
    * the mesh on which Teton will operate. That node will optionally go through
    * the partition() method to repartition the mesh. Various methods from then on
    * will access getMeshBlueprint() when they want to access the original mesh
    * from the host code. The getMeshBlueprintPart() method is used when the
    * partitioned mesh is desired.
    *
    * \param fromRestart Whether we're initializing from a restart, in which case,
    *                    the partitioning should be redone.
    */
   void partition(bool fromRestart);

   /*!
    * \brief Return the names of the topologies that share a coordset with main.
    *        This will include the boundary topo and possibly some surface tally
    *        topos.
    *
    * \param mesh The top-level node that contains the "topologies" node.
    * \return A vector containing topology names that share a coordset with main.
    */
   std::vector<std::string> getPartitionTopologies(const conduit::Node &mesh) const;

   /*!
    * \brief Propagate fields from the blueprint mesh supplied by the host to the
    *        partitioned mesh.
    *
    * \param topoName   The topology that owns the fields.
    * \param fieldNames The names of the fields to repartition and put on the 
    *                   partitiond mesh.
    * \param updateCoords Whether to update the coordinates on the partitioned mesh.
    */
   void sendFieldsOrig2Part(const std::string &topoName, const std::vector<std::string> &fieldNames, bool updateCoords);

   /*!
    * \brief Propagate fields from the partitioned mesh back to the mesh
    *        supplied by the host to the partitioned mesh.
    *
    * \param topoName   The topology that owns the fields.
    * \param fieldNames The names of the fields to repartition and put on the 
    *                   partitiond mesh.
    */
   void sendFieldsPart2Orig(const std::string &topoName, const std::vector<std::string> &fieldNames);

   /*!
    * \brief Optionally repartition the coordinates and then update the mesh
    *        positions for Teton from the part mesh.
    *
    * \param doPartition True if we want to do partitioning (if enabled)
    */
   void updateMeshPositions(bool doPartition);

   /*!
    * \brief Clean up for partitioning. This removes extra fields that were added, etc.
    */
   void partitionCleanup();

   /*!
    * \brief Return the names of the radiation force density fields, taking into
    *        account the mesh dimension.
    * 
    * \return A vector of field names.
    */
   const std::vector<std::string> &radiationForceDensityFields() const;

   /*!
    * \brief Initialize the radiation force density field names. We store them in mRadiationForceDensityFields.
    * 
    */
   void initializeRadiationForceDensityFieldNames();

   /*!
    * \brief Return a vector of double pointers that correspond do the radiation force density fields.
    *
    * \param root The root node in which to search for the fields.
    *
    * \return A vector of double pointers for the fields.
    */
   std::vector<double *> radiationForceDensity(conduit::Node &root) const;

   /*!
    * \brief Create the radiation force density fields on the blueprint mesh if they
    *        do not exist. The fields need to be there in order for getRadiationForceDensity
    *        to work.
    *
    * \param root The root node where the fields will be created if they do not exist.
    * \param elementAssociation Pass true to make element-associated fields or false to make
    *                          veretx-associated fields.
    */
   void createRadiationForceDensity(conduit::Node &root, bool association);

   /*!
    * \brief Create the return radiation temperature field on the blueprint mesh,
    *        if it does not exist.
    */
   void createRadiationTemperature();

   /*!
    * \brief Create the return material temperature field on the blueprint mesh,
    *        if it does not exist.
    */
   void createMaterialTemperature();

   /*!
    * \brief Initialize the radiation flux field names. We store them in mRadiationFluxFields.
    */
   void initializeRadiationFluxFieldNames();

   /*!
    * \brief Return the radiation flux field names for the mesh dimension.
    *
    * \return A vector of field names that contain the radiation flux.
    */
   const std::vector<std::string> &getRadiationFluxFields() const;
   /*!
    * \brief Given a rank and zone id, for the original Teton mesh, return the
    *        rank and zone id in the partitioned mesh, if it exists. If the provided
    *        domain and rank do not exist then return {-1,-1} in the output domain, zone.
    *
    * \param originalDomZone Contains the domain and zone number we're looking for
    *                        in the partitioned mesh. Component [0] contains a 0-origin
    *                        MPI rank number. Component [1] contains a 1-origin Teton
    *                        zone id.
    * \param[out] partDomZone The rank and zone id of the zone in the partitioned mesh.
    *                        Component [0] contains a 0-origin MPI rank number.
    *                        Component [1] contains a 1-origin Teton zone id.
    */
   void zoneLookupOrig2Part(int originalDomZone[2], int partDomZone[2]) const;

   /// Flags that are useful for the test() method.
   enum
   {
      Test_RadiationForceDensity = 1,
      Test_ZonalPsi = 2,
      Test_MaterialTemperature = 4,
      Test_RadiationTemperature = 8,
      Test_RadiationDeposited = 16,
      Test_ReconstructPsi = 32
   };

   /*!
    * \brief Capture the current Blueprint state so we can use it to test
    *        against a baseline.
    *
    * \param[out] n    The node that holds the output node.
    * \param datastore The Teton datastore.
    * \param bp        The Conduit node that contains the Teton mesh.
    * \param options   The Conduit node that contains the Teton options.
    * \param flags     A set of flags (or-ed together) to activate various tests.
    *
    * \return The name of the file to use for current or baseline.
    */
   std::string makeTestNode(conduit::Node &n,
                            conduit::Node &datastore,
                            conduit::Node &bp,
                            conduit::Node &options,
                            int flags);

  private:
   double mDTrad;

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

   std::vector<std::string> mMapBackFields;               //!< Vector of field names to map back during partitioning.
   std::map<std::string, std::string> mMCArrays;          //!< Map of field names to mcarray names
   std::vector<std::string> mRadiationForceDensityFields; //!< Vector of field names for radiation force density
   std::vector<std::string> mRadiationFluxFields;         //!< Vector of field names for radiation flux
};
} //end namespace Teton

#endif // __TETON_CONDUIT_INTERFACE_HH__
