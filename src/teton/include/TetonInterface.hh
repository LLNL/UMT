//--------------------------------------------------------------------------//
// TetonInterface.hh
//
//  TetonInterface  - A set of C/C++ headers to enable
//             host codes to interact with and use the Teton Sn Library
//--------------------------------------------------------------------------//

#ifndef __TETON_INTERFACE_HH__
#define __TETON_INTERFACE_HH__

#include "mpi.h"
#include <iostream>
#include <string>
#include <vector>

#include "TetonVersion.hh"

// Forward declaration so that it can be used in the API definitions:
enum class tetonComptonFlag;

// C versions of the Teton library interface
extern "C"
{
// adjust Teton's outer zone average temperature iteration parameters
void teton_adjusttemperaturegoals(
   const int *maxIters,      // maximum number of outer temperature iterations per time step
   const double *tolerance); // pointwise relative tolerance for zone average electron temperature

// adjust Teton's convergence criteria for radiation energy density
void teton_adjustradenergydensitygoal(
   const double *tolerance); // pointwise relative tolerance for radiation energy density iteration

// adjust grey acceleration sweeping parameters.  There are 2 sweeps per CG iteration
void teton_adjustgreygoals(
   const int *maxIters,      // maximum number of grey acceleration sweeps to perform in each grey acceleration
   const double *tolerance); // pointwise relative tolerance for zone average electron temperature

// adjust Teton's boundary exchange convergence criteria during the linear solve step of each Temperature iteration in a time step
void teton_adjustfluxexchangegoals(
   const int *maxIters, // maximum number of sweep halo exchanges allowed to MPI process boundary incoming fluxes
   const double *
      tolerance); // pointwise relative tolerance on incident angular intensities.  Default internal Teton value is 4e-5

// adjust Teton's nonlinear solver parameters
void teton_adjustnonlinearsolvegoals(
   const int *maxIters,      // maximum number of nonlinear iterations per outer temperature iteration
   const double *tolerance); // pointwise relative tolerance on nonlinear iteration

// Same function as above, but you can set only one of the values, leaving the
// other alone.
void teton_adjust_fluxexchange_reltol(const double *RelTol);
void teton_adjust_fluxexchange_maxits(const int *MaxIts);
void teton_adjust_grey_reltol(const double *RelTol);
void teton_adjust_grey_maxits(const int *MaxIts);
void teton_adjust_nonlinear_reltol(const double *RelTol);
void teton_adjust_nonlinear_maxits(const int *MaxIts);
void teton_adjust_radenergydensity_reltol(const double *RelTol);
// Note: There is no rad energy density max its.  It's the same
// loop as the outer temperature iterations, so we only need one.
void teton_adjust_temperature_reltol(const double *RelTol);
void teton_adjust_temperature_maxits(const int *MaxIts);

void teton_get_fluxexchange_reltol_internal(const double *RelTol);
inline double teton_get_fluxexchange_reltol()
{
   double v{0};
   teton_get_fluxexchange_reltol_internal(&v);
   return v;
}
void teton_get_fluxexchange_maxits_internal(const int *MaxIts);
inline int teton_get_fluxexchange_maxits()
{
   int v{0};
   teton_get_fluxexchange_maxits_internal(&v);
   return v;
}
void teton_get_grey_reltol_internal(const double *RelTol);
inline double teton_get_grey_reltol()
{
   double v{0};
   teton_get_grey_reltol_internal(&v);
   return v;
}
void teton_get_grey_maxits_internal(const int *MaxIts);
inline int teton_get_grey_maxits()
{
   int v{0};
   teton_get_grey_maxits_internal(&v);
   return v;
}
void teton_get_nonlinear_reltol_internal(const double *RelTol);
inline double teton_get_nonlinear_reltol()
{
   double v{0};
   teton_get_nonlinear_reltol_internal(&v);
   return v;
}
void teton_get_nonlinear_maxits_internal(const int *MaxIts);
inline int teton_get_nonlinear_maxits()
{
   int v{0};
   teton_get_nonlinear_maxits_internal(&v);
   return v;
}
void teton_get_radenergydensity_reltol_internal(const double *RelTol);
inline double teton_get_radenergydensity_reltol()
{
   double v{0};
   teton_get_radenergydensity_reltol_internal(&v);
   return v;
}
void teton_get_temperature_reltol_internal(const double *RelTol);
inline double teton_get_temperature_reltol()
{
   double v{0};
   teton_get_temperature_reltol_internal(&v);
   return v;
}
void teton_get_temperature_maxits_internal(const int *MaxIts);
inline int teton_get_temperature_maxits()
{
   int v{0};
   teton_get_temperature_maxits_internal(&v);
   return v;
}

// Internal use only.
void teton_get_default_outer_temp_reltol_internal(const double *var);
void teton_get_default_outer_intensity_reltol_internal(const double *var);
void teton_get_default_grey_reltol_internal(const double *var);
void teton_get_default_incident_flux_reltol_internal(const double *var);
void teton_get_default_inner_nl_reltol_internal(const double *var);
void teton_get_default_outer_max_it_internal(const int *var);
void teton_get_default_grey_max_it_internal(const int *var);
void teton_get_default_incident_flux_max_it_internal(const int *var);
void teton_get_default_inner_nl_max_it_internal(const int *var);

// add a Teton boundary to the Teton boundary list
void teton_addboundary(
   const int *numBcTotal, // scalar, number of bc specified in this call
   const int *bcType,     // array(numBcTotal), of BC type tags integers mapped to the tetonBcTypeFlag enum
   const int *numBdyElem, // array(numBcTotal), specifying number of elements that have a given BC
   const int
      *neighborID); // array(numBcTotal), specifying neighboring process ID's (values used if shared boundary only)

//
//  teton_addprofile_internal
//
//  Don't call this, call one of the teton_addprofile APIs below!
void teton_addprofile_internal(const int *NumTimes,
                               const int *NumValues,
                               const double *Multiplier,
                               const bool *BlackBody,
                               const bool *Isotropic,
                               const double *Times,
                               const double *Values,
                               int *TetonProfileID);

// append source terms to the right hand side for the next time step
//
// the units of sourceValues are such that sourceValues[angle][zone][group]*c*dt has the same units as psi
//   group varies fastest, angle varies slowest.
//
// if numSourceAngles == 1, then all source values will be divided by wtiso (4pi)
// if numSourceAngles == numPolarAngles, then all source values will be divided by half of wtiso (2pi)
//
// Future TODO: variants of this function for isotropic sources?
void teton_appendsourcetopsi(
   const int *numSourceZones,   // scalar
   const int *numSourceAngles,  // scalar, must be 1, # of polar angles, or # of angles
   const int *numSourceGroups,  // scalar, must match # of groups in Teton
   const int *sourceZoneList,   // array of length numSourceZone
   const double *sourceValues); // numAngles x numSourceZones x numGroups, same units as \psi/c/dt

//
//  teton_resetprofile_internal
//
//  Don't call this, call one of the teton_resetprofile APIs below!
void teton_resetprofile_internal(const int *TetonProfileID,
                                 const int *NumTimes,
                                 const int *NumValues,
                                 const double *Multiplier,
                                 const double *Times,
                                 const double *Values);

// verify shared subdomain boundaries with MPI message passing testing
void teton_checksharedboundary();

//
//  teton_checkinputsanity
//
//  Have Teton check all host provided input for sanity (e.g., physicality)
//
//  Notes: verbosityLevel
//           verbose=x0 - No complaints written to stdout
//           verbose=01 - rank 0 complains about the category of bad data
//           verbose=02 - rank 0 complains about individual pieces of bad data
//             (Will print zone, group, and corner indices as appropriate, and the offending value)
//           verbose=11 - all ranks complain about the category of bad data
//           verbose=12 - all ranks complain about individual pieces of bad data
//             (Will print zone, group, and corner indices as appropriate, and the offending value)
//
//          As of 4/25/2018 the following categories can be checked (numbers beside each should be used in listOfCategoriesToCheck
//             1 - scattering cross section (sigs)
//             2 - absorption cross section (siga)
//             3 - effective specific heat  (cve)
//             4 - density                  (rho)
//             5 - Zone avg electron temp   (tez)
//             6 - Electron density         (nez)
//             7 - corner volumes           (Volumes)
//
//  10/7/2020 - sigs and siga checking will now also check for unphysically large opacities (> 1.e40)
//
void teton_checkinputsanity(
   const bool *killIfBad,              // Logical- Call MPI_ABORT on bad input
   const int *verbosityLevel,          // Scalar- Verbosity level of complaining
   const int *numCategoriesToCheck,    // Scalar- number of categories that will be checked
   const int *listOfCategoriesToCheck, // array[numCategoriesToCheck] - integers corresponding to the data
                                       // categories that are to be checked
   int *numberOfBadCategories);        // Scalar- Number of bad input data categories

// create a Teton BoundaryList object
void teton_constructboundary(
   const int *numReflecting, // scalar, number of reflecting BC on this domain
   const int *numVacuum,     // scalar, number of vacuum BC on this domain
   const int *numSource,     // scalar, number of source/profile BC  on this domain
   const int *numShared);    // scalar, number of boundaries this domain has with other processor domains

// set-up Teton's Compton scattering options
// comptonFlag options: (See tetonComptonFlag enum below)
//   none      -- no Compton scattering
//   FP        -- Fokker-Planck
//   Boltzmann -- Boltzmann
//   Thomson   -- Thomson
// To provide a scattering cross section without the energy-exchange physics,
//   use tetonComptonFlag::none and provide the scattering XS in teton_setopacity
// Two wrappers are provided for host code convenience.
void teton_constructcomptoncontrol(const int *comptonFlag); // scalar flag to determine the type of Compton scattering
inline void teton_constructcomptoncontrol_wrapper(tetonComptonFlag comptonFlag_enum)
{
   int comptonFlag = static_cast<int>(comptonFlag_enum);
   teton_constructcomptoncontrol(&comptonFlag);
}

// construct Teton time step monitoring object
void teton_constructdtcontrols(
   const double *dtrad,  // scalar, initial radiation time step size
   const double *dtrmn,  // scalar, minimum radiation time step size
   const double *dtrmx,  // scalar, maximum radiation time step size
   const double *delte,  // scalar, maximum relative change in electron temperature  per time step
   const double *deltr); // scalar, change in Tr^4 per time step

//
//  teton_constructeditor
//
//  construct global radiation tallies
//
//  Notes: PolarSectorPowerEscape stores all bins of a given group consecutively.
//            Ex: PSPE[0] = pow_ab_0_g0 ; PSPE[1] = pow_ab_0_g1 ; ... PSPE[ngr] = pow_ab_1_g0 ...
//          Values input here are used to initialize all polar sector output tallies.
void teton_constructeditor(
   const int *ngr,                           // scalar, number of groups to be output
   const int *nAngleBins,                    // scalar, number of angular bins for polar sector output
   const double *spectrumAngleBinBoundaries, // array[ nAngleBins + 1 ],  ascending cosines required (?)
   const double *RadPowerEscape,             // array[ ngr ],  Values used to initialize escape tallies
   const double *RadPowerIncident,           // array[ ngr ],  Values used to initialize inflow tallies
   const double *PolarSectorPowerEscape); // array[ ngr * nAngleBins] used to initialize angle dependent escape tallies

// Initialize Teton's internal geometry structures
void teton_constructgeometry();

// set the Teton iteration parameters (initialize only)
void teton_constructitercontrols();

// allocate material data storage
void teton_constructmaterial(const bool *nonLTE);

// Construct phase-space (angle, group) decomposition Teton uses for parallelism within this process domain
void teton_constructphasespacesets(const bool *fromRestart);

// prepare phase space sets to store angular flux outputs in bins
// Assuming ngroups energy binning (taken from QuadratureList data structure)
void teton_constructradintensity();

//
//  teton_constructquadrature_new
//
//  build the QuadratureList structure
//
//  Notes: Specify quadrature sets for high-order Sn and GTA sweeps
//         Ask the Teton team what nSetsMaster and nSets should be (related to threading or not).
//
//          Teton will access (in Fortran) QuadDef for the high-order Sn sweep as:
//
//            QuadDef(1,1) = type
//            QuadDef(2,1) = order
//            QuadDef(3,1) = n polar levels
//            QuadDef(4,1) = azimuthal angles
//            QuadDef(5,1) = polar axis [1,2,or 3]
//
//          and output the number of angles as
//
//            QuadDef(6,1)
//
//          Since Fortran is column oriented in memory, in C/C++ define:
//
//            int* quadInfo = new int[6*2]
//            quadInfo[0] = type, 1D(1) ,  rz: (1,2), xyz: (1,2,3)
//            quadInfo[1] = order (valid for level symmetric considerations)
//            quadInfo[2] = polar angles
//            quadInfo[3] = azimuthal angles
//            quadInfo[4] = polar axis
//
//           The second 6 entries in quadInfo is used to define the grey acceleration step angular quadrature
//
//           For example, one can set the grey quadrature as:
//             quaddef[6]  = 1;
//             quaddef[7]  = GTAorder;
//             quaddef[8]  = 1;
//             quaddef[9]  = 1;
//             quaddef[10] = 1;
//
//         The returned value of nSets is NOT the number of sets that Teton will actually use.
//         Host codes should NOT use nSets for anything.  We will probably remove this in the future.
//         Please see the Teton team for more information.
//
void teton_constructquadrature_new(
   const int *nSetsMaster, // scalar, see pnowak for more information
   int *nSets,             // scalar, don't use this, see pnowak for more information
   int *QuadDef,           // array[ 2 x 6] of angular quadrature defining data
   const double *gnu);     // array[ ngrps + 1] that gives frequency group boundaries (lowest to highest)

// Initialize the Size module.  This module knows all there is to know when driving Teton
// 11/19/2020 Update from Branson:
//     - Parameters controlling the behavior of the NLTE source function have changed:
//        - methodNLTE, planckCutoffNLTE1, and planckCutoffNLTE2 have been removed
//        - functionRNLTE controls how the Planck multiplier changes with temperature.
//          It takes values 1 (constant, recommended), 2 (linear), and 3 (nearLTE ansatz).
//     - These changes reflect updates to PhysicsUtils that improve stability of the rad/material coupling.
void teton_constructsize(
   const int *myRankInGroup, // scalar, rank in this process's communication group
   const int *nzones,        // scalar, number of zones this process is responsible for
   const int *ncornr,        // scalar, total number of corners in this domain
   const int *nSides,        // scalar, total number of sides in this domain
   const int *nbelem,        // scalar, total number of boundary elements (whole faces, not corner faces) on this domain
   const int *maxcf,         // scalar, maximum number of corner faces per corner
   const int *maxCorner,     // scalar, maximum number of corners in a zone
   const int *ncomm,         // scalar, number of processes this process communicates to
   const int *ndim,          // scalar, number of spatial dimensions
   const int *ngrp,          // scalar, number of energy groups
   const int *functionRNLTE, // scalar, integer equal to 1, 2, or 3 for controlling NLTE physics (see above comment)
   const double *tFloor,     // scalar, temperature floor
   const double *radForceMultiplier,  // scalar, if set to 0., this would turn off radiation feedback into hydro forces
   const double *betaNLTE,            // scalar, beta used in functionRNLTE=3
   const double *gammaNLTE,           // scalar, gamma for NLTE non-Planckian power law component
   const bool *dopplerShiftOn,        // scalar, account for frequency shifts due to material motion
   const bool *useNewNonLinearSolver, // scalar, this should only be used in 1D and experimentally (?)
   const bool *useNewGTASolver,       // scalar, selects new GTA solver
   const bool *usePWLD,               // scalar, turns on PWLD spatial discretization (rz only)
   const bool *useSurfaceMassLumping, // scalar, turns on surface & mass lumping for PWLD
   const bool *useGPU,                // scalar, offload calculation to the GPU (when available)
   const bool *useCUDASweep,          // scalar, use CUDA solver on the GPU
   const bool *useCUDASolver,         // scalar, use CUDA solver on the GPU
   const int *zoneBatchSize,          // scalar, number of zones to process at a time in NL Compton solver kernel.
   const int *nConcurrentBatches,     // scalar, number of NL zone batches to process concurrently.
   const int *igeom); // scalar, int corresponding to tetonGeometryFlag underlying integer for geometry type

// Set some optional parameters for GTA convergence.
//   enforceHardGTAIterMax defaults to true for the new GTA.  This is needed to avoid buildup of BiCGSTAB roundoff errors that leads to divergence.
//   forceExtraOuter defaults to false in all cases.  This may be needed to prevent the outer iteration from exiting earlier than it should when GTA is struggling.
void teton_setgtaoptions(
   const bool *
      enforceHardGTAIterMax, // bool, whether or not GTA exits at grey_max_sweeps, even if minimum convergence criterion is not met.
   const bool
      *forceExtraOuter); // bool, whether or not the outer is forced to do an extra iteration when GTA isn't converged

// Set verbosity level for controlling output
// verbosity is a two digit number.
//   The first digit determines whether only rank 0 is verbose
//   The second digit determines the verbosity level.
//     verbose=x0 - no verbose output
//     verbose=0x - rank 0 at verbose level x
//     verbose=1x - all ranks at verbose level x
void teton_setverbose(const int *verbose);

// Set sweep version.
// 0 - zone sweep (original)
// 1 - corner sweep (improved parallelism)
void teton_setsweepversion(const int *sweepversion);

// Set number of hyper-domains, for sweep and/or new GTA.
// 0 - Let Teton automatically decide. (default)
// >= 1 - User sets numer of hyper-domains to this number.
void teton_setsweepnumhyperdomains(const int *hyperdomains);
void teton_setgtanumhyperdomains(const int *hyperdomains);

// construct Teton memory allocator.
void teton_constructmemoryallocator(
   int *
      umpire_host_pinned_pool_allocator_id, // scalar, umpire allocator to use for pinned memory allocations on host (cpu).
                                            // value < 0 = have teton create its own allocator, if needed.
   int *umpire_device_pool_allocator_id); // scalar, umpire allocator to use for memory allocations on the accelerator
// device (gpu).  value < 0 = use native OpenMP or CUDA for memory allocation.

// destruct the Teton memory allocator
void teton_destructmemoryallocator();

// destroy mesh data in case it has changed during the simulation
void teton_destructmeshdata(const bool *nonLTE);

//
//  teton_dropvariables
//
//  Open a silo data file in the "/_Teton/" directory (relative to run directory) and drop Teton variables
//
//  Notes: filePtrId is an integer file number for Fortran, chosen by the host, something along the lines:
//           int fortranFileID = DBFortranAllocPointer(restartFile);
//
//         lenTe is the length / number of corner electron temperatures.
//         cornerTe is an array of all corner electron temperatures given to Teton.
//
//         success == -1 , if something goes wrong in this process
//
//         intensity is sorted into directories based on phase space sets
//
//         descriptors of each phase space set are given as integer codes:
//           setDesc(1) = number of groups in set
//           setDesc(2) = number of corners in this set
//           setDesc(3) = number of angles in the set
//           setDesc(4) = lowest energy group number (g0) in this set
//           setDesc(5) = lowest angle number
//           setDesc(6) = quadrature ID
//
//         On 7/28/2017, lenTe was moved to go before cornerTe array for consistency (to self and Fortran standards)
//
void teton_dropvariables(const int *filePtrID, // scalar, integer file index for Fortran to open
                         const int *lenTe, // scalar, total number of electron temperatures = total number of corners
                         const double *cornerTe, // array[lenTe] of electron temperatures
                         int *success);          // scalar, flag for completion or not

//
//  teton_dtnew
//
//  Determine the time step size Teton would suggest based on this cycle's information
//
//  Notes: Input a Compton change parameter here, but calculate dt suggestions of all other factors Teton considers when voting for dt
//         Retrieve time step votes via teton_getrunstats()
void teton_dtnew(const int *maxOSComptonChangeCorner,
                 const double *maxOSComptonChange); //scalar, relative change parameter

// Get angle bins
void teton_getanglebins(const int *numAngleBins, // scalar input, must be equal to what teton_getnumanglebins returns
                        double *angleBinBoundaries); // output array[ numAngleBins + 1 ]

// a function for getting overview information about the most recent cycle
//
// Note: This does NOT actually compute the tallies.  That is done in
//  teton_rtedit.  This only returns stored values from the last
//  teton_rtedit call.  If psi is changed between teton_rtedit and
//  teton_getedits, it will not be reflected in the tally quantities below.
void teton_getedits(int *noutrt,             // scalar, number of thermal [outer] iterations in this cycle
                    int *ninrt,              // scalar, number of inner [transport] iterations in this cycle
                    int *ngdart,             // scalar, number of grey diffusion iterations
                    int *nNLIters,           // scalar, number of nonlinear iterations
                    int *maxNLIters,         // scalar, maximum nonlinear iterations
                    int *TrMaxZone,          // scalar, zone id with maximum T rad
                    int *TeMaxZone,          // scalar, zone id with maximum electron temperature
                    int *TrMaxProcess,       // scalar, MPI process with maximum T rad
                    int *TeMaxProcess,       // scalar, MPI process with maximum electron temperature
                    double *dtused,          // scalar, Teton returns the time step used in cycle just completed
                    double *dtrad,           // scalar, radiation vote for next time step
                    double *TrMax,           // scalar, T rad maximum value
                    double *TeMax,           // scalar, T electron maximum value
                    double *deltaEMat,       // scalar, change in material energy this time step
                    double *EnergyRadiation, // scalar, energy contained in the radiation field
                    double *PowerIncident,   // scalar, power of photons incident
                    double *PowerEscape,     // scalar, power of photons escaping
                    double *PowerAbsorbed,   // scalar, power of energy absorbed
                    double *PowerEmitted,    // scalar, power of photons emitted
                    double *PowerExtSources, // scalar, power of photons from fixed volumetric sources
                    double *PowerCompton,    // scalar, power of Compton scattering photons
                    double *EnergyCheck);    // scalar, total energy not accounted for this cycle.

// same as teton_getedits, but this function stores
// the various edits in to Teton's internal conduit node
// instead of returning them in the arguments
void teton_publishedits(double *dtrad); // scalar, radiation vote for next time step
//
//  teton_getemissionsource
//
//  Retrieve the average energy sourced into every zone and group
//
//  Notes: sourceMatrix is indexed as:
//   sourceMatrix[0] = s_z0_g0
//   sourceMatrix[1] = s_z0_g1 ...
//   sourceMatrix[ngroups] = s_z1_g0 ...
void teton_getemissionsource(double *sourceMatrix); // array[numZonesTotal*numGroupsTotal]

// Get the number of angle bins
void teton_getnumanglebins(int *numAngleBins); //returns a scalar integer equal to the number of polar angle bins

// have Teton re-calculate corner volumes and face areas after the mesh moves or changes
void teton_getvolume();

// given a group-angle phase space index, return its global group and angle number offsets
void teton_getgroupangleoffsets(int *setIdx,  // scalar, input
                                int *g0,      // scalar, returned
                                int *angle0); // scalar, returned

// get the total opacity of a given zone and group
void teton_getopacity(const int *zone, // scalar, zone index
                      const int *grp,  // scalar, group desired
                      double *sigT);   // scalar, returned value

// Zero Psi, PhiTotal, and RadEnergyDensity in void zones; tallies energy as escaped.
void teton_resetpsi(const int *nVoidZones,    // number of void zones
                    const int *voidZoneList); // list of void zone IDs

// Get the number of sets.
void teton_getnumsnsets(int *numSNSets); // number of SN sets

//
//  teton_getpsipointer
//
//  get a pointer to the contiguous memory space of the intensity solution of a given phase space set
//
//  Notes: Specify a set index, and pass a length 3 array specifying the dimensions of that set
//         The dimension specifying array, psiDims is used as:
//           psiDims[0] = number of groups in set, ngrps
//           psiDims[1] = number of corners in this set, ncorners
//           psiDims[2] = number of angles in this set, nAngles
//         Call from C/C++ as:
//           int idx; int psiDims[3]; double* intensity;
//           teton_getpsipointer(&idx , psiDims , &intensity);
void teton_getpsipointer(const int *setIdx,   // scalar, phase space set desired (input)
                         int *psiDims,        // array[3], dimensions of the intensity object pointed at
                         double **intensity); // array[ngrps*ncorners*nAngles], indexed as ( )?

//
//  teton_getzonalpsi
//
//  get the zone-averaged psi for each zone, group, and angle
//
void teton_getzonalpsi(
   const int *numAngles, // scalar, number of angles, must match psi
   double *psi); // array[nzones*ngroups*nAngles], indexed such that zone varies fastest and angle varies slowest

// get the amount of radiation deposited in a zone within a cycle
void teton_getradiationdeposited(
   const int *zoneID,    // scalar, zone index
   double *radDeposited, // scalar, amount of radiation energy deposited
   double *
      radTemperature); // scalar, equivalent black body temperature of radiation energy deposited (for the given zone material)

//
//  teton_getradiationenergydensity
//
//  get a vector of the average radiation energy density  in every zone, for every group
//
//  Notes: RadEnergyDensity is indexed as a (zoneIdx, groupIdx) quantity in Fortran
//         In C/C++ this would be accessed using:
//           radEnergyDensity[0] = radE_z0_g0
//           radEnergyDensity[1] = radE_z1_g0 ....
void teton_getradiationenergydensity(double *radEnergyDensity); // array[nZones * ngroups]

//
//  teton_getradiationflux
//
//  Retrieve a zone's volume averaged or integrated radiation flux
//
//  Notes: radFlux is an array of length (nSpatialDim x nGroups) and is ordered in C/C++ as
//           radFlux[0] = flux_x_g0
//           radFlux[1] = flux_y_g0
//           radFlux[2] = flux_z_g0
//           radFlux[3] = flux_x_g1
void teton_getradiationflux(const int *zoneIndex, // scalar, zone desired
                            double *radFlux);     // array[nSpatialDim*nGroups]

//
//  teton_getbeammetric
//
//  Retrieve a zone's volume averaged, group integrated maximum diagonal of
//  the Eddington tensor.  This value will be 1/3 for isotropic fluxes and 1
//  for a beam.
void teton_getbeammetric(const int *zoneIndex, // scalar, zone desired
                         double *beammetric);  // scalar, double

//
//  teton_getradiationforce
//
//  get the radiation force of every corner in zone zoneIndex
//
//  Notes: matrixRadForce is indexed as
//           matrixRadForce[0] = radF_c0_x
//           matrixRadForce[1] = radF_c0_y
void teton_getradiationforce(const int *zoneIndex,    // scalar, zone index
                             double *matrixRadForce); // array[ nSpatialDim*maxCornersPerZone ]

// get the volumes of the corners in zone zoneIdx
void teton_getcornervolumes(int *zoneIdx,           // scalar, zone index
                            double *cornerVolumes); // array[ maxCornersPerZone ]

// get the radiation temperature of a single zone
void teton_getradiationtemperature(int *zoneIdx,     // scalar, zone index
                                   double *radTemp); // scalar, radiation temperature

// get the material temperature of a single zone
void teton_getmaterialtemperature(int *zoneIdx,     // scalar, zone index
                                  double *matTemp); // scalar, material temperature

// a function for accessing Teton cumulative timers
void teton_getrunstats(double *MatCoupTimeTotal,  // scalar, Cumulative time spent in material coupling
                       double *SweepTimeTotal,    // scalar, Cumulative time spent sweeping
                       double *GPUSweepTimeTotal, // scalar, Cumulative time spent sweeping on the GPU
                       double *GTATimeTotal,      // scalar, Cumulative time spent in grey acceleration
                       double *RadtrTimeTotal,    // scalar, Cumulative time spent in radtr
                       double *InitTimeTotal,     // scalar, Cumulative time spent in initialize sets
                       double *FinalTimeTotal,    // scalar, Cumulative time spent in finalize sets
                       double *timeNonRad,
                       double *timeOther);

// a function for accessing Teton's dt control information in the form of a string:
void teton_getdtmessage(char **dtmessage);

// a function for accessing Teton's dt control information:
void teton_getdtcontrolinfo(
   int *dtControlReason,  // scalar, reason for control, see tetonDtControlFlag below
   int *dtControlProcess, // scalar, controlling process MPI rank (0-based indexing)
   int *dtControlZone);   //scalar, local controlling zone ID on process dtControlProcess (1-based indexing)

//
//  teton_getscalarintensity
//
//  get the multi-group scalar intensity of every corner
//
//  Notes: Indexed as:
//           scalarIntensity[0] =
//           scalarIntensity[1] =
void teton_getscalarintensity(double *scalarIntensity); // array[ngroups*ncornersTotal]

// initialize electron temperature to values in tec array.  Zero out cve, rho, tez, nez, trz of every zone
void teton_initmaterial(const double *tec); //  array[ nCornerTotal ]

// begin initializing Teton opacity structures
void teton_initopacity();

// begin initializing Teton NLTE fields
void teton_initnltefields();

//
//  teton_initteton
//
//  Initialize radiation intensity and corner electron temperatures existing material values
//
//  Notes:  Unknowns are reset to the zone-wise constant values associated with that zone's material values of trz and tez
void teton_initteton(double *radEnergyTotal, // scalar, total radiation energy given the input electron temperatures
                     double *cornerTe); // array[ numCornersTotal ] the assigned electron temperature in each corner

//
//  teton_loadvariables
//
//  Initialize Teton intensity and electron temperature from an existing restart file
//
//  Notes:  filePtrID is the Fortran handle assigned by Silo in the C/C++ calling program
//          Electron temperatures from restart file are loaded by this function into the cornerTe array
void teton_loadvariables(const int *filePtrID, // scalar, corresponding to host code's Silo-based restart file
                         const int *lenTe,     // scalar, numCorners!
                         double *cornerTe,     // array[ lenTe ],
                         int *nSetsMaster,     // scalar,
                         int *success); // scalar, flag to indicate successful opening and reading of restart file

//
//  teton_normalizematerial
//
//  Scale certain material properties
//
//  Notes: This function scales the following ways:
//           te_zone   = te_zone   / cv_e           -- roughly internal energy
//           trad_zone = trad_zone / rho            -- radiation energy per mass
//           cv_e      = max(cv_e_floor, cv_e/rho)  -- Specific heat
//
//         cv_e floor is set at 1.0 x 10^{-8}
void teton_normalizematerial();

//
//  teton_querypsi
//
//  Query a single value of the radiation intensity field.
//
//  Notes: Uses C/C++ (0-based) indexing
//         idxArray is a length three array
//           idxArray[0] = group (global value)
//           idxArray[1] = cornerId (global value)
//           idxArray[2] = angleID (global index value)
void teton_querypsi(const int *idxArray, // array [3]
                    double *value);      // scalar, requested intensity value

// Execute a Teton solve for a given time cycle
void teton_radtr();

//
//  teton_reconstructpsi
//
//  Reconstruct/scale the corner based scalar intensity (Phi) after remap to match the new, zone-based radiation energy density
//
//  Notes: radDensity is treated as a (zone , group ) dimension matrix
//           radDensity = [eRad_z0_g0, ... eRad_zN-1_g0, eRad_z0_g1 ... ]
//         phiTotal is indexed as (group , globalCorner Idx).
//           phiTotal = [phi_g0_c0, phi_g1_c0 ... phi_gN_c0 , phi_g0_c1 ...]
void teton_reconstructpsi(double *EnergyRadiation,   // scalar, total radiation energy density (in/out)
                          const double *radDensity); // array[ nzones x ngroups ], input only

//
//  teton_reconstructpsifrompsi
//
//  Reconstruct/scale the corner based psi after remap to match the new, zone-based remapped psi
//
//  TODO add a note about the ordering of RemappedPsi vs Set%Psi
//
void teton_reconstructpsifrompsi(
   double *EnergyRadiation,    // scalar, total radiation energy density (in/out)
   const int *numAngles,       // scalar, number of angles, must match psi
   const double *RemappedPsi); // array[ nzones x ngroups x numAngles], input only, zone varies fastest

//
//  teton_reconstructpsifromdv
//
//  Reconstruct/scale the corner based scalar intensity (Phi) after Lagrange mesh motion, to match the new, zone-based radiation energy density
void teton_reconstructpsifromdv();

// Change SOME iteration scheme parameters of Teton
// See teton_constructsize for more details
void teton_resetsize(
   const int *functionRNLTE,           // scalar, integer equal to 1, 2, or 3 for controlling NLTE physics
   const double *tFloor,               // scalar, temperature floor
   const double *radForceMultiplier,   // scalar, if set to 0., this would turn off radiation feedback into hydro forces
   const double *betaNLTE,             // scalar, beta used in functionRNLTE=3
   const double *gammaNLTE,            // scalar, gamma power law used for NLTE non-Planckian component
   const bool *dopplerShiftOn,         // scalar, Account for frequency shifts due to material motion
   const bool *useNewNonLinearSolver,  // scalar, enable new non-linear solver
   const bool *useNewGTASolver,        // scalar, enable new GTA solver
   const bool *usePWLD,                // scalar, use PWLD spatial discretization (rz only)
   const bool *useSurfaceMassLumping); // scalar, use surface and mass lumping for PWLD

// reset internal, total duration, Teton timers to the values passed in
void teton_resettimers(const double *radTrTotal,        // scalar, total wall clock time used by Teton
                       const double *matCouplingTotal,  // scalar, total time spent in material coupling
                       const double *sweepTimeTotal,    // scalar, total time spent sweeping
                       const double *gpuSweepTimeTotal, // scalar, total time spent sweeping on the GPU
                       const double *greyAccelTotal,    // scalar, total time spent in grey acceleration step
                       const double *initTotal,         // scalar, total time spent in initialize sets
                       const double *finalTotal);       // scalar, total time spent in finalize sets

//
//  teton_rtedit
//
//  Trigger the calculation of all radiative transfer edits at end of a cycle.
//
//  Notes: Write the new value of corner electron temperatures into tElec array.
//         Other edits must be extracted with different function calls
//
//    This isn't where RadPowerIncident, RadPowerEscape, or
//    PolarSectorPowerEscape are computed.  Those three are computed in
//    BoundaryEdit.F90, which is called at the end of teton_radtr.  For those
//    three tallies, they are only summed over the sets here.  teton_rtedit
//    does compute the remaining tallies (not the aforementioned 3) that can
//    be obtained in teton_getedits.
//
//    In other words, RadPowerIncident, RadPowerEscape, and
//    PolarSectorPowerEscape would not be affected by calling teton_resetpsi
//    immediately before teton_rtedit.  The other tallies in
//    teton_getedits may be affected by calling teton_resetpsi immediately
//    before teton_rtedit.  None of the tallies in teton_getedits would be
//    affected if teton_resetpsi is called after teton_rtedit.
void teton_rtedit(double *tElec); // array[ nCornersTotal ]

//
//  teton_scalepsir
//
//  Isotropically scale radiation intensity to match the specified change in scalar / angle integrated intensity
//
//  Notes: deltaPhi is assumed to have an order of array of
//           deltaPhi[0] = dPhi_g0_c0
//           deltaPhi[1] = dPhi_g1_c0 ...
void teton_scalepsir(const double *deltaPhi); // array[ nGrp * nCornersTotal ]

//
//  teton_setcommunicationgroup
//
//  Set the communicator group that Teton should use
//
//  Notes: Input argument should be the result of the MPI_Comm_c2f( MPI_Comm comm) function
//         MPI_Comm_c2f should return a value that is already of a FORTRAN integer type
void teton_setcommunicationgroup(const MPI_Fint *convertedCommunicatorHandle); // scalar

// Update material in an individual zone by ADDING the passed in values
void teton_setmaterial(const int *zone,                // scalar, zone number
                       const double *cve,              // scalar, "effective" specific heat
                       const double *rho,              // scalar, density
                       const double *tez,              // scalar, electron temperature
                       const double *trz,              // scalar, radiation temperature
                       const double *nez,              // scalar, electron density
                       const double *stimComptonMult); // scalar, stimulated Compton multiplier

// set an external, zone-wise constant, volumetric energy source
void teton_setmaterialsource(const double *sMat); // array[ nZones ]

//
//  teton_setnodeposition
//
//  Set position of a zone's nodes
//
//  Notes: Position matrix should be used as:
//           positionMatrix = new double[maxCorners*nSpatialDim]
//           positionMatrix[0] = x_v0
//           positionMatrix[1] = y_v0
//           positionMatrix[2] = z_v0
//           positionMatrix[3] = x_v1
//         Ordering corresponds to Teton's ordering of nodes/corners of a given zone
void teton_setnodeposition(const int *zone,               // scalar
                           const double *positionMatrix); // array[ nDim * nCornersMax ]

//
//  teton_setnodevelocity
//
//  Set the velocity of mesh nodes
//
//  Notes: Velocity matrix should be used as:
//           velocityMatrix = new double[maxCorners*nSpatialDim]
//           velocityMatrix[0] = v_x_0
//           velocityMatrix[1] = v_y_0
//           velocityMatrix[2] = v_z_0
//           velocityMatrix[3] = v_x_1
//         Ordering corresponds to Teton's ordering of nodes for a given zone
void teton_setnodevelocity(const int *zone,               // scalar, integer of zone to set
                           const double *velocityMatrix); // array[ nDim * nCornersMax ]

// set the absorption and scattering opacity of every group of a zone
//   Note that opacities are clipped at 1.e50.  We may make the clip adjustable by the host-code in the future. See TETON-131.
void teton_setopacity(
   const int *zone,    // scalar, zone number
   const double *siga, // array[ nGrps ] ordered from low to high with absorption opacities
   const double *sigs, // array[ nGrps ] ordered from low to high with scattering opacities
   const bool
      *useTableSigmaS); // if useTableSigmaS==true, then sigs is added to the total opacity here.  else, ignore sigs

// set the opposite face array within the MeshData module
void teton_setoppositeface();
// set the NLTE fields
void teton_setnltefields(
   const int *zone,        // scalar, zone number
   const double *emis,     // array[ nGrps ] ordered from low to high with emissivity
   const double *demisdT); // array[ nGrps ] ordered from low to high with emissivity temperature derivative

// set the NLTE flag
void teton_setnlteflag(const int *zone,     // scalar, zone number
                       const bool *nonLTE); // Boolean, whether the zone is NLTE

// instruct Teton to calculate the radiation flux prior to retrieval
void teton_setradiationflux();

// have Teton calculate the force exhibited by the radiation field on the material
void teton_setradiationforce();

// Have Teton calculate the total out-scatter coefficient for all zones.
//
// comptonFlag should match the value passed to teton_constructcomptoncontrol
// If the flags don't match, a warning is given.
//
// This routine will overwrite the Mat% sigs variable.  If sigs's were
// passed to Teton in teton_setopacity, they will be overwritten by this call.
//
// This routine should probably not be called if setOpacity was called with
// useTableSigmaS == true.
void teton_setscatteringopacity(const int *comptonFlag); // scalar flag to determine the type of Compton scattering
inline void teton_setscatteringopacity_wrapper(tetonComptonFlag comptonFlag_enum)
{
   int comptonFlag = static_cast<int>(comptonFlag_enum);
   teton_setscatteringopacity(&comptonFlag);
}

// add a zone face to a shared (processor subdomain boundary) condition
void teton_setsharedface(const int *bcID,   // scalar, input, index of the boundary we are setting up
                         const int *zoneID, // scalar, global index of zone at this boundary (in C/host code indexing)
                         const int *faceID, // scalar,  face index of this element on the  boundary
                         const int *cornerID); // scalar,  global index of corner at this boundary

//
//  teton_setzone
//
//  Set-up a given zone's topology for Teton
//
//  Notes:  CornerFaces is the total number of corner faces of a given type (zone or inter-corner) faces in the zone
//          Each corner always has an equal number of faces that match with corners external to the given zone (zone corner faces)
//          and faces that separate it from adjacent corners of the zone within it lies (inter-corner corner faces).
//            Example:  Hexahedral conforming mesh.  Each zone has 8 corners.
//                      Each corner has 6 corner faces, 3 corner faces that lie on the zone boundary and
//                      3 corner faces that delineate between corner volumes of adjacent corners in the same zone.
//                      cornerFaces = 8*3 = 24
//
//          zoneNCorner is the number of corners in this zone.
//            1 per corner per node that defines the zone.
//            Hanging nodes are treated as a zone defining node by Teton)
//
//          zoneOpp array indicates if a zone face is a domain boundary (-1) or gives the zone index of zone on the other side of the face
//
//          CornerID an array that associates internal zone corner numbering with local face numbering
//
//          CornerOpp an array that associates a neighboring zone's corner indexing to the  face in question
//
//          nCPerFace array that gives the number of corners that touch a given zone face
//
//          FaceToBcList integers array to associate zone faces that are boundaries to the appropriate BC index
void teton_setzone(
   const int *zoneID,        // scalar, index of zone
   const int *corner0,       // scalar, global index of this zone's local corner zero
   const int *zoneFaces,     // scalar, number of zone faces this zone has that connect to other zones or boundaries
   const int *cornerFaces,   // scalar,
   const int *zoneNCorner,   // scalar, number of corners in this zone
   const int *zoneOpp,       // array[ zoneFaces ] of zones opposite this zone on a given face
   const int *CornerID,      // array[ cornerFaces ]
   int *CornerOpp,           // array[ cornerFaces ]
   const int *nCPerFace,     // array[ zoneFaces ]
   const int *FaceToBCList); // array[ zoneFaces ]

// set up 1D mesh processor layout
void teton_setzone1d(
   const int *zoneID,     // scalar, zone index we are setting up
   const int *numBcTotal, // scalar, total number of BC in the problem (boundaries can be communicator boundaries)
   const int *BcZoneID);  // array[ numBcTotal ] giving communicator ID of shared boundaries

// set time step parameters for Teton
void teton_settimestep(const int *ncycle,      // scalar, current cycle
                       const double *dtRad,    // scalar, time step size for radiation
                       const double *timeRad,  // scalar, physical time  of this cycle evaluation
                       const double *minTemp); // scalar, temperature floor

// Don't call this directly, use one of the two teton_surfaceedit interfaces below
void teton_surfaceedit_internal(const int *nCornerFaces,
                                const bool *labFrame,
                                const int *cornerList,
                                const int *zoneFaceList,
                                const bool *timeShift,
                                const double *centerPoint,
                                const int *numAngleBins,
                                const int *numGroups,
                                const int *numTimeBins,
                                const double *timeBinBoundaries,
                                const bool *computeIncident,
                                const double *scaleTally,
                                const bool *calcErrorMetrics,
                                double *tally,
                                double *tallyIncident,
                                double *err_est_shift,
                                double *err_est_srcsize);
}

// Old API with deprecated verbose for backward compatibility:
inline void teton_getrunstats(const int *verbose,        // deprecated, use teton_setverbose instead
                              double *MatCoupTimeTotal,  // scalar, Cumulative time spent in material coupling
                              double *SweepTimeTotal,    // scalar, Cumulative time spent sweeping
                              double *GPUSweepTimeTotal, // scalar, Cumulative time spent sweeping on the GPU
                              double *GTATimeTotal,      // scalar, Cumulative time spent in grey acceleration
                              double *RadtrTimeTotal,    // scalar, Cumulative time spent in radtr
                              double *InitTimeTotal,     // scalar, Cumulative time spent in initialize sets
                              double *FinalTimeTotal,    // scalar, Cumulative time spent in finalize sets
                              double *timeNonRad,
                              double *timeOther)
{
   (void) verbose;
   teton_getrunstats(MatCoupTimeTotal,
                     SweepTimeTotal,
                     GPUSweepTimeTotal,
                     GTATimeTotal,
                     RadtrTimeTotal,
                     InitTimeTotal,
                     FinalTimeTotal,
                     timeNonRad,
                     timeOther);
}

// Old API for backward compatibility
inline void teton_constructsize(
   const int *myRankInGroup, // scalar, rank in this process's communication group
   const int *nzones,        // scalar, number of zones this process is responsible for
   const int *ncornr,        // scalar, total number of corners in this domain
   const int *nSides,        // scalar, total number of sides in this domain
   const int *nbelem,        // scalar, total number of boundary elements (whole faces, not corner faces) on this domain
   const int *maxcf,         // scalar, maximum number of corner faces per corner
   const int *maxCorner,     // scalar, maximum number of corners in a zone
   const int *ncomm,         // scalar, number of processes this process communicates to
   const int *ndim,          // scalar, number of spatial dimensions
   const int *ngrp,          // scalar, number of energy groups
   const int *maxDynamicIters, // deprecated, will be removed next release.  use teton_adjustfluxexchangegoals instead
   const int *functionRNLTE,   // scalar, integer equal to 1, 2, or 3 for controlling NLTE physics (see above comment)
   const double *tFloor,       // scalar, temperature floor
   const double *radForceMultiplier,  // scalar, if set to 0., this would turn off radiation feedback into hydro forces
   const double *betaNLTE,            // scalar, beta used in functionRNLTE=3
   const double *gammaNLTE,           // scalar, gamma for NLTE non-Planckian power law component
   const bool *dopplerShiftOn,        // scalar, account for frequency shifts due to material motion
   const bool *useNewNonLinearSolver, // scalar, this should only be used in 1D and experimentally (?)
   const bool *useNewGTASolver,       // scalar, selects new GTA solver
   const bool *usePWLD,               // scalar, turns on PWLD spatial discretization (rz only)
   const bool *useSurfaceMassLumping, // scalar, turns on surface & mass lumping for PWLD
   const bool *useGPU,                // scalar, offload calculation to the GPU (when available)
   const bool *useCUDASweep,          // scalar, use CUDA solver on the GPU
   const bool *useCUDASolver,         // scalar, use CUDA solver on the GPU
   const int *zoneBatchSize,          // scalar, number of zones to process at a time in NL Compton solver kernel.
   const int *nConcurrentBatches,     // scalar, number of NL zone batches to process concurrently.
   const int *igeom) // scalar, int corresponding to tetonGeometryFlag underlying integer for geometry type
{
   (void) maxDynamicIters;
   teton_constructsize(myRankInGroup,
                       nzones,
                       ncornr,
                       nSides,
                       nbelem,
                       maxcf,
                       maxCorner,
                       ncomm,
                       ndim,
                       ngrp,
                       functionRNLTE,
                       tFloor,
                       radForceMultiplier,
                       betaNLTE,
                       gammaNLTE,
                       dopplerShiftOn,
                       useNewNonLinearSolver,
                       useNewGTASolver,
                       usePWLD,
                       useSurfaceMassLumping,
                       useGPU,
                       useCUDASweep,
                       useCUDASolver,
                       zoneBatchSize,
                       nConcurrentBatches,
                       igeom);
}

// old version kept for backward compatibility
inline void teton_resetsize(
   const int *maxDynamicIters, // deprecated, will be removed next release.  use teton_adjustfluxexchangegoals instead
   const int *functionRNLTE,   // scalar, integer equal to 1, 2, or 3 for controlling NLTE physics
   const double *tFloor,       // scalar, temperature floor
   const double *radForceMultiplier,  // scalar, if set to 0., this would turn off radiation feedback into hydro forces
   const double *betaNLTE,            // scalar, beta used in functionRNLTE=3
   const double *gammaNLTE,           // scalar, gamma power law used for NLTE non-Planckian component
   const bool *dopplerShiftOn,        // scalar, Account for frequency shifts due to material motion
   const bool *useNewNonLinearSolver, // scalar, enable new non-linear solver
   const bool *useNewGTASolver,       // scalar, enable new GTA solver
   const bool *usePWLD,               // scalar, use PWLD spatial discretization (rz only)
   const bool *useSurfaceMassLumping) // scalar, use surface and mass lumping for PWLD
{
   (void) maxDynamicIters;
   teton_resetsize(functionRNLTE,
                   tFloor,
                   radForceMultiplier,
                   betaNLTE,
                   gammaNLTE,
                   dopplerShiftOn,
                   useNewNonLinearSolver,
                   useNewGTASolver,
                   usePWLD,
                   useSurfaceMassLumping);
}

// Compute time-shifted, angle-resolved radiation tallies on an arbitrary surface
//
// Note: Unlike teton_getedits, this function call performs the tally
//  calculation.  It doesn't just return the values.  This means that the
//  returned quantities will reflect Teton's current stored copy of psi.
//
// tally is incremented with by the energy traveling through the specified surface in this time step.
//   To obtain the energy for this time step only, pass in a zeroed array.
// Groups vary the fastest in the array.  Time varies the slowest.
//   Total leaked energy as a scalar would be sum over the entire length of tally.
//   Should be able to sum over any dimension to get what would have happened if that bin was only one wide.
// tally and tallyIncident are positive unless there are negative fluxes present.
//   The net outgoing/escape energy is given by tally - tallyIncident.
//   Which direction is outgoing/"escape" is determined by the corner ID.  The corner-face is outgoing relative to the corner ID specified.
//
// There are three estimates of the error resulting from the time shifting/approximation of the surface tally as
//   a point source.  For the err_est inputs, the first numTimeBins entries correspond to the estimates for the
//   escape tally, while the latter numTimeBins correspond to the incident tally.  The latter numTimeBins entries
//   are ignored unless computeIncident is true.  These error estimates are accumulated and NOT normalized.
//   To normalize them, divide by group-and-angle integrated tallies.
//   Note that they are NOT filled out if timeShift is False.
//   These arguments are ignored if timeShift is false or calcErrorMetrics is false.
//   If timeShift, the err_est arrays should be of length 2*numTimeBins if computeIncident, length numTimeBins otherwise.
//
inline void teton_surfaceedit(
   const int *nCornerFaces,   // scalar, number of corner-faces that comprise the surface
   const bool *labFrame,      // scalar, transform to lab frame (yes if true)
   const int *cornerList,     // array[ nCornerFaces ], corner IDs on of each corner-face
   const int *zoneFaceList,   // array[ nCornerFaces ], zone-local zone-face ID of each corner-face
   const bool *timeShift,     // scalar, whether to apply a time shift from centerPoint
   const double *centerPoint, // array[ nDim ], coordinates of center point.  ignored if timeShift = False
   const int *numAngleBins,   // either 1 (integrated over all angles) or exactly what teton_getnumanglebins returns
                              //const double* angleBinBoundaries, // Teton already has this
   const int *numGroups,      // scalar.  must be 1 (total energy) or number of groups (spectrum)
                              //const double* groupBinBoundaries, // Teton already has this
   const int *numTimeBins,    // scalar, number of time bins.  Can be 1.
   const double *
      timeBinBoundaries, // array[ numTimeBins+1 ], values must be ascending.  times before 0/after numTimeBins are put into first/last bin
   const bool *
      computeIncident, // scalar, true = append incident energy to tallyIncident, false = tallyIncident is not filled (can be null)
   const double *
      scaleTally, // scalar multiplier on tally.  This can be used to convert from Teton units to your favorite units.  Applied AFTER AllReduction of quantities.
   const bool *
      calcErrorMetrics, // scalar, true = fill out error metrics if timeShift is true, false = err_est_* args below should be NULL
   double *tally, // array[ numGroups * numAngleBins * numTimeBins], escape energy
   double *
      tallyIncident, // array[ numGroups * numAngleBins * numTimeBins], incident energy, ignored if tallyIncident is false
   double *
      err_est_shift, // array[ * ], see above for more info. error metric from distance of time shifting (See <d_{||}> in time-shifted tally notes.)
   double *
      err_est_srcsize) // array[ * ], see above for more info. error metric from size of source inside simulation domain (See <d_{\perp}> in time-shifted tally notes.)
{
   teton_surfaceedit_internal(nCornerFaces,
                              labFrame,
                              cornerList,
                              zoneFaceList,
                              timeShift,
                              centerPoint,
                              numAngleBins,
                              numGroups,
                              numTimeBins,
                              timeBinBoundaries,
                              computeIncident,
                              scaleTally,
                              calcErrorMetrics,
                              tally,
                              tallyIncident,
                              err_est_shift,
                              err_est_srcsize);
}

// Old surfaceedit interface for backward compatibility:
inline void teton_surfaceedit(
   const int *nCornerFaces,   // scalar, number of corner-faces that comprise the surface
   const bool *labFrame,      // scalar, transform to lab frame (yes if true)
   const int *cornerList,     // array[ nCornerFaces ], corner IDs on of each corner-face
   const int *zoneFaceList,   // array[ nCornerFaces ], zone-local zone-face ID of each corner-face
   const bool *timeShift,     // scalar, whether to apply a time shift from centerPoint
   const double *centerPoint, // array[ nDim ], coordinates of center point.  ignored if timeShift = False
   const int *numAngleBins,   // either 1 (integrated over all angles) or exactly what teton_getnumanglebins returns
   const int *numGroups,      // scalar.  must be 1 (total energy) or number of groups (spectrum)
   const int *numTimeBins,    // scalar, number of time bins.  Can be 1.
   const double *
      timeBinBoundaries, // array[ numTimeBins+1 ], values must be ascending.  times before 0/after numTimeBins are put into first/last bin
   const bool *
      computeIncident, // scalar, true = append incident energy to tallyIncident, false = tallyIncident is not filled (can be null)
   double *tally, // array[ numGroups * numAngleBins * numTimeBins], escape energy
   double *
      tallyIncident) // array[ numGroups * numAngleBins * numTimeBins], incident energy, ignored if tallyIncident is false
{
   const double scaleTally = 1.;
   const bool dumpError = false;

   teton_surfaceedit_internal(nCornerFaces,
                              labFrame,
                              cornerList,
                              zoneFaceList,
                              timeShift,
                              centerPoint,
                              numAngleBins,
                              numGroups,
                              numTimeBins,
                              timeBinBoundaries,
                              computeIncident,
                              &scaleTally,
                              &dumpError,
                              tally,
                              tallyIncident,
                              NULL,
                              NULL);
}

//
//  teton_addprofile
//
//  add a time dependent (step function with numtimes steps) boundary source to Teton
//
//  Notes: If not black body emission, must supply n_groups values for each time.
//          Teton expects input values with units of Energy/(Volume*PhotonEnergy) for nonBlackBody sources.
//          Booleans must be 0 or 1 if using int flags instead of BOOL type data.
//          Call after teton_constructboundary().
inline void teton_addprofile(
   const int *NumTimes,      // scalar, number of differnt times for the source
   const int *NumValues,     // scalar, ==numTimes if blackbody BC.  ==numTimes*ngroups if not blackbody BC
   const double *Multiplier, // scalar, multiplier in front of all emission terms
   const bool *BlackBody,    // scalar, isBlackBody emission spectrum?
   const bool *
      Isotropic, // scalar, isIsotropic angular emission?  Must be true unless using Lobatto quadrature.  False means radiation enters along the polar axis only.
   const double *Times,  // array[NumTimes], physical times source changes given in ascending order
   const double *Values, // array[NumValues], group dependent source strength values
   int *TetonProfileID)  // scalar output, integer index corresponding to the profile added (used in teton_resetprofile)
{
   teton_addprofile_internal(NumTimes, NumValues, Multiplier, BlackBody, Isotropic, Times, Values, TetonProfileID);
}
// Overloaded API for setting a constant (in time), isotropic, blackbody source:
inline void teton_addprofile(
   const double Value,  // scalar, blackbody source strength value
   int *TetonProfileID) // scalar output, integer index corresponding to the profile added (used in teton_resetprofile)
{
   const int NumTimes = 2;
   const int NumValues = 2;
   const bool BlackBody = true;
   const bool Isotropic = true;
   const double Multiplier = 1.;
   const double Times[2] = {-1.e50, 1.e50};
   const double Values[2] = {Value, Value};
   teton_addprofile_internal(&NumTimes, &NumValues, &Multiplier, &BlackBody, &Isotropic, Times, Values, TetonProfileID);
}
// Overloaded API for setting a constant (in time), isotropic, frequency-dependent source:
inline void teton_addprofile(
   const double *Multiplier,      // scalar, multiplier in front of all emission terms
   std::vector<double> values_in, // vector of length ngroups, source strength value
   int *TetonProfileID) // scalar output, integer index corresponding to the profile added (used in teton_resetprofile)
{
   const int NumTimes = 2;
   const int ngroups = values_in.size();
   const int NumValues = NumTimes * ngroups;
   const bool BlackBody = false;
   const bool Isotropic = true;
   const double Times[2] = {-1.e50, 1.e50};
   std::vector<double> Values(NumValues);
   for (int t = 0; t < NumTimes; t++)
   {
      for (int g = 0; g < ngroups; g++)
      {
         Values[t * ngroups + g] = values_in[g];
      }
   }

   teton_addprofile_internal(&NumTimes,
                             &NumValues,
                             Multiplier,
                             &BlackBody,
                             &Isotropic,
                             Times,
                             Values.data(),
                             TetonProfileID);
}
// Old API for backward compatibility:
inline void teton_addprofile(const int *NumTimes,
                             const int *NumValues,
                             const int *NumInterpValues, //ignored
                             const double *Multiplier,
                             const bool *BlackBody,
                             const bool *Isotropic,
                             const double *Times,
                             const double *Values)
{
   int dummyInt = -1;
   (void) NumInterpValues;
   teton_addprofile_internal(NumTimes, NumValues, Multiplier, BlackBody, Isotropic, Times, Values, &dummyInt);
}

// Old API for backward compatibility:
inline void teton_constructitercontrols(
   const int *noutmx,    // scalar, maximum number of outer (temperature) iterations before a linear solve
   const int *ninmx,     // scalar, maximum number of inner iterations per coupling update [not used currently]
   const int *ngdamx,    // scalar, maximum number of gray steps before a new temperature iterate
   const double *epstmp, // scalar, temperature solve tolerance
   const double *epsinr, // scalar, linear solve tolerance (should be the smallest)
   const double *epsgda) // scalar, gray diffusion-like tolerance (between epstmp and epsinr)
{
   (void) ninmx;
   // This will be really, really chatty.
   std::cout
      << "WARNING: Calling teton_constructitercontrols with tolerance arguments is deprecated.  Please delete arguments and set tolerances via teton_adjust* functions.\n";
   teton_constructitercontrols();

   teton_adjusttemperaturegoals(noutmx, epstmp);
   teton_adjustradenergydensitygoal(epsinr);
   // The old interface that this supplies took sweeps.
   // The new interface takes iterations, so we convert here
   // so that we don't have to explain it in user manuals.
   int greyIters = (*ngdamx - 1) / 2;
   teton_adjustgreygoals(&greyIters, epsgda);
}

//
//  teton_resetprofile
//
//  reset an existing time dependent (step function with numtimes steps) boundary source in Teton
//
//  Profile must have already been added via teton_addprofile
inline void teton_resetprofile(
   const int *TetonProfileID, // scalar, use the returned value of TetonProfileID from teton_addprofile
   const int *NumTimes,       // scalar, number of differnt times for the source
   const int *NumValues,      // scalar, ==numTimes if blackbody BC.  ==numTimes*ngroups if not blackbody BC
   const double *Multiplier,  // scalar, multiplier in front of all emission terms
   const double *Times,       // array[NumTimes], physical times source changes given in ascending order
   const double *Values)      // array[NumValues], group dependent source strength values
{
   teton_resetprofile_internal(TetonProfileID, NumTimes, NumValues, Multiplier, Times, Values);
}
// Alternative API for updating a constant-in-time blackbody temperature BC
inline void teton_resetprofile(
   const int *TetonProfileID, // scalar, use the returned value of TetonProfileID from teton_addprofile
   const double Value)        // scalar, blackbody source strength value
{
   const int NumTimes = 2;
   const int NumValues = 2;
   const double Multiplier = 1.;
   const double Times[2] = {-1.e50, 1.e50};
   const double Values[2] = {Value, Value};
   teton_resetprofile_internal(TetonProfileID, &NumTimes, &NumValues, &Multiplier, Times, Values);
}
// Alternative API for updating a constant-in-time frequency-dependent isotropic source
inline void teton_resetprofile(
   const int *TetonProfileID,     // scalar, use the returned value of TetonProfileID from teton_addprofile
   const double *Multiplier,      // scalar, multiplier in front of all emission terms
   std::vector<double> values_in) // vector of length ngroups, source strength value
{
   const int NumTimes = 2;
   const int ngroups = values_in.size();
   const int NumValues = NumTimes * ngroups;
   const double Times[2] = {-1.e50, 1.e50};
   std::vector<double> Values(NumValues);
   for (int t = 0; t < NumTimes; t++)
   {
      for (int g = 0; g < ngroups; g++)
      {
         Values[t * ngroups + g] = values_in[g];
      }
   }

   teton_resetprofile_internal(TetonProfileID, &NumTimes, &NumValues, Multiplier, Times, Values.data());
}

// enums for flags for signaling to/from Teton
enum class tetonConvControlFlag : int
{
   FluxIter = 11,
   TempIter = 12,
   Invalid = 0
};

enum class tetonDtControlFlag : int
{
   RadTemp = 21,
   ElecTemp = 22,
   Compton = 23,
   SlowConv = 24,
   NoConv = 25,
   MinDt = 26,
   MaxDt = 27,
   None = 28,
   Invalid = 0
};

enum class tetonBcTypeFlag : int
{
   none = 31,
   refl = 32,
   shared = 33,
   temp = 34,
   vac = 35,
   fds = 36, // used repeatedly in BC.profile.typeIs calls within driver
   invalid = 0
};

enum class tetonGeometryFlag : int
{
   slab = 41,
   sphere = 42,
   cylinder = 43,
   rz = 44,
   xyz = 45,
   invalid = 0
};

enum class tetonBcShapeFlag : int
{
   none = 51, // aka shapeless
   iso = 52,
   fds = 53, // referenced as a shape in TetonProfile comments
   invalid = 0
};

enum class tetonComptonFlag : int
{
   none = 61,
   fp = 62,
   boltzmann = 63,
   thomson = 64,
   external_model = 65,
   invalid = 0
};

enum class tetonIterControlsFlag : int
{
   outerTemp = 71,
   outerPhi = 72,
   grey = 73,
   innernl = 74,
   incidentflux = 75,
   global = 79,
   invalid = 0
};

inline tetonBcShapeFlag translateShape(const std::string &_shape)
{
   tetonBcShapeFlag val = tetonBcShapeFlag::none;
   if (_shape == "iso")
   {
      val = tetonBcShapeFlag::iso;
   }
   else if (_shape == "fds")
   {
      val = tetonBcShapeFlag::fds;
   }
   else
   {
      val = tetonBcShapeFlag::invalid;
   }

   return val;
}

inline tetonBcTypeFlag translateType(const std::string &_shape)
{
   tetonBcTypeFlag val = tetonBcTypeFlag::none;
   if (_shape == "refl")
   {
      val = tetonBcTypeFlag::refl;
   }
   else if (_shape == "shared")
   {
      val = tetonBcTypeFlag::shared;
   }
   else if (_shape == "temp")
   {
      val = tetonBcTypeFlag::temp;
   }
   else if (_shape == "vac")
   {
      val = tetonBcTypeFlag::vac;
   }
   else if (_shape == "fds")
   {
      val = tetonBcTypeFlag::fds;
   }
   else
   {
      val = tetonBcTypeFlag::invalid;
   }

   return val;
}

inline double teton_get_default_outer_temp_reltol()
{
   double var{};
   teton_get_default_outer_temp_reltol_internal(&var);
   return var;
}

inline double teton_get_default_outer_energy_density_reltol()
{
   double var{};
   teton_get_default_outer_intensity_reltol_internal(&var);
   return var;
}

inline double teton_get_default_grey_reltol()
{
   double var{};
   teton_get_default_grey_reltol_internal(&var);
   return var;
}

inline double teton_get_default_incident_flux_reltol()
{
   double var{};
   teton_get_default_incident_flux_reltol_internal(&var);
   return var;
}

inline double teton_get_default_inner_nl_reltol()
{
   double var{};
   teton_get_default_inner_nl_reltol_internal(&var);
   return var;
}

inline int teton_get_default_outer_max_it()
{
   int var{};
   teton_get_default_outer_max_it_internal(&var);
   return var;
}

inline int teton_get_default_grey_max_it()
{
   int var{};
   teton_get_default_grey_max_it_internal(&var);
   return var;
}

inline int teton_get_default_incident_flux_max_it()
{
   int var{};
   teton_get_default_incident_flux_max_it_internal(&var);
   return var;
}

inline int teton_get_default_inner_nl_max_it()
{
   int var{};
   teton_get_default_inner_nl_max_it_internal(&var);
   return var;
}

// Set all tolerances together with one change.  We'll set
// all the tolerances relative to the outer tolerance
inline void teton_adjust_relative_tolerance(const double *reltol)
{
   double localRelTol = std::max(1.0e-14, *reltol);
   // This is the main tolerance that everything else is scaled by.
   teton_adjust_temperature_reltol(&localRelTol);
   teton_adjust_radenergydensity_reltol(&localRelTol);
   double localGreyRelTol = *reltol * teton_get_default_grey_reltol() / teton_get_default_outer_temp_reltol();
   double localFluxRelTol = *reltol * teton_get_default_incident_flux_reltol() / teton_get_default_outer_temp_reltol();
   double localNLRelTol = *reltol * teton_get_default_inner_nl_reltol() / teton_get_default_outer_temp_reltol();
   // These will floor the value to 1e-14, so even if relTol is that
   // small, these won't be smaller.
   teton_adjust_fluxexchange_reltol(&localFluxRelTol);
   teton_adjust_grey_reltol(&localGreyRelTol);
   teton_adjust_nonlinear_reltol(&localNLRelTol);
}

// returns 1 if bad
inline int teton_adjust_control(const tetonIterControlsFlag control, const double reltol, const int maxits)
{
   switch (control)
   {
      case tetonIterControlsFlag::global:
         if (reltol > 0.)
         {
            teton_adjust_relative_tolerance(&reltol);
         }
         if (maxits > 0)
         {
            std::cout << "WARNING: There is no global maxits control.  You must adjust the individual maxit controls."
                      << std::endl;
            return -1;
         }
         break;
      case tetonIterControlsFlag::outerTemp:
         if (reltol > 0.)
         {
            teton_adjust_temperature_reltol(&reltol);
         }
         if (maxits > 0)
         {
            teton_adjust_temperature_maxits(&maxits);
         }
         break;
      case tetonIterControlsFlag::outerPhi:
         if (reltol > 0.)
         {
            teton_adjust_radenergydensity_reltol(&reltol);
         }
         if (maxits > 0)
         {
            std::cout
               << "WARNING: There is no maxits control for phi.  Use the outer temperature control to set a max # of iterations for the outer iteration."
               << std::endl;
            return -1;
         }
         break;
      case tetonIterControlsFlag::grey:
         if (reltol > 0.)
         {
            teton_adjust_grey_reltol(&reltol);
         }
         if (maxits > 0)
         {
            teton_adjust_grey_maxits(&maxits);
         }
         break;
      case tetonIterControlsFlag::innernl:
         if (reltol > 0.)
         {
            teton_adjust_nonlinear_reltol(&reltol);
         }
         if (maxits > 0)
         {
            teton_adjust_nonlinear_maxits(&maxits);
         }
         break;
      case tetonIterControlsFlag::incidentflux:
         if (reltol > 0.)
         {
            teton_adjust_fluxexchange_reltol(&reltol);
         }
         if (maxits > 0)
         {
            teton_adjust_fluxexchange_maxits(&maxits);
         }
         break;
      default:
         std::cout << "WARNING: Bad Teton control flag in teton_adjust_control." << std::endl;
         return 1;
   }

   return 0;
}
inline int teton_adjust_control(const tetonIterControlsFlag control, const int maxits)
{
   return teton_adjust_control(control, -1, maxits);
}
inline int teton_adjust_control(const tetonIterControlsFlag control, const double reltol)
{
   return teton_adjust_control(control, reltol, -1);
}

// Backward compatible function.
inline void teton_getedits(const int *verbose,      // deprecated and ignored, use teton_setverbose
                           int *noutrt,             // scalar, number of thermal [outer] iterations in this cycle
                           int *ninrt,              // scalar, number of inner [transport] iterations in this cycle
                           int *ngdart,             // scalar, number of grey diffusion iterations
                           int *nNLIters,           // scalar, number of nonlinear iterations
                           int *maxNLIters,         // scalar, maximum nonlinear iterations
                           int *TrMaxZone,          // scalar, zone id with maximum T rad
                           int *TeMaxZone,          // scalar, zone id with maximum electron temperature
                           int *TrMaxProcess,       // scalar, MPI process with maximum T rad
                           int *TeMaxProcess,       // scalar, MPI process with maximum electron temperature
                           double *dtused,          // scalar, Teton returns the time step used in cycle just completed
                           double *dtrad,           // scalar, radiation vote for next time step
                           double *TrMax,           // scalar, T rad maximum value
                           double *TeMax,           // scalar, T electron maximum value
                           double *EnergyRadiation, // scalar, energy contained in the radiation field
                           double *PowerIncident,   // scalar, power of photons incident
                           double *PowerEscape,     // scalar, power of photons escaping
                           double *PowerAbsorbed,   // scalar, power of energy absorbed
                           double *PowerEmitted,    // scalar, power of photons emitted
                           double *PowerExtSources, // scalar, power of photons from fixed volumetric sources
                           double *PowerCompton)
{ // scalar, power of Compton scattering photons

   (void) verbose;
   double deltaEMat = 0.0;
   double energyCheck = 0.0;
   teton_getedits(noutrt,
                  ninrt,
                  ngdart,
                  nNLIters,
                  maxNLIters,
                  TrMaxZone,
                  TeMaxZone,
                  TrMaxProcess,
                  TeMaxProcess,
                  dtused,
                  dtrad,
                  TrMax,
                  TeMax,
                  &deltaEMat,
                  EnergyRadiation,
                  PowerIncident,
                  PowerEscape,
                  PowerAbsorbed,
                  PowerEmitted,
                  PowerExtSources,
                  PowerCompton,
                  &energyCheck);
}

// Another interface version:
inline void teton_getedits(int *noutrt,             // scalar, number of thermal [outer] iterations in this cycle
                           int *ninrt,              // scalar, number of inner [transport] iterations in this cycle
                           int *ngdart,             // scalar, number of grey diffusion iterations
                           int *nNLIters,           // scalar, number of nonlinear iterations
                           int *maxNLIters,         // scalar, maximum nonlinear iterations
                           int *TrMaxZone,          // scalar, zone id with maximum T rad
                           int *TeMaxZone,          // scalar, zone id with maximum electron temperature
                           int *TrMaxProcess,       // scalar, MPI process with maximum T rad
                           int *TeMaxProcess,       // scalar, MPI process with maximum electron temperature
                           double *dtused,          // scalar, Teton returns the time step used in cycle just completed
                           double *dtrad,           // scalar, radiation vote for next time step
                           double *TrMax,           // scalar, T rad maximum value
                           double *TeMax,           // scalar, T electron maximum value
                           double *EnergyRadiation, // scalar, energy contained in the radiation field
                           double *PowerIncident,   // scalar, power of photons incident
                           double *PowerEscape,     // scalar, power of photons escaping
                           double *PowerAbsorbed,   // scalar, power of energy absorbed
                           double *PowerEmitted,    // scalar, power of photons emitted
                           double *PowerExtSources, // scalar, power of photons from fixed volumetric sources
                           double *PowerCompton,    // scalar, power of Compton scattering photons
                           double *EnergyCheck)     // scalar, total energy not accounted for this cycle.
{
   double dummyDeltaEMat;
   teton_getedits(noutrt,          // scalar, number of thermal [outer] iterations in this cycle
                  ninrt,           // scalar, number of inner [transport] iterations in this cycle
                  ngdart,          // scalar, number of grey diffusion iterations
                  nNLIters,        // scalar, number of nonlinear iterations
                  maxNLIters,      // scalar, maximum nonlinear iterations
                  TrMaxZone,       // scalar, zone id with maximum T rad
                  TeMaxZone,       // scalar, zone id with maximum electron temperature
                  TrMaxProcess,    // scalar, MPI process with maximum T rad
                  TeMaxProcess,    // scalar, MPI process with maximum electron temperature
                  dtused,          // scalar, Teton returns the time step used in cycle just completed
                  dtrad,           // scalar, radiation vote for next time step
                  TrMax,           // scalar, T rad maximum value
                  TeMax,           // scalar, T electron maximum value
                  &dummyDeltaEMat, // thrown away in old interface
                  EnergyRadiation, // scalar, energy contained in the radiation field
                  PowerIncident,   // scalar, power of photons incident
                  PowerEscape,     // scalar, power of photons escaping
                  PowerAbsorbed,   // scalar, power of energy absorbed
                  PowerEmitted,    // scalar, power of photons emitted
                  PowerExtSources, // scalar, power of photons from fixed volumetric sources
                  PowerCompton,    // scalar, power of Compton scattering photons
                  EnergyCheck);    // scalar, total energy not accounted for this cycle.
}

#endif // __TETON_INTERFACE_HH__
