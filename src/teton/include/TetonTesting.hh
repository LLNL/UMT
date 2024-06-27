//--------------------------------------------------------------------------//
// TetonTesting.hh
//
// This file contains functions that are helpful for testing teton by
// comparing baselines vs current Conduit output.
//
//--------------------------------------------------------------------------//

#ifndef __TETON_TESTING_HH__
#define __TETON_TESTING_HH__

#include "conduit/conduit_node.hpp"

#include <string>

#include <mpi.h>

namespace Teton
{

namespace testing
{

/*!
 * \brief Compare current node against baseline.
 *
 * \param n        The node that holds the data to be tested.
 * \param fileBase The base name of the file to be used for baselines or current data.
 * \param cycle    The current cycle number.
 * \param make     True if we are making a baseline, false if we're comparing to a baseline.
 * \param comm     The MPI communicator to use.
 *
 * \return True if the current results are sufficiently close to baselines.
 *
 * There are some environment variables that affect its operation.
 *   TETON_TESTING_BASELINE_DIR - set path for baseline files.
 *   TETON_TESTING_CURRENT_DIR - set path for current files.
 */
bool test(const conduit::Node &n, const std::string &fileBase, int cycle, bool make, MPI_Comm comm);

} // namespace testing

} // namespace Teton
#endif
