//--------------------------------------------------------------------------//
// TetonUtilities.hh
//
// This file contains functions that are helpful for implementing the Teton
// Conduit interface.
//
//--------------------------------------------------------------------------//

#ifndef __TETON_UTILITIES_HH__
#define __TETON_UTILITIES_HH__

#include "conduit/conduit_blueprint_mesh_utils.hpp"
#include "conduit/conduit_node.hpp"

#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

namespace Teton
{

namespace utilities
{

/**
 \brief Convert paths to int32 storage.

 \param rank The MPI rank
 \param root The root Conduit node to search for the path keys.
 \param keys The vector of path keys that will be modified if they exist.

 \note This is used to ensure that some fields that conduit produces are
       converted to int32 so other algorithms do not end up promoting int
       types to index_t and messing up Teton's int32 assumptions.
 */
void convert_int32(int rank, conduit::Node &root, const std::vector<std::string> &keys);

/**
 \brief Scans through the Conduit tree and returns the paths that have a dtype.

 \param n The node to search.
 \param dtype The dtype we're looking for.
 \param[out] keys The keys with the dtype we're looking for.
 */
void find_dtype(const conduit::Node &n,
                const conduit::DataType &dtype,
                const std::string &path,
                std::vector<std::string> &paths);

/**
 \brief return a vector of keys with dtype int64.
 \param n The node to search.
 \return A vector of paths with dtype int64
 */
std::vector<std::string> find_int64(const conduit::Node &n);

/**
 \brief Scans through field values and returns true if the values do not have errors.
 \param rank The MPI rank
 \param n The Conduit node containing the field values.
 \
 */
bool scan_field_values(int rank, const conduit::Node &n);

/*!
 * \brief Examines a domain and looks for duplicate points with different ids.
 *
 * \param domainId The id of the domain.
 * \param dom The node that contains the domain.
 * \param coordset The node that contains the coordset.
 * \param[out] info A node to contain any findings.
 *
 * \return True if there are duplicate points; False otherwise.
 */
bool find_local_duplicate_points(int domainId,
                                 const conduit::Node &dom,
                                 const conduit::Node &coordset,
                                 conduit::Node &info);

/**
 \brief Gather strings from all ranks and make sure all ranks get that sorted unique vector of strings.

 \param vec The vector of strings on the current MPI rank.
 \param comm The MPI communicator to use.

 \return A vector of strings that includes all of the unique strings across all ranks.
 */
std::vector<std::string> globalizeStringVector(const std::vector<std::string> &vec, MPI_Comm comm);

/**
 \brief Examine blueprint node fields (for provided names, if given) and invoke a function
        on all of the fields that look like they need to be described by an mcarray.

 \param blueprint The blueprint node that contains the topologies and fields.
 \param mainTopologyName The name of the main topology.
 \param options   The options node that contains Teton options.
 \param fieldNames A vector of field names that we want to check. If this is empty, all
                   fields are checked.
 \param func       The function to be invoked when a field looks like it needs to be an
                   mcarray.
*/
template <typename Func>
void iterate_mcarray_candidates(const conduit::Node &blueprint,
                                const std::string &mainTopologyName,
                                const conduit::Node &options,
                                const std::vector<std::string> &fieldNames,
                                Func &&func)
{
   const conduit::Node &fields = blueprint.fetch_existing("fields");
   const conduit::Node &main_topo = blueprint.fetch_existing("topologies/" + mainTopologyName);

   conduit::index_t nzones = conduit::blueprint::mesh::utils::topology::length(main_topo);
   conduit::index_t ngroups = options.fetch_existing("quadrature/num_groups").to_index_t();

   auto looks_like_mcarray = [&](const conduit::Node &f) -> bool
   {
      if (f["association"].as_string() == "element")
      {
         const conduit::Node &v = f.fetch_existing("values");
         // We look for a "scalar" that acts like an mcarray.
         if (v.number_of_children() == 0)
         {
            conduit::index_t nvalues = v.dtype().number_of_elements();
            conduit::index_t values_per_zone = nvalues / nzones;
            if (values_per_zone == ngroups && ngroups > 1 && v.dtype().is_float64())
            {
               // We appear to have a single buffer that contains ngroups values per zone.
               return true;
            }
         }
      }
      return false;
   };

   if (fieldNames.empty())
   {
      for (conduit::index_t i = 0; i < fields.number_of_children(); i++)
      {
         const conduit::Node &f = fields[i];
         if (looks_like_mcarray(f))
            func(f);
      }
   }
   else
   {
      for (const auto &name : fieldNames)
      {
         if (fields.has_path(name))
         {
            const conduit::Node &f = fields.fetch_existing(name);
            if (looks_like_mcarray(f))
               func(f);
         }
      }
   }
}

/**
 * \brief This class prints a simple banner to the console using RAII pattern.
 */
class Banner
{
  public:
   Banner(MPI_Comm c, const std::string &str);
   ~Banner();

  private:
   void printLine(const std::string s) const;
   void emit(char c, int n) const;

   MPI_Comm comm;
   int rank;
   std::string name;
   static int level;
};

} // namespace utilities

} // namespace Teton
#endif
