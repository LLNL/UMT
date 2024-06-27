#include "TetonUtilities.hh"

#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_blueprint_mesh.hpp"
#include "conduit/conduit_relay_mpi.hpp"

#include <cmath>
#include <iostream>
#include <sstream>

namespace Teton
{

namespace utilities
{

void convert_int32(int rank, conduit::Node &root, const std::vector<std::string> &keys)
{
   for (const auto &path : keys)
   {
      if (root.has_path(path))
      {
         conduit::Node &n = root.fetch_existing(path);
         if (!n.dtype().is_int32())
         {
            // Convert the data to int32 and put it back in node n. I tried move/swap
            // to try and steal ifield's data rather than copying it again but that
            // caused n's name to be blank.
            conduit::Node ifield;
            n.to_int32_array(ifield);
            n.set(ifield);
         }
      }
   }
}

void find_dtype(const conduit::Node &n,
                const conduit::DataType &dtype,
                const std::string &path,
                std::vector<std::string> &paths)
{
   auto concat = [](const std::string &path, const std::string &name)
   {
      if (path.empty())
         return name;
      return path + "/" + name;
   };

   if (n.number_of_children() > 0)
   {
      // Make paths be relative to node n at the top level.
      for (conduit::index_t i = 0; i < n.number_of_children(); i++)
      {
         if (path == "")
            find_dtype(n[i], dtype, n[i].name(), paths);
         else
            find_dtype(n[i], dtype, concat(path, n[i].name()), paths);
      }
   }
   else
   {
      if (n.dtype().id() == dtype.id())
         paths.push_back(path);
   }
}

std::vector<std::string> find_int64(const conduit::Node &n)
{
   std::vector<std::string> paths;
   find_dtype(n, conduit::DataType::int64(), "", paths);
   return paths;
}

bool scan_field_values(int rank, const conduit::Node &n)
{
   bool retval = true;
   if (n.dtype().is_float64())
   {
      int count = 0;
      conduit::float64_array arr = n.value();
      for (conduit::index_t i = 0; i < arr.number_of_elements(); i++)
      {
         if (!std::isfinite(arr[i]))
         {
            std::cout << rank << ":" << n.path() << ": elem[" << i << "] is not a number. " << arr[i] << std::endl;
            retval = false;
            count++;
            if (count > 10)
               break;
         }
      }
   }
   return retval;
}

bool find_local_duplicate_points(int domainId,
                                 const conduit::Node &dom,
                                 const conduit::Node &coordset,
                                 conduit::Node &info)
{
   bool retval = false;
   // Make sure Conduit is new enough.
#if (CONDUIT_VERSION_MAJOR == 0 && CONDUIT_VERSION_MINOR >= 9) || (CONDUIT_VERSION_MAJOR > 0)
   using conduit::index_t;

   // See whether any of the points in the local domain are duplicated.
   // If a point's query result does not equal its query index then it
   // must have been defined once before.
   conduit::blueprint::mesh::utils::query::PointQuery localPQ(dom);
   const index_t npts = conduit::blueprint::mesh::coordset::length(coordset);
   for (index_t pi = 0; pi < npts; pi++)
   {
      auto pt = conduit::blueprint::mesh::utils::coordset::_explicit::coords(coordset, pi);
      double pt3[3];
      pt3[0] = pt[0];
      pt3[1] = (pt.size() > 1) ? pt[1] : 0.;
      pt3[2] = (pt.size() > 2) ? pt[2] : 0.;
      localPQ.add(domainId, pt3);
   }
   localPQ.execute(coordset.name());
   for (index_t pi = 0; pi < npts; pi++)
   {
      const auto &res = localPQ.results(static_cast<int>(domainId));
      if (res[pi] != pi)
      {
         const auto pts = localPQ.inputs(domainId);
         double pt3[3]{pts[3 * pi], pts[3 * pi + 1], pts[3 * pi + 2]};
         std::stringstream ss;
         ss << "Domain " << domainId << " duplicated point " << pi << " (" << pt3[0] << ", " << pt3[1] << ", " << pt3[2]
            << ") at " << res[pi] << ".";

         conduit::Node &vn = info.append();
         vn["message"].set(ss.str());
         vn["vertex"] = pi;
         vn["duplicate_vertex"] = res[pi];
         vn["coordinate"].set(pt3, 3);

         retval = true;
      }
   }
#endif
   return retval;
}

std::vector<std::string> globalizeStringVector(const std::vector<std::string> &vec, MPI_Comm comm)
{
   // Make a Conduit node from it.
   conduit::Node send_node;
   for (const auto &value : vec)
      send_node[value] = 1;

   // Send the data to all ranks.
   conduit::Node recv_node;
   conduit::relay::mpi::all_gather_using_schema(send_node, recv_node, comm);

   // Pick through the output and make the output string vector from the node names.
   std::set<std::string> unique;
   for (conduit::index_t i = 0; i < recv_node.number_of_children(); i++)
   {
      const conduit::Node &child = recv_node[i];
      for (conduit::index_t j = 0; j < child.number_of_children(); j++)
         unique.insert(child[j].name());
   }
   std::vector<std::string> retval;
   retval.insert(retval.begin(), unique.begin(), unique.end());

   return retval;
}

//------------------------------------------------------------------------------
int Banner::level = 0;

Banner::Banner(MPI_Comm c, const std::string &str) : comm(c), rank(0), name(str)
{
   MPI_Comm_rank(comm, &rank);
   MPI_Barrier(comm);
   if (rank == 0)
      printLine(name + " (start)");
   MPI_Barrier(comm);
   level++;
}

Banner::~Banner()
{
   level--;
   MPI_Barrier(comm);
   if (rank == 0)
      printLine(name + " (end)");
   MPI_Barrier(comm);
}

void Banner::printLine(const std::string s) const
{
   int n = std::max(2, (80 - 2 - static_cast<int>(s.size())) / 2);
   emit(' ', level * 2);
   emit('=', n);
   std::cout << " " << s << " ";
   emit('=', n);
   std::cout << std::endl;
}

void Banner::emit(char c, int n) const
{
   for (int i = 0; i < n; i++)
      std::cout << c;
}

} // namespace utilities

} // namespace Teton
