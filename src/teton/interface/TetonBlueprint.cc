
#include "TetonBlueprint.hh"
#include "conduit/conduit.hpp"
#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_blueprint_mesh_utils.hpp"
#include "conduit/conduit_blueprint_mpi.hpp"
#include "conduit/conduit_blueprint_o2mrelation.hpp"
#include "conduit/conduit_blueprint_o2mrelation_iterator.hpp"
#include "conduit/conduit_relay.hpp"
#include "conduit/conduit_relay_mpi_io_blueprint.hpp"
#include "dbc_macros.h"

#if defined(TETON_ENABLE_CALIPER) && defined(TETON_ENABLE_BLUEPRINT_GENERATION_CALIPER_TIMERS)
#include "caliper/cali.h"
#else
#define CALI_MARK_BEGIN(label)
#define CALI_MARK_END(label)
#define CALI_CXX_MARK_SCOPE(name)
#define CALI_CXX_MARK_FUNCTION
#endif

#include <algorithm> // for sort
#include <limits>
#include <map>
#include <set>
#include <vector>

#include <string>

using std::cout;
using std::endl;
//using namespace mfem;
// access one-to-many index types
namespace O2MIndex = conduit::blueprint::o2mrelation;

//---------------------------------------------------------------------------
// Utilities
//---------------------------------------------------------------------------
/**
 @brief This class takes a list of ids and makes a unique index from the
        values. The index can then be added as a pair (index,id) that 
        can be looked up later. This is used to do things like build a
        set of unique faces and then look up the id of a face from points
        later to see if it exists in the set.
 */
template <typename T> class unique
{
  public:
   using IdType = T;
   using IndexType = std::uint64_t;

   void beginConstruction(std::int64_t size)
   {
      buffer.reserve(16);
      buffer.clear();
      index_to_id.clear();
      index_to_id.reserve(size);
   }

   IndexType makeIndex(const IdType *ids, size_t n)
   {
      // Copy the ids into the buffer and sort the buffer
      if (buffer.capacity() < n)
         buffer.reserve(n);
      buffer.clear();
      for (size_t i = 0; i < n; i++)
         buffer.push_back(ids[i]);
      std::sort(buffer.begin(), buffer.end());

      return hashBuffer();
   }

   void add(IndexType index, IdType id)
   {
      index_to_id.push_back(std::make_pair(index, id));
   }

   void endConstruction()
   {
      // Sort by the index component.
      std::sort(index_to_id.begin(),
                index_to_id.end(),
                [](const std::pair<IndexType, IdType> &lhs, const std::pair<IndexType, IdType> &rhs)
                { return lhs.first < rhs.first; });
   }

   bool find(IndexType index, IdType &id) const
   {
      std::int64_t i = -1;
      std::int64_t left = 0;
      std::int64_t right = index_to_id.size() - 1;
      while (left <= right)
      {
         std::int64_t m = (left + right) / 2;
         if (index_to_id[m].first < index)
            left = m + 1;
         else if (index_to_id[m].first > index)
            right = m - 1;
         else
         {
            i = m;
            break;
         }
      }
      bool retval = (i != -1) && (index_to_id[i].first == index);
      id = retval ? index_to_id[i].second : 0;
      return retval;
   }

  private:
   IndexType hashBuffer() const
   {
// NOTE: When a new Conduit is out, use conduit::utils::hash.
#ifdef CONDUIT_HAS_INDEX_T_HASH
      return conduit::utils::hash(buffer.data(), buffer.size());
#else
      return hash_uint8(reinterpret_cast<const std::uint8_t *>(buffer.data()), buffer.size() * sizeof(IdType));
#endif
   }

#ifndef CONDUIT_HAS_INDEX_T_HASH
   // TODO: use conduit::utils::hash when it becomes available.
   std::uint64_t hash_uint8(const std::uint8_t *data, size_t n) const
   {
      std::uint32_t hashF = 0;

      // Build the length into the hash so {1} and {0,1} hash to different values.
      const auto ldata = reinterpret_cast<const std::uint8_t *>(&n);
      for (size_t e = 0; e < sizeof(n); e++)
      {
         hashF += ldata[e];
         hashF += hashF << 10;
         hashF ^= hashF >> 6;
      }
      // hash the data forward and backwards.
      std::uint32_t hashB = hashF;
      for (size_t i = 0; i < n; i++)
      {
         hashF += data[i];
         hashF += hashF << 10;
         hashF ^= hashF >> 6;

         hashB += data[n - 1 - i];
         hashB += hashB << 10;
         hashB ^= hashB >> 6;
      }
      hashF += hashF << 3;
      hashF ^= hashF >> 11;
      hashF += hashF << 15;

      hashB += hashB << 3;
      hashB ^= hashB >> 11;
      hashB += hashB << 15;

      // Combine the forward, backward into a uint64_t.
      return (static_cast<std::uint64_t>(hashF) << 32) | (static_cast<std::uint64_t>(hashB));
   }
#endif

   std::vector<IdType> buffer;
   std::vector<std::pair<IndexType, IdType>> index_to_id;
};

//---------------------------------------------------------------------------
/**
 @brief This class copies data from one Conduit node to another, using src_dest
        as a map that determines how elements are copied. If the source is -1
        then a fill value is written into dest. Otherwise the copy happens like
        this: dest[src_dest.second] = src[src_dest.first]. Data are copied
        without using Conduit-style casting for a subset of fast-paths.
 */
class copier
{
  public:
   // Copy src to dest using a src_dest as a guide.
   void operator()(const conduit::Node &src,
                   conduit::Node &dest,
                   const std::vector<std::pair<int, int>> &src_dest,
                   double value)
   {
      bool handled = false;
      // Handle common cases using concrete types.
      if (src.dtype().id() == dest.dtype().id())
      {
         if (src.dtype().is_float())
            handled = copy_type(src.as_float_ptr(), dest.as_float_ptr(), src_dest, static_cast<float>(value));
         else if (src.dtype().is_double())
            handled = copy_type(src.as_double_ptr(), dest.as_double_ptr(), src_dest, value);
         else if (src.dtype().is_int())
            handled = copy_type(src.as_int_ptr(), dest.as_int_ptr(), src_dest, static_cast<int>(value));
         else if (src.dtype().is_int64())
            handled = copy_type(src.as_int64_ptr(), dest.as_int64_ptr(), src_dest, static_cast<conduit::int64>(value));
      }
      else if (src.dtype().is_float() && dest.dtype().is_double())
         handled = copy_type(src.as_float_ptr(), dest.as_double_ptr(), src_dest, static_cast<double>(value));
      else if (src.dtype().is_double() && dest.dtype().is_float())
         handled = copy_type(src.as_double_ptr(), dest.as_float_ptr(), src_dest, static_cast<float>(value));

      // Do the Conduit cast style as a backup case.
      if (!handled)
      {
         size_t n = src_dest.size();
         conduit::Node src_data, dest_data;
         for (size_t i = 0; i < n; i++)
         {
            // Set up src.
            src_data.reset();
            if (src_dest[i].first >= 0)
            {
               src_data.set_external(conduit::DataType(src.dtype().id(), 1),
                                     (void *) src.element_ptr(src_dest[i].first));
            }
            else
            {
               src_data.set(value);
            }
            // Copy src into dest.
            dest_data.set_external(conduit::DataType(dest.dtype().id(), 1),
                                   (void *) dest.element_ptr(src_dest[i].second));
            src_data.to_data_type(dest_data.dtype().id(), dest_data);
         }
      }
   }

  private:
   template <typename Tsrc, typename Tdest>
   bool copy_type(const Tsrc *src, Tdest *dest, const std::vector<std::pair<int, int>> &src_dest, Tdest value)
   {
      size_t n = src_dest.size();
      for (size_t i = 0; i < n; i++)
         dest[src_dest[i].second] = (src_dest[i].first >= 0) ? src[src_dest[i].first] : value;
      return true;
   }
};

//---------------------------------------------------------------------------
/**
 @brief This function iterates over elements in a topology (we use it for faces)
        and calls a functor or lambda function on the current element.

 @param topology The topology to use for iteration.
 @param topo_length The topology length - we pass it in because it can be 
                    calculated externally and used for other things such as
                    sizing arrays.
 @param func The operation to be applied to the current element.
 */
template <typename Func> void iterate_topology(const conduit::Node &topology, int topo_length, Func &&func)
{
   const int *conn = topology.fetch_existing("elements/connectivity").value();
   const int *points = conn;

   if (topology.has_path("elements/sizes"))
   {
      conduit::int_accessor sizes = topology.fetch_existing("elements/sizes").value();
      for (int f = 0; f < topo_length; ++f)
      {
         int size = sizes[f];
         func(f, points, size);
         // Move to the next element.
         points += size;
      }
   }
   else
   {
      // No sizes -- all shapes are the same.
      conduit::blueprint::mesh::utils::ShapeType shape(topology);
      int size = static_cast<int>(shape.indices);
      for (int f = 0; f < topo_length; ++f)
      {
         func(f, points, size);
         // Move to the next element.
         points += size;
      }
   }
}

//---------------------------------------------------------------------------
/**
 @brief This class helps compute the point for a point in the topology. We
        use this as a replacement for Conduit's unstructured::topology so
        we can get more performance by maintaining some internal state and
        just computing the point that we want. Using the internal state
        again and again saves us from having to do a bunch of Conduit
        node lookups, whose costs add up.
 */
class CornerPointHelper
{
  public:
   CornerPointHelper(const conduit::Node &n)
      : topo_shape(n),
        piter(),
        eiter(),
        pidx_node(nullptr),
        eidx_node(nullptr),
        is_polygonal(false)
   {
      is_polygonal = topo_shape.is_polygonal();
      if (is_polygonal)
      {
         const conduit::Node &enode = n["elements"];
         piter = conduit::blueprint::o2mrelation::O2MIterator(enode);
         pidx_node = enode.fetch_ptr("connectivity");
      }
      else
      {
         const conduit::Node &enode = n["subelements"];
         piter = conduit::blueprint::o2mrelation::O2MIterator(enode);
         eiter = conduit::blueprint::o2mrelation::O2MIterator(n["elements"]);
         pidx_node = enode.fetch_ptr("connectivity");
         eidx_node = n.fetch_ptr("elements/connectivity");
      }
   }

   conduit::index_t corner_point(const conduit::index_t ei)
   {
      return is_polygonal ? corner_point_polygonal(ei) : corner_point_polyhedral(ei);
   }

  private:
   conduit::index_t corner_point_polygonal(const conduit::index_t ei)
   {
      conduit::index_t retval = 0;
      conduit::index_t_accessor pidxs_vals = pidx_node->value();

      // Look at the points in polygon ei.
      piter.to(ei, O2MIndex::ONE);
      piter.to_front(O2MIndex::MANY);
      bool compare = false;
      while (piter.has_next(O2MIndex::MANY))
      {
         piter.next(O2MIndex::MANY);
         const conduit::index_t ptidx = pidxs_vals[piter.index(O2MIndex::DATA)];
         if (compare)
            retval = std::min(retval, ptidx);
         else
         {
            retval = ptidx;
            compare = true;
         }
      }
      return retval;
   }

   conduit::index_t corner_point_polyhedral(const conduit::index_t ei)
   {
      conduit::index_t retval = 0;
      conduit::index_t_accessor pidxs_vals = pidx_node->value();
      conduit::index_t_accessor eidxs_vals = eidx_node->value();

      // Iterates the faces for the ei'th PH element. These should be unique.
      eiter.to(ei, O2MIndex::ONE);
      eiter.to_front(O2MIndex::MANY);
      bool compare = false;
      while (eiter.has_next(O2MIndex::MANY))
      {
         eiter.next(O2MIndex::MANY);
         const conduit::index_t faceidx = eidxs_vals[eiter.index(O2MIndex::DATA)];

         // Look at the points in the face.
         piter.to(faceidx, O2MIndex::ONE);
         piter.to_front(O2MIndex::MANY);
         while (piter.has_next(O2MIndex::MANY))
         {
            piter.next(O2MIndex::MANY);
            const conduit::index_t ptidx = pidxs_vals[piter.index(O2MIndex::DATA)];
            if (compare)
               retval = std::min(retval, ptidx);
            else
            {
               retval = ptidx;
               compare = true;
            }
         }
      }
      return retval;
   }

  private:
   conduit::blueprint::mesh::utils::ShapeType topo_shape;
   conduit::blueprint::o2mrelation::O2MIterator piter, eiter;
   const conduit::Node *pidx_node;
   const conduit::Node *eidx_node;
   bool is_polygonal;
};
//---------------------------------------------------------------------------

void TetonBlueprint::verifyInput(conduit::Node &meshNode, MPI_Comm comm)
{
   conduit::Node info;
   bool passed;

   int par_rank, par_size;
   MPI_Comm_rank(comm, &par_rank);
   MPI_Comm_size(comm, &par_size);

   passed = conduit::blueprint::mpi::verify("mesh", meshNode, info, comm);
   if (!passed)
   {
      info.print();
      TETON_FATAL_C(par_rank, "Blueprint mesh input failed verify.");
   }
}

void TetonBlueprint::GenerateElementOffsets(conduit::Node &topo, const std::string &key)
{
   // Generate the offsets only if they do not already exist. This safeguards
   // against the offset type changing for polyhedral topologies in certain versions
   // of Conduit.
   if (!topo.has_path(key))
      conduit::blueprint::mesh::topology::unstructured::generate_offsets(topo, topo[key]);
}

void TetonBlueprint::ComputeCornerOffsets(int nzones)
{
   corner_offsets.resize(nzones);
   int offset = 0;
   for (int zone = 0; zone < nzones; ++zone)
   {
      corner_offsets[zone] = offset;
      offset += zone_to_corners2[zone].size();
   }
}

void TetonBlueprint::ComputeLocalCornersInZone(int nzones, int ncorners)
{
   corner_to_lcorner.resize(ncorners);
   //int offset = 0;
   for (int zone = 0; zone < nzones; ++zone)
   {
      size_t ncorners_local = zone_to_corners2[zone].size();
      for (size_t c = 0; c < ncorners_local; ++c)
      {
         int corner = zone_to_corners2[zone][c];
         corner_to_lcorner[corner] = c;
      }
   }
}

int TetonBlueprint::ComputeMaxCorners(const conduit::Node &topo) const
{
   const int ndim = conduit::blueprint::mesh::utils::topology::dims(topo);
   int nc;

   if (ndim == 3)
   {
      // Assume 8 by default.
      nc = 8;
      if (topo.has_path("elements/shape"))
      {
         std::string shape(topo["elements/shape"].as_string());
         if (shape == "tet")
         {
            nc = 4;
         }
         else if (shape == "wedge")
         {
            nc = 6;
         }
         else if (shape == "pyramid")
         {
            nc = 5;
         }
         else if (shape == "hex")
         {
            nc = 8;
         }
         else if (shape == "polyhedral")
         {
            int maxZoneCorners = 0;

            conduit::int_accessor zoneFaces = topo.fetch_existing("elements/connectivity").value();
            conduit::int_accessor zoneFacesSize = topo.fetch_existing("elements/sizes").value();
            conduit::int_accessor faceConn = topo.fetch_existing("subelements/connectivity").value();
            conduit::int_accessor faceSizes = topo.fetch_existing("subelements/sizes").value();
            conduit::int_accessor faceOffsets = topo.fetch_existing("subelements/offsets").value();

            // Iterate over all zones and find the max number of corners for
            // each zone. Keep the largest one we find.
            conduit::index_t nzone = zoneFacesSize.number_of_elements();
            conduit::index_t zoneOffset = 0;
            for (conduit::index_t ei = 0; ei < nzone; ei++)
            {
               std::set<int> zoneUniqueVertices;
               conduit::index_t nfaces = zoneFacesSize[ei];
               for (conduit::index_t fi = 0; fi < nfaces; fi++)
               {
                  int faceId = zoneFaces[zoneOffset + fi];
                  int faceOffset = faceOffsets[faceId];
                  int faceSize = faceSizes[faceId];
                  // Insert all of the face vertices into the set.
                  for (int vi = 0; vi < faceSize; vi++)
                     zoneUniqueVertices.insert(faceConn[faceOffset + vi]);
               }
               zoneOffset += nfaces;

               const int nZoneCorners = static_cast<int>(zoneUniqueVertices.size());
               maxZoneCorners = std::max(nZoneCorners, maxZoneCorners);
            }

            if (nzone > 0)
               nc = maxZoneCorners;
         }
      }
   }
   else if (ndim == 2)
   {
      if (topo.has_path("elements/shape") && topo["elements/shape"].as_string() == "polygonal")
      {
         conduit::int_accessor counts = topo.fetch_existing("elements/sizes").value();
#if 1
         // For use with earlier Conduit versions that lack conduit::DataAccessor::max().
         int counts_max = std::numeric_limits<int>::lowest();
         for (conduit::index_t ci = 0; ci < counts.number_of_elements(); ci++)
            counts_max = std::max(counts_max, counts[ci]);
         nc = counts_max;
#else
         // For use with later Conduit versions (re-enable later).
         nc = counts.max();
#endif
      }
      else
      {
         nc = 4;
      }
   }
   else // ndim == 1
   {
      nc = 2;
   }

   return nc;
}

void TetonBlueprint::ComputeZoneFaceToHalfFaceDict()
{
   size_t nzones = zone_to_faces2.size();
   int offset = 0;

   for (size_t zone = 0; zone < nzones; ++zone)
   {
      size_t nfaces_local = zone_to_faces2[zone].size();
      int loffset = 0;
      for (size_t f = 0; f < nfaces_local; ++f)
      {
         int face = zone_to_faces2[zone][f];
         std::pair<int, int> zoneface;
         zoneface.first = zone;
         zoneface.second = face;
         zoneface_to_halfface[zoneface] = offset;
         offset += 1;
         zoneface_to_lface[zoneface] = loffset;
         loffset += 1;
      }
   }
}

int TetonBlueprint::GetOppositeZone(int zone, int face)
{
   size_t nzones_on_face = face_to_zones2[face].size();
   if (nzones_on_face == 1)
      return -1;
   int zone1 = face_to_zones2[face][0];
   int zone2 = face_to_zones2[face][1];
   int zone_opp = -1;
   if (zone1 == zone)
   {
      zone_opp = zone2;
      return zone_opp;
   }
   else if (zone2 == zone)
   {
      zone_opp = zone1;
      return zone_opp;
   }
   else
   {
      return zone_opp;
   }
}

int TetonBlueprint::GetOppositeCorner(int zone, int face, int corner, int rank)
{
   int corner1 = -1;
   size_t ncorners = face_to_corners2[face].size();

   const int dim = mParametersNode.fetch_existing("size/ndim").as_int32();
   double node_x = corner_to_node_x[corner];
   double node_y = corner_to_node_y[corner];
   double node_z = (dim > 2) ? corner_to_node_z[corner] : 0.;

   for (size_t c = 0; c < ncorners; ++c)
   {
      corner1 = face_to_corners2[face][c];
      double node1_x = corner_to_node_x[corner1];
      double node1_y = corner_to_node_y[corner1];
      int zone1 = corner_to_zone[corner1];

      // Compilers frequently throw warnings about floating point comparisons being unreliable.
      // In this case, we should have identical binary comparisons because we are just comparing positions.
      if (dim > 2)
      {
         const double node1_z = corner_to_node_z[corner1];
         if (node1_x == node_x && node1_y == node_y && node1_z == node_z && zone1 != zone)
         {
            return corner1;
         }
      }
      else
      {
         if (node1_x == node_x && node1_y == node_y && zone1 != zone)
         {
            return corner1;
         }
      }
   }

   // TODO: change for new blueprint format
   bool is_shared = face_is_shared_bndry[face];

   if (!is_shared)
   {
      return -1;
   }
   else
   {
      std::pair<int, int> zoneface(zone, face);
      try
      {
         const int nbelem = mParametersNode.fetch_existing("size/nbelem").to_int32();
         int halfface = zoneface_to_halfface.at(zoneface);
         return -(nbelem + halfface + 1);
      }
      catch (const std::exception &e)
      {
         TETON_FATAL_C(rank, e.what())
      }
   }
   return -1;
}

void TetonBlueprint::ComputeConnectivityArrays(int zone,
                                               int &corner_offset,
                                               int &nzone_faces,
                                               int &ncorner_faces_total,
                                               int &nzone_corners,
                                               std::vector<int> &zones_opp,
                                               std::vector<int> &corners_local,
                                               std::vector<int> &corners_opp,
                                               std::vector<int> &ncorners_on_zone_face,
                                               std::vector<int> &zone_faces_to_bcids,
                                               int rank)
{
   corner_offset = corner_offsets[zone];
   nzone_faces = zone_to_faces2[zone].size();
   nzone_corners = zone_to_corners2[zone].size();
   ncorner_faces_total = 0;
   ncorners_on_zone_face.resize(nzone_faces);
   zone_faces_to_bcids.resize(nzone_faces);

   for (int f = 0; f < nzone_faces; ++f)
   {
      int face = zone_to_faces2[zone][f];
      int zone_opp = GetOppositeZone(zone, face);
      // NOTE: need to increment by 1 since Fortran indexing starts from 1
      if (zone_opp != -1)
         zone_opp += 1;
      zones_opp.push_back(zone_opp);
      size_t ncorner_faces = face_to_corners2[face].size();
      zone_faces_to_bcids[f] = m_face_to_bcid[face];

      for (size_t c = 0; c < ncorner_faces; ++c)
      {
         int corner = face_to_corners2[face][c];
         int zone1 = corner_to_zone[corner];
         if (zone1 != zone)
         {
            continue;
         }
         int lcorner = corner_to_lcorner[corner];
         // NOTE: this goes from 0-based indexing to 1-based indexing for Teton
         if (lcorner != -1)
            lcorner += 1;
         corners_local.push_back(lcorner);

         int ocorner = GetOppositeCorner(zone, face, corner, rank);
         // NOTE: this goes from 0-based indexing to 1-based indexing for Teton
         if (ocorner > -1)
            ocorner += 1;
         corners_opp.push_back(ocorner);
         ncorner_faces_total += 1;
         ncorners_on_zone_face[f] += 1;
      }
   }
}

void TetonBlueprint::ComputeSharedFaces(int rank)
{
   //conduit::int32 *shared_faces_ptr;
   size_t nfaces = face_to_zones2.size();
   int ndim = mParametersNode.fetch_existing("size/ndim").as_int32();
   face_is_shared_bndry.resize(nfaces);
   for (size_t f = 0; f < nfaces; ++f)
   {
      face_is_shared_bndry[f] = false;
   }

   // Skip all this if we don't have any shared faces.
   if (!mMeshNode.has_path("shared_boundaries/shared_faces"))
   {
      return;
   }

   int sface_offset = 0;
   conduit::int32 *shared_faces_array_ptr = nullptr;
   int nsfaces = 0;
   int sfaces_array_len = mMeshNode.fetch_existing("shared_boundaries/shared_faces").dtype().number_of_elements();
   if (ndim == 2)
   { // for each shared face, {bcid, zone, face, corner1, corner2} is stored
      // in the array shared_boundaries/shared_faces
      nsfaces = sfaces_array_len / 5;
   }
   else if (ndim == 3)
   { // for each shared face, {bcid, zone, face, corner1} is stored
      // in the array shared_boundaries/shared_faces
      nsfaces = sfaces_array_len / 4;
   }
   else
   {
      // TODO: CHANGE !!!!
      std::cout << "the 1D case is not yet implemented!" << std::endl;
   }
   if (nsfaces > 0)
   {
      shared_faces_array_ptr = mMeshNode.fetch_existing("shared_boundaries/shared_faces").as_int32_ptr();
   }
   // We need to tranform the pair (zone, face) to (zone, half-face).
   // The half-face ID uniquely labels a pair (zone,face) (so each interior
   // mesh face has two half-faces associated with it)
   for (int j = 0; j < nsfaces; ++j)
   {
      int zone = shared_faces_array_ptr[sface_offset + 1];
      int face = shared_faces_array_ptr[sface_offset + 2];
      std::pair<int, int> zoneface(zone, face);
      int halfface = 0;
      try
      {
         halfface = zoneface_to_halfface.at(zoneface);
      }
      catch (const std::exception &e)
      {
         TETON_FATAL_C(rank, e.what())
      }
      shared_faces_array_ptr[sface_offset + 2] = halfface;
      if (ndim == 3)
      {
         sface_offset += 4;
      }
      else
      {
         sface_offset += 5;
      }
      face_is_shared_bndry[face] = true;
   }
}

void TetonBlueprint::ComputeTetonMeshConnectivityArray(int rank)
{
   CALI_CXX_MARK_FUNCTION;

   int nzones = mParametersNode.fetch_existing("size/nzones").as_int32();
   int ncorners = mParametersNode.fetch_existing("size/ncornr").as_int32();
   //int ndim = mParametersNode.fetch_existing("size/ndim").as_int32();

   // process mesh topologies in to a more convenient format
   ComputeLocalCornersInZone(nzones, ncorners);
   ComputeCornerOffsets(nzones);
   ComputeZoneFaceToHalfFaceDict();
   ComputeSharedFaces(rank);

   // Compute a connectivity array that standalone Teton uses.
   // For each zone, the following integers are appended:
   //                            zoneID, corner offset, number of corner faces,
   //                            number of corners, zones opposite of each
   //                            zone face, list of local (to the zone) corner
   //                            IDs associated with each corner face, a list
   //                            of global (for each domain decomposition)
   //                            IDs for each corner face, a list of the number
   //                            of corners for each zone face, and a list of
   //                            boundary condition IDs for each zone face
   //int connect_off_set = 0;
   std::vector<int> teton_connectivity;
   for (int zone = 0; zone < nzones; ++zone)
   {
      std::vector<int> zones_opp;
      std::vector<int> corners_local;
      std::vector<int> corners_opp;
      int corner_offset;
      int nzone_faces;
      int nzone_corners;
      int ncorner_faces_total;
      std::vector<int> ncorners_on_zone_face;
      std::vector<int> zone_faces_to_bcids;

      ComputeConnectivityArrays(zone,
                                corner_offset,
                                nzone_faces,
                                ncorner_faces_total,
                                nzone_corners,
                                zones_opp,
                                corners_local,
                                corners_opp,
                                ncorners_on_zone_face,
                                zone_faces_to_bcids,
                                rank);

      // NOTE: here Teton uses indexing starting from 1 (instead of 0 based)
      int zoneID = zone + 1;
      teton_connectivity.push_back(zoneID);
      teton_connectivity.push_back(corner_offset);
      teton_connectivity.push_back(nzone_faces);
      teton_connectivity.push_back(ncorner_faces_total);
      teton_connectivity.push_back(nzone_corners);
      for (int j = 0; j < nzone_faces; ++j)
      {
         teton_connectivity.push_back(zones_opp[j]);
      }
      for (int j = 0; j < ncorner_faces_total; ++j)
      {
         teton_connectivity.push_back(corners_local[j]);
      }
      for (int j = 0; j < ncorner_faces_total; ++j)
      {
         teton_connectivity.push_back(corners_opp[j]);
      }
      for (int j = 0; j < nzone_faces; ++j)
      {
         teton_connectivity.push_back(ncorners_on_zone_face[j]);
      }
      for (int j = 0; j < nzone_faces; ++j)
      {
         teton_connectivity.push_back(zone_faces_to_bcids[j]);
      }
   }

   mMeshNode["teton/arrays/corner_connectivity"].set(teton_connectivity.data(), teton_connectivity.size());
}

void TetonBlueprint::CreateTetonMeshCornerCoords()
{
   CALI_CXX_MARK_FUNCTION;

   const int nzones = mParametersNode.fetch_existing("size/nzones").as_int32();
   const int ndim = mParametersNode.fetch_existing("size/ndim").as_int32();

   // Make an estimate for the vector size.
   const size_t corners_per_zone[] = {1, 2, 4, 8};
   size_t nnodes_est = static_cast<size_t>(nzones) * corners_per_zone[ndim];
   std::vector<double> zone_verts;
   zone_verts.reserve(nnodes_est * ndim);

   // loop over mesh elements
   for (int zone = 0; zone < nzones; ++zone)
   {
      //int zoneID = zone + 1;
      size_t ncorners = zone_to_corners2[zone].size();
      //      zone_ncorners.push_back(ncorners);
      // loop over corners in the element
      for (size_t c = 0; c < ncorners; ++c)
      {
         int corner = zone_to_corners2[zone][c];
         if (ndim == 2)
         {
            // Note: x,y coordinates are swapped in rz geometry
            zone_verts.push_back(corner_to_node_y[corner]);
            zone_verts.push_back(corner_to_node_x[corner]);
         }
         if (ndim == 3)
         {
            zone_verts.push_back(corner_to_node_x[corner]);
            zone_verts.push_back(corner_to_node_y[corner]);
            zone_verts.push_back(corner_to_node_z[corner]);
         }
      }
   }
   mMeshNode["arrays/zone_verts"].set(zone_verts.data(), zone_verts.size());
}

void TetonBlueprint::CreateConnectivityArrays(conduit::Node &meshNode, MPI_Comm comm)
{
   CALI_CXX_MARK_FUNCTION;

   int rank;
   MPI_Comm_rank(comm, &rank);

   conduit::int32 *zone_to_faces;
   conduit::int32 *zone_to_corners;
   conduit::int32 *face_to_zones;

   // zone-to-vertices
   conduit::Node &base_topology = meshNode.fetch_existing("topologies/main");
   conduit::Node &base_coordset = meshNode.fetch_existing("coordsets/coords");
   GenerateElementOffsets(base_topology, "elements/offsets");
   const int nzones = conduit::blueprint::mesh::utils::topology::length(base_topology);
   const int dim = conduit::blueprint::mesh::utils::topology::dims(base_topology);

   TETON_VERIFY_C(rank, (dim == 2 || dim == 3), "Only 2d or 3d blueprint meshes supported.");

   // Create the zone-to-corner and corner-to-zone connectivity info
   conduit::Node &corner_topology = meshNode["topologies/main_corner"];
   conduit::Node &corner_coords = meshNode["coordsets/coords_centroids"];
   conduit::Node &zone2corner_map = meshNode["relations/zone2corner"];
   conduit::Node &corner2zone_map = meshNode["relations/corner2zone"];

   // Call the correct version of generate_corners for a serial vs parallel mesh.
   if (meshNode.has_path("adjsets/main_adjset"))
   {
      CALI_CXX_MARK_SCOPE("generate_corners 1");

      // Store adjset temporarily in this node, until we convert it to pairwise.
      conduit::Node &corner_temp_adjset = meshNode["adjsets/temp_corner"];

      conduit::blueprint::mpi::mesh::generate_corners(meshNode,
                                                      "main_adjset",
                                                      corner_temp_adjset.name(),
                                                      corner_topology.name(),
                                                      corner_coords.name(),
                                                      zone2corner_map,
                                                      corner2zone_map,
                                                      comm);

      conduit::Node &corner_adjset = meshNode["adjsets/main_corner"];

      // Convert adj set to pairwise.
      conduit::blueprint::mesh::adjset::to_pairwise(corner_temp_adjset, corner_adjset);
      // Now get rid of the temp corner adjset.
      meshNode.remove("adjsets/temp_corner");
   }
   else
   {
      CALI_CXX_MARK_SCOPE("generate_corners 2");
      conduit::blueprint::mesh::topology::unstructured::generate_corners(base_topology,
                                                                         corner_topology,
                                                                         corner_coords,
                                                                         zone2corner_map,
                                                                         corner2zone_map);
   }

   GenerateElementOffsets(corner_topology, "elements/offsets");

   // For creating the Teton connectivity array
   zone_to_corners = zone2corner_map["values"].as_int32_ptr();
   conduit::int32 *zone_to_ncorner = zone2corner_map["sizes"].as_int32_ptr();

   // corner-to-verts
   const int ncorners = conduit::blueprint::mesh::utils::topology::length(corner_topology);
   std::vector<int> corner_to_vertex(ncorners);
   CALI_MARK_BEGIN("Making corner_to_node");

   zone_to_corners2.resize(nzones);
   corner_to_node_x.resize(ncorners);
   corner_to_node_y.resize(ncorners);
   corner_to_node_z.resize(ncorners);
   corner_to_zone.resize(ncorners);

   // Get the axis names being used in the coordset and prepare accessors
   // to get at the data. The code already assumes 2d-3d so this should be ok.
   // If 2d, the base_coordset_zc value will alias the x coordinates and won't
   // be used.
   const auto base_axes = conduit::blueprint::mesh::utils::coordset::axes(base_coordset);
   const conduit::Node &bc_values = base_coordset.fetch_existing("values");
   const auto base_coordset_xc = bc_values.fetch_existing(base_axes[0]).as_double_accessor();
   const auto base_coordset_yc = bc_values.fetch_existing(base_axes[1]).as_double_accessor();
   const auto base_coordset_zc = bc_values.fetch_existing(base_axes[dim == 3 ? 2 : 0]).as_double_accessor();

   const auto zone2corner_data = zone2corner_map["values"].as_index_t_accessor();
   conduit::blueprint::o2mrelation::O2MIterator z2c_iter(zone2corner_map);
   CornerPointHelper cph(corner_topology);
   while (z2c_iter.has_next(O2MIndex::DATA))
   {
      z2c_iter.next(O2MIndex::ONE);
      const conduit::index_t zid = z2c_iter.index(O2MIndex::ONE);

      z2c_iter.to_front(O2MIndex::MANY);
      while (z2c_iter.has_next(O2MIndex::MANY))
      {
         z2c_iter.next(O2MIndex::MANY);
         const conduit::index_t cidx = z2c_iter.index(O2MIndex::DATA);
         const conduit::index_t cid = zone2corner_data[cidx];

         zone_to_corners2[zid].push_back(cid);
         corner_to_zone[cid] = zid;

         // Compute the points here.
         conduit::index_t cvid = cph.corner_point(cid);

         corner_to_vertex[cid] = cvid;
         // Set the x,y,[z] values for the corner.
         corner_to_node_x[cid] = base_coordset_xc[cvid];
         corner_to_node_y[cid] = base_coordset_yc[cvid];
         if (dim == 3)
         {
            corner_to_node_z[cid] = base_coordset_zc[cvid];
         }
      }
   }
   CALI_MARK_END("Making corner_to_node");

   // Create the zone-to-face and face-to-zone connectivity info
   conduit::Node &face_topology = meshNode["topologies/main_face"];
   conduit::Node &zone2face_map = meshNode["relations/zone2face"];
   conduit::Node &face2zone_map = meshNode["relations/face2zone"];

   // Call the correct version of generate_faces/lines corners for a serial vs parallel mesh.
   if (meshNode.has_path("adjsets/main_adjset"))
   {
      CALI_CXX_MARK_SCOPE("Generate faces/lines");

      // Store adjset temporarily in this node, until we convert it to pairwise.
      conduit::Node &face_temp_adjset = meshNode["adjsets/temp_face"];

      if (dim == 3)
      {
         conduit::blueprint::mpi::mesh::generate_faces(meshNode,
                                                       "main_adjset",
                                                       face_temp_adjset.name(),
                                                       face_topology.name(),
                                                       zone2face_map,
                                                       face2zone_map,
                                                       comm);
      }
      else if (dim == 2)
      {
         conduit::blueprint::mpi::mesh::generate_lines(meshNode,
                                                       "main_adjset",
                                                       face_temp_adjset.name(),
                                                       face_topology.name(),
                                                       zone2face_map,
                                                       face2zone_map,
                                                       comm);
      }
      conduit::Node &face_adjset = meshNode["adjsets/main_face"];
      // Convert adj set to pairwise.
      conduit::blueprint::mesh::adjset::to_pairwise(face_temp_adjset, face_adjset);
      // Now get rid of the temp adjset.
      meshNode.remove("adjsets/temp_face");
   }
   else
   {
      CALI_CXX_MARK_SCOPE("Generate faces/lines 2");

      if (dim == 3)
      {
         conduit::blueprint::mesh::topology::unstructured::generate_faces(base_topology,
                                                                          face_topology,
                                                                          zone2face_map,
                                                                          face2zone_map);
      }
      else if (dim == 2)
      {
         conduit::blueprint::mesh::topology::unstructured::generate_lines(base_topology,
                                                                          face_topology,
                                                                          zone2face_map,
                                                                          face2zone_map);
      }
   }

   GenerateElementOffsets(face_topology, "elements/offsets");
   if (mMeshNode.has_path("topologies/boundary"))
   {
      CALI_CXX_MARK_SCOPE("Generate offsets for boundaries");

      conduit::Node &bndry_topology = meshNode["topologies/boundary"];
      GenerateElementOffsets(bndry_topology, "elements/offsets");
   }
   // zone-to-faces
   CALI_MARK_BEGIN("zones_to_faces");
   zone_to_faces = zone2face_map["values"].as_int32_ptr();
   conduit::int32 *zone_to_face_offsets = zone2face_map["offsets"].as_int32_ptr();
   conduit::int32 *zone_to_nfaces = zone2face_map["sizes"].as_int32_ptr();

   zone_to_faces2.resize(nzones);
   for (int zone = 0; zone < nzones; ++zone)
   {
      int nfaces_in_zone = zone_to_nfaces[zone];
      int offset = zone_to_face_offsets[zone];
      for (int f = 0; f < nfaces_in_zone; ++f)
      {
         int fid = zone_to_faces[offset + f];
         zone_to_faces2[zone].push_back(fid);
      }
   }
   CALI_MARK_END("zones_to_faces");

   // face-to-zones
   CALI_MARK_BEGIN("faces_to_zones");
   face_to_zones = face2zone_map["values"].as_int32_ptr();
   conduit::int32 *face_to_zone_offsets = face2zone_map["offsets"].as_int32_ptr();
   conduit::int32 *face_to_nzones = face2zone_map["sizes"].as_int32_ptr();
   int nfaces = face2zone_map["sizes"].dtype().number_of_elements();
   face_to_zones2.resize(nfaces);
   for (int face = 0; face < nfaces; ++face)
   {
      int nzones_in_face = face_to_nzones[face];
      int offset = face_to_zone_offsets[face];
      for (int z = 0; z < nzones_in_face; ++z)
      {
         int zid = face_to_zones[offset + z];
         face_to_zones2[face].push_back(zid);
      }
   }
   CALI_MARK_END("faces_to_zones");

   // face-to-vertices
   CALI_MARK_BEGIN("faces_to_vertices");
   std::vector<std::vector<int>> face_to_vertices(nfaces);
   conduit::int32 *face_to_verts_offsets = face_topology["elements/offsets"].as_int32_ptr();
   conduit::int32 *face_to_verts = face_topology["elements/connectivity"].as_int32_ptr();
   int nface2verts = face_topology["elements/connectivity"].dtype().number_of_elements();
   for (int face = 0; face < nfaces; ++face)
   {
      int offset = face_to_verts_offsets[face];
      int nverts_for_face;
      if (face == (nfaces - 1))
      {
         nverts_for_face = nface2verts - offset;
      }
      else
      {
         int offset2 = face_to_verts_offsets[face + 1];
         nverts_for_face = offset2 - offset;
      }
      for (int v = 0; v < nverts_for_face; ++v)
      {
         int vid = face_to_verts[offset + v];
         face_to_vertices[face].push_back(vid);
      }
   }
   CALI_MARK_END("faces_to_vertices");

   // Create face_to_corners
   // Note: we use that the vertex ID for each face is shared
   // by exactly one corner in each zone that has this face
   // TODO: need to reverse order of corners for one zone of each pair of zones
   // sharing a face
   CALI_MARK_BEGIN("faces_to_corners");
   face_to_corners2.resize(nfaces);
   for (int face = 0; face < nfaces; ++face)
   {
      // number of zones with this face
      size_t nzones_for_face = face_to_zones2[face].size();
      // number of vertices on this face
      size_t nverts_for_face = face_to_vertices[face].size();
      // loop over zones that have this face
      for (size_t z = 0; z < nzones_for_face; ++z)
      {
         int zone = face_to_zones2[face][z];
         size_t ncorners_in_zone = zone_to_corners2[zone].size();
         std::vector<int> zoneface_corners;
         for (size_t v = 0; v < nverts_for_face; ++v)
         {
            int face_vertex = face_to_vertices[face][v];
            for (size_t c = 0; c < ncorners_in_zone; ++c)
            {
               int corner = zone_to_corners2[zone][c];
               int corner_vertex = corner_to_vertex[corner];
               if (corner_vertex == face_vertex)
               {
                  zoneface_corners.push_back(corner);
                  break;
               }
            }
         }
         if (dim == 3)
         {
            // NOTE: Teton wants the corner ordering on the zone face
            // ("half face") to be such that, looking at the face from
            // outside the zone, the ordering makes the vector point
            // *into* the zone from the right-hand rule. We have to
            // reverse the orientation for the first zone sharing this mesh face
            // to make this the case.
            // that shares the face (generically there are two zones
            // sharing the mesh face)
            if (z == 0)
            {
               // reverse orientation of corners on this zone face
               size_t ncorners_on_face = zoneface_corners.size();
               face_to_corners2[face].push_back(zoneface_corners[0]);
               for (size_t c = 1; c < ncorners_on_face; c++)
               {
                  int corner = zoneface_corners[ncorners_on_face - c];
                  face_to_corners2[face].push_back(corner);
               }
            }
            else
            {
               size_t ncorners_on_face = zoneface_corners.size();
               for (size_t c = 0; c < ncorners_on_face; c++)
               {
                  int corner = zoneface_corners[c];
                  face_to_corners2[face].push_back(corner);
               }
            }
         }
         else if (dim == 2)
         {
            if (z != 0)
            {
               // reverse orientation of corners on this zone face
               face_to_corners2[face].push_back(zoneface_corners[1]);
               face_to_corners2[face].push_back(zoneface_corners[0]);
            }
            else
            {
               face_to_corners2[face].push_back(zoneface_corners[0]);
               face_to_corners2[face].push_back(zoneface_corners[1]);
            }
         }
      }
   }
   CALI_MARK_END("faces_to_corners");

   // save for mesh motion and forces
   meshNode["arrays/corner_to_vertex"].set(corner_to_vertex.data(), corner_to_vertex.size());
   meshNode["arrays/zone_to_ncorners"].set(zone_to_ncorner, nzones);
   meshNode["arrays/zone_to_corners"].set(zone_to_corners, ncorners);
}

/**
 @brief This method creates a face_attribute field for the main_face topology
        with values that come from boundary_attribute where a face matches a
        boundary face and 0 for other faces.

        The boundary topology is iterated and for each face, we compute a hash
        of the ids and map them to the id of the boundary face. The face topology
        is then iterated and a hash is made for each of its faces. We look up that
        hash to try and map a face in the boundary topology. If we find a matching
        hash then its id is used to copy from the boundary_attributes to the
        face_attribute field.

 @param meshNode The Conduit node that contaisn the meshes.
 @param rank The MPI rank.
 */
void TetonBlueprint::CreateConduitFaceAttributes(conduit::Node &meshNode, int rank)
{
   CALI_CXX_MARK_FUNCTION;
   // Create face boundary condition field.
   conduit::Node src_data, dst_data;

   const conduit::Node &bndry_topo = meshNode["topologies/boundary"];
   const conduit::Node &bndry_field = meshNode["fields/boundary_attribute"];
   const conduit::Node &bndry_values = bndry_field["values"];

   // TODO: clean up how we do this. At this point face_to_zones2 has been
   // created and so has meshNode["topologies/main_face"] from
   // the call to CreateConnectivityArrays
   size_t nfaces = face_to_zones2.size();
   conduit::Node &face_topology = meshNode["topologies/main_face"];
   conduit::Node &face_field = meshNode["fields/face_attribute"];
   face_field["association"].set("element");
   face_field["topology"].set("main_face");
   face_field["values"].set(conduit::DataType::int32(nfaces));
   conduit::Node &face_values = face_field["values"];

   // Go through the boundary faces and make unique faces.
   unique<int> facemap;
   const int bndry_topo_length = conduit::blueprint::mesh::utils::topology::length(bndry_topo);
   facemap.beginConstruction(bndry_topo_length);
   iterate_topology(bndry_topo,
                    bndry_topo_length,
                    [&facemap](int f, const int *face_points, int nface_points)
                    {
      // Make an id for this face.
      auto index = facemap.makeIndex(face_points, nface_points);

      // Add (index,id) to the map,
      facemap.add(index, f);
   });
   facemap.endConstruction();

   // Now, go through the face_topo and use it and the verts_bndryid_map to
   // make the face attributes.
   const int face_topo_length = conduit::blueprint::mesh::utils::topology::length(face_topology);
   std::vector<std::pair<int, int>> src_dest;
   src_dest.reserve(face_topo_length);
   iterate_topology(face_topology,
                    face_topo_length,
                    [&facemap, &src_dest](int f, const int *face_points, int nface_points)
                    {
      // Make an id for this face.
      auto index = facemap.makeIndex(face_points, nface_points);

      // Look up this face to see if we already know it.
      unique<int>::IdType id;
      bool found = facemap.find(index, id);

      // If the face is a boundary face, use the attribute for the associated
      // boundary. Otherwise, give -1 as the source so the copier will assign
      // a default value.
      src_dest.push_back(std::make_pair(found ? id : -1, f));
   });

   // Copy boundary values to face_values. If the face isn't a boundary face
   // (its src is -1), set its attribute value to 0
   copier c;
   c(bndry_values, face_values, src_dest, 0);
}

void TetonBlueprint::GetSurfaceEditZoneFacesAndCorners(int rank,
                                                       const conduit::Node &surf_face_topo,
                                                       const std::map<std::set<int>, int> &verts_face_map,
                                                       std::vector<int> &surf_edits_loczonefaces,
                                                       std::vector<int> &surf_edits_corners) const
{
   const int *corner_to_vertex = mMeshNode.fetch_existing("arrays/corner_to_vertex").value();

   const int surf_face_topo_length = conduit::blueprint::mesh::utils::topology::length(surf_face_topo);
   for (int f = 0; f < surf_face_topo_length; ++f)
   {
      std::vector<conduit::index_t> face_points = conduit::blueprint::mesh::utils::topology::unstructured::points(
         surf_face_topo,
         f);
      std::set<int> face_vertices(face_points.begin(), face_points.end());
      int face = -1;
      if (verts_face_map.find(face_vertices) != verts_face_map.end())
      {
         face = verts_face_map.at(face_vertices);
      }
      else
      {
         TETON_FATAL_C(
            rank,
            "TetonBlueprint::ProcessSurfaceEdits: surface edit face doesn't have vertices that correspond to a mesh face");
      }

      // This mesh face shares (generically) two zones---need to determine the
      // correct zone based on the vertex orientation of the vertices on the face
      // NOTE: this logic will not work for AMR meshes
      // TODO: fix for AMR meshes
      // The code below uses that face_to_corners2[face] has a list of corners
      // associated with the mesh face, face. Generically this list will involves two
      // zones, and the orientation will be according to the left-hand rule. Using the association between each corner and vertex (for non-AMR meshes), we want to find the zone with the same orientation
      int vert1 = face_points[0];
      int vert2 = face_points[1];
      int zone1 = face_to_zones2[face][0];
      int zone_for_surf_face = zone1;
      if (face_to_zones2[face].size() > 1)
      {
         int zone2 = face_to_zones2[face][1];
         zone_for_surf_face = zone2;
         // Loop over corners associated with first zone that
         // shares this face. For non-AMR meshes, each vertex
         // in the zone corresponds to a corner
         int nverts_for_face = face_to_corners2[face].size() / 2;
         for (int c = 0; c < nverts_for_face - 1; ++c)
         {
            int corner1 = face_to_corners2[face][c];
            int corner2 = face_to_corners2[face][c + 1];
            int vert_tail = corner_to_vertex[corner1];
            int vert_head = corner_to_vertex[corner2];
            if ((vert1 == vert_tail) && (vert2 == vert_head))
            {
               zone_for_surf_face = zone1;
               continue;
            }
         }
         if (nverts_for_face > 2) // true if dim == 3
         {
            int corner1 = face_to_corners2[face][nverts_for_face - 1];
            int corner2 = face_to_corners2[face][0];
            int vert_tail = corner_to_vertex[corner1];
            int vert_head = corner_to_vertex[corner2];
            if ((vert1 == vert_tail) && (vert2 == vert_head))
            {
               zone_for_surf_face = zone1;
            }
         }
      }

      // Determine the local zone face ID associated with
      // the pair (zoneID, faceID)
      int nfaces_in_zone = zone_to_faces2[zone_for_surf_face].size();
      int local_face = -1;
      for (int f1 = 0; f1 < nfaces_in_zone; ++f1)
      {
         int face_in_zone = zone_to_faces2[zone_for_surf_face][f1];
         if (face_in_zone == face)
            local_face = f1;
      }
      if (local_face == -1)
      {
         TETON_FATAL_C(
            rank,
            "TetonBlueprint::ProcessSurfaceEdits: surface edit face doesn't have vertices that correspond to a mesh face");
      }

      // Finally, push the global corner IDs and local zone face IDs
      // associated with the pair (zoneID, faceID)
      int ncorner_faces = face_to_corners2[face].size();
      for (int c = 0; c < ncorner_faces; ++c)
      {
         int corner = face_to_corners2[face][c];
         int zone = corner_to_zone[corner];
         if (zone != zone_for_surf_face)
         {
            continue;
         }
         // Note: increment because Teton uses 1-based indexing
         surf_edits_loczonefaces.push_back(local_face + 1);
         surf_edits_corners.push_back(corner + 1);
      }
   }
}

// WORKING ON //
// TODO: fix for AMR meshes
// NOTE: Compare to template <typename Mesh> void Teton<Mesh>::makeCornerLists
// This function outputs, for the jth surface:
// a list surf_edits_loczonefaces[j] of local zone face IDs (i.e., local to the zone) associated with the jth surface
// a list surf_edits_corners[j] of global corner IDs associated with the jth surface

void TetonBlueprint::ProcessSurfaceEdits(int rank)
{
   // Create the map {face vertices => face}. Here the face IDs are those
   // generated by the call to generate_faces. We need to match these face IDs
   // with the vertices associated with the surface edit faces
   conduit::Node &face_topo = mMeshNode["topologies/main_face"];
   std::map<std::set<int>, int> verts_face_map; // {face vertices => face}
   const int face_topo_length = conduit::blueprint::mesh::utils::topology::length(face_topo);
   for (int f = 0; f < face_topo_length; ++f)
   {
      std::vector<conduit::index_t> face_points = conduit::blueprint::mesh::utils::topology::unstructured::points(
         face_topo,
         f);
      std::set<int> face_vertices(face_points.begin(), face_points.end());
      verts_face_map[face_vertices] = f;
   }

   // NOTE: the call to CreateConnectivityArrays creates and stores
   // arrays/corner_to_vertex, so need this function to be called before
   // the call to ProcessSurfaceEdits
   if (!mMeshNode.has_path("arrays/corner_to_vertex"))
   {
      TETON_FATAL_C(
         rank,
         "TetonBlueprint::ProcessSurfaceEdits: field arrays/corner_to_vertex has not yet been created; need to call CreateConnectivityArrays first");
   }

   // Loop over surface_edit topologies
   const conduit::Node &topos = mMeshNode["topologies"];
   const conduit::Node &surface_edits = mParametersNode["surface_edits"];
   conduit::NodeConstIterator surface_edits_it = surface_edits.children();
   while (surface_edits_it.has_next())
   {
      const conduit::Node &param_surface_edit = surface_edits_it.next();
      std::vector<int> surf_edits_loczonefaces, surf_edits_corners;
      // Check whether the surface topology exists on this rank.
      // Only fill in the vectors if the topology exists on this rank. We still add
      // the surface edit since the consumer of it will pass data into teton_surfaceedit,
      // which contains collective communication and all ranks need to participate.
      std::string surface_edit_name_str = param_surface_edit["zone_face_topology_name"].as_string();
      if (topos.has_child(surface_edit_name_str))
      {
         // get the surface topology
         const conduit::Node &surf_face_topo = topos.fetch_existing(surface_edit_name_str);

         // Get the local zone faces and corners.
         GetSurfaceEditZoneFacesAndCorners(rank,
                                           surf_face_topo,
                                           verts_face_map,
                                           surf_edits_loczonefaces,
                                           surf_edits_corners);
      }
      // Append to blueprint mesh node
      conduit::Node &surface_edit = mMeshNode["teton/surface_edits/" + surface_edit_name_str];
      surface_edit["corners"].set(surf_edits_corners);
      surface_edit["local_zone_faces"].set(surf_edits_loczonefaces);
   }
}

/**
 @brief Compute data that lives on the face and associates a face with boundaries,
        if applicable. If the face is not associated with a boundary then it is a
        face that is shared with another domain, figure out which using the face
        adjset. The face adjset is an element-associated adjset that contains a list
        of face ids that touch other domains.

 @param boundaries A map of boundaries to faces where the keys in the map are values
                   present in the fields/boundary_attribute/values field and the
                   values are vectors of face ids in the main_faces topology that
                   use those boundary numbers.
 @param[out] nbelem_corner_faces The number of corner faces that touch boundary faces.
 @param[out] boundaries_types Counts of the various boundary types (how many faces of each)
 @param[out] boundary_conditions An encoding of boundary condition ids, number of corner faces.
 @param[out] face_to_bcid A "field" for the face mesh that indicates for each face
                          the boundary condition id for that face.
 */
void TetonBlueprint::ComputeFaceIDs(std::map<int, std::vector<int>> &boundaries,
                                    int &nbelem_corner_faces,
                                    std::vector<int> &boundaries_types,
                                    std::vector<int> &boundary_conditions,
                                    std::vector<int> &face_to_bcid,
                                    int rank)
{
   CALI_CXX_MARK_FUNCTION;

   std::map<int, int> boundary_id_to_type;
   if (mParametersNode.has_path("boundary_conditions/id_to_type_map"))
   {
      conduit::int_accessor boundary_ids = mParametersNode.fetch_existing("boundary_conditions/id_to_type_map/ids")
                                              .value();
      conduit::int_accessor teton_bc_types = mParametersNode.fetch_existing("boundary_conditions/id_to_type_map/types")
                                                .value();

      TETON_VERIFY_C(rank,
                     boundary_ids.number_of_elements() == teton_bc_types.number_of_elements(),
                     "boundary ids list length != boundary types list length");

      for (int i = 0; i < boundary_ids.number_of_elements(); ++i)
      {
         boundary_id_to_type[boundary_ids[i]] = teton_bc_types[i];
      }
   }

   std::map<int, int> boundary_id_to_ncornerfaces;
   std::map<int, int> number_boundary_types;
   for (int j = 31; j < 37; ++j)
      number_boundary_types[j] = 0;

   const conduit::Node &face_topology = mMeshNode["topologies/main_face"];
   const int nfaces = conduit::blueprint::mesh::utils::topology::length(face_topology);
   const conduit::Node &base_topology = mMeshNode["topologies/main"];
   const int ndim = conduit::blueprint::mesh::utils::topology::dims(base_topology);

   // First tag each face with the original boundary condition ID from the mesh.
   // This will later be changed to the re-enumerated boundary condition ID (so
   // each ID is consecutive and reordered by boundary condition type).
   face_to_bcid.resize(nfaces);
   for (int face = 0; face < nfaces; ++face)
      face_to_bcid[face] = 0;

   for (auto itr = boundaries.begin(); itr != boundaries.end(); ++itr)
   {
      int bc_id = itr->first;
      try
      {
         int bc_type = boundary_id_to_type.at(bc_id);
         number_boundary_types[bc_type] += 1;
      }
      catch (const std::exception &e)
      {
         // This boundary id does not have an associated boundary condition definition.
         // It is either invalid, or the calling code did not set up all the boundary
         // condition definitions needed.
         TETON_FATAL_C(rank,
                       "No boundary condition definition exists for value "
                          << bc_id << ", found in the boundary attributes field.  " << e.what())
      }
      boundary_id_to_ncornerfaces[bc_id] = 0;
      for (size_t j = 0; j < (itr->second).size(); ++j)
      {
         int face = (itr->second)[j];
         size_t ncorner_faces = face_to_corners2[face].size();
         boundary_id_to_ncornerfaces[bc_id] += ncorner_faces;
         face_to_bcid[face] = bc_id;
      }
   }

   int num_face_nbrs_for_teton = 0;

   if (mMeshNode.has_path("adjsets/main_face"))
   {
      num_face_nbrs_for_teton = mMeshNode.fetch_existing("adjsets/main_face/groups").number_of_children();
   }

   boundaries_types[0] = number_boundary_types.at(32);                                // number of reflecting
   boundaries_types[1] = number_boundary_types.at(35);                                // number of vacuum
   boundaries_types[2] = number_boundary_types.at(34) + number_boundary_types.at(36); // number of source
   boundaries_types[3] = num_face_nbrs_for_teton;                                     // number of shared

   // Order the boundaries by: reflecting, vacuum, source, shared.

   // First add the non-shared boundaries
   int num_nonshared_bndrs = boundaries_types[0] + boundaries_types[1] + boundaries_types[2];
   int num_bndrs = num_nonshared_bndrs + boundaries_types[3];
   int nreflec, nvacuum; //, nsource;
   nreflec = boundaries_types[0];
   nvacuum = boundaries_types[1];
   //   nsource = boundaries_types[2];
   //   We don't appear to use nsource for anything -- black27
   std::vector<int> bc_ids_ordered(num_bndrs);
   int jreflect = 0, jvaccuum = 0, jsource = 0;
   for (auto itr = boundaries.begin(); itr != boundaries.end(); ++itr)
   {
      int bc_id = itr->first;
      int bc_type = boundary_id_to_type.at(bc_id);
      if (bc_type == 32) // reflecting
      {
         bc_ids_ordered[jreflect] = bc_id;
         jreflect += 1;
      }
      else if (bc_type == 35) // vaccuum
      {
         int index = jvaccuum + nreflec;
         bc_ids_ordered[index] = bc_id;
         jvaccuum += 1;
      }
      else if (bc_type == 34 || bc_type == 36) // source
      {
         int index = jsource + nreflec + nvacuum;
         bc_ids_ordered[index] = bc_id;
         // TODO: CHANGE THE HARD-CODING HERE BY READING IN SOURCE ID TO PROFILE
         //       DESCRIPTION FROM A FILE
         jsource += 1;
      }
   }

   // Now compute the shared boundaries info
   // First we get the max of the non-shared bc_ids so that we
   // can assign a unique bc_id
   int bc_id_max = 0;
   for (auto itr = boundaries.begin(); itr != boundaries.end(); ++itr)
   {
      bc_id_max = std::max(itr->first, bc_id_max);
   }
   std::map<int, int> bc_id_to_rank;

   if (mMeshNode.has_path("adjsets/main_face"))
   {
      const conduit::Node &face_adjset = mMeshNode.fetch_existing("adjsets/main_face");
      conduit::NodeConstIterator groups_it = face_adjset["groups"].children();
      int fn_counter_teton = 0;
      while (groups_it.has_next())
      {
         const conduit::Node &group = groups_it.next();
         int num_sfaces = group["values"].dtype().number_of_elements();
         int nbr_rank = group["neighbors"].to_int();
         int bc_id = fn_counter_teton + bc_id_max + 1;
         int index = fn_counter_teton + num_nonshared_bndrs;
         // We want to increment this only when num_sfaces > 0
         fn_counter_teton += 1;

         bc_ids_ordered[index] = bc_id;
         boundary_id_to_type[bc_id] = 33; // shared boundary
         bc_id_to_rank[bc_id] = nbr_rank;
         boundary_id_to_ncornerfaces[bc_id] = 0;
         // tag each shared face in this shared boundary with its b.c. ID
         const conduit::Node &sfaces = group["values"];
         for (int f = 0; f < num_sfaces; f++)
         {
            conduit::Node face_data(conduit::DataType(sfaces.dtype().id(), 1), (void *) sfaces.element_ptr(f), true);
            const int face = face_data.to_int();

            // Make sure that the shared face is actually a shared face
            if (face_to_zones2[face].size() > 1)
            {
               std::cout << "WARNING (TetonBlueprint::ComputeFaceIDs): the shared face " << face << " on rank " << rank
                         << " is actually an interior face "
                         << "\n";
               continue;
            }

            size_t ncorners = face_to_corners2[face].size();
            boundary_id_to_ncornerfaces[bc_id] += ncorners;
            face_to_bcid[face] = bc_id;
         }
      }
   }

   // Create boundary conditions array for Teton
   nbelem_corner_faces = 0;
   boundary_conditions.push_back(num_bndrs);
   std::map<int, int> bc_id_to_new_bc_id;
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int bc_type = boundary_id_to_type.at(bc_id);
      boundary_conditions.push_back(bc_type);
      bc_id_to_new_bc_id[bc_id] = j;
   }
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int ncorner_faces = boundary_id_to_ncornerfaces.at(bc_id);
      boundary_conditions.push_back(ncorner_faces);
      nbelem_corner_faces += ncorner_faces;
   }
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int bc_type = boundary_id_to_type.at(bc_id);
      if (bc_type == 33) // shared
      {
         int nbr_rank = bc_id_to_rank[bc_id];
         boundary_conditions.push_back(nbr_rank);
      }
      else
      {
         boundary_conditions.push_back(-1);
      }
   }

   // Finally, fill in face_to_bcid array. Here the bcid (boundary condition ID)
   // corresponds to the re-enumerated boundary condition ID (consecutive and sorted
   // by boundary type)

   for (int face = 0; face < nfaces; ++face)
   {
      int bc_id = face_to_bcid[face];
      if (bc_id > 0) // if this is a boundary face
      {
         face_to_bcid[face] = bc_id_to_new_bc_id[bc_id] + 1;
      }
   }

   // CONDUIT OUTPUT
   // Add to conduit parameters input file
   mParametersNode["boundary_conditions/num_reflecting"] = boundaries_types[0];
   mParametersNode["boundary_conditions/num_vacuum"] = boundaries_types[1];
   mParametersNode["boundary_conditions/num_source"] = boundaries_types[2];
   mParametersNode["boundary_conditions/num_comm"] = boundaries_types[3];
   mParametersNode["boundary_conditions/num_total"] = num_bndrs;

   nbelem_corner_faces = 0;
   // TODO: change these names
   std::vector<int> bc_types, num_bc_cornerfaces, neighbor_ids;
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int bc_type = boundary_id_to_type.at(bc_id);
      bc_types.push_back(bc_type);
   }
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int ncorner_faces = boundary_id_to_ncornerfaces.at(bc_id);
      num_bc_cornerfaces.push_back(ncorner_faces);
      nbelem_corner_faces += ncorner_faces;
   }
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int bc_type = boundary_id_to_type.at(bc_id);
      if (bc_type == 33) // shared
      {
         int nbr_rank = bc_id_to_rank[bc_id];
         neighbor_ids.push_back(nbr_rank);
      }
      else
      {
         neighbor_ids.push_back(-1);
      }
   }

   mParametersNode["boundary_conditions/type"].set(bc_types.data(), bc_types.size());
   mParametersNode["boundary_conditions/corner_face_ids"].set(num_bc_cornerfaces.data(), num_bc_cornerfaces.size());
   mParametersNode["boundary_conditions/neighbor_ids"].set(neighbor_ids.data(), neighbor_ids.size());
}

// Need to populate
//     BCTypeInt = {type1, type2}, where 31 <= type1 <= 37 is Teton's boundary ID type mapping,
//     where type1 corresponds to the first boundary and type2 corresponds to the second boundary.
//     Here BCCornerFaces = {1,1} since always 1 corner per boundary condition (is this ever not true?).
//     BCNeighborID = {rank1, rank2}, where rank1 corresponds to the corresponding rank of a shared boundary.
//     If not a shared boundary, then rank1 = -1 (same with rank2).
//     teton_addboundary(&numBCTotal, &BCTypeInt[0], &BCCornerFaces[0], &BCNeighborID[0]);
// Also need
//         teton_setzone1d(&zoneID, &numBCTotal, &BCZoneID[0]);
// Here numBCTotal = 2 (is this ever not true)? and
//    Teton expects two boundary conditions, one for zone1D == 1 and one for zoneID == nzones
//    Here BCZoneID[0] == 1, BCZoneID[1] == nzones or BCZoneID[0] == nzones and BCZoneID[1] == 1

void TetonBlueprint::ComputeFaceIDs1D(int *boundary_connectivity,         // boundary to vertex mapping
                                      int *boundary_attributes,           // boundary to mesh attribute mapping
                                      std::vector<int> &boundaries_types, // boundary to Teton boundary type mapping
                                      int rank)
{
   // need to fill these out and store them in mMeshNode
   std::vector<int> bc_type_int(2);
   std::vector<int> bc_corner_faces(2);
   std::vector<int> bc_nghbr_id(2);
   std::vector<int> bc_zone_id(2);

   // TODO: check if there is ever a case where this doesn't hold
   bc_corner_faces[0] = 1;
   bc_corner_faces[1] = 1;

   std::map<int, int> boundary_id_to_type;
   // map boundary ID to corresponding Teton boundary type (c.g. 34 for a vaccuum boundary condition)
   if (mParametersNode.has_path("boundary_conditions/id_to_type_map"))
   {
      conduit::int_accessor boundary_ids = mParametersNode.fetch_existing("boundary_conditions/id_to_type_map/ids")
                                              .value();
      conduit::int_accessor teton_bc_types = mParametersNode.fetch_existing("boundary_conditions/id_to_type_map/types")
                                                .value();

      TETON_VERIFY_C(rank,
                     boundary_ids.number_of_elements() == teton_bc_types.number_of_elements(),
                     "boundary ids list length != boundary types list length");

      for (int i = 0; i < boundary_ids.number_of_elements(); ++i)
      {
         boundary_id_to_type[boundary_ids[i]] = teton_bc_types[i];
      }
   }

   std::map<int, int> number_boundary_types;
   for (int j = 31; j < 37; ++j)
      number_boundary_types[j] = 0;

   // First tag each face with the original boundary condition ID from the mesh.
   // This will later be changed to the re-enumerated boundary condition ID (so
   // each ID is consecutive and reordered by boundary condition type).
   //elem_to_bcid.resize(nelem);
   //for (int face = 0; face < nfaces; ++face)
   //   elem_to_bcid[face] = 0;

   for (int j = 0; j < 2; ++j)
   {
      int bc_id = boundary_attributes[j];
      {
         int bc_type = boundary_id_to_type.at(bc_id);
         number_boundary_types[bc_type] += 1;
      }
   }

   int num_face_nbrs_for_teton = 0;

   if (mMeshNode.has_path("adjsets/main_face"))
   {
      num_face_nbrs_for_teton = mMeshNode.fetch_existing("adjsets/main_face/groups").number_of_children();
   }

   boundaries_types[0] = number_boundary_types.at(32);                                // number of reflecting
   boundaries_types[1] = number_boundary_types.at(35);                                // number of vacuum
   boundaries_types[2] = number_boundary_types.at(34) + number_boundary_types.at(36); // number of source
   boundaries_types[3] = num_face_nbrs_for_teton;                                     // number of shared

   // Order the boundaries by: reflecting, vacuum, source, shared.

   // First add the non-shared boundaries
   int num_nonshared_bndrs = boundaries_types[0] + boundaries_types[1] + boundaries_types[2];
   int num_bndrs = num_nonshared_bndrs + boundaries_types[3];
   int nreflec, nvacuum; //, nsource;
   nreflec = boundaries_types[0];
   nvacuum = boundaries_types[1];
   //   nsource = boundaries_types[2];
   //   We don't appear to use nsource for anything -- black27
   std::vector<int> bc_ids_ordered(num_bndrs);
   int jreflect = 0, jvaccuum = 0, jsource = 0;
   for (int j = 0; j < 2; ++j)
   {
      int bc_id = boundary_attributes[j];
      int bc_type = boundary_id_to_type.at(bc_id);
      if (bc_type == BC_Type::bc_reflecting) // reflecting
      {
         bc_ids_ordered[jreflect] = bc_id;
         bc_type_int[jreflect] = bc_type;
         int bc_vertex_id = boundary_connectivity[j];
         if (bc_vertex_id == 0)
         {
            bc_zone_id[jreflect] = 1; // vertexID 0 corresponds to zoneID = 1 in Teton
         }
         else
         {
            bc_zone_id[jreflect] = bc_vertex_id; // vertexID nzones corresponds to zoneID = nzones in Teton
         }
      }
      else if (bc_type == BC_Type::bc_vaccuum) // vaccuum
      {
         int index = jvaccuum + nreflec;
         bc_ids_ordered[index] = bc_id;
         bc_type_int[index] = bc_type;
         int bc_vertex_id = boundary_connectivity[j];
         if (bc_vertex_id == 0)
         {
            bc_zone_id[index] = 1; // vertexID 0 corresponds to zoneID = 1 in Teton
         }
         else
         {
            bc_zone_id[index] = bc_vertex_id; // vertexID nzones corresponds to zoneID = nzones in Teton
         }
         jvaccuum += 1;
      }
      else if (bc_type == BC_Type::bc_source_temp || bc_type == BC_Type::bc_source_fd) // source
      {
         int index = jsource + nreflec + nvacuum;
         bc_ids_ordered[index] = bc_id;
         bc_type_int[index] = bc_type;
         int bc_vertex_id = boundary_connectivity[j];
         if (bc_vertex_id == 0)
         {
            bc_zone_id[index] = 1; // vertexID 0 corresponds to zoneID = 1 in Teton
         }
         else
         {
            bc_zone_id[index] = bc_vertex_id; // vertexID nzones corresponds to zoneID = nzones in Teton
         }
         jsource += 1;
      }
   }

   // Now compute the shared boundaries info

   if (mMeshNode.has_path("adjsets/mesh"))
   {
      const conduit::Node &face_adjset = mMeshNode.fetch_existing("adjsets/mesh");
      conduit::NodeConstIterator groups_it = face_adjset["groups"].children();
      int fn_counter_teton = 0;
      while (groups_it.has_next())
      {
         const conduit::Node &group = groups_it.next();
         int nbr_rank = group["neighbors"].to_int();
         int index = fn_counter_teton + num_nonshared_bndrs;
         // We want to increment this only when num_sfaces > 0
         fn_counter_teton += 1;

         //int *szones_ptr = group["values"].value();
         const conduit::int32 *szones_ptr = group["values"].as_int32_ptr();
         bc_type_int[index] = BC_Type::bc_shared; // shared boundary
         bc_zone_id[index] = szones_ptr[0];       // only one shared zone per face-neighbor in 1d
         bc_nghbr_id[index] = nbr_rank;
         bc_corner_faces[index] = 1;
      }
   }

   // Create boundary conditions array for Teton

   // CONDUIT OUTPUT
   // Add to conduit parameters input file
   mParametersNode["boundary_conditions/num_reflecting"] = boundaries_types[0];
   mParametersNode["boundary_conditions/num_vacuum"] = boundaries_types[1];
   mParametersNode["boundary_conditions/num_source"] = boundaries_types[2];
   mParametersNode["boundary_conditions/num_comm"] = boundaries_types[3];
   mParametersNode["boundary_conditions/num_total"] = num_bndrs;
   mParametersNode["boundary_conditions/type"].set(bc_type_int.data(), 2);
   mParametersNode["boundary_conditions/zone_ids"].set(bc_zone_id.data(), 2);
   mParametersNode["boundary_conditions/neighbor_ids"].set(bc_nghbr_id.data(), 2);
   mParametersNode["boundary_conditions/bc_ncorner_faces"].set(bc_corner_faces.data(), 2);
}

void TetonBlueprint::OutputTetonMesh(int rank, MPI_Comm comm)
{
   const conduit::Node &coords = mMeshNode["coordsets/coords/values"];
   if (!coords.has_path("y") && !coords.has_path("z"))
   {
      int nelem = mMeshNode["topologies/mesh/elements/dims/i"].value();
      // TODO Move to mesh conduit Node
      //  TODO: fix rank
      mParametersNode["size/ndim"] = 1;
      mParametersNode["size/rank"] = rank;
      mParametersNode["size/nzones"] = nelem;
      mParametersNode["size/nverts"] = 2 * nelem;
      mParametersNode["size/ncornr"] = 2 * nelem;
      mParametersNode["size/nsides"] = 2 * nelem;
      mParametersNode["size/nbelem"] = 2;
      mParametersNode["size/maxcf"] = 2;
      mParametersNode["size/maxCorner"] = 2;
      // TODO: enable slab geometry option as well
      mParametersNode["size/geomType"] = 42; //spherical
      //mParametersNode["size/geomType"] = 41; //slab
      int num_face_nbrs_teton = mParametersNode["boundary_conditions/num_comm"].value();
      mParametersNode["size/ncomm"] = num_face_nbrs_teton;
      mMeshNode["shared_boundaries/nsfaces"] = num_face_nbrs_teton;
      mParametersNode["shared_boundaries/nsfaces"] = num_face_nbrs_teton;

      return;
   }

   CALI_CXX_MARK_FUNCTION;
   verifyInput(mMeshNode, comm);

   //  face_to_corners2.
   CreateConnectivityArrays(mMeshNode, comm);

   const conduit::Node &base_topology = mMeshNode["topologies/main"];
   const conduit::Node &base_coordset = mMeshNode["coordsets/coords"];
   const conduit::Node &corner_topology = mMeshNode["topologies/main_corner"];
   const int ndim = conduit::blueprint::mesh::utils::topology::dims(base_topology);
   const int nvert = conduit::blueprint::mesh::utils::coordset::length(base_coordset);
   const int nelem = conduit::blueprint::mesh::utils::topology::length(base_topology);
   const int ncorners_total = conduit::blueprint::mesh::utils::topology::length(corner_topology);

   // Divide up the mesh boundary faces in to groups with
   // the same boundary condition. Note that this also includes
   // shared (processor boundary) faces.
   // NOTE(JRC): boundaries := { bcond_id => list of shared face ids }

   const conduit::Node &face_topology = mMeshNode["topologies/main_face"];
   const int nfaces = conduit::blueprint::mesh::utils::topology::length(face_topology);
   std::map<int, std::vector<int>> boundaries;

   if (mMeshNode.has_path("topologies/boundary"))
   {
      // Create the boundary attributes (to be used for boundary conditions).
      // After this call, mMeshNode["fields/face_attribute/values"]
      // corresponds to the face_to_attr[face] gives the boundary attribute
      // for each face
      CreateConduitFaceAttributes(mMeshNode, rank);

      // Populate boundaries, a list of face ids for each boundary condition id.
      const conduit::int32 *face_to_attr = mMeshNode["fields/face_attribute/values"].as_int32_ptr();

      try
      {
         // Find first the non-shared boundary faces
         for (int face = 0; face < nfaces; ++face)
         {
            bool is_bndry = (face_to_attr[face] > 0);
            // If a physical boundary face, get the boundary condition
            if (is_bndry)
            {
               int bdnry_cond_id = face_to_attr[face];
               if (boundaries.find(bdnry_cond_id) == boundaries.end())
               {
                  boundaries[bdnry_cond_id] = {face};
               }
               else
               {
                  boundaries.at(bdnry_cond_id).push_back(face);
               }
            }
         }
      }
      catch (const std::exception &e)
      {
         TETON_FATAL_C(rank, e.what())
      }
   }
   // Tag each face with a boundary condition ID (i.e., compute
   // m_face_to_bcid[face] = bcID for each mesh face).
   // The BC IDs are reenumerated from those encoded in the original mesh
   // by 1) ordering with respect to reflecting, vaccuum, and source
   //    2) also adding BC IDs corresponding to shared faces
   // For non-boundary faces, m_face_to_bcid[face] = 0.

   std::vector<int> boundaries_types = {0, 0, 0, 0};
   std::vector<int> boundary_conditions;
   int nbelem_corner_faces = 0;
   m_face_to_bcid.resize(nfaces, 0);

   ComputeFaceIDs(boundaries, nbelem_corner_faces, boundaries_types, boundary_conditions, m_face_to_bcid, rank);

   // Compute the shared faces arrays

   // shared_faces: shared face list (ordered by face neighbor) for this domain
   std::vector<int> shared_faces;

   if (mMeshNode.has_path("adjsets/main_face"))
   {
      const conduit::Node &face_adjset = mMeshNode["adjsets/main_face"];
      conduit::NodeConstIterator groups_it = face_adjset["groups"].children();
      while (groups_it.has_next())
      {
         const conduit::Node &fn_face_group = groups_it.next();
         const conduit::Node fn_corner_group = mMeshNode["adjsets/main_corner/groups"][fn_face_group.name()];

         // Store the index of each corner in the corner adjacency list.
         // Precondition: The order of the corners in the corner mesh adjacency set from
         // conduit must match with the order in the neighboring domain's list.
         std::map<int, int> cpoint_to_gindex;
         const auto corner_data = fn_corner_group["values"].as_int_accessor();
         for (int cp = 0; cp < corner_data.number_of_elements(); cp++)
         {
            int cpoint = corner_data[cp];
            cpoint_to_gindex[cpoint] = cp;
         }

         const auto face_data = fn_face_group["values"].as_int_accessor();
         for (int f = 0; f < face_data.number_of_elements(); f++)
         {
            const int face = face_data[f];

            // Makes sure a shared face is really a shared face
            if (face_to_zones2[face].size() > 1)
               continue;

            const int bcid = m_face_to_bcid[face];
            const int zone = face_to_zones2[face].front(); // only 1 b/c on boundary

            std::map<int, int> fpoint_to_corner;
            for (const auto corner : face_to_corners2[face])
            {
               const std::vector<conduit::index_t>
                  corner_points = conduit::blueprint::mesh::utils::topology::unstructured::points(corner_topology,
                                                                                                  corner);
               // NOTE: A corner's vertices are ordered such that the corner vertex
               // is first, which will be the same vertex as is used on the face.
               fpoint_to_corner[corner_points.front()] = corner;
            }

            const std::vector<conduit::index_t>
               face_points = conduit::blueprint::mesh::utils::topology::unstructured::points(face_topology, face);
            std::vector<std::pair<int, int>> fpoint_list;
            for (const auto face_point : face_points)
            {
               fpoint_list.emplace_back(cpoint_to_gindex[face_point], face_point);
            }
            std::sort(fpoint_list.begin(), fpoint_list.end());

            if (ndim == 3)
            {
               for (const auto &fpoint_pair : fpoint_list)
               {
                  const int &fpoint = fpoint_pair.second;
                  const int &fpoint_corner = fpoint_to_corner[fpoint];

                  // 3D: for each shared **corner face**, {bcid, zone, face, corner}
                  shared_faces.push_back(bcid);
                  shared_faces.push_back(zone);
                  shared_faces.push_back(face);
                  shared_faces.push_back(fpoint_corner);
               }
            }
            else if (ndim == 2)
            {
               // 2D: for each shared **face**, {bcid, zone, face, corner1, corner2}
               shared_faces.push_back(bcid);
               shared_faces.push_back(zone);
               shared_faces.push_back(face);
               for (const auto &fpoint_pair : fpoint_list)
               {
                  const int &fpoint = fpoint_pair.second;
                  const int &fpoint_corner = fpoint_to_corner[fpoint];
                  shared_faces.push_back(fpoint_corner);
               }
            }
         }
      }
   }

   // Put mesh info in the conduit mesh
   mParametersNode["boundaries/boundary_conditions"].set(boundary_conditions.data(), boundary_conditions.size());
   mParametersNode["boundaries/boundaries_types"].set(boundaries_types.data(), boundaries_types.size());

   // CONDUIT OUTPUT
   // TODO Move to mesh conduit Node
   mParametersNode["size/ndim"] = ndim;
   mParametersNode["size/rank"] = rank;
   mParametersNode["size/nzones"] = nelem;
   mParametersNode["size/nverts"] = nvert;
   mParametersNode["size/ncornr"] = ncorners_total;
   mParametersNode["size/nsides"] = ncorners_total;
   mParametersNode["size/nbelem"] = nbelem_corner_faces;
   mParametersNode["size/maxCorner"] = ComputeMaxCorners(base_topology);
   if (ndim == 2)
   {
      mParametersNode["size/maxcf"] = 2;
      mParametersNode["size/geomType"] = 44; // rz
   }
   else
   {
      mParametersNode["size/maxcf"] = 3;
      mParametersNode["size/geomType"] = 45; // xyz
   }

   int num_face_nbrs_teton = 0;
   if (mMeshNode.has_path("adjsets/main_face"))
   {
      num_face_nbrs_teton = mMeshNode.fetch_existing("adjsets/main_face/groups").number_of_children();
   }
   mParametersNode["size/ncomm"] = num_face_nbrs_teton;

   // TODO: maybe get a better formula for this
   size_t n_shared_corner_faces = shared_faces.size() / ((ndim == 2) ? 5 : 4);
   mMeshNode["shared_boundaries/nsfaces"] = n_shared_corner_faces;
   // TODO: CHANGE
   mParametersNode["shared_boundaries/nsfaces"] = n_shared_corner_faces;
   std::vector<int> shared_faces_reord;
   if (n_shared_corner_faces > 0)
   {
      mMeshNode["shared_boundaries/shared_faces"].set(shared_faces.data(), shared_faces.size());
   }

   // This transforms the connectivity data in "blueprint" format to
   // the format that standalone Teton expects
   ComputeTetonMeshConnectivityArray(rank);
   CreateTetonMeshCornerCoords();

   //  Transform surface flux tally information to a format that Teton can use
   if (mParametersNode.has_path("surface_edits"))
   {
      ProcessSurfaceEdits(rank);
   }

   verifyInput(mMeshNode, comm);
}
