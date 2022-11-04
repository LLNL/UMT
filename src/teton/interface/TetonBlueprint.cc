
#include "conduit/conduit.hpp"
#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_blueprint_mpi.hpp"
#include "conduit/conduit_blueprint_mesh_utils.hpp"
#include "conduit/conduit_blueprint_o2mrelation.hpp"
#include "conduit/conduit_blueprint_o2mrelation_iterator.hpp"
#include "conduit/conduit_relay.hpp"
#include "conduit/conduit_relay_mpi_io_blueprint.hpp"
#include "TetonBlueprint.hh"
#include "dbc_macros.h"

#if defined(TETON_ENABLE_CALIPER)
#include "caliper/cali.h"
#else
#define CALI_MARK_BEGIN(label)
#define CALI_MARK_END(label)
#endif

#include <map>
#include <vector>
#include <set>
#include <algorithm> // for sort

using std::cout;
using std::endl;
//using namespace mfem;
// access one-to-many index types
namespace O2MIndex = conduit::blueprint::o2mrelation;

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
   nhalffaces = offset;
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

   double node_x = corner_to_node_x[corner];
   double node_y = corner_to_node_y[corner];
   double node_z;
   int dim = mParametersNode.fetch_existing("size/ndim").as_int32();
   if (dim > 2)
      node_z = corner_to_node_z[corner];

   for (size_t c = 0; c < ncorners; ++c)
   {
      corner1 = face_to_corners2[face][c];
      double node1_x = corner_to_node_x[corner1];
      double node1_y = corner_to_node_y[corner1];
      double node1_z;
      if (dim > 2)
         node1_z = corner_to_node_z[corner1];
      int zone1 = corner_to_zone[corner1];

      // Compilers frequently throw warnings about floating point comparisons being unreliable.
      // In this case, we should have identical binary comparisons because we are just comparing positions.
      if (dim > 2)
      {
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
      std::pair<int, int> zoneface;
      zoneface.first = zone;
      zoneface.second = face;
      try
      {
         int halfface = zoneface_to_halfface.at(zoneface);
         return -(nbelem + halfface + 1);
      }
      catch (const std::exception &e)
      {
         TETON_FATAL_C(rank, e.what() )
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
   conduit::int32 *shared_faces_array_ptr;
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
      std::pair<int, int> zoneface;
      zoneface.first = zone;
      zoneface.second = face;
      int halfface;
      try
      {
         halfface = zoneface_to_halfface.at(zoneface);
      }
      catch (const std::exception &e)
      {
         TETON_FATAL_C(rank, e.what() )
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

   // TODO: QUESTION: DO WE NEED TO DO THIS ?
   if (nsfaces > 0)
   {
      mMeshNode["shared_boundaries/shared_faces"].set(shared_faces_array_ptr, sfaces_array_len);
   }
}

void TetonBlueprint::ComputeTetonMeshConnectivityArray(int rank)
{
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
   int nzones = mParametersNode.fetch_existing("size/nzones").as_int32();
   int ndim = mParametersNode.fetch_existing("size/ndim").as_int32();
   std::vector<double> zone_verts;
   //   std::vector<int> zone_ncorners;

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
   int rank;
   MPI_Comm_rank(comm, &rank);

   conduit::int32 *zone_to_faces;
   conduit::int32 *zone_to_corners;
   conduit::int32 *face_to_zones;

   // zone-to-vertices
   conduit::Node &base_topology = meshNode.fetch_existing("topologies/main");
   conduit::Node &base_coordset = meshNode.fetch_existing("coordsets/coords");
   conduit::blueprint::mesh::topology::unstructured::generate_offsets(base_topology, base_topology["elements/offsets"]);
   const int nzones = conduit::blueprint::mesh::utils::topology::length(base_topology);
   const int dim = conduit::blueprint::mesh::utils::topology::dims(base_topology);

   TETON_VERIFY_C(rank, (dim == 2 || dim == 3), "Only 2d or 3d blueprint meshes supported.");

   std::vector<std::vector<int>> zone_to_vertices(nzones);
   conduit::int32 *zone_to_verts_offsets = base_topology.fetch_existing("elements/offsets").as_int32_ptr();
   conduit::int32 *zone_to_verts = base_topology.fetch_existing("elements/connectivity").as_int32_ptr();
   int nzone2verts = base_topology.fetch_existing("elements/connectivity").dtype().number_of_elements();
   for (int zone = 0; zone < nzones; ++zone)
   {
      int offset = zone_to_verts_offsets[zone];
      int nverts_for_zone;
      if (zone == (nzones - 1)) // last zone
      {
         nverts_for_zone = nzone2verts - offset;
      }
      else
      {
         int offset2 = zone_to_verts_offsets[zone + 1];
         nverts_for_zone = offset2 - offset;
      }
      for (int v = 0; v < nverts_for_zone; ++v)
      {
         int vid = zone_to_verts[offset + v];
         zone_to_vertices[zone].push_back(vid);
      }
   }

   // Create the zone-to-corner and corner-to-zone connectivity info
   conduit::Node &corner_topology = meshNode["topologies/main_corner"];
   conduit::Node &corner_coords = meshNode["coordsets/coords_centroids"];
   conduit::Node &zone2corner_map = meshNode["relations/zone2corner"];
   conduit::Node &corner2zone_map = meshNode["relations/corner2zone"];

   // Call the correct version of generate_corners for a serial vs parallel mesh.
   if (meshNode.has_path("adjsets/main_adjset"))
   {
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
      conduit::blueprint::mesh::topology::unstructured::generate_corners(
          base_topology, corner_topology, corner_coords, zone2corner_map, corner2zone_map);
   }

   conduit::blueprint::mesh::topology::unstructured::generate_offsets(corner_topology,
                                                                      corner_topology["elements/offsets"]);
   // For creating the Teton connectivity array
   zone_to_corners = zone2corner_map["values"].as_int32_ptr();
   conduit::int32 *zone_to_ncorner = zone2corner_map["sizes"].as_int32_ptr();

   // corner-to-verts
   const int ncorners = conduit::blueprint::mesh::utils::topology::length(corner_topology);
   std::vector<int> corner_to_vertex(ncorners);
   zone_to_corners2.resize(nzones);
   corner_to_node_x.resize(ncorners);
   corner_to_node_y.resize(ncorners);
   corner_to_node_z.resize(ncorners);
   corner_to_zone.resize(ncorners);

   const conduit::Node &zone2corner_data = zone2corner_map["values"];
   conduit::blueprint::o2mrelation::O2MIterator z2c_iter(zone2corner_map);
   while (z2c_iter.has_next(O2MIndex::DATA))
   {
      z2c_iter.next(O2MIndex::ONE);
      const conduit::index_t zid = z2c_iter.index(O2MIndex::ONE);

      z2c_iter.to_front(O2MIndex::MANY);
      while (z2c_iter.has_next(O2MIndex::MANY))
      {
         z2c_iter.next(O2MIndex::MANY);
         const conduit::index_t cidx = z2c_iter.index(O2MIndex::DATA);

         conduit::Node cdata(
             conduit::DataType(zone2corner_data.dtype().id(), 1), (void *) zone2corner_data.element_ptr(cidx), true);
         const conduit::index_t cid = cdata.to_index_t();

         zone_to_corners2[zid].push_back(cid);
         corner_to_zone[cid] = zid;

         CALI_MARK_BEGIN("Teton_Conduit_Points");
         const std::vector<conduit::index_t> corner_points
             = conduit::blueprint::mesh::utils::topology::unstructured::points(corner_topology, cid);
         CALI_MARK_END("Teton_Conduit_Points");

         const conduit::index_t cvid = corner_points.front();
         const std::vector<conduit::float64> this_corner_coords
             = conduit::blueprint::mesh::utils::coordset::_explicit::coords(base_coordset, cvid);

         corner_to_vertex[cid] = cvid;
         corner_to_node_x[cid] = this_corner_coords[0];
         corner_to_node_y[cid] = this_corner_coords[1];
         if (dim == 3)
         {
            corner_to_node_z[cid] = this_corner_coords[2];
         }
      }
   }

   // Create the zone-to-face and face-to-zone connectivity info
   conduit::Node &face_topology = meshNode["topologies/main_face"];
   conduit::Node &zone2face_map = meshNode["relations/zone2face"];
   conduit::Node &face2zone_map = meshNode["relations/face2zone"];

   // Call the correct version of generate_faces/lines corners for a serial vs parallel mesh.
   if (meshNode.has_path("adjsets/main_adjset"))
   {
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
      if (dim == 3)
      {
         conduit::blueprint::mesh::topology::unstructured::generate_faces(
             base_topology, face_topology, zone2face_map, face2zone_map);
      }
      else if (dim == 2)
      {
         conduit::blueprint::mesh::topology::unstructured::generate_lines(
             base_topology, face_topology, zone2face_map, face2zone_map);
      }
   }

   conduit::blueprint::mesh::topology::unstructured::generate_offsets(face_topology, face_topology["elements/offsets"]);
   if (mMeshNode.has_path("topologies/boundary"))
   {
      conduit::Node &bndry_topology = meshNode["topologies/boundary"];
      conduit::blueprint::mesh::topology::unstructured::generate_offsets(bndry_topology,
                                                                         bndry_topology["elements/offsets"]);
   }
   // zone-to-faces
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
   // face-to-zones
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
   // face-to-vertices
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

   // Create face_to_corners
   // Note: we use that the vertex ID for each face is shared
   // by exactly one corner in each zone that has this face
   // TODO: need to reverse order of corners for one zone of each pair of zones
   // sharing a face
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

   // save for mesh motion and forces
   meshNode["arrays/corner_to_vertex"].set(corner_to_vertex.data(), corner_to_vertex.size());
   meshNode["arrays/zone_to_ncorners"].set(zone_to_ncorner, nzones);
   meshNode["arrays/zone_to_corners"].set(zone_to_corners, ncorners);
}

void TetonBlueprint::CreateConduitFaceAttributes(conduit::Node &meshNode, int rank)
{
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

   std::map<std::set<int>, int> verts_bndryid_map; // {face vertices => face mfem id}
   const conduit::Node *face_topos[] = {&bndry_topo, &face_topology};
   for (int i = 0; i < 2; ++i)
   {
      const conduit::Node &face_topo = *face_topos[i];
      const int face_topo_length = conduit::blueprint::mesh::utils::topology::length(face_topo);
      for (int f = 0; f < face_topo_length; ++f)
      {
         std::vector<conduit::index_t> face_points
             = conduit::blueprint::mesh::utils::topology::unstructured::points(face_topo, f);
         std::set<int> face_vertices(face_points.begin(), face_points.end());

         if (i == 0)
         {
            verts_bndryid_map[face_vertices] = f;
         }
         else
         {
            src_data.reset();
            // if the face isn't a boundary face, set its attribute value to 0
            if (verts_bndryid_map.find(face_vertices) == verts_bndryid_map.end())
            {
               src_data.set(0);
            }
            // if face is a boundary face, then use the attribute for the associated boundary
            else
            {
               try
               {
                  int b = verts_bndryid_map.at(face_vertices);
                  src_data.set_external(conduit::DataType(bndry_values.dtype().id(), 1),
                                        (void *) bndry_values.element_ptr(b));
               }
               catch (const std::exception &e)
               {
                  TETON_FATAL_C(rank, e.what() )
               }
            }
            dst_data.set_external(conduit::DataType(face_values.dtype().id(), 1), (void *) face_values.element_ptr(f));
            src_data.to_data_type(dst_data.dtype().id(), dst_data);
         }
      }
   }
}

void TetonBlueprint::ComputeFaceIDs(std::map<int, std::vector<int>> &boundaries,
                                    int &nbelem_corner_faces,
                                    std::vector<int> &boundaries_types,
                                    std::vector<int> &boundary_conditions,
                                    std::vector<int> &face_to_bcid,
                                    int rank)
{


   std::map<int, int> boundary_id_to_type;
   if (mParametersNode.has_path("boundary_conditions/id_to_type_map"))
   {

      conduit::int_accessor boundary_ids = mParametersNode.fetch_existing("boundary_conditions/id_to_type_map/ids").value();
      conduit::int_accessor teton_bc_types = mParametersNode.fetch_existing("boundary_conditions/id_to_type_map/types").value();

      TETON_VERIFY_C(rank, boundary_ids.number_of_elements() == teton_bc_types.number_of_elements(), "boundary ids list length != boundary types list length");

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
                           << bc_id << ", found in the boundary attributes field.  "  << e.what())
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
   int n_shared_corner_faces = 0;
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
               std::cout << "WARNING (TetonBlueprint::ComputeFaceIDs): the shared face " << face << " on rank " << rank <<  " is actually an interior face " << "\n";
               continue;
            }

            size_t ncorners = face_to_corners2[face].size();
            boundary_id_to_ncornerfaces[bc_id] += ncorners;
            face_to_bcid[face] = bc_id;

            if (ndim == 2)
            {
               n_shared_corner_faces += 1;
            }
            else if (ndim == 3)
            {
               n_shared_corner_faces += ncorners;
            }
            else
            {
               std::cout << "1D is not yet implemented! " << std::endl;
            }
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

void TetonBlueprint::OutputTetonMesh(int rank, MPI_Comm comm)
{
   verifyInput(mMeshNode, comm);

   //  face_to_corners2.
   CreateConnectivityArrays(mMeshNode, comm);

   const conduit::Node &base_coordset = mMeshNode["coordsets/coords"];
   const conduit::Node &base_topology = mMeshNode["topologies/main"];
   const conduit::Node &corner_topology = mMeshNode["topologies/main_corner"];
   const int nvert = conduit::blueprint::mesh::utils::coordset::length(base_coordset);
   const int nelem = conduit::blueprint::mesh::utils::topology::length(base_topology);
   const int ndim = conduit::blueprint::mesh::utils::topology::dims(base_topology);
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
                  //int newface = face_to_newface_mapping[face];
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
         TETON_FATAL_C(rank, e.what() )
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
   m_face_to_bcid.resize(nfaces);

   ComputeFaceIDs(
       boundaries, nbelem_corner_faces, boundaries_types, boundary_conditions, m_face_to_bcid, rank);

   // Compute the shared faces arrays

   // shared_faces: shared face list (ordered by face neighbor) for this domain
   std::vector<int> shared_faces;

   if (mMeshNode.has_path("adjsets/main_face"))
   {
      const conduit::Node &face_adjset = mMeshNode["adjsets/main_face"];
      //const conduit::Node &corner_adjset = mMeshNode["adjsets/main_corner"];
      conduit::NodeConstIterator groups_it = face_adjset["groups"].children();
      while (groups_it.has_next())
      {
         const conduit::Node &fn_face_group = groups_it.next();
         const conduit::Node fn_corner_group = mMeshNode["adjsets/main_corner/groups"][fn_face_group.name()];
         const conduit::Node &fn_faces = fn_face_group["values"];
         const conduit::Node &fn_corners = fn_corner_group["values"];

         // Store the index of each corner in the corner adjacency list.
         // Precondition: The order of the corners in the corner mesh adjacency set from
         // conduit must match with the order in the neighboring domain's list.
         std::map<int, int> cpoint_to_gindex;
         for (int cp = 0; cp < fn_corners.dtype().number_of_elements(); cp++)
         {
            conduit::Node corner_data(
                conduit::DataType(fn_corners.dtype().id(), 1), (void *) fn_corners.element_ptr(cp), true);
            cpoint_to_gindex[corner_data.to_int()] = cp;
         }

         for (int f = 0; f < fn_faces.dtype().number_of_elements(); f++)
         {
            conduit::Node face_data(
                conduit::DataType(fn_faces.dtype().id(), 1), (void *) fn_faces.element_ptr(f), true);
            const int face = face_data.to_int();

            // Makes sure a shared face is really a shared face
            if (face_to_zones2[face].size() > 1) continue;

            const int bcid = m_face_to_bcid[face];
            const int zone = face_to_zones2[face].front(); // only 1 b/c on boundary
            const std::vector<int> corners = face_to_corners2[face];

            const std::vector<conduit::index_t> face_points
                = conduit::blueprint::mesh::utils::topology::unstructured::points(face_topology, face);

            std::map<int, int> fpoint_to_corner;
            for (size_t c = 0; c < corners.size(); c++)
            {
               const int corner = corners[c];
               const std::vector<conduit::index_t> corner_points
                   = conduit::blueprint::mesh::utils::topology::unstructured::points(corner_topology, corner);
               // NOTE: A corner's vertices are ordered such that the corner vertex
               // is first, which will be the same vertex as is used on the face.
               fpoint_to_corner[corner_points.front()] = corner;
            }

            std::vector<std::pair<int, int>> fpoint_list;
            for (size_t fp = 0; fp < face_points.size(); fp++)
            {
               const int face_point = face_points[fp];
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
   // TODO: CHANGE THIS!!!!
   nbelem = nbelem_corner_faces;
   // TODO: CHANGE THIS!!!!
   if (ndim == 2)
   {
      mParametersNode["size/maxcf"] = 2;
      mParametersNode["size/maxCorner"] = 4;
      mParametersNode["size/geomType"] = 44; // rz
   }
   else
   {
      mParametersNode["size/maxcf"] = 3;
      mParametersNode["size/maxCorner"] = 8;
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

   verifyInput(mMeshNode, comm);
}
