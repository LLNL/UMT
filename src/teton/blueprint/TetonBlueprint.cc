#if defined(TETON_ENABLE_MFEM)

#include "mfem.hpp"
#include "conduit/conduit.hpp"
#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_relay.hpp"
#include "TetonBlueprint.hh"

#include <vector>


bool CompareCoords3D(std::vector<double> &a, std::vector<double> &b)
{
   double tol = 1e-12;
   if (a[0] < (b[0] - tol))
   {
      return true;
   }
   else
   {
      if ((std::abs(a[0] - b[0]) < tol) && (a[1] < (b[1] - tol)))
      {
         return true;
      }
      else
      {
         if ((std::abs(a[0] - b[0]) < tol) && (std::abs(a[1] - b[1]) < tol) && (a[2] < (b[2] - tol)))
         {
            return true;
         }
      }
   }
   return false;
}

bool CompareCoords2D(std::vector<double> &a, std::vector<double> &b)
{
   if (a[0] < b[0])
   {
      return true;
   }
   else
   {
      if ((a[0] == b[0]) && (a[1] <= b[1]))
      {
         return true;
      }
   }
   return false;
}

bool CompareCoords(std::vector<double> &a, std::vector<double> &b)
{
   if (a.size() == 2)
   {
      return CompareCoords2D(a, b);
   }
   else
   {
      return CompareCoords3D(a, b);
   }
}

void TetonBlueprint::ComputeFaceToBndryElement(mfem::Mesh *mesh, mfem::Array<int> &face_to_bndry)
{
   size_t nbndry = mesh->GetNBE();
   face_to_bndry.SetSize(mesh->GetNumFaces());
   int dim = mesh->SpaceDimension();

   face_to_bndry = -1;

   for (size_t i = 0; i < nbndry; ++i)
   {
      int face, orient;
      if (dim == 3)
      {
         mesh->GetBdrElementFace(i, &face, &orient);
      }
      else
      {
         face = mesh->GetBdrElementEdgeIndex(i);
      }
      face_to_bndry[face] = i;
   }
}

void TetonBlueprint::ComputeLfaceToSface(mfem::ParMesh *mesh, mfem::Array<int> &lface_to_sface)
{
   int n_sfaces = mesh->GetNSharedFaces();
   lface_to_sface.SetSize(mesh->GetNumFaces());
   lface_to_sface = -1;
   for (int sface = 0; sface < n_sfaces; ++sface)
   {
      int lface = mesh->GetSharedFace(sface);
      lface_to_sface[lface] = sface;
   }
}

//// TODO: MODIFY FOR NON-CONFORMING MESH

void TetonBlueprint::OutputConduitMesh(mfem::ParMesh *pmesh,
                                       mfem::Vector &density,
                                       mfem::Vector &heat_capacity,
                                       mfem::Vector &rad_temp,
                                       mfem::Vector &material_temp,
                                       const mfem::Vector &gr_bounds,
                                       mfem::Vector &abs_opacity,
                                       mfem::Vector &scat_opacity,
                                       mfem::Vector &electron_density)
{
   int nelems = pmesh->GetNE();
   std::vector<int> zone_to_faces;
   std::vector<int> zone_to_corners;
   std::vector<int> corner_to_vertex;
   std::vector<int> corner_to_zone;
   std::vector<double> corner_to_node_x;
   std::vector<double> corner_to_node_y;
   std::vector<double> corner_to_node_z;
   std::vector<int> face_to_corners;
   std::vector<int> face_to_zones;
   std::vector<int> face_to_bcid;

   int nfaces = pmesh->GetNumFaces();
   std::vector<std::vector<int>> face_to_corners2(nfaces);
   std::vector<std::vector<int>> face_to_zones2(nfaces);
   std::vector<std::vector<int>> face_to_corners2_reord(nfaces);
   std::vector<std::vector<int>> face_to_zones2_reord(nfaces);
   std::vector<int> zone_to_faces_reord;
   std::vector<int> face_to_corners_reord;
   std::vector<int> face_to_zones_reord;
   std::vector<int> face_to_bcid_reord(nfaces);

   // zone ID and (local-to-a-zone) vertex ID to corner ID
   std::vector<std::vector<int>> zone_and_lvertex_to_corner;
   std::map<std::pair<int, int>, int> zone_and_gvertex_to_lvertex;

   int dim = pmesh->Dimension();
   int corner_id = 0;
   zone_and_lvertex_to_corner.resize(nelems);

   // re-enumerate mesh faces

   mfem::Array<int> lface_to_sface(nfaces);
   mfem::Array<int> face_to_bndry(nfaces);
   ComputeLfaceToSface(pmesh, lface_to_sface);
   face_to_bndry = -1;
   ComputeFaceToBndryElement(pmesh, face_to_bndry);
   std::vector<int> face_to_newface_mapping(nfaces);
   mfem::Array<int> face_added(nfaces);
   face_added = -1;
   int face_counter = 0;
   vert_prev_global = 0;
   global_vert_counter = 0;
   int nverts_all = pmesh->GetNV();
   vert_to_counters.resize(nverts_all);
   for (int elem = 0; elem < nelems; ++elem)
   {
      std::vector<int> faces_ordered;
      std::vector<std::vector<int>> face_verts_ordered;
      std::vector<int> vertices_ordered;
      if (dim == 2)
      {
         OrderVertices(pmesh, elem, vertices_ordered, faces_ordered);
      }
      else
      {
         OrderVertices3D(pmesh, elem, vertices_ordered, faces_ordered, face_verts_ordered);
      }
      int nfaces_loc = faces_ordered.size();
      for (int j = 0; j < nfaces_loc; ++j)
      {
         int face = faces_ordered[j];
         if (face_added[face] < 0)
         {
            face_to_newface_mapping[face] = face_counter;
            face_counter++;
            face_added[face] = 1;
         }
      }
   }

   // Compute zone_to, face_to, and corner_to arrays

   vert_prev_global = 0;
   global_vert_counter = 0;
   for (int v = 0; v < nverts_all; ++v)
   {
      vert_to_counters[v].clear();
   }
   for (int elem = 0; elem < nelems; ++elem)
   {
      // Compute corner ordering and connectivity info
      std::vector<int> faces_ordered;
      std::vector<std::vector<int>> face_verts_ordered;
      std::vector<int> vertices_ordered;
      if (dim == 2)
      {
         OrderVertices(pmesh, elem, vertices_ordered, faces_ordered);
      }
      else
      {
         OrderVertices3D(pmesh, elem, vertices_ordered, faces_ordered, face_verts_ordered);
      }
      int nfaces_loc = faces_ordered.size();
      zone_to_faces.push_back(nfaces_loc);
      zone_to_faces_reord.push_back(nfaces_loc);

      for (int j = 0; j < nfaces_loc; ++j)
      {
         int face = faces_ordered[j];
         int newface = face_to_newface_mapping[face];
         zone_to_faces.push_back(face);
         face_to_zones2[face].push_back(elem);
         face_to_zones2_reord[newface].push_back(elem);
         zone_to_faces_reord.push_back(newface);
      }
      int nverts = vertices_ordered.size();
      zone_to_corners.push_back(nverts);
      zone_and_lvertex_to_corner[elem].resize(nverts);
      for (int j = 0; j < nverts; ++j)
      {
         int v = vertices_ordered[j];
         corner_to_vertex.push_back(v);
         corner_to_zone.push_back(elem);
         double vert_coords[3];
         pmesh->GetNode(v, vert_coords);
         zone_and_gvertex_to_lvertex[{elem, v}] = j;
         corner_to_node_x.push_back(vert_coords[0]);
         corner_to_node_y.push_back(vert_coords[1]);
         if (dim == 3)
         {
            corner_to_node_z.push_back(vert_coords[2]);
         }
         zone_to_corners.push_back(corner_id);
         zone_and_lvertex_to_corner[elem][j] = corner_id;
         corner_id += 1;
      }
      if (dim == 2)
      {
         // Note: faces correspond to {v1,v2}, {v2,v3}, {v3,v4}, and {v4,v1}
         for (int j = 0; j < nfaces_loc - 1; ++j)
         {
            int face = faces_ordered[j];
            int newface = face_to_newface_mapping[face];
            int v1 = vertices_ordered[j];
            int v2 = vertices_ordered[j + 1];
            int i1 = zone_and_gvertex_to_lvertex[{elem, v1}];
            int c1 = zone_and_lvertex_to_corner[elem][i1];
            int i2 = zone_and_gvertex_to_lvertex[{elem, v2}];
            int c2 = zone_and_lvertex_to_corner[elem][i2];
            face_to_corners2[face].push_back(c1);
            face_to_corners2[face].push_back(c2);
            face_to_corners2_reord[newface].push_back(c1);
            face_to_corners2_reord[newface].push_back(c2);
         }
         int face = faces_ordered[nfaces_loc - 1];
         int newface = face_to_newface_mapping[face];
         int v1 = vertices_ordered[nfaces_loc - 1];
         int v2 = vertices_ordered[0];
         int i1 = zone_and_gvertex_to_lvertex[{elem, v1}];
         int c1 = zone_and_lvertex_to_corner[elem][i1];
         int i2 = zone_and_gvertex_to_lvertex[{elem, v2}];
         int c2 = zone_and_lvertex_to_corner[elem][i2];
         face_to_corners2[face].push_back(c1);
         face_to_corners2[face].push_back(c2);
         face_to_corners2_reord[newface].push_back(c1);
         face_to_corners2_reord[newface].push_back(c2);
         // TODO: In addition to element vertices,
         //       need to also add "hanging nodes" for AMR meshes
      }
      else
      {
         for (int j = 0; j < nfaces_loc; ++j)
         {
            int face = faces_ordered[j];
            int newface = face_to_newface_mapping[face];
            for (int k = 0; k < face_verts_ordered[j].size(); ++k)
            {
               int v1 = face_verts_ordered[j][k];
               int i1 = zone_and_gvertex_to_lvertex[{elem, v1}];
               int c1 = zone_and_lvertex_to_corner[elem][i1];
               face_to_corners2[face].push_back(c1);
               face_to_corners2_reord[newface].push_back(c1);
            }
         }
      }
   }
   int ncorners_total = corner_id;

   // Output zone and face connectivity in to blueprint format
   std::vector<int> face_to_ncorners(nfaces);
   for (int face = 0; face < nfaces; ++face)
   {
      int ncorners = face_to_corners2[face].size();
      face_to_ncorners[face] = ncorners;
      face_to_corners.push_back(ncorners);
      for (int j = 0; j < ncorners; ++j)
      {
         int c = face_to_corners2[face][j];
         face_to_corners.push_back(c);
      }
      int nelem_on_face = face_to_zones2[face].size();
      face_to_zones.push_back(nelem_on_face);
      for (int j = 0; j < nelem_on_face; ++j)
      {
         int elem = face_to_zones2[face][j];
         face_to_zones.push_back(elem);
      }
   }
   // Output reordered zone and face connectivity in to blueprint format
   for (int newface = 0; newface < nfaces; ++newface)
   {
      // reordered
      int ncorners = face_to_corners2_reord[newface].size();
      //face_to_ncorners[face] = ncorners;
      face_to_corners_reord.push_back(ncorners);
      for (int j = 0; j < ncorners; ++j)
      {
         int c = face_to_corners2_reord[newface][j];
         face_to_corners_reord.push_back(c);
      }
      int nelem_on_face = face_to_zones2_reord[newface].size();
      face_to_zones_reord.push_back(nelem_on_face);
      for (int j = 0; j < nelem_on_face; ++j)
      {
         int elem = face_to_zones2_reord[newface][j];
         face_to_zones_reord.push_back(elem);
      }
   }

   std::vector<std::vector<double>> shared_face_verts;
   int nsfaces = pmesh->GetNSharedFaces();
   shared_face_verts.resize(nsfaces);

   // Divide up the mesh boundary faces in to groups with
   // the same boundary condition. Note that this also includes
   // shared (processor boundary) faces.
   std::map<int, std::vector<int>> boundaries;

   // Find first the non-shared boundary faces
   for (int face = 0; face < nfaces; ++face)
   {
      bool is_bndry = (face_to_bndry[face] >= 0);
      // If a physical boundary face, get the boundary condition
      if (is_bndry)
      {
         int i = face_to_bndry[face];
         int bdnry_cond_id = pmesh->GetBdrAttribute(i);
         if (boundaries.find(bdnry_cond_id) == boundaries.end())
         {
            //int newface = face_to_newface_mapping[face];
            boundaries[bdnry_cond_id] = {face};
         }
         else
         {
            boundaries[bdnry_cond_id].push_back(face);
         }
      }
   }

   // Tag each face with a boundary condition ID (i.e., compute
   //  face_to_bcid[face] = bcID for each mesh face).
   // The BC IDs are reenumerated from those encoded in the original mesh
   // by 1) ordering with respect to reflecting, vaccuum, and source
   //    2) also adding BC IDs corresponding to shared faces
   // For non-boundary faces, face_to_bcid[face] = 0.

   std::vector<int> boundaries_types = {0, 0, 0, 0};
   std::vector<int> boundary_conditions;
   std::vector<std::vector<int>> gr_sfaces;
   int nbelem_corner_faces = 0;
   ComputeFaceIDs(pmesh,
                  face_to_ncorners,
                  boundaries,
                  nbelem_corner_faces,
                  gr_sfaces,
                  boundaries_types,
                  boundary_conditions,
                  face_to_bcid);
   // reorder face_to_bcid
   for (int face = 0; face < nfaces; ++face)
   {
      int newface = face_to_newface_mapping[face];
      face_to_bcid_reord[newface] = face_to_bcid[face];
   }

   // Compute the shared faces arrays

   std::vector<int> shared_faces;
   std::vector<int> shared_zone_and_face_pairs;
   std::vector<std::vector<int>> shared_faces_all_groups;
   std::vector<std::vector<int>> shared_zone_and_face_pairs_all_groups;

   int num_face_nbrs = pmesh->GetNFaceNeighbors();
   for (int fn = 0; fn < num_face_nbrs; fn++)
   {
      int nbr_group = pmesh->face_nbr_group[fn];
      std::vector<int> shared_faces_in_group;
      std::vector<int> shared_zone_and_face_pairs_in_group;
      SortSharedFaces(pmesh,
                      zone_and_gvertex_to_lvertex,
                      zone_and_lvertex_to_corner,
                      face_to_bcid,
                      fn,
                      gr_sfaces,
                      shared_faces_in_group,
                      shared_zone_and_face_pairs_in_group);
      for (int j = 0; j < shared_faces_in_group.size(); ++j)
      {
         int val = shared_faces_in_group[j];
         shared_faces.push_back(val);
      }
      int npairs = shared_zone_and_face_pairs_in_group.size() / 2;
      for (int j = 0; j < npairs; ++j)
      {
         int elem = shared_zone_and_face_pairs_in_group[2 * j];
         int face = shared_zone_and_face_pairs_in_group[2 * j + 1];
         shared_zone_and_face_pairs.push_back(elem);
         shared_zone_and_face_pairs.push_back(face);
      }
      shared_zone_and_face_pairs_all_groups.push_back(shared_zone_and_face_pairs_in_group);
      shared_faces_all_groups.push_back(shared_faces_in_group);
   }
   CheckSharedFaceData(pmesh,
                       corner_to_node_x,
                       corner_to_node_y,
                       corner_to_node_z,
                       shared_zone_and_face_pairs_all_groups,
                       shared_faces_all_groups);

   // Put mesh info in the conduit mesh
   mParametersNode["boundaries/boundary_conditions"].set(boundary_conditions.data(), boundary_conditions.size());
   mParametersNode["boundaries/boundaries_types"].set(boundaries_types.data(), boundaries_types.size());

   int ndim = pmesh->Dimension();

   // CONDUIT OUTPUT
   mParametersNode["size/ndim"] = ndim;
   mParametersNode["size/rank"] = pmesh->GetMyRank();
   mParametersNode["size/nzones"] = nelems;
   mParametersNode["size/ncornr"] = ncorners_total;
   mParametersNode["size/nverts"] = pmesh->GetNV();
   mParametersNode["size/nsides"] = ncorners_total;
   mParametersNode["size/nbelem"] = nbelem_corner_faces;
   nbelem = nbelem_corner_faces;
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

   // This accounts for the fact that some "MFEM face neighbors"
   // might have shared vertices but no shared faces
   int num_face_nbrs_teton = 0;
   int nsfaces_teton = 0;
   for (int fn = 0; fn < gr_sfaces.size(); ++fn)
   {
      if (gr_sfaces[fn].size() > 0)
      {
         num_face_nbrs_teton += 1;
         nsfaces_teton += gr_sfaces[fn].size();
      }
   }
   mParametersNode["size/ncomm"] = num_face_nbrs_teton;

   // TODO change hard-coding of this
   mParametersNode["size/nangsn"] = 8;
   mParametersNode["size/enableNLTE"] = 0;
   mParametersNode["size/methodNLTE"] = 1;
   mParametersNode["size/functionRNLTE"] = 2;
   mParametersNode["size/tfloor"] = 2.5e-05;
   mParametersNode["size/radForceMultiplier"] = 1.0;
   mParametersNode["size/betaNLTE"] = 0.2;
   mParametersNode["size/gammaNLTE"] = 4.0;
   mParametersNode["size/planckCutoffNLTE1"] = 20.0;
   mParametersNode["size/planckCutoffNLTE2"] = 30.0;
   mParametersNode["size/DopplerShiftOn"] = 1;
   mParametersNode["size/useNewNonLinearSolver"] = 0;
   mParametersNode["size/useNewGTASolver"] = 0;
   mParametersNode["size/usePWLD"] = 0;
   mParametersNode["size/useSurfaceMassLumping"] = 0;
   mParametersNode["size/useGPU"] = 0;
   mParametersNode["size/useCUDASolver"] = 0;
   mParametersNode["size/useCUDASweep"] = 0;
   mParametersNode["size/zoneBatchSize"] = 500;
   mParametersNode["size/nConcurrentBatches"] = 3;

   // GPU related

   // CONDUIT OUTPUT

   mMeshNode["shared_boundaries/nsfaces"] = n_shared_corner_faces;
   mParametersNode["shared_boundaries/nsfaces"] = n_shared_corner_faces;
   std::vector<int> shared_faces_reord, shared_zonefacepairs_reord;
   // Relabel the face indices using the new face ordering
   ReorderSharedFaces(ndim, face_to_newface_mapping, shared_faces, shared_faces_reord, shared_zonefacepairs_reord);
   if (n_shared_corner_faces > 0)
   {
      mMeshNode["shared_boundaries/shared_zone_and_face_pairs"].set(shared_zone_and_face_pairs.data(),
                                                                    shared_zone_and_face_pairs.size());
      mMeshNode["shared_boundaries/shared_faces"].set(shared_faces.data(), shared_faces.size());
   }

   mMeshNode["coordsets"]["coords"]["type"] = "explicit";
   conduit::Node &mesh_faces = mMeshNode["topologies/mesh/faces"];
   conduit::Node &mesh_elements = mMeshNode["topologies/mesh/elements"];

   mesh_faces["face_to_corners"].set(face_to_corners.data(), face_to_corners.size());
   mesh_faces["face_to_zones"].set(face_to_zones.data(), face_to_zones.size());
   mesh_faces["face_to_bcid"].set(face_to_bcid.data(), face_to_bcid.size());
   mesh_elements["zone_to_faces"].set(zone_to_faces.data(), zone_to_faces.size());

   // set topologies
   mMeshNode["topologies"]["mesh"]["type"] = "unstructured";
   mMeshNode["topologies"]["mesh"]["coordset"] = "coords";
   mMeshNode["topologies"]["mesh"]["elements"]["shape"] = "polygonal";

   mesh_elements["zone_to_corners"].set(zone_to_corners.data(), zone_to_corners.size());

   // corner to vertex connectivity array
   if (ndim == 2)
   {
      swap_axis_ordering = true;
   }
   else
   {
      swap_axis_ordering = false;
   }
   if (swap_axis_ordering)
   {
      mMeshNode["topologies"]["mesh"]["corners"]["corner_to_node_x"].set(corner_to_node_y.data(),
                                                                         corner_to_node_y.size());
      mMeshNode["topologies"]["mesh"]["corners"]["corner_to_node_y"].set(corner_to_node_x.data(),
                                                                         corner_to_node_x.size());
   }
   else
   {
      mMeshNode["topologies"]["mesh"]["corners"]["corner_to_node_x"].set(corner_to_node_x.data(),
                                                                         corner_to_node_x.size());
      mMeshNode["topologies"]["mesh"]["corners"]["corner_to_node_y"].set(corner_to_node_y.data(),
                                                                         corner_to_node_y.size());
   }

   if (ndim > 2)
   {
      mMeshNode["topologies"]["mesh"]["corners"]["corner_to_node_z"].set(corner_to_node_z.data(),
                                                                         corner_to_node_z.size());
   }

   // corner to zone connectivity array
   mMeshNode["topologies"]["mesh"]["corners"]["corner_to_zone"].set(corner_to_zone.data(), corner_to_zone.size());
   mMeshNode["topologies"]["mesh"]["corners"]["corner_to_vertex"].set(corner_to_vertex.data(), corner_to_vertex.size());

   OutputConduitParameters(
       pmesh, density, heat_capacity, rad_temp, material_temp, gr_bounds, abs_opacity, scat_opacity, electron_density);

   // This transforms the connectivity data in "blueprint" format to
   // the format that standalone Teton expects
   ComputeTetonMeshConnectivityArray();
   CreateTetonMeshCornerCoords();

   int rank = pmesh->GetMyRank();
   std::string file_name = "mesh_rank" + std::to_string(rank) + ".json";
   conduit::relay::io::save(mMeshNode, file_name);
   file_name = "parameters_rank" + std::to_string(rank) + ".json";
   conduit::relay::io::save(mParametersNode, file_name);

}

// Assumes 2D for now
// Note: faces correspond to {v1,v2}, {v2,v3}, {v3,v4}, and {v4,v1}
void TetonBlueprint::OrderVertices(mfem::ParMesh *pmesh,
                                   int elem,
                                   std::vector<int> &vertices_ordered,
                                   std::vector<int> &faces_ordered)
{
   int dim = pmesh->Dimension();
   mfem::Array<int> fcs, cor;

   mfem::Array<int> elem_verts;
   pmesh->GetElementVertices(elem, elem_verts);
   int nverts = elem_verts.Size();

   pmesh->GetElementEdges(elem, fcs, cor);
   int nfaces = fcs.Size();
   std::set<int> faces_left;
   for (int j = 0; j < nfaces; ++j)
   {
      faces_left.insert(fcs[j]);
   }

   // Order the vertices consecutively
   // Get one vertex making up the first edge
   int face = fcs[0];
   faces_ordered.push_back(face);
   mfem::Array<int> vert;
   pmesh->GetEdgeVertices(face, vert);
   int v = vert[0];
   int v_next = vert[1];
   vertices_ordered.push_back(v);
   faces_left.erase(face);
   while (vertices_ordered.size() < nverts)
   {
      std::set<int>::iterator itr;
      for (itr = faces_left.begin(); itr != faces_left.end(); ++itr)
      {
         int face = *itr;
         mfem::Array<int> vert;
         pmesh->GetEdgeVertices(face, vert);
         int v1 = vert[0];
         int v2 = vert[1];
         if (v1 == v_next)
         {
            vertices_ordered.push_back(v1);
            v_next = v2;
            faces_left.erase(face);
            faces_ordered.push_back(face);
            break;
         }
         else if (v2 == v_next)
         {
            vertices_ordered.push_back(v2);
            v_next = v1;
            faces_left.erase(face);
            faces_ordered.push_back(face);
            break;
         }
      }
   }

   // Fix orientation to be counter-clockwise
   // Look at the three points {v1,v2}, {v2,v3}, {v3,v4}
   // See https://www.geeksforgeeks.org/orientation-3-ordered-points/

   // First check whether the vertices are ordered clockwise or counter-clockwise
   std::vector<std::vector<double>> v_coords(nverts);
   for (int j = 0; j < 3; ++j)
   {
      int v = vertices_ordered[j];
      double vert_coords[3];
      pmesh->GetNode(v, vert_coords);
      v_coords[j].push_back(vert_coords[0]);
      v_coords[j].push_back(vert_coords[1]);
   }
   double p1x = v_coords[0][0];
   double p1y = v_coords[0][1];
   double p2x = v_coords[1][0];
   double p2y = v_coords[1][1];
   double p3x = v_coords[2][0];
   double p3y = v_coords[2][1];
   double val = (p2y - p1y) * (p3x - p2x) - (p2x - p1x) * (p3y - p2y);
   bool is_clockwise;
   if (val == 0)
   {
      std::cerr << "degenerate quad!" << std::endl;
      std::cerr << "  elem: " << elem << std::endl;
      exit(1);
   }
   else if (val > 0)
   {
      is_clockwise = true;
   }
   else
   {
      is_clockwise = false;
   }

   if (is_clockwise)
   {
      // If oriented counter-clockwise, reverse orientation
      // Want faces {v1,v2}, {v2,v3}, {v3,v4}, {v4,v1} to go to faces
      // {v1,v4}, {v4,v3},{v3,v2}, {v2,v1} and vertices {v1,v2,v3,v4} to go to
      // {v1,v4,v3,v2}
      std::vector<int> vertices_ordered2(nverts);
      vertices_ordered2[0] = vertices_ordered[0];
      for (int j = 0; j < nverts - 1; ++j)
      {
         vertices_ordered2[j + 1] = vertices_ordered[nverts - 1 - j];
      }
      std::vector<int> faces_ordered2(nfaces);
      for (int j = 0; j < nfaces; ++j)
      {
         faces_ordered2[j] = faces_ordered[nfaces - 1 - j];
      }
      // Overwrite old data
      for (int j = 0; j < nverts; ++j)
      {
         vertices_ordered[j] = vertices_ordered2[j];
      }
      for (int j = 0; j < nfaces; ++j)
      {
         faces_ordered[j] = faces_ordered2[j];
      }
   }

   // Finally, choose the starting vertex to be the one with the smallest
   // x coordinate. In case of a tie, choose the one with the largest y coordinate.
   // Then start counter clock-wise ordering at this node.

   // First find minimum x coordinate
   std::vector<int> vertices_ordered2(nverts);
   std::vector<int> faces_ordered2(nfaces);
   double minx = 1e10;
   std::vector<int> vertids_minxs;
   for (int j = 0; j < nverts; ++j)
   {
      int v = vertices_ordered[j];
      double vert_coords[3];
      pmesh->GetNode(v, vert_coords);
      double px = vert_coords[0];
      if (px < minx)
      {
         minx = px;
      }
   }
   // Next find vertex IDs with this minimum x coordinate
   for (int j = 0; j < nverts; ++j)
   {
      int v = vertices_ordered[j];
      double vert_coords[3];
      pmesh->GetNode(v, vert_coords);
      double px = vert_coords[0];
      if (px == minx)
      {
         vertids_minxs.push_back(j);
      }
   }

   // Now break ties by choosing largest y coordinate
   double maxy = -1e10;
   int vertid_first; // local (to element) vertex ID
   int vfirst;
   for (int j = 0; j < vertids_minxs.size(); ++j)
   {
      int j1 = vertids_minxs[j];
      int v = vertices_ordered[j1];
      double vert_coords[3];
      pmesh->GetNode(v, vert_coords);
      double py = vert_coords[1];
      if (py > maxy)
      {
         maxy = py;
         vfirst = v;
      }
   }
   for (int j = 0; j < nverts; ++j)
   {
      int v = vertices_ordered[j];
      if (v == vfirst)
      {
         vertid_first = j;
         break;
      }
   }

   // Arrange the vertex IDs counter-clockwise starting at the vertex
   // selected above. For example, say that the first coordinate is v3.
   // Then we want {v1,v2,v3,v4} -> {v3,v4,v1,v2} and
   // {v1,v2}, {v2,v3}, {v3,v4}, {v4,v1} ->
   // {v3,v4}, {v4,v1}, {v1,v2}, {v2,v3}

   for (int j = vertid_first; j < nverts; ++j)
   {
      vertices_ordered2[j - vertid_first] = vertices_ordered[j];
      faces_ordered2[j - vertid_first] = faces_ordered[j];
   }
   for (int j = 0; j < vertid_first; ++j)
   {
      vertices_ordered2[nverts - vertid_first + j] = vertices_ordered[j];
      faces_ordered2[nverts - vertid_first + j] = faces_ordered[j];
   }
   for (int j = 0; j < nverts; ++j)
   {
      vertices_ordered[j] = vertices_ordered2[j];
   }
   for (int j = 0; j < nfaces; ++j)
   {
      faces_ordered[j] = faces_ordered2[j];
   }

   // Finally, update the data on previous vertices encountered,
   // as well as what the last vertex is

   vert_prev_global = vertices_ordered[nverts - 1];
   for (int j = 0; j < nverts; ++j)
   {
      int v = vertices_ordered[j];
      vert_to_counters[v].push_back(global_vert_counter);
      global_vert_counter += 1;
   }
}

// Note: edges correspond to {v1,v2}, {v2,v3}, {v3,v4}, and {v4,v1}
void TetonBlueprint::OrderVertices2D(mfem::ParMesh *pmesh, int elem, int elem_face, std::vector<int> &vertices_ordered)
{
   mfem::Array<int> fcs, cor;
   mfem::Array<int> face_verts;
   pmesh->GetFaceVertices(elem_face, face_verts);
   int nverts = face_verts.Size();

   pmesh->GetFaceEdges(elem_face, fcs, cor);

   // Confusingly, I'm using "faces" in OrderVertices2D for the elem_face edges
   int nfaces = fcs.Size();
   std::set<int> faces_left;
   for (int j = 0; j < nfaces; ++j)
   {
      faces_left.insert(fcs[j]);
   }

   std::vector<int> faces_ordered;
   // Order the vertices consecutively
   // Get one vertex making up the first edge
   int face = fcs[0];
   faces_ordered.push_back(face);
   mfem::Array<int> vert;
   pmesh->GetEdgeVertices(face, vert);
   int v = vert[0];
   int v_next = vert[1];
   vertices_ordered.push_back(v);
   faces_left.erase(face);
   while (vertices_ordered.size() < nverts)
   {
      std::set<int>::iterator itr;
      for (itr = faces_left.begin(); itr != faces_left.end(); ++itr)
      {
         int face = *itr;
         mfem::Array<int> vert;
         pmesh->GetEdgeVertices(face, vert);
         int v1 = vert[0];
         int v2 = vert[1];
         if (v1 == v_next)
         {
            vertices_ordered.push_back(v1);
            v_next = v2;
            faces_left.erase(face);
            faces_ordered.push_back(face);
            break;
         }
         else if (v2 == v_next)
         {
            vertices_ordered.push_back(v2);
            v_next = v1;
            faces_left.erase(face);
            faces_ordered.push_back(face);
            break;
         }
      }
   }

   // Fix orientation to be clockwise

   // compute a normal to the surface
   int v1 = vertices_ordered[0];
   double p1_coords[3];
   pmesh->GetNode(v1, p1_coords);
   int v2 = vertices_ordered[1];
   double p2_coords[3];
   pmesh->GetNode(v2, p2_coords);
   int v4 = vertices_ordered[3];
   double p4_coords[3];
   pmesh->GetNode(v4, p4_coords);
   double a1 = p2_coords[0] - p1_coords[0];
   double a2 = p2_coords[1] - p1_coords[1];
   double a3 = p2_coords[2] - p1_coords[2];
   double b1 = p4_coords[0] - p1_coords[0];
   double b2 = p4_coords[1] - p1_coords[1];
   double b3 = p4_coords[2] - p1_coords[2];
   mfem::Vector norm_vec(3);
   norm_vec(0) = (a2 * b3 - a3 * b2);
   norm_vec(1) = (a3 * b1 - a1 * b3);
   norm_vec(2) = (a1 * b2 - a2 * b1);
   // Compute a vector pointing from the element centroid to
   // the face centroid.
   // To do so, first find the centroid of the element
   mfem::Vector pc(3);
   mfem::Array<int> elem_verts;
   pmesh->GetElementVertices(elem, elem_verts);
   double pc1 = 0, pc2 = 0, pc3 = 0;
   int nverts_elem = elem_verts.Size();
   for (int j = 0; j < nverts_elem; ++j)
   {
      int v = elem_verts[j];
      double p_coords[3];
      pmesh->GetNode(v, p_coords);
      pc1 += p_coords[0] / nverts_elem;
      pc2 += p_coords[1] / nverts_elem;
      pc3 += p_coords[2] / nverts_elem;
   }
   // Find the centroid of the face
   double pa1 = 0, pa2 = 0, pa3 = 0;
   int nverts_face = vertices_ordered.size();
   for (int j = 0; j < nverts_face; ++j)
   {
      int v = vertices_ordered[j];
      double p_coords[3];
      pmesh->GetNode(v, p_coords);
      pa1 += p_coords[0] / nverts_face;
      pa2 += p_coords[1] / nverts_face;
      pa3 += p_coords[2] / nverts_face;
   }
   // Find (p_a - p_c) \cdot n
   mfem::Vector outward_vec(3);
   outward_vec(0) = pa1 - pc1;
   outward_vec(1) = pa2 - pc2;
   outward_vec(2) = pa3 - pc3;

   // If (p_a - p_c) \cdot n < 0, reverse orientation
   bool is_counter_clockwise;
   if (outward_vec * norm_vec > 0)
   {
      is_counter_clockwise = true;
   }
   else
   {
      is_counter_clockwise = false;
   }
   if (is_counter_clockwise)
   {
      // If oriented counter-clockwise, reverse orientation
      // Want faces {v1,v2}, {v2,v3}, {v3,v4}, {v4,v1} to go to faces
      // {v1,v4}, {v4,v3},{v3,v2}, {v2,v1} and vertices {v1,v2,v3,v4} to go to
      // {v1,v4,v3,v2}
      int nverts = vertices_ordered.size();
      std::vector<int> vertices_ordered2(nverts);
      vertices_ordered2[0] = vertices_ordered[0];
      for (int j = 0; j < nverts - 1; ++j)
      {
         vertices_ordered2[j + 1] = vertices_ordered[nverts - 1 - j];
      }
      // Overwrite old data
      for (int j = 0; j < nverts; ++j)
      {
         vertices_ordered[j] = vertices_ordered2[j];
      }
   }

}

void TetonBlueprint::OrderFaces3D(mfem::ParMesh *pmesh, int elem, std::vector<int> &faces_ordered)
{
   mfem::Array<int> fcs, cor;
   pmesh->GetElementFaces(elem, fcs, cor);

   faces_ordered.resize(fcs.Size());
   for (int j = 0; j < fcs.Size(); j++)
   {
      faces_ordered[j] = fcs[j];
   }

   // Get starting vertex in first face.
   // First get min x coordinate
   mfem::Array<int> face_verts;
   pmesh->GetFaceVertices(faces_ordered[0], face_verts);
   v_start_glob = face_verts[0];
   double vert_coords[3];
   pmesh->GetNode(v_start_glob, vert_coords);
   // First get min x coordinate
   double x_min = vert_coords[0];
   for (int j = 1; j < face_verts.Size(); j++)
   {
      int v = face_verts[j];
      double vert_coords[3];
      pmesh->GetNode(v, vert_coords);
      if (x_min > vert_coords[0])
      {
         x_min = vert_coords[0];
      }
   }
   // Now get the vertex with the minimum x coordinate and
   // (in case of a tie) the maximum y coordinate
   double y_max = -1.0e10;
   for (int j = 0; j < face_verts.Size(); j++)
   {
      int v = face_verts[j];
      double vert_coords[3];
      pmesh->GetNode(v, vert_coords);
      if ((std::abs(x_min - vert_coords[0]) < 1.0e-15) && (y_max < vert_coords[1]))
      {
         y_max = vert_coords[1];
         v_start_glob = v;
      }
      else
      {
      }
   }
}

void TetonBlueprint::OrderVertices3D(mfem::ParMesh *pmesh,
                                     int elem,
                                     std::vector<int> &vertices_ordered,
                                     std::vector<int> &faces_ordered,
                                     std::vector<std::vector<int>> &face_verts_ordered)
{
   mfem::Array<int> fcs, cor;
   pmesh->GetElementFaces(elem, fcs, cor);
   face_verts_ordered.resize(fcs.Size());

   // Set counter on all edges to zero
   v1_last_glob = -1;
   v2_last_glob = -1;
   std::map<std::set<int>, int>::iterator it = elem_edges_count.begin();
   while (it != elem_edges_count.end())
   {
      elem_edges_count.erase(it);
      it++;
   }
   elem_edges_count.clear();
   mfem::Array<int> edges;
   pmesh->GetElementEdges(elem, edges, cor);
   for (int j = 0; j < fcs.Size(); j++)
   {
      mfem::Array<int> vert;
      pmesh->GetEdgeVertices(edges[j], vert);
      std::set<int> edge;
      edge.insert(vert[0]);
      edge.insert(vert[1]);
      elem_edges_count[edge] = 0;
   }

   // Get face ordering and starting vertex, v_start_glob
   OrderFaces3D(pmesh, elem, faces_ordered);

   // Get the vertex ordering for each face
   for (int j = 0; j < faces_ordered.size(); j++)
   {
      int face = faces_ordered[j];
      nfaces_count_glob = j + 1;
      OrderVertices2D(pmesh, elem, face, face_verts_ordered[j]);
   }

   // Vertices are ordered to agree with other mesh sources (to the extent possible)
   int nzonefaces = fcs.Size();
   for (int k = 0; k < 4; ++k)
   {
      int v = face_verts_ordered[0][k];
      vertices_ordered.push_back(v);
   }
   int v = face_verts_ordered[nzonefaces - 1][0];
   vertices_ordered.push_back(v);
   for (int k = 3; k > 0; --k)
   {
      int v = face_verts_ordered[nzonefaces - 1][k];
      vertices_ordered.push_back(v);
   }
}

void TetonBlueprint::ComputeFaceIDs(mfem::ParMesh *pmesh,
                                    std::vector<int> &face_to_ncorners,
                                    std::map<int, std::vector<int>> &boundaries,
                                    int &nbelem_corner_faces,
                                    std::vector<std::vector<int>> &gr_sface,
                                    std::vector<int> &boundaries_types,
                                    std::vector<int> &boundary_conditions,
                                    std::vector<int> &face_to_bcid)
{
   // create boundary conditions info
   // This is Teton's scheme for associating an
   // integer with a B.C. type
   // bcType_none    = 31
   // bcType_refl    = 32
   // bcType_shared  = 33
   // bcType_temp    = 34
   // bcType_vac     = 35
   // bcType_fds     = 36
   // bcType_invalid = 0
   //std::map<int, int> boundary_id_to_type;
   std::map<int, int> boundary_id_to_ncornerfaces;
   std::map<int, int> number_boundary_types;
   for (int j = 31; j < 37; ++j)
      number_boundary_types[j] = 0;
   //FormBoundaryIDToBCType(boundary_id_to_type);
   std::map<int, std::vector<int>>::iterator itr;

   int nfaces = pmesh->GetNumFaces();
   int ndim = pmesh->Dimension();

   // First tag each face with the original boundary condition ID from the mesh.
   // This will later be changed to the re-enumerated boundary condition ID (so
   // each ID is consecutive and reordered by boundary condition type).
   face_to_bcid.resize(nfaces);
   for (int face = 0; face < nfaces; ++face)
      face_to_bcid[face] = 0;

   for (itr = boundaries.begin(); itr != boundaries.end(); ++itr)
   {
      int bc_id = itr->first;
      int bc_type = (boundary_id_to_type)[bc_id];
      number_boundary_types[bc_type] += 1;
      boundary_id_to_ncornerfaces[bc_id] = 0;
      for (int j = 0; j < (itr->second).size(); ++j)
      {
         int face = (itr->second)[j];
         int ncorner_faces = face_to_ncorners[face];
         boundary_id_to_ncornerfaces[bc_id] += ncorner_faces;
         face_to_bcid[face] = bc_id;
      }
   }

   pmesh->ExchangeFaceNbrData();

   // Determine the number of "MFEM face neighbors" that have non-zero shared faces
   // Note that (IS THIS TRUE???) two subdomains (each associated with an MPI rank) can share
   // vertices but no faces, and so might be "MFEM face neighbors" without shared faces
   int num_face_nbrs = pmesh->GetNFaceNeighbors();
   int ngroups = pmesh->GetNGroups();
   CreateSharedFacesTable(pmesh, gr_sface);
   int num_face_nbrs_for_teton = 0;
   for (int g = 0; g < gr_sface.size(); g++)
   {
      int num_sfaces = gr_sface[g].size();
      if (num_sfaces > 0)
         num_face_nbrs_for_teton += 1;
   }

   boundaries_types[0] = number_boundary_types[32];                             // number of reflecting
   boundaries_types[1] = number_boundary_types[35];                             // number of vaccuum
   boundaries_types[2] = number_boundary_types[34] + number_boundary_types[36]; // number of source
   boundaries_types[3] = num_face_nbrs_for_teton;                               // number of shared

   // Order the boundaries by: reflecting, vaccuum, source, shared.

   // First add the non-shared boundaries
   int num_nonshared_bndrs = boundaries_types[0] + boundaries_types[1] + boundaries_types[2];
   int num_bndrs = num_nonshared_bndrs + boundaries_types[3];
   int nreflec, nvaccuum, nsource;
   nreflec = boundaries_types[0];
   nvaccuum = boundaries_types[1];
   nsource = boundaries_types[2];
   std::vector<int> bc_ids_ordered(num_bndrs);
   int jreflect = 0, jvaccuum = 0, jsource = 0, jshared = 0;
   for (itr = boundaries.begin(); itr != boundaries.end(); ++itr)
   {
      int bc_id = itr->first;
      int bc_type = (boundary_id_to_type)[bc_id];
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
         int index = jsource + nreflec + nvaccuum;
         bc_ids_ordered[index] = bc_id;
         jsource += 1;
      }
   }

   // Now compute the shared boundaries info
   // First we get the max of the non-shared bc_ids so that we
   // can assign a unique bc_id
   int bc_id_max = 0;
   for (itr = boundaries.begin(); itr != boundaries.end(); ++itr)
   {
      bc_id_max = std::max(itr->first, bc_id_max);
   }
   std::map<int, int> bc_id_to_rank;

   // Tag each shared face in the same group with its boundary ID
   int fn_counter_teton = 0;
   for (int fn = 0; fn < num_face_nbrs; fn++)
   {
      int num_sfaces = gr_sface[fn].size();
      // If there aren't any shared faces for this "MFEM face neighbor", don't
      // do anything
      if (num_sfaces < 1)
      {
         continue;
      }

      // Store rank of MPI process that shares this face
      int nbr_rank = pmesh->GetFaceNbrRank(fn);
      // this gives a unique boundary id for each shared boundary
      int bc_id = fn_counter_teton + bc_id_max + 1;
      int index = fn_counter_teton + num_nonshared_bndrs;
      // We want to increment this only when num_sfaces > 0
      fn_counter_teton += 1;

      bc_ids_ordered[index] = bc_id;
      boundary_id_to_type[bc_id] = 33; // shared boundary
      bc_id_to_rank[bc_id] = nbr_rank;
      boundary_id_to_ncornerfaces[bc_id] = 0;

      // tag each shared face in this shared boundary with its b.c. ID
      for (int i = 0; i < num_sfaces; i++)
      {
         //shared_face_bcid[sface[i]] = fn + num_nonshared_bndrs;
         //int face = pmesh -> GetSharedFace(sface[i]);
         int face = pmesh->GetSharedFace(gr_sface[fn][i]);
         int ncorners = face_to_ncorners[face];
         boundary_id_to_ncornerfaces[bc_id] += ncorners;
         face_to_bcid[face] = bc_id;

         // TODO: figure out the logic here for n_shared_corner_faces
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
            std::cerr << "1D is not yet implemented! " << std::endl;
         }
         // TODO: figure out the logic here for n_shared_corner_faces
      }
   }

   // Create boundary conditions array for Teton
   nbelem_corner_faces = 0;
   boundary_conditions.push_back(num_bndrs);
   std::map<int, int> bc_id_to_new_bc_id;
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int bc_type = boundary_id_to_type[bc_id];
      boundary_conditions.push_back(bc_type);
      bc_id_to_new_bc_id[bc_id] = j;
   }
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int ncorner_faces = boundary_id_to_ncornerfaces[bc_id];
      boundary_conditions.push_back(ncorner_faces);
      nbelem_corner_faces += ncorner_faces;
   }
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int bc_type = boundary_id_to_type[bc_id];
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
      int bc_type = (boundary_id_to_type)[bc_id];
      bc_types.push_back(bc_type);
   }
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int ncorner_faces = boundary_id_to_ncornerfaces[bc_id];
      num_bc_cornerfaces.push_back(ncorner_faces);
      nbelem_corner_faces += ncorner_faces;
   }
   for (int j = 0; j < num_bndrs; ++j)
   {
      int bc_id = bc_ids_ordered[j];
      int bc_type = (boundary_id_to_type)[bc_id];
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
   // CONDUIT OUTPUT
}

void TetonBlueprint::SortSharedFaces(mfem::ParMesh *pmesh,
                                     std::map<std::pair<int, int>, int> &zone_and_gvertex_to_lvertex,
                                     std::vector<std::vector<int>> &zone_and_lvertex_to_corner,
                                     std::vector<int> &face_to_bcid,
                                     int fn,
                                     std::vector<std::vector<int>> &gr_sfaces,
                                     std::vector<int> &shared_faces,
                                     std::vector<int> &shared_zone_and_face_pairs)
{
   // Compute the centroid for each shared faces
   int nsfaces = gr_sfaces[fn].size();
   std::vector<std::vector<double>> shared_face_verts(nsfaces);
   std::map<std::vector<double>, int> coords_to_sface;
   int dim = pmesh->Dimension();
   for (int j = 0; j < nsfaces; ++j)
   {
      int sface = gr_sfaces[fn][j];
      // average the face vertex coordinates
      mfem::Array<int> vert;
      int face = pmesh->GetSharedFace(sface);
      pmesh->GetFaceVertices(face, vert);
      double vx = 0, vy = 0, vz = 0;
      int nvert = vert.Size();
      for (int ivert = 0; ivert < nvert; ++ivert)
      {
         double v[3];
         pmesh->GetNode(vert[ivert], v);
         vx += v[0];
         vy += v[1];
         if (dim == 3)
         {
            vz += v[2];
         }
      }
      vx /= nvert;
      vy /= nvert;
      vz /= nvert;
      // get the shared face id and append "face centroid"
      shared_face_verts[j].resize(dim);
      shared_face_verts[j][0] = vx;
      shared_face_verts[j][1] = vy;
      if (dim == 3)
         shared_face_verts[j][2] = vz;
      coords_to_sface[shared_face_verts[j]] = sface;
   }

   // Sort the shared faces based on their centroids
   std::sort(shared_face_verts.begin(), shared_face_verts.end(), CompareCoords);

   // Sort the corners on each shared face by their corresponding vertices
   int nsfaces_total = 0;
   for (int j_teton = 0; j_teton < nsfaces; ++j_teton)
   {
      // index for next shared face of sorted shared faces
      std::vector<double> sface_coords(dim);
      sface_coords[0] = shared_face_verts[j_teton][0];
      sface_coords[1] = shared_face_verts[j_teton][1];
      if (dim == 3)
         sface_coords[2] = shared_face_verts[j_teton][2];
      int sface = coords_to_sface[sface_coords];
      int face = pmesh->GetSharedFace(sface);
      int elem1, elem2;
      pmesh->GetFaceElements(face, &elem1, &elem2);
      int bcID = face_to_bcid[face];
      // get the face vertex coordinates
      mfem::Array<int> vert;
      pmesh->GetFaceVertices(face, vert);
      int nvert = vert.Size();
      // shared_face_corners[cID] = {corner_to_vertex_coords[cID], cID}
      std::vector<std::vector<double>> shared_face_corners(nvert);
      std::map<std::vector<double>, int> coords_to_corner;
      for (int ivert = 0; ivert < nvert; ++ivert)
      {
         int v = vert[ivert];
         int j = zone_and_gvertex_to_lvertex[{elem1, v}];
         int c1 = zone_and_lvertex_to_corner[elem1][j];
         shared_face_corners[ivert].resize(dim);
         double coords[3];
         pmesh->GetNode(vert[ivert], coords);
         shared_face_corners[ivert][0] = coords[0];
         shared_face_corners[ivert][1] = coords[1];
         if (dim == 3)
         {
            shared_face_corners[ivert][2] = coords[2];
         }
         // store the relationship (corner coordinates) <--> (local vertex ID)
         std::vector<double> coords2(dim);
         coords2[0] = coords[0];
         coords2[1] = coords[1];
         if (dim == 3)
            coords2[2] = coords[2];
         coords_to_corner[coords2] = c1;
      }
      // Sort the corners based on the corresponding vertex coords.
      // In 3D, each corner-face is represented via {bcID, elem, face, c},
      // where bcID == boundary condition ID, elem = element the corner-face belongs to,
      // face == face the corner-face belongs to, and c == corner the corner face
      // belongs to. In 2D, the corner faces are stored as {bcID, elem, face, c1, c2} ,
      // where c1 and c2 are the two corners on each face.
      std::sort(shared_face_corners.begin(), shared_face_corners.end(), CompareCoords);
      for (int ivert_new = 0; ivert_new < nvert; ++ivert_new)
      {
         // get the corner ID associated with this vertex
         std::vector<double> coords(dim);
         coords[0] = shared_face_corners[ivert_new][0];
         coords[1] = shared_face_corners[ivert_new][1];
         if (dim == 3)
            coords[2] = shared_face_corners[ivert_new][2];
         int c = coords_to_corner[coords];
         // NOTE: the face, zone, and corner IDs are increment by
         // one in "teton_setsharedface"
         shared_faces.push_back(bcID);
         shared_faces.push_back(elem1);
         shared_faces.push_back(face);
         shared_faces.push_back(c);
         nsfaces_total += 1;
         if (dim == 2) // store pairs of corners for each shared face
         {
            ivert_new++;
            // get corner ID for next vertex for this shared face
            std::vector<double> coords(dim);
            coords[0] = shared_face_corners[ivert_new][0];
            coords[1] = shared_face_corners[ivert_new][1];
            int c_next = coords_to_corner[coords];
            shared_faces.push_back(c_next);
            nsfaces_total += 1;
         }
         shared_zone_and_face_pairs.push_back(elem1);
         shared_zone_and_face_pairs.push_back(face);
      }
   }
}


void TetonBlueprint::OutputConduitParameters(mfem::ParMesh *pmesh,
                                             mfem::Vector &density,
                                             mfem::Vector &heat_capacity,
                                             mfem::Vector &rad_temp,
                                             mfem::Vector &material_temp,
                                             const mfem::Vector &gr_bounds,
                                             mfem::Vector &abs_opacity,
                                             mfem::Vector &scat_opacity,
                                             mfem::Vector &electron_density)
{
   double dtrad = 3e-06;
   double dtrmn = 1e-40;
   double dtrmx = 0.1;
   double delte = 0.4;
   double deltr = 0.4;
   double tfloor = 2.5e-5;
   int ngr = gr_bounds.Size() - 1;

   // Iteration controls, but use defaults.
   //mParametersNode["iteration/outerMaxIt"] = teton_get_temperature_maxits();
   //mParametersNode["iteration/incidentFluxMaxIt"] = teton_get_fluxexchange_maxits();
   //mParametersNode["iteration/greyMaxIt"] = teton_get_grey_maxits();
   //mParametersNode["iteration/innerNLMaxIt"] = teton_get_nonlinear_maxits();
   //mParametersNode["iteration/outerTempRelTol"] = teton_get_temperature_reltol();
   //mParametersNode["iteration/outerPhiRelTol"] = teton_get_radenergydensity_reltol();
   //mParametersNode["iteration/incidentFluxRelTol"] = teton_get_fluxexchange_reltol();
   //mParametersNode["iteration/greyRelTol"] = teton_get_fluxexchange_reltol();
   //mParametersNode["iteration/innerNLRelTol"] = teton_get_fluxexchange_reltol();
   // Time step controls
   mParametersNode["iteration/dtrad"] = dtrad;
   mParametersNode["iteration/dtrmn"] = dtrmn;
   mParametersNode["iteration/dtrmx"] = dtrmx;
   mParametersNode["iteration/delte"] = delte;
   mParametersNode["iteration/deltr"] = deltr;
   mParametersNode["iteration/tfloor"] = tfloor;

   // Compton scattering settings.
   mParametersNode["compton/use_internal_compton"] = 1;
   mParametersNode["compton/use_thomas_compton"] = 0;
   mParametersNode["compton/use_boltzmann_compton"] = 1;
   mParametersNode["compton/use_internal_sigma_s"] = 0;
   mParametersNode["compton/use_nka"] = 0;

   // Initial conditions
   int nelem = pmesh->GetNE();
   std::vector<double> scattering_multiplier(nelem);
   std::vector<double> corner_temperature;
   int dim = pmesh->Dimension();
   int ncorner_in_elem;
   if (dim == 2)
   {
      ncorner_in_elem = 4;
   }
   else
   {
      ncorner_in_elem = 8;
   }
   for (int elem = 0; elem < nelem; ++elem)
   {
      int elem_attr = pmesh->GetAttribute(elem);
      scattering_multiplier[elem] = 1.0;
      for (int c = 0; c < ncorner_in_elem; ++c)
      {
         double rhoCvTe = material_temp(elem);
         double rhoCv = heat_capacity(elem);
         double Te = rhoCvTe / rhoCv;
         // Note: Confusingly, Teton calls "normalizeMaterial" which
         // divides the zone-based material temperature
         // (set from fields/electron_temp) by rho*Cv,
         // but the corner_temperature field is not normalized.
         corner_temperature.push_back(Te);
      }
   }
   double radiation_energy = 7.54305739119688e-06;
   mMeshNode["fields/density"].set(density.GetData(), nelem);
   mMeshNode["fields/heat_capacity"].set(heat_capacity.GetData(), nelem);
   mMeshNode["fields/electron_temp"].set(material_temp.GetData(), nelem);
   mMeshNode["fields/rad_temp"].set(rad_temp.GetData(), nelem);
   mMeshNode["fields/electron_density"].set(electron_density.GetData(), nelem);
   mMeshNode["fields/corner_temperature"].set(corner_temperature.data(), corner_temperature.size());
   mMeshNode["fields/scattering_multiplier"].set(scattering_multiplier.data(), nelem);
   mMeshNode["fields/radiation_energy"] = radiation_energy;
   // TODO: change
   mParametersNode["radiation_energy"] = radiation_energy;

   // Opacity
   mMeshNode["fields/absorption_opacity"].set(abs_opacity.GetData(), abs_opacity.Size());
   mMeshNode["fields/scattering_opacity"].set(scat_opacity.GetData(), scat_opacity.Size());

   //Energy groups and SN quadrature info
   int qtype = 1;
   int qorder = 10;
   // NOTE: npolar and nazimu need to have a particular relationship
   // based on the number of threads and other fine-grained parallelism.
   // TODO: DOCUMENT THIS FOR UMT
   int npolar = 4;
   int nazimu = 3;
   int paxis = 1;
   std::string top = "quadrature/";
   mParametersNode[top + "gnu"].set(gr_bounds.GetData(), gr_bounds.Size());
   mParametersNode[top + "qtype"] = qtype;
   mParametersNode[top + "qorder"] = qorder;
   mParametersNode[top + "npolar"] = npolar;
   mParametersNode[top + "nazimu"] = nazimu;
   mParametersNode[top + "paxis"] = paxis;
   mParametersNode[top + "num_groups"] = ngr;
   mParametersNode[top + "gtaorder"] = 2;
   mParametersNode[top + "nSetsMaster"] = -1;
   mParametersNode[top + "nSets"] = 1;

   // TODO: change
   mParametersNode["boundary_edits/numSpectrumAngleBins"] = 1;
   std::vector<double> spectrumAngleBinBoundaryList = {-1.0, 1.0};
   mParametersNode["boundary_edits/spectrumAngleBinBoundaryList"].set(spectrumAngleBinBoundaryList.data(),
                                                                      spectrumAngleBinBoundaryList.size());

   std::vector<double> rad_edits(ngr);
   for (int g = 0; g < ngr; ++g)
   {
      rad_edits[g] = 0.0;
   }
   mParametersNode["boundary_edits/RadPowerEscape"].set(rad_edits.data(), rad_edits.size());
   mParametersNode["boundary_edits/RadPowerIncident"].set(rad_edits.data(), rad_edits.size());
   mParametersNode["boundary_edits/PolarSectorPowerEscape"].set(rad_edits.data(), rad_edits.size());

   mParametersNode["memory_allocator/umpire_host_allocator_id"] = -1;
   mParametersNode["memory_allocator/umpire_device_allocator_id"] = -1;
}

void TetonBlueprint::CreateSharedFacesTable(mfem::ParMesh *pmesh, std::vector<std::vector<int>> &group_sfaces)
{
   mfem::Array<int> lface_to_sface;
   ComputeLfaceToSface(pmesh, lface_to_sface);
   int num_face_nbrs = pmesh->GetNFaceNeighbors();

   int dim = pmesh->Dimension();
   group_sfaces.resize(num_face_nbrs);
   for (int fn = 0; fn < num_face_nbrs; fn++)
   {
      int group = pmesh->face_nbr_group[fn];
      if (dim == 2)
      {
         int nedges = pmesh->GroupNEdges(group);
         for (int i = 0; i < nedges; ++i)
         {
            int edge, o;
            pmesh->GroupEdge(group, i, edge, o);
            int sface = lface_to_sface[edge];
            if (sface > -1)
               group_sfaces[fn].push_back(sface);
         }
      }
      else
      {
         int nsfaces_trias = pmesh->GroupNTriangles(group);
         int nsfaces_quads = pmesh->GroupNQuadrilaterals(group);
         for (int i = 0; i < nsfaces_trias; ++i)
         {
            int face, o;
            pmesh->GroupTriangle(group, i, face, o);
            int sface = lface_to_sface[face];
            if (sface > -1)
               group_sfaces[fn].push_back(sface);
         }
         for (int i = 0; i < nsfaces_quads; ++i)
         {
            int face, o;
            pmesh->GroupQuadrilateral(group, i, face, o);
            int sface = lface_to_sface[face];
            if (sface > -1)
               group_sfaces[fn].push_back(sface);
         }
      }
   }
}

void TetonBlueprint::ReorderSharedFaces(int ndim,
                                        std::vector<int> &face_to_newface_mapping,
                                        std::vector<int> &shared_faces,
                                        std::vector<int> &shared_faces_reord,
                                        std::vector<int> &shared_zonefacepairs_reord)
{
   int nsfaces = shared_faces.size();
   if (ndim == 2)
   { // for each shared face, {zone, face, bcid, corner1, corner2} is stored
      // in the array shared_boundaries/shared_faces
      nsfaces = nsfaces / 5;
   }
   else if (ndim == 3)
   { // for each shared face, {zone, face, bcid, corner1} is stored
      // in the array shared_boundaries/shared_faces
      nsfaces = nsfaces / 4;
   }
   else
   {
      std::cerr << "the 1D case is not yet implemented!" << std::endl;
      exit(1);
   }

   // We need to transform the pair (zone, face) to (zone, half-face).
   // The half-face ID uniquely labels a pair (zone,face) (so each interior
   // mesh face has two half-faces associated with it)
   int sface_offset = 0;
   for (int j = 0; j < nsfaces; ++j)
   {
      int bcid = shared_faces[sface_offset];
      int zone = shared_faces[sface_offset + 1];
      int face = shared_faces[sface_offset + 2];
      int corner1 = shared_faces[sface_offset + 3];
      int newface = face_to_newface_mapping[face];
      shared_faces_reord.push_back(bcid);
      shared_faces_reord.push_back(zone);
      shared_faces_reord.push_back(newface);
      shared_faces_reord.push_back(corner1);
      int corner2;
      if (ndim == 2)
      {
         corner2 = shared_faces[sface_offset + 4];
         shared_faces_reord.push_back(corner2);
         sface_offset += 5;
      }
      else
      {
         sface_offset += 4;
      }
      shared_zonefacepairs_reord.push_back(zone);
      shared_zonefacepairs_reord.push_back(newface);
   }
}

//////////// Added for blueprint format standardization ////////////
//
//           FOR TRANSFORMING BLUEPRINT MESH DESCRIPTION TO TETON'S
//           INTERNAL MESH DESCRIPTION. EVENTUALLY WILL NOT NEED TO
//           DO THIS
//
//////////// Added for blueprint format standardization ////////////

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
   int offset = 0;
   for (int zone = 0; zone < nzones; ++zone)
   {
      int ncorners_local = zone_to_corners2[zone].size();
      for (int c = 0; c < ncorners_local; ++c)
      {
         int corner = zone_to_corners2[zone][c];
         corner_to_lcorner[corner] = c;
      }
   }
}

void TetonBlueprint::ComputeZoneFaceToHalfFaceDict()
{
   int nzones = zone_to_faces2.size();
   int offset = 0;

   for (int zone = 0; zone < nzones; ++zone)
   {
      int nfaces_local = zone_to_faces2[zone].size();
      int loffset = 0;
      for (int f = 0; f < nfaces_local; ++f)
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

void TetonBlueprint::ProcessTopologyConnectivity(int nzones, int ncorners, conduit::Node &mesh_faces)
{
   int nvals = mesh_faces["face_to_zones"].dtype().number_of_elements();
   int offset = 0;
   while (offset < nvals)
   {
      int nzones_for_face = face_to_zones[offset];
      offset += 1;
      std::vector<int> zones;
      zones.push_back(face_to_zones[offset]);
      offset += 1;
      if (nzones_for_face == 2)
      {
         zones.push_back(face_to_zones[offset]);
         offset += 1;
      }
      face_to_zones2.push_back(zones);
   }

   offset = 0;
   zone_to_faces2.resize(nzones);
   for (int zone = 0; zone < nzones; ++zone)
   {
      int nfaces = zone_to_faces[offset];
      offset += 1;
      for (int k = 0; k < nfaces; ++k)
      {
         int face = zone_to_faces[offset];
         zone_to_faces2[zone].push_back(face);
         offset += 1;
      }
   }

   offset = 0;
   zone_to_corners2.resize(nzones);
   for (int zone = 0; zone < nzones; ++zone)
   {
      int ncorners = zone_to_corners[offset];
      offset += 1;

      for (int k = 0; k < ncorners; ++k)
      {
         int corner = zone_to_corners[offset];
         zone_to_corners2[zone].push_back(corner);
         offset += 1;
      }
   }

   nvals = mesh_faces["face_to_corners"].dtype().number_of_elements();
   offset = 0;
   while (offset < nvals)
   {
      int ncorners_on_face = face_to_corners[offset];
      offset += 1;
      std::vector<int> corners;
      for (int k = 0; k < ncorners_on_face; ++k)
      {
         int corner = face_to_corners[offset];
         corners.push_back(corner);
         offset += 1;
      }
      face_to_corners2.push_back(corners);
   }

   ComputeLocalCornersInZone(nzones, ncorners);
   ComputeCornerOffsets(nzones);
   ComputeZoneFaceToHalfFaceDict();
   //ComputeZoneFacesToFaces();
   // tag shared faces
   ComputeSharedFaces();
}

int TetonBlueprint::GetOppositeZone(int zone, int face)
{
   int nzones_on_face = face_to_zones2[face].size();
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

int TetonBlueprint::GetOppositeCorner(int zone, int face, int corner)
{
   int corner1 = -1;
   int ncorners = face_to_corners2[face].size();

   double node_x = corner_to_node_x[corner];
   double node_y = corner_to_node_y[corner];
   double node_z;
   int dim = mParametersNode["size/ndim"].as_int32();
   if (dim > 2)
      node_z = corner_to_node_z[corner];

   for (int c = 0; c < ncorners; ++c)
   {
      corner1 = face_to_corners2[face][c];
      //int node1 = corner_to_node[corner1];
      double node1_x = corner_to_node_x[corner1];
      double node1_y = corner_to_node_y[corner1];
      double node1_z;
      if (dim > 2)
         node1_z = corner_to_node_z[corner1];
      int zone1 = corner_to_zone[corner1];

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
      int halfface = zoneface_to_halfface[zoneface];
      return -(nbelem + halfface + 1);
   }
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
                                               std::vector<int> &zone_faces_to_bcids)
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
      int ncorner_faces = face_to_corners2[face].size();
      zone_faces_to_bcids[f] = face_to_bcid[face];

      for (int c = 0; c < ncorner_faces; ++c)
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

         int ocorner = GetOppositeCorner(zone, face, corner);
         // NOTE: this goes from 0-based indexing to 1-based indexing for Teton
         if (ocorner > -1)
            ocorner += 1;
         corners_opp.push_back(ocorner);
         ncorner_faces_total += 1;
         ncorners_on_zone_face[f] += 1;
      }
   }
}

void TetonBlueprint::ComputeSharedFaces()
{
   conduit::int32 *shared_faces_ptr;
   conduit::int32 *shared_zoneface_pairs_ptr;
   int nfaces = face_to_zones2.size();
   int ndim = mParametersNode["size/ndim"].as_int32();
   face_is_shared_bndry.resize(nfaces);
   for (int f = 0; f < nfaces; ++f)
   {
      face_is_shared_bndry[f] = false;
   }

   int sface_offset = 0;
   conduit::int32 *shared_faces_array_ptr;
   conduit::int32 *shared_zonefacepairs_array_ptr;
   int nsfaces;
   int sfaces_array_len = mMeshNode["shared_boundaries/shared_faces"].dtype().number_of_elements();
   if (ndim == 2)
   { // for each shared face, {zone, face, bcid, corner1, corner2} is stored
      // in the array shared_boundaries/shared_faces
      nsfaces = sfaces_array_len / 5;
   }
   else if (ndim == 3)
   { // for each shared face, {zone, face, bcid, corner1} is stored
      // in the array shared_boundaries/shared_faces
      nsfaces = sfaces_array_len / 4;
   }
   else
   {
      // TODO: CHANGE !!!!
      std::cerr << "the 1D case is not yet implemented!" << std::endl;
   }
   if (nsfaces > 0)
   {
      shared_faces_array_ptr = mMeshNode["shared_boundaries/shared_faces"].as_int32_ptr();
      shared_zonefacepairs_array_ptr = mMeshNode["shared_boundaries/shared_zone_and_face_pairs"].as_int32_ptr();
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
      int halfface = zoneface_to_halfface[zoneface];
      shared_faces_array_ptr[sface_offset + 2] = halfface;
      shared_zonefacepairs_array_ptr[2 * j + 1] = halfface;
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

   mMeshNode["shared_boundaries/shared_zone_and_face_pairs"].set(shared_zonefacepairs_array_ptr, 2 * nsfaces);
   mMeshNode["shared_boundaries/shared_faces"].set(shared_faces_array_ptr, sfaces_array_len);
}

void TetonBlueprint::ComputeTetonMeshConnectivityArray()
{
   int nzones = mParametersNode["size/nzones"].as_int32();
   int ncornr = mParametersNode["size/ncornr"].as_int32();
   int ndim = mParametersNode["size/ndim"].as_int32();

   zone_to_corners = mMeshNode["topologies"]["mesh"]["elements"]["zone_to_corners"].as_int32_ptr();
   corner_to_node_x = mMeshNode["topologies"]["mesh"]["corners"]["corner_to_node_x"].as_float64_ptr();
   corner_to_node_y = mMeshNode["topologies"]["mesh"]["corners"]["corner_to_node_y"].as_float64_ptr();
   if (ndim > 2)
   {
      corner_to_node_z = mMeshNode["topologies"]["mesh"]["corners"]["corner_to_node_z"].as_float64_ptr();
   }
   corner_to_zone = mMeshNode["topologies"]["mesh"]["corners"]["corner_to_zone"].as_int32_ptr();
   conduit::Node &mesh_faces = mMeshNode["topologies/mesh/faces"];
   conduit::Node &mesh_elements = mMeshNode["topologies/mesh/elements"];
   face_to_corners = mesh_faces["face_to_corners"].as_int32_ptr();
   face_to_zones = mesh_faces["face_to_zones"].as_int32_ptr();
   zone_to_faces = mesh_elements["zone_to_faces"].as_int32_ptr();
   face_to_bcid = mesh_faces["face_to_bcid"].as_int32_ptr();
   int nfaces = mesh_faces["face_to_corners"].dtype().number_of_elements();

   // process mesh topologies in to a more convenient format
   ProcessTopologyConnectivity(nzones, ncornr, mesh_faces);

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
   int connect_off_set = 0;
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
                                zone_faces_to_bcids);

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

      // Used as arguments to:
      //teton_setzone(&zoneID, &corner_offset, &nzone_faces,
      //              &ncorner_faces_total, &nzone_corners,
      //              &zones_opp[0], &corners_local[0], &corners_opp[0],
      //              &ncorners_on_zone_face[0], &zone_faces_to_bcids[0]);
   }

   mMeshNode["fields/corner_connectivity"].set(teton_connectivity.data(), teton_connectivity.size());
}

void TetonBlueprint::CreateTetonMeshCornerCoords()
{
   int nzones = mParametersNode["size/nzones"].as_int32();
   int ndim = mParametersNode["size/ndim"].as_int32();
   std::vector<double> zone_verts;
   std::vector<int> zone_ncorners;

   // loop over mesh elements
   for (int zone = 0; zone < nzones; ++zone)
   {
      int zoneID = zone + 1;
      int ncorners = zone_to_corners2[zone].size();
      zone_ncorners.push_back(ncorners);
      // loop over corners in the element
      for (int c = 0; c < ncorners; ++c)
      {
         int corner = zone_to_corners2[zone][c];
         zone_verts.push_back(corner_to_node_x[corner]);
         zone_verts.push_back(corner_to_node_y[corner]);
         if (ndim == 3)
         {
            zone_verts.push_back(corner_to_node_z[corner]);
         }
      }
      // used in this function in Teton
      //teton_setnodeposition(&zoneID, &zoneCoordinates[0]);
   }
   mMeshNode["fields/zone_verts"].set(zone_verts.data(), zone_verts.size());
   mMeshNode["fields/ncorners"].set(zone_ncorners.data(), zone_ncorners.size());
}

void TetonBlueprint::CheckSharedFaceData(mfem::ParMesh *pmesh,
                                         std::vector<double> &corner_to_node_x,
                                         std::vector<double> &corner_to_node_y,
                                         std::vector<double> &corner_to_node_z,
                                         std::vector<std::vector<int>> &shared_zone_and_face_pairs_all_groups,
                                         std::vector<std::vector<int>> &shared_faces_all_groups)
{
   int my_rank = pmesh->GetMyRank();
   std::ofstream outfile;
   outfile.open("shared_faces" + std::to_string(my_rank) + ".txt");

   int dim = pmesh->Dimension();
   int num_face_nbrs = pmesh->GetNFaceNeighbors();
   for (int fn = 0; fn < num_face_nbrs; fn++)
   {
      int nbr_group = pmesh->face_nbr_group[fn];
      int nbr_rank = pmesh->GetFaceNbrRank(fn);
      int offset = 0;
      outfile << "group:   " << nbr_group << std::endl;
      outfile << "nbr rank:   " << nbr_rank << std::endl;
      int corner_face_counter = 0;
      while (offset < shared_faces_all_groups[fn].size())
      {
         int bcID = shared_faces_all_groups[fn][offset];
         int elem = shared_faces_all_groups[fn][offset + 1];
         int face = shared_faces_all_groups[fn][offset + 2];
         int corner = shared_faces_all_groups[fn][offset + 3];
         int corner2;
         if (dim == 2)
            corner2 = shared_faces_all_groups[fn][offset + 4];
         double px = corner_to_node_x[corner];
         double py = corner_to_node_y[corner];
         double pz;
         if (dim == 3)
            pz = corner_to_node_z[corner];
         // get corner coordinates and face coordinates
         mfem::Array<int> vert;
         pmesh->GetFaceVertices(face, vert);
         int nvert = vert.Size();
         if (dim == 2)
         {
            offset += 5; // store pairs of corners for each shared face
         }
         else
         {
            offset += 4;
         }
         if (dim == 2)
         {
            px = corner_to_node_x[corner];
            py = corner_to_node_y[corner];
            outfile << "    corner-face count:   " << corner_face_counter << std::endl;
            corner_face_counter++;
            //for (int ivert = 0; ivert < nvert; ++ivert)
            //{
            //   const double *v = pmesh->GetVertex(vert[ivert]);
            //   //if (dim == 3) { shared_face_corners[ivert][2] = v[2]; }
            //   outfile << "        face vert coords:   "
            //           << "(" << v[1] << " , " << v[0] << ")" << std::endl;
            //}
            outfile << "        corner coords:   "
                    << "(" << py << " , " << px << ")" << std::endl;
            px = corner_to_node_x[corner2];
            py = corner_to_node_y[corner2];
            outfile << "    corner-face count:   " << corner_face_counter << std::endl;
            corner_face_counter++;
            //for (int ivert = 0; ivert < nvert; ++ivert)
            //{
            //   const double *v = pmesh->GetVertex(vert[ivert]);
            //   //if (dim == 3) { shared_face_corners[ivert][2] = v[2]; }
            //   outfile << "        face vert coords:   "
            //           << "(" << v[1] << " , " << v[0] << ")" << std::endl;
            //}
            //outfile << "        corner coords:   "
            //        << "(" << py << " , " << px << ")" << std::endl;
         }
         if (dim == 3)
         {
            px = corner_to_node_x[corner];
            py = corner_to_node_y[corner];
            pz = corner_to_node_z[corner];
            outfile << "    corner-face count:   " << corner_face_counter << std::endl;
            corner_face_counter++;
            //for (int ivert = 0; ivert < nvert; ++ivert)
            //{
            //   const double *v = pmesh->GetVertex(vert[ivert]);
            //   //if (dim == 3) { shared_face_corners[ivert][2] = v[2]; }
            //   outfile << "        face vert coords:   "
            //           << "(" << v[1] << " , " << v[0] << ")" << std::endl;
            //}
            outfile << "        corner coords:   "
                    << "(" << px << " , " << py << " , " << pz << ")" << std::endl;
         }
      }
   }
   outfile.close();
}

void TetonBlueprint::UpdateConduitMeshCoords(mfem::ParMesh *pmesh,
                                             int *corner_to_vertex,
                                             conduit::Node &mesh_conduit_node)

{
   mfem::Array<double> zone_verts;
   // loop over mesh elements
   int nzones = pmesh->GetNE();
   int ndim = pmesh->Dimension();
   int corner = 0;
   int zone_verts_size = 0;
   for (int zone = 0; zone < nzones; ++zone)
   {
      int ncorners = zone_to_corners2[zone].size();
      zone_verts_size += ncorners * ndim;
   }
   zone_verts.Reserve(zone_verts_size);
   for (int zone = 0; zone < nzones; ++zone)
   {
      int zoneID = zone + 1;
      int ncorners = zone_to_corners2[zone].size();
      // loop over corners in the element
      for (int c = 0; c < ncorners; ++c)
      {
         //int corner = zone_to_corners2[zone][c];
         // get the MFEM vertex ID corresponding to the corner
         int v = corner_to_vertex[corner];
         //double *vert_coords = pmesh->GetVertex(v);
         double vert_coords[3];
         pmesh->GetNode(v, vert_coords);
         // store the vertex coordinates
         if (ndim == 2)
         {
            // NOTE: need to swap the x-y ordering in rz geometry
            zone_verts.Append(vert_coords[1]);
            zone_verts.Append(vert_coords[0]);
         }
         else if (ndim == 3)
         {
            zone_verts.Append(vert_coords[0]);
            zone_verts.Append(vert_coords[1]);
            zone_verts.Append(vert_coords[2]);
         }
         else
         {
            std::cerr << "1D not implemented yet! " << std::endl;
            exit(1);
         }
         corner += 1;
      }
      // used in this function in Teton
      //teton_setnodeposition(&zoneID, &zoneCoordinates[0]);
   }
   mesh_conduit_node["fields/zone_verts"].set(zone_verts.GetData(), zone_verts.Size());
}

#endif
