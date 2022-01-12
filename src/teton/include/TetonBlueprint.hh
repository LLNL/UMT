// Blueprint
// This class provides functions for converting a MFEM mesh to a teton conduit mesh node.

#if defined(TETON_ENABLE_MFEM)

#ifndef TETONBLUEPRINT_HPP__
#define TETONBLUEPRINT_HPP__

#include "mfem.hpp"
#include "conduit/conduit.hpp"
#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_relay.hpp"
#include "TetonConduitInterface.hh"

class TetonBlueprint
{
private:
   conduit::Node &mParametersNode;
   conduit::Node &mMeshNode;

   int vert_prev_global, global_vert_counter;
   std::vector<std::vector<int>> vert_to_counters;
   // SWAP ORDERING OF corner_to_node_x AND corner_to_node_x
   // SINCE APPARENTLY THIS IS DONE FOR rz GEOMETRY
   bool swap_axis_ordering = true;
   std::map<int, int> boundary_id_to_type;
   std::vector<double> source_bndry_temperatures;
   // for determining vertex ordering
   int v1_last_glob, v2_last_glob, v_start_glob;
   std::map<std::set<int>, int> elem_edges_count;
   // REMOVE, FOR DEBUGGING ONLY //
   int nfaces_count_glob = 0;
   // REMOVE, FOR DEBUGGING ONLY //
   // TODO: make this more sane
   int number_of_source_profiles = 0;

   // REMOVE //
   int myrank_ = 0;
   // REMOVE //

   // For converting from the Blueprint Conduit Node format
   // to Teton's standalone Conduit Node format
   void ProcessTopologyConnectivity(int nzones, int ncorners, conduit::Node &mesh_faces);
   int GetOppositeZone(int zone, int face);
   int GetOppositeCorner(int zone, int face, int corner);
   void ComputeCornerOffsets(int nzones);
   void ComputeLocalCornersInZone(int nzones, int ncorners);
   void ComputeZoneAndFaceToHalfFaceArray();
   void ComputeZoneFaceToHalfFaceDict();
   void ComputeZoneFacesToFaces();
   void ComputeSharedFaces();
   void ComputeConnectivityArrays(int zone,
                                  int &corner_offset,
                                  int &nzone_faces,
                                  int &ncorner_faces,
                                  int &ncorners,
                                  std::vector<int> &zones_opp,
                                  std::vector<int> &corners_local,
                                  std::vector<int> &corners_opp,
                                  std::vector<int> &ncorners_on_zone_face,
                                  std::vector<int> &zone_faces_to_bcids);

   // For blueprint standardization
   int dim;
   int nbelem = 0;
   int n_shared_corner_faces = 0;
   conduit::int32 *zone_to_faces;
   conduit::int32 *zone_to_nodes;
   conduit::int32 *zone_to_corners;
   conduit::int32 *face_to_zones;
   conduit::int32 *face_to_corners;
   conduit::int32 *face_to_bcid;
   conduit::float64 *corner_to_node_x;
   conduit::float64 *corner_to_node_y;
   conduit::float64 *corner_to_node_z;
   conduit::int32 *corner_to_zone;
   int nhalffaces = 0;
   std::vector<std::vector<int>> zone_and_face_to_halfface;
   std::vector<bool> face_is_shared_bndry;
   std::map<std::pair<int, int>, int> zoneface_to_halfface;
   std::map<std::pair<int, int>, int> zoneface_to_halfface2;
   std::map<std::pair<int, int>, int> zoneface_to_lface;
   std::vector<int> zoneface_to_face;
   std::vector<std::vector<int>> zone_to_faces2;
   std::vector<std::vector<int>> zone_to_corners2;
   std::vector<std::vector<int>> zone_to_nodes2;
   std::vector<std::vector<int>> face_to_zones2;
   std::vector<std::vector<int>> face_to_corners2;
   std::vector<int> corner_to_lcorner;
   std::vector<int> corner_offsets;
   std::vector<int> unique_zone_ordering;

public:
   TetonBlueprint(Teton::Teton &teton) : mMeshNode(teton.getMeshBlueprint()), mParametersNode(teton.getOptions())
   {
   }
   ~TetonBlueprint()
   {
   }
   void SetTimeStep(double dt)
   {
      mParametersNode["iteration/dtrad"] = dt;
   }
   conduit::Node &GetConduitMeshNode()
   {
      return mMeshNode;
   }
   conduit::Node &GetConduitInputNode()
   {
      return mParametersNode;
   }
   void SetBoundaryIDs(std::map<int, int> &boundary_id_to_type0)
   {
      boundary_id_to_type = boundary_id_to_type0;
   }
   void SetSourceBoundaryTemperature(double Tval)
   {
      source_bndry_temperatures.push_back(Tval);
   }
   void ComputeFaceToBndryElement(mfem::Mesh *mesh, mfem::Array<int> &face_to_bndry);
   void ComputeLfaceToSface(mfem::ParMesh *mesh, mfem::Array<int> &lface_to_sface);
   void OutputConduitMesh(mfem::ParMesh *pmesh,
                          mfem::Vector &density,
                          mfem::Vector &heat_capacity,
                          mfem::Vector &rad_temp,
                          mfem::Vector &material_temp,
                          const mfem::Vector &gr_bounds,
                          mfem::Vector &abs_opacity,
                          mfem::Vector &scat_opacity,
                          mfem::Vector &electron_density);
   void SortSharedFaces(mfem::ParMesh *pmesh,
                        std::vector<std::vector<int>> &zone_and_vertex_to_corner,
                        std::vector<int> &face_to_bcid,
                        std::vector<int> &shared_faces,
                        std::vector<int> &shared_zone_and_face_pairs);
   void SortSharedFaces(mfem::ParMesh *pmesh,
                        std::map<std::pair<int, int>, int> &zone_and_gvertex_to_lvertex,
                        std::vector<std::vector<int>> &zone_and_lvertex_to_corner,
                        std::vector<int> &face_to_bcid,
                        int group,
                        std::vector<std::vector<int>> &gr_sfaces,
                        std::vector<int> &shared_faces,
                        std::vector<int> &shared_zone_and_face_pairs);
   void ReorderSharedFaces(int ndim,
                           std::vector<int> &face_to_newface_mapping,
                           std::vector<int> &shared_faces,
                           std::vector<int> &shared_faces_reord,
                           std::vector<int> &shared_zonefacepairs_reord);
   void ComputeFaceIDs(mfem::ParMesh *pmesh,
                       std::vector<int> &face_to_ncorners,
                       std::map<int, std::vector<int>> &boundaries,
                       int &nbelem_corner_faces,
                       std::vector<std::vector<int>> &gr_sface,
                       std::vector<int> &boundaries_types,
                       std::vector<int> &boundary_conditions,
                       std::vector<int> &face_to_bcid);
   void CreateSharedFacesTable(mfem::ParMesh *pmesh, mfem::Array<mfem::Connection> &table_list);
   void CreateSharedFacesTable(mfem::ParMesh *pmesh, std::vector<std::vector<int>> &group_sfaces);
   void CheckSharedFaceData(mfem::ParMesh *pmesh,
                            std::vector<double> &corner_to_node_x,
                            std::vector<double> &corner_to_node_y,
                            std::vector<double> &corner_to_node_z,
                            std::vector<std::vector<int>> &shared_zone_and_face_pairs_all_groups,
                            std::vector<std::vector<int>> &shared_faces_all_groups);
   void FormBoundaryIDToBCType(std::map<int, int> &boundary_id_to_type);
   void OutputConduitParameters(mfem::ParMesh *pmesh,
                                mfem::Vector &density,
                                mfem::Vector &heat_capacity,
                                mfem::Vector &rad_temp,
                                mfem::Vector &material_temp,
                                const mfem::Vector &gr_bounds,
                                mfem::Vector &abs_opacity,
                                mfem::Vector &scat_opacity,
                                mfem::Vector &electron_density);

   void ComputeTetonMeshConnectivityArray();
   void CreateTetonMeshCornerCoords();
   void OrderVertices(mfem::ParMesh *pmesh,
                      int elem,
                      std::vector<int> &vertices_ordered,
                      std::vector<int> &faces_ordered);
   void OrderVertices2D(mfem::ParMesh *pmesh, int elem, int elem_face, std::vector<int> &vertices_ordered);
   void OrderFaces3D(mfem::ParMesh *pmesh, int elem, std::vector<int> &faces_ordered);
   void OrderVertices3D(mfem::ParMesh *pmesh,
                        int elem,
                        std::vector<int> &vertices_ordered,
                        std::vector<int> &faces_ordered,
                        std::vector<std::vector<int>> &face_verts_ordered);
   bool CompareVectors(std::vector<double> &vec1, std::vector<double> &vec2);
   bool ComparePairs(const std::pair<double, double> &x1, const std::pair<double, double> &x2);
   void CreateZoneOrdering(mfem::Mesh *mesh);
   void UpdateMaterialProperties(mfem::Vector &heat_capacity, mfem::Vector &abs_opacity, mfem::Vector &scat_opacity);
   void UpdateConduitMeshCoords(mfem::ParMesh *pmesh, int *corner_to_vertex, conduit::Node &mesh_conduit_node);

   // REMOVE //
   // For new blueprint
   void ComputeConnectivityViaConduitLib(mfem::ParMesh *pmesh);
   // REMOVE //
};

#endif // TETONBLUEPRINT_HPP__
#endif // TETON_ENABLE_MFEM
