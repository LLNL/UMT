#include "mfem.hpp"
#include <fstream>
#include <iostream>

mfem::Mesh makeUnstructBoxMesh(int zoneSplits)
{
   int red = 1;

   int top = 1;
   int bottom = 2;
   int left = 3;
   int right = 4;
   int front = 5;
   int back = 6;

   /*
     This makes a mesh that looks like this.
     On the outside, all nodes are spaced evenly.
     The g node is meant to make box zone 5 a square.
     Point f is somewhere.



   i      j k l
   *------*-*-*
   |      |3|5|
   |  2   | *-* h
   |      |/ g|
   *------* 4 |
   |e    f|\  |
   |      | \ |
   |      |  \|
   |  0   | 1 *d
   |      |   |
   *------*---*
   a      b   c
   */

   int dim = 3;
   int num2DVert = 12;
   int numVert = num2DVert * 3;
   int numElem = 6 * 2;
   int numBdrElem = 6 * 2;
   mfem::Mesh mesh(dim, numVert, numElem, numBdrElem, dim);

   double width = 1.0;

   // TODO: Fix to be "width/2" so that this is a cube.  But
   // maybe it doesn't matter.
   for (int i = 0; i < 3; ++i)
   {
      // a
      mesh.AddVertex(0.0, 0.0, i * width / 3.0);
      // b
      mesh.AddVertex(width / 2, 0.0, i * width / 3.0);
      // c
      mesh.AddVertex(width, 0.0, i * width / 3.0);
      // d
      mesh.AddVertex(width, width / 3.0, i * width / 3.0);
      // e
      mesh.AddVertex(0.0, width / 2.0, i * width / 3.0);
      // f
      mesh.AddVertex(width / 3.0, width / 2.0, i * width / 3.0);
      // g
      mesh.AddVertex(2.0 * width / 3.0, 2.0 * width / 3.0, i * width / 3.0);
      // h
      mesh.AddVertex(width, 2.0 * width / 3.0, i * width / 3.0);
      // i
      mesh.AddVertex(0.0, width, i * width / 3.0);
      // j
      mesh.AddVertex(width / 3.0, width, i * width / 3.0);
      // k
      mesh.AddVertex(2.0 * width / 3.0, width, i * width / 3.0);
      // l
      mesh.AddVertex(width, width, i * width / 3.0);
   }

   int a = 0;
   int b = 1;
   int c = 2;
   int d = 3;
   int e = 4;
   int f = 5;
   int g = 6;
   int h = 7;
   int i = 8;
   int j = 9;
   int k = 10;
   int l = 11;

   // MFEM for hexes and quads (and a few others) shares the VTK file format.  See
   // https://kitware.github.io/vtk-examples/site/VTKFileFormats/
   int B, T;
   B = 0 * num2DVert; // bottom layer
   T = 1 * num2DVert; // top layer
   mesh.AddHex(a + B, b + B, f + B, e + B, a + T, b + T, f + T, e + T, red);
   mesh.AddHex(b + B, c + B, d + B, f + B, b + T, c + T, d + T, f + T, red);
   mesh.AddHex(e + B, f + B, j + B, i + B, e + T, f + T, j + T, i + T, red);
   mesh.AddHex(f + B, g + B, k + B, j + B, f + T, g + T, k + T, j + T, red);
   mesh.AddHex(f + B, d + B, h + B, g + B, f + T, d + T, h + T, g + T, red);
   mesh.AddHex(g + B, h + B, l + B, k + B, g + T, h + T, l + T, k + T, red);
   B = 1 * num2DVert;
   T = 2 * num2DVert;
   mesh.AddHex(a + B, b + B, f + B, e + B, a + T, b + T, f + T, e + T, red);
   mesh.AddHex(b + B, c + B, d + B, f + B, b + T, c + T, d + T, f + T, red);
   mesh.AddHex(e + B, f + B, j + B, i + B, e + T, f + T, j + T, i + T, red);
   mesh.AddHex(f + B, g + B, k + B, j + B, f + T, g + T, k + T, j + T, red);
   mesh.AddHex(f + B, d + B, h + B, g + B, f + T, d + T, h + T, g + T, red);
   mesh.AddHex(g + B, h + B, l + B, k + B, g + T, h + T, l + T, k + T, red);

   // Should be counter-clockwise from the outside of the mesh.
   B = 0 * num2DVert;
   mesh.AddBdrQuad(a + B, e + B, f + B, b + B, bottom);
   mesh.AddBdrQuad(b + B, f + B, d + B, c + B, bottom);
   mesh.AddBdrQuad(e + B, i + B, j + B, f + B, bottom);
   mesh.AddBdrQuad(f + B, j + B, k + B, g + B, bottom);
   mesh.AddBdrQuad(f + B, g + B, h + B, d + B, bottom);
   mesh.AddBdrQuad(g + B, k + B, l + B, h + B, bottom);

   B = 2 * num2DVert;
   mesh.AddBdrQuad(a + B, b + B, f + B, e + B, top);
   mesh.AddBdrQuad(b + B, c + B, d + B, f + B, top);
   mesh.AddBdrQuad(e + B, f + B, j + B, i + B, top);
   mesh.AddBdrQuad(f + B, g + B, k + B, j + B, top);
   mesh.AddBdrQuad(f + B, d + B, h + B, g + B, top);
   mesh.AddBdrQuad(g + B, h + B, l + B, k + B, top);

   T = 1 * num2DVert;
   B = 0 * num2DVert;
   mesh.AddBdrQuad(a + B, a + T, e + T, e + B, left);
   mesh.AddBdrQuad(e + B, e + T, i + T, i + B, left);
   T = 2 * num2DVert;
   B = 1 * num2DVert;
   mesh.AddBdrQuad(a + B, a + T, e + T, e + B, left);
   mesh.AddBdrQuad(e + B, e + T, i + T, i + B, left);

   T = 1 * num2DVert;
   B = 0 * num2DVert;
   mesh.AddBdrQuad(b + B, b + T, a + T, a + B, front);
   mesh.AddBdrQuad(c + B, c + T, b + T, b + B, front);
   T = 2 * num2DVert;
   B = 1 * num2DVert;
   mesh.AddBdrQuad(b + B, b + T, a + T, a + B, front);
   mesh.AddBdrQuad(c + B, c + T, b + T, b + B, front);

   T = 1 * num2DVert;
   B = 0 * num2DVert;
   mesh.AddBdrQuad(d + B, d + T, c + T, c + B, right);
   mesh.AddBdrQuad(h + B, h + T, d + T, d + B, right);
   mesh.AddBdrQuad(l + B, l + T, h + T, h + B, right);
   T = 2 * num2DVert;
   B = 1 * num2DVert;
   mesh.AddBdrQuad(d + B, d + T, c + T, c + B, right);
   mesh.AddBdrQuad(h + B, h + T, d + T, d + B, right);
   mesh.AddBdrQuad(l + B, l + T, h + T, h + B, right);

   T = 1 * num2DVert;
   B = 0 * num2DVert;
   mesh.AddBdrQuad(k + B, k + T, l + T, l + B, right);
   mesh.AddBdrQuad(j + B, j + T, k + T, k + B, right);
   mesh.AddBdrQuad(i + B, i + T, j + T, j + B, right);
   T = 2 * num2DVert;
   B = 1 * num2DVert;
   mesh.AddBdrQuad(k + B, k + T, l + T, l + B, right);
   mesh.AddBdrQuad(j + B, j + T, k + T, k + B, right);
   mesh.AddBdrQuad(i + B, i + T, j + T, j + B, right);

   // Build internal structures
   mesh.FinalizeTopology();

   if (zoneSplits > 0)
   {
      // Refine the initial mesh with this many zones in each
      // dimension in each element above.
      int ref_type = mfem::BasisType::ClosedUniform;
      mesh = mfem::Mesh::MakeRefined(mesh, zoneSplits, ref_type);
   }

   if (int wrong = mesh.CheckElementOrientation(true) > 0)
   {
      std::cout << "There were " << wrong << " 3D mesh elements with the wrong orientation.\n";
   }
   if (int wrong = mesh.CheckBdrElementOrientation(true) > 0)
   {
      std::cout << "There were " << wrong << " 3D mesh boundary elements with the wrong orientation.\n";
   }

   // Sort the grid for better locality
   mfem::Array<int> ordering;
   mesh.GetHilbertElementOrdering(ordering);
   mesh.ReorderElements(ordering);

   return mesh;
}

void twist(const mfem::Vector &in, mfem::Vector &p)
{
   const double x = in[0] - 0.5;
   const double y = in[1] - 0.5;
   const double z = in[2] - 0.5;
   const double r = std::hypot(x, y);
   const double theta = std::atan2(y, x);
   p[0] = r * std::cos(theta + 0.2 * M_PI * z) + 0.5;
   p[1] = r * std::sin(theta + 0.2 * M_PI * z) + 0.5;
   p[2] = z + 0.5;
}

int main(int argc, char *argv[])
{
   MPI_Comm comm = MPI_COMM_WORLD;
   int request = MPI_THREAD_SINGLE;
   int provided = 0;
   int claimed = 0;
   MPI_Init_thread(&argc, &argv, request, &provided);
   int myRank = 0;
   int mySize = 0;

   MPI_Comm_rank(comm, &myRank);
   MPI_Comm_size(comm, &mySize);

   //int zoneSplit = 10;
   int zoneSplit = 0;

   mfem::Mesh mesh3D = makeUnstructBoxMesh(zoneSplit);

   // Make it harder.  Only do once on the fully refined mesh.
   //mesh3D.Transform(twist);

   if (myRank == 0)
   {
      if (int wrong = mesh3D.CheckElementOrientation(true) > 0)
      {
         std::cout << "There were " << wrong << " 3D mesh elements with the wrong orientation after reordering.\n";
      }
      if (int wrong = mesh3D.CheckBdrElementOrientation(true) > 0)
      {
         std::cout << "There were " << wrong
                   << " 3D mesh boundary elements with the wrong orientation after reordering.\n";
      }

      mfem::VisItDataCollection vdc("unstructBox3D", &mesh3D);
      vdc.Save();

      // Save mesh to file
      std::ofstream mesh_ofs("unstructBox3D.mesh");
      mesh_ofs.precision(16);
      mesh3D.Print(mesh_ofs);
      mesh_ofs.close();

      // Show info about the mesh.
      mesh3D.PrintInfo();
   }

   /*
   mesh3D.SetCurvature(1);
   mfem::ParMesh pmesh(comm, mesh3D);

   int ref_type = mfem::BasisType::ClosedUniform;
   pmesh = mfem::ParMesh::MakeRefined(pmesh, 2, ref_type);
   pmesh.PrintInfo();

   if (int wrong = pmesh.CheckElementOrientation(true) > 0)
   {
      std::cout << "There were " << wrong << " 3D mesh elements with the wrong orientation.\n";
   }
   if (int wrong = pmesh.CheckBdrElementOrientation(true) > 0)
   {
      std::cout << "There were " << wrong << " 3D mesh boundary elements with the wrong orientation.\n";
   }
*/
   return 0;
}
