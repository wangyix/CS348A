#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "curvature.h"

using namespace OpenMesh;
using namespace Eigen;
using namespace std;

void computeCurvature(Mesh& mesh, OpenMesh::VPropHandleT<CurvatureInfo>& curvature) {
  // WRITE CODE HERE TO COMPUTE THE CURVATURE AT THE CURRENT VERTEX ----------------------------------------------
  // -------------------------------------------------------------------------------------------------------------
}

void computeViewCurvature(Mesh& mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo>& curvature,
  OpenMesh::VPropHandleT<double>& viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f>& viewCurvatureDerivative) {
  // WRITE CODE HERE TO COMPUTE CURVATURE IN THE VIEW PROJECTION PROJECTED ON THE TANGENT PLANE ------------------
  // Compute vector to viewer and project onto tangent plane, then use components in principal directions to find curvature
  // -------------------------------------------------------------------------------------------------------------

  // We'll use the finite elements piecewise hat method to find per-face gradients of the view curvature
  // CS 348a doesn't cover how to differentiate functions on a mesh (Take CS 468!) so we provide code here

  for (Mesh::ConstFaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
    const Mesh::FaceHandle fh = (*f_it);

    Vec3f p[3];
    double c[3];
    Mesh::ConstFaceVertexIter fv_it = mesh.cfv_begin(fh);
    for (int i = 0; i < 3 && fv_it != mesh.cfv_end(fh); ++i, ++fv_it) {
      const Mesh::VertexHandle n_vh = (*fv_it);
      p[i] = mesh.point(n_vh);
      c[i] = mesh.property(viewCurvature, n_vh);
    }

    const Vec3f n = mesh.normal(fh);
    double area = mesh.calc_sector_area(mesh.halfedge_handle(fh));

    mesh.property(viewCurvatureDerivative, fh) =
      cross(n, (p[0] - p[2])) * (c[1] - c[0]) / (2 * area) +
      cross(n, (p[1] - p[0])) * (c[2] - c[0]) / (2 * area);
  }
}
