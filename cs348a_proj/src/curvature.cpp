#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "curvature.h"

using namespace OpenMesh;
using namespace Eigen;
using namespace std;

void computeCurvature(Mesh& mesh, OpenMesh::VPropHandleT<CurvatureInfo>& curvature) {
  // WRITE CODE HERE TO COMPUTE THE CURVATURE AT EACH VERTEX ----------------------------------------------
  Mesh::VertexIter v_it, v_it_end = mesh.vertices_end();
  for (v_it = mesh.vertices_begin(); v_it != v_it_end; ++v_it) {
    const Mesh::VertexHandle vh = (*v_it);
    const Vec3f normal = mesh.normal(vh);
    const Vector3d n(normal[0], normal[1], normal[2]);

    Vec3f vi_temp = mesh.point(vh);
    Vector3d vi(vi_temp[0], vi_temp[1], vi_temp[2]);

    Matrix3d Mvi = Matrix3d::Zero();
    double wSum = 0.0;
    for (Mesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(vh); voh_it; ++voh_it) {
      Mesh::VertexHandle vh_j = mesh.to_vertex_handle(*voh_it);
      Vec3f vj_temp = mesh.point(vh_j);
      Vector3d vj(vj_temp[0], vj_temp[1], vj_temp[2]);

      Vector3d Tij = ((Matrix3d::Identity() - n * n.transpose()) * (vi - vj)).normalized();
      double kij = (2.0 * n.dot(vj - vi)) / (vj - vi).squaredNorm();

      double leftFaceArea = mesh.calc_sector_area(*voh_it);
      double rightFaceArea = mesh.calc_sector_area(mesh.opposite_halfedge_handle(*voh_it));
      double wij = leftFaceArea + rightFaceArea;

      Mvi += (wij * kij * Tij * Tij.transpose());
      wSum += wij;
    }
    Mvi /= wSum;

    // Find T1,T2, which are the eigenvectors of Mvi with nonzero eigenvalues
    EigenSolver<Matrix3d> es(Mvi);
    Vector3d evalues = es.eigenvalues().real();
    Matrix3d evectors = es.eigenvectors().real();
    // find min eigenvalue; we'll assume that's the 0 eigenvalue
    double minAbsEvalue = std::numeric_limits<double>::max();
    int zeroEvalueAt = -1;
    for (int i = 0; i < 3; i++) {
      double absEvalue = abs(evalues[i]);
      if (absEvalue < minAbsEvalue) {
        minAbsEvalue = absEvalue;
        zeroEvalueAt = i;
      }
    }
    // extract eigenvectors T1,T2 and compute k1,k2 from eigenvalues
    int index1 = (zeroEvalueAt == 0) ? 1 : 0;
    int index2 = (zeroEvalueAt == 2) ? 1 : 2;
    double m11 = evalues[index1];
    double m22 = evalues[index2];
    Vector3d T1 = evectors.col(index1).normalized();
    Vector3d T2 = evectors.col(index2).normalized();
    double k1 = 3.0*m11 - m22;
    double k2 = 3.0*m22 - m11;

    // In the end you need to fill in this struct
    CurvatureInfo info;
    info.curvatures[0] = -k1;   // DeCarlo's kr defined to have opposite sign of Taubin's kw?
    info.curvatures[1] = -k2;
    info.directions[0] = Vec3f(T1[0], T1[1], T1[2]);
    info.directions[1] = Vec3f(T2[0], T2[1], T2[2]);

    mesh.property(curvature, vh) = info;
  }
  // -------------------------------------------------------------------------------------------------------------
}

void computeViewCurvature(Mesh& mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo>& curvature,
  OpenMesh::VPropHandleT<double>& viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f>& viewCurvatureDerivative) {
  // WRITE CODE HERE TO COMPUTE CURVATURE IN THE VIEW PROJECTION PROJECTED ON THE TANGENT PLANE ------------------
  // Compute vector to viewer and project onto tangent plane, then use components in principal directions to find curvature
  Mesh::VertexIter v_it, v_it_end = mesh.vertices_end();
  for (v_it = mesh.vertices_begin(); v_it != v_it_end; ++v_it) {
    Mesh::VertexHandle vh = *v_it;
    Vec3f toCamera = camPos - mesh.point(vh);
    Vec3f n = mesh.normal(vh);
    Vec3f w = (toCamera - (toCamera | n)*n).normalized();
    
    CurvatureInfo info = mesh.property(curvature, vh);
    float cosTheta = w | info.directions[0];
    float cosThetaSq = cosTheta * cosTheta;
    mesh.property(viewCurvature, vh) = cosThetaSq * info.curvatures[0] + (1.f - cosThetaSq) * info.curvatures[1];
  }
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
