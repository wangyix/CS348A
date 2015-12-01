#include <Eigen/Dense>
#include "mesh_features.h"

using namespace OpenMesh;
using namespace Eigen;

bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos)  {
  // CHECK IF e IS A SILHOUETTE HERE -----------------------------------------------------------------------------
  Mesh::HalfedgeHandle heh = mesh.halfedge_handle(e, 0);

  Vec3f v0 = mesh.point(mesh.from_vertex_handle(heh));
  Vec3f v1 = mesh.point(mesh.to_vertex_handle(heh));
  Vec3f midpoint = 0.5f * (v0 + v1);
  Vec3f toCamera = cameraPos - midpoint;

  Vec3f n0 = mesh.calc_face_normal(mesh.face_handle(heh));
  Vec3f n1 = mesh.calc_face_normal(mesh.opposite_face_handle(heh));

  bool f0FacesCamera = (n0 | toCamera) > 0.0;
  bool f1FacesCamera = (n1 | toCamera) > 0.0;
  return f0FacesCamera != f1FacesCamera;
  // -------------------------------------------------------------------------------------------------------------
}

bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e) {
  // CHECK IF e IS SHARP HERE ------------------------------------------------------------------------------------
  Mesh::HalfedgeHandle heh = mesh.halfedge_handle(e, 0);
  Vec3f n0 = mesh.calc_face_normal(mesh.face_handle(heh));
  Vec3f n1 = mesh.calc_face_normal(mesh.opposite_face_handle(heh));
  return (n0 | n1) < 0.5f;
  // -------------------------------------------------------------------------------------------------------------
}

bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos) {
  return mesh.is_boundary(e) || isSilhouette(mesh, e, cameraPos) || isSharpEdge(mesh, e);
}

bool isContourLineAcceptable(Mesh& mesh, Mesh::FaceHandle fh,
                             const Vec3f& s, const Vec3f& t, const Vec3f& camPos,
                             float nDotViewMax, float DwkrMin,
                             FPropHandleT<Vec3f>& viewCurvatureDerivative) {
  Vec3f toCamera = camPos - 0.5f*(s + t);   // v
  Vec3f n = mesh.calc_face_normal(fh);
  if ((n | toCamera.normalized()) > nDotViewMax) {
    // angle between face normal and is small; bail.
    return false;
  }
  Vec3f w = toCamera - (toCamera | n)*n;    // w is v projected onto tangent plane
  Vec3f kwGradient = mesh.property(viewCurvatureDerivative, fh);
  return ((kwGradient | w.normalized()) > DwkrMin);   // reject if Dwkr is negative or too small of a positive; bail.
}


bool isSuggestiveContourFace(Mesh& mesh, Mesh::FaceHandle fh, const Vec3f& actualCamPos, 
                             VPropHandleT<double>& viewCurvature,
                             FPropHandleT<Vec3f>& viewCurvatureDerivative,
                             float nDotViewMax, float DwkrMin,
                             Vec3f* s, Vec3f* t,
                             Mesh::HalfedgeHandle* s_edge, Mesh::HalfedgeHandle* t_edge) {
  Mesh::HalfedgeHandle heh[3];
  double kw[3];
  int numKwPositive = 0, posKwIndex = -1, negKwIndex = -1;
  {
    int i = 0;
    Mesh::FaceHalfedgeCCWIter fh_it, fh_it_end = mesh.fh_ccwend(fh);
    for (fh_it = mesh.fh_ccwbegin(fh); fh_it != fh_it_end; ++fh_it) {
      heh[i] = *fh_it;
      Mesh::VertexHandle to_vh = mesh.to_vertex_handle(*fh_it);
      kw[i] = mesh.property(viewCurvature, to_vh);
      if (kw[i] > 0.0) {
        numKwPositive++;
        posKwIndex = i;
      } else {
        negKwIndex = i;
      }
      i++;
    }
    assert(i == 3);
  }
  if (numKwPositive == 0 || numKwPositive == 3) {
    // face does not have kw zero-crossing; bail.
    return false;
  }

  // assign v0 to be the odd-man-out vertex, and v1,v2 to be the other two vertices
  // i0, i1, i2 are in CCW order
  int i0 = (numKwPositive == 1) ? posKwIndex : negKwIndex;
  assert(0 <= i0 && i0 < 3);
  int i1 = (i0 + 1) % 3;
  int i2 = (i0 + 2) % 3;
  Vec3f v0 = mesh.point(mesh.to_vertex_handle(heh[i0]));
  Vec3f v1 = mesh.point(mesh.to_vertex_handle(heh[i1]));
  Vec3f v2 = mesh.point(mesh.to_vertex_handle(heh[i2]));
  double kw0 = kw[i0];
  double kw1 = kw[i1];
  double kw2 = kw[i2];
  // compute the two edge points where kw = 0
  double a1 = kw0 / (kw0 - kw1);
  Vec3f p1 = (1.0 - a1)*v0 + a1*v1;
  double a2 = kw0 / (kw0 - kw2);
  Vec3f p2 = (1.0 - a2)*v0 + a2*v2;

  if (!isContourLineAcceptable(mesh, fh, p1, p2, actualCamPos, nDotViewMax, DwkrMin,
                               viewCurvatureDerivative)) {
    return false;
  }

  *s = p2;
  *t = p1;
  *s_edge = heh[i0];
  *t_edge = heh[i1];
  return true;
}
