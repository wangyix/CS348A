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
