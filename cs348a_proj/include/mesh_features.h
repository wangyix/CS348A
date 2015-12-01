#ifndef MESH_FEATURES_H
#define MESH_FEATURES_H

#include "mesh_definitions.h"

bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, OpenMesh::Vec3f cameraPos);
bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e);
bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, OpenMesh::Vec3f cameraPos);

bool isContourLineAcceptable(Mesh& mesh, Mesh::FaceHandle fh,
                             const OpenMesh::Vec3f& s, const OpenMesh::Vec3f& t, const OpenMesh::Vec3f& camPos,
                             float nDotViewMax, float DwkrMin,
                             OpenMesh::FPropHandleT<OpenMesh::Vec3f>& viewCurvatureDerivative);

bool isSuggestiveContourFace(Mesh& mesh, Mesh::FaceHandle fh, const OpenMesh::Vec3f& actualCamPos,
                             OpenMesh::VPropHandleT<double>& viewCurvature,
                             OpenMesh::FPropHandleT<OpenMesh::Vec3f>& viewCurvatureDerivative,
                             float nDotViewMax, float DwkrMin,
                             OpenMesh::Vec3f* s, OpenMesh::Vec3f* t,
                             Mesh::HalfedgeHandle* s_edge, Mesh::HalfedgeHandle* t_edge);

#endif
