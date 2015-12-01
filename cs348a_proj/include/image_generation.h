#ifndef IMAGE_GENERATION_H
#define IMAGE_GENERATION_H

#include "mesh_definitions.h"
#include <string>

void writeImage(Mesh &mesh, int width, int height, const std::string& filename, 
                const OpenMesh::Vec3f& camPos, const OpenMesh::Vec3f& camLookDir,
                float nDotViewMax, float DwkrMin,
                OpenMesh::EPropHandleT<bool>& edgeVisited, OpenMesh::FPropHandleT<bool>& faceVisited,
                OpenMesh::VPropHandleT<double>& viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f>& viewCurvatureDerivative);

#endif
