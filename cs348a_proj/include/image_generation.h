#ifndef IMAGE_GENERATION_H
#define IMAGE_GENERATION_H

#include "mesh_definitions.h"
#include <string>

//bool isVisible(OpenMesh::Vec3f point);
void writeImage(Mesh &mesh, int width, int height, std::string filename, 
                const OpenMesh::Vec3f& camPos, const OpenMesh::Vec3f& camLookDir,
                OpenMesh::EPropHandleT<bool>& edgeVisited);

#endif
