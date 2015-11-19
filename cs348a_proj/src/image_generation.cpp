#include "image_generation.h"
#include "mesh_features.h"
#include <GL/glut.h>
#include <fstream>
//#include <set>
//#include <map>
//#include <list>
#include <deque>

using namespace OpenMesh;
using namespace std;

Vec3f toImagePlane(Vec3f point) {
  GLdouble point3DX = point[0], point3DY = point[1], point3DZ = point[2];

  GLdouble modelMatrix[16], projMatrix[16];
  GLint viewport[4];
  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glGetIntegerv(GL_VIEWPORT, viewport);

  GLdouble point2DX, point2DY, point2DZ;
  gluProject(point3DX, point3DY, point3DZ, modelMatrix, projMatrix, viewport, &point2DX, &point2DY, &point2DZ);

  return Vec3f(point2DX,point2DY,point2DZ);
}

/*// Adapted from
// http://stackoverflow.com/questions/1311869/opengl-how-to-determine-if-a-3d-rendered-point-is-occluded-by-other-3d-rende
bool isVisible(Vec3f point) {
  Vec3f projected = toImagePlane(point);

  GLfloat bufDepth = 0.0;
  glReadPixels(static_cast<GLint>( projected[0] ), static_cast<GLint>( projected[1] ),
      1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &bufDepth);

  GLdouble EPSILON = 0.01;
  return (bufDepth - projected[2]) > -EPSILON; // check sign!
}*/

// Returns true if edge is wholly or partially visible. If only partially visible, the occluded endpoint
// will be moved along the edge to the boundary between the visible and occluded portions of the edge.
bool isEdgePartialVisible(const Vec3f& s, const Vec3f& t, const GLfloat* depthBuffer, int width,
                          Vec3f* s_vis, Vec3f* t_vis) {
  GLdouble EPSILON = 0.0075;

  Vec3f s_proj = toImagePlane(s);
  Vec3f t_proj = toImagePlane(t);
  float s_depth = s_proj[2];
  float t_depth = t_proj[2];  
  bool sIsVisible = (depthBuffer[(GLint)s_proj[1]*width + (GLint)s_proj[0]] - s_depth) > -EPSILON;
  bool tIsVisible = (depthBuffer[(GLint)t_proj[1]*width + (GLint)t_proj[0]] - t_depth) > -EPSILON;
  
  Vec3f a, b;   // a will be visible bound, b will be occluded bound
  Vec3f a_proj, b_proj;
  Vec3f** vis;  // ptr to the endpoint that will be set to the result of the binary search
  if (sIsVisible) {
    if (tIsVisible) {
      *s_vis = s;
      *t_vis = t;
      return true;
    } else {
      *s_vis = s;
      vis = &t_vis;
      a = s;
      a_proj = s_proj;
      b = t;
      b_proj = t_proj;
    }
  } else {
    if (tIsVisible) {
      *t_vis = t;
      vis = &s_vis;
      a = t;
      a_proj = t_proj;
      b = s;
      b_proj = s_proj;
    } else {
      return false;
    }
  }

  // binary search with starting range [a,b].  Search until a, b are within same pixel
  while ((abs(a_proj[0] - b_proj[0]) + abs(a_proj[1] - b_proj[1])) >= 1.f) {
    Vec3f c = 0.5f * (a + b);
    Vec3f c_proj = toImagePlane(c);
    float c_depth = c_proj[2];
    bool cIsVisible = (depthBuffer[(GLint)c_proj[1]*width + (GLint)c_proj[0]] - c_depth) > -EPSILON;
    if (cIsVisible) {
      a = c;
      a_proj = c_proj;
    } else {
      b = c;
      b_proj = c_proj;
    }
  }

  **vis = a;
  return true;
}


void writeImage(Mesh &mesh, int width, int height, string filename, Vec3f camPos,
                EPropHandleT<void*>& edgeData) {
  // copy entire depth buffer
  vector<GLfloat> depthBuffer(width * height);
  glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, &depthBuffer[0]);

  ofstream outfile(filename.c_str());
  outfile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
  outfile << "<svg width=\"5in\" height=\"5in\" viewBox=\"0 0 " << width << ' ' << height << "\">\n";
  outfile << "<g stroke=\"black\" fill=\"black\">\n";

  /*// Sample code for generating image of the entire triangle mesh:
  for (Mesh::ConstEdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
    Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(it,0);
    Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(it,1);
    Vec3f source(mesh.point(mesh.from_vertex_handle(h0)));
    Vec3f target(mesh.point(mesh.from_vertex_handle(h1)));

    //if (!isVisible(source, &depthBuffer[0], width) || !isVisible(target, &depthBuffer[0], width)) continue;
    Vec3f source_vis, target_vis;
    if (!isEdgePartialVisible(source, target, &depthBuffer[0], width, &source_vis, &target_vis)) continue;

    Vec3f p1 = toImagePlane(source_vis);
    Vec3f p2 = toImagePlane(target_vis);
    outfile << "<line ";
    outfile << "x1=\"" << p1[0] << "\" ";
    outfile << "y1=\"" << height-p1[1] << "\" ";
    outfile << "x2=\"" << p2[0] << "\" ";
    outfile << "y2=\"" << height-p2[1] << "\" stroke-width=\"1\" />\n";
  }*/

  struct SilhouetteLink {
    deque<Mesh::HalfedgeHandle> halfEdges;  // halfedges oriented and ordered from start to end
    Mesh::VertexHandle start, end;
  };

  vector<SilhouetteLink> silhouetteLinks;
  Mesh::EdgeIter e_it, e_it_end = mesh.edges_end();
  for (e_it = mesh.edges_begin(); e_it != e_it_end; ++e_it) {
    Mesh::EdgeHandle eh = *e_it;
    if (!isSilhouette(mesh, eh, camPos)) {
      continue;
    }
    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
    Mesh::VertexHandle s = mesh.from_vertex_handle(heh);
    Mesh::VertexHandle t = mesh.to_vertex_handle(heh);

    SilhouetteLink* firstAttachedLink = NULL;    // This will be updated to the first link the edge is attached to (if any)
    Mesh::VertexHandle unattachedVertex;  // This will be updated to the unattached vertex of the edge after the edge is first attached to a link (if any)
    bool unattachedIsStart = false;       // This will be updated at the same time as unattachedVertex; will indicate if unattachedVertex is the start or end of the link the edge attached to.
    
    auto links_it = silhouetteLinks.begin();
    while(links_it != silhouetteLinks.end()) {
      SilhouetteLink& link = *links_it;
      if (link.start == link.end) {   // skip any links that are loops
        continue;
      }
      // It's possible that s=start and t=end, or s=end and t=start.  In those cases, the edge will
      // be attached at s, and no further action is required.
      if (s == link.start) {
        link.halfEdges.push_front(mesh.opposite_halfedge_handle(heh));
        link.start = t;
        unattachedVertex = t;
        unattachedIsStart = true;
        firstAttachedLink = &link;
        break;
      } else if (s == link.end) {
        link.halfEdges.push_back(heh);
        link.end = t;
        unattachedVertex = t;
        unattachedIsStart = false;
        firstAttachedLink = &link;
        break;
      } else if (t == link.start) {
        link.halfEdges.push_front(heh);
        link.start = s;
        unattachedVertex = s;
        unattachedIsStart = true;
        firstAttachedLink = &link;
        break;
      } else if (t == link.end) {
        link.halfEdges.push_back(mesh.opposite_halfedge_handle(heh));
        link.end = s;
        unattachedVertex = s;
        unattachedIsStart = false;
        firstAttachedLink = &link;
        break;
      }
      ++links_it;
    }
    ++links_it;

    // At this point, if i < silhouetteLinks.size(), then the edge has been attached to a link.  We'll
    // now check if the unattached end of the edge connects to any of the remaining links.
    while (links_it != silhouetteLinks.end()) {
      SilhouetteLink& link = *links_it;
      if (link.start == link.end) {   // skip any links that are loops
        continue;
      }
      if (unattachedVertex == link.start) {
        if (unattachedIsStart) {
          // firstAttachedLink <-<-<- unattachedVertex ->->->-> link
          for (int i = 0; i < link.halfEdges.size(); i++) {
            Mesh::HalfedgeHandle& link_heh = link.halfEdges[i];
            mesh.property(edgeData, mesh.edge_handle(link_heh)) = reinterpret_cast<void*>(firstAttachedLink);
            firstAttachedLink->halfEdges.push_front(mesh.opposite_halfedge_handle(link_heh));
          }
          firstAttachedLink->start = link.end;
          silhouetteLinks.erase(links_it);
        } else {
          // firstAttachedLink ->->-> unattachedVertex ->->->-> link
          for (int i = 0; i < link.halfEdges.size(); i++) {
            Mesh::HalfedgeHandle& link_heh = link.halfEdges[i];
            mesh.property(edgeData, mesh.edge_handle(link_heh)) = reinterpret_cast<void*>(firstAttachedLink);
            firstAttachedLink->halfEdges.push_back(link_heh);
          }
          firstAttachedLink->end = link.end;
          silhouetteLinks.erase(links_it);
        }
      } else if (unattachedVertex == link.end) {
        if (unattachedIsStart) {
          // firstAttachedLink <-<-<- unattachedVertex <-<-<-<- link
          for (int i = link.halfEdges.size() - 1; i >= 0; i--) {
            Mesh::HalfedgeHandle& link_heh = link.halfEdges[i];
            mesh.property(edgeData, mesh.edge_handle(link_heh)) = reinterpret_cast<void*>(firstAttachedLink);
            firstAttachedLink->halfEdges.push_front(link_heh);
          }
          firstAttachedLink->start = link.start;
          silhouetteLinks.erase(links_it);
        } else {
          // firstAttachedLink ->->-> unattachedVertex <-<-<-<- link
          for (int i = link.halfEdges.size() - 1; i >= 0; i--) {
            Mesh::HalfedgeHandle& link_heh = link.halfEdges[i];
            mesh.property(edgeData, mesh.edge_handle(link_heh)) = reinterpret_cast<void*>(firstAttachedLink);
            firstAttachedLink->halfEdges.push_back(mesh.opposite_halfedge_handle(link_heh));
          }
          firstAttachedLink->end = link.start;
          silhouetteLinks.erase(links_it);
        }
      }
      else {
        ++links_it;
      }
    }

    // if this edge doesn't touch any existing link, create a new link from it
    if (firstAttachedLink == NULL) {
      silhouetteLinks.push_back(SilhouetteLink());
      SilhouetteLink& newLink = silhouetteLinks.back();
      newLink.halfEdges.push_back(heh);
      newLink.start = s;
      newLink.end = t;
      firstAttachedLink = &newLink;
    }

    mesh.property(edgeData, eh) = reinterpret_cast<void*>(firstAttachedLink);
  }

  // WRITE CODE HERE TO GENERATE A .SVG OF THE MESH --------------------------------------------------------------

  // -------------------------------------------------------------------------------------------------------------

  outfile << "</g>\n";
  outfile << "</svg>\n";
  outfile.close();
}
