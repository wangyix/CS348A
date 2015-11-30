#include "image_generation.h"
#include "mesh_features.h"
#include <GL/glut.h>
#include <fstream>
#include <iostream>
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

Vec2f toImagePlane(Vec3f point, float* depth) {
  GLdouble point3DX = point[0], point3DY = point[1], point3DZ = point[2];

  GLdouble modelMatrix[16], projMatrix[16];
  GLint viewport[4];
  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glGetIntegerv(GL_VIEWPORT, viewport);

  GLdouble point2DX, point2DY, point2DZ;
  gluProject(point3DX, point3DY, point3DZ, modelMatrix, projMatrix, viewport, &point2DX, &point2DY, &point2DZ);

  *depth = point2DZ;
  return Vec2f(point2DX, point2DY);
}

const GLdouble EPSILON = 0.0075;

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

/*bool isVisible(const Vec3f& point, const GLfloat* depthBuffer, int width, int height, Vec3f* projected) {
Vec3f proj = toImagePlane(point);
*projected = proj;
if (!(0 <= proj[0] && proj[0] < width && 0 <= proj[1] && proj[1] < height)) return false;
return (depthBuffer[((int)proj[1]) * width + (int)proj[0]] - proj[2]) > -EPSILON;
}*/
/*// Returns true if edge is wholly or partially visible. If only partially visible, the occluded endpoint
// will be adjusted along the edge to the boundary between the visible and occluded portions of the edge.
bool isEdgePartialVisible(const Vec3f& s, const Vec3f& t, const GLfloat* depthBuffer, int width, int height,
                          bool* s_visible, bool* t_visible, Vec3f* s_adj, Vec3f* t_adj) {
  Vec3f s_proj, t_proj;
  *s_visible = isVisible(s, &depthBuffer[0], width, height, &s_proj);
  *t_visible = isVisible(t, &depthBuffer[0], width, height, &t_proj);
  
  Vec3f a, b;   // a will be visible bound, b will be occluded bound
  Vec3f a_proj, b_proj;
  Vec3f** adj;  // ptr to the endpoint that will be adjusted to the result of the binary search
  if (*s_visible) {
    if (*t_visible) {
      *s_adj = s;
      *t_adj = t;
      return true;
    } else {
      *s_adj = s;
      adj = &t_adj;
      a = s;
      a_proj = s_proj;
      b = t;
      b_proj = t_proj;
    }
  } else {
    if (*t_visible) {
      *t_adj = t;
      adj = &s_adj;
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
    Vec3f c_proj;
    bool cIsVisible = isVisible(c, &depthBuffer[0], width, height, &c_proj);
    if (cIsVisible) {
      a = c;
      a_proj = c_proj;
    } else {
      b = c;
      b_proj = c_proj;
    }
  }

  **adj = a;
  return true;
}*/

bool isVisible(const Vec3f& point_proj, const GLfloat* depthBuffer, int width, int height) {
  if (!(0 <= point_proj[0] && point_proj[0] < width && 0 <= point_proj[1] && point_proj[1] < height)) return false;
  return (depthBuffer[((int)point_proj[1]) * width + (int)point_proj[0]] - point_proj[2]) > -EPSILON;
}

bool isVisible(const Vec2f& point_proj_xy, float point_depth, const GLfloat* depthBuffer, int width, int height) {
  if (!(0 <= point_proj_xy[0] && point_proj_xy[0] < width && 0 <= point_proj_xy[1] && point_proj_xy[1] < height)) return false;
  return (depthBuffer[((int)point_proj_xy[1]) * width + (int)point_proj_xy[0]] - point_depth) > -EPSILON;
}

Vec3f findVisibilityBoundary(const Vec3f& visible, const Vec3f& occluded, const GLfloat* depthBuffer, int width, int height) {
  Vec3f a = visible;
  Vec3f b = occluded;
  float unused;
  Vec2f a_proj = toImagePlane(a, &unused);
  Vec2f b_proj = toImagePlane(b, &unused);
  while ((b_proj - a_proj).sqrnorm() >= 1.f) {
    Vec3f c = 0.5f * (a + b);
    float c_depth;
    Vec2f c_proj = toImagePlane(c, &c_depth);
    bool cIsVisible = isVisible(c_proj, c_depth,  &depthBuffer[0], width, height);
    if (cIsVisible) {
      a = c;
      a_proj = c_proj;
    } else {
      b = c;
      b_proj = c_proj;
    }
  }
  return 0.5f * (a + b);
}


void writeImage(Mesh &mesh, int width, int height, string filename,
                const Vec3f& camPos, const Vec3f& camLookDir,
                EPropHandleT<bool>& edgeVisited) {

  // copy entire depth buffer
  vector<GLfloat> depthBuffer(width * height);
  glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, &depthBuffer[0]);


  // We want to group visible feature edges into links.  All edges in a visible link must be fully
  // visible except the links on the ends: those may be partially visible. Specifically, the mesh
  // vertices that the link starts/ends on may not be visible, but all other vertices must be visible.
  //
  // o: vertex
  // ==: visible portion of edge
  // --: occluded portion of edge
  //
  // Both ends occluded:
  //
  // start  o---|==>o=====>o=====>o....o=====>o===|-->o  end
  //            ^                                 ^
  //           pos                               pos
  //
  // One end occluded:
  //
  // start  o======>o=====>o=====>o....o=====>o===|-->o  end
  //        ^                                     ^
  //       pos                                   pos
  //
  // The LinkEndpoint struct describes one end of a link: vh points to the mesh vertex that this end
  // is attached to. If the vertex is not visible, then pos is the position on the link edge that's on
  // the boundary between the visible and occluded portions of the edge. If vertex is visible, then pos
  // is simply the position of the mesh vertex.
  // The VisibleLink struct describes a visible link. It's considered complete when no more edges can
  // be added to either end, which can occur in two ways: the link forms a loop, or both ends of the link
  // are occluded.

  struct LinkEndpoint {
    Mesh::VertexHandle vh;
    bool isVisible;
    Vec3f pos;
  };
  struct Link {
    deque<Mesh::HalfedgeHandle> halfEdges;  // halfedges oriented and ordered from start to end
    LinkEndpoint start, end;
  };


  vector<Link> silhouetteLinks;

  // clear all edgeVisited flags
  Mesh::EdgeIter e_it, e_it_end = mesh.edges_end();
  for (e_it = mesh.edges_begin(); e_it != e_it_end; ++e_it) {
    mesh.property(edgeVisited, *e_it) = false;
  }

  // Visit all edges and group them into visible links.
  for (e_it = mesh.edges_begin(); e_it != e_it_end; ++e_it) {
    Mesh::EdgeHandle eh = *e_it;

    if (!isSilhouette(mesh, eh, camPos)) {
      continue;
    }
    if (mesh.property(edgeVisited, eh)) {
      continue;
    }

    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);

    LinkEndpoint s, t;
    s.vh = mesh.from_vertex_handle(heh);
    t.vh = mesh.to_vertex_handle(heh);
    // isVisible, pos are not used yet

    silhouetteLinks.push_back(Link());
    Link& newLink = silhouetteLinks.back();
    newLink.start = s;
    newLink.end = t;
    newLink.halfEdges.push_back(heh);
    mesh.property(edgeVisited, eh) = true;

    // extend the end of the link as much as possible
    while (true) {
      // Try to find an outgoing halfedge from the current link end that's a silhouette edge
      Mesh::VertexOHalfedgeCCWIter voh_it, voh_it_end = mesh.voh_ccwend(newLink.end.vh);
      for (voh_it = mesh.voh_ccwbegin(newLink.end.vh); voh_it != voh_it_end; ++voh_it) {
        Mesh::HalfedgeHandle heh = *voh_it;
        Mesh::EdgeHandle eh = mesh.edge_handle(heh);
        if (!mesh.property(edgeVisited, eh) && isSilhouette(mesh, eh, camPos)) {
          newLink.end.vh = mesh.to_vertex_handle(heh);
          newLink.halfEdges.push_back(heh);
          mesh.property(edgeVisited, eh) = true;
          break;
        }
      }
      if (voh_it == voh_it_end) {
        // No such halfedge was found; the link end cannot be extended further.
        break;
      }
    }
    // extend the start of the link as much as possible
    while (true) {
      // Try to find an incoming halfedge from the current link start that's a silhouette edge
      Mesh::VertexIHalfedgeCCWIter vih_it, vih_it_end = mesh.vih_ccwend(newLink.start.vh);
      for (vih_it = mesh.vih_ccwbegin(newLink.start.vh); vih_it != vih_it_end; ++vih_it) {
        Mesh::HalfedgeHandle heh = *vih_it;
        Mesh::EdgeHandle eh = mesh.edge_handle(heh);
        if (!mesh.property(edgeVisited, eh) && isSilhouette(mesh, eh, camPos)) {
          newLink.start.vh = mesh.from_vertex_handle(heh);
          newLink.halfEdges.push_front(heh);
          mesh.property(edgeVisited, eh) = true;
          break;
        }
      }
      if (vih_it == vih_it_end) {
        // No such halfedge was found; the link start cannot be extended further
        break;
      }
    }
  }


  /*// sanity check!!!
  for (Link& link : visibleLinks) {
    assert(link.isComplete);
    if (link.start.isVisible) assert(link.start.pos == mesh.point(link.start.vh));
    if (link.end.isVisible) assert(link.end.pos == mesh.point(link.end.vh));
    Mesh::VertexHandle prevT = link.start.vh;
    for (Mesh::HalfedgeHandle& heh : link.halfEdges) {
      assert(prevT == mesh.from_vertex_handle(heh));
      prevT = mesh.to_vertex_handle(heh);
    }
    assert(prevT == link.end.vh);
  }*/

  ofstream outfile(filename.c_str());
  outfile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
  outfile << "<svg width=\"5in\" height=\"5in\" viewBox=\"0 0 " << width << ' ' << height << "\">\n";

  string strokeStrings[8] = { "<g stroke=\"black\" fill=\"black\">\n",
                              "<g stroke=\"red\" fill=\"black\">\n",
                              "<g stroke=\"blue\" fill=\"black\">\n",
                              "<g stroke=\"magenta\" fill=\"black\">\n",
                              "<g stroke=\"orange\" fill=\"black\">\n",
                              "<g stroke=\"deeppink\" fill=\"black\">\n",
                              "<g stroke=\"olivedrab\" fill=\"black\">\n"
                              "<g stroke=\"deepskyblue\" fill=\"black\">\n" };

  for (int i = 0; i < silhouetteLinks.size(); i++) {
    Link& link = silhouetteLinks[i];

    outfile << strokeStrings[i % 8];

    // Sample code for generating image of the entire triangle mesh:
    /*for (Mesh::ConstEdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
      Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(it, 0);
      Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(it, 1);
      Vec3f source(mesh.point(mesh.from_vertex_handle(h0)));
      Vec3f target(mesh.point(mesh.from_vertex_handle(h1)));*/

    //for (Mesh::HalfedgeHandle& heh : link.halfEdges) {
    for (int i = 0; i < link.halfEdges.size(); i++) {
      Mesh::HalfedgeHandle heh = link.halfEdges[i];
      //Vec3f source = (i == 0 ? link.start.pos : mesh.point(mesh.from_vertex_handle(heh)));
      //Vec3f target = (i == link.halfEdges.size() - 1 ? link.end.pos : mesh.point(mesh.to_vertex_handle(heh)));
      Vec3f source = mesh.point(mesh.from_vertex_handle(heh));
      Vec3f target = mesh.point(mesh.to_vertex_handle(heh));

      //if (!isVisible(source, &depthBuffer[0], width) || !isVisible(target, &depthBuffer[0], width)) continue;
      //Vec3f source_vis, target_vis;
      //if (!isEdgePartialVisible(source, target, &depthBuffer[0], width, &source_vis, &target_vis)) continue;
      Vec3f source_vis = source, target_vis = target;

      Vec3f p1 = toImagePlane(source_vis);
      Vec3f p2 = toImagePlane(target_vis);
      outfile << "<line ";
      outfile << "x1=\"" << p1[0] << "\" ";
      outfile << "y1=\"" << height - p1[1] << "\" ";
      outfile << "x2=\"" << p2[0] << "\" ";
      outfile << "y2=\"" << height - p2[1] << "\" stroke-width=\"1\" />\n";
    }

    outfile << "</g>\n";
  }
  
  outfile << "</svg>\n";
  outfile.close();

  std::cout << filename << " written" << std::endl;
}
