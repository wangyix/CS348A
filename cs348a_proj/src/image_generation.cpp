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

const GLdouble EPSILON = 0.01;

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

bool isVisible(const Vec3f& point_proj, const GLfloat* depthBuffer, int width, int height) {
  if (!(0 <= point_proj[0] && point_proj[0] < width && 0 <= point_proj[1] && point_proj[1] < height)) return false;
  return (depthBuffer[((int)point_proj[1]) * width + (int)point_proj[0]] - point_proj[2]) > -EPSILON;

  /*GLfloat bufDepth = 0.0;
  glReadPixels(static_cast<GLint>(point_proj[0]), static_cast<GLint>(point_proj[1]),
    1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &bufDepth);
  return bufDepth - point_proj[2] > -EPSILON;*/
}

bool isVisible(const Vec2f& point_proj_xy, float point_depth, const GLfloat* depthBuffer, int width, int height) {
  if (!(0 <= point_proj_xy[0] && point_proj_xy[0] < width && 0 <= point_proj_xy[1] && point_proj_xy[1] < height)) return false;
  return (depthBuffer[((int)point_proj_xy[1]) * width + (int)point_proj_xy[0]] - point_depth) > -EPSILON;

  /*GLfloat bufDepth = 0.0;
  glReadPixels(static_cast<GLint>(point_proj_xy[0]), static_cast<GLint>(point_proj_xy[1]),
    1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &bufDepth);
  return bufDepth - point_depth > -EPSILON;*/
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




struct Link {
  deque<Mesh::HalfedgeHandle> halfEdges;  // halfedges oriented and ordered from start to end
  Vec3f start, end;
};


void writeImage(Mesh &mesh, int width, int height, string filename,
                const Vec3f& camPos, const Vec3f& camLookDir,
                EPropHandleT<bool>& edgeVisited) {

  // copy entire depth buffer
  vector<GLfloat> depthBuffer(width * height);
  glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, &depthBuffer[0]);

  vector<Link> silhouetteLinks;

  // clear all edgeVisited flags
  Mesh::EdgeIter e_it, e_it_end = mesh.edges_end();
  for (e_it = mesh.edges_begin(); e_it != e_it_end; ++e_it) {
    mesh.property(edgeVisited, *e_it) = false;
  }

  // Visit all edges and group them into visible links.
  for (e_it = mesh.edges_begin(); e_it != e_it_end; ++e_it) {
    Mesh::EdgeHandle eh = *e_it;
    if (!isFeatureEdge(mesh, eh, camPos) || mesh.property(edgeVisited, eh)) {
      continue;
    }

    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
    Mesh::VertexHandle s_vh = mesh.from_vertex_handle(heh);
    Mesh::VertexHandle t_vh = mesh.to_vertex_handle(heh);
    Vec3f s = mesh.point(s_vh);
    Vec3f t = mesh.point(t_vh);

    silhouetteLinks.push_back(Link());
    Link& newLink = silhouetteLinks.back();
    newLink.start = s;
    newLink.end = t;
    newLink.halfEdges.push_back(heh);
    mesh.property(edgeVisited, eh) = true;

    // extend the end of the link as much as possible
    Mesh::VertexHandle newLinkEndVh = t_vh;
    while (true) {
      // Try to find an outgoing halfedge from the current link end that's a silhouette edge
      Mesh::VertexOHalfedgeCCWIter voh_it, voh_it_end = mesh.voh_ccwend(newLinkEndVh);
      for (voh_it = mesh.voh_ccwbegin(newLinkEndVh); voh_it != voh_it_end; ++voh_it) {
        Mesh::HalfedgeHandle heh = *voh_it;
        Mesh::EdgeHandle eh = mesh.edge_handle(heh);
        if (!mesh.property(edgeVisited, eh) && isFeatureEdge(mesh, eh, camPos)) {
          newLinkEndVh = mesh.to_vertex_handle(heh);
          newLink.end = mesh.point(newLinkEndVh);
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
    Mesh::VertexHandle newLinkStartVh = s_vh;
    while (true) {
      // Try to find an incoming halfedge from the current link start that's a silhouette edge
      Mesh::VertexIHalfedgeCCWIter vih_it, vih_it_end = mesh.vih_ccwend(newLinkStartVh);
      for (vih_it = mesh.vih_ccwbegin(newLinkStartVh); vih_it != vih_it_end; ++vih_it) {
        Mesh::HalfedgeHandle heh = *vih_it;
        Mesh::EdgeHandle eh = mesh.edge_handle(heh);
        if (!mesh.property(edgeVisited, eh) && isFeatureEdge(mesh, eh, camPos)) {
          newLinkStartVh = mesh.from_vertex_handle(heh);
          newLink.start = mesh.point(newLinkStartVh);
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

  vector<Link> visibleLinks;

  for (Link& silhouetteLink : silhouetteLinks) {

    bool startVisible = isVisible(toImagePlane(silhouetteLink.start), &depthBuffer[0], width, height);
    int firstLinkAt = -1;
    if (startVisible) {
      firstLinkAt = visibleLinks.size();
      visibleLinks.push_back(Link());
      Link& newLink = visibleLinks.back();
      newLink.start = silhouetteLink.start;
      newLink.end = newLink.start;
      newLink.halfEdges.push_back(silhouetteLink.halfEdges.front());
    }

    Vec3f prevVertex = silhouetteLink.start;
    bool prevVertexVisible = startVisible;

    // Walk through link edges, recording visible portions as we encounter them
    for (Mesh::HalfedgeHandle heh : silhouetteLink.halfEdges) {
      Mesh::VertexHandle s_vh = mesh.from_vertex_handle(heh);
      Mesh::VertexHandle t_vh = mesh.to_vertex_handle(heh);
      Vec3f s = mesh.point(s_vh);
      Vec3f t = mesh.point(t_vh);
      float s_depth, t_depth;
      Vec2f s_proj = toImagePlane(s, &s_depth);
      Vec2f t_proj = toImagePlane(t, &t_depth);
      //bool sVisible = isVisible(s_proj, s_depth, &depthBuffer[0], width, height);
      //bool tVisible = isVisible(t_proj, t_depth, &depthBuffer[0], width, height);
      float w_s = (s - camPos) | camLookDir;
      float w_t = (t - camPos) | camLookDir;

      // split edge into segments of equal length in screen-space.
      const float MAX_EDGE_SEG_LENGTH = 4.f;  // measured in pixels in screenspace
      int numSegments = ceilf((t_proj - s_proj).length() / MAX_EDGE_SEG_LENGTH);
      for (int i = 1; i <= numSegments; i++) {
        float b = i / (float)(numSegments);
        float a = b*w_s / ((1.f - b)*w_t + b*w_s);
        Vec3f vertex = (1.f - a)*s + a*t;
        bool vertexVisible = isVisible(toImagePlane(vertex), &depthBuffer[0], width, height);

        /*float unused;               // debug code to check screen-space lengths of segments
        Vec2f vertexProj = toImagePlane(vertex, &unused);
        Vec2f prevVertexProj = toImagePlane(prevVertex, &unused);
        float d = (vertexProj - prevVertexProj).norm();*/

        if (prevVertexVisible) {
          Vec3f newEnd;
          if (vertexVisible) {
            newEnd = vertex;
          } else {
            newEnd = findVisibilityBoundary(prevVertex, vertex, &depthBuffer[0], width, height);
          }
          // Extend the current link
          Link& currentLink = visibleLinks.back();
          currentLink.end = vertex;
          if (heh != currentLink.halfEdges.back()) {
            currentLink.halfEdges.push_back(heh);
          }
          /*visibleLinks.push_back(Link());     // debug code; for visiualizing individual segments
          Link& newLink = visibleLinks.back();
          newLink.start = prevVertex;
          newLink.end = newEnd;
          newLink.halfEdges.push_back(heh);*/
        } else {
          if (vertexVisible) {
            // Start a new link
            visibleLinks.push_back(Link());
            Link& newLink = visibleLinks.back();
            newLink.start = findVisibilityBoundary(vertex, prevVertex, &depthBuffer[0], width, height);;
            newLink.end = vertex;
            newLink.halfEdges.push_back(heh);
          }
        }
        prevVertex = vertex;
        prevVertexVisible = vertexVisible;
      }
    }

    if (silhouetteLink.start == silhouetteLink.end && startVisible) {
      // The silhouette link is a loop; the first and last visible links need to be joined into one link
      // (unless the whole loop is visible, in which case the first and last links are already the same link).
      // Append the edges in the last visible link to the front of the edges in the first visible link,
      // then remove the last visible link from visibleLinks.
      assert(firstLinkAt >= 0);
      if (firstLinkAt != visibleLinks.size() - 1) {
        Link& firstLink = visibleLinks[firstLinkAt];
        Link& lastLink = visibleLinks.back();
        deque<Mesh::HalfedgeHandle>& firstLinkEdges = firstLink.halfEdges;
        deque<Mesh::HalfedgeHandle>& lastLinkEdges = lastLink.halfEdges;
        firstLinkEdges.insert(firstLinkEdges.begin(), lastLinkEdges.begin(), lastLinkEdges.end());
        firstLink.start = lastLink.start;
        visibleLinks.pop_back();
      }
    }
  }

  // sanity check!!!
  for (Link& link : visibleLinks) {
    for (int i = 1; i < link.halfEdges.size(); ++i) {
      assert(mesh.from_vertex_handle(link.halfEdges[i]) == mesh.to_vertex_handle(link.halfEdges[i - 1]));
    }
  }

  ofstream outfile(filename.c_str());
  outfile << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
  outfile << "<svg width=\"5in\" height=\"5in\" viewBox=\"0 0 " << width << ' ' << height << "\">" << endl;

  char* strokeStrings[8] = { "<g stroke=\"black\" fill=\"black\">",
                              "<g stroke=\"red\" fill=\"black\">",
                              "<g stroke=\"blue\" fill=\"black\">",
                              "<g stroke=\"magenta\" fill=\"black\">",
                              "<g stroke=\"orange\" fill=\"black\">",
                              "<g stroke=\"deeppink\" fill=\"black\">",
                              "<g stroke=\"olivedrab\" fill=\"black\">",
                              "<g stroke=\"deepskyblue\" fill=\"black\">" };

  for (int i = 0; i < visibleLinks.size(); i++) {
    Link& link = visibleLinks[i];

    outfile << strokeStrings[i % 8] << endl;

    for (int i = 0; i < link.halfEdges.size(); i++) {
      Mesh::HalfedgeHandle heh = link.halfEdges[i];
      Vec3f source = (i == 0 ? link.start : mesh.point(mesh.from_vertex_handle(heh)));
      Vec3f target = (i == link.halfEdges.size() - 1 ? link.end : mesh.point(mesh.to_vertex_handle(heh)));
      //Vec3f source = mesh.point(mesh.from_vertex_handle(heh));
      //Vec3f target = mesh.point(mesh.to_vertex_handle(heh));

      Vec3f source_vis = source, target_vis = target;

      Vec3f p1 = toImagePlane(source_vis);
      Vec3f p2 = toImagePlane(target_vis);
      outfile << "\t<line ";
      outfile << "x1=\"" << p1[0] << "\" ";
      outfile << "y1=\"" << height - p1[1] << "\" ";
      outfile << "x2=\"" << p2[0] << "\" ";
      outfile << "y2=\"" << height - p2[1] << "\" stroke-width=\"1\" />" << endl;
    }

    outfile << "</g>" << endl;
  }
  
  outfile << "</svg>" << endl;
  outfile.close();

  std::cout << filename << " written" << std::endl;
}
