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
  deque<Vec3f> vertices;
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
    newLink.vertices.push_back(s);
    newLink.vertices.push_back(t);
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
          newLink.vertices.push_back(mesh.point(newLinkEndVh));
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
          newLink.vertices.push_front(mesh.point(newLinkStartVh));
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

    Vec3f start = silhouetteLink.vertices.front();
    bool startVisible = isVisible(toImagePlane(start), &depthBuffer[0], width, height);
    int firstLinkAt = -1;
    if (startVisible) {
      firstLinkAt = visibleLinks.size();
      visibleLinks.push_back(Link());
      visibleLinks.back().vertices.push_back(start);
    }

    Vec3f prevSamplePos = start;
    bool prevSamplePosVisible = startVisible;

    // Walk through link edges, recording visible portions as we encounter them
    for (int tIndex = 1; tIndex < silhouetteLink.vertices.size(); ++tIndex) {
      Vec3f s = silhouetteLink.vertices[tIndex - 1];
      Vec3f t = silhouetteLink.vertices[tIndex];
      float s_depth, t_depth;
      Vec2f s_proj = toImagePlane(s, &s_depth);
      Vec2f t_proj = toImagePlane(t, &t_depth);
      float w_s = (s - camPos) | camLookDir;
      float w_t = (t - camPos) | camLookDir;

      // split edge into segments of equal length in screen-space.
      const float MAX_EDGE_SEG_LENGTH = 4.f;  // measured in pixels in screenspace
      int numSegments = ceilf((t_proj - s_proj).length() / MAX_EDGE_SEG_LENGTH);
      bool linkEndsOnThisEdge = false;
      for (int i = 1; i <= numSegments; i++) {
        float b = i / (float)(numSegments);
        float a = b*w_s / ((1.f - b)*w_t + b*w_s);
        Vec3f samplePos = (1.f - a)*s + a*t;
        bool samplePosVisible = isVisible(toImagePlane(samplePos), &depthBuffer[0], width, height);

        /*float unused;               // debug code to check screen-space lengths of segments
        Vec2f vertexProj = toImagePlane(vertex, &unused);
        Vec2f prevVertexProj = toImagePlane(prevVertex, &unused);
        float d = (vertexProj - prevVertexProj).norm();*/

        if (prevSamplePosVisible) {
          Vec3f newEnd;
          if (samplePosVisible) {
            newEnd = samplePos;
          } else {
            newEnd = findVisibilityBoundary(prevSamplePos, samplePos, &depthBuffer[0], width, height);
          }
          // Extend the current link
          if (linkEndsOnThisEdge) {
            visibleLinks.back().vertices.back() = newEnd;
          } else {
            visibleLinks.back().vertices.push_back(newEnd);
            linkEndsOnThisEdge = true;
          }
          /*visibleLinks.push_back(Link());       // debug code for visualizing the segmented edges
          visibleLinks.back().vertices.push_back(prevSamplePos);
          visibleLinks.back().vertices.push_back(newEnd);
          linkEndsOnThisEdge = true;*/
        } else {
          if (samplePosVisible) {
            // Start a new link
            visibleLinks.push_back(Link());
            Vec3f newStart = findVisibilityBoundary(samplePos, prevSamplePos, &depthBuffer[0], width, height);
            visibleLinks.back().vertices.push_back(newStart);
            visibleLinks.back().vertices.push_back(samplePos);
            linkEndsOnThisEdge = true;
          }
        }
        prevSamplePos = samplePos;
        prevSamplePosVisible = samplePosVisible;
      }
    }

    if (start == silhouetteLink.vertices.back() && startVisible) {
      // The silhouette link is a loop; the first and last visible links need to be joined into one link
      // (unless the whole loop is visible, in which case the first and last links are already the same link).
      // Append the edges in the last visible link to the front of the edges in the first visible link,
      // then remove the last visible link from visibleLinks.
      assert(0 <= firstLinkAt && firstLinkAt < visibleLinks.size());
      if (firstLinkAt != visibleLinks.size() - 1) {
        Link& firstLink = visibleLinks[firstLinkAt];
        Link& lastLink = visibleLinks.back();
        deque<Vec3f>& firstLinkVertices = firstLink.vertices;
        deque<Vec3f>& lastLinkVertices = lastLink.vertices;
        firstLinkVertices.insert(firstLinkVertices.begin(), lastLinkVertices.begin(), lastLinkVertices.end() - 1);
        visibleLinks.pop_back();
      }
    }
  }

  // sanity check!!!
  for (Link& link : visibleLinks) {
    for (int i = 1; i < link.vertices.size(); ++i) {
      assert(link.vertices[i - 1] != link.vertices[i]);
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

    for (int i = 1; i < link.vertices.size(); i++) {
      Vec3f source = link.vertices[i - 1];
      Vec3f dest = link.vertices[i];
      Vec3f p1 = toImagePlane(source);
      Vec3f p2 = toImagePlane(dest);
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
