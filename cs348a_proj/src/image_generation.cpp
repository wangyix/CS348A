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

static Vec3f toImagePlane(Vec3f point) {
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

static Vec2f toImagePlane(Vec3f point, float* depth) {
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

static void toImagePlaneInvertHeight(const deque<Vec3f>& points, int height, vector<Vec2f>* projPoints) {
  GLdouble modelMatrix[16], projMatrix[16];
  GLint viewport[4];
  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glGetIntegerv(GL_VIEWPORT, viewport);

  projPoints->resize(points.size());
  for (int i = 0; i < points.size(); ++i) {
    GLdouble point3DX = points[i][0], point3DY = points[i][1], point3DZ = points[i][2];
    GLdouble point2DX, point2DY, point2DZ;
    gluProject(point3DX, point3DY, point3DZ, modelMatrix, projMatrix, viewport, &point2DX, &point2DY, &point2DZ);
    projPoints->at(i) = Vec2f(point2DX, height - 1 - point2DY);
  }
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

static bool isVisible(const Vec3f& point_proj, const GLfloat* depthBuffer, int width, int height, float epsilon) {
  if (!(0 <= point_proj[0] && point_proj[0] < width && 0 <= point_proj[1] && point_proj[1] < height)) return false;
  return (depthBuffer[((int)point_proj[1]) * width + (int)point_proj[0]] - point_proj[2]) > -epsilon;

  /*GLfloat bufDepth = 0.0;
  glReadPixels(static_cast<GLint>(point_proj[0]), static_cast<GLint>(point_proj[1]),
    1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &bufDepth);
  return bufDepth - point_proj[2] > -EPSILON;*/
}

static bool isVisible(const Vec2f& point_proj_xy, float point_depth, const GLfloat* depthBuffer, int width, int height, float epsilon) {
  if (!(0 <= point_proj_xy[0] && point_proj_xy[0] < width && 0 <= point_proj_xy[1] && point_proj_xy[1] < height)) return false;
  return (depthBuffer[((int)point_proj_xy[1]) * width + (int)point_proj_xy[0]] - point_depth) > -epsilon;

  /*GLfloat bufDepth = 0.0;
  glReadPixels(static_cast<GLint>(point_proj_xy[0]), static_cast<GLint>(point_proj_xy[1]),
    1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &bufDepth);
  return bufDepth - point_depth > -EPSILON;*/
}

static Vec3f findVisibilityBoundary(const Vec3f& visible, const Vec3f& occluded, const GLfloat* depthBuffer, int width, int height, float epsilon) {
  Vec3f a = visible;
  Vec3f b = occluded;
  float unused;
  Vec2f a_proj = toImagePlane(a, &unused);
  Vec2f b_proj = toImagePlane(b, &unused);
  while ((b_proj - a_proj).sqrnorm() >= 1.f) {
    Vec3f c = 0.5f * (a + b);
    float c_depth;
    Vec2f c_proj = toImagePlane(c, &c_depth);
    bool cIsVisible = isVisible(c_proj, c_depth, &depthBuffer[0], width, height, epsilon);
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


static Mesh::VertexHandle extendFeatureLinkEnd(Mesh& mesh, Mesh::VertexHandle linkEnd, deque<Vec3f>* linkVertices,
                                              void (deque<Vec3f>::*dequePushFunc)(const Vec3f&), const Vec3f& camPos, 
                                              EPropHandleT<bool>& edgeVisited) {
  while (true) {
    // Try to find an outgoing halfedge from the current link end that's a silhouette edge
    Mesh::VertexOHalfedgeCCWIter voh_it, voh_it_end = mesh.voh_ccwend(linkEnd);
    for (voh_it = mesh.voh_ccwbegin(linkEnd); voh_it != voh_it_end; ++voh_it) {
      Mesh::HalfedgeHandle heh = *voh_it;
      Mesh::EdgeHandle eh = mesh.edge_handle(heh);
      if (mesh.property(edgeVisited, eh)) {
        continue;
      }
      mesh.property(edgeVisited, eh) = true;

      if (isFeatureEdge(mesh, eh, camPos)) {
        linkEnd = mesh.to_vertex_handle(heh);
        (linkVertices->*dequePushFunc)(mesh.point(linkEnd));
        break;
      }
    }
    if (voh_it == voh_it_end) {
      // No such halfedge was found; the link end cannot be extended further.
      break;
    }
  }
  return linkEnd;
}

static void generateFeatureLinks(Mesh& mesh, const Vec3f& camPos, EPropHandleT<bool>& edgeVisited,
                                 vector<Link>* featureLinks) {
  featureLinks->clear();

  // clear all edgeVisited flags
  Mesh::EdgeIter e_it, e_it_end = mesh.edges_end();
  for (e_it = mesh.edges_begin(); e_it != e_it_end; ++e_it) {
    mesh.property(edgeVisited, *e_it) = false;
  }

  // Visit all edges and group them into visible links.
  for (e_it = mesh.edges_begin(); e_it != e_it_end; ++e_it) {
    Mesh::EdgeHandle eh = *e_it;
    if (mesh.property(edgeVisited, eh)) {
      continue;
    }
    mesh.property(edgeVisited, eh) = true;

    if (!isFeatureEdge(mesh, eh, camPos)) {
      continue;
    }

    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
    Mesh::VertexHandle s_vh = mesh.from_vertex_handle(heh);
    Mesh::VertexHandle t_vh = mesh.to_vertex_handle(heh);
    Vec3f s = mesh.point(s_vh);
    Vec3f t = mesh.point(t_vh);

    featureLinks->push_back(Link());
    Link& newLink = featureLinks->back();
    newLink.vertices.push_back(s);
    newLink.vertices.push_back(t);

    // extend the two ends of the link as much as possible
    extendFeatureLinkEnd(mesh, t_vh, &newLink.vertices, &deque<Vec3f>::push_back, camPos, edgeVisited);
    extendFeatureLinkEnd(mesh, s_vh, &newLink.vertices, &deque<Vec3f>::push_front, camPos, edgeVisited);
  }
}

static Mesh::HalfedgeHandle extendContourLinkEnd(Mesh& mesh, Vec3f linkEnd, Mesh::HalfedgeHandle linkEnd_heh,
                                                 const Vec3f& camPos, float nDotViewMax, float DwkrMin,
                                                 deque<Vec3f>* linkVertices, void(deque<Vec3f>::*dequePushFunc)(const Vec3f&),
                                                 FPropHandleT<bool>& faceVisited,
                                                 VPropHandleT<double>& viewCurvature, FPropHandleT<Vec3f>& viewCurvatureDerivative) {
  // extend the t end of the link as much as possible
  while (true) {
    Vec3f newLinkEnd;
    Mesh::HalfedgeHandle newLinkEnd_heh;
    // Find a line of kw=0 on this face, if possible
    Mesh::HalfedgeHandle linkEnd_opp_heh = mesh.opposite_halfedge_handle(linkEnd_heh);
    Mesh::FaceHandle linkEnd_opp_fh = mesh.face_handle(linkEnd_opp_heh);
    if (mesh.property(faceVisited, linkEnd_opp_fh)) {
      break;
    }
    mesh.property(faceVisited, linkEnd_opp_fh) = true;
    Mesh::HalfedgeHandle linkEnd_opp_next_heh = mesh.next_halfedge_handle(linkEnd_opp_heh);
    Mesh::VertexHandle vh0 = mesh.to_vertex_handle(linkEnd_opp_heh);
    Mesh::VertexHandle vh1 = mesh.to_vertex_handle(linkEnd_opp_next_heh);
    double kw0 = mesh.property(viewCurvature, vh0);
    double kw1 = mesh.property(viewCurvature, vh1);
    if (kw0 > 0.0 != kw1 > 0.0) {
      double a = kw0 / (kw0 - kw1);
      newLinkEnd = (1.0 - a)*mesh.point(vh0) + a*mesh.point(vh1);
      newLinkEnd_heh = linkEnd_opp_next_heh;
    } else {
      Mesh::HalfedgeHandle linkEnd_opp_prev_heh = mesh.prev_halfedge_handle(linkEnd_opp_heh);
      Mesh::VertexHandle vh2 = mesh.to_vertex_handle(linkEnd_opp_prev_heh);
      double kw2 = mesh.property(viewCurvature, vh2);
      if (kw1 > 0.0 != kw2 > 0.0) {
        double a = kw1 / (kw1 - kw2);
        newLinkEnd = (1.0 - a)*mesh.point(vh1) + a*mesh.point(vh2);
        newLinkEnd_heh = linkEnd_opp_prev_heh;
      } else {
        break;
      }
    }
    if (!isContourLineAcceptable(mesh, linkEnd_opp_fh, linkEnd, newLinkEnd, camPos, nDotViewMax, DwkrMin, viewCurvatureDerivative)) {
      break;
    }
    (linkVertices->*dequePushFunc)(newLinkEnd);
    linkEnd = newLinkEnd;
    linkEnd_heh = newLinkEnd_heh;
  }
  return linkEnd_heh;
}

static void generateContourLinks(Mesh &mesh, const Vec3f& camPos,
                                 float nDotViewMax, float DwkrMin,
                                 vector<Link>* contourLinks,
                                 FPropHandleT<bool>& faceVisited,
                                 VPropHandleT<double>& viewCurvature, FPropHandleT<Vec3f>& viewCurvatureDerivative) {
  contourLinks->clear();

  // clear all faceVisited flags
  Mesh::FaceIter f_it, f_it_end = mesh.faces_end();
  for (f_it = mesh.faces_begin(); f_it != f_it_end; ++f_it) {
    mesh.property(faceVisited, *f_it) = false;
  }

  for (f_it = mesh.faces_begin(); f_it != f_it_end; ++f_it) {
    Mesh::FaceHandle fh = *f_it;
    if (mesh.property(faceVisited, fh)) {
      continue;
    }
    mesh.property(faceVisited, fh) = true;

    Vec3f s, t;
    Mesh::HalfedgeHandle s_heh, t_heh;
    if (!isSuggestiveContourFace(mesh, fh, camPos, viewCurvature, viewCurvatureDerivative, nDotViewMax, DwkrMin,
                                 &s, &t, &s_heh, &t_heh)) {
      continue;
    }

    contourLinks->push_back(Link());
    Link& newLink = contourLinks->back();
    newLink.vertices.push_back(s);
    newLink.vertices.push_back(t);

    // extend the ends of the link as much as possible
    extendContourLinkEnd(mesh, t, t_heh, camPos, nDotViewMax, DwkrMin, &newLink.vertices, &deque<Vec3f>::push_back,
                         faceVisited, viewCurvature, viewCurvatureDerivative);
    extendContourLinkEnd(mesh, s, s_heh, camPos, nDotViewMax, DwkrMin, &newLink.vertices, &deque<Vec3f>::push_front,
                         faceVisited, viewCurvature, viewCurvatureDerivative);
  }
}

static void generateVisibleLinks(const vector<Link> links,
                                 const GLfloat* depthBuffer, int width, int height, float epsilon,
                                 const Vec3f& camPos, const Vec3f& camLookDir,
                                 vector<Link>* visibleLinks) {
  visibleLinks->clear();

  for (const Link& link : links) {

    Vec3f start = link.vertices.front();
    bool startVisible = isVisible(toImagePlane(start), &depthBuffer[0], width, height, epsilon);
    int firstLinkAt = -1;
    if (startVisible) {
      firstLinkAt = visibleLinks->size();
      visibleLinks->push_back(Link());
      visibleLinks->back().vertices.push_back(start);
    }

    Vec3f prevSamplePos = start;
    bool prevSamplePosVisible = startVisible;

    // Walk through link edges, recording visible portions as we encounter them
    for (int tIndex = 1; tIndex < link.vertices.size(); ++tIndex) {
      Vec3f s = link.vertices[tIndex - 1];
      Vec3f t = link.vertices[tIndex];
      float s_depth, t_depth;
      Vec2f s_proj = toImagePlane(s, &s_depth);
      Vec2f t_proj = toImagePlane(t, &t_depth);
      float w_s = (s - camPos) | camLookDir;
      float w_t = (t - camPos) | camLookDir;

      // split edge into segments of equal length in screen-space.
      const float MAX_EDGE_SEG_LENGTH = 4.f;  // measured in pixels in screenspace
      int numSegments = 2;            //ceilf((t_proj - s_proj).length() / MAX_EDGE_SEG_LENGTH);
      bool linkEndsOnThisEdge = false;
      for (int i = 1; i <= numSegments; i++) {
        float b = i / (float)(numSegments);
        float a = b*w_s / ((1.f - b)*w_t + b*w_s);
        Vec3f samplePos = (1.f - a)*s + a*t;
        bool samplePosVisible = isVisible(toImagePlane(samplePos), &depthBuffer[0], width, height, epsilon);

        /*float unused;               // debug code to check screen-space lengths of segments
        Vec2f vertexProj = toImagePlane(vertex, &unused);
        Vec2f prevVertexProj = toImagePlane(prevVertex, &unused);
        float d = (vertexProj - prevVertexProj).norm();*/

        if (prevSamplePosVisible) {
          Vec3f newEnd;
          if (samplePosVisible) {
            newEnd = samplePos;
          } else {
            newEnd = findVisibilityBoundary(prevSamplePos, samplePos, &depthBuffer[0], width, height, epsilon);
          }
          // Extend the current link
          if (linkEndsOnThisEdge) {
            visibleLinks->back().vertices.back() = newEnd;
          } else {
            visibleLinks->back().vertices.push_back(newEnd);
            linkEndsOnThisEdge = true;
          }
          /*visibleLinks->push_back(Link());       // debug code for visualizing the segmented edges
          visibleLinks->back().vertices.push_back(prevSamplePos);
          visibleLinks->back().vertices.push_back(newEnd);
          linkEndsOnThisEdge = true;*/
        } else {
          if (samplePosVisible) {
            // Start a new link
            visibleLinks->push_back(Link());
            Vec3f newStart = findVisibilityBoundary(samplePos, prevSamplePos, &depthBuffer[0], width, height, epsilon);
            visibleLinks->back().vertices.push_back(newStart);
            visibleLinks->back().vertices.push_back(samplePos);
            linkEndsOnThisEdge = true;
          }
        }
        prevSamplePos = samplePos;
        prevSamplePosVisible = samplePosVisible;
      }
    }

    if (start == link.vertices.back() && startVisible) {
      // The silhouette link is a loop; the first and last visible links need to be joined into one link
      // (unless the whole loop is visible, in which case the first and last links are already the same link).
      // Append the edges in the last visible link to the front of the edges in the first visible link,
      // then remove the last visible link from visibleLinks->
      assert(0 <= firstLinkAt && firstLinkAt < visibleLinks->size());
      if (firstLinkAt != visibleLinks->size() - 1) {
        Link& firstLink = visibleLinks->at(firstLinkAt);
        Link& lastLink = visibleLinks->back();
        deque<Vec3f>& firstLinkVertices = firstLink.vertices;
        deque<Vec3f>& lastLinkVertices = lastLink.vertices;
        firstLinkVertices.insert(firstLinkVertices.begin(), lastLinkVertices.begin(), lastLinkVertices.end() - 1);
        visibleLinks->pop_back();
      }
    }
  }

  // sanity check!!!
  for (Link& link : *visibleLinks) {
    for (int i = 1; i < link.vertices.size(); ++i) {
      assert(link.vertices[i - 1] != link.vertices[i]);
    }
  }
}


static void writeCatmullRomSpline(ofstream& outfile, const vector<Vec2f>& points) {
  if (points.size() == 2) {
    outfile << "\t<line x1=\"" << points[0][0] << "\" y1=\"" << points[0][1]
            << "\" x2=\"" << points[1][0] << "\" y2=\"" << points[1][1] << "\" stroke-width=\"1\" />" << endl;
  } else if (points.size() > 2) {
    Vec2f mPrev;
    for (int i = 1; i < points.size() - 1; ++i) {
      const Vec2f& pPrev = points[i - 1];
      const Vec2f& p = points[i];
      const Vec2f& pNext = points[i + 1];
      Vec2f m = (pNext - pPrev) / 2.f;//(sqrtf((pNext - p).norm()) + sqrtf((p - pPrev).norm()));  // centripetal (alpha = 0.5)
      if (i == 1) {
        // Quadratic curve to start: Bezier ctrl pts are pPrev, p - m/2, p
        Vec2f q = p - m / 2.f;
        outfile << "\t<path d=\"M " << pPrev[0] << " " << pPrev[1]
          << " Q " << q[0] << " " << q[1] << ", " << p[0] << " " << p[1] << "\" stroke-width=\"1\" />" << endl;
      } else {
        // Bezier ctrl pts are pPrev, pPrev + mPrev/3, p - m/3, p
        Vec2f c2 = p - m / 3.f;
        if (i == 2) {
          Vec2f c1 = pPrev + mPrev / 3.f;
          outfile << "\t<path d=\"M " << pPrev[0] << " " << pPrev[1];
          outfile << endl << "\t\tC " << c1[0] << " " << c1[1] << ", " << c2[0] << " " << c2[1] << ", " << p[0] << " " << p[1];
        } else {
          outfile << endl << "\t\tS " << c2[0] << " " << c2[1] << ", " << p[0] << " " << p[1];
        }
      }
      mPrev = m;
    }
    if (points.size() > 3) {
      outfile << "\"" << endl << "\tstroke-width=\"1\" />" << endl;
    }
    // Quadratic curve to end: Bezier ctrl pts are pPrev, pPrev + mPrev/2, p
    const Vec2f& pPrev = points[points.size() - 2];
    const Vec2f& p = points[points.size() - 1];
    const Vec2f& q = pPrev + mPrev / 2.f;
    outfile << "\t<path d=\"M " << pPrev[0] << " " << pPrev[1]
      << " Q " << q[0] << " " << q[1] << ", " << p[0] << " " << p[1] << "\" stroke-width=\"1\" />" << endl;
  }
}

void writeImage(Mesh &mesh, int width, int height, const string& filename,
                const Vec3f& camPos, const Vec3f& camLookDir,
                float nDotViewMax, float DwkrMin, float epsilon,
                EPropHandleT<bool>& edgeVisited, FPropHandleT<bool>& faceVisited,
                VPropHandleT<double>& viewCurvature, FPropHandleT<Vec3f>& viewCurvatureDerivative) {

  // copy entire depth buffer
  vector<GLfloat> depthBuffer(width * height);
  glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, &depthBuffer[0]);

  vector<Link> featureLinks;
  generateFeatureLinks(mesh, camPos, edgeVisited, &featureLinks);
  vector<Link> visibleFeatureLinks;
  generateVisibleLinks(featureLinks, &depthBuffer[0], width, height, epsilon, camPos, camLookDir, &visibleFeatureLinks);
  
  vector<Link> contourLinks;
  generateContourLinks(mesh, camPos, nDotViewMax, DwkrMin, &contourLinks,
                       faceVisited, viewCurvature, viewCurvatureDerivative);
  vector<Link> visibleContourLinks;
  generateVisibleLinks(contourLinks, &depthBuffer[0], width, height, epsilon, camPos, camLookDir, &visibleContourLinks);
  



  ofstream outfile(filename.c_str());
  outfile << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
  outfile << "<svg width=\"5in\" height=\"5in\" viewBox=\"0 0 " << width << ' ' << height << "\">" << endl;

  char* strokeStrings[8] = { "<g stroke=\"black\" fill=\"none\">",
                              "<g stroke=\"red\" fill=\"none\">",
                              "<g stroke=\"blue\" fill=\"none\">",
                              "<g stroke=\"magenta\" fill=\"none\">",
                              "<g stroke=\"orange\" fill=\"none\">",
                              "<g stroke=\"deeppink\" fill=\"none\">",
                              "<g stroke=\"olivedrab\" fill=\"none\">",
                              "<g stroke=\"deepskyblue\" fill=\"none\">" };


  vector<Vec2f> imagePoints;
  for (int i = 0; i < visibleFeatureLinks.size(); i++) {
    Link& link = visibleFeatureLinks[i];

    outfile << strokeStrings[i % 8] << endl;

    toImagePlaneInvertHeight(link.vertices, height, &imagePoints);
    writeCatmullRomSpline(outfile, imagePoints);

    /*for (int i = 1; i < link.vertices.size(); i++) {
      Vec3f source = link.vertices[i - 1];
      Vec3f dest = link.vertices[i];
      Vec3f p1 = toImagePlane(source);
      Vec3f p2 = toImagePlane(dest);
      outfile << "\t<line ";
      outfile << "x1=\"" << p1[0] << "\" ";
      outfile << "y1=\"" << height - p1[1] << "\" ";
      outfile << "x2=\"" << p2[0] << "\" ";
      outfile << "y2=\"" << height - p2[1] << "\" stroke-width=\"1\" />" << endl;
    }*/
    outfile << "</g>" << endl;
  }


  /*for (int i = 0; i < visibleFeatureLinks.size(); i++) {
    Link& link = visibleFeatureLinks[i];
    outfile << strokeStrings[0] << endl;
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
  }*/
  
  outfile << "</svg>" << endl;
  outfile.close();

  std::cout << filename << " written" << std::endl;
}
