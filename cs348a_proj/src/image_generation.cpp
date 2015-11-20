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

bool isVisible(const Vec3f& point, const GLfloat* depthBuffer, int width, int height) {
  Vec3f proj = toImagePlane(point);
  if (!(0 <= proj[0] && proj[0] < width && 0 <= proj[1] && proj[1] < height)) return false;
  return (depthBuffer[((int)proj[1]) * width + (int)proj[0]] - proj[2]) > -EPSILON;
}

bool isVisible(const Vec3f& point, const GLfloat* depthBuffer, int width, int height, Vec3f* projected) {
  Vec3f proj = toImagePlane(point);
  *projected = proj;
  if (!(0 <= proj[0] && proj[0] < width && 0 <= proj[1] && proj[1] < height)) return false;
  return (depthBuffer[((int)proj[1]) * width + (int)proj[0]] - proj[2]) > -EPSILON;
}

// Returns true if edge is wholly or partially visible. If only partially visible, the occluded endpoint
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
}


void writeImage(Mesh &mesh, int width, int height, string filename, Vec3f camPos) {

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
  //           adj                               adj
  //
  // One end occluded:
  //
  // start  o======>o=====>o=====>o....o=====>o===|-->o  end
  //        ^                                     ^
  //       adj                                   adj
  //
  // The LinkEndpoint struct describes one end of a link: vh points to the mesh vertex that this end
  // is attached to. If the vertex is not visible, then adj is the position on the link edge that's on
  // the boundary between the visible and occluded portions of the edge. If vertex is visible, then adj
  // is simply the position of the mesh vertex.
  // The VisibleLink struct describes a visible link. It's considered complete when no more edges can
  // be added to either end, which can occur in two ways: the link forms a loop, or both ends of the link
  // are occluded.

  struct LinkEndpoint {
    Mesh::VertexHandle vh;
    bool isVisible;
    Vec3f adj;
  };
  struct VisibleLink {
    bool checkIsComplete() {
      isComplete = (start.vh == end.vh || (!start.isVisible && !end.isVisible));
      return isComplete;
    }
    deque<Mesh::HalfedgeHandle> halfEdges;  // halfedges oriented and ordered from start to end
    LinkEndpoint start, end;
    bool isComplete;
  };

  vector<VisibleLink> silhouetteLinks;

  // For each edge, check if it's connected to any existing silhouette link(s) at either endpoint.  If so, add the edge
  // to that link.  If the edge is connected to two different links, the two links need to be merged.
  Mesh::EdgeIter e_it, e_it_end = mesh.edges_end();
  for (e_it = mesh.edges_begin(); e_it != e_it_end; ++e_it) {
    Mesh::EdgeHandle eh = *e_it;
    if (!isSilhouette(mesh, eh, camPos)) {
      continue;
    }
    Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);

    LinkEndpoint s, t;
    s.vh = mesh.from_vertex_handle(heh);
    t.vh = mesh.to_vertex_handle(heh);

    if (!isEdgePartialVisible(mesh.point(s.vh), mesh.point(t.vh), &depthBuffer[0], width, height, &s.isVisible, &t.isVisible, &s.adj, &t.adj)) {
      continue;
    }

    VisibleLink* firstAttachedLink = NULL;    // This will be updated to the first link the edge is attached to (if any)
    LinkEndpoint unattachedEndpoint;      // This will be updated to be the unattached end of the edge after the edge is first attached to a link (if any)
    bool unattachedIsStart = false;       // This will be updated at the same time as unattachedVertex; will indicate if unattachedVertex is the start or end of the link the edge attached to.
    
    // Iterate over the existing links until we find one this edge connects to, or we run out of links.
    auto links_it = silhouetteLinks.begin();
    while(links_it != silhouetteLinks.end()) {
      VisibleLink& link = *links_it;
      if (link.isComplete) {
        ++links_it;
        continue;
      }
      // Check if edge vertex s connects to any links
      if (s.isVisible) {
        if (s.vh == link.start.vh) {
          link.halfEdges.push_front(mesh.opposite_halfedge_handle(heh));
          link.start = t;
          unattachedEndpoint = t;
          unattachedIsStart = true;
          firstAttachedLink = &link;
          ++links_it;
          break;
        } else if (s.vh == link.end.vh) {
          link.halfEdges.push_back(heh);
          link.end = t;
          unattachedEndpoint = t;
          unattachedIsStart = false;
          firstAttachedLink = &link;
          ++links_it;
          break;
        }
      }
      // Check if edge vertex t connects to any links
      if (t.isVisible) {
        if (t.vh == link.start.vh) {
          link.halfEdges.push_front(heh);
          link.start = s;
          unattachedEndpoint = s;
          unattachedIsStart = true;
          firstAttachedLink = &link;
          ++links_it;
          break;
        } else if (t.vh == link.end.vh) {
          link.halfEdges.push_back(mesh.opposite_halfedge_handle(heh));
          link.end = s;
          unattachedEndpoint = s;
          unattachedIsStart = false;
          firstAttachedLink = &link;
          ++links_it;
          break;
        }
      }
      // It's possible that s=start and t=end, or s=end and t=start, which menas this edge made a link
      // into a loop. In those cases, the edge will be attached at s, and no further action is required.
      ++links_it;
    }

    
    if (firstAttachedLink != NULL) {
      // If this edge was attached to a link and turned it into a loop, then that link is complete.
      // Additionally, if both ends of this link are occluded, then no more links can attach at either end, so it's complete.
      // We don't need to check this edge against the remaining links; move on to next edge.
      if (firstAttachedLink->checkIsComplete()) {
        continue;
      }
      // If the unattached end of the edge is not visible, there's no need to check if it connects to any other links since nothing
      // further cound be connected to the edge on that end.
      if (!unattachedEndpoint.isVisible) {
        continue;
      }
      
      // At this point, if the edge has been attached to a link, we'll check if the unattached end of the
      // edge connects to any of the remaining links. If so, the two links connected by the edge need to be merged.
      bool linksMerged = false;
      while (links_it != silhouetteLinks.end()) {
        VisibleLink& link = *links_it;
        if (link.isComplete) {
          ++links_it;
          continue;
        }
        if (unattachedEndpoint.vh == link.start.vh) {
          if (unattachedIsStart) {
            // firstAttachedLink <-<-<- unattachedVertex ->->->-> link
            for (int i = 0; i < link.halfEdges.size(); i++) {
              Mesh::HalfedgeHandle& link_heh = link.halfEdges[i];
              firstAttachedLink->halfEdges.push_front(mesh.opposite_halfedge_handle(link_heh));
            }
            firstAttachedLink->start = link.end;
            links_it = silhouetteLinks.erase(links_it);
            linksMerged = true;
            break;
          } else {
            // firstAttachedLink ->->-> unattachedVertex ->->->-> link
            for (int i = 0; i < link.halfEdges.size(); i++) {
              Mesh::HalfedgeHandle& link_heh = link.halfEdges[i];
              firstAttachedLink->halfEdges.push_back(link_heh);
            }
            firstAttachedLink->end = link.end;
            links_it = silhouetteLinks.erase(links_it);
            linksMerged = true;
            break;
          }
        } else if (unattachedEndpoint.vh == link.end.vh) {
          if (unattachedIsStart) {
            // firstAttachedLink <-<-<- unattachedVertex <-<-<-<- link
            for (int i = link.halfEdges.size() - 1; i >= 0; i--) {
              Mesh::HalfedgeHandle& link_heh = link.halfEdges[i];
              firstAttachedLink->halfEdges.push_front(link_heh);
            }
            firstAttachedLink->start = link.start;
            links_it = silhouetteLinks.erase(links_it);
            linksMerged = true;
            break;
          } else {
            // firstAttachedLink ->->-> unattachedVertex <-<-<-<- link
            for (int i = link.halfEdges.size() - 1; i >= 0; i--) {
              Mesh::HalfedgeHandle& link_heh = link.halfEdges[i];
              firstAttachedLink->halfEdges.push_back(mesh.opposite_halfedge_handle(link_heh));
            }
            firstAttachedLink->end = link.start;
            links_it = silhouetteLinks.erase(links_it);
            linksMerged = true;
            break;
          }
        }
        ++links_it;
      }
      if (linksMerged) {
        if (firstAttachedLink->checkIsComplete()) {
          continue;
        }
      }

    } else {
      // This edge doesn't touch any existing link, create a new link from it
      silhouetteLinks.push_back(VisibleLink());
      VisibleLink& newLink = silhouetteLinks.back();
      newLink.halfEdges.push_back(heh);
      newLink.start = s;
      newLink.end = t;
      newLink.isComplete = false;   // at least one of s,t is visible, so this can't be true initially
      //firstAttachedLink = &newLink;
    }
  }


  // sanity check!!!
  for (VisibleLink& link : silhouetteLinks) {
    assert(link.isComplete);
    if (link.start.isVisible) assert(link.start.adj == mesh.point(link.start.vh));
    if (link.end.isVisible) assert(link.end.adj == mesh.point(link.end.vh));
    Mesh::VertexHandle prevT = link.start.vh;
    for (Mesh::HalfedgeHandle& heh : link.halfEdges) {
      assert(prevT == mesh.from_vertex_handle(heh));
      prevT = mesh.to_vertex_handle(heh);
    }
    assert(prevT == link.end.vh);
  }

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
    VisibleLink& link = silhouetteLinks[i];

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
      Vec3f source = (i == 0 ? link.start.adj : mesh.point(mesh.from_vertex_handle(heh)));
      Vec3f target = (i == link.halfEdges.size() - 1 ? link.end.adj : mesh.point(mesh.to_vertex_handle(heh)));

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

  cout << filename << " written" << endl;
}
