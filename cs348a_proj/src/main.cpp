#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <GL/glut.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
#include "mesh_features.h"
#include "image_generation.h"
#include "decimate.h"

using namespace std;
using namespace OpenMesh;
using namespace Eigen;

VPropHandleT<double> viewCurvature;
FPropHandleT<Vec3f> viewCurvatureDerivative;
VPropHandleT<CurvatureInfo> curvature;
EPropHandleT<bool> edgeVisited;
FPropHandleT<bool> faceVisited;
Mesh mesh;


bool leftDown = false, rightDown = false, middleDown = false;
int lastPos[2];
float cameraPos[4] = { 0, 0, 4, 1 };
Vec3f up, pan;
int windowWidth = 640, windowHeight = 480;
bool showSurface = true, showAxes = true, showCurvature = false, showNormals = false, showEdges = false;
float nDotViewMax = 0.6f, DwkrMin = 1.f;

float specular[] = { 1.0, 1.0, 1.0, 1.0 };
float shininess[] = { 50.0 };

void renderSuggestiveContours(Vec3f actualCamPos) { // use this camera position to account for panning etc.
  glColor3f(0.5, 0.5, 0.5);

  // RENDER SUGGESTIVE CONTOURS HERE -----------------------------------------------------------------------------
  //glColor3f(1, 0, 0);
  glBegin(GL_LINES);
  Mesh::ConstFaceIter f_it, f_it_end = mesh.faces_end();
  for (f_it = mesh.faces_begin(); f_it != f_it_end; ++f_it) {
    Vec3f p1, p2;
    Mesh::HalfedgeHandle unused;
    if (isSuggestiveContourFace(mesh, *f_it, actualCamPos, viewCurvature, viewCurvatureDerivative, nDotViewMax, DwkrMin, &p1, &p2, &unused, &unused)) {
      glVertex3f(p1[0], p1[1], p1[2]);
      glVertex3f(p2[0], p2[1], p2[2]);
    }
  }
  glEnd();
  // -------------------------------------------------------------------------------------------------------------
}

void renderMesh() {
  if (!showSurface) glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // render regardless to remove hidden lines

  glEnable(GL_LIGHTING);
  glLightfv(GL_LIGHT0, GL_POSITION, cameraPos);

  glDepthRange(0.001, 1);
  glEnable(GL_NORMALIZE);

  // WRITE CODE HERE TO RENDER THE TRIANGLES OF THE MESH ---------------------------------------------------------
  glBegin(GL_TRIANGLES);
  Mesh::ConstFaceIter f_it, f_it_end = mesh.faces_end();
  for (f_it = mesh.faces_begin(); f_it != f_it_end; ++f_it) {
    Mesh::ConstFaceVertexIter cfv_it = mesh.cfv_iter(*f_it);
    Vec3f points[3];
    Vec3f normals[3];
    for (int i = 0; i < 3; i++) {
      points[i] = mesh.point(*cfv_it);
      normals[i] = mesh.normal(*cfv_it);
      ++cfv_it;
    }
    glNormal3f(normals[0][0], normals[0][1], normals[0][2]);
    glVertex3f(points[0][0], points[0][1], points[0][2]);
    glNormal3f(normals[1][0], normals[1][1], normals[1][2]);
    glVertex3f(points[1][0], points[1][1], points[1][2]);
    glNormal3f(normals[2][0], normals[2][1], normals[2][2]);
    glVertex3f(points[2][0], points[2][1], points[2][2]);
  }
  glEnd();  
  // -------------------------------------------------------------------------------------------------------------

  if (!showSurface) glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

  glDisable(GL_LIGHTING);
  glDepthRange(0, 0.999);

  Vec3f actualCamPos(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2]);
  renderSuggestiveContours(actualCamPos);

  // We'll be nice and provide you with code to render feature edges below
  glBegin(GL_LINES);
  glColor3f(0, 0, 0);
  glLineWidth(2.0f);
  for (Mesh::ConstEdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
    const Mesh::EdgeHandle eh = (*e_it);
    if (isFeatureEdge(mesh, eh, actualCamPos)) {
      const Mesh::HalfedgeHandle heh_0 = mesh.halfedge_handle(eh, 0);
      const Mesh::HalfedgeHandle heh_1 = mesh.halfedge_handle(eh, 1);
      const Vec3f source(mesh.point(mesh.from_vertex_handle(heh_0)));
      const Vec3f target(mesh.point(mesh.from_vertex_handle(heh_1)));
      glVertex3f(source[0], source[1], source[2]);
      glVertex3f(target[0], target[1], target[2]);
    }
  }
  glEnd();

  if (showCurvature) {
    // WRITE CODE HERE TO RENDER THE PRINCIPAL DIRECTIONS YOU COMPUTED ---------------------------------------------
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    Mesh::VertexIter v_it, v_it_end = mesh.vertices_end();
    for (v_it = mesh.vertices_begin(); v_it != v_it_end; ++v_it) {
      CurvatureInfo info;
      info = mesh.property(curvature, *v_it);
      const Vec3f p = mesh.point(*v_it);
      const Vec3f pT1 = p + 0.02 * info.directions[0];
      const Vec3f pT2 = p + 0.02 * info.directions[1];
      glVertex3f(p[0], p[1], p[2]);
      glVertex3f(pT1[0], pT1[1], pT1[2]);
      glVertex3f(p[0], p[1], p[2]);
      glVertex3f(pT2[0], pT2[1], pT2[2]);
    }
    glEnd();
    // -------------------------------------------------------------------------------------------------------------
  }

  if (showNormals) {
    glBegin(GL_LINES);
    glColor3f(0, 1, 0);
    for (Mesh::ConstVertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
      const Mesh::VertexHandle vh = (*v_it);
      const Vec3f n = mesh.normal(vh);
      const Vec3f p = mesh.point(vh);
      const Vec3f d = p + n * 0.01;
      glVertex3f(p[0], p[1], p[2]);
      glVertex3f(d[0], d[1], d[2]);
    }
    glEnd();
  }

  if (showEdges) {
    glBegin(GL_LINES);
    glColor3f(0, 0, 0);
    Mesh::ConstEdgeIter e_it, e_it_end = mesh.edges_end();
    for (e_it = mesh.edges_begin(); e_it != e_it_end; ++e_it) {
      Mesh::HalfedgeHandle heh = mesh.halfedge_handle(*e_it, 0);
      Mesh::VertexHandle vh_s = mesh.from_vertex_handle(heh);
      Mesh::VertexHandle vh_t = mesh.to_vertex_handle(heh);
      Vec3f s = mesh.point(vh_s);
      Vec3f t = mesh.point(vh_t);
      Vec3f sn = mesh.normal(vh_s);
      Vec3f tn = mesh.normal(vh_t);
      s += 0.001f * sn;
      t += 0.001f * tn;
      glVertex3f(s[0], s[1], s[2]);
      glVertex3f(t[0], t[1], t[2]);
    }
    glEnd();
  }

  glDepthRange(0, 1);
}

void display() {
  glClearColor(1, 1, 1, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glShadeModel(GL_SMOOTH);
  glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
  glEnable(GL_LIGHT0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0, 0, windowWidth, windowHeight);

  float ratio = (float)windowWidth / (float)windowHeight;
  gluPerspective(50, ratio, 1, 1000); // 50 degree vertical viewing angle, zNear = 1, zFar = 1000

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2],
    pan[0], pan[1], pan[2], up[0], up[1], up[2]);

  // Draw mesh
  renderMesh();

  // Draw axes
  if (showAxes) {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    glLineWidth(1);
    glColor3f(1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(1, 0, 0); // x axis
    glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, 1, 0); // y axis
    glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, 1); // z axis
    glEnd(/*GL_LINES*/);
  }

  glutSwapBuffers();
}

void mouse(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) leftDown = (state == GLUT_DOWN);
  else if (button == GLUT_RIGHT_BUTTON) rightDown = (state == GLUT_DOWN);
  else if (button == GLUT_MIDDLE_BUTTON) middleDown = (state == GLUT_DOWN);

  lastPos[0] = x;
  lastPos[1] = y;
}

void mouseMoved(int x, int y) {
  const float speed = 30.0f;

  int dx = x - lastPos[0];
  int dy = y - lastPos[1];

  Vec3f curCamera(cameraPos[0], cameraPos[1], cameraPos[2]);
  Vec3f curCameraNormalized = curCamera.normalized();
  Vec3f right = up % curCameraNormalized;

  if (middleDown || (leftDown && rightDown)) {
    pan += -speed * (float)((float)dx / (float)windowWidth) * right +
      speed * (float)((float)dy / (float)windowHeight) * up;
  }
  else if (leftDown) {
    // Assume here that up vector is (0,1,0)
    Vec3f newPos = curCamera - speed * (float)((float)dx / (float)windowWidth) * right +
      speed * (float)((float)dy / (float)windowHeight) * up;
    newPos = newPos.normalized() * curCamera.length();

    up = up - (up | newPos) * newPos / newPos.sqrnorm();
    up.normalize();

    for (int i = 0; i < 3; i++) cameraPos[i] = newPos[i];
  }
  else if (rightDown) {
    for (int i = 0; i < 3; i++) cameraPos[i] *= pow(1.1, dy * 0.1);
  }


  lastPos[0] = x;
  lastPos[1] = y;

  Vec3f actualCamPos(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2]);
  computeViewCurvature(mesh, actualCamPos, curvature, viewCurvature, viewCurvatureDerivative);

  glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
  Vec3f actualCamPos(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2]);

  if (key == 's' || key == 'S') showSurface = !showSurface;
  else if (key == 'a' || key == 'A') showAxes = !showAxes;
  else if (key == 'c' || key == 'C') showCurvature = !showCurvature;
  else if (key == 'n' || key == 'N') showNormals = !showNormals;
  else if (key == 'e' || key == 'E') showEdges = !showEdges;
  else if (key == '1') {
    nDotViewMax -= 0.01f;
    printf("nDotViewMax = %f\n", nDotViewMax);
  } else if (key == '2') {
    nDotViewMax += 0.01f;
    printf("nDotViewMax = %f\n", nDotViewMax);
  } else if (key == '3') {
    DwkrMin /= 1.05f;
    printf("DwkrMin = %f\n", DwkrMin);
  } else if (key == '4') {
    DwkrMin *= 1.05f;
    printf("DwkrMin = %f\n", DwkrMin);
  }
  else if (key == 'd' || key == 'D') {
    float percentage = 1.0f;
    while (percentage <= 0.0f || percentage >= 1.0f) {
      std::cout << "Type percentage of vertices: ";
      std::cin >> percentage;
    }
    simplify(mesh, percentage, "output.off");
  }
  else if (key == 'w' || key == 'W') {
    Vec3f cameraLookDir(-cameraPos[0], -cameraPos[1], -cameraPos[2]);
    cameraLookDir.normalize();
    writeImage(mesh, windowWidth, windowHeight, "renderedImage.svg", actualCamPos, cameraLookDir,
               nDotViewMax, DwkrMin, edgeVisited, faceVisited, viewCurvature, viewCurvatureDerivative);
  }
  else if (key == 'q' || key == 'Q') exit(0);
  glutPostRedisplay();
}

void reshape(int width, int height) {
  windowWidth = width;
  windowHeight = height;
  glutPostRedisplay();
}

int main(int argc, char** argv) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " mesh_filename\n";
    exit(0);
  }

  IO::Options opt;
  opt += IO::Options::VertexNormal;
  opt += IO::Options::FaceNormal;

  mesh.request_face_normals();
  mesh.request_vertex_normals();

  cout << "Reading from file " << argv[1] << "...\n";
  if (!IO::read_mesh(mesh, argv[1], opt)) {
    cout << "Read failed.\n";
    exit(0);
  }

  cout << "Mesh stats:\n";
  cout << '\t' << mesh.n_vertices() << " vertices.\n";
  cout << '\t' << mesh.n_edges() << " edges.\n";
  cout << '\t' << mesh.n_faces() << " faces.\n";

  mesh.update_normals();

  mesh.add_property(viewCurvature);
  mesh.add_property(viewCurvatureDerivative);
  mesh.add_property(curvature);
  mesh.add_property(edgeVisited);

  // Move center of mass to origin
  Vec3f center(0, 0, 0);
  for (Mesh::ConstVertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    center += mesh.point(*v_it);
  center /= mesh.n_vertices();

  for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    mesh.point(*v_it) -= center;

  // Fit in the unit sphere
  float maxLength = 0;
  for (Mesh::ConstVertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    maxLength = max(maxLength, mesh.point(*v_it).length());

  if (maxLength > 0) {
    for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
      mesh.point(*v_it) /= maxLength;
  }
                                                            simplify(mesh, 0.5, "output.off");    // TESTOESTESTSTS!!!!!!!!
  computeCurvature(mesh, curvature);

  up = Vec3f(0, 1, 0);
  pan = Vec3f(0, 0, 0);

  Vec3f actualCamPos(cameraPos[0] + pan[0], cameraPos[1] + pan[1], cameraPos[2] + pan[2]);
  computeViewCurvature(mesh, actualCamPos, curvature, viewCurvature, viewCurvatureDerivative);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(windowWidth, windowHeight);
  glutCreateWindow(argv[0]);

  glutDisplayFunc(display);
  glutMotionFunc(mouseMoved);
  glutMouseFunc(mouse);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);

  GLint depth;
  glGetIntegerv(GL_DEPTH_BITS, &depth);
  printf("%d bits in depth buffer\n", depth);

  glutMainLoop();

  return 0;
}
