/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "crv.h"
#include "crvAdapt.h"
#include "crvBezier.h"
#include "crvShape.h"
#include "crvSnap.h"

#include <lionPrint.h>
#include <pcu_util.h>
#include <PCU.h>
#include <apfDynamicVector.h>
#include <apfElement.h>

#include <apf.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <pcu_util.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <cassert>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#ifdef HAVE_SIMMETRIX
//#include <gmi_sim.h>
//#include <SimUtil.h>
//#include <MeshSim.h>
//#include <SimModel.h>
#endif

#include <reel.h>
#include <sys/types.h> 
#include <sys/stat.h> 
#include <errno.h> 

static void safe_mkdir(
    const char* path)
{
  mode_t const mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
  {
    reel_fail("Err: could not create directory \"%s\"\n", path);
  }
}

static apf::Vector3 getEdgeCenter(
    apf::Mesh2* m,
    apf::MeshEntity* e)
{
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::EDGE);
  apf::MeshEntity* dv[2];
  m->getDownward(e, 0, dv);
  apf::Vector3 center(0., 0., 0.);
  for (int i = 0; i < 2; i++) {
    apf::Vector3 p;
    m->getPoint(dv[i], 0, p);
    center = center + p;
  }
  center = center * (0.5);
  return center;
}

static apf::Mesh2*  makePoint(
    apf::Mesh2* m,
    apf::MeshEntity* e)
{
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::VERTEX);
  apf::Mesh2* sphMesh = apf::makeEmptyMdsMesh(gmi_load(".null"), 1, false);
  double xrange[2] = {1.e16, -1.e16};
  double yrange[2] = {1.e16, -1.e16};
  double zrange[2] = {1.e16, -1.e16};

  apf::Adjacent adj;
  m->getAdjacent(e, 1, adj);

  for (int i = 0; i < (int)adj.getSize(); i++) {
    apf::Vector3 p = getEdgeCenter(m, adj[i]);
    if (p[0] < xrange[0]) xrange[0] = p[0];
    if (p[0] > xrange[1]) xrange[1] = p[0];

    if (p[1] < yrange[0]) yrange[0] = p[1];
    if (p[1] > yrange[1]) yrange[1] = p[1];

    if (p[2] < zrange[0]) zrange[0] = p[2];
    if (p[2] > zrange[1]) zrange[1] = p[2];
  }

  double minsize = 1.e16;
  if (xrange[1]-xrange[0] < minsize) minsize = xrange[1] - xrange[0];
  if (yrange[1]-yrange[0] < minsize) minsize = yrange[1] - yrange[0];
  if (zrange[1]-zrange[0] < minsize) minsize = zrange[1] - zrange[0];

  double radius = minsize / 20.;
  apf::Vector3 center;
  m->getPoint(e, 0, center);

  int n = 20;
  const double pi = 3.141595;
  const apf::Vector3 param(0., 0., 0.);
  std::vector<apf::MeshEntity*> vs;
  vs.clear();
  for (int i = 0; i < n; i++) {
    apf::Vector3 p(0., 0., 0.);
    p[0] = center[0] + radius * std::cos(2.*i*pi/n);
    p[1] = center[1] + radius * std::sin(2.*i*pi/n);
    p[2] = center[2];
    apf::MeshEntity* newV = sphMesh->createVertex(m->toModel(e), p, param);
    vs.push_back(newV);
  }

  for (int i = 0; i < n; i++) {
    apf::MeshEntity* dv[2];
    dv[0] = vs[i];
    dv[1] = vs[(i+1)%20];
    sphMesh->createEntity(apf::Mesh::EDGE, m->toModel(e), dv);
  }

  vs.clear();
  for (int i = 0; i < n; i++) {
    apf::Vector3 p(0., 0., 0.);
    p[0] = center[0] + radius * std::cos(2.*i*pi/n);
    p[1] = center[1];
    p[2] = center[2] + radius * std::sin(2.*i*pi/n);
    apf::MeshEntity* newV = sphMesh->createVertex(m->toModel(e), p, param);
    vs.push_back(newV);
  }

  for (int i = 0; i < n; i++) {
    apf::MeshEntity* dv[2];
    dv[0] = vs[i];
    dv[1] = vs[(i+1)%20];
    sphMesh->createEntity(apf::Mesh::EDGE, m->toModel(e), dv);
  }

  vs.clear();
  for (int i = 0; i < n; i++) {
    apf::Vector3 p(0., 0., 0.);
    p[0] = center[0];
    p[1] = center[1] + radius * std::cos(2.*i*pi/n);
    p[2] = center[2] + radius * std::sin(2.*i*pi/n);
    apf::MeshEntity* newV = sphMesh->createVertex(m->toModel(e), p, param);
    vs.push_back(newV);
  }

  for (int i = 0; i < n; i++) {
    apf::MeshEntity* dv[2];
    dv[0] = vs[i];
    dv[1] = vs[(i+1)%20];
    sphMesh->createEntity(apf::Mesh::EDGE, m->toModel(e), dv);
  }
  sphMesh->acceptChanges();
  apf::deriveMdsModel(sphMesh);
  sphMesh->verify();
  return sphMesh;
}

static void makeEntMeshes(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Mesh2* &entMeshLinear,
    apf::Mesh2* &entMeshCurved)
{
  PCU_ALWAYS_ASSERT(!entMeshLinear);
  PCU_ALWAYS_ASSERT(!entMeshCurved);

  const apf::Vector3 param(0., 0., 0.);

  if (m->getType(e) == apf::Mesh::VERTEX)
  {
    entMeshLinear = makePoint(m, e);
    entMeshCurved = makePoint(m, e);
    return;
  }
  if (m->getType(e) == apf::Mesh::EDGE)
  {
    entMeshLinear = apf::makeEmptyMdsMesh(gmi_load(".null"), 1, false);
    entMeshCurved = apf::makeEmptyMdsMesh(gmi_load(".null"), 1, false);
    apf::MeshEntity* vs[2];
    m->getDownward(e, 0, vs);
    apf::Vector3 p[2];
    m->getPoint(vs[0], 0, p[0]);
    m->getPoint(vs[1], 0, p[1]);
    apf::MeshEntity* newVs[2];
    newVs[0] = entMeshLinear->createVertex(0, p[0], param);
    newVs[1] = entMeshLinear->createVertex(0, p[1], param);
    entMeshLinear->createEntity(apf::Mesh::EDGE, 0, newVs);

    apf::MeshEntity* newVsc[2];
    newVsc[0] = entMeshCurved->createVertex(0, p[0], param);
    newVsc[1] = entMeshCurved->createVertex(0, p[1], param);
    apf::MeshEntity* edge = entMeshCurved->createEntity(apf::Mesh::EDGE, 0, newVsc);

    entMeshLinear->acceptChanges();
    apf::deriveMdsModel(entMeshLinear);
    entMeshLinear->verify();

    entMeshCurved->acceptChanges();
    apf::deriveMdsModel(entMeshCurved);
    entMeshCurved->verify();

    apf::FieldShape* fs = m->getShape();
    entMeshCurved->changeShape(fs, true);
    if (fs->countNodesOn(apf::Mesh::EDGE))
    {
      for (int i = 0; i < fs->countNodesOn(apf::Mesh::EDGE); i++) {
	apf::Vector3 p;
	m->getPoint(e, i, p);
	entMeshCurved->setPoint(edge, i, p);
      }
    }
    entMeshCurved->acceptChanges();
    return;
  }
  if (m->getType(e) == apf::Mesh::TRIANGLE)
  {
    entMeshLinear = apf::makeEmptyMdsMesh(gmi_load(".null"), 2, false);
    entMeshCurved = apf::makeEmptyMdsMesh(gmi_load(".null"), 2, false);
    apf::MeshEntity* downverts[3];
    apf::MeshEntity* downedges[3];
    m->getDownward(e, 0, downverts);
    m->getDownward(e, 1, downedges);
    int edge_vert[3][2];
    for (int i = 0; i < 3; i++) {
      apf::MeshEntity* dv[2];
      m->getDownward(downedges[i], 0, dv);
      int i0 = apf::findIn(downverts, 3, dv[0]);
      int i1 = apf::findIn(downverts, 3, dv[1]);
      PCU_ALWAYS_ASSERT(i0 != -1);
      PCU_ALWAYS_ASSERT(i1 != -1);
      edge_vert[i][0] = i0;
      edge_vert[i][1] = i1;
    }

    apf::MeshEntity* newvertsLinear[3];
    apf::MeshEntity* newedgesLinear[3];
    apf::MeshEntity* newvertsCurved[3];
    apf::MeshEntity* newedgesCurved[3];
    apf::Vector3 param(0.,0.,0.);
    for (int i = 0; i < 3; i++) {
      apf::Vector3 p;
      m->getPoint(downverts[i], 0, p);
      newvertsLinear[i] = entMeshLinear->createVertex(0, p, param);
      newvertsCurved[i] = entMeshCurved->createVertex(0, p, param);
    }

    for (int i = 0; i < 3; i++) {
      apf::MeshEntity* evLinear[2] = {
      	newvertsLinear[edge_vert[i][0]],
      	newvertsLinear[edge_vert[i][1]]
      };
      newedgesLinear[i] = entMeshLinear->createEntity(
    	apf::Mesh::EDGE, 0, evLinear);

      apf::MeshEntity* evCurved[2] = {
      	newvertsCurved[edge_vert[i][0]],
      	newvertsCurved[edge_vert[i][1]]
      };
      newedgesCurved[i] = entMeshCurved->createEntity(
    	apf::Mesh::EDGE, 0, evCurved);
    }


    entMeshLinear->createEntity(
    	apf::Mesh::TRIANGLE, 0, newedgesLinear);

    apf::MeshEntity* face =
    entMeshCurved->createEntity(
    	apf::Mesh::TRIANGLE, 0, newedgesCurved);


    entMeshLinear->acceptChanges();
    apf::deriveMdsModel(entMeshLinear);
    entMeshLinear->verify();

    entMeshCurved->acceptChanges();
    apf::deriveMdsModel(entMeshCurved);
    entMeshCurved->verify();

    apf::FieldShape* fs = m->getShape();
    entMeshCurved->changeShape(fs, true);

    int nnodes = fs->countNodesOn(apf::Mesh::EDGE);
    if (nnodes) {
      for (int i = 0; i < 3; i++) {
	for (int n = 0; n < nnodes; n++) {
	  apf::Vector3 p;
	  m->getPoint(downedges[i], n, p);
	  entMeshCurved->setPoint(newedgesCurved[i], n, p);
	}
      }
    }

    nnodes = fs->countNodesOn(apf::Mesh::TRIANGLE);
    if (nnodes) {
      for (int n = 0; n < nnodes; n++) {
	apf::Vector3 p;
	m->getPoint(e, n, p);
	entMeshCurved->setPoint(face, n, p);
      }
    }

    entMeshCurved->acceptChanges();
    return;
  }
}

static void writeMeshes(
    apf::Mesh2* m,
    const char* prefix0,
    const char* prefix1,
    const char* prefix2,
    const char* name,
    int res)
{
  PCU_ALWAYS_ASSERT(prefix0);
  PCU_ALWAYS_ASSERT(name);
  if (!prefix1)
    PCU_ALWAYS_ASSERT(!prefix2);

  int order = m->getShape()->getOrder();
  std::stringstream ss;
  ss << prefix0 << "/";
  if (prefix1) {
    ss << prefix1;
    safe_mkdir(ss.str().c_str());
    ss << "/";
  }
  if (prefix2) {
    ss << prefix2;
    safe_mkdir(ss.str().c_str());
    ss << "/";
  }
  ss << name;

  if (order == 1) {
    apf::writeVtkFiles(ss.str().c_str(), m);
  }
  else {
    crv::writeCurvedVtuFiles(m, apf::Mesh::TRIANGLE, res, ss.str().c_str());
    crv::writeCurvedWireFrame(m, res, ss.str().c_str());
  }
  ss << ".smb";
  //m->writeNative(ss.str().c_str());
}

static void getEntIds(
    const std::vector<std::string>& ids,
    std::vector<int>& vids,
    std::vector<int>& eids,
    std::vector<int>& fids)
{
  vids.clear();
  eids.clear();
  fids.clear();

  for (std::size_t i = 0; i < ids.size(); i++) {
    std::string key = ids[i].substr(0,1);
    int value = atoi(ids[i].substr(1).c_str());
    if (key.compare(std::string("v")) == 0)
      vids.push_back(value);
    if (key.compare(std::string("e")) == 0)
      eids.push_back(value);
    if (key.compare(std::string("f")) == 0)
      fids.push_back(value);
  }

}

static void makeSurfMesh(
    apf::Mesh2* m, const char* prefix, const int res)
{
  safe_mkdir(prefix);
  writeMeshes(m, prefix, "mesh", NULL, "curved", res);

  apf::Mesh2* cavityMeshCurved = 0;
  apf::Mesh2* cavityMeshLinear = 0;

  typedef std::vector<apf::MeshEntity*> Cavity;
  typedef std::vector<apf::MeshEntity*>::iterator CavityIter;
  PCU_ALWAYS_ASSERT(!cavityMeshLinear);
  PCU_ALWAYS_ASSERT(!cavityMeshCurved);

  int dim = m->getDimension();

  Cavity icavity3;
  Cavity icavity2;
  Cavity icavity1;
  Cavity icavity0;
  icavity3.clear();
  icavity2.clear();
  icavity1.clear();
  icavity0.clear();

  std::vector<std::vector<int>> face_edge;
  std::vector<std::vector<int>> edge_vert;


  apf::MeshEntity* e;
  apf::MeshIterator* it;
  it = m->begin(2);
  int n_rc_faces = 0;
  while ( (e = m->iterate(it)) ) {
    auto gent = m->toModel(e);
    auto gdim = m->getModelType(gent);
    auto gid = m->getModelTag(gent);
    apf::Vector3 coords;
    apf::MeshEntity* dv[3];
    m->getDownward(e, 0, dv);
    double r_min=1e16;
    for (int v=0; v<3; ++v) {
      m->getPoint(dv[v], 0, coords);
      double r = std::sqrt((coords[0]*coords[0] + coords[1]*coords[1]));
      if (r < r_min) r_min = r;
    }
    if ((gdim == 2) && (gid == 3))
    //if ((gdim == 2) && (r_min > 0.8) && (std::abs(coords[0]) < 0.5) && (gid != 3) && (gid != 13))
      icavity2.push_back(e);
  }
  m->end(it);

  for (int i = 0; i < (int)icavity2.size(); i++) {
    apf::Downward dents;
    int nents = m->getDownward(icavity2[i], 1, dents);
    for (int j = 0; j < nents; j++) {
      if (std::find(icavity1.begin(), icavity1.end(), dents[j]) == icavity1.end())
      	icavity1.push_back(dents[j]);
    }
    std::vector<int> conn;
    for (int j = 0; j < nents; j++) {
      CavityIter it = std::find(icavity1.begin(), icavity1.end(), dents[j]);
      PCU_ALWAYS_ASSERT(it != icavity1.end());
      conn.push_back(std::distance(icavity1.begin(), it));
    }
    PCU_ALWAYS_ASSERT((int)conn.size() == nents);
    face_edge.push_back(conn);
  }
  PCU_ALWAYS_ASSERT(icavity2.size() == face_edge.size());

  for (int i = 0; i < (int)icavity1.size(); i++) {
    apf::Downward dents;
    int nents = m->getDownward(icavity1[i], 0, dents);
    for (int j = 0; j < nents; j++) {
      if (std::find(icavity0.begin(), icavity0.end(), dents[j]) == icavity0.end())
      	icavity0.push_back(dents[j]);
    }
    std::vector<int> conn;
    for (int j = 0; j < nents; j++) {
      CavityIter it = std::find(icavity0.begin(), icavity0.end(), dents[j]);
      PCU_ALWAYS_ASSERT(it != icavity0.end());
      conn.push_back(std::distance(icavity0.begin(), it));
    }
    PCU_ALWAYS_ASSERT((int)conn.size() == nents);
    edge_vert.push_back(conn);
  }
  PCU_ALWAYS_ASSERT(icavity1.size() == edge_vert.size());

  Cavity ocavity2linear;
  Cavity ocavity1linear;
  Cavity ocavity0linear;
  Cavity ocavity2curved;
  Cavity ocavity1curved;
  Cavity ocavity0curved;
  ocavity2linear.clear();
  ocavity1linear.clear();
  ocavity0linear.clear();
  ocavity2curved.clear();
  ocavity1curved.clear();
  ocavity0curved.clear();

  cavityMeshLinear = apf::makeEmptyMdsMesh(gmi_load(".null"), dim, false);
  cavityMeshCurved = apf::makeEmptyMdsMesh(gmi_load(".null"), dim, false);

  for (int i = 0; i < (int) icavity0.size(); i++) {
    apf::MeshEntity* ent = icavity0[i];
    apf::ModelEntity* c = m->toModel(ent);
    apf::Vector3 coords;
    apf::Vector3 params;
    m->getPoint(ent, 0, coords);
    m->getParam(ent, params);

    apf::MeshEntity* newEnt = cavityMeshLinear->createVertex(c, coords, params);
    ocavity0linear.push_back(newEnt);

    apf::MeshEntity* newEntc = cavityMeshCurved->createVertex(c, coords, params);
    ocavity0curved.push_back(newEntc);
  }
  PCU_ALWAYS_ASSERT(icavity0.size() == ocavity0linear.size());
  PCU_ALWAYS_ASSERT(icavity0.size() == ocavity0curved.size());

  for (int i = 0; i < (int) icavity1.size(); i++) {
    apf::MeshEntity* ent = icavity1[i];
    apf::ModelEntity* c = m->toModel(ent);
    apf::MeshEntity* downv[2];
    downv[0] = ocavity0linear[edge_vert[i][0]];
    downv[1] = ocavity0linear[edge_vert[i][1]];
    apf::MeshEntity* newEnt = cavityMeshLinear->createEntity(
    	apf::Mesh::EDGE, c, downv);
    ocavity1linear.push_back(newEnt);

    downv[0] = ocavity0curved[edge_vert[i][0]];
    downv[1] = ocavity0curved[edge_vert[i][1]];
    apf::MeshEntity* newEntc = cavityMeshCurved->createEntity(
    	apf::Mesh::EDGE, c, downv);
    ocavity1curved.push_back(newEntc);
  }
  PCU_ALWAYS_ASSERT(icavity1.size() == ocavity1linear.size());
  PCU_ALWAYS_ASSERT(icavity1.size() == ocavity1curved.size());

  for (int i = 0; i < (int) icavity2.size(); i++) {
    apf::MeshEntity* ent = icavity2[i];
    apf::ModelEntity* c = m->toModel(ent);
    apf::MeshEntity* downe[3];
    downe[0] = ocavity1linear[face_edge[i][0]];
    downe[1] = ocavity1linear[face_edge[i][1]];
    downe[2] = ocavity1linear[face_edge[i][2]];
    apf::MeshEntity* newEnt = cavityMeshLinear->createEntity(
    	apf::Mesh::TRIANGLE, c, downe);
    ocavity2linear.push_back(newEnt);

    downe[0] = ocavity1curved[face_edge[i][0]];
    downe[1] = ocavity1curved[face_edge[i][1]];
    downe[2] = ocavity1curved[face_edge[i][2]];
    apf::MeshEntity* newEntc = cavityMeshCurved->createEntity(
    	apf::Mesh::TRIANGLE, c, downe);
    ocavity2curved.push_back(newEntc);
  }
  PCU_ALWAYS_ASSERT(icavity2.size() == ocavity2linear.size());
  PCU_ALWAYS_ASSERT(icavity2.size() == ocavity2curved.size());

  cavityMeshLinear->acceptChanges();
  apf::deriveMdsModel(cavityMeshLinear);

  cavityMeshCurved->acceptChanges();
  apf::deriveMdsModel(cavityMeshCurved);

  cavityMeshCurved->changeShape(m->getShape(), true);
  apf::FieldShape* fs = cavityMeshCurved->getShape();

  int nnodes = fs->countNodesOn(apf::Mesh::TRIANGLE);
  if (nnodes) {
    for (int i = 0; i < (int)icavity2.size(); i++) {
      apf::MeshEntity* fromEnt = icavity2[i];
      apf::MeshEntity* toEnt   = ocavity2curved[i];
      for (int j = 0; j < nnodes; j++) {
	apf::Vector3 p;
	m->getPoint(fromEnt, j, p);
	cavityMeshCurved->setPoint(toEnt, j, p);
      }
    }
  }

  nnodes = fs->countNodesOn(apf::Mesh::EDGE);
  if (nnodes) {
    for (int i = 0; i < (int)icavity1.size(); i++) {
      apf::MeshEntity* fromEnt = icavity1[i];
      apf::MeshEntity* toEnt   = ocavity1curved[i];
      for (int j = 0; j < nnodes; j++) {
	apf::Vector3 p;
	m->getPoint(fromEnt, j, p);
	cavityMeshCurved->setPoint(toEnt, j, p);
      }
    }
  }

  cavityMeshCurved->acceptChanges();

  char cavityFolderName[128];
  char cavityFileNameLinear[128];
  char cavityFileNameCurved[128];
  char entityFileNameLinear[128];
  char entityFileNameCurved[128];
  sprintf(cavityFolderName, "%s_%05d", "triangle", 0);
  sprintf(cavityFileNameLinear, "%s", "cavity_linear");
  sprintf(cavityFileNameCurved, "%s", "cavity_curved");
  sprintf(entityFileNameLinear, "%s", "entity_linear");
  sprintf(entityFileNameCurved, "%s", "entity_curved");
  writeMeshes(cavityMeshLinear, prefix, "cavities",
      cavityFolderName, cavityFileNameLinear, res);
  writeMeshes(cavityMeshCurved, prefix, "cavities",
      cavityFolderName, cavityFileNameCurved, res);

  cavityMeshLinear->destroyNative();
  cavityMeshCurved->destroyNative();
  apf::destroyMesh(cavityMeshLinear);
  apf::destroyMesh(cavityMeshCurved);

}

static void makeInvalidMesh(
    apf::Mesh2* m, const int res)
{
  std::string prefix_s = "invalid_tetp2";
  const char * prefix = prefix_s.c_str();
  safe_mkdir(prefix);
  //writeMeshes(m, prefix, "mesh", NULL, "curved", res);
  char cavityFolderName[128];
  char cavityFileNameLinear[128];
  char cavityFileNameCurved[128];
  char entityFileNameLinear[128];
  char entityFileNameCurved[128];
  sprintf(cavityFolderName, "%s_%05d", "triangle", 0);
  sprintf(cavityFileNameLinear, "%s", "cavity_linear");
  sprintf(cavityFileNameCurved, "%s", "cavity_curved");
  sprintf(entityFileNameLinear, "%s", "entity_linear");
  sprintf(entityFileNameCurved, "%s", "entity_curved");
 
  
  apf::MeshIterator* it2 = m->begin(m->getDimension());
  apf::MeshEntity* e2;
  long n = 0;
  while ((e2 = m->iterate(it2)))
  {
    if (!isSimplex(m->getType(e2)))
      continue;
    double v = measure(m,e2);
    auto ce = getLinearCentroid(m,e2);
    if ((std::abs(ce[0]-0.779547)<.00001) && (std::abs(ce[1]-0.340564)<.00001) && (std::abs(ce[2]-0.103453)<.00001)){
    //if (v < 0)
      std::stringstream ss;
      ss << "printing invalid element, volume " << v
	<< " at " << getLinearCentroid(m, e2) << '\n';
      std::string s = ss.str();
      lion_oprint(1, "%s", s.c_str());
      fflush(stdout);
      ++n;

      apf::Mesh2* cavityMeshCurved = 0;
      apf::Mesh2* cavityMeshLinear = 0;

      typedef std::vector<apf::MeshEntity*> Cavity;
      typedef std::vector<apf::MeshEntity*>::iterator CavityIter;
      PCU_ALWAYS_ASSERT(!cavityMeshLinear);
      PCU_ALWAYS_ASSERT(!cavityMeshCurved);

      int dim = m->getDimension();

      Cavity icavity3;
      Cavity icavity2;
      Cavity icavity1;
      Cavity icavity0;
      icavity3.clear();
      icavity2.clear();
      icavity1.clear();
      icavity0.clear();

      std::vector<std::vector<int>> face_edge;
      std::vector<std::vector<int>> edge_vert;

      printf("0\n");

      apf::MeshEntity* e;
      apf::MeshIterator* it;
      apf::MeshEntity* f[4];
      m->getDownward(e2,2,f);
      int n_rc_faces = 0;
      for (int i=0; i<4; ++i) {
	e = f[i];
	auto gent = m->toModel(e);
	auto gid = m->getModelTag(gent);
	auto gdim = m->getModelType(gent);
	printf("face %d classified on Gent %d Gdim %d\n", i, gid, gdim);
	icavity2.push_back(e);
      }
      printf("1\n");

      for (int i = 0; i < (int)icavity2.size(); i++) {
	apf::Downward dents;
	int nents = m->getDownward(icavity2[i], 1, dents);
	for (int j = 0; j < nents; j++) {
	  if (std::find(icavity1.begin(), icavity1.end(), dents[j]) == icavity1.end())
	    icavity1.push_back(dents[j]);
	}
	std::vector<int> conn;
	for (int j = 0; j < nents; j++) {
	  CavityIter it = std::find(icavity1.begin(), icavity1.end(), dents[j]);
	  PCU_ALWAYS_ASSERT(it != icavity1.end());
	  conn.push_back(std::distance(icavity1.begin(), it));
	}
	PCU_ALWAYS_ASSERT((int)conn.size() == nents);
	face_edge.push_back(conn);
      }
      PCU_ALWAYS_ASSERT(icavity2.size() == face_edge.size());
      printf("2\n");

      for (int i = 0; i < (int)icavity1.size(); i++) {
	apf::Downward dents;
	int nents = m->getDownward(icavity1[i], 0, dents);
	for (int j = 0; j < nents; j++) {
	  if (std::find(icavity0.begin(), icavity0.end(), dents[j]) == icavity0.end())
	    icavity0.push_back(dents[j]);
	}
	std::vector<int> conn;
	for (int j = 0; j < nents; j++) {
	  CavityIter it = std::find(icavity0.begin(), icavity0.end(), dents[j]);
	  PCU_ALWAYS_ASSERT(it != icavity0.end());
	  conn.push_back(std::distance(icavity0.begin(), it));
	}
	PCU_ALWAYS_ASSERT((int)conn.size() == nents);
	edge_vert.push_back(conn);
      }
      PCU_ALWAYS_ASSERT(icavity1.size() == edge_vert.size());
      printf("3\n");

      Cavity ocavity2linear;
      Cavity ocavity1linear;
      Cavity ocavity0linear;
      Cavity ocavity2curved;
      Cavity ocavity1curved;
      Cavity ocavity0curved;
      ocavity2linear.clear();
      ocavity1linear.clear();
      ocavity0linear.clear();
      ocavity2curved.clear();
      ocavity1curved.clear();
      ocavity0curved.clear();

      cavityMeshLinear = apf::makeEmptyMdsMesh(gmi_load(".null"), dim, false);
      cavityMeshCurved = apf::makeEmptyMdsMesh(gmi_load(".null"), dim, false);

      printf("4\n");
      for (int i = 0; i < (int) icavity0.size(); i++) {
	apf::MeshEntity* ent = icavity0[i];
	apf::ModelEntity* c = m->toModel(ent);
	apf::Vector3 coords;
	apf::Vector3 params;
	m->getPoint(ent, 0, coords);
	m->getParam(ent, params);

	apf::MeshEntity* newEnt = cavityMeshLinear->createVertex(c, coords, params);
	ocavity0linear.push_back(newEnt);

	apf::MeshEntity* newEntc = cavityMeshCurved->createVertex(c, coords, params);
	ocavity0curved.push_back(newEntc);
      }
      PCU_ALWAYS_ASSERT(icavity0.size() == ocavity0linear.size());
      PCU_ALWAYS_ASSERT(icavity0.size() == ocavity0curved.size());

      printf("5\n");
      for (int i = 0; i < (int) icavity1.size(); i++) {
	apf::MeshEntity* ent = icavity1[i];
	apf::ModelEntity* c = m->toModel(ent);
	apf::MeshEntity* downv[2];
	downv[0] = ocavity0linear[edge_vert[i][0]];
	downv[1] = ocavity0linear[edge_vert[i][1]];
	apf::MeshEntity* newEnt = cavityMeshLinear->createEntity(
	    apf::Mesh::EDGE, c, downv);
	ocavity1linear.push_back(newEnt);

	downv[0] = ocavity0curved[edge_vert[i][0]];
	downv[1] = ocavity0curved[edge_vert[i][1]];
	apf::MeshEntity* newEntc = cavityMeshCurved->createEntity(
	    apf::Mesh::EDGE, c, downv);
	ocavity1curved.push_back(newEntc);
      }
      PCU_ALWAYS_ASSERT(icavity1.size() == ocavity1linear.size());
      PCU_ALWAYS_ASSERT(icavity1.size() == ocavity1curved.size());
      printf("6\n");

      for (int i = 0; i < (int) icavity2.size(); i++) {
	apf::MeshEntity* ent = icavity2[i];
	apf::ModelEntity* c = m->toModel(ent);
	apf::MeshEntity* downe[3];
	downe[0] = ocavity1linear[face_edge[i][0]];
	downe[1] = ocavity1linear[face_edge[i][1]];
	downe[2] = ocavity1linear[face_edge[i][2]];
	apf::MeshEntity* newEnt = cavityMeshLinear->createEntity(
	    apf::Mesh::TRIANGLE, c, downe);
	ocavity2linear.push_back(newEnt);

	downe[0] = ocavity1curved[face_edge[i][0]];
	downe[1] = ocavity1curved[face_edge[i][1]];
	downe[2] = ocavity1curved[face_edge[i][2]];
	apf::MeshEntity* newEntc = cavityMeshCurved->createEntity(
	    apf::Mesh::TRIANGLE, c, downe);
	ocavity2curved.push_back(newEntc);
      }
      PCU_ALWAYS_ASSERT(icavity2.size() == ocavity2linear.size());
      PCU_ALWAYS_ASSERT(icavity2.size() == ocavity2curved.size());
      printf("7\n");

      //cavityMeshLinear->acceptChanges();
      //apf::deriveMdsModel(cavityMeshLinear);

      //cavityMeshCurved->acceptChanges();
      //apf::deriveMdsModel(cavityMeshCurved);

      cavityMeshCurved->changeShape(m->getShape(), true);
      apf::FieldShape* fs = cavityMeshCurved->getShape();

      int nnodes = fs->countNodesOn(apf::Mesh::TRIANGLE);
      printf("7.5\n");
      if (nnodes) {
	for (int i = 0; i < (int)icavity2.size(); i++) {
	  apf::MeshEntity* fromEnt = icavity2[i];
	  apf::MeshEntity* toEnt   = ocavity2curved[i];
	  for (int j = 0; j < nnodes; j++) {
	    apf::Vector3 p;
	    m->getPoint(fromEnt, j, p);
	    cavityMeshCurved->setPoint(toEnt, j, p);
	  }
	}
      }

      printf("8\n");
      nnodes = fs->countNodesOn(apf::Mesh::EDGE);
      if (nnodes) {
	for (int i = 0; i < (int)icavity1.size(); i++) {
	  apf::MeshEntity* fromEnt = icavity1[i];
	  apf::MeshEntity* toEnt   = ocavity1curved[i];
	  for (int j = 0; j < nnodes; j++) {
	    apf::Vector3 p;
	    m->getPoint(fromEnt, j, p);
	    cavityMeshCurved->setPoint(toEnt, j, p);
	  }
	}
      }

      printf("9\n");
      //cavityMeshCurved->acceptChanges();

      printf("10\n");
      //writeMeshes(cavityMeshLinear, prefix, "cavities",
	  //cavityFolderName, cavityFileNameLinear, res);
      writeMeshes(cavityMeshCurved, prefix, "cavities",
	  NULL, "cavity_curved", res);
      printf("11\n");

      cavityMeshLinear->destroyNative();
      cavityMeshCurved->destroyNative();
      apf::destroyMesh(cavityMeshLinear);
      apf::destroyMesh(cavityMeshCurved);

      printf("printed invalid mesh\n");
    }
  }
  m->end(it2);
  PCU_Barrier();

}

namespace crv {

void convertInterpolationPoints(int n, int ne,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<double>& c,
    apf::NewArray<apf::Vector3>& newNodes){

  for(int i = 0; i < ne; ++i)
    newNodes[i].zero();

  for( int i = 0; i < ne; ++i)
    for( int j = 0; j < n; ++j)
      newNodes[i] += nodes[j]*c[i*n+j];

}

void convertInterpolationPoints(apf::Mesh2* m, apf::MeshEntity* e,
    int n, int ne, apf::NewArray<double>& c){

  apf::NewArray<apf::Vector3> l, b(ne);
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),e);
  apf::getVectorNodes(elem,l);

  crv::convertInterpolationPoints(n,ne,l,c,b);

  for(int i = 0; i < ne; ++i)
    m->setPoint(e,i,b[i]);

  apf::destroyElement(elem);
}

void interpolatingToBezier(apf::Mesh2* m)
{
  apf::FieldShape * fs = m->getShape();
  int order = fs->getOrder();

  int md = m->getDimension();
  int blendingOrder = getBlendingOrder(apf::Mesh::simplexTypes[md]);
  // go downward, and convert interpolating to control points
  int startDim = md - (blendingOrder > 0);

  for(int d = startDim; d >= 1; --d){
    if(!fs->hasNodesIn(d)) continue;
    int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
    apf::NewArray<double> c;
    getBezierTransformationCoefficients(order,
        apf::Mesh::simplexTypes[d],c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(d);
    while ((e = m->iterate(it))){
      if(m->isOwned(e))
        convertInterpolationPoints(m,e,n,ne,c);
    }
    m->end(it);
  }
  // if we have a full representation, we need to place internal nodes on
  // triangles and tetrahedra
  for(int d = 2; d <= md; ++d){
    if(!fs->hasNodesIn(d) ||
        getBlendingOrder(apf::Mesh::simplexTypes[d])) continue;
    int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
    apf::NewArray<double> c;
    getInternalBezierTransformationCoefficients(m,order,1,
        apf::Mesh::simplexTypes[d],c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(d);
    while ((e = m->iterate(it))){
      if(!isBoundaryEntity(m,e) && m->isOwned(e))
        convertInterpolationPoints(m,e,n-ne,ne,c);
    }
    m->end(it);
  }
  apf::synchronize(m->getCoordinateField());

}

void snapToInterpolate(apf::Mesh2* m, apf::MeshEntity* e, bool isNew)
{
  PCU_ALWAYS_ASSERT(m->canSnap());
  int type = m->getType(e);
  if(type == apf::Mesh::VERTEX){
    apf::Vector3 p, pt(0,0,0);
    apf::ModelEntity* g = m->toModel(e);
    m->getParamOn(g,e,p);
    m->snapToModel(g,p,pt);
    m->setPoint(e,0,pt);
    return;
  }
  // e is an edge or a face
  // either way, get a length-scale by computing
  // the distance b/w first two downward verts
  apf::MeshEntity* down[12];
  m->getDownward(e, 0, down);
  apf::Vector3 p0, p1;
  m->getPoint(down[0], 0, p0);
  m->getPoint(down[1], 0, p1);
  double lengthScale = (p1 - p0).getLength();
  apf::FieldShape * fs = m->getShape();
  int non = fs->countNodesOn(type);
  apf::Vector3 p, xi, pt0, pt(0,0,0);
  if (type == 2) {
    if (m->getModelType(m->toModel(e)) != 2) {
      printf("g id %d, class dim\n",m->getModelTag(m->toModel(e)),m->getModelType(m->toModel(e)));
      fail("error: cannot interior face\n");
    }
  }
  for(int i = 0; i < non; ++i){
    apf::ModelEntity* g = m->toModel(e);
    fs->getNodeXi(type,i,xi);
    if(type == apf::Mesh::EDGE)
      transferParametricOnEdgeSplit(m,e,0.5*(xi[0]+1.),p);
    else
      transferParametricOnTriSplit(m,e,xi,p);
    m->snapToModel(g,p,pt);
    if (isNew || !m->canGetClosestPoint()) {
      m->setPoint(e,i,pt);
      continue;
    }
    m->getPoint(e,i,pt0);
    if (!m->isOnModel(g, pt0, lengthScale))
      m->setPoint(e,i,pt);
  }

  // sample model entity into m points for fitting
  // consider the 2d gauss/circle bump cases
  {
    const double m20_n2_Binv [3][21] = {
      0.356295878035009,0.289666854884246,0.228684359119142,0.173348390739695,0.123658949745906,0.0796160361377752,0.041219649915302,0.00846979107848662,-0.0186335403726709,-0.0400903444381706,-0.0559006211180125,-0.0660643704121966,-0.0705815923207228,-0.0694522868435912,-0.0626764539808018,-0.0502540937323546,-0.0321852060982496,-0.00846979107848674,0.020892151326934,0.0559006211180125,0.0965556182947488,
      -0.338226990400904,-0.222473178994918,-0.118903979315879,-0.027519391363785,0.0516805848613629,0.118695949359565,0.173526702130821,0.216172843175132,0.246634372492496,0.264911290082915,0.271003595946388,0.264911290082915,0.246634372492496,0.216172843175131,0.173526702130821,0.118695949359565,0.0516805848613628,-0.027519391363785,-0.118903979315879,-0.222473178994918,-0.338226990400904,
      0.0965556182947488,0.0559006211180123,0.0208921513269339,-0.00846979107848674,-0.0321852060982496,-0.0502540937323546,-0.0626764539808018,-0.0694522868435912,-0.0705815923207228,-0.0660643704121966,-0.0559006211180124,-0.0400903444381705,-0.0186335403726708,0.00846979107848673,0.0412196499153021,0.0796160361377753,0.123658949745906,0.173348390739695,0.228684359119142,0.289666854884246,0.356295878035009
    };
    const double m20_n3_Binv [4][21] = {
      0.544042913608131,0.364765669113495,0.220779220779221,0.108789760963674,0.0255034820252211,-0.0323734236777715,-0.0681347637869376,-0.0850743459439112,-0.0864859777903257,-0.0756634669678148,-0.0559006211180124,-0.0304912478825521,-0.00272915490306781,0.0240918501788068,0.0466779597214381,0.0617353660831923,0.0659702616224356,0.0560888386975344,0.0287972896668549,-0.0191981931112365,-0.0911914172783739,
      -0.776094276094276,-0.319509797770667,0.0251430068363706,0.2696075571743,0.425627272690583,0.504945572832682,0.519305877048058,0.480451604784175,0.400126175488495,0.290073008608478,0.162035523591588,0.027757139885286,-0.101018723062964,-0.212548645805701,-0.295089208895464,-0.336896992884789,-0.326228578326214,-0.251340545772277,-0.100489475775517,0.138068051111529,0.476075454336324,
      0.476075454336324,0.138068051111529,-0.100489475775517,-0.251340545772277,-0.326228578326214,-0.336896992884789,-0.295089208895464,-0.212548645805702,-0.101018723062965,0.0277571398852864,0.162035523591588,0.290073008608478,0.400126175488495,0.480451604784175,0.519305877048059,0.504945572832682,0.425627272690583,0.2696075571743,0.0251430068363704,-0.319509797770667,-0.776094276094276,
      -0.0911914172783736,-0.0191981931112364,0.028797289666855,0.0560888386975342,0.0659702616224355,0.061735366083192,0.0466779597214378,0.0240918501788066,-0.00272915490306797,-0.0304912478825524,-0.0559006211180125,-0.0756634669678148,-0.0864859777903257,-0.0850743459439111,-0.0681347637869377,-0.0323734236777714,0.0255034820252213,0.108789760963674,0.220779220779221,0.364765669113495,0.544042913608131
    };
    const double m20_n4_Binv [5][21] = {
      0.708187473409206,0.36476565709375,0.134387357502348,-0.0063994094316919,-0.0786749402606679,-0.101148135835206,-0.0901562094915379,-0.0596649788642247,-0.0212685740734334,0.0158102704623412,0.0447205078389975,0.0609824868541059,0.0624882438196302,0.0495012107492064,0.0246565071708638,-0.00703935168569621,-0.0382081630763789,-0.0591003285554174,-0.0575945621626516,-0.0191981822362491,0.0729531804000156,
      -1.33358116312781,-0.148440854630516,0.541168380325945,0.852225028526053,0.88888366789942,0.742477180569769,0.491514539805945,0.20168302307091,-0.0741540010272527,-0.295954859686346,-0.43650178895905,-0.481398720703929,-0.429073495634426,-0.290775650269871,-0.0905786299844687,0.134622424041689,0.335108461773633,0.448338737370508,0.398948596136583,0.098751687570239,-0.553262247685018,
      1.24696035039086,-0.0907210565663788,-0.772920278932157,-0.971196238290646,-0.836922106520793,-0.501290293637357,-0.0753074599498626,0.350200496096352,0.704593425761185,0.937411942168717,1.01838240814827,0.93741194839334,0.704593437302672,0.35020051113919,-0.0753074441289471,-0.501290280669394,-0.836922100944569,-0.971196245552705,-0.772920305386802,-0.0907211094756714,1.24696026285711,
      -0.553262494554797,0.0987516827957529,0.398948437634559,0.448338808500357,0.335108364817624,0.134622541220603,-0.090578696770731,-0.29077551842663,-0.429073562131614,-0.481398601060433,-0.436501877502103,-0.295954768535871,-0.0741541203552493,0.201683086056018,0.491514400599922,0.742477238388223,0.888883545418415,0.85222513289776,0.541168342919258,-0.148440617214326,-1.33358100931449,
      0.0729531804000153,-0.0191982998242995,-0.0575944938682612,-0.0591004031331327,-0.0382080713448737,-0.00703941730433604,0.0246565978629322,0.0495011379564905,0.062488324451202,0.0609824037470375,0.0447205799192716,0.0158101819682861,-0.0212685034302329,-0.059665064675492,-0.0901561324913942,-0.101148212676735,-0.0786748533550052,-0.00639947772458944,0.134387448691432,0.364765585295087,0.708187549163706
    };
    const double m20_n5_Binv [6][21] = {
      0.831646510782998,0.303036225676515,0.0239238620797841,-0.0897882595217121,-0.103774504008318,-0.0672895004752876,-0.0151761057174141,0.0301262787455872,0.0565782793374601,0.0605307533056197,0.0447205986603458,0.0162620784667895,-0.0153584906276285,-0.0402900906661246,-0.0503236372902417,-0.040898261106781,-0.0131088114781731,0.0242887309932694,0.0528686797690083,0.0425315533251177,-0.050505822979151,
      -1.91395372457395,0.448562952190189,1.34446356098562,1.34832246795317,0.896382518890465,0.302594171607156,-0.225301060225096,-0.569683238259049,-0.687016549243884,-0.591747141705648,-0.340257753768803,-0.0147769080648868,0.292677818674148,0.496376999642731,0.532949689432892,0.377449830541904,0.0594351395263026,-0.320975041173248,-0.577010059558813,-0.419203703314889,0.560709084929318,
      2.38361933676342,-1.19824510370117,-2.18786173975884,-1.70676589637661,-0.58753772019111,0.591029203808873,1.46864144994903,1.86820785546364,1.76067582984413,1.22969462534779,0.436429944595957,-0.415742521499404,-1.11645275355536,-1.48360603033217,-1.39863542401224,-0.841742801406701,0.0728266233994717,1.061565345207,1.63640038964965,1.0695112048131,-1.64201035227704,
      -1.64200357680825,1.06950322472808,1.63640422278886,1.06155453872438,0.0728280000842051,-0.841735561980836,-1.39863625472873,-1.48360338046489,-1.11645529001363,-0.415741350915918,0.436426674576716,1.22969477934522,1.76067321536832,1.86820922710839,1.46864097193397,0.591034755170837,-0.587535085638718,-1.7067749095215,-2.18785638947945,-1.19825265773246,2.38362449333281,
      0.560705898529053,-0.419199548677787,-0.577011849903559,-0.320969384408559,0.0594346119978461,0.377444988531286,0.532950355577265,0.496375363147522,0.292679384320605,-0.0147777104823279,-0.340255912546467,-0.591746527252507,-0.687015306460646,-0.569683471742696,-0.22530127592314,0.302591833068326,0.896380447098189,1.34832623753409,1.34446033076671,0.448566777071761,-1.91395550302515,
      -0.0505058229791512,0.0425309625445233,0.0528692721106524,0.0242881878460195,-0.0131086410926881,-0.040897247682369,-0.0503239478473577,-0.0402901197429378,-0.0153589478604045,0.0162622129671662,0.0447203340041208,0.060530520695455,0.0565783370170216,0.0301264590589881,-0.0151758461241843,-0.0672890631941011,-0.103774280003404,-0.0897891999511792,0.0239240028058856,0.303035516498323,0.831647191032501
    };
    const double m20_n6_Binv [7][21] = {
      0.912705227609145,0.21387192666163,-0.056281650401418,-0.101308035289009,-0.0536333025501318,0.00893865840360135,0.0486673095456289,0.0552342998998059,0.0350463486506416,0.00246008858985797,-0.0270528501534785,-0.0418097313478337,-0.0368905775854432,-0.0151824271240998,0.0135196056666986,0.0353299416538685,0.0370323787247057,0.0127695722080586,-0.0273366757525609,-0.0466331299670386,0.030552746106313,
      -2.42872478680291,1.49391350251802,2.0865106489112,1.24682688185614,0.128199333492804,-0.673477798945145,-0.956139658223423,-0.770905707451444,-0.304788033342446,0.213576101193128,0.584892767485453,0.687020441940358,0.499633134357767,0.105741313618309,-0.330122631688456,-0.606703599495369,-0.554143722586158,-0.125245033092819,0.490110288750932,0.72736464818203,-0.513536923683858,
      3.73189186859999,-3.70817543422675,-3.76195143698865,-1.08359970251796,1.6272530048478,3.1099642977715,3.09420955324071,1.9169231734604,0.20610128206057,-1.36962477832278,-2.28472764988548,-2.27430393831166,-1.38541616781499,0.0377376737175841,1.43544301090974,2.17975207885321,1.78850612950278,0.205531921013922,-1.85294539468312,-2.48561538396767,1.87304470424711,
      -3.45474462249615,4.14372048572353,3.50956668398318,0.221068803719308,-2.62377514749855,-3.72292515136169,-2.97809019384233,-0.992650479726901,1.33830821017565,3.1476117232658,3.82377110424121,3.14768173354648,1.33830737491896,-0.992657168795306,-2.97808693185153,-3.72295850413978,-2.62376345818023,0.221053377384503,3.50957945842844,4.14374944729894,-3.45476638219011,
      1.87302790598009,-2.48559216935262,-1.85293655735757,0.205516971758034,1.78851780907937,2.1797294840911,1.43544860156137,0.0377274971268636,-1.38541580479723,-2.27425216071165,-2.28472799343434,-1.36967656378274,0.206102653954808,1.91692388147164,3.09421001540408,3.10999159506315,1.62724718435426,-1.0835932294722,-3.76196167697789,-3.70819369581154,3.73190641346092,
      -0.513532638583127,0.727357419262641,0.490108572433856,-0.125238593472584,-0.554149395400495,-0.606696223827648,-0.330126932673499,0.105746800960618,0.499632233414001,0.687001319671445,0.58489357152371,0.213595599090282,-0.304788087821029,-0.770904587911242,-0.956140752713292,-0.673489223796926,0.128200163182265,1.246826236147,2.08651490719719,1.49391878225167,-2.42872991225011,
      0.0305527461063126,-0.0466329218274189,-0.0273368405973939,0.0127688176041216,0.0370332550981113,0.0353294847033453,0.013520537725021,-0.0151831761960567,-0.0368903455079547,-0.0418078836675619,-0.0270531815580689,0.00245821694319088,0.0350460782723201,0.0552339847926492,0.0486674348199609,0.0089401265266058,-0.0536331957917287,-0.101307940437496,-0.0562821723417539,0.213871579067842,0.912705527962718
    };
    const int m_dataPts = 20;
    double data_x[m_dataPts+1];
    double data_y[m_dataPts+1];
    double data_z[m_dataPts+1];
    //double data_y[m_dataPts+1];
    //FILE *fp;
    //fp = fopen("/users/joshia5/Documents/curved/matlab-files/cubic_bezier_pInv/m20_n4_B.txt","w");
    if (m->getModelType(m->toModel(e)) == 1) {
      if ((m->getModelTag(m->toModel(e))) == 19) {   
        for(int i = 0; i <= m_dataPts; ++i){ // ignore the pts on vtx
          apf::ModelEntity* g = m->toModel(e);
          //fs->getNodeXi(type,i,xi);
          if(type == apf::Mesh::EDGE) {
            xi[0] = i/(1.*m_dataPts);
            transferParametricOnEdgeSplit(m,e,xi[0],p);
          }
          //else
          //transferParametricOnTriSplit(m,e,xi,p);
          m->snapToModel(g,p,pt);
          data_x[i] = pt[0];
          data_y[i] = pt[1];
          data_z[i] = pt[2];
          //printf("%.16f,%.16f,%.16f\n", std::pow((1.-xi[0]),2), 2*xi[0]*(1.-xi[0]), std::pow(xi[0],2)); // p2
          //printf("%.16f,%.16f,%.16f,%.16f\n", std::pow((1.-xi[0]),3), 3*xi[0]*std::pow((1.-xi[0]),2), 3*(1.-xi[0])*std::pow(xi[0],2), std::pow(xi[0],3)); // p3
          //fprintf(stderr,"%.16f,%.16f,%.16f,%.16f,%16f\n", std::pow((1.-xi[0]),4.), 4.*xi[0]*std::pow((1.-xi[0]),3.), 6.*std::pow(xi[0],2.)*std::pow((1.-xi[0]),2.), 4.*std::pow(xi[0],3.)*(1.-xi[0]), 1.*std::pow(xi[0],4.)); // p4
          //fprintf(stderr,"%.16f,%.16f,%.16f,%.16f,%16f,%16f\n", std::pow((1.-xi[0]),5.), 5.*xi[0]*std::pow((1.-xi[0]),4.), 10.*std::pow(xi[0],2.)*std::pow((1.-xi[0]),3.), 10.*std::pow(xi[0],3.)*std::pow((1.-xi[0]),2.), 5.*std::pow(xi[0],4.)*(1.-xi[0]), std::pow(xi[0],5.)); // p5
          fprintf(stderr,"%.16f,%.16f,%.16f,%.16f,%16f,%16f,%16f\n", std::pow((1.-xi[0]),6.), 6.*xi[0]*std::pow((1.-xi[0]),5.), 15.*std::pow(xi[0],2.)*std::pow((1.-xi[0]),4.), 20.*std::pow(xi[0],3.)*std::pow((1.-xi[0]),3.), 15.*std::pow(xi[0],4.)*std::pow((1.-xi[0]),2.), 6.*std::pow(xi[0],5.)*(1.-xi[0]), std::pow(xi[0],6.)); // p6
          //printf("%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n", pt[0], pt[1], pt[2],
              //std::pow((1.-xi[0]),3), 3*xi[0]*std::pow((1.-xi[0]),2), 3*(1.-xi[0])*std::pow(xi[0],2), std::pow(xi[0],3), xi[0]);
          //if (isNew || !m->canGetClosestPoint()) {
          //m->setPoint(e,i,pt);
          //continue;
          //}
          //m->getPoint(e,i,pt0);
          //if (!m->isOnModel(g, pt0, lengthScale))
          //m->setPoint(e,i,pt);
        }
        double c_x[2], c_y[2], c_z[2];
        const int order=3; // cubic for now
        for (int i=0; i<order-1; ++i) {
          c_x[i] = 0.;
          c_y[i] = 0.;
          c_z[i] = 0.;
          for (int j=0; j<=m_dataPts; ++j) {
            c_x[i] += m20_n3_Binv[i+1][j]*data_x[j];
            c_y[i] += m20_n3_Binv[i+1][j]*data_y[j];
            c_z[i] += m20_n3_Binv[i+1][j]*data_z[j];
          }
        }
        printf("computed control points x: %.16f, %.16f\n",c_x[0],c_x[1]);
        printf("computed control points y: %.16f, %.16f\n",c_y[0],c_y[1]);
        printf("computed control points z: %.16f, %.16f\n",c_z[0],c_z[1]);
      }
    }
    //fclose(fp);
  }

}

void MeshCurver::synchronize()
{
  // this causes the matched entities to collapse onto each other.
  // In other words, the matched vertexes will have the same coords
  // after the following call. This is not a desired behavior.
  // TODO: fix this!
  apf::synchronize(m_mesh->getCoordinateField());
}

void MeshCurver::snapToInterpolate(int dim)
{
  PCU_ALWAYS_ASSERT(m_mesh->canSnap());
  apf::MeshEntity* e;
  apf::MeshIterator* it = m_mesh->begin(dim);
  while ((e = m_mesh->iterate(it))) {
    if(isBoundaryEntity(m_mesh,e) && m_mesh->isOwned(e))
      crv::snapToInterpolate(m_mesh,e);
  }
  m_mesh->end(it);
}

bool InterpolatingCurver::run()
{
  if (!m_mesh->canSnap())
    fail("Cannot snap to geometry, "
        "this operation is pointless.\n");
  // interpolate points in each dimension
  for(int d = 1; d < 2; ++d)
    snapToInterpolate(d);

  synchronize();

  m_mesh->acceptChanges();
  m_mesh->verify();
  return true;
}

void BezierCurver::convertInterpolatingToBezier()
{
  interpolatingToBezier(m_mesh);
}

bool BezierCurver::run()
{
  std::string name = m_mesh->getShape()->getName();
  if(m_order < 1 || m_order > 6){
    fail("trying to convert to unimplemented Bezier order\n");
  }
  // if its already bezier, check what needs to be done, if anything
  if(name == std::string("Bezier")){
    changeMeshOrder(m_mesh,m_order);
    return true;
  } else {
    // project the new mesh onto the old, with interpolating shapes
    apf::changeMeshShape(m_mesh, getBezier(m_order),true);
  }

  if (m_mesh->canSnap()){
    for(int d = 1; d <= 2; ++d)
      snapToInterpolate(d);
    synchronize();
  }

  convertInterpolatingToBezier();

  if( m_mesh->getDimension() >= 2 && m_order == 2){
    //ma::Input* shapeFixer = configureShapeCorrection(m_mesh);
    //crv::adapt(shapeFixer);
  }

  m_mesh->acceptChanges();
  m_mesh->verify();
  if (m_mesh->canSnap()) {
    crv::computeMeanDist(m_mesh);
  }

  makeSurfMesh(m_mesh, "SurfMesh_outerfluxsurf", 15);
  makeInvalidMesh(m_mesh, 15);
  synchronize();
  PCU_Barrier();

  return true;
}

void GregoryCurver::setCubicEdgePointsUsingNormals()
{
  apf::MeshEntity* e;
  apf::Vector3 p, xi, pt;
  apf::Vector3 points[4];

  apf::MeshIterator* it = m_mesh->begin(1);

  while ((e = m_mesh->iterate(it))) {
    apf::ModelEntity* g = m_mesh->toModel(e);

    if(m_mesh->getModelType(g) == 3) continue;
    // set edges using normals
    if(m_mesh->getModelType(g) == 1) {

      apf::Vector3 t[2];
      apf::MeshEntity* v[2];
      m_mesh->getDownward(e,0,v);
      for(int i = 0; i < 2; ++i){
        m_mesh->getPoint(v[i],0,points[i*3]);
        m_mesh->getParamOn(g,v[i],p);
        m_mesh->getFirstDerivative(g,p,t[i],t[i]);
        t[i] = t[i].normalize();
      }
      double d = (points[3]-points[0]).getLength();
      apf::Vector3 l = (points[3]-points[0])/d;
      points[1] = points[0] + t[0]*(l*t[0])/fabs(l*t[0])*d/3.;
      points[2] = points[3] - t[1]*(l*t[1])/fabs(l*t[1])*d/3.;
      for(int i = 0; i < 2; ++i)
        m_mesh->setPoint(e,i,points[i+1]);

    } else {
      // set edges using tangents
      apf::Vector3 n[2];
      apf::MeshEntity* v[2];
      m_mesh->getDownward(e,0,v);
      for(int i = 0; i < 2; ++i){
        m_mesh->getPoint(v[i],0,points[i*3]);
        m_mesh->getParamOn(g,v[i],p);
        m_mesh->getNormal(g,p,n[i]);
      }
      double d = (points[3]-points[0]).getLength();
      apf::Vector3 l = (points[3]-points[0])/d;
      double a[3] = {n[0]*l,n[1]*l,n[0]*n[1]};

      double r = 6.*(2.*a[0]+a[2]*a[1])/(4.-a[2]*a[2]);
      double s = 6.*(2.*a[1]+a[2]*a[0])/(4.-a[2]*a[2]);
      points[1] = points[0] + (l*6. - n[0]*2.*r+ n[1]*s)   *d/18.;
      points[2] = points[3] - (l*6. + n[0]*r   - n[1]*2.*s)*d/18.;
      for(int i = 0; i < 2; ++i)
        m_mesh->setPoint(e,i,points[i+1]);

    }
  }
  m_mesh->end(it);
}

static void elevateBezierCurves(apf::Mesh2* m)
{

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(1);
  while ((e = m->iterate(it))) {
    if(isBoundaryEntity(m,e))
      elevateBezierCurve(m,e,3,1);
  }
  m->end(it);
}

void GregoryCurver::setInternalPointsLocally()
{
  apf::Vector3 D[3][4];
  apf::Vector3 W[3][3];
  apf::Vector3 A[3][3];

  double lam[3][2];
  double mu[3][2];
  apf::Vector3 G[6];

  apf::MeshEntity* e;
  apf::MeshIterator* it = m_mesh->begin(2);
  while ((e = m_mesh->iterate(it))) {
    apf::ModelEntity* g = m_mesh->toModel(e);
    if(!m_mesh->isOwned(e) || m_mesh->getModelType(g) != 2) continue;

    apf::Vector3 n[3];
    apf::MeshEntity* verts[3];
    apf::MeshEntity* edges[3];
    m_mesh->getDownward(e,0,verts);
    m_mesh->getDownward(e,1,edges);

    // elevated edges
    apf::NewArray<apf::Vector3> q(12);

    for(int i = 0; i < 3; ++i){
      apf::Vector3 param;
      m_mesh->getPoint(verts[i],0,q[i]);
      m_mesh->getParamOn(g,verts[i],param);
      m_mesh->getNormal(g,param,n[i]);
    }

    // elevate the edge points without formally setting them to q
    // compute tangent vectors, W
    for(int i = 0; i < 3; ++i){
      apf::Element* edge =
          apf::createElement(m_mesh->getCoordinateField(),edges[i]);
      apf::NewArray<apf::Vector3> ep;
      apf::getVectorNodes(edge,ep);

      bool flip;
      int which, rotate;
      apf::getAlignment(m_mesh,e,edges[i],which,flip,rotate);

      if(flip){
        W[i][0] = ep[3]-ep[1];
        W[i][1] = ep[2]-ep[3];
        W[i][2] = ep[0]-ep[2];
        q[i*3+3] = ep[1]*0.25+ep[3]*0.75;
        q[i*3+4] = ep[3]*0.5+ep[2]*0.5;
        q[i*3+5] = ep[2]*0.75+ep[0]*0.25;
      } else {
        W[i][0] = ep[2]-ep[0];
        W[i][1] = ep[3]-ep[2];
        W[i][2] = ep[1]-ep[3];
        q[i*3+3] = ep[0]*0.25+ep[2]*0.75;
        q[i*3+4] = ep[2]*0.5+ep[3]*0.5;
        q[i*3+5] = ep[3]*0.75+ep[1]*0.25;
      }

      apf::destroyElement(edge);
    }
    int const (*tev)[2] = apf::tri_edge_verts;

    for(int i = 0; i < 3; ++i){
      A[i][0] = apf::cross(n[tev[i][0]],W[i][0].normalize());
      A[i][2] = apf::cross(n[tev[i][1]],W[i][2].normalize());
      A[i][1] = (A[i][0]+A[i][2]).normalize();
    }

    D[0][0] = q[11] - (q[0]+q[3] )*0.5;
    D[0][3] = q[6]  - (q[1]+q[5] )*0.5;

    D[1][0] = q[5]  - (q[1]+q[6] )*0.5;
    D[1][3] = q[9]  - (q[2]+q[8] )*0.5;

    D[2][0] = q[8]  - (q[2]+q[9] )*0.5;
    D[2][3] = q[3]  - (q[0]+q[11])*0.5;

    for(int i = 0; i < 3; ++i){
      lam[i][0] = D[i][0]*W[i][0]/(W[i][0]*W[i][0]);
      lam[i][1] = D[i][3]*W[i][2]/(W[i][2]*W[i][2]);

      mu[i][0]  = D[i][0]*A[i][0];
      mu[i][1]  = D[i][3]*A[i][2];
    }

    for(int i = 0; i < 3; ++i){
      G[i] = (q[i*3+3]+q[i*3+4])*0.5
          + W[i][1]*2./3.*lam[i][0] + W[i][0]*1./3.*lam[i][1]
          + A[i][1]*2./3.*mu[i][0] + A[i][0]*1./3.*mu[i][1];
      G[3+i] = (q[i*3+4]+q[i*3+5])*0.5
          + W[i][2]*1./3.*lam[i][0] + W[i][1]*2./3.*lam[i][1]
          + A[i][2]*1./3.*mu[i][0] + A[i][1]*2./3.*mu[i][1];
    }
    for(int i = 0; i < 6; ++i)
      m_mesh->setPoint(e,i,G[i]);

  }
  m_mesh->end(it);
}

bool GregoryCurver::run()
{
  if(m_order != 4){
    fail("cannot only convert to G1 of order 4\n");
  }
  if(m_mesh->getDimension() != 3){
    fail("can only convert 3D mesh\n");
  }
  if (!m_mesh->canSnap()){
     fail("Cannot snap to geometry, "
         "cannot convert mesh to G1.\n");
  }

  apf::changeMeshShape(m_mesh, getGregory(),true);
  int md = m_mesh->getDimension();
  apf::FieldShape * fs = m_mesh->getShape();

  // interpolate points in each dimension
  for(int d = 1; d < 2; ++d)
    snapToInterpolate(d);

  synchronize();

  // go downward, and convert interpolating to control points
  for(int d = md; d >= 1; --d){
    if(!fs->hasNodesIn(d)) continue;

    int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
    apf::NewArray<apf::Vector3> l, b(ne);

    apf::NewArray<double> c;

    getGregoryTransformationCoefficients(apf::Mesh::simplexTypes[d],c);

    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);

    while ((e = m_mesh->iterate(it))) {
      if(m_mesh->isOwned(e))
        convertInterpolationPoints(m_mesh,e,n,ne,c);
    }
    m_mesh->end(it);
  }

  setCubicEdgePointsUsingNormals();
  setInternalPointsLocally();

  elevateBezierCurves(m_mesh);

  for(int d = 2; d <= md; ++d){
    if(!fs->hasNodesIn(d) ||
        getBlendingOrder(apf::Mesh::simplexTypes[d])) continue;
    int type = apf::Mesh::simplexTypes[d];
    int n = fs->getEntityShape(type)->countNodes();
    int ne = fs->countNodesOn(type);
    apf::NewArray<double> c;
    getGregoryBlendedTransformationCoefficients(1,type,c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);
    while ((e = m_mesh->iterate(it))){
      if(!isBoundaryEntity(m_mesh,e) && m_mesh->isOwned(e))
        convertInterpolationPoints(m_mesh,e,n-ne,ne,c);
    }
    m_mesh->end(it);
  }

  synchronize();

  m_mesh->acceptChanges();
  m_mesh->verify();
  return true;
}

void computeMeanDist(apf::Mesh2* m) {
  PCU_ALWAYS_ASSERT(m->canSnap());
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(2);
  double mean_dist = 0.;
  double max_dist = -1.e32;
  int rc_faces = 0;

  int dim = m->getDimension();
  apf::Field* crd_field = m->getCoordinateField();
  apf::FieldShape *fs = m->getShape();
  int s_order = fs->getOrder();
  printf("order of coordinate field %d\n", s_order);
      
  double const xi_p6tri[28*2] = {
  0.,0.,  0.0805979,0.,  0.2650895,0.,  0.5,0., 0.7349105,0.,  0.9194021,0.,  1.,0.,
  0.,0.0805979,  0.0805979,0.0805979,  0.02650895,0.0805979,  0.5,0.0805979,  
    0.7349105,0.0805979, 0.9194021,0.0805979,
  0.,0.2650895,  0.0805979,0.2650895,  0.2650895,0.2650895,  0.5,0.2650895,
    0.7349105,0.2650895,
  0.,0.5,  0.0805979,0.5,  0.2650895,0.5,  0.5,0.5,
  0., 0.7349105,  0.0805979,0.7349105,  0.2650895,0.7349105,
  0.,0.9194021,  0.0805979,0.9194021,
  0.,1.
  };

  while ((e = m->iterate(it))){
    int type = m->getType(e);
    auto gent = m->toModel(e);
    auto gdim = m->getModelType(gent);
    if (gdim == (dim-1)) {
      apf::MeshElement* mE = apf::createMeshElement(m, e);
      apf::Element* elem = apf::createElement(crd_field, mE); 
      int bad_disp = -1;

      apf::MeshEntity* down[3];
      m->getDownward(e, 0, down);
      apf::Vector3 p0, p1, p2;
      m->getPoint(down[0], 0, p0);
      m->getPoint(down[1], 0, p1);
      m->getPoint(down[2], 0, p2);
      double const lengthScale1 = (p1 - p0).getLength();
      double const lengthScale2 = (p2 - p1).getLength();
      double const lengthScale3 = (p0 - p2).getLength();
      double minlength = 1.e16;
      if (minlength > lengthScale1) minlength = lengthScale1;
      if (minlength > lengthScale2) minlength = lengthScale2;
      if (minlength > lengthScale3) minlength = lengthScale3;
      double const lengthScale = (lengthScale1+lengthScale2+lengthScale3)/3.;
      double const minlength_allowed = 0.0004; // measured min length 401 micron
      if (minlength < 0.0004) {
	printf("error: very small element min length %f, lengthScale %f\n", minlength, lengthScale);
	continue;
      }
      int non = 28; //order 6 BC
      //int non = 5; // for equidistant sampling
      //int non = getNumControlPoints(type, s_order);
      //int non = fs->countNodesOn(type);
      apf::Vector3 p, xi, pt(0,0,0), pt0_v;
      apf::DynamicVector pt0(dim);
      double haus_dist = -1.e32;
      apf::NewArray<double> nodeData;
      elem->getElementNodeData(nodeData);
      //printf("element node data size %d\n");

      //for(int i=0; i<non; ++i)

	//for (int d=0; d<dim; ++d) xi[d] = nodeData[i*dim+d]; // reference coords
      /*
      for (int j = 0; j <= non; ++j)
	xi[1] = 1.*j/non;
	for (int i = 0; i <= non-j; ++i)
	  xi[0] = 1.*i/non;
      */
      
      for(int i=0; i<non; ++i) {
	xi[0] = xi_p6tri[i*2+0];
	xi[1] = xi_p6tri[i*2+1];

	apf::getComponents(elem, xi, &pt0[0]);
	for (int d=0; d<dim; ++d) pt0_v[d] = pt0[d];
	m->getClosestPoint(gent,pt0_v,pt,p);

	double disp = 0.;
	for (int d=0; d<dim; ++d)
	  disp += (pt[d] - pt0[d])*(pt[d] - pt0[d]);
	disp = sqrt(disp);
	if ((isnan(disp))||(disp < 0)) {
	  printf("error computing disp : %f\n",disp);
	  continue;
	}
	//if (disp > lengthScale) 
	if (!m->isOnModel(gent, pt0_v, lengthScale, 1.0)) {
	  if (disp > lengthScale) {
	    printf("mesh face not on gface %d, disp %f, lscale %f\n",m->getModelTag(gent), disp, lengthScale);
	    printf("mesh pt {%f,%f,%f}, gpt {%f,%f,%f}\n",pt0[0], pt0[1], pt0[2], pt[0], pt[1], pt[2]);
	    printf("mesh pt xi {%f,%f,%f}\n",xi[0], xi[1], xi[2]);
	  }
	  printf("error, skipping this pt..., lengthScale %f\n", lengthScale);
	  bad_disp = 1;
	  break;
	}
	if (disp > haus_dist) haus_dist = disp;
	// curly bracket for double i,j for loop
      }
      if (bad_disp > 0) continue;
      if (haus_dist > max_dist) max_dist = haus_dist;
      mean_dist += haus_dist;

      apf::destroyElement(elem);
      apf::destroyMeshElement(mE);
      ++rc_faces;
    }
  }
  m->end(it);
  double mean_dist_g, max_dist_g = 0.;
  int rc_faces_g =0;
  MPI_Allreduce(&mean_dist, &mean_dist_g, 1, MPI_DOUBLE, MPI_SUM, PCU_Get_Comm());
  MPI_Allreduce(&rc_faces, &rc_faces_g, 1, MPI_INT, MPI_SUM, PCU_Get_Comm());
  mean_dist_g = mean_dist_g/(rc_faces_g*1.);
  MPI_Allreduce(&max_dist, &max_dist_g, 1, MPI_DOUBLE, MPI_MAX, PCU_Get_Comm());
  if (!PCU_Comm_Self()) {
    lion_eprint(1, "mean distance %f max distance %f nrc faces %d\n", mean_dist_g, 
      max_dist_g, rc_faces_g);
  }
  return;
}

}
// namespace crv

