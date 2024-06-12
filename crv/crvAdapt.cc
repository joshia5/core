/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvAdapt.h"
#include "crvShape.h"
#include <apf.h>
#include <apfMesh.h>
#include <maBalance.h>
#include <maCoarsen.h>
#include <maShape.h>
#include <maSnap.h>
#include <maStats.h>
#include <maLayer.h>
#include <PCU.h>
#include <pcu_util.h>

//for printing crvvtk
#include <reel.h>
#include <sys/types.h> 
#include <sys/stat.h> 
#include <errno.h> 
#include <algorithm>
#include <gmi_null.h>
#include <apfMDS.h>

namespace crv {

Adapt::Adapt(ma::Input* in)
: ma::Adapt(in)
{
  validityTag = mesh->createIntTag("crv_tags",1);
}

// rather than use the destructor to delete validityTag,
// this function takes care of it (since ~ma::Adapt() isn't virtual)
static void clearTags(Adapt* a)
{
  ma::Mesh* m = a->mesh;
  ma::Entity* e;
  for (int d=0; d <= 3; ++d)
  {
    ma::Iterator* it = m->begin(d);
    while ((e = m->iterate(it)))
      if (m->hasTag(e,a->validityTag))
        m->removeTag(e,a->validityTag);
    m->end(it);
  }
  m->destroyTag(a->validityTag);
}

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
    //if ((gdim == 2) && (gid == 3))
    if ((gdim == 2) && (r_min > 0.8) && (std::abs(coords[0]) < 0.5) && (gid != 3) && (gid != 13))
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

static int getTags(Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;
  if ( ! m->hasTag(e,a->validityTag))
    return 0; //we assume 0 is the default (unset) value for all tags
  int tags;
  m->getIntTag(e,a->validityTag,&tags);
  return tags;
}

static void setTags(Adapt* a, ma::Entity* e, int tags)
{
  a->mesh->setIntTag(e,a->validityTag,&tags);
}

void splitEdges(ma::Adapt* a)
{
  PCU_ALWAYS_ASSERT(ma::checkFlagConsistency(a,1,ma::SPLIT));
  //printf("split edges 1\n");
  ma::Refine* r = a->refine;
  //printf("split edges 2\n");
  ma::resetCollection(r);
  //printf("split edges 3\n");
  ma::collectForTransfer(r);
  printf("split edges 4\n");
  ma::addAllMarkedEdges(r);
  printf("split edges 5\n");
  ma::splitElements(r);
  printf("split edges 6\n");
  ma::processNewElements(r);
  printf("split edges 7\n");
  ma::destroySplitElements(r);
  printf("split edges 8\n");
  ma::forgetNewEntities(r);
  printf("split edges 9\n");
}

static void refine(ma::Adapt* a)
{
  double t0 = PCU_Time();
  --(a->refinesLeft);
  long count = ma::markEdgesToSplit(a);
  if ( ! count) {
    return;
  }
  splitEdges(a);
  double t1 = PCU_Time();
  ma::print("split %li edges in %f seconds",count,t1-t0);
}

int getValidityTag(ma::Mesh* m, ma::Entity* e,
    ma::Entity* bdry)
{
  if (bdry == e) return 18;
  m->getType(bdry);
  int dim = apf::getDimension(m,bdry);
  apf::Downward down;
  int n = m->getDownward(e,dim,down);
  int index = apf::findIn(down,n,bdry);
  // set up the tag here;
  switch (dim) {
    case 0:
      return index+2;
    case 1:
      return index+8;
    case 2:
      return index+14;
    default:
      fail("invalid lower entity in quality check\n");
      break;
  }
  return -1;
}

int markInvalidEntities(Adapt* a)
{
  ma::Entity* e;
  int count = 0;
  ma::Mesh* m = a->mesh;
  int dimension = m->getDimension();
  ma::Iterator* it = m->begin(dimension);
  Quality* qual = makeQuality(m,2);
  while ((e = m->iterate(it)))
  {
    /* this skip conditional is powerful: it affords us a
       3X speedup of the entire adaptation in some cases */
    int qualityTag = crv::getTag(a,e);
    if (qualityTag) continue;
    qualityTag = qual->checkValidity(e);
    if (qualityTag >= 2)
    {
      crv::setTag(a,e,qualityTag);
      if (m->isOwned(e))
        ++count;
    }
  }
  m->end(it);
  delete qual;
  return PCU_Add_Int(count);
}

int getTag(Adapt* a, ma::Entity* e)
{
  return getTags(a,e);
}

void setTag(Adapt* a, ma::Entity* e, int tag)
{
  setTags(a,e,tag);
}

void clearTag(Adapt* a, ma::Entity* e)
{
  setTags(a,e,0);
}
// use an identity configuration but with default fixing values
ma::Input* configureShapeCorrection(
    ma::Mesh* m, ma::SizeField* f,
    ma::SolutionTransfer* s)
{
  ma::Input* in = ma::makeAdvanced(ma::configureIdentity(m,f,s));
  in->shouldFixShape = true;
  in->shouldSnap = in->mesh->canSnap();
  in->shouldTransferParametric = in->mesh->canSnap();
  return in;
}

static int fixInvalidElements(crv::Adapt* a)
{
  a->input->shouldForceAdaptation = true;
  int count = crv::fixLargeBoundaryAngles(a)
            + crv::fixInvalidEdges(a);
  int originalCount = count;
  int prev_count;
  do {
    if ( ! count)
      break;
    prev_count = count;
    count = crv::fixLargeBoundaryAngles(a)
          + crv::fixInvalidEdges(a);
  } while(count < prev_count);

  crv::fixLargeBoundaryAngles(a);
  ma::clearFlagFromDimension(a,ma::COLLAPSE | ma::BAD_QUALITY,1);
  a->input->shouldForceAdaptation = false;
  return originalCount - count;
}

static void flagCleaner(crv::Adapt* a)
{
  int dim = a->mesh->getDimension();

  for (int d = 0; d <= dim; d++) {
    ma::clearFlagFromDimension(a, ma::BAD_QUALITY, d);
    ma::clearFlagFromDimension(a, ma::OK_QUALITY, d);
  }
}

void adapt(ma::Input* in)
{
  std::string name = in->mesh->getShape()->getName();
  if(name != std::string("Bezier"))
    fail("mesh must be bezier to adapt\n");

  in->shapeHandler = crv::getShapeHandler;
  ma::print("Curved Adaptation Version 2.0 ! : test modify for build");
  double t0 = PCU_Time();
  printf("ok1\n");
  ma::validateInput(in);
  printf("ok2\n");
  Adapt* a = new Adapt(in);
  printf("ok3\n");
  ma::preBalance(a);
  printf("ok4\n");

  if (in->shouldCoarsen) {
    //printf("fixing invalid elems\n");
    fixInvalidElements(a);
  }
  printf("ok4.5\n");

  for (int i=0; i < in->maximumIterations; ++i)
  {
    ma::print("iteration %d",i);
    if (in->shouldCoarsen) {
      ma::coarsen(a);
    }
    ma::midBalance(a);
    printf("ok5\n");
    //apf::synchronize(in->sizeField->sizes);
    crv::refine(a);
    if (in->shouldCoarsen) {
      allowSplitCollapseOutsideLayer(a);
    }
  printf("ok6\n");
    flagCleaner(a); // all true-flags must be false before using markEntities
    if (in->shouldFixShape) {
      fixCrvElementShapes(a);
    }
  }

  if (in->shouldCoarsen) {
    allowSplitCollapseOutsideLayer(a);
  }

  if (in->maximumIterations > 0 && in->shouldCoarsen) {
    fixInvalidElements(a);
    flagCleaner(a); // all true-flags must be false before using markEntities
    if (in->shouldFixShape) {
      fixCrvElementShapes(a);
    }
  }
  cleanupLayer(a);
  ma::printQuality(a);
  ma::postBalance(a);
  double t1 = PCU_Time();
  ma::print("mesh adapted in %f seconds",t1-t0);
  apf::printStats(a->mesh);
  crv::clearTags(a);
  
//makeSurfMesh(a->mesh, "108kp2uniref_vis", 15);

  delete a;
  // cleanup input object and associated sizefield and solutiontransfer objects
  if (in->ownsSizeField)
    delete in->sizeField;
  if (in->ownsSolutionTransfer)
    delete in->solutionTransfer;
  delete in;
}


void adapt(const ma::Input* in)
{
  crv::adapt(ma::makeAdvanced(in));
}

/** \brief Measures entity related quantities for a given mesh
  \details  quantities include normalized edge length, linear quality
  and curved quality. The values can be computed in both metric (if
  inMetric = true) and physical (if inMetric = false) spaces.*/
void stats(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities,
    std::vector<double> &curvedQualities,
    bool inMetric)
{
  ma::stats(m, sf, edgeLengths, linearQualities, inMetric);


  /* curved qualities are approximately the same in both
     metric and physical spaces
     metric quality = min(QJ) / max(QJ) ~ min(J) / max(J)
   */
  curvedQualities.clear();
  if (m->getShape()->getOrder() == 1)
    curvedQualities = std::vector<double>(linearQualities.size(), 0.0);
  else {
    crv::Quality* qual = makeQuality(m, 2);
    ma::Entity* e;
    ma::Iterator* it = m->begin(m->getDimension());
    while( (e = m->iterate(it)) ) {
      if (! m->isOwned(e))
	continue;
      if (! apf::isSimplex(m->getType(e))) // ignore non-simplex elements
        continue;
      curvedQualities.push_back(qual->getQuality(e));
    }
    m->end(it);
  }
}

}
