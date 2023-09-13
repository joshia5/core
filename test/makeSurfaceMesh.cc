#include <apf.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <PCU.h>
#include <pcu_util.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <crv.h>
#include <cassert>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
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
  m->writeNative(ss.str().c_str());
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

  printf("0\n");

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
    if ((gdim == 2) && (r_min > 0.8) && (std::abs(coords[0]) < 0.5) && (gid != 3) && (gid != 13))
      icavity2.push_back(e);
  }
  m->end(it);
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

int main(int argc, char** argv) {

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if (PCU_Comm_Peers() > 1) {
    printf("%s should only be used for serial (single part) meshes!\n", argv[0]);
    printf("use the serialize utility to get a serial mesh, and retry!\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  if (argc != 5) {
    printf("USAGE: %s <model> <mesh> <prefix> <resolution> \n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif

  gmi_register_null();

  const char* modelFile = argv[1];
  const char* meshFile  = argv[2];
  const char* prefix    = argv[3];
  int         res       = atoi(argv[4]);

  std::string gmi_native_path = argv[1];
  gmi_native_path = gmi_native_path.substr(0, gmi_native_path.length() - 4);
  gmi_native_path += "_nat.x_t";
  gmi_model* mdl;
  mdl = gmi_sim_load(gmi_native_path.c_str(),argv[1]);
 
  // load the mesh and check if the tag exists on the mesh
  apf::Mesh2* m = apf::loadMdsMesh(mdl,argv[2]);
  //apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);

  // make the root directory to save the cavity info
  int order = 3;
  crv::BezierCurver bc(m,order,0);
  bc.run();
  //m->changeShape(crv::getBezier(2), true);
  makeSurfMesh(m, prefix, res);

  // rest of the clean up
  m->destroyNative();
  apf::destroyMesh(m);

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif

  PCU_Comm_Free();
  MPI_Finalize();
}
