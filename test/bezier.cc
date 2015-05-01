#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <PCU.h>
#include <ma.h>
#include <maCurveMesh.h>
#include <apfField.h>

void vert0(double const p[2], double x[3], void*)
{
  (void)p;
  (void)x;
}
// edges go counter clockwise
void edge0(double const p[2], double x[3], void*)
{
  x[0] = p[0];
  x[1] = 0.;
  x[2] = (2.*p[0]-1.)*(2.*p[0]-1.)-1.;
}
void edge1(double const p[2], double x[3], void*)
{
  x[0] = 1;
  x[1] = p[0];
  x[2] = 1.-(2.*p[0]-1.)*(2.*p[0]-1.);
}
void edge2(double const p[2], double x[3], void*)
{
  double u = 1.-p[0];
  x[0] = u;
  x[1] = 1.;
  x[2] = (2.*u-1.)*(2.*u-1.)-1.;
}
void edge3(double const p[2], double x[3], void*)
{
  double v = 1.-p[0];
  x[0] = 0;
  x[1] = v;
  x[2] = 1.-(2.*v-1.)*(2.*v-1.);
}
void face0(double const p[2], double x[3], void*)
{
  x[0] = p[0];
  x[1] = p[1];
  x[2] = (2.*p[0]-1.)*(2.*p[0]-1.)-(2.*p[1]-1.)*(2.*p[1]-1.);
}
void reparam_zero(double const from[2], double to[2], void*)
{
  (void)from;
  to[0] = 0;
  to[1] = 0;
}

void reparam_one(double const from[2], double to[2], void*)
{
  (void)from;
  to[0] = 1;
  to[1] = 0;
}
agm_bdry add_bdry(gmi_model* m, gmi_ent* e)
{
  return agm_add_bdry(gmi_analytic_topo(m), agm_from_gmi(e));
}

agm_use add_adj(gmi_model* m, agm_bdry b, int tag)
{
  agm* topo = gmi_analytic_topo(m);
  int dim = agm_dim_from_type(agm_bounds(topo, b).type);
  gmi_ent* de = gmi_find(m, dim - 1, tag);
  return agm_add_use(topo, b, agm_from_gmi(de));
}

void make_edge_topo(gmi_model* m, gmi_ent* e, int v0tag, int v1tag)
{
  agm_bdry b = add_bdry(m, e);
  agm_use u0 = add_adj(m, b, v0tag);
  gmi_add_analytic_reparam(m, u0, reparam_zero, 0);
  agm_use u1 = add_adj(m, b, v1tag);
  gmi_add_analytic_reparam(m, u1, reparam_one, 0);
}


gmi_model* makeModel()
{
  gmi_model* model = gmi_make_analytic();
  int edPer = 0;
  double edRan[2] = {0, 1};
  for(int i = 0; i < 4; ++i)
    gmi_add_analytic(model, 0, i, vert0, NULL,NULL,NULL);
  gmi_ent* eds[4];
  eds[0] = gmi_add_analytic(model, 1, 0, edge0, &edPer, &edRan, 0);
  eds[1] = gmi_add_analytic(model, 1, 1, edge1, &edPer, &edRan, 0);
  eds[2] = gmi_add_analytic(model, 1, 2, edge2, &edPer, &edRan, 0);
  eds[3] = gmi_add_analytic(model, 1, 3, edge3, &edPer, &edRan, 0);
  for(int i = 0; i < 4; ++i)
    make_edge_topo(model, eds[i], i, (i+1) % 4);
  int faPer[2] = {0, 0};
  double faRan[2][2] = {{0,1},{0,1}};
  gmi_add_analytic(model, 2, 0, face0, faPer, faRan, 0);

  return model;
}


apf::Mesh2* createMesh2D()
{
  gmi_model* model = makeModel();
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, true);
  apf::MeshEntity* v[4];
  apf::Vector3 points2D[4] =
  {apf::Vector3(0,0,0),
   apf::Vector3(1,0,0),
   apf::Vector3(1,1,0),
   apf::Vector3(0,1,0)};

  for (int i = 0; i < 4; ++i){
    v[i] = m->createVertex(m->findModelEntity(0,i),points2D[i],points2D[i]);
  }
  for (int i = 0; i < 4; ++i){
    apf::ModelEntity* edge = m->findModelEntity(1,i);
    apf::MeshEntity* ved[2] = {v[i],v[(i+1) % 4]};
    apf::buildElement(m, edge, apf::Mesh::EDGE, ved);
  }

  apf::MeshEntity* ved[2] = {v[0],v[2]};
  apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::EDGE, ved);

  apf::MeshEntity* vf0[3] = {v[0],v[1],v[2]};
  apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, vf0);
  apf::MeshEntity* vf1[3] = {v[0],v[2],v[3]};
  apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, vf1);
  m->acceptChanges();
  m->verify();
  return m;
}
apf::Mesh2* createMesh3D()
{
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, true);
  apf::MeshEntity* v[4];
  apf::Vector3 points3D[4] =
  {apf::Vector3(0,0,0),
   apf::Vector3(1,0,0),
   apf::Vector3(0,1,0),
   apf::Vector3(0,0,1)};

  int const (*ttv)[3] = apf::tet_tri_verts;
  int const (*tev)[2] = apf::tet_edge_verts;

  for (int i = 0; i < 4; ++i){
    v[i] = m->createVertex(0,points3D[i],points3D[i]);
  }
  for (int i = 0; i < 6; ++i){
    apf::MeshEntity* ved[2] = {v[tev[i][0]],v[tev[i][1]]};
    apf::buildElement(m, 0, apf::Mesh::EDGE, ved);
  }
  apf::MeshEntity* ved0[2] = {v[1],v[3]};
  apf::buildElement(m, 0, apf::Mesh::EDGE, ved0);
  apf::MeshEntity* ved1[2] = {v[0],v[2]};
  apf::buildElement(m, 0, apf::Mesh::EDGE, ved1);

  for(int i = 0; i < 4; ++i){
    // 0,1,3 and 1,2,3
    apf::MeshEntity* vf[3] = {v[ttv[i][0]],v[ttv[i][1]],v[ttv[i][2]]};
    if(i == 1 || i == 2)
      apf::buildElement(m, 0, apf::Mesh::TRIANGLE, vf);
    else
      apf::buildElement(m, 0, apf::Mesh::TRIANGLE, vf);

  }
  apf::buildElement(m,0,apf::Mesh::TET,v);
  apf::deriveMdsModel(m);
  m->acceptChanges();
  m->verify();
  return m;
}
void testInterpolatedPoints2D(apf::Mesh2* m){
  apf::FieldShape * fs = m->getCoordinateField()->getShape();
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    apf::ModelEntity* g = m->toModel(e);
    if (m->getModelType(g) == m->getDimension())
      continue;
    ma::Vector pt,pa(0.,0.,0.),cpt,cpa;
    apf::Element* elem =
        apf::createElement(m->getCoordinateField(),e);

    for(int i = 0; i < fs->countNodesOn(1); ++i){
      fs->getNodeXi(1,i,pa);
      apf::getVector(elem,pa,pt);
      pa[0] = 0.5*(pa[0]+1.);
      m->snapToModel(g,pa,cpt);
      double error = (cpt-pt).getLength();
      if (error > 1.e-8) {
        std::cout << apf::Mesh::typeName[m->getType(e)] <<
          " is not being interpolated correctly by " <<
          m->getCoordinateField()->getShape()->getName() <<
          " with error " << error << "\n";
        abort();
      }
    }
    apf::destroyElement(elem);
  }
  m->end(it);
}
void testSize2D(apf::Mesh2* m)
{
  double sizes[2] = {2.323391881216468,1.625569984255547};
  int order = m->getCoordinateField()->getShape()->getOrder();
  for(int d = 1; d <= 2; ++d){
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::ModelEntity* g = m->toModel(e);
      apf::MeshElement* me = apf::createMeshElement(m,e);
      double v = apf::measure(me);
      // these are checks for the middle edge and first order
      if (m->getModelType(g) == 2 && d == 1){
        assert(fabs(v-sqrt(2))< 1e-10);
      } else if(order == 1 && d == 1){
        assert(fabs(v-1.0 )< 1e-10);
      } else if(order == 1 && d == 2){
        assert(fabs(v-0.5)< 1e-10);
      } else if(fabs(v-sizes[d-1])/sizes[d-1] > (1.-0.5*exp(-order))){
        std::stringstream ss;
        ss << "error: " << apf::Mesh::typeName[m->getType(e)]
           << " size " << v
           << " at " << getLinearCentroid(m, e) << '\n';
        std::string s = ss.str();
        fprintf(stderr, "%s", s.c_str());
        abort();
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }
}
void test2D()
{
  for(int order = 1; order < 7; ++order){
    apf::Mesh2* m = createMesh2D();
    ma::curveMeshToBezier(m,order);
    testInterpolatedPoints2D(m);
    testSize2D(m);
    m->destroyNative();
    apf::destroyMesh(m);
  }
}
void testSize3D(apf::Mesh2* m)
{
  int order = m->getCoordinateField()->getShape()->getOrder();
  printf("order %d\n",order);
  for(int d = 1; d <= 2; ++d){
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
//      apf::ModelEntity* g = m->toModel(e);
      apf::MeshElement* me = apf::createMeshElement(m,e);
      double v = apf::measure(me);
//      printf("element type %s has size %f\n",
//          apf::Mesh::typeName[m->getType(e)],v);
      if((d == 1 && !(fabs(v-1.) < 1e-8 || fabs(v-sqrt(2)) < 1e-8 ))
          || (d == 2 && !(fabs(v-0.5) < 1e-8 || fabs(v-sqrt(3)/2.) < 1e-8))){
        std::stringstream ss;
        ss << "error: " << apf::Mesh::typeName[m->getType(e)]
           << " size " << v
           << " at " << getLinearCentroid(m, e) << '\n';
        std::string s = ss.str();
        fprintf(stderr, "%s", s.c_str());
//        abort();
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }
}
void test3D()
{
  for(int order = 1; order < 7; ++order){

    gmi_register_null();

    apf::Mesh2* m = createMesh3D();
    apf::changeMeshShape(m, apf::getBezier(3,order),true);

    apf::FieldShape * fs = m->getCoordinateField()->getShape();

    apf::MeshEntity* e;

    // go downward, and convert interpolating to control points
    for(int d = 2; d >= 1; --d){
      int n = (d == 2)? (order+1)*(order+2)/2 : order+1;
      int ne = fs->countNodesOn(d);

      apf::NewArray<double> c;
      apf::getTransformationCoefficients(order,3,d,c);

      apf::MeshIterator* it = m->begin(d);
      while ((e = m->iterate(it))) {
        apf::ModelEntity* g = m->toModel(e);
        if(m->getModelType(g) == m->getDimension()) continue;
        ma::convertInterpToBezier(m,e,n,ne,c);
      }
      m->end(it);
    }

    m->acceptChanges();
    m->verify();
//    testSize3D(m);
    m->destroyNative();
    apf::destroyMesh(m);
  }
}
int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  test2D();
  test3D();
  PCU_Comm_Free();
  MPI_Finalize();
}