#include <crv.h>
#include <crvBezier.h>
#include <crvTables.h>
#include <crvSnap.h>
#include <crvMath.h>
#include <crvBezierShapes.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <PCU.h>
#include <lionPrint.h>
#include <mth.h>
#include <mth_def.h>
#include <pcu_util.h>
#include <ostream>

/* This file contains miscellaneous tests relating to bezier fitting
 */

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
 
  if ( argc != 3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];

  gmi_model* g = 0;
  g = gmi_load(modelFile);
  apf::Mesh2* m = 0;
  m = apf::loadMdsMesh(g, meshFile);

  // set the order
  apf::FieldShape* bezierShape = crv::getBezier(order);
  int non = bezierShape->getEntityShape(type)->countNodes();
  apf::Vector3 xi(xi0, xi1, xi2);
  apf::NewArray<double> vals(non);

  crv::setBlendingOrder(type, b);
  crv::BlendedTetGetValues(m,ent,xi,vals);

  PCU_Comm_Free();
  MPI_Finalize();
}
