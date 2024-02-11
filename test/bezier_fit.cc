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

#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#include <cstdlib>


/* This file contains miscellaneous tests relating to bezier fitting
 */

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
  gmi_register_mesh();
 
  if (argc != 5) {
    if (!PCU_Comm_Self())
      printf("Usage: %s <nat-model> <model> <mesh> order\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  gmi_model* g = 0;
  g = gmi_sim_load(argv[1], argv[2]);
  apf::Mesh2* m = 0;
  m = apf::loadMdsMesh(g, argv[3]);
  int order = atoi(argv[4]);
  crv::BezierCurver bc(m,order,0);
  bc.run();

  m->destroyNative();
  apf::destroyMesh(m);

  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit(); 
  PCU_Comm_Free();
  MPI_Finalize();
}
