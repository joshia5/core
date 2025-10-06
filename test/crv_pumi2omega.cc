#include <PCU.h>
#include <lionPrint.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include "SimParasolidKrnl.h"
#include <SimAdvMeshing.h>
#include <SimUtil.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <gmi.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cassert>
#include <getopt.h>
#include <string.h>
#include <stdio.h>

#include <crv.h>

#include <apfOmega_h.h>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_adapt.hpp>

using namespace std;

apf::Field* convert_my_tag(apf::Mesh* m, apf::MeshTag* t) {
  apf::MeshEntity* vtx;
  apf::MeshIterator* it = m->begin(0);
  apf::Field* f = apf::createFieldOn(m, "fathers2D_field", apf::SCALAR);
  int vals[1];
  double vals_d;
  while ((vtx = m->iterate(it))) {
    m->getIntTag(vtx, t, vals);
    vals_d = vals[0];
    apf::setScalar(f, vtx, 0, vals_d);
  }
  m->end(it);
  return f;
}

static void attachOrder(apf::Mesh* m) {
  apf::numberOverlapDimension(m, "sim_order", m->getDimension());
}

const char* gmi_path = NULL;
const char* gmi_native_path = NULL;
const char* sms_path = NULL;
const char* smb_path = NULL;
int should_log = 0;
int should_fix_pyramids = 1;
int should_attach_order = 0;
const char* extruRootPath = NULL;
int ExtruRootId =0;
bool found_bad_arg = false;

void getConfig(int argc, char** argv) {

  opterr = 0;

  static struct option long_opts[] = {
    {"no-pyramid-fix", no_argument, &should_fix_pyramids, 0},
    {"attach-order", no_argument, &should_attach_order, 1},
    {"enable-log", no_argument, &should_log, 2},
    {"model-face-root", required_argument, 0, 'e'},
    {"native-model", required_argument, 0, 'n'},
    {0, 0, 0, 0}  // terminate the option array
  };

  const char* usage=""
    "--native-model=/path/to/model <model file> <simmetrix mesh> <scorec mesh>\n";

  int option_index = 0;
  while(1) {
    int c = getopt_long(argc, argv, "", long_opts, &option_index);
    if (c == -1) break; //end of options
    switch (c) {
      case 0: // pyramid fix flag
      case 1: // attach order flag
      case 2: // enable simmetrix logging
        break;
      case 'e':
        extruRootPath = optarg;
        break;
      case 'n':
        gmi_native_path = optarg;
        break;
      case '?':
        if (!PCU_Comm_Self())
          printf ("warning: skipping unrecognized option \'%s\'\n", argv[optind-1]);
        break;
      default:
        if (!PCU_Comm_Self())
          printf("Usage %s %s", argv[0], usage);
        exit(EXIT_FAILURE);
    }
  }

  if(argc-optind != 3) {
    if (!PCU_Comm_Self())
      printf("Usage %s %s", argv[0], usage);
    exit(EXIT_FAILURE);
  }
  int i=optind;
  gmi_path = argv[i++];
  sms_path = argv[i++];
  smb_path = argv[i++];
  if (!PCU_Comm_Self()) {
    printf ("fix_pyramids %d attach_order %d enable_log %d extruRootPath %s\n",
            should_fix_pyramids, should_attach_order, should_log, extruRootPath);
    printf ("native-model \'%s\' model \'%s\' simmetrix mesh \'%s\' output mesh \'%s\'\n",
      gmi_native_path, gmi_path, sms_path, smb_path);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);

  getConfig(argc, argv);
  if( should_log )
    Sim_logOn("crv_pumi2omega.sim.log");

  if (should_attach_order && should_fix_pyramids) {
    if (!PCU_Comm_Self())
      std::cout << "disabling pyramid fix because --attach-order was given\n";
    should_fix_pyramids = false;
  }

  gmi_sim_start();
  gmi_register_sim();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  gmi_model* mdl;
  if( gmi_native_path ) {
    if (!PCU_Comm_Self())
      fprintf(stderr, "loading native model %s\n", gmi_native_path);
    mdl = gmi_sim_load(gmi_native_path,gmi_path);
  } else {
    mdl = gmi_load(gmi_path);
  }

  pGModel simModel = gmi_export_sim(mdl);
/*
  pParasolidNativeModel nModel = ParasolidNM_createFromFile(gmi_native_path,0);
  pGModel    Amodel = GAM_createFromNativeModel(nModel,progress); 
*/

  double t0 = PCU_Time();
  pParMesh sim_mesh = PM_load(sms_path, simModel, progress);
  double t1 = PCU_Time();
  if(!PCU_Comm_Self())
    fprintf(stderr, "read and created the simmetrix mesh in %f seconds\n", t1-t0);

  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);

  double t2 = PCU_Time();
  if(!PCU_Comm_Self())
    fprintf(stderr, "created the apf_sim mesh in %f seconds\n", t2-t1);
  if (should_attach_order) attachOrder(simApfMesh);

  apf::Mesh2* a_mesh = apf::createMdsMesh(mdl, simApfMesh);
  double t3 = PCU_Time();
  if(!PCU_Comm_Self())
    fprintf(stderr, "created the apf_mds mesh in %f seconds\n", t3-t2);

  apf::printStats(a_mesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  fixMatches(a_mesh);
  if (should_fix_pyramids) fixPyramids(a_mesh);
  a_mesh->verify();
  //a_mesh->writeNative(smb_path);

  int order = 3;
  //int order = atoi(argv[]);
  crv::BezierCurver bc(a_mesh,order,0);
  bc.run();

  auto o_lib = Omega_h::Library(&argc, &argv);
  Omega_h::Mesh o_mesh(&o_lib);
  apf::to_omega_h(&o_mesh, a_mesh);

  auto opts = Omega_h::AdaptOpts(&o_mesh);
  opts.verbosity = Omega_h::EXTRA_STATS;
  opts.length_histogram_max = 2.0;
  opts.max_length_allowed = opts.max_length_desired*2.0;
  opts.should_smooth_snap = 0;
  opts.should_coarsen = 0;
  opts.should_swap = 0;
  opts.should_coarsen_slivers = 0;
  opts.check_crv_qual = 0;
  opts.min_quality_allowed = 0.1;
  opts.min_quality_desired = 0.25;
  int desired_group_nelems = 2000;
  (approach_metric(&o_mesh, opts));
  int nelems = o_mesh.nelems();
  if (nelems < 8000) Omega_h::adapt(&o_mesh, opts);

  /*
  Omega_h::vtk::FullWriter writer;
  writer = Omega_h::vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus3d-24_crvsmb2osh.vtk",
      &o_mesh);
  writer.write();
  auto wireframe_mesh = Omega_h::Mesh(&o_lib);
  wireframe_mesh.set_comm(o_mesh.comm());
  Omega_h::build_cubic_wireframe_3d(&o_mesh, &wireframe_mesh, 10);
  std::string vtuPath =
    "/lore/joshia5/Meshes/curved/annulus3d-24-p2o_wire.vtu";
  Omega_h::vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Omega_h::Mesh(&o_lib);
  cubic_curveVtk_mesh.set_comm(o_mesh.comm());
  Omega_h::build_cubic_curveVtk_3d(&o_mesh, &cubic_curveVtk_mesh, 10);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus3d-24-p2o.vtu";
  Omega_h::vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  */

  a_mesh->destroyNative();
  apf::destroyMesh(a_mesh);

  apf::Mesh2* am2 = apf::makeEmptyMdsMesh(mdl, o_mesh.dim(), false);
  printf("ok0\n");
  apf::from_omega_h(am2, &o_mesh);
  //am2->writeNative(argv[3]);
  am2->destroyNative();
  apf::destroyMesh(am2);

  Progress_delete(progress);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
  if( should_log )
    Sim_logOff();
  PCU_Comm_Free();
  MPI_Finalize();
}
