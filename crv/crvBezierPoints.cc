/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvBezier.h"
#include "crvTables.h"

namespace crv {

void getBezierCurveNodeXi(int type, int P, int node, apf::Vector3& xi)
{
  static double eP2[1] = {0.0};
  static double eP3[2] = {-0.4306648,0.4306648};
  static double eP4[3] = {-0.6363260,0.0,0.6363260};
  static double eP5[4] = {-0.7485748,-0.2765187,0.2765187,0.7485748};
  static double eP6[5] = {-0.8161268,-0.4568660,0.0,
      0.4568660,0.8161268};

  static double* edgePoints[6] =
  {eP2, eP2, eP3, eP4, eP5, eP6 };

  if(type == apf::Mesh::EDGE && P > 1){
    xi[0] = edgePoints[P-1][node];
  } else {
    getBezierNodeXi(type,P,node,xi);
  }
}

void getBezierNodeXi(int type, int P, int node, apf::Vector3& xi)
{
  static double eP2[1] = {0.0};
  static double eP3[2] = {-0.4503914,0.4503914};
  static double eP4[3] = {-0.6612048,0.0,0.6612048};
  static double eP5[4] = {-0.7732854,-0.2863522,0.2863522,0.7732854};
  static double eP6[5] = {-0.8388042,-0.469821,0.0,
      0.469821,0.8388042};
  static double* edgePoints[6] =
  {eP2, eP2, eP3, eP4, eP5, eP6 };
  switch (type) {
  case apf::Mesh::EDGE:
    xi[0] = edgePoints[P-1][node];
    break;
  case apf::Mesh::TRIANGLE:
    xi = apf::Vector3(1./3.,1./3.,1./3.);
    if(node == curved_face_internal[BEZIER][P-1]-1 && P % 3 == 0){
      return;
    } else { // technically only two of these numbers are needed
      switch (P) {
      case 1:
      case 2:
      case 3:
        fail("expected P >= 4");
      case 4:
        xi[(node+2) % 3] = 0.5582239;
        xi[(node+0) % 3] = 0.22088805;
        xi[(node+1) % 3] = 0.22088805;
        break;
      case 5:
        if(node % 2 == 0) {
          xi[(node/2+2) % 3] = 0.6949657;
          xi[(node/2+0) % 3] = 0.15251715;
          xi[(node/2+1) % 3] = 0.15251715;
        } else {
          xi[((node-1)/2+2) % 3] = 0.4168658;
          xi[((node-1)/2+0) % 3] = 0.4168658;
          xi[((node-1)/2+1) % 3] = 0.1662684;
        }
        break;
      case 6:
        if (node % 3 == 0) {
          xi[(node/3+2) % 3] = 0.7805723;
          xi[(node/3+0) % 3] = 0.10971385;
          xi[(node/3+1) % 3] = 0.10971385;
        } else if ((node-1) % 3 == 0) {
          xi[((node-1)/3+2) % 3] = 0.5586077;
          xi[((node-1)/3+0) % 3] = 0.3157892;
          xi[((node-1)/3+1) % 3] = 0.1256031;
        } else if ((node-2) % 3 == 0) {
          xi[((node-2)/3+2) % 3] = 0.3157892;
          xi[((node-2)/3+0) % 3] = 0.5586077;
          xi[((node-2)/3+1) % 3] = 0.1256031;
        }
        break;
      }
    }
    break;
  case apf::Mesh::TET:
    switch (P) {
    case 1:
    case 2:
    case 3:
      fail("expected P == 4");
    case 4:
      xi = apf::Vector3(0.25,0.25,0.25);
      break;
    case 5:
    case 6:
      fail("expected P == 4");
    }
    break;
    default:
      xi.zero();
      break;
  }
}

static void getBezierCurveTransform(int P, apf::NewArray<double> & c)
{
  double e2[3] = {-0.5,-0.5,2};
  double e3[8] = {
      -0.970273514083553,0.333333333333333,2.71895067382449,-1.08201049307427,
      0.333333333333333,-0.970273514083553,-1.08201049307427,2.71895067382449};
  double e4[15] = {
      -1.4304202857228,-0.25,3.39545839723405,-1.46967987431139,
      0.754641762800137,0.953613523815196,0.953613523815197,-2.76673344002279,4.62623983241519,
      -2.76673344002279,-0.25,-1.4304202857228,0.754641762800137,-1.46967987431139,
      3.39545839723405};
  double e5[24] = {
      -1.88592269024942,0.2,4.05614415979432,-1.81653638123435,
      1.0295423816296,-0.583227469940158,1.85476912333284,-0.942961345124708,-5.01939997635205,6.96205913930752,
      -4.56234099538677,2.70787405422317,-0.942961345124708,1.85476912333285,2.70787405422317,-4.56234099538677,
      6.96205913930752,-5.01939997635206,0.2,-1.88592269024942,-0.583227469940158,1.0295423816296,
      -1.81653638123435,4.05614415979432};
  double e6[35] = {
      -2.33890800235808,-0.166666666666667,4.70907763497668,-2.14695478588352,
      1.2670886004356,-0.80040589915343,0.476769118649422,3.03457283388393,
      0.935563200943235,-7.82909978199834,9.74813267975089,
      -6.60581336123903,4.37362214799981,-2.65697771934049,-2.27592962541295,
      -2.27592962541295,6.3088040999163,-9.70710791530195,
      12.3484668815972,-9.70710791530194,6.3088040999163,0.935563200943235,
      3.03457283388393,-2.65697771934049,4.37362214799981,
      -6.60581336123903,9.74813267975088,-7.82909978199834,-0.166666666666667,
      -2.33890800235809,0.476769118649422,-0.80040589915343,
      1.2670886004356,-2.14695478588352,4.70907763497668};
  double* table[5] = {
      e2,e3,e4,e5,e6};
  int nb = P-1;
  int ni = P+1;
  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = table[P-2][i*ni+j];

}
static void getBezierEdgeTransform(int P, apf::NewArray<double> & c)
{
  double e2[3] = {-0.5,-0.5,2};
  double e3[8] = {
      -1.00596379148431,0.333333333333333,2.69317845753742,-1.02054799938644,
      0.333333333333333,-1.00596379148431,-1.02054799938644,2.69317845753742};
  double e4[15] = {
      -1.52680420766155,-0.25,3.37567603341243,-1.28732567375034,
      0.688453847999458,1.01786947177436,1.01786947177436,-2.70941992094126,4.38310089833379,
      -2.70941992094126,-0.25,-1.52680420766155,0.688453847999459,-1.28732567375034,
      3.37567603341243};
  double e5[24] = {
      -2.06136018481524,0.2,4.05936188849008,-1.52513018373641,
      0.846118038541145,-0.518989558479574,2.07392858296555,-1.03068009240762,
      -5.04403528329688,6.4850808350761,-4.1256786562572,2.64138461392004,
      -1.03068009240762,2.07392858296555,2.64138461392004,-4.1256786562572,
      6.48508083507611,-5.04403528329688,0.2,-2.06136018481524,-0.518989558479574,
      0.846118038541145,-1.52513018373641,4.05936188849008};
  double e6[35] = {
      -2.60465921875445,-0.166666666666667,4.74317259573776,-1.74847211573074,
      0.99151406014263,-0.630691218757935,0.415802564029398,3.50945493261743,
      1.04186368750178,-8.00777336834505,8.99109434126308,
      -5.80152934540333,3.85156704410589,-2.5846772917398,-2.63209119946307,
      -2.63209119946307,6.39664544713349,-8.91824703868012,
      11.3073855820194,-8.91824703868012,6.39664544713349,1.04186368750178,
      3.50945493261743,-2.5846772917398,3.85156704410589,
      -5.80152934540333,8.99109434126308,-8.00777336834505,-0.166666666666667,
      -2.60465921875445,0.415802564029398,-0.630691218757935,
      0.991514060142631,-1.74847211573074,4.74317259573777};

  double* table[5] = {
      e2,e3,e4,e5,e6};

  int nb = P-1;
  int ni = P+1;
  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = table[P-2][i*ni+j];
}

static void getBezierTriangleTransform(int P, apf::NewArray<double> & c)
{
  double f3[10] = {
      0.5059637914843115,0.5059637914843116,0.5059637914843119,-0.8363152290754889,
      -0.8363152290754891,-0.8363152290754891,-0.8363152290754897,-0.8363152290754891,
      -0.8363152290754898,4.5};
  double f4[45] = {
      1.473866405971784,-0.4873245458152482,-0.4873245458152483,-2.157170326583087,
      -0.7895371825786641,1.002268609355065,0.4569922360188257,0.4160673296503333,
      0.4569922360188257,1.002268609355065,-0.7895371825786643,-2.157170326583087,
      7.066481928037422,-2.00343662222666,-2.00343662222666,
      -0.4873245458152483,1.473866405971784,-0.4873245458152482,1.002268609355065,
      -0.7895371825786646,-2.157170326583087,-2.157170326583087,-0.7895371825786643,
      1.002268609355065,0.4569922360188255,0.4160673296503333,0.4569922360188255,
      -2.00343662222666,7.066481928037421,-2.00343662222666,
      -0.4873245458152483,-0.4873245458152481,1.473866405971784,0.4569922360188258,
      0.4160673296503332,0.4569922360188255,1.002268609355065,-0.7895371825786646,
      -2.157170326583087,-2.157170326583088,-0.7895371825786643,1.002268609355065,
      -2.00343662222666,-2.00343662222666,7.066481928037422};
  double f5[126] = {
      2.955509375394112,0.4850633218816851,0.4850633218816851,-4.015153971419451,
      -0.6339185304220427,1.140973028173454,-1.041278651590119,-0.3192965376118003,
      -0.2178043751582093,-0.2178043751582093,-0.3192965376118008,-1.041278651590118,
      1.140973028173453,-0.633918530422042,-4.015153971419451,10.16896081482358,
      -3.086925777444234,1.177507607441262,0.8971975820812148,1.177507607441263,
      -3.086925777444236,
      -1.87982956469617,-1.879829564696171,0.5695483481139566,3.558924203285685,
      -2.437779912677446,-2.437779912677446,3.558924203285686,1.772590081128005,
      0.9751122788728792,0.1149068856949735,-0.833741865064019,-0.8337418650640189,
      0.1149068856949739,0.9751122788728801,1.772590081128005,-6.137619043989807,
      12.47852594090006,-6.137619043989807,-2.280366067778754,2.247531721435293,
      -2.280366067778758,
      0.4850633218816852,2.955509375394111,0.4850633218816852,-1.041278651590119,
      1.140973028173453,-0.6339185304220422,-4.01515397141945,-4.015153971419449,
      -0.6339185304220425,1.140973028173453,-1.041278651590118,-0.3192965376118005,
      -0.2178043751582091,-0.2178043751582096,-0.3192965376118004,1.177507607441262,
      -3.086925777444234,10.16896081482358,-3.086925777444233,1.177507607441262,
      0.8971975820812153,
      0.569548348113957,-1.87982956469617,-1.879829564696171,-0.8337418650640194,
      0.114906885694974,0.97511227887288,1.772590081128004,3.558924203285685,
      -2.437779912677446,-2.437779912677443,3.558924203285685,1.772590081128006,
      0.9751122788728793,0.1149068856949746,-0.8337418650640199,2.247531721435295,
      -2.280366067778758,-6.137619043989805,12.47852594090005,-6.137619043989806,
      -2.280366067778759,
      0.4850633218816858,0.4850633218816853,2.955509375394112,-0.3192965376118007,
      -0.2178043751582097,-0.2178043751582098,-0.3192965376118001,-1.041278651590119,
      1.140973028173454,-0.6339185304220429,-4.015153971419452,-4.01515397141945,
      -0.6339185304220436,1.140973028173455,-1.04127865159012,1.177507607441263,
      0.8971975820812174,1.177507607441262,-3.086925777444234,10.16896081482359,
      -3.086925777444237,
      -1.879829564696172,0.5695483481139572,-1.879829564696172,1.772590081128005,
      0.9751122788728803,0.114906885694975,-0.8337418650640201,-0.8337418650640197,
      0.114906885694974,0.9751122788728801,1.772590081128007,3.558924203285685,
      -2.437779912677443,-2.437779912677449,3.558924203285688,-6.137619043989808,
      -2.280366067778761,2.247531721435295,-2.280366067778758,-6.137619043989812,
      12.4785259409000};
  double f6[280] = {
      4.990795106388393,-0.4798837984147291,-0.4798837984147287,-6.458245423578343,
      -0.3639631629214212,1.22336563258444,-1.237281132017878,1.048637951819922,
      0.2399023612066681,0.1476405704559505,0.09371675625975406,0.1476405704559514,
      0.2399023612066667,1.048637951819923,-1.23728113201788,1.223365632584442,
      -0.3639631629214235,-6.458245423578343,13.89283735202935,-4.276280903149708,
      1.808916611395355,-0.8152389231864423,-0.4900563063939803,-0.4900563063939829,
      -0.8152389231864404,1.808916611395354,-4.276280903149705,1.327623829722839,
      -4.689458981957318,2.303409358216713,-0.6686820957216264,8.401259643873185,
      -5.327640035897593,-2.380733253738293,4.596086431435174,-4.67688716273589,
      -1.613859854887045,-0.8313011663168295,-0.3532695050699889,0.01692167342423888,
      0.7838078136607817,1.115432125134806,-0.6396258416621832,-0.4063403163660048,
      1.633507742829366,4.420179509817743,-13.10100871536221,20.06968509845678,
      -10.70127780752898,5.23777942309397,2.406990467529688,0.7247013815918657,
      -2.116258714513936,3.071067981878374,-2.093285349716039,-4.18119984946874,
      2.303409358216714,-4.689458981957322,-0.6686820957216258,-4.676887162735897,
      4.596086431435186,-2.380733253738301,-5.327640035897589,8.401259643873187,
      4.420179509817745,1.633507742829365,-0.4063403163660045,-0.6396258416621829,
      1.115432125134805,0.7838078136607807,0.01692167342423913,-0.353269505069989,
      -0.8313011663168312,-1.613859854887043,5.23777942309397,-10.70127780752899,
      20.06968509845678,-13.10100871536222,-2.093285349716035,3.071067981878372,
      -2.116258714513934,0.7247013815918641,2.406990467529694,-4.181199849468742,
      -0.4798837984147298,4.990795106388394,-0.4798837984147283,1.048637951819924,
      -1.237281132017881,1.223365632584441,-0.3639631629214217,-6.458245423578346,
      -6.458245423578344,-0.3639631629214246,1.223365632584442,-1.237281132017879,
      1.048637951819922,0.2399023612066666,0.1476405704559507,0.09371675625975429,
      0.1476405704559513,0.2399023612066678,-0.8152389231864428,1.808916611395358,
      -4.27628090314971,13.89283735202935,-4.276280903149705,1.808916611395352,
      -0.8152389231864391,-0.4900563063939816,-0.490056306393983,1.32762382972284,
      -0.6686820957216246,-4.689458981957318,2.303409358216711,1.115432125134803,
      -0.6396258416621817,-0.4063403163660055,1.633507742829363,4.420179509817745,
      8.40125964387318,-5.32764003589758,-2.3807332537383,4.596086431435179,
      -4.676887162735889,-1.613859854887039,-0.8313011663168305,-0.3532695050699883,
      0.016921673424238,0.78380781366078,-2.116258714513931,3.071067981878368,
      -2.093285349716028,-13.10100871536221,20.06968509845676,-10.70127780752896,
      5.23777942309396,2.406990467529691,0.724701381591865,-4.181199849468744,
      -0.6686820957216236,2.30340935821671,-4.689458981957321,0.7838078136607789,
      0.01692167342423904,-0.353269505069988,-0.8313011663168306,-1.61385985488704,
      -4.676887162735886,4.596086431435173,-2.380733253738294,-5.327640035897591,
      8.401259643873187,4.420179509817743,1.633507742829364,-0.4063403163660066,
      -0.6396258416621801,1.115432125134801,-2.116258714513929,0.7247013815918617,
      2.406990467529691,5.237779423093964,-10.70127780752897,20.06968509845676,
      -13.10100871536221,-2.093285349716029,3.071067981878367,-4.18119984946874,
      -0.4798837984147294,-0.4798837984147293,4.990795106388396,0.2399023612066681,
      0.1476405704559513,0.09371675625975392,0.147640570455951,0.2399023612066676,
      1.048637951819923,-1.23728113201788,1.223365632584443,-0.3639631629214256,
      -6.458245423578345,-6.458245423578348,-0.3639631629214194,1.223365632584439,
      -1.237281132017879,1.048637951819923,-0.8152389231864435,-0.4900563063939818,
      -0.4900563063939812,-0.8152389231864419,1.808916611395355,-4.276280903149707,
      13.89283735202936,-4.27628090314971,1.808916611395358,1.327623829722839,
      2.303409358216711,-0.6686820957216234,-4.689458981957324,-1.613859854887043,
      -0.8313011663168315,-0.353269505069988,0.01692167342423824,0.7838078136607798,
      1.115432125134801,-0.639625841662179,-0.4063403163660082,1.633507742829368,
      4.420179509817745,8.401259643873196,-5.327640035897597,-2.380733253738291,
      4.596086431435175,-4.676887162735889,5.23777942309397,2.406990467529692,
      0.7247013815918628,-2.11625871451393,3.071067981878368,-2.093285349716032,
      -13.10100871536222,20.06968509845678,-10.70127780752898,-4.181199849468741,
      -4.689458981957318,-0.6686820957216246,2.303409358216713,4.420179509817742,
      1.633507742829366,-0.4063403163660054,-0.6396258416621825,1.115432125134803,
      0.7838078136607793,0.01692167342423931,-0.3532695050699873,-0.831301166316834,
      -1.61385985488704,-4.676887162735893,4.596086431435182,-2.380733253738301,
      -5.327640035897584,8.401259643873182,-13.10100871536221,-2.093285349716037,
      3.071067981878372,-2.11625871451393,0.7247013815918595,2.406990467529697,
      5.237779423093967,-10.70127780752898,20.06968509845677,-4.181199849468742,
      2.740410900541109,2.740410900541107,2.740410900541116,-3.896719679725134,
      0.8525687056196629,2.629196981978395,0.8525687056196668,-3.896719679725134,
      -3.89671967972513,0.8525687056196584,2.6291969819784,0.8525687056196664,
      -3.896719679725144,-3.896719679725143,0.8525687056196644,2.629196981978399,
      0.8525687056196618,-3.896719679725134,9.218530840490745,-7.999447648758338,
      -7.999447648758347,9.218530840490741,-7.99944764875833,-7.999447648758361,
      9.218530840490766,-7.999447648758353,-7.99944764875834,23.49717556815212
  };
  double* table[4] = {f3,f4,f5,f6};
  int nb = curved_face_internal[BEZIER][P-1];
  int ni = curved_face_total[BEZIER][P-1];
  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = table[P-3][i*ni+j];
}

static void getBezierTetTransform(int P, apf::NewArray<double> & c)
{
  assert(P == 4 && getBlendingOrder() == 0);
  double t4[35] = {
      -0.665492638178598,-0.665492638178598,-0.665492638178598,-0.665492638178598,
      0.697909481209196,0.496340368840329,0.697909481209197,0.697909481209196,
      0.496340368840329,0.697909481209196,0.697909481209196,0.49634036884033,
      0.697909481209196,0.697909481209196,0.496340368840329,0.697909481209196,
      0.697909481209196,0.496340368840329,0.697909481209196,0.697909481209196,
      0.496340368840329,0.697909481209196,-1.52980434179205,-1.52980434179205,
      -1.52980434179205,-1.52980434179205,-1.52980434179205,-1.52980434179205,
      -1.52980434179205,-1.52980434179205,-1.52980434179205,-1.52980434179205,
      -1.52980434179205,-1.52980434179205,10.6666666666667,
  };

  int nb = curved_tet_internal[BEZIER][P-1];
  int ni = curved_tet_total[BEZIER][P-1];

  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = t4[i*ni+j];
}
void getTransformationCoefficients(int dim, int P, int type,
    apf::NewArray<double>& c){
  if(dim == 2 && getBlendingOrder() > 0)
    getBezierCurveTransform(P,c);
  else if(type == apf::Mesh::EDGE)
    getBezierEdgeTransform(P,c);
  else if(type == apf::Mesh::TRIANGLE)
    getBezierTriangleTransform(P,c);
  else if(type == apf::Mesh::TET)
    getBezierTetTransform(P,c);
}

static void getGregoryTriangleTransform(int P, apf::NewArray<double> & c)
{
  apf::NewArray<double> d;
  getBezierTriangleTransform(P,d);

  int nbBezier = curved_face_internal[BEZIER][P-1];
  int niBezier = curved_face_total[BEZIER][P-1];

  int nb = curved_face_internal[GREGORY][P-1];
  int ni = curved_face_total[GREGORY][P-1];
  c.allocate(ni*nb);

  int map[3] = {1,2,0};
  // copy the bezier point locations
  for(int i = 0; i < nbBezier; ++i){
    for(int j = 0; j < niBezier; ++j)
      c[i*ni+j] = d[i*niBezier+j];
    for(int j = niBezier; j < ni; ++j)
      c[i*ni+j] = 0.;
  }
  if(P == 3){
    for(int i = nbBezier; i < nb; ++i){
      for(int j = 0; j < niBezier; ++j)
        c[i*ni+j] = d[j];
      for(int j = niBezier; j < ni; ++j)
        c[i*ni+j] = 0.;
    }
  }
  if(P == 4){
    for(int i = nbBezier; i < nb; ++i){
      for(int j = 0; j < niBezier; ++j)
        c[i*ni+j] = d[map[i-nbBezier]*niBezier+j];
      for(int j = niBezier; j < ni; ++j)
        c[i*ni+j] = 0.;
    }
  }
}

static void getGregoryTetTransform(int P, apf::NewArray<double> & c)
{
  assert(P == 4 && getBlendingOrder() == 0);
  double t4[47] = {
      -0.665492638178598,-0.665492638178598,-0.665492638178598,-0.665492638178598,
      0.697909481209196,0.496340368840329,0.697909481209197,0.697909481209196,
      0.496340368840329,0.697909481209196,0.697909481209196,0.49634036884033,
      0.697909481209196,0.697909481209196,0.496340368840329,0.697909481209196,
      0.697909481209196,0.496340368840329,0.697909481209196,0.697909481209196,
      0.496340368840329,0.697909481209196,
      -1.52980434179205,-1.52980434179205,-1.52980434179205,0.,0.,0.,
      -1.52980434179205,-1.52980434179205,-1.52980434179205,0.,0.,0.,
      -1.52980434179205,-1.52980434179205,-1.52980434179205,0.,0.,0.,
      -1.52980434179205,-1.52980434179205,-1.52980434179205,0.,0.,0.,
      10.6666666666667,
  };

  int nb = curved_tet_internal[GREGORY][P-1];
  int ni = curved_tet_total[GREGORY][P-1];

  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = t4[i*ni+j];
}

void getGregoryTransformationCoefficients(int /*dim*/, int P, int type,
    apf::NewArray<double>& c){
  assert(P == 3 || P == 4);
  if(type == apf::Mesh::EDGE)
    getBezierEdgeTransform(P,c);
  else if(type == apf::Mesh::TRIANGLE)
    getGregoryTriangleTransform(P,c);
  else if(type == apf::Mesh::TET)
    getGregoryTetTransform(P,c);
}

void getTransformationMatrix(apf::Mesh* m, int type, apf::DynamicMatrix& A)
{
  setBlendingOrder(0); // makes sure blending is turned off.
  apf::FieldShape* fs = m->getShape();
  apf::EntityShape* es = fs->getEntityShape(type);
  int n = es->countNodes();
  int typeDim = apf::Mesh::typeDimension[type];
  apf::Downward down;
  apf::Vector3 xi, exi;
  apf::NewArray<double> values;

  apf::MeshIterator* it = m->begin(typeDim);
  apf::MeshEntity* e = m->iterate(it); // take the first one
  m->end(it);

  A.setSize(n,n);

  int row = 0;
  for(int d = 0; d <= typeDim; ++d){
    int nDown = m->getDownward(e,d,down);
    for(int j = 0; j < nDown; ++j){
      int bt = m->getType(down[j]);
      for(int x = 0; x < fs->countNodesOn(bt); ++x){
        fs->getNodeXi(bt,x,xi);
        exi = apf::boundaryToElementXi(m,down[j],e,xi);
        es->getValues(m,e,exi,values);
        for(int i = 0; i < n; ++i){
          A(row,i) = values[i];
        }
        ++row;
      }
    }
  }
}

} // namespace crv
