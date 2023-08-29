#include <omp.h>

#include "qphix/clover.h"
#include "qphix/qdp_packer.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#include "qphix/inv_richardson_multiprec.h"
#include "dslashm_w.h"

#include "setup.h"
#include "params.h"

struct CliArgs {
  int nrow_in[4];
  int By;
  int Bz;
  int PadXY;
  int PadXYZ;
  int NCores;
  int Sy;
  int Sz;
  int MinCt;
};

CliArgs args;

// WARNING: This is not thread safe
// Wrap calls to this in #pragma omp critical if needed
void * allocate_spinor(void * params) {
  void * spinor;
  spinor = ((Params *) params)->geom->allocCBFourSpinor();
  return spinor;
}

// also not thread safe
void free_spinor(void * params, void * spinor) {
  ((Params *) params)->geom->free(
    (typename QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock *) spinor);
}

void init_dims(int nx, int ny, int nz, int nt) {
  args.nrow_in[0] = nx;
  args.nrow_in[1] = ny;
  args.nrow_in[2] = nz;
  args.nrow_in[3] = nt;
}

void setup_QDP(int * argc, char *** argv) {
  args.By = 4;
  args.Bz = 4;
  args.PadXY= 0;
  args.PadXYZ = 0;
  args.NCores = 6;
  args.Sy = 1;
  args.Sz = 1;
  args.MinCt = 1;

  printf("Initializing QDP\n");
  QDP_initialize(argc, argv);
  printf("QDP initialized\n");

  multi1d<int> nrow(Nd);
  for (int i = 0; i < Nd; ++i) {
    nrow[i] = args.nrow_in[i];
  }
  Layout::setLattSize(nrow);
  printf("Creating layout\n");
  Layout::create();
  printf("Finished creating layout\n");

  omp_set_num_threads(args.NCores * args.Sy * args.Sz);
  printf("Test setup done\n");
}

using namespace std;
using namespace QPhiX;

void * create_geometry(int inner) {
  if (inner)
    return (void*) (new Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS> 
                  (Layout::subgridLattSize().slice(),
                   args.By,
                   args.Bz,
                   args.NCores,
                   args.Sy,
                   args.Sz,
                   args.PadXY,
                   args.PadXYZ,
                   args.MinCt,
                   true));
  else
    return (void*) (new Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> 
                  (Layout::subgridLattSize().slice(),
                   args.By,
                   args.Bz,
                   args.NCores,
                   args.Sy,
                   args.Sz,
                   args.PadXY,
                   args.PadXYZ,
                   args.MinCt,
                   true));
}

double * load_gauge(char * filename_char) {
  QDP::multi1d<GAUGE_TYPE> * u = new QDP::multi1d<GAUGE_TYPE>(4);
  std::string filename (filename_char);
  QDPIO::cout << "Inititalizing QDP++ gauge field" << endl;

  // Make a random gauge field
  // Start up the gauge field somehow
  // We choose: u = reunit(  1 + factor*gaussian(u) );
  // Adjust factor to get reasonable inversion times in invert test.
  // Bug gauge field is always unitary


  QDPIO::cout << "Reading field from file" << filename << endl;
  XMLReader file_xml, record_xml;
  QDPFileReader from(file_xml, filename, QDPIO_PARALLEL);
  read(from, record_xml, *u);
  close(from);
  return (double*) u;
}

void pack_gauge(double * gauge,
                double ** packed_gauge,
                float ** packed_gauge_inner,
                void * geometry,
                void * geometry_inner) {
  Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> * geom 
    = (Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> *)
                    geometry;
  Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS> * geom_inner
    = (Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS> *)
                    geometry_inner;
  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::SU3MatrixBlock Gauge;
  typedef typename Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS>::SU3MatrixBlock GaugeInner;
  QDPIO::cout << "Allocating packed gauge fields" << endl;
  packed_gauge[0] = (double*) geom->allocCBGauge();
  packed_gauge[1] = (double*) geom->allocCBGauge();

  packed_gauge_inner[0] = (float*) geom_inner->allocCBGauge();
  packed_gauge_inner[1] = (float*) geom_inner->allocCBGauge();

  // Pack the gauge field
  QDPIO::cout << "Packing gauge field...";
  qdp_pack_gauge<>(*(QDP::multi1d<GAUGE_TYPE>*) gauge, 
                   (Gauge*) packed_gauge[0], 
                   (Gauge*) packed_gauge[1], *geom);
  qdp_pack_gauge<>(*(QDP::multi1d<GAUGE_TYPE>*) gauge, 
                   (GaugeInner*) packed_gauge_inner[0], 
                   (GaugeInner*) packed_gauge_inner[1], *geom_inner);

  QDPIO::cout << "done" << endl;
}

void * create_solver(double mass,
                     double clov_coeff,
                     void * geometry,
                     void * geometry_inner,
                     double * u,
                     double ** packed_gauge,
                     float ** packed_gauge_inner) {
  Params * params = (Params *) malloc(sizeof(Params));
  params->geom = (Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> *)
                    geometry;
  params->u = (QDP::multi1d<GAUGE_TYPE> *) u;

  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock Spinor;
  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::SU3MatrixBlock Gauge;
  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::CloverBlock Clover;

  typedef typename Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS>::FourSpinorBlock SpinorInner;
  typedef typename Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS>::SU3MatrixBlock GaugeInner;
  typedef typename Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS>::CloverBlock CloverInner;
  QDPIO::cout << endl << "ENTERING CLOVER DSLASH TEST" << endl;

  // Diagnostic information:
  const multi1d<int> &lattSize = Layout::subgridLattSize();
  int Nx = lattSize[0];
  int Ny = lattSize[1];
  int Nz = lattSize[2];
  int Nt = lattSize[3];

  // Set anisotropy parameters -- pick some random numbers
  double xi_0_f = 1.0;
  double nu_f = 1.0;

  // This is what makeFermCoeffs does under the hood.
  // Get the direct spatial and temporal anisotropy factors
  double aniso_fac_s = (double)(1);
  double aniso_fac_t = (double)(1);

  QDPIO::cout << "Setting Clover Term Parameters" << endl;

  CloverFermActParams clparam;
  AnisoParam_t aniso;

  // Aniso prarams
  aniso.anisoP = false;
  aniso.xi_0 = 1;
  aniso.nu = 1;
  aniso.t_dir = 3;

  // Set up the Clover params
  clparam.anisoParam = aniso;

  // Some mass
  // Now use the real mass

  clparam.Mass = Real(mass);

  // Some random clover coeffs
  clparam.clovCoeffR = Real(clov_coeff);
  clparam.clovCoeffT = Real(clov_coeff);

  double t_boundary = double(-1);
  // Create Dslash
  Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS> * geom_inner
    = (Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS> *) geometry_inner;

  QDPIO::cout << "Allocate Packed Clover Term" << endl;
  Clover *A_cb0 = params->geom->allocCBClov();
  Clover *A_cb1 = params->geom->allocCBClov();
  Clover *A_inv_cb0 = params->geom->allocCBClov();
  Clover *A_inv_cb1 = params->geom->allocCBClov();
  Clover *invclov_packed[2] = {A_inv_cb0, A_inv_cb1};
  Clover *clov_packed[2] = {A_cb0, A_cb1};

  CloverInner *A_cb0_i = geom_inner->allocCBClov();
  CloverInner *A_cb1_i = geom_inner->allocCBClov();
  CloverInner *A_inv_cb0_i = geom_inner->allocCBClov();
  CloverInner *A_inv_cb1_i = geom_inner->allocCBClov();
  CloverInner *invclov_packed_i[2] = {A_inv_cb0_i, A_inv_cb1_i};
  CloverInner *clov_packed_i[2] = {A_cb0_i, A_cb1_i};

  QDPIO::cout << "Fields allocated" << endl;

  Gauge **u_packed = (Gauge**) packed_gauge;
  GaugeInner **u_packed_i = (GaugeInner**) packed_gauge_inner;

  QDPIO::cout << "done" << endl;

  QDPIO::cout << " Packing fermions...";


  QDPIO::cout << "done" << endl;

  QDPIO::cout << "Creating the Clover Term " << endl;

  // Clover term deals with anisotropy internally -- so use original u field.
  QDPCloverTermT<FERM_TYPE, GAUGE_TYPE> clov_qdp;

  QDPIO::cout << "Adding on boundary field" << endl;
  // Modify u (antiperiodic BC's)
  (*params->u)[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3] - 1),
                Real(t_boundary),
                Real(1));

  clov_qdp.create(*params->u, clparam);
  QDPIO::cout << "Inverting Clover Term" << endl;
  params->invclov_qdp = new QDPCloverTermT<FERM_TYPE, GAUGE_TYPE>(clov_qdp);
  for (int cb = 0; cb < 2; cb++) {
    params->invclov_qdp->choles(cb);
  }
  QDPIO::cout << "Done" << endl;

  QDPIO::cout << "Packing Clover term..." << endl;
  for (int cb = 0; cb < 2; cb++) {
    qdp_pack_clover<>(*params->invclov_qdp, invclov_packed[cb], *params->geom, cb);
    qdp_pack_clover<>(*params->invclov_qdp, invclov_packed_i[cb], *geom_inner, cb);
    qdp_pack_clover<>(clov_qdp, clov_packed[cb], *params->geom, cb);
    qdp_pack_clover<>(clov_qdp, clov_packed_i[cb], *geom_inner, cb);
  }
  QDPIO::cout << "Done" << endl;

  int max_iters = 5000;

  EvenOddCloverOperator<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> * M
    = new EvenOddCloverOperator<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> (u_packed,
                                              clov_packed[1],
                                              invclov_packed[0],
                                              params->geom,
                                              t_boundary,
                                              aniso_fac_s,
                                              aniso_fac_t);

  EvenOddCloverOperator<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS> * M_inner
    = new EvenOddCloverOperator<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS>(u_packed_i,
                                                         clov_packed_i[1],
                                                         invclov_packed_i[0],
                                                         geom_inner,
                                                         t_boundary,
                                                         aniso_fac_s,
                                                         aniso_fac_t);
  InvBiCGStab<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS> *solver_inner
    = new InvBiCGStab<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS>(*M_inner, max_iters);
  params->solver = new InvRichardsonMultiPrec<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS, INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS>
                      (*M, *solver_inner, 0.1, 5000);

  // free the clover pieces we no longer need
  free(clov_packed[0]);
  free(clov_packed_i[0]);
  // NOTE: This leads to problems if we try to free invclov_i.  Why?
  free(invclov_packed[1]);

  return (void *) params;
}
