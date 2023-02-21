#include <omp.h>

#include "qphix/clover.h"
#include "qphix/qdp_packer.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#include "qphix/inv_richardson_multiprec.h"
#include "dslashm_w.h"

#include "setup.h"
#include "params.h"

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

void setup_QDP(int * argc, char *** argv) {
  CliArgs args;
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
void * create_solver(double mass,
                     double clov_coeff,
                     char * filename_char,
                     int num_solvers)
{
  std::string filename (filename_char);
  CliArgs args;
  Params * params = (Params *) malloc(sizeof(Params));
  params->geom = new Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> 
                  (Layout::subgridLattSize().slice(),
                   args.By,
                   args.Bz,
                   args.NCores,
                   args.Sy,
                   args.Sz,
                   args.PadXY,
                   args.PadXYZ,
                   args.MinCt,
                   true);
  params->u = new multi1d<GAUGE_TYPE>(4);

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

  // What we consider to be small enough...
  QDPIO::cout << "Inititalizing QDP++ gauge field" << endl;

  // Make a random gauge field
  // Start up the gauge field somehow
  // We choose: u = reunit(  1 + factor*gaussian(u) );
  // Adjust factor to get reasonable inversion times in invert test.
  // Bug gauge field is always unitary


  QDPIO::cout << "Reading field from file" << filename << endl;
  XMLReader file_xml, record_xml;
  QDPFileReader from(file_xml, filename, QDPIO_PARALLEL);
  read(from, record_xml, *params->u);
  close(from);

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
    = new Geometry<INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS>(Layout::subgridLattSize().slice(),
                                               args.By,
                                               args.Bz,
                                               args.NCores,
                                               args.Sy,
                                               args.Sz,
                                               args.PadXY,
                                               args.PadXYZ,
                                               args.MinCt,
                                               true);
  QDPIO::cout << "Finished creating geom_inner" << endl;


  QDPIO::cout << "Allocating packged gauge fields" << endl;
  Gauge *packed_gauge_cb0 = (Gauge *)params->geom->allocCBGauge();
  Gauge *packed_gauge_cb1 = (Gauge *)params->geom->allocCBGauge();

  GaugeInner *packed_gauge_cb0_i = (GaugeInner *)geom_inner->allocCBGauge();
  GaugeInner *packed_gauge_cb1_i = (GaugeInner *)geom_inner->allocCBGauge();

  QDPIO::cout << "Allocate Packed Clover Term" << endl;
  Clover *A_cb0 = (Clover *)params->geom->allocCBClov();
  Clover *A_cb1 = (Clover *)params->geom->allocCBClov();
  Clover *A_inv_cb0 = (Clover *)params->geom->allocCBClov();
  Clover *A_inv_cb1 = (Clover *)params->geom->allocCBClov();
  Clover *invclov_packed[2] = {A_inv_cb0, A_inv_cb1};
  Clover *clov_packed[2] = {A_cb0, A_cb1};

  CloverInner *A_cb0_i = (CloverInner *)geom_inner->allocCBClov();
  CloverInner *A_cb1_i = (CloverInner *)geom_inner->allocCBClov();
  CloverInner *A_inv_cb0_i = (CloverInner *)geom_inner->allocCBClov();
  CloverInner *A_inv_cb1_i = (CloverInner *)geom_inner->allocCBClov();
  CloverInner *invclov_packed_i[2] = {A_inv_cb0_i, A_inv_cb1_i};
  CloverInner *clov_packed_i[2] = {A_cb0_i, A_cb1_i};

  QDPIO::cout << "Fields allocated" << endl;

  // Pack the gauge field
  QDPIO::cout << "Packing gauge field...";
  qdp_pack_gauge<>(*params->u, packed_gauge_cb0, packed_gauge_cb1, *params->geom);
  qdp_pack_gauge<>(*params->u, packed_gauge_cb0_i, packed_gauge_cb1_i, *geom_inner);

  Gauge *u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  GaugeInner *u_packed_i[2];
  u_packed_i[0] = packed_gauge_cb0_i;
  u_packed_i[1] = packed_gauge_cb1_i;

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


  params->solver = (AbstractSolver<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> **) 
          malloc(sizeof(AbstractSolver<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> *) * num_solvers);

  for (int solver_num = 0; solver_num < num_solvers; solver_num ++) {
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
    params->solver[solver_num] = new InvRichardsonMultiPrec<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS, INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS>
                        (*M, *solver_inner, 0.1, 5000);
  }
  return (void *) params;
}
