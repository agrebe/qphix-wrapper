#include "sources.h"
#include "setup.h"
#include "global.h"
#include "invert.h"
#include "params.h"

#include "qphix/clover_dslash_def.h"
#include "qphix/clover_dslash_body.h"

using namespace QDP;


#include "qphix/clover.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#include "qphix/inv_richardson_multiprec.h"

#include <omp.h>

using namespace std;
using namespace QPhiX;



void pion_correlator(double * correlator,
    double ** qphix_ferm,
    const Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> &s)
{
    typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock **psi_s = (Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock**) qphix_ferm;
  // Get the subgrid latt size.
  int Nt = s.Nt();
  int Nz = s.Nz();
  int Ny = s.Ny();
  int Nxh = s.Nxh();
  int nvecs = s.nVecs();
  int Pxy = s.getPxy();
  int Pxyz = s.getPxyz();

  // adapted from qdp_packer_parscalar.h
  for (int cb = 0; cb < 2; cb ++) {
    for (int64_t t = 0; t < Nt; t++) {
      for (int64_t z = 0; z < Nz; z++) {
        for (int64_t y = 0; y < Ny; y++) {
          for (int64_t s = 0; s < nvecs; s++) {
            for (int col = 0; col < 3; col++) {
              for (int spin = 0; spin < 4; spin++) {
                for (int x = 0; x < OUTER_SOALEN; x++) {

                  int ind = t * Pxyz + z * Pxy + y * nvecs + s; //((t*Nz+z)*Ny+y)*nvecs+s;
                  int x_coord = s * OUTER_SOALEN + x;

                  int qdp_ind = ((t * Nz + z) * Ny + y) * Nxh + x_coord;
                  for (int i = 0; i < 2; i ++) { // real versus imag
                    double temp = psi_s[cb][ind][col][spin][i][x];
                    correlator[t] += temp * temp;
                  }
                }
              }
            }
          } // s
        } // y
      } // z
    } // t
  } // cb
}

void run_test(Params * params_pointer, int solve_num) {
  Params params = *params_pointer;
  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock Spinor;
  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::SU3MatrixBlock Gauge;
  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::CloverBlock Clover;
  double start, end;

  const multi1d<int> &lattSize = Layout::subgridLattSize();
  int Nx = lattSize[0];
  int Nt = lattSize[3];

  start = omp_get_wtime();

  QDPIO::cout << "Allocating packed spinor fields" << endl;

  int vol = Nx * Nx * Nx * Nt;
  double ** psi_s = create_ferm(vol);
  double ** chi_s = create_ferm(vol);
  double * correlator = (double *) malloc(sizeof(double) * Nt * 2);

  // Make a random source
  QDPIO::cout << "Initializing QDP++ input spinor" << endl;

  //point_source(psi_s, 0, 0, 0, 0, 0, 0, Nx, Nt);
  int color = solve_num / 4;
  int spin  = solve_num % 4;
  wall_source(psi_s, 0, spin, color, Nx, Nt);
  project_positive(psi_s, Nx, Nt);
  //project_negative(psi_s, Nx, Nt);
  end = omp_get_wtime();
  printf("Pre-inversion setup time: %f sec\n", end - start);

  invert(chi_s, psi_s, (void *) &params);

  // check that we get the same answers after converting to SpinMat format
  double * spin_mat;
  spin_mat = (double *) malloc(sizeof(double) * Nt * Nx * Nx * Nx * 144 * 2);

  start = omp_get_wtime();
  to_spin_mat(spin_mat, chi_s, 0, 0, Nx, Nt);


  // pion correlator
  pion_correlator(correlator, chi_s, *params.geom);
  //pion_correlator_wall_sink(correlator, chi_s, *params.geom);

  for (int t = 0; t < 8; t ++)
    printf("t = %d: %e\n", t, correlator[t]);
  destroy_ferm(psi_s);
  destroy_ferm(chi_s);
  free(spin_mat);
  end = omp_get_wtime();
  printf("Post-inversion time: %f sec\n", end - start);
}


int main(int argc, char **argv) {
  //init_dims(32, 32, 32, 64);
  init_dims(16, 16, 16, 32);
  setup_QDP(&argc, &argv);

  char filename [] = "su3_16_32_b5p87793.lime1000";
  //char filename [] = "su3_32_64_b6p30168.lime1";
  double kappa = 0.109;
  double mass=1.0/(2*kappa) - 4;
  double clov_coeff=1.87567;
  double * u = load_gauge((char*) filename);
  Params * params = (Params*) create_solver(mass, clov_coeff, u);
  for (int j = 0; j < 6; j ++) {
    run_test(params, j);
  }
}
