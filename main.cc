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
    typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock **psi_s,
    const Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> &s)
{
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

// FIXME: I think this method is wrong but don't know why
void pion_correlator_wall_sink(double * correlator,
    typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock **psi_s,
    const Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> &s)
{
  // Get the subgrid latt size.
  int Nt = s.Nt();
  int Nz = s.Nz();
  int Ny = s.Ny();
  int Nxh = s.Nxh();
  int nvecs = s.nVecs();
  int Pxy = s.getPxy();
  int Pxyz = s.getPxyz();
  printf("Pxy = %d\n", Pxy);
  printf("Pxyz = %d\n", Pxyz);
  int nt = Nt;
  int nx = Nz;
  for (int t = 0; t < nt; t ++) {
    double wall_sink_re, wall_sink_im;
    for (int z = 0; z < nx; z ++) {
      for (int y = 0; y < nx; y ++) {
        for (int x_full = 0; x_full < nx; x_full ++) {
          int x = x_full;
          // find index of position
          // first find checkerboard
          int cb = (x + y + z + t) % 2;
          x /= 2;
          // then find vector index
          int v = x % OUTER_SOALEN;
          x /= OUTER_SOALEN;
          // finally get index in lattice
          // x runs fastest, t is slowest
          int index = ((t * nx + z) * nx + y) * nx / (OUTER_SOALEN * 2) + x;
          for (int col = 0; col < 3; col++) {
            for (int spin = 0; spin < 4; spin++) {
              wall_sink_re += psi_s[cb][index][col][spin][0][v];
              wall_sink_im += psi_s[cb][index][col][spin][1][v];
            }
          }
        }
      }
    }
    correlator[t] += wall_sink_re * wall_sink_re;
                   + wall_sink_im * wall_sink_im;
  }
}

void run_test(Params * params_pointer, int solve_num) {
  Params params = *params_pointer;
  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock Spinor;
  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::SU3MatrixBlock Gauge;
  typedef typename Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::CloverBlock Clover;

  const multi1d<int> &lattSize = Layout::subgridLattSize();
  int Nx = lattSize[0];
  int Nt = lattSize[3];

  QDPIO::cout << "Allocating packed spinor fields" << endl;
  Spinor * psi_even, * psi_odd, * chi_even, * chi_odd;
  #pragma omp critical
  {
    psi_even = (Spinor *)params.geom->allocCBFourSpinor();
    psi_odd = (Spinor *)params.geom->allocCBFourSpinor();
    chi_even = (Spinor *)params.geom->allocCBFourSpinor();
    chi_odd = (Spinor *)params.geom->allocCBFourSpinor();
  }
  Spinor *psi_s[2] = {psi_even, psi_odd};
  Spinor *chi_s[2] = {chi_even, chi_odd};

  // Make a random source
  QDPIO::cout << "Initializing QDP++ input spinor" << endl;

  // pion correlator
  double * correlator = (double*) calloc(Nt, sizeof(double));

  //point_source(psi_s, 0, 0, 0, 0, 0, 0, Nx, Nt);
  int color = solve_num / 4;
  int spin  = solve_num % 4;
  wall_source(psi_s, 0, spin, color, Nx, Nt);
  project_positive(psi_s, Nx, Nt);
  //project_negative(psi_s, Nx, Nt);

  invert(chi_s, psi_s, (void *) &params, solve_num);

  printf("psi_s[0][0][0][0][0][0] = %f\n", psi_s[0][0][0][0][0][0]);
  printf("psi_s[1][0][0][0][0][0] = %f\n", psi_s[1][0][0][0][0][0]);
  // real and imaginary parts at 0
  printf("Re[chi_s(0,0,0,0)] = %f\n", chi_s[0][0][0][0][0][0]);
  printf("Im[chi_s(0,0,0,0)] = %f\n", chi_s[0][0][0][0][1][0]);
  // real part at (1, 0, 0, 0)
  printf("Re[chi_s(1,0,0,0)] = %f\n", chi_s[1][0][0][0][0][0]);
  // real part at (2, 0, 0, 0)
  printf("Re[chi_s(2,0,0,0)] = %f\n", chi_s[0][0][0][0][0][1]);
  printf("Re[chi_s(8,0,0,0)] = %f\n", chi_s[0][0][0][0][0][4]);
  printf("Re[chi_s(16,0,0,0)] = %f\n", chi_s[0][1][0][0][0][0]);
  printf("Re[chi_s(0,1,0,0)] = %f\n", chi_s[1][2][0][0][0][0]);

  // check that we get the same answers after converting to SpinMat format
  double * spin_mat;
  #pragma omp critical
  {
    spin_mat = (double *) malloc(sizeof(double) * Nt * Nx * Nx * Nx * 144 * 2);
  }

  double start = omp_get_wtime();
  to_spin_mat(spin_mat, chi_s, 0, 0, Nx, Nt);
  double end = omp_get_wtime();

  printf("Conversion to SpinMat format: %f sec\n", end - start);

  // real and imaginary parts at 0
  printf("Re[chi_s(0,0,0,0)] = %f\n", spin_mat[0]);
  printf("Im[chi_s(0,0,0,0)] = %f\n", spin_mat[1]);
  // real part at (1, 0, 0, 0)
  printf("Re[chi_s(1,0,0,0)] = %f\n", spin_mat[288]);
  // real part at (2, 0, 0, 0)
  printf("Re[chi_s(2,0,0,0)] = %f\n", spin_mat[288*2]);
  printf("Re[chi_s(8,0,0,0)] = %f\n", spin_mat[288*8]);
  printf("Re[chi_s(16,0,0,0)] = %f\n", spin_mat[288*16]);
  printf("Re[chi_s(0,1,0,0)] = %f\n", spin_mat[288*32]);

  // pion correlator
  pion_correlator(correlator, chi_s, *params.geom);
  //pion_correlator_wall_sink(correlator, chi_s, *params.geom);

  for (int t = 0; t < 8; t ++)
    printf("t = %d: %e\n", t, correlator[t]);
  #pragma omp critical
  {
    params.geom->free(psi_even);
    params.geom->free(psi_odd);
    params.geom->free(chi_even);
    params.geom->free(chi_odd);
    free(spin_mat);
  }
}


int main(int argc, char **argv) {
  setup_QDP(&argc, &argv);

  char filename [] = "cl3_32_48_b6p1_m0p2450-sgf.lime";
  double mass=-0.245;
  double clov_coeff=1.24930970916466;
  int num_solvers = 6;
  Params * params = (Params*) create_solver(mass, clov_coeff, (char*) filename, num_solvers);
  omp_set_nested(1);
  for (int j = 0; j < 1; j ++) {
    #pragma omp parallel for
    for (int i = 0; i < num_solvers; i ++) {
      run_test(params, i);
    }
  }
}
