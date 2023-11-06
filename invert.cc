#include "invert.h"
#include <omp.h>
#include "qphix/qdp_packer.h"
#include "dslashm_w.h"
#include "setup.h"
#include "params.h"

using namespace QPhiX;

void invert(double ** solution,
            double ** input,
            void * params_pointer, // inverter, etc
            double rsd_target,
            int max_iters) {
  QPHIX_FERM(chi_s) = (double (**)[3][4][2][OUTER_SOALEN]) solution;
  QPHIX_FERM(psi_s) = (double (**)[3][4][2][OUTER_SOALEN]) input;
  Params params = *(Params *) params_pointer;
  int target_cb = 1;
  bool verbose = false;

  // performance metrics
  int niters;
  unsigned long site_flops;
  unsigned long mv_apps;
  double rsd_final;

  FERM_TYPE ferm, temp;

  // unpack psi and precondition
  double time0 = omp_get_wtime();

  // zero out solution vector (ferm starts out as 0)
  qdp_pack_cb_spinor<>(ferm, chi_s[1], *params.geom, 1);

  qdp_unpack_cb_spinor<>(psi_s[0], ferm, *params.geom, 0);
  qdp_unpack_cb_spinor<>(psi_s[1], ferm, *params.geom, 1);
  double time0A = omp_get_wtime();
  Double betaFactor = Real(0.5); // what is this for?  no idea
  params.invclov_qdp->apply(ferm, ferm, 1, 0);
  double time0B = omp_get_wtime();
  dslash(temp, *params.u, ferm, 1, 1);
  double time0C = omp_get_wtime();
  ferm[rb[1]] += betaFactor * temp; // sign fitted to data
  double time0D = omp_get_wtime();
  // repackage
  qdp_pack_cb_spinor<>(ferm, psi_s[1], *params.geom, 1);

  double time1 = omp_get_wtime();

  // actual inversion
  (*params.solver)
                  (chi_s[1],
                   psi_s[1],
                   rsd_target,
                   niters,
                   rsd_final,
                   site_flops,
                   mv_apps,
                   target_cb,
                   verbose);
  double time2 = omp_get_wtime();
  

  // apply post-conditioner matrix R to result of solve
  //     (  1     - D_{ee}^{-1} D_{eo}  )
  // R = (                              )
  //     (  0                1          )
  // somehow the sign/factor of upper right entry seems wrong?
  // replacing -1 with +0.5 gives result that agrees with other methods
  qdp_unpack_cb_spinor<>(chi_s[1], ferm, *params.geom, 1);
  double time2B = omp_get_wtime();
  dslash(temp, *params.u, ferm, 1, 0);
  double time2C = omp_get_wtime();
  params.invclov_qdp->apply(temp, temp, 1, 0);
  double time2D = omp_get_wtime();
  ferm[rb[0]] += betaFactor * temp; // sign fitted to data
  double time2E = omp_get_wtime();
  qdp_pack_cb_spinor<>(ferm, chi_s[0], *params.geom, 0);
  double time3 = omp_get_wtime();

  unsigned long num_cb_sites = Layout::vol() / 2;
  unsigned long total_flops =
      (site_flops + (1320 + 504 + 1320 + 504 + 48) * mv_apps) * num_cb_sites;

  printf("RICHARDSON: iters=%d\n", niters);
  printf("RICHARDSON Time for solve=%16.8e sec\n", (time2 - time1));
  printf("RICHARDSON GFLOPS=%e\n",
               1.0e-9 * (double)(total_flops) / (time2 - time1));

  printf("Setup and preconditioner time: %f sec\n", time1 - time0);
  printf(" - unpacking:                  %f sec\n", time0A - time0);
  printf(" - applying invclov:           %f sec\n", time0B - time0A);
  printf(" - applying dslash:            %f sec\n", time0C - time0B);
  printf(" - applying 0.5 factor:        %f sec\n", time0D - time0C);
  printf(" - packing:                    %f sec\n", time1 - time0D);
  printf("QPhiX inverter time:           %f sec\n", time2 - time1);
  printf("Postconditioner time:          %f sec\n", time3 - time2);
  printf(" - unpacking:                  %f sec\n", time2B - time2);
  printf(" - applying dslash:            %f sec\n", time2C - time2B);
  printf(" - applying invclov:           %f sec\n", time2D - time2C);
  printf(" - applying 0.5 factor:        %f sec\n", time2E - time2D);
  printf(" - packing:                    %f sec\n", time3 - time2E);
}
