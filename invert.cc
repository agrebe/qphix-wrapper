#include "invert.h"
#include <omp.h>
#include "qphix/qdp_packer.h"
#include "dslashm_w.h"
#include "setup.h"
#include "params.h"

using namespace QPhiX;

void invert(QPHIX_FERM(chi_s), // solution
            QPHIX_FERM(psi_s), // input
            void * params_pointer,   // inverter, etc.
            int solver_num) {
  Params params = *(Params *) params_pointer;
  // parameters for inversion
  double rsd_target = 1.0e-12;
  int max_iters = 5000;
  int target_cb = 1;
  bool verbose = false;

  // performance metrics
  int niters;
  unsigned long site_flops;
  unsigned long mv_apps;
  double rsd_final;

  FERM_TYPE psi, psi2, clov_psi,
      chi, chi2, clov_chi;

  // unpack psi and precondition
  double time0 = omp_get_wtime();
  qdp_unpack_cb_spinor<>(psi_s[0], psi, *params.geom, 0);
  qdp_unpack_cb_spinor<>(psi_s[1], psi, *params.geom, 1);
  Double betaFactor = Real(0.5); // what is this for?  no idea
  params.invclov_qdp->apply(clov_psi, psi, 1, 0);
  dslash(psi2, *params.u, clov_psi, 1, 1);
  psi[rb[1]] += betaFactor * psi2; // sign fitted to data
  // repackage
  qdp_pack_spinor<>(psi, psi_s[0], psi_s[1], *params.geom);

  qdp_pack_spinor<>(chi, chi_s[0], chi_s[1], *params.geom);
  double time1 = omp_get_wtime();

  // actual inversion
  (*params.solver[solver_num])
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

  // apply D_{ee}^{-1} to even piece
  params.invclov_qdp->apply(chi, psi, 1, 0);
  

  // apply post-conditioner matrix R to result of solve
  //     (  1     - D_{ee}^{-1} D_{eo}  )
  // R = (                              )
  //     (  0                1          )
  // somehow the sign/factor of upper right entry seems wrong?
  // replacing -1 with +0.5 gives result that agrees with other methods
  qdp_unpack_cb_spinor<>(chi_s[1], chi, *params.geom, 1);
  dslash(chi2, *params.u, chi, 1, 0);
  params.invclov_qdp->apply(clov_chi, chi2, 1, 0);
  chi[rb[0]] += betaFactor * clov_chi; // sign fitted to data
  qdp_pack_spinor<>(chi, chi_s[0], chi_s[1], *params.geom);
  double time3 = omp_get_wtime();

  unsigned long num_cb_sites = Layout::vol() / 2;
  unsigned long total_flops =
      (site_flops + (1320 + 504 + 1320 + 504 + 48) * mv_apps) * num_cb_sites;

  printf("RICHARDSON: iters=%d\n", niters);
  printf("RICHARDSON Time for solve=%16.8e sec\n", (time2 - time1));
  printf("RICHARDSON GFLOPS=%e\n",
               1.0e-9 * (double)(total_flops) / (time2 - time1));

  printf("Setup and preconditioner time: %f sec\n", time1 - time0);
  printf("QPhiX inverter time:           %f sec\n", time2 - time1);
  printf("Postconditioner time:          %f sec\n", time3 - time2);
}
