#include "global.h"

void point_source(QPHIX_FERM(psi_s),
                  int x, int y, int z, int t,
                  int spin, int col,
                  int nx, int nt);

void wall_source(QPHIX_FERM(psi_s),
                 int t,
                 int spin, int col,
                 int nx, int nt);

void project_positive(QPHIX_FERM(psi_s),
                      int nx, int nt);

void project_negative(QPHIX_FERM(psi_s),
                      int nx, int nt);

void to_spin_mat(double * output,
    QPHIX_FERM(psi_s),
    int s1, int c1,
    int nx, int nt);
