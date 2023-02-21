#include "global.h"
#include "qphix/abs_solver.h"

void point_source(typename QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock ** psi_s,
                  int x, int y, int z, int t,
                  int spin, int col,
                  int nx, int nt);

void wall_source(typename QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock ** psi_s,
                 int t,
                 int spin, int col,
                 int nx, int nt);

void project_positive(typename QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock ** psi_s,
                      int nx, int nt);

void project_negative(typename QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock ** psi_s,
                      int nx, int nt);

void to_spin_mat(double * output,
    typename QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock ** phi_s,
    int s1, int c1,
    int nx, int nt);
