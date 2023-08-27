#include "global.h"

EXTERN_C double ** create_ferm(int vol);
EXTERN_C void destroy_ferm(double** ferm);

EXTERN_C void point_source(double ** qphix_ferm,
                           int x, int y, int z, int t,
                           int spin, int col,
                           int nx, int nt);
EXTERN_C void wall_source(double ** qphix_ferm,
                          int t,
                          int spin, int col,
                          int nx, int nt);

EXTERN_C void project_positive(double ** qphix_ferm,
                               int nx, int nt);
EXTERN_C void project_negative(double ** qphix_ferm,
                               int nx, int nt);

EXTERN_C void to_spin_mat(double * output,
             double ** qphix_ferm,
             int s1, int c1,
             int nx, int nt);

EXTERN_C void to_qc_single(float * output,
    double ** qphix_ferm,
    int nx, int nt);
EXTERN_C void from_qc_single(float * input,
    double ** qphix_ferm,
    int nx, int nt);
EXTERN_C void to_qc_double(double * output,
    double ** qphix_ferm,
    int nx, int nt);
EXTERN_C void from_qc_double(double * input,
    double ** qphix_ferm,
    int nx, int nt);
