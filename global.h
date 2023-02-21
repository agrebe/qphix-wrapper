#pragma once

#define OUTER_PREC   double
#define OUTER_VECLEN 8
#define OUTER_SOALEN 8
#define COMPRESS     true
#define INNER_PREC   float
#define INNER_VECLEN 16
#define INNER_SOALEN 8
#define GAUGE_TYPE   QDP::LatticeColorMatrixD
#define FERM_TYPE    QDP::LatticeDiracFermionD

#define QPHIX_FERM(x) double (** x) [3][4][2][OUTER_SOALEN]

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

struct CliArgs {
  int nrow_in[4] = {32, 32, 32, 48};

  int By = 4;
  int Bz = 4;
  int PadXY = 0;
  int PadXYZ = 0;
  int NCores = 16;
  int Sy = 1;
  int Sz = 1;
  int MinCt = 1;
};
