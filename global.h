#pragma once

#define OUTER_PREC   double
#define OUTER_VECLEN 4
#define OUTER_SOALEN 4
#define COMPRESS     true
#define INNER_PREC   float
#define INNER_VECLEN 8
#define INNER_SOALEN 4
#define GAUGE_TYPE   QDP::LatticeColorMatrixD
#define FERM_TYPE    QDP::LatticeDiracFermionD

#define QPHIX_FERM(x) double (** x) [3][4][2][OUTER_SOALEN]

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif
