#pragma once

#include "qdp.h"
#include "global.h"

#include "clover_term_qdp_w.h"
#include "qphix/abs_solver.h"

void setup_QDP(int * argc, char *** argv);

struct Params {
  QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> * geom;
  QDP::multi1d<GAUGE_TYPE> * u;
  //InvRichardsonMultiPrec<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS, INNER_PREC, INNER_VECLEN, INNER_SOALEN, COMPRESS> * solver;
  QPhiX::AbstractSolver<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> * solver;
  QPhiX::QDPCloverTermT<FERM_TYPE, GAUGE_TYPE> * invclov_qdp;
};

Params create_solver(double mass,
                     double clov_coeff,
                     const std::string &filename);
