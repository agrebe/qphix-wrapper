#pragma once

#include "qdp.h"
#include "global.h"

#include "clover_term_qdp_w.h"
#include "qphix/abs_solver.h"

struct Params {
  QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> * geom;
  QDP::multi1d<GAUGE_TYPE> * u;
  QPhiX::AbstractSolver<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS> ** solver;
  QPhiX::QDPCloverTermT<FERM_TYPE, GAUGE_TYPE> * invclov_qdp;
};
