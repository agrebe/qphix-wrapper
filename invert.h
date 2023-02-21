#include "global.h"
#include "setup.h"

void invert(typename QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock ** chi_s, // solution
            typename QPhiX::Geometry<OUTER_PREC, OUTER_VECLEN, OUTER_SOALEN, COMPRESS>::FourSpinorBlock ** psi_s, // input
            Params params);
