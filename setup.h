#include "global.h"

EXTERN_C void setup_QDP(int * argc, char *** argv);

EXTERN_C void * create_solver(double mass,
                              double clov_coeff,
                              char * filename_char);

EXTERN_C void * allocate_spinor(void * params);
EXTERN_C void free_spinor(void * params, void * spinor);
