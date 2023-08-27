#include "global.h"

EXTERN_C void init_dims(int nx, int ny, int nz, int nt);
EXTERN_C void setup_QDP(int * argc, char *** argv);

EXTERN_C double * read_gauge(char * filename_char);
EXTERN_C void * create_solver(double mass,
                              double clov_coeff,
                              double * u);

EXTERN_C void * allocate_spinor(void * params);
EXTERN_C void free_spinor(void * params, void * spinor);
