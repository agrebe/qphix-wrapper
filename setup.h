#include "global.h"

EXTERN_C void init_dims(int nx, int ny, int nz, int nt);
EXTERN_C void setup_QDP(int * argc, char *** argv);

EXTERN_C void * create_geometry(int inner);

EXTERN_C double * load_gauge(char * filename_char);
EXTERN_C void apply_t_boundary(double * gauge, double t_boundary);
EXTERN_C void pack_gauge(double * gauge,
                         double ** packed_gauge,
                         float ** packed_gauge_inner,
                         void * geometry,
                         void * geometry_inner);

EXTERN_C void * create_solver(double mass,
                              double clov_coeff,
                              void * geom,
                              void * geom_inner,
                              double * u,
                              double ** packed_gauge,
                              float ** packed_gauge_inner);

EXTERN_C void * allocate_spinor(void * params);
EXTERN_C void free_spinor(void * params, void * spinor);
