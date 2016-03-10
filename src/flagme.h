#include <R.h>
#include <Rinternals.h>

void cos_ndp_lowmem(double *score, long *nr_, long *nc1_, long *nc2_, double *x1, double *x2, double *t1, double *t2, double *D_, long *df_);
void cos_ndp_himem(double *score, long *nr_, long *nc1_, long *nc2_, double *x1, double *x2, double *D, long *df_, double *timedf);
void dp( double *D, double *M, double *gap_, int *phi, long *nr_, long *nc_, int *match, int *nmatch);
void pearson(int* size, double* x, double* y, double* result);
