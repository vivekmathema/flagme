
#include "R_ext/Rdynload.h"
#include "flagme.h"

static const R_CMethodDef cMethods[]  = {
  {"dp",(DL_FUNC)&dp,8},
  {"cos_ndp_himem", (DL_FUNC)&cos_ndp_himem,9},
  {"cos_ndp_lowmem", (DL_FUNC)&cos_ndp_lowmem,10},
  {NULL, NULL, 0}
};

void R_init_flagme(DllInfo *info){
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
