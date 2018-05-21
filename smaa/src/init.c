#include "smaa.h"
#include <R_ext/Visibility.h>

static const R_CallMethodDef callMethods[] = {
  { "smaa_countRankings", (DL_FUNC) &smaa_countRankings, 1 },
  { "smaa_pvf", (DL_FUNC) &smaa_pvf, 3 },
  { "smaa_ranks", (DL_FUNC) &smaa_ranks, 1 },
  { "smaa_values", (DL_FUNC) &smaa_values, 3 },
  { "smaa_smaa", (DL_FUNC) &smaa, 3 },
  { NULL, NULL, 0 }
};

void attribute_visible R_init_smaa(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
