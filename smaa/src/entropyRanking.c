#include "smaa.h"

#include <stdlib.h>
#include <string.h>

// Compare 0-terminated arrays of ints
static int smaa_nAlt;
static int smaa_cmpintarr(void const *p1, void const *p2) {
  return memcmp(*(void const **)p1, *(void const **)p2, sizeof(int) * smaa_nAlt);
}

/*
 * Count unique rankings
 * @param m * N matrix of rankings
 * @return N-vector of counts of unique rankings
 */
SEXP smaa_countRankings(SEXP _r) {
  int const nAlt = nrows(_r);
  int const nIter = ncols(_r);

  _r = PROTECT(coerceVector(_r, INTSXP));
  int const *r = INTEGER(_r);

  SEXP _counts = PROTECT(allocVector(INTSXP, nIter));
  int *counts = INTEGER(_counts);
  memset(counts, 0, sizeof(int) * nIter);

  // Sort the rankings
  int const **sorted = malloc(nIter * sizeof(int *));
  for (int k = 0; k < nIter; ++k) {
    sorted[k] = r;
    r += nAlt;
  }
  smaa_nAlt = nAlt;
  qsort(sorted, nIter, sizeof(int *), smaa_cmpintarr);

  // Count the unique rankings
  int i = 0;
  ++counts[i];
  for (int k = 1; k < nIter; ++k) {
    if (memcmp(sorted[k], sorted[k-1], sizeof(int) * nAlt)) ++i;
    ++counts[i];
  }

  free(sorted);

  UNPROTECT(2);
  return _counts;
}
