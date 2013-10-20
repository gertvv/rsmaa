#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// Compare 0-terminated arrays of ints
static int smaa_nAlt;
static int smaa_cmpintarr(void const *p1, void const *p2) {
  return memcmp(*(void const **)p1, *(void const **)p2, sizeof(int) * smaa_nAlt);
}

void smaa_countRankings(
    int const *r,
    int const *nIter, int const *nAlt,
    int *counts) {
  // Sort the rankings
  int const **sorted = malloc(*nIter * sizeof(int *));
  for (int k = 0; k < *nIter; ++k) {
    sorted[k] = r;
    r += *nAlt;
  }
  smaa_nAlt = *nAlt;
  qsort(sorted, *nIter, sizeof(int *), smaa_cmpintarr);

  // Count the unique rankings
  int i = 0;
  ++counts[i];
  for (int k = 1; k < *nIter; ++k) {
    if (memcmp(sorted[k], sorted[k-1], sizeof(int) * *nAlt)) ++i;
    ++counts[i];
  }

  free(sorted);
}
