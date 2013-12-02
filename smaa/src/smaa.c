#include <R.h>
#include <R_ext/BLAS.h>

#include <math.h>
#include <stdlib.h>

#include <stdio.h>

typedef struct Matrix {
  double * const data;
  int const nRow;
  int const nCol;
} Matrix;

/**
 * @param i Row index.
 * @param j Column index.
 */
static inline double *get(Matrix *m, int i, int j) {
  return m->data + j * (m->nRow) + i;
}

/**
 * Rank the n-array of doubles t, writing results to r.
 * The rank of t[i] is the number of elements greater than t[i].
 * Assumes that n is small (complexity O(n^2)).
 */
static inline void smaa_rank(double const *t, int *r, int n) {
  for (int i = 0; i < n; ++i) {
    r[i] = 0;
    for (int j = 0; j < n; ++j) {
      if (t[j] > t[i]) {
        ++r[i];
      }
    }
  }
}

void smaa_ranks(
    double const *v,
    int const *nIter, int const *nAlt,
    int *r) {
  for (int k = 0; k < *nIter; ++k) {
    smaa_rank(v, r, *nAlt);

    v += *nAlt;
    r += *nAlt;
  }
}

void smaa_values(
    double const *meas, double const *pref,
    int const *nIter, int const *nAlt, int const *nCrit,
    int const *singleWeight,
    double *v) {
  const int inc1 = 1;
  const double one = 1.0, zero = 0.0; // for BLAS
  const char trans = 'N';

  for (int k = 0; k < *nIter; ++k) {
    // calculate value of each alternative
    F77_CALL(dgemv)(&trans, nAlt, nCrit,
      &one, meas, nAlt, pref, &inc1,
      &zero, v, &inc1); // t := 1Aw + 0t

    // advance measurement and weight pointers
    meas += *nAlt * *nCrit;
    if (!*singleWeight) pref += *nCrit;
    // advance alternative value pointer
    v += *nAlt;
  }
}

#include <stdio.h>

void smaa(
    double const *meas, double const *pref,
    int const *nIter, int const *nAlt, int const *nCrit,
    int const *singleWeight,
    double *hData, double *wcData) {
  const int inc1 = 1;
  const double one = 1.0, zero = 0.0; // for BLAS
  const char trans = 'N';

  Matrix h = { hData, *nAlt, *nAlt };
  Matrix wc = { wcData, *nAlt, *nCrit };

  double t[*nAlt];
  int r[*nAlt]; // alternative ranks
  for (int k = 0; k < *nIter; ++k) {
    // calculate value of each alternative
    F77_CALL(dgemv)(&trans, nAlt, nCrit,
      &one, meas, nAlt, pref, &inc1,
      &zero, t, &inc1); // t := 1Aw + 0t

    // rank the alternatives
    smaa_rank(t, r, *nAlt);
    for (int i = 0; i < *nAlt; ++i) {
      *get(&h, i, r[i]) = *get(&h, i, r[i]) + 1; // update rank counts
      if (!*singleWeight && r[i] == 0) { // update central weights
        for (int j = 0; j < *nCrit; ++j) {
          *get(&wc, i, j) = *get(&wc, i, j) + pref[j];
        }
      }
    }

    // advance measurement and weight pointers
    meas += *nAlt * *nCrit;
    if (!*singleWeight) pref += *nCrit;
  }

  // normalize central weights
  if (!*singleWeight) {
    for (int i = 0; i < *nAlt; ++i) {
      double const r1 = *get(&h, i, 0);
      if (r1 > 0.0) {
        for (int j = 0; j < *nCrit; ++j) {
          *get(&wc, i, j) = *get(&wc, i, j) / r1;
        }
      }
    }
  }
}
