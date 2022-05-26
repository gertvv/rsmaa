#include "smaa.h"
#include <R_ext/BLAS.h>
#ifndef FCONE
#define FCONE
#endif

#include <math.h>
#include <stdlib.h>

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

SEXP smaa_ranks(SEXP _v) {
  int const nAlt = nrows(_v);
  int const nIter = ncols(_v);

  _v = PROTECT(coerceVector(_v, REALSXP));
  double const *v = REAL(_v);

  SEXP _r = PROTECT(allocMatrix(INTSXP, nAlt, nIter));
  int *r = INTEGER(_r);

  for (int k = 0; k < nIter; ++k) {
    smaa_rank(v, r, nAlt);

    v += nAlt;
    r += nAlt;
  }

  UNPROTECT(2);
  return _r;
}

SEXP smaa_values(SEXP _meas, SEXP _pref, SEXP _singleWeight) {
  SEXP dim = getAttrib(_meas, R_DimSymbol);
  int const nAlt  = INTEGER(dim)[0];
  int const nCrit = INTEGER(dim)[1];
  int const nIter = INTEGER(dim)[2];
  int const singleWeight = asLogical(_singleWeight);
  const int inc1 = 1;
  const double one = 1.0, zero = 0.0; // for BLAS
  const char trans = 'N';

  _meas = PROTECT(coerceVector(_meas, REALSXP));
  _pref = PROTECT(coerceVector(_pref, REALSXP));
  double const *meas = REAL(_meas);
  double const *pref = REAL(_pref);

  SEXP _v = PROTECT(allocMatrix(REALSXP, nAlt, nIter));
  double *v = REAL(_v);

  for (int k = 0; k < nIter; ++k) {
    // calculate value of each alternative
    F77_CALL(dgemv)(&trans, &nAlt, &nCrit,
      &one, meas, &nAlt, pref, &inc1,
      &zero, v, &inc1 FCONE); // t := 1Aw + 0t

    // advance measurement and weight pointers
    meas += nAlt * nCrit;
    if (!singleWeight) pref += nCrit;
    // advance alternative value pointer
    v += nAlt;
  }

  UNPROTECT(3);
  return _v;
}

SEXP smaa(SEXP _meas, SEXP _pref, SEXP _singleWeight) {
  SEXP dim = getAttrib(_meas, R_DimSymbol);
  int const nAlt  = INTEGER(dim)[0];
  int const nCrit = INTEGER(dim)[1];
  int const nIter = INTEGER(dim)[2];
  int const singleWeight = asLogical(_singleWeight);
  const int inc1 = 1;
  const double one = 1.0, zero = 0.0; // for BLAS
  const char trans = 'N';

  _meas = PROTECT(coerceVector(_meas, REALSXP));
  _pref = PROTECT(coerceVector(_pref, REALSXP));
  double const *meas = REAL(_meas);
  double const *pref = REAL(_pref);

  SEXP _h = PROTECT(allocMatrix(REALSXP, nAlt, nAlt));
  SEXP _wc = PROTECT(allocMatrix(REALSXP, nAlt, nCrit));
  Matrix h = { REAL(_h), nAlt, nAlt };
  memset(h.data, 0, nAlt * nAlt * sizeof(double));
  Matrix wc = { REAL(_wc), nAlt, nCrit };
  memset(wc.data, 0, nAlt * nCrit * sizeof(double));

  double t[nAlt];
  int r[nAlt]; // alternative ranks
  for (int k = 0; k < nIter; ++k) {
    // calculate value of each alternative
    F77_CALL(dgemv)(&trans, &nAlt, &nCrit,
      &one, meas, &nAlt, pref, &inc1,
      &zero, t, &inc1 FCONE); // t := 1Aw + 0t

    // rank the alternatives
    smaa_rank(t, r, nAlt);
    for (int i = 0; i < nAlt; ++i) {
      *get(&h, i, r[i]) = *get(&h, i, r[i]) + 1; // update rank counts
      if (!singleWeight && r[i] == 0) { // update central weights
        for (int j = 0; j < nCrit; ++j) {
          *get(&wc, i, j) = *get(&wc, i, j) + pref[j];
        }
      }
    }

    // advance measurement and weight pointers
    meas += nAlt * nCrit;
    if (!singleWeight) pref += nCrit;
  }

  // normalize central weights
  if (!singleWeight) {
    for (int i = 0; i < nAlt; ++i) {
      double const r1 = *get(&h, i, 0);
      if (r1 > 0.0) {
        for (int j = 0; j < nCrit; ++j) {
          *get(&wc, i, j) = *get(&wc, i, j) / r1;
        }
      }
    }
  }

  SEXP ans = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ans, 0, _h);
  SET_VECTOR_ELT(ans, 1, _wc);

  UNPROTECT(5);
  return ans;
}
