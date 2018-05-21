#include <R.h>
#include <Rinternals.h>

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

/*
 * Calculate ranks from values
 * @param _v: m * N matrix of values
 * @return m * N array of ranks
 */
SEXP smaa_ranks(SEXP _v);

/*
 * Calculate aggregate values from partial values
 * @param _meas: m * n * N array of partial values
 * @param _pref: n * N matrix of weights, or a single n-vector
 * @param _singleWeight: TRUE iff _pref is a single n-vector
 * @return m * N matrix of values
 */
SEXP smaa_values(SEXP _meas, SEXP _pref, SEXP _singleWeight);

/*
 * Calculate SMAA metrics from partial values
 * @param _meas: m * n * N array of partial values
 * @param _pref: n * N matrix of weights, or a single n-vector
 * @param _singleWeight: TRUE iff _pref is a single n-vector
 * @return A list of (1) m * m matrix of rank acceptabilities; (2) m * n matrix of central weights
 */
SEXP smaa(SEXP _meas, SEXP _pref, SEXP _singleWeight);

/*
 * Calculate piece-wise linear partial value function
 * @param x: vector of measurements (length N)
 * @param y: vector of cut-offs (length n) defining the PVF
 * @param v: vector of values corresponding to the cut-offs (length n)
 * @return vector of values (length N)
 */
SEXP smaa_pvf(SEXP _x, SEXP _y, SEXP _v);

/*
 * Count unique rankings
 * @param m * N matrix of rankings
 * @return N-vector of counts of unique rankings
 */
SEXP smaa_countRankings(SEXP _r);
