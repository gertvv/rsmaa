#include <R.h>
#include <Rinternals.h>

/*
 * Calculate piece-wise linear partial value function
 * @param x: vector of measurements (length N)
 * @param y: vector of cut-offs (length n) defining the PVF
 * @param v: vector of values corresponding to the cut-offs (length n)
 * @return vector of values (length N)
 */
SEXP smaa_pvf(SEXP _x, SEXP _y, SEXP _v) {
	int const N = length(_x);
	int const n = length(_y);
	
	_x = PROTECT(coerceVector(_x, REALSXP));
	_y = PROTECT(coerceVector(_y, REALSXP));
	_v = PROTECT(coerceVector(_v, REALSXP));
	double const *x = REAL(_x);
	double const *y = REAL(_y);
	double const *v = REAL(_v);

	SEXP _out = PROTECT(allocVector(REALSXP, N));
	double *out = REAL(_out);

	for (unsigned i = 0; i < N; ++i) {
		unsigned j;
		for (j = 1; j < (n - 1) && y[j] < x[i]; ++j) ;
		out[i] = v[j - 1] + (x[i] - y[j - 1]) * ((v[j] - v[j - 1]) / (y[j] - y[j - 1]));
	}

	UNPROTECT(4);
	return _out;
}
