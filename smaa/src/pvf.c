#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void smaa_pvf(
    double const *x, int const *N,
    double const *y, double const *v, int const *n,
    double *out) {
  for (unsigned i = 0; i < *N; ++i) {
    unsigned j;
    for (j = 1; j < (*n - 1) && y[j] < x[i]; ++j) ;
    out[i] = v[j - 1] + (x[i] - y[j - 1]) * ((v[j] - v[j - 1]) / (y[j] - y[j - 1]));
  }
}
