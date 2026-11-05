/* conf_c.c: C side of conf.f90. Complex values cross the language boundary by value both
 * ways, and an array arrives by reference: the calling conventions sander's C helpers use. */
#include <complex.h>

double complex f_cadd(double complex a, double complex b);

double complex c_cmul(double complex a, double complex b) { return a * b; }

double c_dot3(const double *p, const double *q) { return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]; }

void c_calls_fortran(double complex *r) { *r = f_cadd(1.5 + 2.0 * I, 2.5 - 3.0 * I); }
