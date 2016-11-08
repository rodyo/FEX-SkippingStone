#ifndef GRAVITY_HIGHPRECISION_H
#define GRAVITY_HIGHPRECISION_H

void
gravity_highprecision(const double pos_ECEF[3],
                      double JD,
                      unsigned int nmax,
                      int compute_gradient,   double (*gradient)[3],
                      int compute_hessian,    double (*hessian)[3][3]);

#endif