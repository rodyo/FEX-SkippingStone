#include <stdio.h>

#define GRAVITY_EARTH_EGM96
#include "gravity_highprecision.h"

int
main(int argc,
     char** argv)
{
    const double pos_ECEF[3] = {6378136.31, 0.0, 0.0};
    double JD = 2446431.5;
    unsigned int nmax = 12;

    double gradient[3] = 0.0,
           hessian [9] = 0.0;


    gravity_highprecision(pos, JD, nmax,
                          &gradient,
                          &hessian);

    fprintf(stdout,
            "Gradient = {%e, %e, %e}",
            gradient[0], gradient[1], gradient[2]);

}



