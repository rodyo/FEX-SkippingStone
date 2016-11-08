#include <math.h>
#include "gravity_highprecision.h"
#include "constants.h"


#if   defined(GRAVITY_EARTH_EGM96)
    #include "EGM96_parameters.h"

#elif defined(GRAVITY_EARTH_EGM2008)
    #include "EGM2008_parameters.h"

#elif defined(GRAVITY_MOON_LP100K)
    #include "LP100K_parameters.h"

#elif defined(GRAVITY_MOON_LP165P)
    #include "LP165P_parameters.h"

#elif defined(GRAVITY_MARS_GMM2B)
    #include "GMM2B_parameters.h"

#else
    #error Missing gravity model definition.
#endif


/*
 * References:
 * [1] "Fast Gravity, Gravity Partials, Normalized Gravity, Gravity Gradient
 *      Torque and Magnetic Field: Derivation, Code and Data". 1993.
 *      Robert G. Gottlieb,  McDonnell Douglas Space Systems, Houston
 *      Devision, Houston Texas. Prepared for Lyndon B. Johnson Space Center,
 *      contract NAS9-17885.
*/


static double nn[GRAVITY_COEFF_MAXDEGREE] = 0.0,
              mn[GRAVITY_COEFF_MAXDEGREE] = 0.0;


/* Fully normalized associated Legendre polynomials */
static
void
ComputeLegendrePolys(double phi,
                     unsigned int nmax,
                     double (*Pnm[GRAVITY_COEFF_MAXDEGREE+3u][GRAVITY_COEFF_MAXDEGREE+3u]))
{
    /* Computation based on equations 7-9 and 7-10 from [1] */

    static unsigned int Pnm_initialized = 0u;

    static double  xi[GRAVITY_COEFF_MAXDEGREE][GRAVITY_COEFF_MAXDEGREE] = 0.0,
                  eta[GRAVITY_COEFF_MAXDEGREE][GRAVITY_COEFF_MAXDEGREE] = 0.0;

    unsigned int n,
                 m;

    double cphi = cos(PI_OVER_2 - phi),
           sphi = sin(PI_OVER_2 - phi);

    if (fabs(cphi) <= EPS) cphi = 0.0;
    if (fabs(sphi) <= EPS) sphi = 0.0;

    if (!Pnm_initialized)
    {
        /* TODO: implement eqs. 7-10 and 7-12 */
    }



    (*Pnm)[0][0] = 1.0;
    (*Pnm)[1][0] = SQUAREROOT3 * cphi;
    (*Pnm)[1][1] = SQUAREROOT3 * sphi;

    for (n=2u; n<nmax+3u; ++n)
    {
        unsigned int nm1 = n - 1u;

        double nm1d  = nn[n] - 1.0,
               n2    = nn[n] * 2.0,
               sn2   = sqrt(n2+0.0),  sn2p1 = sqrt(n2+1.0),
               sn2m3 = sqrt(n2-3.0),  sn2m1 = sqrt(n2-1.0);

        (*Pnm)[n][0] = sn2p1/nn[n] *
                         (sn2m1*cphi*(*Pnm)[nm1][0u] - nm1d/sn2m3 * (*Pnm)[n-2u][0u]);

        for (m=1u; m<n+1u; ++m) {
            (*Pnm)[n][m] = sn2p1 / ( sqrt(nn[n]+mm[m]) * sqrt(nn[n]-mm[m]) ) *
                             (sn2m1*cphi*(*Pnm)[nm1][m] - sqrt(nm1d+mm[m]) * sqrt(nm1d-mm[m]) * (*Pnm)[n-2u][m]/sn2m3);
        }

        (*Pnm)[n][n] = sn2p1/sn2 * sphi*(*Pnm)[nm1][nm1];

    }

}


void
gravity_highprecision(const double pos_ECEF[3],
                      double JD,
                      unsigned int nmax,
                      int compute_gradient,   double (*gradient)[3],
                      int compute_hessian,    double (*hessian)[3][3])
{
    static unsigned int grav_initialized = 0u;

    double Pnm[GRAVITY_COEFF_MAXDEGREE+3u][GRAVITY_COEFF_MAXDEGREE+3u] = 0.0,
           cmlambda[GRAVITY_COEFF_MAXDEGREE+1u] = 0.0,
           smlambda[GRAVITY_COEFF_MAXDEGREE+1u] = 0.0;

    unsigned int i,
                 n,
                 m;

    double deltaT                 = (JD - Gravity_reference_epoch)/DAYS_PER_YEAR,
           dUdr,      dUdr_n      = 1.0,
           dUdphi,    dUdphi_n    = 0.0,
           dUdlambda, dUdlambda_n = 0.0;

    double x = pos_ECEF[0],
           y = pos_ECEF[1],
           z = pos_ECEF[2],

           r      = sqrt (x*x + y*y + z*z),
           phi    = asin (z/r),
           lambda = atan2(y,x),

           slambda = sin(lambda),
           clambda = cos(lambda),

           rm1 = 1.0/r,
           rm2 = rm1*rm1,

           roR   = Gravity_Re*rm1,  xy2  = x*x + y*y,
           roR_n = roR,             nxy  = sqrt(xy2),
           oxy2  = 1/xy2,           onxy = 1/nxy;

    /* Casting the integers in a loop is quite expensive; so cast them once at
     * program start into memory, and reuse them */
    if (!grav_initialized)
    {
        for (n=0u; n<GRAVITY_COEFF_MAXDEGREE; ++n) {
            nn[n] = (double)n;
            mm[n] = nn[n];
        }
        grav_initialized = 1u;
    }

    /* Initialize output */
    for (i=0u; i<3u; ++i) (*gradient)[i] = 0.0;
    for (i=0u; i<9u; ++i) (*hessian)[i]  = 0.0;

    /* Fully normalized associated Legendre polynomials */
    ComputeLegendrePolys(phi, nmax, &Pnm);


    /* Compute partial derivatives of potential
     * ================================================================================ */

    /* pre-compute trig factors (recursion relations) */
    smlambda[0] = slambda;  smlambda[1] = 2.0*clambda*slambda;
    cmlambda[0] = clambda;  cmlambda[1] = 2.0*clambda*clambda - 1.0;

    for (m=2u; m<=nmax; ++m)
    {
        smlambda[m] = 2.0*clambda*smlambda[m-1u] - smlambda[m-2u];
        cmlambda[m] = 2.0*clambda*cmlambda[m-1u] - cmlambda[m-2u];
    }

    /* Compute double sum */
    for (n=2u; n<=nmax; ++n)
    {
        double dUdr_m      = 0.0,
               dUdphi_m    = 0.0,
               dUdlambda_m = 0.0;

        /* summation over m */
        for (m=0u; m<=n; ++m)
        {
            double fac1, Cnm = Gravity_Cnm[n][m],
                   fac2, Snm = Gravity_Snm[n][m];

            /* Take into account coefficient drifts */
            if (n < GRAVITY_DRIFT_MAXDEGREE && m < GRAVITY_DRIFT_MAXDEGREE) {
                Cnm += DeltaT * Gravity_Cnm_dot[n][m];
                Snm += DeltaT * Gravity_Snm_dot[n][m];
            }

            fac1 = Cnm * cmlambda[m] + Snm * smlambda[m],
            fac2 = Snm * cmlambda[m] - Cnm * smlambda[m];

            dUdr_m      += Pnm[n][m] * fac1;
            dUdphi_m    += (Pnm[n][m+1u] - z*onxy * mm[m]*Pnm[n][m]) * fac1;
            dUdlambda_m += mm[m]*Pnm[n][m] * fac2;
        }

        /* summation over n */
        roR_n       *= roR;
        dUdr_n      += roR_n * dUdr_m * (nn[n] + 1.0);
        dUdphi_n    += roR_n * dUdphi_m;
        dUdlambda_n += roR_n * dUdlambda_m;

    }

    /* Acceleration in spherical coordinates */
    dUdr      = -Gravity_GM*rm2 * dUdr_n;
    dUdphi    =  Gravity_GM*rm1 * dUdphi_n;
    dUdlambda =  Gravity_GM*rm1 * dUdlambda_n;

    /**/
    if (compute_gradient)
    {
        /* Special case for positions near the poles */
        if (fabs(atan2(z,nxy)) == PI_OVER_2)
        {
            (*gradient)[0] = 0.0;
            (*gradient)[1] = 0.0;
            (*gradient)[2] = z*dUdr*rm1;
        }

        /* All other locations: */
        else
        {
            /* Convert back to Planet centered/planet fixed */
            double f = dUdr*rm1,
                   g = f - z*dUdphi*rm2*onxy,
                   h = dUdlambda*oxy2;

            (*gradient)[0] = g*x - h*y;
            (*gradient)[1] = g*y + h*x;
            (*gradient)[2] = f*z + nxy*dUdphi*rm2;
        }
    }

    /* Compute second derivatives (Hessian) */
    if (compute_hessian)
    {
        /* TODO: eq 5-15, 8-30*/
    }
}







