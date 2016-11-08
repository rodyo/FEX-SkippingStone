#include <stddef.h>
#include "Ephemerides_wrapper.h"
#include "jpleph.h"

/* data struct malloc'ed by JPL routines */
static void *ephem = NULL;


void
Ephem_initialize()
{
    char   *nams = NULL;
    double *vals = NULL;

    /* initialize ephem data structure */
    ephem = jpl_init_ephemeris(nams, vals);
}


int
Ephem_step(const double JD,
           unsigned char select_bodies[JPL_NUMDEFS],
           unsigned int origin,
           double (*statevec)[JPL_NUMDEFS][6])
{
    /* FIXME: (Rody Oldenhuis) put this in global struct or something */
    #define AU              (149597870.700)
    #define SECONDS_PER_DAY (86400.0)

    int err_code = 0;

    unsigned int i = 0u,
                 j = 0u;

    /* Compute relevant ephemeris */
    for (i=0u; i<JPL_NUMDEFS && err_code==0; ++i)
    {
        if (select_bodies[i])
        {
            err_code = jpl_pleph(ephem, JD, i+1u, origin,
                                 (*statevec)[i], 1);

            /* Convert units from AU to km, and AU/day to km/s */
            for (j=0u; j<3u; ++j) (*statevec)[i][j] *= AU;
            for (j=3u; j<6u; ++j) (*statevec)[i][j] *= AU/SECONDS_PER_DAY;
        }
        else {
            for (j=0u; j<6u; ++j)
                (*statevec)[i][j] = 0.0;
        }
    }

    /* error condition */
    if (err_code != 0)
    {
        /* Output error, and explicitly assign all-zero output values */
        for (i=0u; i<JPL_NUMDEFS; ++i)
            for (j=0u; j<6u; ++j)
                (*statevec)[i][j] = 0.0;
    }

    return err_code;
}


void
Ephem_cleanup()
{
    if (ephem != NULL)
        jpl_close_ephemeris(&ephem);
}

