#include "Ephemerides_wrapper.h"

#include "Utils.h"
#include "jpleph.h"

/* Setup logging for this S-function */
#if EXECUTING_BINARY

    #include "logger.h"
    #define DO_LOGGING (1u)
    static LOGTOKEN logToken = LOGTOKEN_NULL;

#else
    #define DO_LOGGING (0u)
    #include "simulatorControl.h"
#endif

/* data struct malloc'ed by JPL routines */
static void *ephem = NULL;

/* NOTE: (Rody Oldenhuis) 3 elements will be inserted into the output
 * signal for each celestial  body added; make sure the signal size
 * matches.*/
static const unsigned int
    bodies_of_interest[] = {JPL_MOON, JPL_SUN},
    num_bodies           = sizeof(bodies_of_interest)/sizeof(unsigned int);



void
Ephem_check_params(const unsigned int signalWidth)
{
    /* Check consistency */
    if (3u*num_bodies != signalWidth)
    {
        #if EXECUTING_BINARY
        {
            logMessage(logToken, LOG_ERROR,
                       "%s; line %d: number of bodies requested inconsistent "
                       "with size of output signal.",
                       __FILE__, __LINE__);
        }
        #else
        {
            sim_error("number of bodies requested inconsistent with size of "
                      "output signal.");
        }
        #endif
    }
}

void
Ephem_initialize(const char* fcnname)
{
    char   *nams = NULL;
    double *vals = NULL;

    /* initialize logging */
    #if DO_LOGGING
    {
        logToken = initialize_eventlog(fcnname);
        #if DEBUGGING
        {
            logMessage(logToken, LOG_DEBUG,
                       "Entered Ephemerides initialization function");
        }
        #endif
    }
    #else
    {
        UNREFERENCED_FORMAL_PARAMETER(fcnname);
    }
    #endif

    /* initialize ephem data structure */
    ephem = jpl_init_ephemeris(nams, vals);
    if (ephem == NULL)
    {
        #if DO_LOGGING
        {
            logMessage(logToken, LOG_ERROR,
                       "%s; line %d: Error reading ephemerides data.",
                       __FILE__, __LINE__);
        }
        #else
        {
            sim_error("Error reading ephemerides data.");
        }
        #endif
    }

    #if DO_LOGGING
    #if DEBUGGING
    {
        logMessage(logToken, LOG_DEBUG,
                   "Completed Ephemerides initialization");
    }
    #endif
    #endif
}


void
Ephem_step(const double JD,
           double *positions)
{
    /* FIXME: (Rody Oldenhuis) put this in global struct or something */
    static const double AU = 149597870.700;

    /* NOTE: (Rody Oldenhuis) must be 6; it should be able to hold
     *       velocities, even though they are not used */
    double statevec[JPL_NUMDEFS][6];

    int err_code = 0;

    unsigned int
        i = 0u;

    #if DO_LOGGING
    #if DEBUGGING
    {
        logMessage(logToken, LOG_DEBUG,
                   "Entered Ephemerides main function; starting computations");
    }
    #endif
    #endif

    /* Compute relevant ephemeris */
    for (i=0u; i<num_bodies && err_code==0; ++i) {
        err_code = jpl_pleph(ephem, JD, bodies_of_interest[i], JPL_EARTH,
                             statevec[bodies_of_interest[i]], 0);
    }

    if (err_code != 0)
    {
        /* Output error, and explicitly assign all-zero output values */
        #if DO_LOGGING
        {
            logMessage(logToken,LOG_ERROR,
                       "%s; line %d: Received error code %d from 'jpl_ephem'; "
                       "possibly out of range Julian date (%f)",
                       __FILE__, __LINE__, err_code, JD);
        }
        #else
        {
            sim_error("Received error code %d from 'jpl_ephem'; possibly out "
                      "of range Julian date (%f)",
                      err_code, JD);
        }
        #endif

        for (i=0u; i<3u*num_bodies; ++i)
            positions[i] = 0.0;

    }
    else
    {
        /* Convert units from AU to km and assign output signal */
        for (i=0u; i<num_bodies; ++i) {
            unsigned int j;
            for (j=0u; j<3u; ++j)
                positions[3u*i + j] = AU * statevec[bodies_of_interest[i]][j];
        }
    }

    #if DO_LOGGING
    #if DEBUGGING
    {
        logMessage(logToken, LOG_DEBUG,
                   "Computed all requested ephemerides; leaving "
                   "Ephemerides main function");
    }
    #endif
    #endif
}


void
Ephem_cleanup()
{
    #if DO_LOGGING
    #if DEBUGGING
    {
        logMessage(logToken, LOG_DEBUG,
                   "Cleaning up Ephemerides data");
    }
    #endif
    #endif

    if (ephem != NULL)
        jpl_close_ephemeris(&ephem);

    #if DO_LOGGING
    {
        #if DEBUGGING
        {
            logMessage(logToken, LOG_DEBUG,
                       "Cleaned up Ephemerides data; terminating");
        }
        #endif

        if (logToken != LOGTOKEN_NULL)
            terminate_eventlog(logToken);
    }
    #endif
}

