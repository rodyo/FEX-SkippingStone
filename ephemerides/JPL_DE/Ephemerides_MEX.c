/*
 * JPL developement ephemerides
 *
 *
 * Authors:
 *     Rody Oldenhuis
 *
 *
 */

#include <stddef.h>
#include "mex.h"
#include "Ephemerides_wrapper.h"



static unsigned char initialized = 0u;



/* Gateway function */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    #define MEXNAME "JPL_DE405_Ephemerides"

    double *JD     = NULL,
           *bodies = NULL;

    double *statevectors = NULL;
    unsigned int output_size[3];

    unsigned int origin = JPL_SOLAR_SYSTEM_BARYCENTER;

    unsigned char select_bodies[JPL_NUMDEFS] = {0u};

    double statevecs[JPL_NUMDEFS][6] = {0u};

    unsigned int i,j,k,
                 bodyindex = 0u,
                 numtimes,
                 numbodies;


    /* Get inputs */
    JD        = mxGetPr(prhs[0]);
    numtimes  = mxGetNumberOfElements(prhs[0]);

    bodies    = mxGetPr(prhs[1]);
    numbodies = mxGetNumberOfElements(prhs[1]);

    if (nrhs == 3)
    {
        origin = *(unsigned int*)mxGetPr(prhs[2]);

        if (mxGetNumberOfElements(prhs[2]) != 1 || origin > JPL_NUMDEFS)
            mexErrMsgIdAndTxt(MEXNAME ":invalid_origin",
                              "Input argument \"origin\" must be given as a scalar integer between 1 and %d.",
                              JPL_NUMDEFS);
    }

    /* Check argument counts */
    if (nrhs < 2 || nrhs > 3) {
        mexErrMsgIdAndTxt(MEXNAME ":argcount",
                          MEXNAME " requires between two and three input arguments.");
    }
    if(nlhs > 1) {
        mexErrMsgIdAndTxt(MEXNAME ":argcount",
                          MEXNAME " produces only one output argument.");
    }

    /* Initialize output arguments */
    output_size[0] = numtimes;
    output_size[1] = 6u;
    output_size[2] = numbodies;

    plhs[0] = mxCreateNumericArray(3u,
                                   output_size,
                                   mxDOUBLE_CLASS,
                                   mxREAL);

    statevectors = mxGetPr(plhs[0]);


    /* Initialize the JPL routines */
    if (!initialized) {
        Ephem_initialize();
        initialized = 1u;
    }

    /* Initialize celestial body request vector */
    for (i=0u; i<numbodies; ++i)
    {
        unsigned int body = (unsigned int)bodies[i];

        if (body > JPL_NUMDEFS)
            mexErrMsgIdAndTxt(MEXNAME ":invalid_body",
                              "Target body %d not defined (maximum %d).",
                              i+1u, JPL_NUMDEFS);

        select_bodies[body-1u] = 1u;
    }

    /* Compute all requested ephemerides */
    for (i=0u; i<numtimes; ++i)
    {
        unsigned int errorcode = 0u;

        /* Get ephemerides for this epoch */
        errorcode = Ephem_step(JD[i],
                               select_bodies,
                               origin,
                               &statevecs);

        /* Handle any error */
        if (errorcode != 0)
            mexErrMsgIdAndTxt(MEXNAME ":jpl_error",
                              "An error occurred during evaluation of the ephemerides "
                              "at time %f; the error code was %d.",
                              JD[i], errorcode);

        /* Assign to outputs */
        bodyindex = 0u;
        for (k=0u; k<JPL_NUMDEFS; ++k)
        {
            if (select_bodies[k])
            {
                for (j=0u; j<6u; ++j)
                    statevectors[i + numtimes*(j + 6u*bodyindex)] = statevecs[k][j];

                ++bodyindex;
            }
        }
    }

}


