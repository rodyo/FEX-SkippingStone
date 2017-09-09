#include <stddef.h>
#include <stdlib.h>  /* malloc */

#include "jpleph.h"
#include "jpl_int.h"
#include "mex.h"

#define PROGRAM_NAME "jpl_ephemerides_mex"

/* data struct malloc'ed by JPL routines */
static void *ephem = NULL;

/* Wrapper for JPL routines */
static
int
compute_ephemerides(const double         JD,
                    const unsigned int   origin,
                    const unsigned int   select_bodies[],
                    const unsigned int   num_bodies,                    
                          double       (*statevec)[][6])
{
    /* FIXME: (Rody Oldenhuis) put this in global struct or something */
    #define AU              (149597870.700)
    #define SECONDS_PER_DAY (86400.0)

    int err_code = 0;

    unsigned int i = 0u,
                 j = 0u;

    /* Compute relevant ephemeris */
    for (i=0u; i<num_bodies && err_code==0; ++i)
    {        
        err_code = jpl_pleph(ephem, 
                             JD,
                             select_bodies[i], 
                             origin,
                             (*statevec)[i],
                             1);

        /* Convert units from AU to km, and AU/day to km/s */
        for (j=0u; j<3u; ++j) (*statevec)[i][j] *= AU;
        for (j=3u; j<6u; ++j) (*statevec)[i][j] *= AU/SECONDS_PER_DAY;        
    }

    /* error condition */
    if (err_code != 0)
    {
        /* Output error, and explicitly assign all-zero output values */
        for (i=0u; i<num_bodies; ++i)
            for (j=0u; j<6u; ++j)
                (*statevec)[i][j] = 0.0;
    }

    return err_code;
}


/* Clean-up; when MEX file is cleared */
static
void 
close_fcn() 
{    
	if (ephem != NULL)
        jpl_close_ephemeris(&ephem);
}


/* Gateway function */
void 
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    
    /* Initialize
     * ==================================================================*/
    
    double JD = 0.0;
    
    unsigned int  origin          = 0u,
                  num_bodies      = 0u, 
                 *selected_bodies = NULL;
    
    /* first call: initialize ephem data structure */
    if (ephem == NULL)
    {
        char   *nams = NULL;
        double *vals = NULL;
        
        ephem = jpl_init_ephemeris(nams, vals);
    }
    
    /* last call; allow JPL to clean up after itself */
    mexAtExit( close_fcn );
        
    /* Check inputs & assign data 
     * ==================================================================*/
    
    /* Argin count */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt(PROGRAM_NAME ":argin_error", 
                          PROGRAM_NAME " takes exactly 3 inputs.");
    }

    if (nlhs > 2) {
        mexErrMsgIdAndTxt(PROGRAM_NAME ":argout_error", 
                          PROGRAM_NAME " returns at most 2 outputs.");
    }

    /* JD */
    if (!mxIsNumeric(prhs[0]) ||
         mxIsComplex(prhs[0]) || 
        !mxIsDouble(prhs[0])  ||    
        !mxIsScalar(prhs[0]))
    {
        mexErrMsgIdAndTxt(PROGRAM_NAME ":invalid_JD", 
                          "First argument to " PROGRAM_NAME " (Julian "
                          "date) must be a real scalar double.");
    }

    JD = mxGetScalar(prhs[0]);
    {
        struct jpl_eph_data *eph = (struct jpl_eph_data *)ephem;
        if (JD < eph->ephem_start || JD > eph->ephem_end)
            mexErrMsgIdAndTxt(PROGRAM_NAME ":JD_out_of_range", 
                              "Given JD lies outside of data range.");
    }
    
    /* Origin */
    if (!mxIsNumeric(prhs[1]) ||
         mxIsComplex(prhs[1]) || 
        !mxIsScalar(prhs[1]) || 
        !mxIsDouble(prhs[1]))
    {
        mexErrMsgIdAndTxt(PROGRAM_NAME ":invalid_origin", 
                          "Second argument to " PROGRAM_NAME " (origin"
                          "body) must be a real scalar integer.");            
    }

    {
        double origin_D = 0.0;

        origin_D = *(double*)mxGetData(prhs[1]);
        origin   = (unsigned int)origin_D;

        if ((double)origin != origin_D)             
            mexErrMsgIdAndTxt(PROGRAM_NAME ":noninteger_origin", 
                              "Origin must be integer.");

        if (origin == 0u || origin > (unsigned)JPL_NUMDEFS)             
            mexErrMsgIdAndTxt(PROGRAM_NAME ":invalid_origin", 
                              "Given origin (%u) not defined in JPL routines.", 
                              origin);
    }

    /* Body selection vector */            
    if (!mxIsNumeric(prhs[2]) ||
         mxIsComplex(prhs[2]) || 
        !mxIsDouble(prhs[2])  || 
        mxGetNumberOfDimensions(prhs[2]) != 2 || 
        (mxGetM(prhs[2]) != 1 && mxGetN(prhs[2]) != 1))
    {
        mexErrMsgIdAndTxt(PROGRAM_NAME ":invalid_body_selected", 
                          "Third argument to " PROGRAM_NAME " (body"
                          "selection vector) must be a vector of real "
                          "integers.");            

    }

    {
        double *selected_bodies_D = NULL;

        num_bodies        = mxGetNumberOfElements(prhs[2]);
        selected_bodies_D = mxGetData(prhs[2]);
        selected_bodies   = malloc(num_bodies * sizeof(unsigned int));

        unsigned int i = 0u;
        for (i=0u; i<num_bodies; ++i)
        {
            selected_bodies[i] = (unsigned int)selected_bodies_D[i];
            
            if ((double)selected_bodies[i] != selected_bodies_D[i]) 
            {
                free(selected_bodies);                
                mexErrMsgIdAndTxt(PROGRAM_NAME ":noninteger_origin", 
                                  "Body %u must be integer.", 
                                  i);
            }
            
            if (selected_bodies[i] == 0u || 
                selected_bodies[i] > (unsigned)JPL_NUMDEFS) 
            {
                free(selected_bodies);
                mexErrMsgIdAndTxt(PROGRAM_NAME ":invalid_body_selected", 
                                  "Selected body %i is not defined in "
                                  "JPL routines.", 
                                  i);
            }
        }   
    }    
    
    /* Initialize outputs 
     * ==================================================================*/
    
    if (nlhs == 1) {
        plhs[0] = mxCreateDoubleMatrix(num_bodies, 6, mxREAL);
    }
    else {
        plhs[0] = mxCreateDoubleMatrix(num_bodies, 3, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(num_bodies, 3, mxREAL);
    }
    
    /* Do the computation and assign outputs
     * ==================================================================*/
    { 
        int exit_code = 0;
        double (*statevec)[][6] = NULL;
        
        /* initialize temporary data container */
        statevec = malloc(num_bodies * sizeof(double[6]));
        if (statevec == NULL) {
            free(selected_bodies);
            mexErrMsgIdAndTxt(PROGRAM_NAME ":malloc_error", 
                              "Could not allocate sufficient memory.");
        }
        
        /* Call the wrapper */
        exit_code = compute_ephemerides(JD,
                                        origin,
                                        selected_bodies, 
                                        num_bodies,                    
                                        statevec);
        if (exit_code != 0) {
            free(statevec);
            free(selected_bodies);
            mexErrMsgIdAndTxt(PROGRAM_NAME ":jpl_error", 
                              "JPL routine returned error code %i.",
                              exit_code);
        }
        
        /* Assign outputs */
        {
            unsigned int i = 0u, 
                         j = 0u; 
            
            /* 1 output: 6-element statevector per body */            
            if (nlhs == 1)
            {
                double *sv = mxGetPr(plhs[0]);                
                for (i=0u; i<num_bodies; ++i)
                    for (j=0u; j<6u; ++j)
                        sv[i+num_bodies*j] = (*statevec)[i][j];
            }
            
            /* 2 outputs: position and velocity */
            else 
            {
                double *r = mxGetPr(plhs[0]),  
                       *v = mxGetPr(plhs[1]);
                
                for (i=0u; i<num_bodies; ++i) {
                    for (j=0u; j<3u; ++j)
                        r[i+num_bodies*j] = (*statevec)[i][j];
                    for (j=0u; j<3u; ++j)
                        v[i+num_bodies*j] = (*statevec)[i][j+3u];
                }                
            }
        }
        
        free(statevec);
        free(selected_bodies);
    }   
    
}


