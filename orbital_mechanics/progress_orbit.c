#include <math.h>
#include "mex.h"
#include <matrix.h>

/* gateway function */
void mexFunction(int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[])
{ 
	double x,y,z;                        /* initial positions */
    double xdot,ydot,zdot;               /* initial velocities */
    double *dts, *statevec;               /* time steps & statevectors */
    double muC;                          /* std. grav. parameter Sun */
    char *tunits = "days";               /* time units */
	unsigned int numdts, i;              /* number of time steps, loop indicess */
    const int *dimst;                    /* dimensions */
    const double pi = 3.141592653589793115997963; /* pi  */
    double *state, *exitflag;            /* output arguments */
    double *xm,*ym,*zm,*xdotm,*ydotm,*zdotm; /* output arguments, different format */
    mxArray *xx,*yy,*zz,*xxdot,*yydot,*zzdot;/* output arguments, different format */
    double nu0,r1m,beta,P,DeltaU,t,dt,bu,q;  /* used in time loop */
    double f,g,G,F,A,B,cFracPrev,cFrac;      /* used in time loop & cont. fraction */
    int n,k,d,l,iter,q_corrected;            /* used in time loop & cont. fraction */
    double U0w2,U1w2,U,U0,U1,U2,U3,r,u,DeltaUFactor,deltaT;/* used in time loop */
    bool qisbad = false, t_in_days = true;   /* escape clause and time units */
    mxArray *vecin[7], *vecout[7];            /* used for MATLAB functions (in case of failure) */
            
    /* narg check */
    if ((nrhs > 9) || (nrhs < 3))
        mexErrMsgTxt("PROGRESS_ORBIT() requires more than 3, but less than 9 input arguments.");
    
    /* parse input */
    dts    = mxGetPr(prhs[0]);
    dimst  = mxGetDimensions(prhs[0]);
    numdts = dimst[0];
    if (nrhs <= 4){
        statevec = mxGetPr(prhs[1]);        
        x = statevec[0];    xdot = statevec[3];        
        y = statevec[1];    ydot = statevec[4];
        z = statevec[2];    zdot = statevec[5];
        muC = *mxGetPr(prhs[2]);        
        if (nrhs == 4)
            tunits = mxArrayToString(prhs[3]);
    }
    if (nrhs >= 8){        
        x = *mxGetPr(prhs[1]);  xdot = *mxGetPr(prhs[4]);
        y = *mxGetPr(prhs[2]);  ydot = *mxGetPr(prhs[5]);  
        z = *mxGetPr(prhs[3]);  zdot = *mxGetPr(prhs[6]);
        muC = *mxGetPr(prhs[7]);
        if (nrhs == 9)
            tunits = mxArrayToString(prhs[8]);
    }
    
    /* initial output is empty */
    plhs[0] = mxCreateDoubleMatrix(numdts, 6, mxREAL);
    state = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(numdts, 1, mxREAL); /* exitflag */
    exitflag = mxGetPr(plhs[1]); 
        
    /* process time unit; 
       times are given in days, unless specified otherwise */
    if (!strcmp(tunits, "seconds"))
        t_in_days = false;
    else if (strcmp(tunits, "days")>0){
        mexErrMsgTxt("Time units must be specified as either 'seconds' or 'days'.");
        return;
    }
    
    /* initialize */
    nu0  = x*xdot+y*ydot+z*zdot;
    r1m  = sqrt(x*x+y*y+z*z);
    beta = 2*muC/r1m - (xdot*xdot+ydot*ydot+zdot*zdot);
    /* initialize period effects */
    if (beta > 0){
        P = 2*pi*muC*pow(beta, -3.0/2.0);
        DeltaUFactor = 2*pi*pow(beta, -5.0/2.0);
    }        
    
    /* main call: loop through all elements in [dts] */
    for(i = 0; i < numdts; i++){
        
        /* extract relevant time step (and convert to seconds) */
        dt = dts[i]; if (t_in_days) dt = dt * 86400;
        
        /* possibly quick exit */
        if (dt == 0){
            state[0*numdts+i] = x;   state[3*numdts+i] = xdot;
            state[1*numdts+i] = y;   state[4*numdts+i] = ydot;
            state[2*numdts+i] = z;   state[5*numdts+i] = zdot;
            exitflag[i] = 1;   continue;
        }
        
        /* period effects */
        DeltaU = 0;
        if (beta > 0){            
            n = floor((dt + P/2.0 - 2.0*nu0/beta)/P);
            DeltaU = n*DeltaUFactor;
        }
        
        /* loop until convergence of the time step */
        u = 0; t = 0; iter = 0; qisbad = false; deltaT = t-dt;
        while (fabs(t-dt) > 1){ /* one second tolerance seems fine
            
            /* increase iterations */
            iter += 1;
            
            /* compute q 
               NOTE: [q] may not exceed 1/2. In principle, this will never 
               occur, but the iterative nature of the procedure can bring 
               it above 1/2 for some iterations. */
            bu = beta*u*u;
            q  = bu/(1 + bu); 
            
            /* escape clause;
               The value for [q] will almost always stabilize to a value less 
               than 1/2 after a few iterations, but NOT always. In those
               cases, just use repeated coordinate transformations */
            if ((iter > 25) || (q >= 1)) {qisbad = true; break;}
            
            /* evaluate continued fraction
               (when q < 1, it *always* converges) */
            A = 1.0; B = 1.0; n = 0; k = -9; d = 15; l = 3;
            cFrac = 1.0, cFracPrev = 2.0;
            while (fabs(cFrac-cFracPrev) > 1e-14){
                k = -k;                 l += 2;
                d += 4*l;               n += (1+k)*l;
                A = d/(d - n*A*q);      B *= (A-1);
                cFracPrev = cFrac;      cFrac += B;                
            }
            
            /* continue time loop */
            U0w2   = 1.0 - 2.0*q;
            U1w2   = 2.0*(1-q)*u;
            U      = 16.0/15.0*pow(U1w2, 5.0)*cFrac + DeltaU;
            U0     = 2.0*U0w2*U0w2-1;
            U1     = 2.0*U0w2*U1w2;
            U2     = 2.0*U1w2*U1w2;
            U3     = beta*U + U1*U2/3;
            r      = r1m*U0 + nu0*U1 + muC*U2;
            t      = r1m*U1 + nu0*U2 + muC*U3;
            deltaT = t-dt;
            /* Newton-Raphson method works most of the time, but is 
               not too stable; the method fails far too often for my
               liking...
               u -= deltaT/4.0/(1-q)/r; 
               Halley's method is much better in that respect. Working 
               out all substitutions and collecting terms gives the 
               following simplification: */
            u -= deltaT/(1-q)/(4.0*r + deltaT*beta*u);
        } /* time loop */
        
        /* do it the slow way if state transition matrix fails for some q
           (repeated coordinate transformations) */
        if (qisbad){   
            vecin[0] = mxCreateDoubleScalar(x);     vecin[3] = mxCreateDoubleScalar(xdot);
            vecin[1] = mxCreateDoubleScalar(y);     vecin[4] = mxCreateDoubleScalar(ydot);
            vecin[2] = mxCreateDoubleScalar(z);     vecin[5] = mxCreateDoubleScalar(zdot);           
            vecin[6] = mxCreateDoubleScalar(muC);   vecout[6] = mxCreateDoubleScalar(muC);          
            /* call CART2KEP() */
            mexCallMATLAB(6, vecout, 7, vecin, "cart2kep");
            /* adjust M             */
            *mxGetPr(vecout[5]) += dt * sqrt(muC/pow(fabs(*mxGetPr(vecout[0])),3.0));
            /* call KEP2CART() */
            mexCallMATLAB(6, vecin, 7, vecout, "kep2cart");
            /* insert results in state[] */
            state[0*numdts+i] = *mxGetPr(vecin[0]);    state[3*numdts+i] = *mxGetPr(vecin[3]);
            state[1*numdts+i] = *mxGetPr(vecin[1]);    state[4*numdts+i] = *mxGetPr(vecin[4]);
            state[2*numdts+i] = *mxGetPr(vecin[2]);    state[5*numdts+i] = *mxGetPr(vecin[5]);          
            /* this DOES mean failure */
            exitflag[i] = -1;
            
        /* use transition matrix if all went well */
        } else {
            /* Kepler solution */
            f = 1.0 - muC/r1m*U2;   F = -muC*U1/r/r1m;
            g = r1m*U1 + nu0*U2;    G = 1.0 - muC/r*U2;
            /* carry out matrix multiplication by hand */
            state[0*numdts+i] = x*f + xdot*g;    state[3*numdts+i] = x*F + xdot*G;
            state[1*numdts+i] = y*f + ydot*g;    state[4*numdts+i] = y*F + ydot*G;
            state[2*numdts+i] = z*f + zdot*g;    state[5*numdts+i] = z*F + zdot*G;
            /* all went fine */
            exitflag[i] = 1;  
        }   
        
    } /* loop through all time steps */
    
    /* assign output arguments */
    if (nlhs >= 3){
        xx    = mxCreateDoubleMatrix(numdts,1,mxREAL);   xm = mxGetPr(xx);        
        yy    = mxCreateDoubleMatrix(numdts,1,mxREAL);   ym = mxGetPr(yy);        
        zz    = mxCreateDoubleMatrix(numdts,1,mxREAL);   zm = mxGetPr(zz);
        xxdot = mxCreateDoubleMatrix(numdts,1,mxREAL);   xdotm = mxGetPr(xxdot);
        yydot = mxCreateDoubleMatrix(numdts,1,mxREAL);   ydotm = mxGetPr(yydot);
        zzdot = mxCreateDoubleMatrix(numdts,1,mxREAL);   zdotm = mxGetPr(zzdot);                
        for (i = 0; i < numdts; i++){
            xm[i] = state[i+0*numdts];    xdotm[i] = state[i+3*numdts];
            ym[i] = state[i+1*numdts];    ydotm[i] = state[i+4*numdts];
            zm[i] = state[i+2*numdts];    zdotm[i] = state[i+5*numdts];
        }
        /* exitflag first */
        if (nlhs == 7)
            plhs[6] = plhs[1];
        /* the rest later */
        plhs[0] = xx;  plhs[1] = yy;  plhs[2] = zz;  
        if (nlhs >= 4){
            plhs[3] = xxdot;
            if (nlhs >= 5){
                plhs[4] = yydot;
                if (nlhs >= 6){
                    plhs[5] = zzdot;
        }}}
    } /* different output format */
        
} /* end of gateway function */
