/*
 EM2E        Convert eccentricity and mean anomaly to eccentric anomaly

 [E] = EM2E(e, M) converts the eccentricity [e] and mean anomaly [M] to 
 the corresponding eccentric anomaly [E]. The mean anomaly [M] should be 
 given in radians. This conversion is essentially a black-box function 
 that solves Kepler's equation for any arbitrary conic section; it works 
 correctly for all types of conic sections, and handles scalar/vector/
 matrix input intuitively. Note that the eccentric anomaly has no
 definition in the parabolic case; for parabolic cases, [E] = [theta] 
 (the true anomaly) is returned. 

 EM2E(e, M, tol) (with a third argument) uses a maximum error-tolerance 
 given by [tol]. The default is 1e-12.

 Algorithm: To calculate [theta] from [M], EM2E uses a carefully selected 
 first approximate root of Kepler's equation, followed by a Newton-Raphson 
 iteration scheme if the first estimate is not within the limits set by 
 [tol]. See [Seppo Mikkola, "A cubic approximation for Kepler's equation", 
 1987] for more details.

 See also eM2theta, eE2theta, eE2M.

 Authors
 .�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.�`�.
 Name       : Rody P.S. Oldenhuis
 E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
 Affiliation: Delft University of Technology

 Last edited 11/Nov/2009
*/

#include <math.h>
#include "mex.h"

/* define simple signum function for clarity */
double sign(double x){return(x>0)-(x<0);}
/* round towards zero */
double fix(double x){return((x<0)?(ceil(x)):(floor(x)));}

/* gateway function */
void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{     
    /* declarations */
    double *E;     /* initialize solution */
	double *Ms, M; /* mean anomalies */   
    double *es, e; /* eccentricities */
    double tol;    /* tolerance on Newton-Raphson iterations */
	unsigned int row, i; /* dimensions, loop index */
    const int *dimsM, *dimse; /* dimensions */
    /* variables in main loop */
    double alpha,beta,gamma,z,s,s2,Ep,x,y;
    unsigned int iters;
    unsigned int maxiters = 250;
    double pi = 3.141592653589793;
    
    /* empty initial output */
    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    
    /* check on number of inputs */
	if (nrhs < 2 || nrhs > 3){
	    mexErrMsgTxt("eM2E() requires arguments [e] and [M], and optional argument [tol].");         
        return;
	}
	
    /* get pointers, dimensions */
    es    = mxGetPr(prhs[0]);    
	Ms    = mxGetPr(prhs[1]);    
    dimse = mxGetDimensions(prhs[0]);
	dimsM = mxGetDimensions(prhs[1]);    
	row   = dimsM[0];
    
    /* error traps */
    if (row != dimse[0]){
        mexErrMsgTxt("Dimensions of [e] do not agree with those of [M].");  return;
    }
    if (dimsM[1] != 1){
        mexErrMsgTxt("[M] must be a columnvector.");  return;
    }
    if (dimse[1] != 1){
        mexErrMsgTxt("[e] must be a columnvector.");  return;
    }
    
    /* get tolerance or set default value */
    (nrhs > 2) ? (tol = *mxGetPr(prhs[2])) : (tol = 1e-12);
	
	/* initialize output */
	plhs[0] = mxCreateDoubleMatrix(row, 1, mxREAL); 
	E   = (double *) mxGetPr(plhs[0]);
    
    /* main procedure */
    for (i=0; i<row; i++){
        
        /* reset iteration counter */
        iters = 0;
        
        /* extract current eccentricity & mean anomaly */
        e = es[i];  M = Ms[i];
        
        /* quick exit */
        if (M == 0){E[i] = 0; continue;}  
        
        /* spend a few flops on finding a very accurate approximation 
         * to the root, as described in Mikkola's paper */
        gamma = (4.0*e + 0.5);
        alpha = fabs((1 - e)/gamma);
        
        /* elliptic case */
        if (e < 1.0){   
            
            /* make sure the M used for the calculation obeys  
             * -pi <= M <= pi */
            if ((M > pi) || (M < -pi)){
                M = fmod(M, 2.0*pi);
                if (M > +pi) M -= 2.0*pi;
                if (M < -pi) M += 2.0*pi;
            }
            
            /* quick exit */
            if (M == pi){E[i] = pi; continue;}
            if (M == -pi){E[i] = -pi; continue;}
            
            /* Continuation of Mikkola's scheme, elliptic case */
            beta  = M/2.0/gamma;
            z     = (beta + sign(beta)*sqrt(beta*beta + pow(alpha, 3.0)));
            z     = sign(z)*pow(fabs(z), 1.0/3.0);
            if (z == 0) z = 1e-100; /* prevent division by zero */
            s     = z - alpha/z;
            s    -= (0.078*pow(s,5.0))/(1 + e);
            E[i]  = M + e*(3.0*s - 4.0*pow(s,3.0));
            Ep    = 0; 
            
            /* Newton-Raphson iterations */
            while (fabs(Ep-E[i]) >= tol){
                iters += 1;  Ep = E[i];  
                E[i] += (M + e*sin(E[i]) - E[i]) / (1 - e*cos(E[i]));
                /* escape clause */
                if (iters >= maxiters){
                    mexWarnMsgTxt("Maximum iterations exceeded -- exiting.");
                    break;
                }
            } /* Newton-Raphson while loop */
            
            /* also put result in correct quadrant */
            E[i] += sign(Ms[i])*fix((fabs(Ms[i])+pi)/2.0/pi)*2.0*pi;
            
        }
        
        /* parabolic case */
        if (e == 1.0){
            /* [E] is not defined; return [theta]
             * (analytic solution to Barker's equation)
             * compute [th] */
            y = fabs(atan(1.0/3.0/M));
            x = atan(tan(pow(y/2.0, 1.0/3.0)));
            E[i] = 2*atan(2.0/tan(2.0*x));  
            /* make sure the quadrant is also correct */
            E[i] *= sign(M);
        }
        
        /* hyperbolic case */
        if (e > 1.0){            
            /* Continuation of Mikkola's scheme, hyperbolic case */
            beta = M/2.0/gamma;
            z    = (beta + sign(beta)*sqrt(beta*beta + pow(alpha, 3.0)));
            z    = sign(z)*pow(fabs(z), 1.0/3.0);
            if (z == 0) z = 1e-100; /* prevent division by zero */
            s    = z - alpha/z;
            s2   = s*s;
            s    += 0.071*pow(s,5.0) / (1 + 0.45*s2)/(1 + 4.0*s2)/e;
            E[i] = 3.0*log(s + sqrt(1 + s2));            
            Ep   = 0; 
            
            /* Newton-Raphson iterations */
            while (fabs(E[i]-Ep) >= tol){
                iters += 1;  Ep = E[i]; 
                E[i] += (E[i]- e*sinh(E[i]) + M) / (e*cosh(E[i]) - 1);
                /* escape clause */
                if (iters >= maxiters){
                    mexWarnMsgTxt("Maximum iterations exceeded -- exiting."); 
                    break;
                }
            } /* Newton-Raphson while loop */
            
        } /* type of conic section selector */
    } /* main loop */
} /* gateway routine */
