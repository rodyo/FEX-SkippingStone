#include <math.h>
#include "mex.h"

/*
    paretofront MEX-file
    original version by Yi Cao 
    (see MATLAB File Exchange, file 17251) 
 
    PARETO returns the logical Pareto membership of a set of points.

    synopsis:  front = paretofront(objMat)

    created by Yi Cao    
    y.cao@cranfield.ac.uk
    
    for compiling type 
        mex paretofront.c   
    in MATLAB command window
    
*/

/* definitions */
void paretofront(bool * front, double * M, unsigned int row, unsigned int col);

/* gateway function */
void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
    bool * front;              /* logical front membership array */
	double * M;                /* matrix of objective function values */
    double * N;                /* matrix of constraint violations */
	unsigned int rowM, colM;   /* matrix dimensions for M */
    unsigned int rowN, colN;   /* matrix dimensions for N */
	const int *dimsM, *dimsN;  /* simple array containing dimensions */
    int i, j;                  /* loop indices */
    
    /* check on number of inputs */
	if(nrhs == 0 || nlhs > 2){
	    printf("synopsis:   front = paretofront(objfun_vals) or \n paretofront(objfun_vals, con_vals)\n");
	    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	    return;
	}
	
    /* get pointers, dimensions */
	M     = mxGetPr(prhs[0]);    
	dimsM = mxGetDimensions(prhs[0]);
	rowM  = dimsM[0]; colM  = dimsM[1];
    if (nrhs == 2){
        N      = mxGetPr(prhs[1]);  
        dimsN  = mxGetDimensions(prhs[1]);
        rowN   = dimsN[0]; colN  = dimsN[1];
        /* only use the SUM of the constraint violations */        
        if (colN > 1){
            j=rowN;
            for (i=0; i<rowN; i++){
                j+=rowN; N[i] += N[j]; 
            }
        }
    }
	
	/* initialize output */
	plhs[0] = mxCreateLogicalMatrix(rowM, 1); /* initialize all as true */
	front   = (bool *) mxGetPr(plhs[0]);      /* cast output to boolean */
    
	/* main call */
    /* if (nrhs == 1) */
        /* unconstrained version by Yi Chao */
        paretofront(front, M, rowM, colM);
    /* else 
        /* constrained version partially by Yi Chao, partially by Rody Oldenhuis 
        paretofront_constrained(front, M, rowM, colM, N); */
}

/* unconstrained Pareto membership test */
void paretofront(bool * front, double * M, unsigned int row, unsigned int col)
{
    /* declarations */
    unsigned int t, s, i, j , j1, j2;  /* all loop indices */
    bool *checklist, coldominatedflag; /* working variables */
    
    /* initialize checklist */
    checklist = (bool *)mxMalloc(row*sizeof(bool));
    for(t = 0; t<row; t++) checklist[t] = true;
    
    /* loop through all rows */
    for(s = 0; s<row; s++) {
        t=s;
        if (!checklist[t]) continue;
        checklist[t]=false;
        coldominatedflag=true;
        for(i=t+1;i<row;i++) {
            if (!checklist[i]) continue;
            checklist[i]=false;
            for (j=0,j1=i,j2=t;j<col;j++,j1+=row,j2+=row) {
                if (M[j1] < M[j2]) {
                    checklist[i]=true;
                    break;
                }
            }
            if (!checklist[i]) continue;
            coldominatedflag=false;
            for (j=0,j1=i,j2=t;j<col;j++,j1+=row,j2+=row) {
                if (M[j1] > M[j2]) {
                    coldominatedflag=true;
                    break;
                }
            }
            if (!coldominatedflag) {     /*swap active index continue checking */
                front[t]=false;
                checklist[i]=false;
                coldominatedflag=true;
                t=i;
            }
        }
        front[t]=coldominatedflag;
        if (t>s) {
            for (i=s+1; i<t; i++) {
                if (!checklist[i]) continue;
                checklist[i]=false;
                for (j=0,j1=i,j2=t;j<col;j++,j1+=row,j2+=row) {
                    if (M[j1] < M[j2]) {
                        checklist[i]=true;
                        break;
                    }
                }
            }
        }
    }
    
    /* clear working variable */
    mxFree(checklist); 
}
