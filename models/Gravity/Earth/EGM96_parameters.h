#ifndef EGM96_PARAMETERS_H
#define EGM96_PARAMETERS_H

 
#define GRAVITY_COEFF_MAXDEGREE (361u)
#define GRAVITY_DRIFT_MAXDEGREE (3u)
 
double Gravity_P[GRAVITY_COEFF_MAXDEGREE+3u][GRAVITY_COEFF_MAXDEGREE+3u],
       Gravity_cmlambda[GRAVITY_COEFF_MAXDEGREE+1u],
       Gravity_smlambda[GRAVITY_COEFF_MAXDEGREE+1u];
 
extern const double Gravity_GM;
extern const double Gravity_Re;
 
extern const double Gravity_Cnm[GRAVITY_COEFF_MAXDEGREE][GRAVITY_COEFF_MAXDEGREE];
extern const double Gravity_Snm[GRAVITY_COEFF_MAXDEGREE][GRAVITY_COEFF_MAXDEGREE];
 
extern const double Gravity_reference_epoch;
extern const double Gravity_Cnm_dot[GRAVITY_DRIFT_MAXDEGREE][GRAVITY_DRIFT_MAXDEGREE];
extern const double Gravity_Snm_dot[GRAVITY_DRIFT_MAXDEGREE][GRAVITY_DRIFT_MAXDEGREE];
 
#endif
