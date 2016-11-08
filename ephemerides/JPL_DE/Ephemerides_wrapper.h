#ifndef EPHEMERIDES_WRAPPER_H
#define EPHEMERIDES_WRAPPER_H

void 
Ephem_check_params(const unsigned int signalWidth);

void 
Ephem_initialize(const char* fcnname);

void 
Ephem_step(const double JD, 
           double *positions);

void 
Ephem_cleanup(void);

#endif