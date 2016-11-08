#ifndef EPHEMERIDES_WRAPPER_H
#define EPHEMERIDES_WRAPPER_H

#include "jpleph.h"

void
Ephem_initialize(void);

int
Ephem_step(const double JD,
           unsigned char select_bodies[JPL_NUMDEFS],
           unsigned int origin,
           double (*statevec)[JPL_NUMDEFS][6]);

void
Ephem_cleanup(void);

#endif