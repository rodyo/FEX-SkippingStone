/* Rody Oldenhuis - Summary of changes to original
 *
 * - removed dependency on external file; the DE405 definition file is
 *   comes in C-include form (binary), and is read here as such
 * - changed function signatures to comply with coding standards
 * - removed dependency on structure alignment (indexing out-of-bounds
 *   was used as a trick to quickly swap two adjacent structure members)
 * - Added high level/function level #ifdef/#endif mechanism to
 *   include/exclude specific functionalities, to easily improve
 *   code coverage
 */


 /* Functionality level */
/* #define JPL_USE_BYTE_SWAP */
/* #define JPL_USE_NUTATION */
/* #define JPL_USE_LIBRATION */





/***************************************************************************
*******                  JPLEPH.H                                  *********
****************************************************************************
**  This header file is used both by ASC2EPH and TESTEPH programs.        **
****************************************************************************
**  Written: May 28, 1997 by PAD   **  Last modified: June 23,1997 by PAD **
**  Modified further by Bill Gray,  Jun-Aug 2001                          **
****************************************************************************
**  PAD: dr. Piotr A. Dybczynski,          e-mail: dybol@phys.amu.edu.pl  **
**   Astronomical Observatory of the A.Mickiewicz Univ., Poznan, Poland   **
***************************************************************************/

/* By default,  in Windoze 32,  the JPL ephemeris functions are compiled
   into a DLL.  This is not really all that helpful at present, but may
   be useful to people who want to use the functions from languages other
   than C. */

#ifdef __cplusplus
extern "C" {
#endif

void* jpl_init_ephemeris(char *nam, double *val);

void jpl_close_ephemeris(void **ephem);

int  jpl_state(void *ephem, const double et, const int list[12],
               double pv[][6], double nut[4], const int bary);

int  jpl_pleph(void *ephem, const double et, const int ntarg,
               const int ncent, double rrd[], const int calc_velocity);

int make_sub_ephem(const void *ephem, const char *sub_filename,
                   const double start_jd, const double end_jd);

#ifdef __cplusplus
}
#endif

/* Following are constants used in        */
/* jpl_get_double( ) and jpl_get_long( ): */

#define JPL_EPHEM_START_JD               0
#define JPL_EPHEM_END_JD                 8
#define JPL_EPHEM_STEP                  16
#define JPL_EPHEM_N_CONSTANTS           24
#define JPL_EPHEM_AU_IN_KM              28
#define JPL_EPHEM_EARTH_MOON_RATIO      36
#define JPL_EPHEM_EPHEMERIS_VERSION    200
#define JPL_EPHEM_KERNEL_SIZE          204
#define JPL_EPHEM_KERNEL_RECORD_SIZE   208
#define JPL_EPHEM_KERNEL_NCOEFF        212
#define JPL_EPHEM_KERNEL_SWAP_BYTES    216

/* Some useful constants */
#define JPL_SUN                        11
#define JPL_MERCURY                    1
#define JPL_VENUS                      2
#define JPL_EARTH                      3
#define JPL_MOON                       10
#define JPL_MARS                       4
#define JPL_JUPITER                    5
#define JPL_SATURN                     6
#define JPL_URANUS                     7
#define JPL_NEPTUNE                    8
#define JPL_PLUTO                      9
#define JPL_SOLAR_SYSTEM_BARYCENTER    12
#define JPL_EARTH_MOON_BARYCENTER      13
#define JPL_NUTATIONS                  14
#define JPL_LIBRATIONS                 15

#define JPL_NUMDEFS                    15



