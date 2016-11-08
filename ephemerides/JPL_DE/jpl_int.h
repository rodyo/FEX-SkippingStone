/* Right now,DEs 403 and 405 have the maximum kernel size of 2036.    */
/* This value may need to be updated the next time JPL releases a new DE: */

#define MAX_KERNEL_SIZE 2036

/***** THERE IS NO NEED TO MODIFY THE REST OF THIS SOURCE *********/


/* A JPL binary ephemeris header contains five doubles and */
/* (up to) 41 long integers,  so:                          */
#define JPL_HEADER_SIZE (5*sizeof(double) + 41*sizeof(long))

#pragma pack(1)

struct jpl_eph_data 
{
   double ephem_start, ephem_end, ephem_step;
   long   ncon;
   double au;
   double emrat;
   long   ipt[13][3];
   long   ephemeris_version;
   long   kernel_size, recsize, ncoeff;
   long   swap_bytes;
   long   curr_cache_loc;
   double pvsun[6];
   double *cache;
   void   *iinfo;
   
   unsigned char *fileData;
   unsigned int  fileDataLength;
};

struct interpolation_info
{
   double pc[18], vc[18], twot;
   int np, nv;
};

#pragma pack()
