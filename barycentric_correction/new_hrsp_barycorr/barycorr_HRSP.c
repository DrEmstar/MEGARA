
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "vector.h"
#include "numeric.h"
#include "astro.h"
#include "config.h"
#include "fits.h"
#include "process.h"


#define OBJSKY    1
#define OBJSTAR   0
 
char *object_description[] = {"Object is a star","Object is the blue sky"};

int main(int argc,char **argv)
{
   int hipnum,usrnum;
   double day,month,year,uthr,utmin,utsec;
   epochtype epoch;
   topotype topo;
   startype star;
   double drv,djd;

 /*  if (argc != 3) get_error("Usage:  barycorr_HRSP  STARNAME EXPDATE MIDTIME"); */
/* 
   inform("Barycentric correction");
   inform("----------------------");
*/
   astro_initialize();
   load_earth();

/* from identify object -> getting the starname parameter (to be supplied as inputs in this version)
get_keyword_textual(&fits,&fits.head,"STARNAME",starname,MAX_STARNAME); */

/* setting the starname, expdate & midtime, the three required inputs */

/* 
starname = argv[1];
expdate = argv[2];
midtime = argv[3];
STARNAME= 'HD182640'           / Star name                                      
EXPDATE = '2009-08-18'         / Exposure date                                  
MIDTIME = '13:11:18.2'         / Universal time at middle of exposure           
*/

/* from identify object -> getting object info needed for get_barycentric_correction */
get_object_info(argv[1],&hipnum,&usrnum,&star);
/* hipnum = hipnum */
/* usrnum = usrnum */

/*
eq = jd_to_julian_year(star.equinox.tdt.JD);
ep = jd_to_julian_year(star.epoch.tdt.JD);
ra = star.ra.theta.frac;
dec = star.dec.delta.frac;

   star.equinox = julian_epoch(eq);
   star.epoch = julian_epoch(ep);
   star.ra = form_angle(getrad(ra));
   star.dec = form_angle(getrad(dec));
*/


/* 
possibly nicer to go like this after get_object_info:
*/                                                                                                                                     
   star.equinox = julian_epoch(jd_to_julian_year(star.equinox.tdt.JD));
   star.epoch = julian_epoch(jd_to_julian_year(star.epoch.tdt.JD));
   star.ra = form_angle(getrad(star.ra.theta.frac));
   star.dec = form_angle(getrad(star.dec.delta.frac));

   sscanf(argv[2],"%lf-%lf-%lf",&year,&month,&day);
   sscanf(argv[3],"%lf:%lf:%lf",&uthr,&utmin,&utsec);
 
   epoch = set_ut(day,month,year,uthr,utmin,utsec);
   topo = set_topocentre('+',11.0,21.0,51.6,'-',43.0,59.0,12.0,1027.0,&epoch);
 
   get_barycentric_correction(epoch,star,topo,&drv,&djd);

/*
   inform("RV correction: %8.3f km/s",drv);
   inform("JD correction: %12.5e day",djd);
*/ 

   inform("%8.3f",drv);                                                                                                                             

   free_earth();

   return(0);
}
