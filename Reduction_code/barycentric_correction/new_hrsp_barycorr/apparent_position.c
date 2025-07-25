/* ------------------------------------------------------------------------

   Program:  apparent_position
   Purpose:  Compute the apparent position of a star.

   ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "datetime.h"
#include "vector.h"
#include "numeric.h"
#include "nutation.h"
#include "astro.h"
#include "config.h"
#include "fits.h"
#include "process.h"

int main(int argc,char **argv)
{
   starlabel_type lab;
   int usrnum;
   char starstr[256];
   char ratxt[40],dectxt[40];
   star_type astar,bstar;
   int day,month,year,tthr,ttmin;
   double ttsec;
   epoch_type epoch;
   int ref,pos;

   if (argc != 10) get_error("Usage:  apparent_position  <star> "
     "<day> <month> <year> <hr> <min> <sec> <ref> <pos>");

   load_astro();
   load_cats();

   clear_star(&astar);
   clear_star(&bstar);
   memset(ratxt,0,sizeof(ratxt));
   memset(dectxt,0,sizeof(dectxt));

   inform("Star name: %s",argv[1]);
   parse_star_name_err(argv[1],&lab);

   if (lab.kind == USR_STAR)
   {
      usrnum = lab.catnum;
      inform("User catalogue number: USR %d",usrnum);
   }
   else
   {
     any_to_hip(&lab);
     print_star_label(&lab,starstr);
     inform("Hipparcos catalog number: %s",starstr);
   }

   load_star(&lab,&astar);

   day = atoi(argv[2]);
   month = atoi(argv[3]);
   year = atoi(argv[4]);

   tthr = atoi(argv[5]);
   ttmin = atoi(argv[6]);
   ttsec = atof(argv[7]);

   set_epoch_tdt(&epoch,day,month,year,tthr,ttmin,ttsec);

   inform("Epoch: %d/%02d/%04d %d:%02d:%07.4f TT",
     epoch.tdt.date.day,epoch.tdt.date.month,epoch.tdt.date.year,
     epoch.tdt.time.hour,epoch.tdt.time.min,epoch.tdt.time.sec);
   inform("JD %17.9f",epoch.tdt.jd);

   ref = ICRS;
   if (strcasecmp(argv[8],"MEAN") == 0) ref = MEAN_EQ;
   if (strcasecmp(argv[8],"TRUE") == 0) ref = TRUE_EQ;

   if (strcasecmp(argv[9],"APP") == 0)
     pos = APPARENT;
   else
     pos = GEOMETRIC;
   
   switch(ref)
   {
     case MEAN_EQ:
       inform("Mean equator and equinox of date");
       break;
     case TRUE_EQ:
       inform("True equator and equinox of date");
       break;
     default:
       inform("International Celestial Reference System (ICRS)");
   }

   if (pos == APPARENT)
     inform("Apparent position:");
   else
     inform("Geometric position:");

   move_star(&astar,&bstar,ref,GEO_CENT,pos,&epoch,&epoch,NULL);

   print_alpha(&bstar.ra.alpha,ratxt,4);
   print_delta(&bstar.dec.delta,dectxt,3);

   inform("RA =   %s",ratxt);
   inform("DEC = %s",dectxt);

   free_astro();
   free_cats();

   check_memory();

   return(0);
}
