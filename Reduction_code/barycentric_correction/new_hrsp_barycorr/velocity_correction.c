/* ------------------------------------------------------------------------

   Program:  velocity_correction
   Purpose:  Compute the velocity correction for a given star.

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
   obs_type obs;
   starlabel_type lab;
   int day_arg;
   int usrnum;
   char starstr[256];
   char ratxt[40],dectxt[40];
   star_type astar,bstar;
   int rah,ram;
   double ras;
   char decz;
   int decd,decm;
   double decs;
   int ref,cen;
   int day,month,year,tthr,ttmin;
   double ttsec;
   epoch_type epoch;
   topo_type topo;
   double vorb,vrot,dv,dt;

   if (argc != 9 && argc != 16)
     get_error("Usage:  velocity_correction  "
       "<spec> <star> <day> <month> <year> <hr> <min> <sec>\n"
               "        velocity_correction  "
       "<spec> <rah> <ram> <ras> <decd> <decm> <decs> <ref> <cen> "
       "<day> <month> <year> <hr> <min> <sec>\n");

   load_astro();
   load_cats();

   get_observatory(&obs,argv[1]);

   clear_star(&astar);
   clear_star(&bstar);
   memset(ratxt,0,sizeof(ratxt));
   memset(dectxt,0,sizeof(dectxt));

   if (argc == 9)
   {
     day_arg = 3;
     inform("Star name: %s",argv[2]);
     parse_star_name_err(argv[2],&lab);

     if (lab.kind == USR_STAR)
     {
        usrnum = lab.catnum;
        inform("User catalogue number: USR %d",usrnum);
     }
     else
     {
       any_to_hip(&lab);
       print_star_label(&lab,starstr);
       inform("Hipparcos catalogue number: %s",starstr);
     }

     load_star(&lab,&astar);
   }
   else
   {
     day_arg = 10;
     rah = atoi(argv[2]);
     ram = atoi(argv[3]);
     ras = atof(argv[4]);
     if (*argv[5] == '-')
       decz = '-';
     else
       decz = '+';
     decd = abs(atoi(argv[5]));
     decm = atoi(argv[6]);
     decs = atof(argv[7]);

     ref = ICRS;
     if (strcasecmp(argv[8],"MEAN") == 0) ref = MEAN_EQ;
     if (strcasecmp(argv[8],"TRUE") == 0) ref = TRUE_EQ;

     cen = BARY_CENT;
     if (strcasecmp(argv[9],"GEO") == 0) cen = GEO_CENT;
     if (strcasecmp(argv[9],"TOPO") == 0) cen = TOPO_CENT;

     astar.ref = ref;
     astar.cent = cen;
     astar.pos = GEOMETRIC;
     set_angle_alpha(&astar.ra,rah,ram,ras);
     set_angle_delta(&astar.dec,decz,decd,decm,decs);
     astar.par = 0.0;
     astar.mura = 0.0;
     astar.mudec = 0.0;
     astar.rv = 0.0;
     astar.mass_ratio = 0.0;
     astar.component = CENTER_OF_MASS;
   }

   day = atoi(argv[day_arg]);
   month = atoi(argv[day_arg+1]);
   year = atoi(argv[day_arg+2]);

   tthr = atoi(argv[day_arg+3]);
   ttmin = atoi(argv[day_arg+4]);
   ttsec = atof(argv[day_arg+5]);

   set_epoch_tdt(&epoch,day,month,year,tthr,ttmin,ttsec);

   if (argc == 15)
   {
     astar.equinox = epoch;
     astar.epoch = epoch;
   }

   inform("Epoch: %d/%02d/%04d %d:%02d:%07.4f TT",
     epoch.tdt.date.day,epoch.tdt.date.month,epoch.tdt.date.year,
     epoch.tdt.time.hour,epoch.tdt.time.min,epoch.tdt.time.sec);
   inform("JD %17.9f",epoch.tdt.jd);

   inform("Location: %s",obs.name);

   set_topo(&topo,&obs.lon,&obs.lat,obs.height,&epoch);

   astar.topo = topo;

   full_corr(&astar,&bstar,&epoch,&topo,&vorb,&vrot,&dv,&dt);

   print_alpha(&bstar.ra.alpha,ratxt,4);
   print_delta(&bstar.dec.delta,dectxt,3);

   inform("ICRS topocentric position:");
   inform("RA =   %s",ratxt);
   inform("DEC = %s",dectxt);

   inform("Orbital correction: %10.5f km/s",vorb);
   inform("Diurnal correction: %10.5f km/s",vrot);
   inform("Total correction:   %10.5f km/s",dv);
   inform("Time correction:    %12.9f days",dt);

   free_astro();
   free_cats();

   check_memory();

   return(0);
}
