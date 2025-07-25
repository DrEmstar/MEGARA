/* ------------------------------------------------------------------------

   Program:  barycorr_hrsp
   Purpose:  Compute the velocity correction for a given star.
             This is the same as "velocity_correction", but with different
             input and output, to be used with the University of Canterbury
             Matlab code.

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
   char starstr[256];
   star_type astar,bstar;
   int day,month,year,tthr,ttmin;
   double ttsec;
   epoch_type epoch;
   topo_type topo;
   double vorb,vrot,dv,dt;

   if (argc != 4) get_error("Usage:  barycorr_hrsp  <star> <date> <time>");

   load_astro();
   load_cats();

   get_observatory(&obs,"hercules");

   clear_star(&astar);
   clear_star(&bstar);

   parse_star_name_err(argv[1],&lab);

   if (lab.kind != USR_STAR)
   {
      any_to_hip(&lab);
      print_star_label(&lab,starstr);
   }

   load_star(&lab,&astar);

   if (sscanf(argv[2],"%d-%d-%d",&year,&month,&day) != 3)
     get_error("Invalid date format!");
   if (sscanf(argv[3],"%d:%d:%lf",&tthr,&ttmin,&ttsec) != 3)
     get_error("Invalid time format!");

   set_epoch_tdt(&epoch,day,month,year,tthr,ttmin,ttsec);

   set_topo(&topo,&obs.lon,&obs.lat,obs.height,&epoch);

   astar.topo = topo;

   full_corr(&astar,&bstar,&epoch,&topo,&vorb,&vrot,&dv,&dt);

   printf("%8.3f\n",dv);

   free_astro();
   free_cats();

   check_memory();

   return(0);
}
