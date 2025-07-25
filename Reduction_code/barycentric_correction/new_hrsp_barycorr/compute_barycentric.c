/* ------------------------------------------------------------------------

   Program:  compute_barycentric
   Purpose:  Compute the barycentric RV correction

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
   char spectrograph[MAX_SPECNAME+1];
   obs_type obs;
   fits_file xfits;
   int hipnum,usrnum;
   char objtype[MAX_OBJTYPE+1],hipcomp[4];
   char expdate[MAX_DATE+1],midtime[MAX_TIME+1];
   int day,month,year,uthr,utmin;
   double utsec;
   epoch_type epoch;
   topo_type topo;
   star_type astar,bstar;
   starlabel_type lab;
   double drv,djd,jd_bary;
   char labtxt[40],ratxt[40],dectxt[40];

   if (argc != 2) get_error("Usage:  compute_barycentric  <img>");

   inform_title("Barycentric correction");

   load_astro();
   load_cats();

   read_raw_image(&xfits,argv[1]);
   get_spectrograph(&xfits,&xfits.head,spectrograph);
   get_observatory(&obs,spectrograph);

   get_keyword_textual(&xfits,&xfits.head,"EXPDATE",expdate,MAX_DATE);
   get_keyword_textual(&xfits,&xfits.head,"MIDTIME",midtime,MAX_TIME);
   get_keyword_textual(&xfits,&xfits.head,"OBJTYPE",objtype,MAX_OBJTYPE);
   hipnum = get_keyword_integer(&xfits,&xfits.head,"HIPNUM");
   usrnum = get_keyword_integer(&xfits,&xfits.head,"USRNUM");
   get_keyword_textual(&xfits,&xfits.head,"HIPCOMP",hipcomp,3);

   inform("Object type: %s",objtype);

   if (strcmp(objtype,"STAR")!=0 && strcmp(objtype,"SKY")!=0)
     get_error("Invalid object type!");

   clear_star(&astar);
   clear_star(&bstar);
   memset(ratxt,0,sizeof(ratxt));
   memset(dectxt,0,sizeof(dectxt));
   drv = djd = 0.0;

   if (strcmp(objtype,"STAR") == 0)
   {
      clear_label(&lab);
      if (hipnum > 0)
        set_hip_label(&lab,hipnum,hipcomp);
      else
         set_usr_label(&lab,usrnum);
      if (lab.catnum < 1) get_error("Missing catalogue number (HIP or USR)!");
      inform("Star recognized: %s",print_star_label(&lab,labtxt));
      load_star(&lab,&astar);
   }

   inform("Date: %s",expdate);
   inform("Time: %s UTC",midtime);

   sscanf(expdate,"%d-%d-%d",&year,&month,&day);
   sscanf(midtime,"%d:%d:%lf",&uthr,&utmin,&utsec);
 
   set_epoch_utc(&epoch,day,month,year,uthr,utmin,utsec);
   set_topo(&topo,&obs.lon,&obs.lat,obs.height,&epoch);
 
   if (strcmp(objtype,"STAR") == 0)
   {
      star_corr(&astar,&bstar,&epoch,&topo,&drv,&djd);
      print_alpha(&bstar.ra.alpha,ratxt,4);
      print_delta(&bstar.dec.delta,dectxt,3);
   }
   else
      sun_corr(&epoch,&topo,&drv,&djd);

   jd_bary = epoch.utc.jd0 + (epoch.utc.jdf + djd);
 
   write_keyword_textual(&xfits,&xfits.head,"RA_OBS",ratxt,
      "Apparent right ascension (hour,min,sec)");
   write_keyword_textual(&xfits,&xfits.head,"DEC_OBS",dectxt,
      "Apparent declination (deg,min,sec)");
   write_keyword_double(&xfits,&xfits.head,"RVCORR",drv,
      "Barycentric RV correction (km/s)");
   write_keyword_double(&xfits,&xfits.head,"JDCORR",djd,
      "Barycentric JD correction (days)");
   write_keyword_double(&xfits,&xfits.head,"JD_BARY",jd_bary,
      "Corrected barycentric JD at mid-exposure");
 
   write_raw_image(&xfits,xfits.file.path);

   inform("Apparent RA: %s",ratxt); 
   inform("Apparent DEC: %s",dectxt); 
   inform("RV correction: %12.6f km/s",drv);
   inform("JD correction: %12.6f d",djd);

   free_astro();
   free_cats();

   check_memory();

   return(0);
}
