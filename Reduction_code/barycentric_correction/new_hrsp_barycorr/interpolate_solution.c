/* ------------------------------------------------------------------------

   Program:  interpolate_solution
   Purpose:  Interpolate the dispersion solution between two thorium images.

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
   fits_file afits,bfits,cfits;
   regre_block areg,breg;
   int i,h,m;
   double s,q;
   char midtime[MAX_TIME+1];
   time_type autc,butc,cutc;

   if (argc != 4)
      get_error("Usage:  interpolate_solution  <atab> <btab> <img>");

   inform_title("Dispersion solution interpolation");

   read_table_header(&afits,argv[1]);
   read_table_header(&bfits,argv[2]);
   read_raw_image(&cfits,argv[3]);

   read_regression_block(&afits,&afits.extend,"XREGRE",&areg);
   read_regression_block(&bfits,&bfits.extend,"XREGRE",&breg);

   if (areg.ma != breg.ma)
      get_error("Same number of coefficients expected!");

   get_keyword_textual(&afits,&afits.extend,"MIDTIME",midtime,MAX_TIME);
   sscanf(midtime,"%d:%d:%lf",&h,&m,&s);
   set_time_hms(&autc,h,m,s);

   get_keyword_textual(&bfits,&bfits.extend,"MIDTIME",midtime,MAX_TIME);
   sscanf(midtime,"%d:%d:%lf",&h,&m,&s);
   set_time_hms(&butc,h,m,s);

   get_keyword_textual(&cfits,&cfits.head,"MIDTIME",midtime,MAX_TIME);
   sscanf(midtime,"%d:%d:%lf",&h,&m,&s);
   set_time_hms(&cutc,h,m,s);

   inform("First solution at  %2d:%02d:%04.1f UTC",
     autc.hour,autc.min,autc.sec);
   inform("Second solution at %2d:%02d:%04.1f UTC",
     butc.hour,butc.min,butc.sec);
   inform("Stellar image at   %2d:%02d:%04.1f UTC",
     cutc.hour,cutc.min,cutc.sec);

   if (butc.frac < autc.frac) get_error("Thorium solutions out of sequence!");

   if (fabs(butc.frac-autc.frac) < 1.0e-12)
      q = 0;
   else
      q = (cutc.frac-autc.frac)/(butc.frac-autc.frac);


   inform("Interpolation argument q = %12.8f",q);

   for (i=1;i<=areg.ma;i++)
     areg.a[i] =  linint(autc.frac,areg.a[i],butc.frac,breg.a[i],cutc.frac);

   write_regression_block(&cfits,&cfits.head,"XREGRE",&areg);

   write_raw_image(&cfits,argv[3]);

   check_memory();

   return(0);
}
