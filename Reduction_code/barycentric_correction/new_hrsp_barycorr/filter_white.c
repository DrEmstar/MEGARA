/* ------------------------------------------------------------------------

   Program:  filter_white
   Purpose:  Filter a white lamp image

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
   file_name_type tmp_med,tmp_bkg;
   fits_file imgfits;
   char spectrograph[MAX_SPECNAME+1];
   param_type par;

   if (argc < 3 || argc > 4)
      get_error("Usage:  filter_white  <in> <out>\n"
                "        filter_white  <in> <out> <bkg>");

   read_image_header(&imgfits,argv[1]);
   get_spectrograph(&imgfits,&imgfits.head,spectrograph);

   get_parameters(&par,spectrograph);

   get_system("filter_median %s %s %d %d 0 ABS",
               argv[1],next_tmp_file_name(tmp_med),
               par.medfilt_rx,par.medfilt_ry);
   get_system("background_echelle %s %s",tmp_med,next_tmp_file_name(tmp_bkg));

   if (argc == 4)
   {
      get_system("cp %s.fit %s.fit",tmp_bkg,argv[3]);
      inform("Image './%s.fit' created.",argv[3]);
   }

   subtract_image_file(argv[1],tmp_bkg,argv[2],YES_INFO);

   remove_tmp_file("%s.fit",tmp_med);
   remove_tmp_file("%s.fit",tmp_bkg);

   check_memory();

   return(0);
}
