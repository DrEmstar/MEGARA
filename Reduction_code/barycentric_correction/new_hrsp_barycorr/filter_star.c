/* ------------------------------------------------------------------------

   Program:  filter_star
   Purpose:  Filter a stellar image

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
   file_name_type tmp_med,tmp_bkg,tmp_dif;
   fits_file imgfits;
   char spectrograph[MAX_SPECNAME+1];
   param_type par;

   if (argc < 3 || argc > 5)
      get_error("Usage:  filter_star  <in> <out>\n"
                "        filter_star  <in> <out> <cr>\n"
                "        filter_star  <in> <out> <cr> <bkg>");

   inform("");
   inform("Background and cosmic ray subtraction");
   inform("=====================================");

   read_image_header(&imgfits,argv[1]);
   get_spectrograph(&imgfits,&imgfits.head,spectrograph);

   get_parameters(&par,spectrograph);

   get_system("filter_median %s %s %d %d 0 ABS",
               argv[1],next_tmp_file_name(tmp_med),
               par.medfilt_rx,par.medfilt_ry);
   get_system("background_echelle %s %s",tmp_med,next_tmp_file_name(tmp_bkg));

   if (argc == 5)
   {
      get_system("cp %s.fit %s.fit",tmp_bkg,argv[4]);
      inform("Image './%s.fit' created.",argv[4]);
   }

   get_system("cosmic_rays %s %s %s %s",argv[1],tmp_med,tmp_bkg,argv[2]);

   if (argc >= 4)
   {
      inform_title("Cosmic ray image creation");

      subtract_image_file(argv[1],tmp_bkg,next_tmp_file_name(tmp_dif),NO_INFO);
      subtract_image_file(tmp_dif,argv[2],argv[3],YES_INFO);
      remove_tmp_file("%s.fit",tmp_dif);
   }

   remove_tmp_file("%s.fit",tmp_med);
   remove_tmp_file("%s.fit",tmp_bkg);

   check_memory();

   return(0);
}
