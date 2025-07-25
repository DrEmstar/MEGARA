/* ------------------------------------------------------------------------

   Program:  form_rflat
   Purpose:  Create a full-size (two-dimensional) flat-field image from
             a given white-lamp image.

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
   fits_file imgfits;
   char spectrograph[MAX_SPECNAME+1];
   param_type par;
   float *fptr,fmax;
   int seq;

   if (argc != 2) get_error("Usage:  form_rflat  <img>");

   setbuf(stderr,NULL);

   inform_title("Full-size flat-field image creation");

   read_image(&imgfits,argv[1]);
   check_real_image(&imgfits);
   get_spectrograph(&imgfits,&imgfits.head,spectrograph);
   get_parameters(&par,spectrograph);

   for (seq=1,fptr=imgfits.pix,fmax=0.0;seq<=imgfits.totpix;seq++,fptr++)
   {
     if (*fptr < par.flat_thres_full) *fptr = 0.0;
     if (*fptr > fmax) fmax = *fptr;
   }

   if (fmax == 0.0) get_error("All pixels below the threshold!");

   for (seq=1,fptr=imgfits.pix;seq<=imgfits.totpix;seq++,fptr++)
     *fptr = *fptr / fmax;

   write_keyword_textual(&imgfits,&imgfits.head,"WHITFILE",
                         imgfits.file.name,
                         "Image used for flat field creation");

   update_fits_minmax(&imgfits);

   write_image(&imgfits,"rflat");

   inform("Image '%s' created.",imgfits.file.path);

   check_memory();

   return(0);
}
