/* ------------------------------------------------------------------------

   Program:  flat_full
   Purpose:  Divide a given echelle image by the white-lamp image
             (two-dimensional flat-fielding). This will correct the CCD
             image for any pixel-to-pixel variations.

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
   param_type par;
   fits_file xfits,flatfits;
   int seq;
   float *xptr,*fptr;

   if (argc != 3) get_error("Usage:  flat_full  <in> <out>");

   inform_title("Flat-field correction (2-D)");

   read_image_header(&xfits,argv[1]);
   get_spectrograph(&xfits,&xfits.head,spectrograph);
   get_parameters(&par,spectrograph);

   if (par.flat_method.id != FLAT_FULL)
   {
     get_system("cp %s.fit %s.fit",argv[1],argv[2]);
     inform("No flat-fielding required. File copy performed instead.");
     exit(0);
   }

   read_image(&xfits,argv[1]);
   read_image(&flatfits,"rflat");

   check_image_2d(&xfits);
   check_image_2d(&flatfits);

   if (xfits.npix[1] != flatfits.npix[1])
      get_error("Same number of pixels along X-axis expected "
                "for '%s' and '%s'!",xfits.file.path,flatfits.file.path);

   if (xfits.npix[2] != flatfits.npix[2])
      get_error("Same number of pixels along Y-axis expected "
                "for '%s' and '%s'!",xfits.file.path,flatfits.file.path);

   xptr=xfits.pix;
   fptr=flatfits.pix;
   for (seq=1;seq<=xfits.totpix;seq++,xptr++,fptr++)
   {
      if (*fptr == 0.0)
        *xptr = 0.0;
      else
        *xptr /= *fptr;
   }

   inform("%dx%d pixels processed.",xfits.npix[1],xfits.npix[2]);

   copy_fits_keyword(&flatfits.head,"WHITFILE",&xfits,&xfits.head);
   copy_fits_keyword(&flatfits.head,"WHITEJD",&xfits,&xfits.head);
   copy_fits_keyword(&flatfits.head,"WHITEIMG",&xfits,&xfits.head);

   update_fits_minmax(&xfits);
   write_image(&xfits,argv[2]);
   inform("Image '%s' created.",xfits.file.path);

   free_binary_data(&flatfits);

   check_memory();

   return(0);
}
