/* ------------------------------------------------------------------------

   Program:  flat_extracted
   Purpose:  Divide a given extracted stellar spectrum by the extracted
             white-lamp spectrum which has optionally been smoothed to
             eliminate any high frequency variations between pixels.
             This will also correct the curvature of the continuum in the
             stellar spectrum due to the echelle blaze function.

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
   int row,order,flatrow,i;
   float *xptr,*flatptr;

   if (argc != 3) get_error("Usage:  flat_extracted  <in> <out>");

   inform_title("Flat-field correction (1-D)");

   read_image_header(&xfits,argv[1]);
   get_spectrograph(&xfits,&xfits.head,spectrograph);
   get_parameters(&par,spectrograph);

   if (par.flat_method.id != FLAT_EXTRACTED)
   {
     get_system("cp %s.fit %s.fit",argv[1],argv[2]);
     inform("No flat-fielding required. File copy performed instead.");
     exit(0);
   }

   read_image(&xfits,argv[1]);
   read_image(&flatfits,"eflat");

   check_image_2d(&xfits);
   check_image_2d(&flatfits);

   if (xfits.npix[1] != flatfits.npix[1])
      get_error("Same number of pixels along X-axis expected "
                "for '%s' and '%s'!",xfits.file.path,flatfits.file.path);

   for (row=1;row<=xfits.npix[2];row++)
   {
      order = (int)world_coordinate(&xfits,2,(double)row);
      flatrow = pixel_coordinate(&flatfits,2,(double)order);
      if (flatrow < 1 || flatrow > flatfits.npix[2])
      {
         inform("Order %d cannot be flat-fielded "
                "(order not found in '%s')!",order,flatfits.file.path);
      }
      else
      {
         xptr = xfits.pix + (row-1)*xfits.npix[1];
         flatptr = flatfits.pix + (flatrow-1)*flatfits.npix[1];
         for (i=0;i<xfits.npix[1];i++)
         {
            if (flatptr[i] == 0.0)
              xptr[i] = 0.0;
            else
              xptr[i] /= flatptr[i];
         }
      }
   }

   inform("%d orders processed.",xfits.npix[2]);

   copy_fits_keyword(&flatfits.head,"WHITFILE",&xfits,&xfits.head);
   copy_fits_keyword(&flatfits.head,"WHITEJD",&xfits,&xfits.head);
   copy_fits_keyword(&flatfits.head,"WHITEIMG",&xfits,&xfits.head);

   write_image(&xfits,argv[2]);
   inform("Image '%s' created.",xfits.file.path);

   free_binary_data(&flatfits);

   check_memory();

   return(0);
}
