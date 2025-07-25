/* ------------------------------------------------------------------------

   Program:  frame_image
   Purpose:  Add a frame around a FITS image.

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
   fits_file xfits,yfits;
   int left,right,bottom,top;
   int nx,ny;
   float *src,*dest;
   int row;

   if (argc != 7) get_error
      ("Usage:  frame_image  <in> <out> <left> <right> <bottom> <top>");

   inform("Image framing:");

   left = atoi(argv[3]);
   right = atoi(argv[4]);
   bottom = atoi(argv[5]);
   top = atoi(argv[6]);

   inform("   left=%d right=%d bottom=%d top=%d",left,right,bottom,top);

   read_image(&xfits,argv[1]);
   check_image_2d(&xfits);

   nx = xfits.npix[1] + left + right;
   ny = xfits.npix[2] + bottom + top;

   create_image(argv[2],NULL,-32,2,nx,ny,1.0,1.0,1.0,1.0,1.0,1.0);
   read_image(&yfits,argv[2]);

   src = xfits.pix;
   dest = yfits.pix + bottom * nx + left;

   for (row=1;row<=xfits.npix[2];row++)
   {
      memmove(dest,src,xfits.npix[1]*sizeof(float));
      src += xfits.npix[1];
      dest += nx;
   }

   write_image(&yfits,argv[2]);
   inform("Image '%s' created.",yfits.file.path);
   free_binary_data(&xfits);

   check_memory();

   return(0);
}
