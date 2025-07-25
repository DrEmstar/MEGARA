/* ------------------------------------------------------------------------

   Program:  normalize_spectrum
   Purpose:  Normalize a given stellar spectrum.

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
   fits_file xfits;
   int row,col;
   float *fptr,fmax;

   if (argc != 3) get_error("Usage:  normalize_spectrum  <in> <out>");

   inform_title("Echelle order normalization");

   read_image(&xfits,argv[1]);
   check_image_2d(&xfits);

   fptr = xfits.pix;
   for (row=1;row<=xfits.npix[2];row++,fptr+=xfits.npix[1])
   {
      for (col=0,fmax=*fptr;col<xfits.npix[1];col++)
         if (fptr[col] > fmax) fmax = fptr[col];
      if (fmax == 0.0) fmax=1.0;
      for (col=0;col<xfits.npix[1];col++) fptr[col] /= fmax;
   }

   inform("%d orders processed.",xfits.npix[2]);

   write_image(&xfits,argv[2]);
   inform("Image '%s' created.",xfits.file.path);

   check_memory();

   return(0);
}
