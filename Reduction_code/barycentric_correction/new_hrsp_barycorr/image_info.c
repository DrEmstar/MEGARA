/* ------------------------------------------------------------------------

   Program:  image_info
   Purpose:  Display some basic information about a given FITS image.

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
   fits_file fits;
   int nx,ny;
   double xa,xb,ya,yb;

   if (argc != 2) get_error("Usage:  image_info  <name>");

   read_image_header(&fits,argv[1]);
   check_image_2d(&fits);

   nx = fits.npix[1];
   ny = fits.npix[2];

   xa = world_coordinate(&fits,1,1.0);
   xb = world_coordinate(&fits,1,(double)nx);

   ya = world_coordinate(&fits,2,1.0);
   yb = world_coordinate(&fits,2,(double)ny);

   inform("Image size: %d x %d",fits.npix[1],fits.npix[2]);

   if (floor(xa)==xa && floor(xb)==xb)
     inform("X-axis: %1.0f..%1.0f",xa,xb);
   else
     inform("X-axis: %5.3f..%5.3f",xa,xb);

   if (floor(ya)==ya && floor(yb)==yb)
     inform("Y-axis: %1.0f..%1.0f",ya,yb);
   else
     inform("Y-axis: %5.3f..%5.3f",ya,yb);

   check_memory();

   return(0);
}
