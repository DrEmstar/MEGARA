/* ------------------------------------------------------------------------

   Program:  rebin_image
   Purpose:  Rebin a FITS image by combining more pixels into one. Only
             an integer number of pixels can be specified. The input image
             can be cropped as well.

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
   int binsize,arow,acol,brow,bcol,wx,hx,wy,hy;
   int i,j,iy,jy,ix,jx;
   float pixsum,*ypix,*xpix;
   float gamma;

   if (argc != 8) get_error("Usage:  "
     "rebin_image  <in> <out> <binsize> <arow> <acol> <brow> <bcol>");

   inform("Image rebinning:");

   binsize = atoi(argv[3]);
   arow = atoi(argv[4]);
   acol = atoi(argv[5]);
   brow = atoi(argv[6]);
   bcol = atoi(argv[7]);

   inform("   binsize=%d arow=%d acol=%d brow=%d bcol=%d",
     binsize,arow,acol,brow,bcol);

   read_image(&xfits,argv[1]);
   check_image_2d(&xfits);

   if (acol < 1 || acol > xfits.npix[1]) get_error("Outside image area!");
   if (arow < 1 || arow > xfits.npix[2]) get_error("Outside image area!");
   if (bcol < 1 || bcol > xfits.npix[1]) get_error("Outside image area!");
   if (brow < 1 || brow > xfits.npix[2]) get_error("Outside image area!");

   wx = bcol - acol + 1;
   hx = brow - arow + 1;

   if (wx < 1 || hx < 1) get_error("Bad crop rectangle specification!");

   wy = wx / binsize;
   hy = hx / binsize;

   if (wy * binsize != wx || hy * binsize != hx) get_error
      ("The crop rectangle must contain an integer number of bins!");

   gamma = 1.0 / (float)(binsize * binsize);

   create_image(argv[2],NULL,-32,2,wy,hy,1.0,1.0,1.0,1.0,1.0,1.0);
   read_image(&yfits,argv[2]);

   ypix = yfits.pix;
   for (iy=1;iy<=hy;iy++)
   {
      for (jy=1;jy<=wy;jy++,ypix++)
      {
         pixsum = 0;
         for (i=1;i<=binsize;i++)
         {
            ix = arow + (iy - 1) * binsize + i - 1;
            for (j=1;j<=binsize;j++)
            {
               jx = acol + (jy - 1) * binsize + j - 1;
               xpix = xfits.pix + (ix - 1) * xfits.npix[1]  + jx - 1;
               pixsum += *xpix;
            }
         }
         *ypix = pixsum * gamma;
      }
   }

   write_image(&yfits,argv[2]);
   inform("Image '%s' created.",yfits.file.path);
   free_binary_data(&xfits);

   check_memory();

   return(0);
}
