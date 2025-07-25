/* ------------------------------------------------------------------------

   Program:  filter_median
   Purpose:  Apply a median filter to a given CCD image

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
   int rx,ry,dx,dy,npix,nmed;
   double thres;
   float *box,*dest,*src,*xpix,*ypix,med;
   int row,col,i,pold,pnew;
   int u,v,usize,vsize;

   if (argc != 7)
      get_error("Usage:  filter_median  <in> <out> <rx> <ry> <thres> <mthd>");

   setbuf(stderr,NULL);

   inform_title("Median filter (2-D)");

   rx = atoi(argv[3]);
   ry = atoi(argv[4]);
   thres = atof(argv[5]);
   if (strcasecmp(argv[6],"ABS") && strcasecmp(argv[6],"REL"))
      get_error("filter_median: Bad method specification!");

   read_image(&xfits,argv[1]);
   check_image_2d(&xfits);

   dx = 2*rx+1;
   dy = 2*ry+1;
   npix = dx*dy;
   nmed = (npix+1)/2;

   inform("Box size: %d x %d",dx,dy);

   usize = xfits.npix[1] - 2*rx;
   vsize = xfits.npix[2] - 2*ry;

   if (rx < 0 || usize < 1)
      get_error("filter_median: Bad filter specification along X-axis!");

   if (ry < 0 || vsize < 1)
      get_error("filter_median: Bad filter specification along Y-axis!");

   if (dx == 0 && dy == 0)
      get_error("filter_median: Bad filter specification!");

   yfits = xfits;
   set_fits_name(&yfits,argv[2]);
   allocate_whole_binary(&yfits);
   memmove(yfits.bin,xfits.bin,xfits.totpix*sizeof(float));

   box = fvector(npix);

   fprintf(stderr,"Image filter:   0%%");
   for (row=ry+1,v=1;v<=vsize;row++,v++)
   {
      pnew = nint((v-1)*100.0/(double)(vsize-1));
      if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
      xpix = xfits.pix + (row-1)*xfits.npix[1]+rx;
      ypix = yfits.pix + (row-1)*xfits.npix[1]+rx;
      for (col=rx+1,u=1;u<=usize;col++,xpix++,ypix++,u++)
      {
         src = xfits.pix + (row-ry-1)*xfits.npix[1] + col-rx-1;
         for (i=row-ry,dest=&box[1];i<=row+ry;i++,dest+=dx,src+=xfits.npix[1])
            memmove(dest,src,dx*sizeof(float));
         med = fselect(nmed,npix,box);
         if (strcasecmp(argv[6],"ABS") == 0)
         {
            if (fabs(*xpix-med) > thres) *ypix = med;
         }
         else
         {
            if (fabs(*xpix-med) > thres*sqrt(*xpix)) *ypix = med;
         }
         if (u == 1)
            for (i=1;i<=rx;i++) ypix[-i] = *ypix;
         else
            if (u == usize) for (i=1;i<=rx;i++) ypix[i] = *ypix;
      }
      if (v == 1)
      {
         src = yfits.pix + ry*yfits.npix[1];
         for (i=1,dest=src-yfits.npix[1];i<=ry;i++,dest-=yfits.npix[1])
            memmove(dest,src,yfits.npix[1]*sizeof(float));
      }
      else
      {
         if (v == vsize)
         {
            src = yfits.pix + (yfits.npix[2]-ry-1)*yfits.npix[1];
            for (i=1,dest=src+yfits.npix[1];i<=ry;i++,dest+=yfits.npix[1])
               memmove(dest,src,yfits.npix[1]*sizeof(float));
         }
      }
   }
   fprintf(stderr,"\n");

   free_fvector(box);

   write_image(&yfits,argv[2]);
   inform("Image '%s' created.",yfits.file.path);

   free_binary_data(&xfits);

   check_memory();

   return(0);
}
