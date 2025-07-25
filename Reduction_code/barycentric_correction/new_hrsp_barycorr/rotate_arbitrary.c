/* ------------------------------------------------------------------------

   Program:  rotate_arbitrary
   Purpose:  Rotate a two-dimensional image by an arbitrary angle

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

double real_coordinate(int i,int n)
{
   int k;

   k = n >> 1;
   if ((k << 1) == n)
      return((double)(i-k) - 0.5);
   else
      return((double)(i-k-1));
}

double pix_coordinate(double x,int n)
{
   int k;

   k = n >> 1;
   if ((k << 1) == n)
      return(x + (double)k + 0.5);
   else
      return(x + (double)(k + 1));
}

double linear(double a,double b,double q)
{
   return(a + q *(b - a));
}

int main(int argc,char **argv)
{
   fits_file xfits,yfits;
   int row,col;
   float *src,*dest;
   double theta,ct,st,x,y,u,v,up,vp;
   int h,w,r,c,i,j;
   double rq,cq;
   float p[2][2];

   if (argc != 4)
      get_error("Usage:  rotate_arbitrary  <in> <out> <angle>");

   inform("Image rotation.");

   read_image(&xfits,argv[1]);
   check_image_2d(&xfits);

   w = xfits.npix[1];
   h = xfits.npix[2];

   yfits = xfits;
   set_fits_name(&yfits,argv[2]);
   allocate_whole_binary(&yfits);

   theta = deg2rad(atof(argv[3]));
   inform("   theta=%s degrees",argv[3]);

   ct = cos(theta);
   st = sin(theta);

   dest = yfits.pix;

   for (row=1;row<=h;row++)
   {
      y = real_coordinate(row,h);
      for (col=1;col<=w;col++)
      {
         x = real_coordinate(col,w);
         u = x * ct + y * st;
         v = y * ct - x * st;
         up = pix_coordinate(u,w);
         vp = pix_coordinate(v,h);
         c = (int)floor(up);
         r = (int)floor(vp);
         cq = up - (double)c;
         rq = vp - (double)r;
         src = xfits.pix + (r-1) * w + (c-1);
         for (i=0;i<2;i++)
         {
            for (j=0;j<2;j++)
            {
               if ((r+i>=1) && (r+i<=h) && (c+j>=1) && (c+j<=w))
                  p[i][j] = *(src + j);
               else
                  p[i][j] = 0.0;
            }
            src += w;
         }
         *dest = linear
            (linear(p[0][0],p[0][1],cq),linear(p[1][0],p[1][1],cq),rq);
         dest++;
      }
   }

   write_image(&yfits,argv[2]);
   inform("Image '%s' created.",yfits.file.path);
   free_binary_data(&xfits);

   check_memory();

   return(0);
}
