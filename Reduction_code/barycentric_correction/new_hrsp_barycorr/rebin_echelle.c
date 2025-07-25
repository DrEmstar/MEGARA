/* ------------------------------------------------------------------------

   Program:  rebin_echelle
   Purpose:  Rebin a given spectrum from pixel space to wavelength space

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

double xord;
double *b;
int mb;
spec_type spec;

double func_lin(double w)
{
   double mlnorm;

   mlnorm = normalized_mlambda(&spec,xord*w);
   return(compute_polynomial(b,mb,mlnorm));
}

double func_log(double lw)
{
   double mlnorm;

   mlnorm = normalized_mlambda(&spec,xord*exp(lw));
   return(compute_polynomial(b,mb,mlnorm));
}

int main(int argc,char **argv)
{
   ccd_type ccd;
   fits_file xfits,wfits,gfits;
   int axis;
   fitskey_type key;
   regre_block xreg;
   double mnorm,*xval,*wval;
   int row,ord;
   double start,step,*wavstart,*wavstep,*logstart,*logstep;
   double lamcen[1024],disp[1024];
   int midord;
   double log_step;
   int k;
   double ml[66],xc[66],ml2[66];
   double mla,mlb;

   if (argc != 4)
      get_error("Usage:  rebin_echelle  <in> <linout> <logout>");

   inform_title("Wavelength rebinning");

   read_image(&xfits,argv[1]);
   check_image_2d(&xfits);
   get_specinfo_fits(&xfits,&xfits.head,&spec);
   get_ccdinfo_fits(&xfits,&xfits.head,&ccd);

   get_dispersion_list(&spec,lamcen,disp);
   midord = (spec.firstord + spec.lastord) / 2;

   wfits = xfits;
   set_fits_name(&wfits,argv[2]);
   wfits.npix[1]= (xfits.npix[1] * 3) / 2;
   write_keyword_integer
     (&wfits,&wfits.head,"NAXIS1",wfits.npix[1],"Number of bins");

   wfits.crpix[1] = 1.0;
   wfits.crval[1] = 1.0;
   wfits.cdelt[1] = 1.0;

   for (axis=1;axis<=2;axis++)
   {
      sprintf(key,"CRPIX%d",axis);
      write_keyword_double(&wfits,&wfits.head,key,wfits.crpix[axis],
                                           "Reference pixel");
      sprintf(key,"CRVAL%d",axis);
      write_keyword_double(&wfits,&wfits.head,key,wfits.crval[axis],
                                           "Coordinate at reference pixel");
      sprintf(key,"CDELT%d",axis);
      write_keyword_double(&wfits,&wfits.head,key,wfits.cdelt[axis],
                                           "Coordinate increment per pixel");
   }
   wfits.totpix = wfits.npix[1]*wfits.npix[2];
   wfits.binstart = wfits.head.nrec * RECORD_BYTES;
   set_img_totbinrec(&wfits);

   gfits = wfits;
   set_fits_name(&gfits,argv[3]);

   allocate_whole_binary(&wfits);
   allocate_whole_binary(&gfits);

   read_regression_block(&xfits,&xfits.head,"XREGRE",&xreg);
   mb = xreg.deg[1] + 1;

   b=dvector(mb);
   xval=dvector(xfits.npix[1]);
   wval=dvector(wfits.npix[1]);
   wavstart=dvector(xfits.npix[2]);
   wavstep=dvector(xfits.npix[2]);
   logstart=dvector(xfits.npix[2]);
   logstep=dvector(xfits.npix[2]);

   log_step = roundx(disp[midord] * ccd.pixsize[1] / lamcen[midord],4);
   for (row=1;row<=xfits.npix[2];row++) logstep[row] = log_step;

   inform("%d rows to process.",xfits.npix[2]);

   fprintf(stderr,"Processing row number:   0");
   for (row=1;row<=xfits.npix[2];row++)
   {
      fprintf(stderr,"\b\b\b%3d",row);
      xord = world_coordinate(&xfits,2,(double)row);
      ord = nint(xord);
      mnorm = normalized_order(&spec,xord);
      accumulate_coefficients(xreg.a,xreg.ma,b,mb,mnorm,2);

      for (k=1;k<=65;k++)
      {
        ml[k] = (double)(k-1) / 64.0;
        xc[k] = compute_polynomial(b,mb,ml[k]);
      }

      spline(xc,ml,65,1e30,1e30,ml2);
      splint(xc,ml,ml2,65,1.0,&mla);
      splint(xc,ml,ml2,65,(double)xfits.npix[1],&mlb);

      start = denormalized_mlambda(&spec,(5.0*mla-mlb)*0.25) / xord;
      step = disp[ord] * ccd.pixsize[1];
      wavstart[row] = roundx(start,8);
      wavstep[row] = roundx(step,4);
      logstart[row] = roundx(log(start),8);

      collect_pixel_row(&xfits,row,&xval[1]);
      rebin_spectrum(xval,xfits.npix[1],1.0,1.0,
                wval,wfits.npix[1],wavstart[row],wavstep[row],func_lin);

      deposit_pixel_row(&wfits,row,&wval[1]);

      rebin_spectrum(xval,xfits.npix[1],1.0,1.0,
                wval,wfits.npix[1],logstart[row],logstep[row],func_log);

      deposit_pixel_row(&gfits,row,&wval[1]);
   }
   fprintf(stderr,"\n");

   write_double_array(&wfits,&wfits.head,"STARTVAL",
                      &wavstart[1],xfits.npix[2]);
   write_double_array(&wfits,&wfits.head,"STEPVAL",
                      &wavstep[1],xfits.npix[2]);

   write_double_array(&gfits,&gfits.head,"STARTVAL",
                      &logstart[1],xfits.npix[2]);
   write_double_array(&gfits,&gfits.head,"STEPVAL",
                      &logstep[1],xfits.npix[2]);

   inform("Linear and logarithmic rebinning done.");

   write_image(&wfits,argv[2]);
   write_image(&gfits,argv[3]);

   inform("Images '%s' and '%s' created.",wfits.file.path,gfits.file.path);

   free_binary_data(&xfits);

   free_dvector(b);
   free_dvector(xval);
   free_dvector(wval);
   free_dvector(wavstart);
   free_dvector(wavstep);
   free_dvector(logstart);
   free_dvector(logstep);

   check_memory();

   return(0);
}
