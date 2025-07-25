/* ------------------------------------------------------------------------

   Program:  splice_echelle
   Purpose:  Splice together the individual echelle orders into one long
             spectrum.

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

double fx(double u)
{
   return(u);
}

int main(int argc,char **argv)
{
   fits_file xfits,yfits;
   int numpix,numord,firstord,lastord;
   double *startval,*stepval;
   double *xval,*yval,*cval,*rval;
   file_name_type tmp;
   double step;
   int row,nbins,i,j,k;
   double t,rstart,cstart;
   int first,last,nmid;

   if (argc != 2) get_error("Usage:  splice_echelle  <img>");

   setbuf(stderr,NULL);

   inform_title("Splicing echelle orders");

   read_image(&xfits,argv[1]);
   check_image_2d(&xfits);
   check_real_image(&xfits);

   numpix = xfits.npix[1];
   if (numpix < 16) get_error
     ("Spectral orders too short in '%s'!",xfits.file.path);
   numord = xfits.npix[2];
   if (numord < 2) get_error
     ("At least two echelle orders expected in '%s'!",xfits.file.path);

   if (xfits.crpix[2] != 1.0)
     get_error("Invalid reference pixel CRPIX2 in '%s'!",xfits.file.path);
   if (xfits.cdelt[2] != -1.0)
     get_error("Decreasing order numbers expected in '%s'!",xfits.file.path);
   if (xfits.crval[2] != floor(xfits.crval[2]))
     get_error("Fractional order number detected in '%s'!",xfits.file.path);

   lastord = (int)xfits.crval[2];
   firstord = lastord - numord + 1;
   if (firstord < 50 || lastord > 200)
     get_error("Order number out of range in '%s'!",xfits.file.path);

   startval=dvector(numord);
   stepval=dvector(numord);

   get_double_array(&xfits,&xfits.head,"STARTVAL",&startval[1],numord);
   get_double_array(&xfits,&xfits.head,"STEPVAL",&stepval[1],numord);

   for (row=1;row<=numord;row++)
     if (startval[row] < 7.0 || startval[row] > 10.0)
       get_error("Unexpected start value detected in '%s'!",xfits.file.path);

   for (row=1;row<=numord;row++)
     if (stepval[row] < 0.0 || stepval[row] > 1.0e-3)
       get_error("Unexpected step value detected in '%s'!",xfits.file.path);

   for (row=1,step=stepval[1];row<=numord;row++)
     if (stepval[row] != step)
       get_error("Unequal step detected in '%s'!",xfits.file.path);

   nbins = (int)ceil((startval[numord]-startval[1])/step) + 2*numpix;

   create_image(next_tmp_file_name(tmp),&xfits.head,xfits.bitpix,
                1,nbins,1.0,startval[1],step);

   read_image(&yfits,tmp);

   rval = dvector(numpix);
   xval = dvector(numpix);
   yval = dvector(nbins);

   clear_dvector(yval,nbins);

   inform("%d rows to process (%d..%d).",numord,firstord,lastord);

   fprintf(stderr,"Processing row number:   0");
   for (row=1;row<=numord;row++)
   {
      fprintf(stderr,"\b\b\b%3d",row);
      collect_pixel_row(&xfits,row,&xval[1]);
      for (first=1;first<=numpix;first++) if (xval[first] != 0.0) break;
      if (first > numpix) continue;
      for (last=numpix;last>=1;last--) if (xval[last] != 0.0) break;
      if (last < 1) continue;
      nmid = last - first + 1;
      if (nmid < 2) continue;
      cval = &xval[first-1];
      cstart = startval[row] + (first - 1) * step;
      t = (cstart - startval[1]) / step;
      k = (int)ceil(t);
      rstart = startval[1] + (k - 1) * step;
      linear_rebinning(cval,nmid,cstart,step,
        rval,numpix,rstart,step,fx);
      for (i=1,j=k;i<=numpix;i++,j++)
      {
        if (rval[i] == 0.0 || yval[j] == 0.0)
          yval[j] += rval[i];
        else
          yval[j] = (yval[j] + rval[i]) * 0.5;
      }
   }
   fprintf(stderr,"\n");

   deposit_pixel_row(&yfits,1,&yval[1]);
   update_fits_minmax(&yfits);
   write_keyword_textual(&yfits,&yfits.head,"CTYPE1","LOG_LAM","Units");
   write_image(&yfits,"%s-long",argv[1]);
   inform("Image '%s' created",yfits.file.path);

   free_dvector(rval);
   free_dvector(xval);
   free_dvector(yval);

   free_dvector(startval);
   free_dvector(stepval);

   free_binary_data(&xfits);

   check_memory();

   return(0);
}
