/* ------------------------------------------------------------------------

   Program:  form_eflat
   Purpose:  Create an extracted (one-dimensional) flat-field image from
             an extracted white-lamp image.

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
   fits_file imgfits,flagfits;
   char spectrograph[MAX_SPECNAME+1];
   param_type par;
   float *imgptr,pixmax;
   int row,col,left,right,first,last,n,nsel,i,pold,pnew,seq;
   byte *flagptr;
   double *pixval,*pixno,*xnorm,*fit,*resid,rms,*c;
   int *sel;

   if (argc != 2) get_error("Usage:  form_eflat  <img>");

   setbuf(stderr,NULL);

   inform_title("Extracted flat-field image creation");

   read_image(&imgfits,argv[1]);
   read_image(&flagfits,"%s#",argv[1]);
   get_spectrograph(&imgfits,&imgfits.head,spectrograph);
   get_parameters(&par,spectrograph);

   inform("%d orders to process.",imgfits.npix[2]);

   pixval = dvector(imgfits.npix[1]);
   pixno = dvector(imgfits.npix[1]);
   xnorm = dvector(2*par.flat_rfilt+1);
   fit = dvector(2*par.flat_rfilt+1);
   resid = dvector(2*par.flat_rfilt+1);
   sel = ivector(2*par.flat_rfilt+1);
   c = dvector(par.flat.deg[1]+1);

   for (col=1;col<=imgfits.npix[1];col++) pixno[col] = (double)col;

   if (par.flat_smooth == 1)
   {
     imgptr=imgfits.pix - 1;
     flagptr=flagfits.bin - 1;
     fprintf(stderr,"Image filter:   0%%");
     for (row=1,seq=0;row<=imgfits.npix[2];row++)
     {
        collect_pixel_row(&imgfits,row,&pixval[1]);
        for (col=1;col<=imgfits.npix[1];col++,seq++)
        {
           pnew = nint(seq*100.0/(double)(imgfits.totpix-1));
           if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
           if (flagptr[col] >= 0x80) continue;

           for (left=col-1;left>=1;left--)
           if (flagptr[left] >= 0x80) break;

           for (right=col+1;right<=imgfits.npix[1];right++)
           if (flagptr[right] >= 0x80) break;

           if (2*col < left+right)
           {
              first = imax(col-par.flat_rfilt,left+1);
              last = imin(first+2*par.flat_rfilt,right-1);
           }
           else
           {
              last = imin(col+par.flat_rfilt,right-1);
              first = imax(last-2*par.flat_rfilt,left+1);
           }

           n = last - first + 1;

           if (n < 2*par.flat_rfilt)
           {
              fit[col-first+1] = pixval[col];
           }
           else
           {
              clear_ivector(sel,n);
              for (i=1;i<=n;i++) sel[i] = 1;
              for (i=1;i<=n;i++) xnorm[i] = (i-1)/(double)(n-1);
              best_regression_polynomial_1d(xnorm,&pixval[first-1],
                            n,sel,par.flat.deg[1],c,
                            &rms,&nsel,fit,resid,par.flat.kappa,NO_INFO);
           }

           imgptr[col] = (float)fit[col-first+1];
        }
        imgptr += imgfits.npix[1];
        flagptr += flagfits.npix[1];
     }
     fprintf(stderr,"\n");
   }

   for (i=0,imgptr=imgfits.pix,pixmax=0;i<imgfits.totpix;i++,imgptr++)
   {
     if (*imgptr < par.flat_thres_ext) *imgptr = 0.0;
     if (*imgptr > pixmax) pixmax = *imgptr;
   }

   if (pixmax <= 0.0) get_error("A positive maximum pixel value expected!");

   for (i=0,imgptr=imgfits.pix;i<imgfits.totpix;i++,imgptr++)
     *imgptr /= pixmax;

   write_keyword_textual(&imgfits,&imgfits.head,"WHITFILE",
                         imgfits.file.name,
                         "Image used for flat field creation");

   update_fits_minmax(&imgfits);
   write_image(&imgfits,"eflat");

   free_dvector(pixval);
   free_dvector(pixno);
   free_dvector(xnorm);
   free_dvector(fit);
   free_dvector(resid);
   free_ivector(sel);
   free_dvector(c);

   free_binary_data(&flagfits);

   inform("Image '%s' created.",imgfits.file.path);

   check_memory();

   return(0);
}
