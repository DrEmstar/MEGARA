/* ------------------------------------------------------------------------

   Program:  extract_echelle
   Purpose:  Extract the one-dimensional spectrum from an echelle image

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

spec_type spec;
ccd_type ccd;

int main(int argc,char **argv)
{
   fits_file xfits,yfits,zfits,tblfits;
   param_type par;
   file_name_type tmp;
   double voffset,bias,maxval,minord,maxord,xord;
   double **p;
   int i,firstord,lastord,midord,nord,ord,ka,kb,kc,k,pnew,pold,seq;
   float *pix,val;
   byte *flag;
   double ycen,yfwhm;
   regre_block mreg,yreg,wreg;
   double *startval,*stepval;
   int rslit,slitsize;

   if (argc != 3) get_error("Usage:  extract_echelle  <in> <out>");

   setbuf(stderr,NULL);

   inform_title("Echelle order extraction");

   read_image(&xfits,argv[1]);
   get_specinfo_fits(&xfits,&xfits.head,&spec);
   get_ccdinfo_fits(&xfits,&xfits.head,&ccd);
   get_parameters(&par,spec.specname);
   read_table_header(&tblfits,"order");

   voffset = get_keyword_double(&xfits,&xfits.head,"VOFFSET");
   bias = get_keyword_double(&xfits,&xfits.head,"ADU_BIAS");
   check_keyword_logical(&xfits,&xfits.head,"BIASDONE",'T');

   read_regression_block(&tblfits,&tblfits.extend,"MREGRE",&mreg);

   p = dmatrix(1,2);

   p[1][2] = normalized_pixel(&xfits,2,1.0+voffset);
   for (i=1,maxord=0.0;i<=xfits.npix[1];i++)
   {
      p[1][1] = normalized_pixel(&xfits,1,(double)i);
      xord = polynomial_val(p,1,2,mreg.deg,mreg.a,mreg.ma);
      if (xord > maxord) maxord = xord;
   }

   p[1][2] = normalized_pixel(&xfits,2,(double)xfits.npix[2]+voffset);
   for (i=1,minord=maxord;i<=xfits.npix[1];i++)
   {
      p[1][1] = normalized_pixel(&xfits,1,(double)i);
      xord = polynomial_val(p,1,2,mreg.deg,mreg.a,mreg.ma);
      if (xord < minord) minord = xord;
   }

   firstord = (int)floor(minord) + 1;
   lastord = (int)floor(maxord);
   nord = lastord - firstord + 1;
   midord = nint((double)(firstord+lastord)*0.5);

   inform("%d orders to extract (%d-%d).",nord,firstord,lastord);

   yfits = xfits;

   set_fits_name(&yfits,argv[2]);

   yfits.npix[2] = nord;
   yfits.totpix = yfits.npix[1]*yfits.npix[2];
   yfits.crval[2] = (double)lastord;
   yfits.cdelt[2] = -1.0;
   set_img_totbinrec(&yfits);
   allocate_whole_binary(&yfits);

   write_keyword_integer(&yfits,&yfits.head,"NAXIS2",yfits.npix[2],
                                         "Number of rows");
   write_keyword_double(&yfits,&yfits.head,"CRVAL2",yfits.crval[2],
                                         "Coordinate at reference pixel");
   write_keyword_double(&yfits,&yfits.head,"CDELT2",yfits.cdelt[2],
                                         "Coordinate increment per pixel");

   read_regression_block(&tblfits,&tblfits.extend,"YREGRE",&yreg);
   read_regression_block(&tblfits,&tblfits.extend,"WREGRE",&wreg);

   create_image(next_tmp_file_name(tmp),NULL,8,2,yfits.npix[1],yfits.npix[2],
                yfits.crpix[1],yfits.crval[1],yfits.cdelt[1],
                yfits.crpix[2],yfits.crval[2],yfits.cdelt[2]);
   read_image(&zfits,tmp);
   remove_tmp_file("%s.fit",tmp);

   maxval = ccd.maxpix - bias;

   p[1][1] = 0.5;
   p[1][2] = normalized_order(&spec,(double)midord);
   yfwhm = polynomial_val(p,1,2,wreg.deg,wreg.a,wreg.ma);
   rslit = nint(yfwhm*par.extract_slit/2) + 1;
   slitsize = 2 * rslit + 1;
   inform("Order profile FWHM: %4.2f pix",yfwhm);
   inform("Extraction slit height: %d pix",slitsize);

   fprintf(stderr,"Order extraction:   0%%");
   pix=yfits.pix,flag=zfits.bin;
   for (ord=lastord,seq=0,pold=-1;ord>=firstord;ord--)
   {
      p[1][2] = normalized_order(&spec,(double)ord);
      for (i=1;i<=yfits.npix[1];i++,pix++,flag++,seq++)
      {
         pnew = nint(seq*100.0/(double)(yfits.totpix-1));
         if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
         p[1][1] = normalized_pixel(&xfits,1,(double)i);
         ycen = polynomial_val(p,1,2,yreg.deg,yreg.a,yreg.ma)
                                                                 + voffset;
         kc = nint(ycen);
         ka = kc - rslit;
         kb = kc + rslit;

         if (ka >= 1 && kb <= xfits.npix[2])
         {
            for (k=ka,*pix=0.0,*flag='\x00';k<=kb;k++)
            {
               val = collect_pixel_value(&xfits,i,k);
               if (val >= maxval) *flag |= '\x01';
               *pix += val;
            }
         }
         else
         {
            *flag = '\x80';
            *pix = 0.0;
         }
      }
   }
   fprintf(stderr,"\n");

   startval = dvector(nord);
   stepval = dvector(nord);

   for (ord=1;ord<=nord;ord++) startval[ord]=stepval[ord]=1.0;

   write_double_array(&yfits,&yfits.head,"STARTVAL",&startval[1],nord);
   write_double_array(&yfits,&yfits.head,"STEPVAL",&stepval[1],nord);

   write_double_array(&zfits,&zfits.head,"STARTVAL",&startval[1],nord);
   write_double_array(&zfits,&zfits.head,"STEPVAL",&stepval[1],nord);

   free_dvector(startval);
   free_dvector(stepval);

   copy_fits_keyword(&tblfits.extend,"TRANSIMG",&yfits,&yfits.head);
   copy_fits_keyword(&tblfits.extend,"ORDERIMG",&yfits,&yfits.head);
   copy_fits_keyword(&tblfits.extend,"WHITEJD",&yfits,&yfits.head);
   copy_fits_keyword(&tblfits.extend,"WHITEIMG",&yfits,&yfits.head);

   write_keyword_integer(&yfits,&yfits.head,"SLITSIZE",slitsize,
                                         "Extraction slit size (in pixels)");

   write_image(&yfits,argv[2]);
   inform("Output image '%s' created.",yfits.file.path);
   write_image(&zfits,"%s#",argv[2]);
   inform("Extraction status image '%s' created.",zfits.file.path);

   free_dmatrix(p);

   free_binary_data(&xfits);

   check_memory();

   return(0);
}
