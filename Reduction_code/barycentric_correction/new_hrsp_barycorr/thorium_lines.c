/* ------------------------------------------------------------------------

   Program:  thorium_lines
   Purpose:  Locate the emission lines in a given thorium image.

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
   spec_type spec;
   fits_file imgfits,linfits,transfits;
   trans_table transtbl;
   lin_table lintbl;
   file_name_type tmp;
   double **p;
   int seq,pold,pnew,row,ndone;
   gaussian g;
   float *pix;
   regre_block xreg;
   char midtime[MAX_TIME+1];
   int fib;
   ccd_type ccd;
   double fwhm;

   if (argc != 3) get_error("Usage:  thorium_lines  <img> <tbl>");

   setbuf(stderr,NULL);

   inform_title("Thorium line detection");

   read_image(&imgfits,argv[1]);
   get_ccdinfo_fits(&imgfits,&imgfits.head,&ccd);
   get_specinfo_fits(&imgfits,&imgfits.head,&spec);

   read_trans_table(&transfits,&transtbl,"transform");

   create_lin_table(transfits.nrows,next_tmp_file_name(tmp));
   read_lin_table(&linfits,&lintbl,tmp);
   remove_tmp_file("%s.fit",tmp);

   get_keyword_textual(&imgfits,&imgfits.head,"MIDTIME",midtime,MAX_TIME);

   fib = get_keyword_integer(&imgfits,&imgfits.head,"FIBRENO");

   copy_integer_array(transtbl.order,transfits.nrows,lintbl.order);
   copy_double_array(transtbl.airlam,transfits.nrows,lintbl.airlam);
   copy_double_array(transtbl.vaclam,transfits.nrows,lintbl.vaclam);
   copy_double_array(transtbl.waveno,transfits.nrows,lintbl.waveno);
   copy_string_array(transtbl.species,transfits.nrows,8,lintbl.species);
   copy_double_array(transtbl.ustart,transfits.nrows,lintbl.ustart);
   copy_double_array(transtbl.ucen,transfits.nrows,lintbl.ucen);
   copy_double_array(transtbl.uend,transfits.nrows,lintbl.uend);
   copy_double_array(transtbl.vstart,transfits.nrows,lintbl.vstart);
   copy_double_array(transtbl.vcen,transfits.nrows,lintbl.vcen);
   copy_double_array(transtbl.vend,transfits.nrows,lintbl.vend);

   read_regression_block(&transfits,&transfits.extend,"XREGRE",&xreg);

   p = dmatrix(1,2);

   for (seq=0;seq<linfits.nrows;seq++)
   {
      p[1][2] = normalized_v(&spec,lintbl.vcen[seq]);
      p[1][1] = normalized_u(&spec,lintbl.ustart[seq]);
      lintbl.xstart[seq] = polynomial_val(p,1,2,xreg.deg,xreg.a,xreg.ma);
      p[1][1] = normalized_u(&spec,lintbl.ucen[seq]);
      lintbl.xcenter[seq] = polynomial_val(p,1,2,xreg.deg,xreg.a,xreg.ma);
      p[1][1] = normalized_u(&spec,lintbl.uend[seq]);
      lintbl.xend[seq] = polynomial_val(p,1,2,xreg.deg,xreg.a,xreg.ma);
      lintbl.status[seq] = -1;
   }

   fwhm = spec.fib[fib].dx / ccd.pixsize[1];

   fprintf(stderr,"Locating spectral lines:   0%%");
   for (seq=ndone=0,pold=-1;seq<linfits.nrows;seq++)
   {
      pnew = nint(seq*100.0/(double)(linfits.nrows-1));
      if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
      row = (int)(world_coordinate(&imgfits,2,1.0)-lintbl.order[seq]) + 1;
      if (row < 1 || row > imgfits.npix[2]) continue;
      pix = imgfits.pix + (row-1)*imgfits.npix[1];
      locate_gaussian_1d(pix,imgfits.npix[1],1.0,1.0,lintbl.xstart[seq],
                         lintbl.xcenter[seq],lintbl.xend[seq],fwhm,&g);
      lintbl.base[seq] = g.base;
      lintbl.icen[seq] = g.icen;
      lintbl.xcen[seq] = g.xcen;
      lintbl.xfwhm[seq] = g.xfwhm;
      lintbl.chisq[seq] = g.chisq;
      lintbl.rms[seq] = g.rms;
      lintbl.niter[seq] = g.niter;
      lintbl.status[seq] = g.status;
      if (g.status == 0) ndone++;
   }
   fprintf(stderr,"\n");

   inform("%d lines detected.",ndone);

   write_keyword_textual(&linfits,&linfits.extend,"SPECNAME",spec.specname,
                         "Spectrograph name");
   write_keyword_textual(&linfits,&linfits.extend,"MIDTIME",midtime,
                         "Universal time at mid-exposure");

   write_table(&linfits,argv[2]);
   inform("Table '%s' created.",linfits.file.path);

   free_dmatrix(p);

   free_binary_data(&imgfits);
   free_binary_data(&transfits);

   check_memory();

   return(0);
}
