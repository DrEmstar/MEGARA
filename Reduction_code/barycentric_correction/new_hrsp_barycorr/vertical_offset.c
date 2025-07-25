/* ------------------------------------------------------------------------

   Program:  vertical_offset
   Purpose:  Determine the vertical offset for a given stellar image

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
   fits_file imgfits,ordfits,rotfits;
   file_name_type tmp;
   param_type par;
   double *voffset,*column,*fit,*del;
   int *sel;
   regre_block yreg,mreg;
   double **p,ystart,ycenter,yend,centval,alpha,rms;
   int midord,pold,pnew,point,x,nsel;
   gaussian g;
   float *pix;
   double a[4];
   ccd_type ccd;
   double fwhm;

   if (argc != 2) get_error("Usage:  vertical_offset  <img>");

   setbuf(stderr,NULL);

   inform_title("Vertical offset determination");

   read_image_header(&imgfits,argv[1]);
   read_table_header(&ordfits,"order");

   get_ccdinfo_fits(&imgfits,&imgfits.head,&ccd);
   get_specinfo_fits(&imgfits,&imgfits.head,&spec);
   get_parameters(&par,spec.specname);

   rotate_image_proc(argv[1],next_tmp_file_name(tmp),-90,YES_INFO);
   read_image(&rotfits,tmp);
   remove_tmp_file("%s.fit",tmp);

   read_regression_block(&ordfits,&ordfits.extend,"YREGRE",&yreg);
   read_regression_block(&ordfits,&ordfits.extend,"MREGRE",&mreg);

   p = dmatrix(1,2);
   p[1][2] = p[1][1] = 0.5;
   midord = nint(polynomial_val(p,1,2,mreg.deg,mreg.a,mreg.ma));

   voffset = dvector(par.echdef_ncols);
   column = dvector(par.echdef_ncols);
   fit = dvector(par.echdef_ncols);
   del = dvector(par.echdef_ncols);
   sel = ivector(par.echdef_ncols);

   alpha = (double)(imgfits.npix[1]-1)/(double)(par.echdef_ncols-1);

   fwhm = spec.fib[get_fibre_number(&imgfits)].dy / ccd.pixsize[2];

   fprintf(stderr,"Tracing order No. %d:   0%%",midord);
   for (point=1,pold=-1;point<=par.echdef_ncols;point++)
   {
      pnew = nint((double)(point-1)*100.0/(double)(par.echdef_ncols-1));
      if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
      x = nint((double)(point-1)*alpha)+1;
      pix = rotfits.pix + (imgfits.npix[1]-x)*imgfits.npix[2];
      p[1][1] = normalized_pixel(&imgfits,1,(double)x);
      p[1][2] = normalized_order(&spec,(double)midord+0.5);
      ystart = polynomial_val(p,1,2,yreg.deg,yreg.a,yreg.ma);
      p[1][2] = normalized_order(&spec,(double)midord);
      ycenter = polynomial_val(p,1,2,yreg.deg,yreg.a,yreg.ma);
      p[1][2] = normalized_order(&spec,(double)midord-0.5);
      yend = polynomial_val(p,1,2,yreg.deg,yreg.a,yreg.ma);
      locate_gaussian_1d(pix,imgfits.npix[2],1.0,1.0,
                         ystart,ycenter,yend,fwhm,&g);
      voffset[point] = g.xcen - ycenter;
      column[point] = (double)x;
      sel[point] = 1;
   }
   fprintf(stderr,"\n");

   best_regression_polynomial_1d(column,voffset,par.echdef_ncols,sel,2,
                                 a,&rms,&nsel,fit,del,3.0,NO_INFO);

   centval = fit[(1+par.echdef_ncols)/2];
   inform("Vertical offset: %8.3f",centval);
   inform("Total number of data points used: %d. Rejected: %d.",
           par.echdef_ncols,par.echdef_ncols-nsel);
   inform("Parabolic fit r.m.s. error: %8.3f",rms);

   free_dvector(voffset);
   free_dvector(column);
   free_dvector(fit);
   free_dvector(del);
   free_ivector(sel);

   free_dmatrix(p);

   free_binary_data(&rotfits);

   write_image_descriptor_double(argv[1],"VOFFSET",centval,"Vertical offset");

   check_memory();

   return(0);
}
