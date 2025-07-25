/* ------------------------------------------------------------------------

   Program:  background_echelle
   Purpose:  Estimate the background scattered light for a given echelle
             spectrum

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
param_type par;
fits_file xfits,ordfits;
regre_block mreg,yreg;
int firstord,lastord;
double voffset;

void order_limits(void)
{
   double **p;
   double minord,maxord,xord;
   int i;

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
   lastord = (int)floor(maxord) + 1;

   free_dmatrix(p);
}

void method_cubic(void)
{
  int col,row,m,n;
  double order;
  double **p,*x,*y;
  int k,ka,kb,pold,pnew,seq;
  double bkg;

  inform("Using third-order polynomial interpolation.");

  p = dmatrix(1,2);
  x = dvector(xfits.npix[2]);
  y = dvector(xfits.npix[2]);

  fprintf(stderr,"Background image creation:   0%%");
  for (col=1,pold=-1,seq=1;col<=xfits.npix[1];col++)
  {
    p[1][1] = normalized_pixel(&xfits,1,(double)col);
    for (m=lastord,n=0;m>=firstord;m--)
    {
      order = (double)m-0.5;
      p[1][2] = normalized_order(&spec,order);
      row = nint(polynomial_val(p,1,2,yreg.deg,yreg.a,yreg.ma)+voffset);
      if (row < 1 || row > xfits.npix[2]) continue;
      n++;
      x[n] = (double)row;
      y[n] = (double)collect_pixel_value(&xfits,col,row);
    }
    if (n < 4) get_error
      ("Column %d: Not enough data points for interpolation!",col);
    for (row=1;row<=xfits.npix[2];row++,seq++)
    {
      k = locate_node(x,n,(double)row);
      if (k < n/2)
      {
        ka = k - 1;
        if (ka < 1) ka = 1;
      }
      else
      {
        kb = k + 2;
        if (kb > n) kb = n;
        ka = kb - 3;
      }
      bkg = interpolation(&x[ka],&y[ka],4,(double)row);
      deposit_pixel_value(&xfits,col,row,bkg);
      pnew = nint(seq*100.0/(double)xfits.totpix);
      if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
    }
  }
  fprintf(stderr,"\n");

  free_dmatrix(p);
  free_dvector(x);
  free_dvector(y);
}

void method_spline(void)
{
  int col,row,m,n;
  double order;
  double **p,*x,*y,*y2;
  int pold,pnew,seq;
  double bkg;

  inform("Using spline interpolation.");
  p = dmatrix(1,2);
  x = dvector(xfits.npix[2]);
  y = dvector(xfits.npix[2]);
  y2 = dvector(xfits.npix[2]);

  fprintf(stderr,"Background image creation:   0%%");
  for (col=1,pold=-1,seq=1;col<=xfits.npix[1];col++)
  {
    p[1][1] = normalized_pixel(&xfits,1,(double)col);
    for (m=lastord,n=0;m>=firstord;m--)
    {
      order = (double)m-0.5;
      p[1][2] = normalized_order(&spec,order);
      row = nint(polynomial_val(p,1,2,yreg.deg,yreg.a,yreg.ma)+voffset);
      if (row < 1 || row > xfits.npix[2]) continue;
      n++;
      x[n] = (double)row;
      y[n] = (double)collect_pixel_value(&xfits,col,row);
    }
    if (n < 4) get_error
      ("Column %d: Not enough data points for interpolation!",col);
    spline(x,y,n,1e30,1e30,y2);
    for (row=1;row<=xfits.npix[2];row++,seq++)
    {
      splint(x,y,y2,n,(double)row,&bkg);
      deposit_pixel_value(&xfits,col,row,bkg);
      pnew = nint(seq*100.0/(double)xfits.totpix);
      if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
    }
  }
  fprintf(stderr,"\n");

  free_dmatrix(p);
  free_dvector(x);
  free_dvector(y);
  free_dvector(y2);
}

void method_polyfit(void)
{
   fits_file backfits;
   back_table backtbl;
   file_name_type tmp;
   int nord,nrows;
   double alpha,order,xnorm,ynorm;
   int m,seq,pold,pnew,point,x,y,mb;
   double **p,*b;
   float *pix;
   regre_block reg;

   inform("Using two-dimensional polynomial fit.");

   nord = lastord - firstord + 1;
   nrows = nord*par.echback_ncols;
   create_back_table(nrows,next_tmp_file_name(tmp));
   read_back_table(&backfits,&backtbl,tmp);
   remove_tmp_file("%s.fit",tmp);

   alpha = (double)(xfits.npix[1]-1)/(double)(par.echback_ncols-1);

   p = dmatrix(1,2);

   fprintf(stderr,"Collecting pixel values:   0%%");
   for (m=firstord,seq=0,pold=-1;m<=lastord;m++)
   {
      order = (double)m-0.5;
      p[1][2] = normalized_order(&spec,order);
      for (point=1;point<=par.echback_ncols;point++,seq++)
      {
         pnew = nint(seq*100.0/(double)(nrows-1));
         if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
         x = nint((double)(point-1)*alpha)+1;
         p[1][1] = normalized_pixel(&xfits,1,(double)x);
         y = nint(polynomial_val(p,1,2,yreg.deg,yreg.a,yreg.ma)+voffset);
         backtbl.x[seq] = x;
         backtbl.y[seq] = y;
         backtbl.order[seq] = order;
         if (pixel_inside(&xfits,x,y))
         {
            backtbl.sel[seq] = 1;
            backtbl.pixval[seq] = (double)collect_pixel_value(&xfits,x,y);
         }
         else
            backtbl.sel[seq] = 0;
      }
   }
   fprintf(stderr,"\n");

   for (seq=0;seq<nrows;seq++)
   {
      backtbl.xnorm[seq] = normalized_pixel(&xfits,1,(double)backtbl.x[seq]);
      backtbl.ynorm[seq] = normalized_pixel(&xfits,2,(double)backtbl.y[seq]);
   }

   write_table(&backfits,"back");
   inform("Table '%s' created.",backfits.file.path);

   inform("");
   inform("Two-dimensional fit: Pixval = f(xnorm,ynorm)");

   regression_polynomial_proc("back","Pixval","Xnorm,Ynorm",
                    par.echback.degarg,
                    "Sel","Fit","Resid","REGRE",par.echback.kappa);

   read_table_header(&backfits,"back");
   read_regression_block(&backfits,&backfits.extend,"REGRE",&reg);

   mb = reg.deg[1]+1;
   b = dvector(mb);

   inform("");
   fprintf(stderr,"Background image creation:   0%%");
   for (y=1,pix=xfits.pix,seq=0,pold=-1;y<=xfits.npix[2];y++)
   {
      ynorm = normalized_pixel(&xfits,2,(double)y);
      accumulate_coefficients(reg.a,reg.ma,b,mb,ynorm,2);
      for (x=1;x<=xfits.npix[1];x++,pix++,seq++)
      {
         pnew = nint(seq*100.0/(double)(xfits.totpix-1));
         if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
         xnorm = normalized_pixel(&xfits,1,(double)x);
         *pix = (float)compute_polynomial(b,mb,xnorm);
      }
   }
   fprintf(stderr,"\n");

   free_dvector(b);
   free_dmatrix(p);
}

int main(int argc,char **argv)
{
   if (argc != 3) get_error("Usage:  background_echelle  <in> <out>");

   setbuf(stderr,NULL);

   inform_title("Background estimation");

   read_image(&xfits,argv[1]);
   get_specinfo_fits(&xfits,&xfits.head,&spec);
   read_table_header(&ordfits,"order");

   set_fits_name(&xfits,argv[2]);

   get_parameters(&par,spec.specname);

   voffset = get_keyword_double(&xfits,&xfits.head,"VOFFSET");

   read_regression_block(&ordfits,&ordfits.extend,"MREGRE",&mreg);
   read_regression_block(&ordfits,&ordfits.extend,"YREGRE",&yreg);

   order_limits();

   switch (par.back_method.id)
   {
     case BACK_CUBIC: method_cubic();
                      break;
     case BACK_SPLINE: method_spline();
                      break;
     default: method_polyfit();
   }

   write_image(&xfits,argv[2]);
   inform("Image '%s' created.",xfits.file.path);

   check_memory();

   return(0);
}
