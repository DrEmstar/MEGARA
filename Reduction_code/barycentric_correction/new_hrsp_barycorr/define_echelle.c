/* ------------------------------------------------------------------------

   Program:  define_echelle
   Purpose:  Locate the echelle orders in a CCD image.

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
   param_type par;
   ccd_type ccd;
   fits_file tharfits,transfits,imgfits,ordfits,rotfits;
   ord_table ordtbl;
   file_name_type tmp;
   int fibno;
   double **p,**q;
   int i,nrows;
   double minord,maxord,u,v,xord,alpha;
   int firstord,lastord,nord,pold,pnew,point,x,order,seq;
   double *mval,*pval,*pval2;
   gaussian g;
   float *pix;
   regre_block ureg,vreg,mreg,yreg,wreg;
   double fwhm;
   int ntot;

   if (argc != 2) get_error("Usage:  define_echelle  <img>");

   setbuf(stderr,NULL);

   inform_title("Echelle order definition");

   read_image_header(&imgfits,argv[1]);

   get_ccdinfo_fits(&imgfits,&imgfits.head,&ccd);
   get_specinfo_fits(&imgfits,&imgfits.head,&spec);
   get_parameters(&par,spec.specname);
   fibno = get_fibre_number(&imgfits);

   read_table_header(&tharfits,"%s/spec/%s/tbl/thar",hrsp_dir(),spec.specname);
   read_table_header(&transfits,"transform");

   read_regression_block(&transfits,&transfits.extend,"UREGRE",&ureg);
   read_regression_block(&transfits,&transfits.extend,"VREGRE",&vreg);
   read_regression_block(&tharfits,&tharfits.extend,"MREGRE",&mreg);

   p = dmatrix(1,2);
   q = dmatrix(1,2);

   p[1][2] = 0.0;
   for (i=1,maxord=0.0;i<=imgfits.npix[1];i++)
   {
      p[1][1] = normalized_pixel(&imgfits,1,(double)i);
      u = polynomial_val(p,1,2,ureg.deg,ureg.a,ureg.ma);
      v = polynomial_val(p,1,2,vreg.deg,vreg.a,vreg.ma);
      q[1][2] = normalized_v(&spec,v);
      q[1][1] = normalized_u(&spec,u);
      xord = polynomial_val(q,1,2,mreg.deg,mreg.a,mreg.ma);
      if (xord > maxord) maxord = xord;
   }

   p[1][2] = 1.0;
   for (i=1,minord=maxord;i<=imgfits.npix[1];i++)
   {
      p[1][1] = normalized_pixel(&imgfits,1,(double)i);
      u = polynomial_val(p,1,2,ureg.deg,ureg.a,ureg.ma);
      v = polynomial_val(p,1,2,vreg.deg,vreg.a,vreg.ma);
      q[1][2] = normalized_v(&spec,v);
      q[1][1] = normalized_u(&spec,u);
      xord = polynomial_val(q,1,2,mreg.deg,mreg.a,mreg.ma);
      if (xord < minord) minord = xord;
   }

   firstord = (int)floor(minord) + 1;
   lastord = (int)floor(maxord);
   nord = lastord - firstord + 1;

   inform("%d orders to trace (%d-%d).",nord,firstord,lastord);

   nrows = nord*par.echdef_ncols;
   create_ord_table(nrows,next_tmp_file_name(tmp));
   read_ord_table(&ordfits,&ordtbl,tmp);
   remove_tmp_file("%s.fit",tmp);

   write_keyword_textual(&ordfits,&ordfits.extend,"ORDERIMG",
                         imgfits.file.name,
                         "Image used for order tracing");
   copy_fits_keyword(&transfits.extend,"TRANSIMG",&ordfits,&ordfits.extend);
   copy_fits_keyword(&imgfits.head,"WHITEJD",&ordfits,&ordfits.extend);
   copy_fits_keyword(&imgfits.head,"WHITEIMG",&ordfits,&ordfits.extend);

   rotate_image_proc(argv[1],next_tmp_file_name(tmp),-90,YES_INFO);
   read_image(&rotfits,tmp);
   remove_tmp_file("%s.fit",tmp);

   alpha = (double)(imgfits.npix[1]-1)/(double)(par.echdef_ncols-1);

   mval = dvector(imgfits.npix[2]);
   pval = dvector(imgfits.npix[2]);
   pval2 = dvector(imgfits.npix[2]);

   for (i=1;i<=imgfits.npix[2];i++) pval[i] = (double)(imgfits.npix[2]+1-i);

   fwhm = spec.fib[fibno].dy / ccd.pixsize[2];

   fprintf(stderr,"Echelle order detection:   0%%");
   ntot = 0;
   for (point=1,seq=0,pold=-1;point<=par.echdef_ncols;point++)
   {
      x = nint((double)(point-1)*alpha)+1;
      p[1][1] = normalized_pixel(&imgfits,1,(double)x);
      for (i=1;i<=imgfits.npix[2];i++)
      {
         p[1][2] = normalized_pixel(&imgfits,2,pval[i]);
         u = polynomial_val(p,1,2,ureg.deg,ureg.a,ureg.ma);
         v = polynomial_val(p,1,2,vreg.deg,vreg.a,vreg.ma);
         q[1][2] = normalized_v(&spec,v);
         q[1][1] = normalized_u(&spec,u);
         mval[i] = polynomial_val(q,1,2,mreg.deg,mreg.a,mreg.ma);
      }
      spline(mval,pval,imgfits.npix[2],1e30,1e30,pval2);
      pix = rotfits.pix + (imgfits.npix[1]-x)*imgfits.npix[2];
      for (order=lastord;order>=firstord;order--,seq++)
      {
         pnew = nint(seq*100.0/(double)(nrows-1));
         if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
         ordtbl.order[seq] = order;
         ordtbl.x[seq] = x;

         splint(mval,pval,pval2,imgfits.npix[2],(double)order+0.5,
                                                 &ordtbl.ystart[seq]);
         splint(mval,pval,pval2,imgfits.npix[2],(double)order,
                                                 &ordtbl.ycenter[seq]);
         splint(mval,pval,pval2,imgfits.npix[2],(double)order-0.5,
                                                 &ordtbl.yend[seq]);

         if (nint(ordtbl.ystart[seq]) < 1 && 
             nint(ordtbl.ycenter[seq]) > 1)
                ordtbl.ystart[seq] = 1.0;
         if (nint(ordtbl.yend[seq]) > imgfits.npix[2] && 
             nint(ordtbl.ycenter[seq]) < imgfits.npix[2])
                ordtbl.yend[seq] =  (double)imgfits.npix[2];

         locate_gaussian_1d(pix,imgfits.npix[2],1.0,1.0,ordtbl.ystart[seq],
            ordtbl.ycenter[seq],ordtbl.yend[seq],fwhm,&g);

         ordtbl.base[seq] = g.base;
         ordtbl.icen[seq] = g.icen;
         ordtbl.ycen[seq] = g.xcen;
         ordtbl.yfwhm[seq] = g.xfwhm;
         ordtbl.chisq[seq] = g.chisq;
         ordtbl.rms[seq] = g.rms;
         ordtbl.niter[seq] = g.niter;
         ordtbl.status[seq] = g.status;

         if (g.status == 0) ntot++;
      }
   }
   fprintf(stderr,"\n");

   for (seq=0;seq<nrows;seq++)
   {
      ordtbl.xnorm[seq] = normalized_pixel(&imgfits,1,(double)ordtbl.x[seq]);
      ordtbl.ynorm[seq] = normalized_pixel(&imgfits,2,ordtbl.ycen[seq]);
      ordtbl.mnorm[seq] = normalized_order(&spec,(double)ordtbl.order[seq]);
      ordtbl.ysel[seq] = ordtbl.wsel[seq] =
                         ordtbl.msel[seq] = !ordtbl.status[seq];
   }

   free_dvector(mval);
   free_dvector(pval);
   free_dvector(pval2);

   free_dmatrix(p);
   free_dmatrix(q);

   free_binary_data(&rotfits);

   write_image_descriptor_integer(argv[1],"VOFFSET",0,"Vertical offset");

   write_table(&ordfits,"order");
   inform("Table '%s' created.",ordfits.file.path);

   inform("");
   inform("Two-dimensional fit: ycen = f(xnorm,mnorm)");
   regression_polynomial_proc("order","Ycen","Xnorm,Mnorm",
                    par.echdef_ycen.degarg,
                    "Ysel","Yfit","Yresid","YREGRE",par.echdef_ycen.kappa);
 
   inform("");
   inform("Two-dimensional fit: yfwhm = f(xnorm,mnorm)");
   regression_polynomial_proc("order","Yfwhm","Xnorm,Mnorm",
                    par.echdef_fwhm.degarg,
                    "Wsel","Wfit","Wresid","WREGRE",par.echdef_fwhm.kappa);
 
   inform("");
   inform("Two-dimensional fit: order = f(xnorm,ynorm)");
   regression_polynomial_proc("order","Order","Xnorm,Ynorm",
                    par.echdef_ord.degarg,
                    "Msel","Mfit","Mresid","MREGRE",par.echdef_ord.kappa);

   read_table_header(&ordfits,"order");

   read_regression_block(&ordfits,&ordfits.extend,"YREGRE",&yreg);
   read_regression_block(&ordfits,&ordfits.extend,"WREGRE",&wreg);
   read_regression_block(&ordfits,&ordfits.extend,"MREGRE",&mreg);

   inform("");
   inform("                     SUMMARY");

   inform_table("   Fit        N_tot     N_sel      R.M.S. Error");

   inform("y = f(x,m)    %4d      %4d     %12.3e",ntot,yreg.nsel,yreg.rms);
   inform("w = f(x,m)    %4d      %4d     %12.3e",ntot,wreg.nsel,wreg.rms);
   inform("m = f(x,y)    %4d      %4d     %12.3e",ntot,mreg.nsel,mreg.rms);
   inform("");

   check_memory();

   return(0);
}
