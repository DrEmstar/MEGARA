/* ------------------------------------------------------------------------

   Program:  transform_echelle
   Purpose:  Compute the coordinate transformation coefficients between
             the fixed 2-D spectrum and a given thorium image.

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

#define MAX_TCOUNT 32

fits_file imgfits,tharfits;
ccd_type ccd;
spec_type spec;
param_type par;
thar_table thartbl;
int fibno;
double imgxc,imgyc;
int zerpix,tcount;
double rmm;

int binning_factor(void)
{
   int n;

   if (ccd.pixsize[1] != ccd.pixsize[2])
      get_error("Square pixels expected!");
   n = (int)floor(spec.binsize / ccd.pixsize[1] + 0.5);
   if (fabs(n * ccd.pixsize[1] - spec.binsize) > 1.0e-12)
      get_error("Bad pixel size!");
   return(n);
}

void rebin_boundaries(int n,int m,int *a,int *b)
{
   int h,k,e;

   h = n >> 1;

   k = h / m;
   while (k+k>spec.bincount) k--;
   if (k < 1) get_error("Image too small, or bin too large!");
   e = k * m;
   *a = h - e + 1;
   *b = h + e; 
}

void extract_boundaries(int n,int *a,int *b)
{
   int h,e;

   if (n < spec.bincount) get_error("Image too small!");

   h = n >> 1;
   e = spec.bincount >> 1;

   *a = h - e + 1;
   *b = h + e; 
}

void unwrap(int *k)
{
   if (*k < zerpix)
      *k = zerpix + *k - 1;
   else
      *k = *k - zerpix + 1;
}

void init_regression_blocks(double ucen,double vcen,double theta,
   regre_block *xreg,regre_block *yreg)
{
   double **uv,*x,*y,*fit,*del;
   int *sel;
   int i,j,k;
   double ct,st,u,v,du,dv,ur,vr;
   regre_block reg;
   int n,m;

   reg.ndim = 2;
   reg.deg[1] = 1;
   reg.deg[2] = 1;
   strcpy(reg.arg,"1,1");
   reg.ma = 4;

   n = 9;
   m = 3;

   uv = dmatrix(n,2);
   x = dvector(n);
   y = dvector(n);
   fit = dvector(n);
   del = dvector(n);
   sel = ivector(n);

   for (i=1;i<=n;i++) sel[i] = 1;

   ct = cos(theta);
   st = sin(theta);

   for (i=0,k=1;i<m;i++)
   {
      for (j=0;j<m;j++,k++)
      {
         u = i * spec.usize/2;
         v = j * spec.vsize/2;
         du = u - ucen;
         dv = v - vcen;
         ur = du * ct + dv * st;
         vr = dv * ct - du * st;
         uv[k][1] = normalized_u(&spec,u);
         uv[k][2] = normalized_v(&spec,v);
         x[k] = imgxc + ur / ccd.pixsize[1];
         y[k] = imgyc + vr / ccd.pixsize[2];
      }
   }

   regression_polynomial
      (uv,n,2,x,sel,reg.deg,reg.a,reg.ma,&reg.rms,&reg.nsel,fit,del,YES_INFO);

   *xreg = reg;

   regression_polynomial
      (uv,n,2,y,sel,reg.deg,reg.a,reg.ma,&reg.rms,&reg.nsel,fit,del,YES_INFO);

   *yreg = reg;

   free_dmatrix(uv);
   free_dvector(x);
   free_dvector(y);
   free_dvector(fit);
   free_dvector(del);
   free_ivector(sel);
}

void align_image(char *imgname,regre_block *xreg,regre_block *yreg)
{
   file_name_type reb,frm,ref,ext,rot,ccf;
   fits_file afits,areal,aimag;
   fits_file bfits,breal,bimag;
   fits_file cfits,creal,cimag;
   int nbin,acol,arow,bcol,brow,ccol,crow;
   int pcol,prow,qcol,qrow,mcol,mrow;
   int seq,maxseq;
   int width,height,hpad,vpad;
   float maxpix[MAX_TCOUNT+1],maxccf;
   double theta[MAX_TCOUNT+1];
   int maxrow[MAX_TCOUNT+1],maxcol[MAX_TCOUNT+1];
   double dx,dy,ct,st;
   double ucen,vcen;

   nbin = binning_factor();
   rebin_boundaries(imgfits.npix[1],nbin,&acol,&bcol);
   rebin_boundaries(imgfits.npix[2],nbin,&arow,&brow);

   next_tmp_file_name(reb);
   get_system("rebin_image %s %s %d %d %d %d %d",
      imgname,reb,nbin,arow,acol,brow,bcol);

   read_image_header(&afits,reb);
   width = afits.npix[1];
   height = afits.npix[2];
   if (width > spec.bincount || height > spec.bincount)
      get_error("Image too large!");
   hpad = (spec.bincount - width) >> 1;
   vpad = (spec.bincount - height) >> 1;
   next_tmp_file_name(frm);
   get_system
      ("frame_image %s %s %d %d %d %d",reb,frm,hpad,hpad,vpad,vpad);

   read_image(&afits,frm);

   areal = afits;
   allocate_whole_binary(&areal);
   aimag = afits;
   allocate_whole_binary(&aimag);

   cfits = afits;
   allocate_whole_binary(&cfits);
   creal = afits;
   allocate_whole_binary(&creal);
   cimag = afits;
   allocate_whole_binary(&cimag);

   normalize_pixel_sum(&afits);
   fft_fits(&afits,&areal,&aimag);

   inform("\nStart coarse alignment with reference thorium image...\n");
   for (seq=1;seq<=tcount;seq++)
   {
      inform("Reference image number: %d",seq);
      read_image_header
        (&bfits,"%s/spec/%s/img/Thorium%02d",hrsp_dir(),spec.specname,seq);
      if (bfits.npix[1] != spec.bincount || bfits.npix[2] != spec.bincount)
         get_error("Bad image size!");
      theta[seq] = get_keyword_double(&bfits,&bfits.head,"ROTATION");
      inform("Rotation angle: %5.2f degrees.",theta[seq]);
      read_image
        (&breal,"%s/spec/%s/img/Real%02d",hrsp_dir(),spec.specname,seq);
      read_image
        (&bimag,"%s/spec/%s/img/Imag%02d",hrsp_dir(),spec.specname,seq);
      complex_product_fits(&areal,&aimag,&breal,&bimag,&creal,&cimag,1);
      ifft_fits(&cfits,&creal,&cimag);
      locate_fits_max(&cfits,&maxpix[seq],&maxcol[seq],&maxrow[seq]);
      unwrap(&maxcol[seq]);
      unwrap(&maxrow[seq]);
      free_binary_data(&breal);
      free_binary_data(&bimag);
      inform("");
   }

   for (seq=1,maxseq=1,maxccf=maxpix[1];seq<=tcount;seq++)
   {
      if (maxpix[seq] > maxccf)
      {
         maxccf = maxpix[seq];
         maxseq = seq;
      }
   }

   dx = -(maxcol[maxseq] - zerpix) * spec.binsize;
   dy = -(maxrow[maxseq] - zerpix) * spec.binsize;

   ct = cos(deg2rad(theta[maxseq]));
   st = sin(deg2rad(theta[maxseq]));

   ucen = spec.uref + dx * ct + dy * st;
   vcen = spec.vref + dy * ct - dx * st;

   inform("Image rotation: %5.2f degrees.",theta[maxseq]);
   inform("Image centre at: u=%4.1f mm, v=%4.1f mm",ucen,vcen);

   free_binary_data(&afits);
   free_binary_data(&areal);
   free_binary_data(&aimag);
   free_binary_data(&cfits);
   free_binary_data(&creal);
   free_binary_data(&cimag);

   inform("\nStart fine alignment with reference thorium image...\n");

   ccol = floor((rmm + ucen - spec.uref)/ccd.pixsize[1] + 1.0);
   crow = floor((rmm + vcen - spec.vref)/ccd.pixsize[2] + 1.0);
   acol = ccol - (spec.bincount >> 1);
   arow = crow - (spec.bincount >> 1);
   bcol = acol + spec.bincount - 1;
   brow = arow + spec.bincount - 1; 
   next_tmp_file_name(ref);
   get_system("extract_image %s/spec/%s/img/Th-%s %s %d %d %d %d",
      hrsp_dir(),spec.specname,ccd.ccdname,ref,acol,arow,bcol,brow);

   extract_boundaries(imgfits.npix[1],&pcol,&qcol);
   extract_boundaries(imgfits.npix[2],&prow,&qrow);
   next_tmp_file_name(ext);
   get_system("extract_image %s %s %d %d %d %d",
      imgname,ext,pcol,prow,qcol,qrow);

   next_tmp_file_name(rot);
   get_system("rotate_arbitrary %s %s %4.2f",ext,rot,-theta[maxseq]);

   next_tmp_file_name(ccf);
   get_system("correlate_images %s %s %s",ref,rot,ccf);
   read_image_header(&cfits,ccf);
   mcol = get_keyword_integer(&cfits,&cfits.head,"MAXCOL");
   mrow = get_keyword_integer(&cfits,&cfits.head,"MAXROW");
   ucen = spec.uref + (acol + mcol - 2) * ccd.pixsize[1] - rmm;
   vcen = spec.vref + (arow + mrow - 2) * ccd.pixsize[2] - rmm;
   inform("Image centre at: u=%5.2f mm, v=%5.2f mm\n",ucen,vcen);

   init_regression_blocks(ucen,vcen,-deg2rad(theta[maxseq]),xreg,yreg);

   remove_tmp_file("%s.fit",reb);
   remove_tmp_file("%s.fit",frm);
   remove_tmp_file("%s.fit",ref);
   remove_tmp_file("%s.fit",ext);
   remove_tmp_file("%s.fit",rot);
   remove_tmp_file("%s.fit",ccf);
}

void transformation(double radius,char *deg,double kappa,
   regre_block *xreg,regre_block *yreg,int *ntot)
{
   fits_file afits,transfits;
   trans_table transtbl;
   file_name_type tbl;
   double xrad,yrad;
   int seq;
   double xfwhm,yfwhm;
   double percent;
   int pnew,pold;
   gaussian g;
   double dx,dy,rho;
   double **uv;

   inform("\nStart transformation inside a radius of R=%3.1f pix\n",radius);

   read_image(&afits,imgfits.file.path);

   xrad = par.trans_rad[1] * 0.001 / ccd.pixsize[1];
   yrad = par.trans_rad[2] * 0.001 / ccd.pixsize[2];

   create_trans_table(tharfits.nrows,next_tmp_file_name(tbl));
   read_trans_table(&transfits,&transtbl,tbl);
   remove_tmp_file("%s.fit",tbl);

   write_keyword_textual(&transfits,&transfits.extend,"TRANSIMG",
            imgfits.file.name,"Image used for coordinate transformation");

   copy_integer_array(thartbl.order,tharfits.nrows,transtbl.order);
   copy_double_array(thartbl.airlam,tharfits.nrows,transtbl.airlam);
   copy_double_array(thartbl.vaclam,tharfits.nrows,transtbl.vaclam);
   copy_double_array(thartbl.waveno,tharfits.nrows,transtbl.waveno);
   copy_string_array(thartbl.species,tharfits.nrows,8,transtbl.species);
   copy_double_array(thartbl.ustart,tharfits.nrows,transtbl.ustart);
   copy_double_array(thartbl.ucen,tharfits.nrows,transtbl.ucen);
   copy_double_array(thartbl.uend,tharfits.nrows,transtbl.uend);
   copy_double_array(thartbl.vstart,tharfits.nrows,transtbl.vstart);
   copy_double_array(thartbl.vcen,tharfits.nrows,transtbl.vcen);
   copy_double_array(thartbl.vend,tharfits.nrows,transtbl.vend);

   uv = dmatrix(1,2);

   for (seq=0;seq<transfits.nrows;seq++)
   {
      uv[1][1] = normalized_u(&spec,transtbl.ucen[seq]);
      uv[1][2] = normalized_v(&spec,transtbl.vcen[seq]);

      transtbl.xcenter[seq] =
         polynomial_val(uv,1,2,xreg->deg,xreg->a,xreg->ma);
      transtbl.xstart[seq] = transtbl.xcenter[seq] - xrad;
      transtbl.xend[seq] = transtbl.xcenter[seq] + xrad;

      transtbl.ycenter[seq] =
         polynomial_val(uv,1,2,yreg->deg,yreg->a,yreg->ma);
      transtbl.ystart[seq] = transtbl.ycenter[seq] - yrad;
      transtbl.yend[seq] = transtbl.ycenter[seq] + yrad;
   }

   free_dmatrix(uv);

   xfwhm = spec.fib[fibno].dx / ccd.pixsize[1];
   yfwhm = spec.fib[fibno].dy / ccd.pixsize[2];

   percent = 100.0/(double)(transfits.nrows-1);
   fprintf(stderr,"Locating thorium lines:   0%%");
   for (seq=0,pold=-1,*ntot=0;seq<transfits.nrows;seq++)
   {
      dx = transtbl.xcenter[seq] - imgxc;
      dy = transtbl.ycenter[seq] - imgyc;
      rho = sqrt(dx*dx + dy*dy);
      if (rho < radius)
      {
         locate_gaussian_2d(afits.pix,afits.npix[1],afits.npix[2],
            get_axis_start(&afits,1),get_axis_start(&afits,2),
            afits.cdelt[1],afits.cdelt[2],
            transtbl.xstart[seq],transtbl.ystart[seq],
            transtbl.xcenter[seq],transtbl.ycenter[seq],
            transtbl.xend[seq],transtbl.yend[seq],xfwhm,yfwhm,&g);
      }
      else
      {
         clear_gaussian(&g);
         g.status = 5000;
      }
      transtbl.base[seq] = g.base;
      transtbl.icen[seq] = g.icen;
      transtbl.xcen[seq] = g.xcen;
      transtbl.ycen[seq] = g.ycen;
      transtbl.xfwhm[seq] = g.xfwhm;
      transtbl.yfwhm[seq] = g.yfwhm;
      transtbl.xyfactor[seq] = g.xyfactor;
      transtbl.chisq[seq] = g.chisq;
      transtbl.rms[seq] = g.rms;
      transtbl.niter[seq] = g.niter;
      transtbl.status[seq] = g.status;
      if (g.status == 0) (*ntot)++;
      pnew = nint((double)seq*percent);
      if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
   }
   fprintf(stderr,"\n");
   inform("%d lines detected.",*ntot);

   for (seq=0;seq<transfits.nrows;seq++)
   {
      transtbl.unorm[seq] = normalized_u(&spec,transtbl.ucen[seq]);
      transtbl.vnorm[seq] = normalized_v(&spec,transtbl.vcen[seq]);
      transtbl.xnorm[seq] = normalized_pixel(&imgfits,1,transtbl.xcen[seq]);
      transtbl.ynorm[seq] = normalized_pixel(&imgfits,2,transtbl.ycen[seq]);
      transtbl.xsel[seq] = transtbl.ysel[seq] = 
         transtbl.usel[seq] = transtbl.vsel[seq] = !transtbl.status[seq];
   }

   write_table(&transfits,"transform");
   inform("Table '%s' created.",transfits.file.path);

   inform("\nTwo-dimensional fit: xcen = f(unorm,vnorm)");
   regression_polynomial_proc("transform","Xcen","Unorm,Vnorm",
      deg,"Xsel","Xfit","Xresid","XREGRE",kappa);

   inform("\nTwo-dimensional fit: ycen = f(unorm,vnorm)");
   regression_polynomial_proc("transform","Ycen","Unorm,Vnorm",
      deg,"Ysel","Yfit","Yresid","YREGRE",kappa);

   inform("\nTwo-dimensional fit: ucen = f(xnorm,ynorm)");
   regression_polynomial_proc("transform","Ucen","Xnorm,Ynorm",
      deg,"Usel","Ufit","Uresid","UREGRE",kappa);

   inform("\nTwo-dimensional fit: vcen = f(xnorm,ynorm)");
   regression_polynomial_proc("transform","Vcen","Xnorm,Ynorm",
      deg,"Vsel","Vfit","Vresid","VREGRE",kappa);

   read_table_header(&transfits,"transform");
   read_regression_block(&transfits,&transfits.extend,"XREGRE",xreg);
   read_regression_block(&transfits,&transfits.extend,"YREGRE",yreg);

   free_binary_data(&afits);

   inform("\nTransformation complete.");
}

int main(int argc,char **argv)
{
   regre_block xreg,yreg,ureg,vreg;
   fits_file transfits;
   double imgr;
   int ntot;

   if (argc != 2) get_error("Usage:  transform_echelle  <img>");

   setbuf(stderr,NULL);

   inform_title("Coordinate transformation");

   read_image_header(&imgfits,argv[1]);
   inform("Thorium image: %s",imgfits.file.name);

   get_ccdinfo_fits(&imgfits,&imgfits.head,&ccd);
   get_specinfo_fits(&imgfits,&imgfits.head,&spec);

   inform("CCD: %s",ccd.ccdname);
   inform("Spectrograph: %s",spec.specname);

   fibno = get_fibre_number(&imgfits);
   inform("Fibre: %d",fibno);

   get_parameters(&par,spec.specname);

   read_thar_table
     (&tharfits,&thartbl,"%s/spec/%s/tbl/thar",hrsp_dir(),spec.specname);

   zerpix = spec.bincount / 2 + 1;
   rmm = spec.bincount * spec.binsize * 0.5;
   tcount = count_thorium_ref(spec.specname);
   if (tcount > MAX_TCOUNT) get_error("Too many reference thorium images!");

   imgxc = (1.0 + imgfits.npix[1]) * 0.5;
   imgyc = (1.0 + imgfits.npix[2]) * 0.5;
   imgr = sqrt(imgxc*imgxc + imgyc*imgyc);

   align_image(argv[1],&xreg,&yreg);

   if (imgr > 256.0) transformation
      (256.0,"1,1",par.trans.kappa,&xreg,&yreg,&ntot);
   if (imgr > 512.0 ) transformation
      (512.0,"1,1",par.trans.kappa,&xreg,&yreg,&ntot);
   if (imgr > 1024.0 ) transformation
      (1024.0,"1,1",par.trans.kappa,&xreg,&yreg,&ntot);

   transformation
      (imgr,par.trans.degarg,par.trans.kappa,&xreg,&yreg,&ntot);

   free_binary_data(&tharfits);

   read_table_header(&transfits,"transform");

   read_regression_block(&transfits,&transfits.extend,"XREGRE",&xreg);
   read_regression_block(&transfits,&transfits.extend,"YREGRE",&yreg);
   read_regression_block(&transfits,&transfits.extend,"UREGRE",&ureg);
   read_regression_block(&transfits,&transfits.extend,"VREGRE",&vreg);

   inform("");
   inform("                     SUMMARY");

   inform_table("   Fit        N_tot     N_sel      R.M.S. Error");

   inform("x = f(u,v)    %4d      %4d     %12.3e",ntot,xreg.nsel,xreg.rms);
   inform("y = f(u,v)    %4d      %4d     %12.3e",ntot,yreg.nsel,yreg.rms);
   inform("u = f(x,y)    %4d      %4d     %12.3e",ntot,ureg.nsel,ureg.rms);
   inform("v = f(x,y)    %4d      %4d     %12.3e",ntot,vreg.nsel,vreg.rms);
   inform("");

   check_memory();

   return(0);
}
