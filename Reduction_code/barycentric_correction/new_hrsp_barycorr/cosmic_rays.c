/* ------------------------------------------------------------------------

   Program:  cosmic_rays
   Purpose:  Eliminate any cosmic rays from a given stellar image

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

int u,v,vm;
double *area;
double **norm,**bkg,**raw;

void norm_profile(void)
{
   for (v=1,area[u]=0;v<=vm;v++)
   {
      norm[v][u] = raw[v][u] - bkg[v][u];
      area[u] += norm[v][u];
   }

   if (area[u] == 0.0) area[u] = 1.0;

   for (v=1;v<=vm;v++) norm[v][u] /= area[u];
}

int main(int argc,char **argv)
{
   spec_type spec;
   fits_file rfits,bfits,cfits,ordfits;
   param_type par;
   ccd_type ccd;
   int gainset,gain;
   double voffset,bsize;
   int umax,um,xa,xb;
   int *kcen,*nsel,**sel,**inside,*whole;
   double *bin,*ycen,*shift,*pix,*rms,**fit,**del;
   int i,x,y,profile_adjusted;
   float *rptr,*cptr,*bptr;
   regre_block mreg,yreg;
   double **p,maxord,minord,xord;
   int firstord,lastord,nord,pold,pnew,order,box,seq,totbox,pass;
   double *b,*c;
   int mb,mc;
   double mnorm,ymin,ymax,xnorm,value,sigma;

   if (argc != 5)
      get_error("Usage:  cosmic_rays  <raw> <med> <bkg> <out>");

   setbuf(stderr,NULL);

   inform_title("Cosmic ray subtraction");

   read_image(&rfits,argv[1]);
   read_image(&cfits,argv[2]);
   read_image(&bfits,argv[3]);
   read_table_header(&ordfits,"order");

   get_specinfo_fits(&rfits,&rfits.head,&spec);
   get_ccdinfo_fits(&rfits,&rfits.head,&ccd);
   get_parameters(&par,spec.specname);

   voffset = get_keyword_double(&rfits,&rfits.head,"VOFFSET");
   gainset = get_keyword_integer(&rfits,&rfits.head,"GAINSET");

   for (gain=1;gain<=ccd.ngain && ccd.gain[gain].gainset!=gainset;gain++);
   if (gain > ccd.ngain) get_error("Bad GAINSET value (%d)!",gainset);

   bsize = (double)(rfits.npix[1])/(double)par.crfilt_nboxes;
   umax =  (int)floor(bsize);
   if (bsize > umax) umax++;

   vm = 2*par.crfilt_rslit + 1;

   read_regression_block(&ordfits,&ordfits.extend,"MREGRE",&mreg);
   read_regression_block(&ordfits,&ordfits.extend,"YREGRE",&yreg);

   p = dmatrix(1,2);
 
   p[1][2] = normalized_pixel(&rfits,2,1.0+voffset);
   for (i=1,maxord=0.0;i<=rfits.npix[1];i++)
   {
      p[1][1] = normalized_pixel(&rfits,1,(double)i);
      xord = polynomial_val(p,1,2,mreg.deg,mreg.a,mreg.ma);
      if (xord > maxord) maxord = xord;
   }
 
   p[1][2] = normalized_pixel(&rfits,2,(double)rfits.npix[2]+voffset);
   for (i=1,minord=maxord;i<=rfits.npix[1];i++)
   {
      p[1][1] = normalized_pixel(&rfits,1,(double)i);
      xord = polynomial_val(p,1,2,mreg.deg,mreg.a,mreg.ma);
      if (xord < minord) minord = xord;
   }
 
   free_dmatrix(p);

   firstord = (int)floor(minord) + 1;
   lastord = (int)floor(maxord);
   nord = lastord - firstord + 1;

   bin = dvector(umax);
   ycen = dvector(umax);
   shift = dvector(umax);
   area = dvector(umax);
   kcen = ivector(umax);
   whole = ivector(umax);
   pix = dvector(vm);
   rms = dvector(vm);
   nsel = ivector(vm);
   norm = dmatrix(vm,umax);
   bkg = dmatrix(vm,umax);
   raw = dmatrix(vm,umax);
   fit = dmatrix(vm,umax);
   del = dmatrix(vm,umax);
   sel = imatrix(vm,umax);
   inside = imatrix(vm,umax);

   for (v=1;v<=vm;v++) pix[v] = v;

   rptr = rfits.pix;
   cptr = cfits.pix;
   bptr = bfits.pix;
 
   for (i=1;i<=rfits.totpix;i++,rptr++,cptr++,bptr++)
   {
      if (*rptr - *cptr > par.crfilt_nsigma*sqrt(*cptr/ccd.gain[gain].invgain))
         *cptr -= *bptr;
      else
         *cptr = *rptr - *bptr;
   }

   totbox = nord*par.crfilt_nboxes;

   mb = yreg.deg[1] + 1;
   mc = par.crfilt.deg[1] + 1;
   b = dvector(mb);
   c = dvector(mc);
   
   fprintf(stderr,"Cosmic ray detection:   0%%");
   for (order=firstord,seq=0;order<=lastord;order++)
   {
      mnorm = normalized_order(&spec,(double)order);
      accumulate_coefficients(yreg.a,yreg.ma,b,mb,mnorm,2);
      for (box=1;box<=par.crfilt_nboxes;box++,seq++)
      {
         pnew = nint(seq*100.0/(double)(totbox-1));
         if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
         xa = 1 + nint((box-1)*bsize);
         xb = nint(box*bsize);
         um = xb - xa + 1;

         xnorm = normalized_pixel(&rfits,1,(double)xa);
         ymin = ymax = compute_polynomial(b,mb,xnorm) + voffset;
         for (u=1,x=xa;u<=um;u++,x++)
         {
            xnorm = normalized_pixel(&rfits,1,(double)x);
            ycen[u] = compute_polynomial(b,mb,xnorm) + voffset;
            if (ycen[u] < ymin)
               ymin = ycen[u];
            else
               if (ycen[u] > ymax) ymax = ycen[u];
         }

         if (ymax-ymin > 1.0)
         {
            for (u=1;u<=um;u++)
            {
               kcen[u] = nint(ycen[u]);
               shift[u] = ycen[u] - kcen[u];
            }
         }
         else
         {
            for (u=1;u<=um;u++)
            {
               kcen[u] = nint(ycen[1]);
               shift[u] = (double)(u-1)/(double)(um-1) - 0.5;
            }
         }

         for (u=1,x=xa;u<=um;u++,x++)
         {
            whole[u] = 1;
            for (v=1,y=kcen[u]-par.crfilt_rslit;v<=vm;v++,y++)
            {
               inside[v][u] = sel[v][u] = pixel_inside(&rfits,x,y);
               if (inside[v][u])
               {
                  raw[v][u] = (double)collect_pixel_value(&rfits,x,y);
                  bkg[v][u] = (double)collect_pixel_value(&bfits,x,y);
               }
               else
               {
                  raw[v][u] = 0;
                  bkg[v][u] = 0;
                  whole[u] = 0;
               }
            }
            norm_profile();
         }

         for (v=1;v<=vm;v++)
           for (u=1,nsel[v]=0;u<=um;u++) nsel[v] += sel[v][u];

         for (v=1;v<=vm;v++)
         {
            if (nsel[v] > par.crfilt.deg[1]+1)
            {
               best_regression_polynomial_1d(shift,norm[v],um,sel[v],
                               par.crfilt.deg[1],c,
                               &rms[v],&nsel[v],fit[v],del[v],
                               par.crfilt.kappa,NO_INFO);
            }
            else
            {
               for (u=1;u<=um;u++)
               {
                  fit[v][u] = norm[v][u];
                  del[v][u] = 0;
               }
               rms[v] = 0;
            }
         }

         for (u=1,x=xa;u<=um;u++,x++)
         {
            if (!whole[u]) continue;
            for (pass=1;pass<=par.crfilt_npasses;pass++)
            {
               for (v=1,profile_adjusted=0;v<=vm;v++)
               {
                  value = fit[v][u]*area[u]+bkg[v][u];
                  sigma = sqrt(value/ccd.gain[gain].invgain);
                  if (raw[v][u] - value > par.crfilt_nsigma*sigma)
                     raw[v][u] = value, profile_adjusted=1;
               }
               if (!profile_adjusted) break;
               norm_profile();
            }
         }

         for (u=1,x=xa;u<=um;u++,x++)
         {
            for (v=1,y=kcen[u]-par.crfilt_rslit;v<=vm;v++,y++)
            {
               if (!inside[v][u]) continue;
               deposit_pixel_value(&cfits,x,y,(float)(raw[v][u]-bkg[v][u]));
            }
         }
      }
   }
   fprintf(stderr,"\n");

   free_dvector(c);
   free_dvector(b);

   free_imatrix(inside);
   free_imatrix(sel);
   free_dmatrix(del);
   free_dmatrix(fit);
   free_dmatrix(raw);
   free_dmatrix(bkg);
   free_dmatrix(norm);
   free_ivector(nsel);
   free_dvector(rms);
   free_dvector(pix);
   free_ivector(whole);
   free_ivector(kcen);
   free_dvector(area);
   free_dvector(shift);
   free_dvector(ycen);
   free_dvector(bin);

   free_binary_data(&rfits);
   free_binary_data(&bfits);

   write_image(&cfits,argv[4]);

   inform("Image '%s' created.",cfits.file.path);

   check_memory();

   return(0);
}
