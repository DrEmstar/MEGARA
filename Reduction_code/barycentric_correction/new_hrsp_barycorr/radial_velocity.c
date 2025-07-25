/* ------------------------------------------------------------------------

   Program:  radial_velocity
   Purpose:  Determine the radial velocity between two spectra using the
             cross-correlation technique.

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

radvel_type rvcfg;

int core_bins(float *rptr,int n,int *first,int *last,int *flag)
{
   float *fptr;
   int a,b,m,s;

   for (a=1,fptr=rptr;a<=n;a++,fptr++) if (*fptr != 0.0) break;
   for (b=n,fptr=rptr+n-1;b>=1;b--,fptr--) if (*fptr != 0.0) break;
   m = b - a + 1;
   if (m < 16)
      s = 0x04;
   else
      s = 0;
   *first = a;
   *last = b;
   *flag |= s;
   return(m);
}

void copy_order(fits_file *afits,fits_file *cfits,int ord,
  double *astart,double *astep,double *cstart,double *cstep,int *flag)
{
   int an,am,cn,cm,arow,crow;
   float *aptr,*cptr;

   an = afits->npix[1];
   am = afits->npix[2];

   cn = cfits->npix[1];
   cm = cfits->npix[2];

   if (an > cn) get_error("copy_order: Destination image too small!");

   crow = nint(pixel_coordinate(cfits,2,(double)ord));
   if (crow < 1 || crow > cm)
     get_error("copy_order: Order number out of range!");
   cptr = cfits->pix + (crow-1)*cn;
   memset(cptr,0,cn*sizeof(float));
   cstart[crow] = 0.0;
   cstep[crow] = 0.0;
   flag[crow] = 0;
   if (!rvcfg.sel[ord]) flag[crow] |= 0x01;

   arow = nint(pixel_coordinate(afits,2,(double)ord));
   if (arow >= 1 && arow <= am)
   {
     aptr = afits->pix + (arow-1)*an;
     memmove(cptr,aptr,an*sizeof(float));
     cstart[crow] = astart[arow];
     cstep[crow] = astep[arow];

   }
   else
     flag[crow] |= 0x02;
}

void prepare_order(fits_file *fits,int row,double start,double step,
  double vbar,double *aml,double *bml,int *flag)
{
   int nbins,ord;
   double left,right;
   float *rowptr,*binptr,ftot,fmean,bell;
   int first,last,chunk,ntot,i;
   double rvcorr,gshift;

   nbins = fits->npix[1];
   ord = nint(world_coordinate(fits,2,(double)row));

   rowptr = fits->pix + (row-1)*nbins;

   rvcorr = get_keyword_double(fits,&fits->head,"RVCORR");
   gshift = log_doppler(vbar - rvcorr);

   if ((flag[row] & 0x02) == 0)
   {
     core_bins(rowptr,nbins,&first,&last,&flag[row]);
     aml[row] = exp(start + (first-1)*step - gshift) * (double)ord;
     bml[row] = exp(start + (last-1)*step - gshift) * (double)ord;
   }

   if (flag[row] != 0) return;

   if (rvcfg.ordprep.id == PREP_NONE) return;

   left = log(rvcfg.mlamstart/(double)ord) + gshift;
   right = log(rvcfg.mlamend/(double)ord) + gshift;
   first = nint((left-start)/step) + 1;
   last = nint((right-start)/step) + 1;
   chunk = first-1;
   if (chunk > 0) memset(rowptr,0,chunk*sizeof(float));
   chunk = nbins-last;
   if (chunk > 0) memset(rowptr+last,0,chunk*sizeof(float));

   ntot = core_bins(rowptr,nbins,&first,&last,&flag[row]);
   if (flag[row] != 0) return;

   if (rvcfg.ordprep.id == PREP_FLIP && ntot > 0)
   {
      for (i=first,binptr=rowptr+first-1;i<=last;i++,binptr++)
         *binptr = 1.0 - *binptr;
   }

   if (rvcfg.ordprep.id == PREP_MEAN && ntot > 0)
   {
      for (i=first,binptr=rowptr+first-1,ftot=0.0;i<=last;i++,binptr++)
         ftot += *binptr;
      fmean = ftot/(double)ntot;
      for (i=first,binptr=rowptr+first-1;i<=last;i++,binptr++)
         *binptr -= fmean;
   }

   if (rvcfg.cosbell > 0 && ntot > 2*rvcfg.cosbell)
   {
      for (i=1;i<=rvcfg.cosbell;i++)
      {
         bell = (float)(cos((i-1)/(double)(rvcfg.cosbell-1)*PI)+1.0)*0.5;
         rowptr[last-rvcfg.cosbell+i] *= bell;
         rowptr[first+rvcfg.cosbell-i] *= bell;
      }
   }
}

int main(int argc,char **argv)
{
   fits_file afits,bfits,fitsvar,fitsref,fitsccf,fitstbl;
   char spectrograph[MAX_SPECNAME+1];
   file_name_type tmp;
   char outname[81];
   int afirstord,alastord,anumord,bfirstord,blastord,bnumord;
   int cnumbin,cnumord,cfirstord,clastord;
   int row,ord,nrv,seq;
   double *astart,*astep,*bstart,*bstep;
   double *varstart,*varstep,*refstart,*refstep,*ccfstart,*ccfstep;
   int *varflag,*refflag;
   double *vardat,*refdat,*ccfdat,*radvel;
   double rvcorr,rv_raw,rv,shift;
   double *aml,*bml;
   double maxaml,minbml;
   gaussian gauss;
   ccf_table ccftbl;
   int first_bin,last_bin;
   double wsum;
   double vref,vrel,vvar;

   if (argc != 6)
      get_error("Usage:  radial_velocity  <var> <ref> <out> <vref> <vrel>");

   setbuf(stderr,NULL);

   inform_title("Radial velocity determination");

   vref = atof(argv[4]);
   vrel = atof(argv[5]);
   vvar = vref + vrel;

   read_image(&afits,argv[1]);
   check_image_2d(&afits);

   read_image(&bfits,argv[2]);
   check_image_2d(&bfits);

   get_spectrograph(&afits,&afits.head,spectrograph);

   strcpy(outname,argv[3]);

   anumord = afits.npix[2];
   bnumord = bfits.npix[2];

   astart = dvector(anumord);
   astep = dvector(anumord);
   bstart = dvector(bnumord);
   bstep = dvector(bnumord);

   get_double_array(&afits,&afits.head,"STARTVAL",&astart[1],anumord);
   get_double_array(&afits,&afits.head,"STEPVAL",&astep[1],anumord);
   get_double_array(&bfits,&bfits.head,"STARTVAL",&bstart[1],bnumord);
   get_double_array(&bfits,&bfits.head,"STEPVAL",&bstep[1],bnumord);

   afirstord = nint(world_coordinate(&afits,2,(double)anumord));
   alastord  = nint(world_coordinate(&afits,2,1.0));
   if (alastord-afirstord != anumord-1)
     get_error("Order number mismatch in '%s'!",afits.file.path);

   bfirstord = nint(world_coordinate(&bfits,2,(double)bnumord));
   blastord  = nint(world_coordinate(&bfits,2,1.0));
   if (blastord-bfirstord != bnumord-1)
     get_error("Order number mismatch in '%s'!",bfits.file.path);

   cfirstord = imin(afirstord,bfirstord);
   clastord = imax(alastord,blastord);

   cnumord = clastord-cfirstord+1;
   cnumbin = imax(afits.npix[1],bfits.npix[1]);

   inform("%d orders to process (%d-%d).",cnumord,cfirstord,clastord);

   varstart = dvector(cnumord);
   varstep = dvector(cnumord);
   refstart = dvector(cnumord);
   refstep = dvector(cnumord);
   ccfstart = dvector(cnumord);
   ccfstep = dvector(cnumord);
   varflag = ivector(cnumord);
   refflag = ivector(cnumord);

   clear_dvector(varstart,cnumord);
   clear_dvector(varstep,cnumord);
   clear_dvector(refstart,cnumord);
   clear_dvector(refstep,cnumord);
   clear_dvector(ccfstart,cnumord);
   clear_dvector(ccfstep,cnumord);
   clear_ivector(varflag,cnumord);
   clear_ivector(refflag,cnumord);

   get_radvel_info(&rvcfg,spectrograph);

   create_image(next_tmp_file_name(tmp),&afits.head,-32,2,
                cnumbin,cnumord,1.0,1.0,1.0,1.0,(double)clastord,-1.0);
   read_image(&fitsvar,tmp);
   read_image(&fitsccf,tmp);
   remove_tmp_file("%s.fit",tmp);

   create_image(next_tmp_file_name(tmp),&bfits.head,-32,2,
                cnumbin,cnumord,1.0,1.0,1.0,1.0,(double)clastord,-1.0);
   read_image(&fitsref,tmp);
   remove_tmp_file("%s.fit",tmp);

   for (ord=cfirstord;ord<=clastord;ord++)
   {
      copy_order(&afits,&fitsvar,ord,astart,astep,varstart,varstep,varflag);
      copy_order(&bfits,&fitsref,ord,bstart,bstep,refstart,refstep,refflag);
   }

   write_double_array(&fitsvar,&fitsvar.head,"STARTVAL",&varstart[1],cnumord);
   write_double_array(&fitsvar,&fitsvar.head,"STEPVAL",&varstep[1],cnumord);

   write_double_array(&fitsref,&fitsref.head,"STARTVAL",&refstart[1],cnumord);
   write_double_array(&fitsref,&fitsref.head,"STEPVAL",&refstep[1],cnumord);

   aml = dvector(cnumord);
   bml = dvector(cnumord);

   clear_dvector(aml,cnumord);
   clear_dvector(bml,cnumord);

   for (row=1;row<=cnumord;row++)
   {
      prepare_order
        (&fitsvar,row,varstart[row],varstep[row],vvar,aml,bml,varflag);
      prepare_order
        (&fitsref,row,refstart[row],refstep[row],vref,aml,bml,refflag);
   }

   rvcorr = get_keyword_double(&afits,&afits.head,"RVCORR") -
            get_keyword_double(&bfits,&bfits.head,"RVCORR");

   vardat = dvector(cnumbin);
   refdat = dvector(cnumbin);
   ccfdat = dvector(cnumbin);
   radvel = dvector(cnumord);

   clear_dvector(vardat,cnumbin);
   clear_dvector(refdat,cnumbin);
   clear_dvector(radvel,cnumord);

   create_ccf_table(cnumord,next_tmp_file_name(tmp));
   read_ccf_table(&fitstbl,&ccftbl,tmp);
   remove_tmp_file("%s.fit",tmp);

   inform("Using method: %s",rvcfg.method.txt);

   fprintf(stderr,"Processing row number:   0");
   maxaml = 0.0;
   minbml = 1.0e6;
   for (row=1,seq=0,rv_raw=0.0,wsum=0.0,nrv=0,ord=clastord;
        row<=cnumord;row++,seq++,ord--)
   {
      fprintf(stderr,"\b\b\b%3d",row);

      ccftbl.order[seq] = ord;
      ccftbl.sel[seq] = rvcfg.sel[ord];
      ccftbl.status[seq] = 4000 + 16 * varflag[row] + refflag[row];

      if (ccftbl.status[seq] != 4000) continue;

      ccftbl.status[seq] = 0;

      if (aml[row] > maxaml) maxaml = aml[row];
      if (bml[row] < minbml) minbml = bml[row];

      if (varstep[row] != refstep[row])
        get_error("Same step expected in order %d!",ord);
      ccfstep[row] = refstep[row];

      collect_pixel_row(&fitsvar,row,&vardat[1]);
      collect_pixel_row(&fitsref,row,&refdat[1]);
      correlate_arrays(vardat,refdat,ccfdat,cnumbin,
            varstart[row],refstart[row],ccfstep[row],&ccfstart[row]);
      deposit_pixel_row(&fitsccf,row,&ccfdat[1]);

      ccftbl.radius[seq] = rvcfg.r[ord];
      ccftbl.weight[seq] = rvcfg.w[ord];

      switch(rvcfg.method.id)
      {
         case METHOD_GAUSS:
              shift = locate_gauss_maximum(ccfstart[row],ccfstep[row],
                        ccfdat,cnumbin,rvcfg.checkrad,rvcfg.r[ord],NO_INFO,
                        &first_bin,&last_bin,&gauss);
              break;
         case METHOD_PARABOLA:
              shift = locate_parab_maximum(ccfstart[row],ccfstep[row],
                        ccfdat,cnumbin,rvcfg.checkrad,rvcfg.r[ord],NO_INFO,
                        &first_bin,&last_bin,&gauss);
              break;
         case METHOD_SPLINE:
              shift = locate_spline_maximum(ccfstart[row],ccfstep[row],
                        ccfdat,cnumbin,rvcfg.checkrad,1.0e-6,NO_INFO,&gauss);
              break;
         case METHOD_PEAK:
              shift = locate_peak_maximum(ccfstart[row],ccfstep[row],
                                  ccfdat,cnumbin,rvcfg.checkrad,NO_INFO,
                                  &first_bin,&last_bin,&gauss);
              break;
         case METHOD_CCFMAX:
              shift = locate_ccf_maximum(ccfstart[row],ccfstep[row],
                                          ccfdat,cnumbin,&gauss);
              break;
         default: get_error("Bad RV method (%d)!",rvcfg.method);
      }

      if (gauss.status == 0) radvel[row] =  doppler_vel(shift);

      ccftbl.icen[seq] = gauss.icen;
      ccftbl.xcen[seq] = gauss.xcen;
      ccftbl.status[seq] = gauss.status;

      switch(rvcfg.method.id)
      {
         case METHOD_GAUSS:
            ccftbl.xstart[seq] = first_bin;
            ccftbl.xend[seq] = last_bin;
            ccftbl.base[seq] = gauss.base;
            ccftbl.xfwhm[seq] = gauss.xfwhm;
            ccftbl.chisq[seq] = gauss.chisq;
            ccftbl.rms[seq] = gauss.rms;
            ccftbl.niter[seq] = gauss.niter;
            break;
         case METHOD_PARABOLA:
            ccftbl.xstart[seq] = first_bin;
            ccftbl.xend[seq] = last_bin;
            ccftbl.rms[seq] = gauss.rms;
            break;
         case METHOD_PEAK:
            ccftbl.xstart[seq] = first_bin;
            ccftbl.xend[seq] = last_bin;
            break;
      }

      if (gauss.status != 0) continue;

      ccftbl.rv_raw[seq] = radvel[row];
      ccftbl.rv[seq] = radvel[row] + rvcorr;

      rv_raw += radvel[row] * rvcfg.w[ord];
      wsum += rvcfg.w[ord];
      nrv++;
   }
   fprintf(stderr,"\n");

   if (maxaml > rvcfg.mlamstart || minbml < rvcfg.mlamend)
   {
     inform("* MLAM limits too wide in 'radvel.cfg'");
     inform("* Suggested MLAM limits: %6.0f..%6.0f",
       ceil(maxaml),floor(minbml));
   }

   for (seq=0;seq<cnumord;seq++)
   {
      if (ccftbl.status[seq]<4000 && ccftbl.status[seq]!=0)
         inform("Order %d: Radial velocity determination failed!",
                 ccftbl.order[seq]);
   }

   if (nrv == 0) get_error("No usable orders found!");

   inform("%d orders used.",nrv);

   rv_raw /= wsum;
   rv = rv_raw + rvcorr;

   inform("Corrected radial velocity: %10.5f km/s",rv);

   write_double_array(&fitsccf,&fitsccf.head,"STARTVAL",&ccfstart[1],cnumord);
   write_double_array(&fitsccf,&fitsccf.head,"STEPVAL",&ccfstep[1],cnumord);

   write_double_array(&fitsccf,&fitsccf.head,"RADVEL",&radvel[1],cnumord);
   write_keyword_double(&fitsccf,&fitsccf.head,"RV_RAW",rv_raw,
                        "Raw radial velocity (km/s)");
   write_keyword_double(&fitsccf,&fitsccf.head,"RV",rv,
                        "Corrected radial velocity (km/s)");
   write_keyword_textual(&fitsccf,&fitsccf.head,"RVREF",bfits.file.name,
                        "RV reference image name");

   write_image(&fitsvar,"%s-var",outname);
   inform("Image '%s' created.",fitsvar.file.path);

   write_image(&fitsref,"%s-ref",outname);
   inform("Image '%s' created.",fitsref.file.path);

   write_image(&fitsccf,"%s-ccf",outname);
   inform("Image '%s' created.",fitsccf.file.path);

   write_table(&fitstbl,"%s-tbl",outname);
   inform("Table '%s' created.",fitstbl.file.path);

   free_dvector(vardat);
   free_dvector(refdat);
   free_dvector(ccfdat);
   free_dvector(radvel);
   free_dvector(aml);
   free_dvector(bml);

   free_dvector(astart);
   free_dvector(astep);
   free_dvector(bstart);
   free_dvector(bstep);
   free_dvector(varstart);
   free_dvector(varstep);
   free_dvector(refstart);
   free_dvector(refstep);
   free_dvector(ccfstart);
   free_dvector(ccfstep);
   free_ivector(varflag);
   free_ivector(refflag);

   free_binary_data(&afits);
   free_binary_data(&bfits);

   check_memory();

   return(0);
}
