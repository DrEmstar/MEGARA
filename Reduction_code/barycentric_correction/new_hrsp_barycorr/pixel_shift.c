/* ------------------------------------------------------------------------

   Program:  pixel_shift
   Purpose:  Determine the shift in pixel units between two extracted
             spectra by cross-correlating individual echelle orders.

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

shift_type sft;

void copy_order(fits_file *afits,fits_file *cfits,int ord,int *flag)
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
   flag[crow] = 0;
   if (!sft.sel[ord]) flag[crow] |= 0x01;

   arow = nint(pixel_coordinate(afits,2,(double)ord));
   if (arow >= 1 && arow <= am)
   {
     aptr = afits->pix + (arow-1)*an;
     memmove(cptr,aptr,an*sizeof(float));
   }
   else
     flag[crow] |= 0x02;
}

void prepare_order(fits_file *fits,int row,int *flag)
{
   int nbins,ntot;
   float *rowptr,*binptr,ftot,fmean,bell;
   int first,last,chunk,i;

   nbins = fits->npix[1];

   rowptr = fits->pix + (row-1)*nbins;

   if (flag[row] != 0) return;

   if (sft.ordprep.id == PREP_NONE) return;

   first = imax(1,sft.xleft);
   last = imin(sft.xright,nbins);

   chunk = first-1;
   if (chunk > 0) memset(rowptr,0,chunk*sizeof(float));
   chunk = nbins-last;
   if (chunk > 0) memset(rowptr+last,0,chunk*sizeof(float));

   ntot = last - first + 1;
   if (ntot < 16) flag[row] |= 0x04;
   if (flag[row] != 0) return;

   if (sft.ordprep.id == PREP_FLIP && ntot > 0)
   {
      for (i=first,binptr=rowptr+first-1;i<=last;i++,binptr++)
         *binptr = 1.0 - *binptr;
   }

   if (sft.ordprep.id == PREP_MEAN && ntot > 0)
   {
      for (i=first,binptr=rowptr+first-1,ftot=0.0;i<=last;i++,binptr++)
         ftot += *binptr;
      fmean = ftot/(double)ntot;
      for (i=first,binptr=rowptr+first-1;i<=last;i++,binptr++)
         *binptr -= fmean;
   }

   if (sft.cosbell > 0 && ntot > 2*sft.cosbell)
   {
      for (i=1;i<=sft.cosbell;i++)
      {
         bell = (float)(cos((i-1)/(double)(sft.cosbell-1)*PI)+1.0)*0.5;
         rowptr[last-sft.cosbell+i] *= bell;
         rowptr[first+sft.cosbell-i] *= bell;
      }
   }
}

int main(int argc,char **argv)
{
   fits_file afits,bfits,fitsvar,fitsref,fitsccf,fitstbl;
   char spectrograph[MAX_SPECNAME+1];
   int anumord,bnumord,cnumord;
   int afirstord,alastord,bfirstord,blastord,cfirstord,clastord;
   int cnumbin;
   int *varflag,*refflag;
   file_name_type tmp;
   gaussian gauss;
   sft_table stbl;
   int row,seq,ord;
   double *vardat,*refdat,*ccfdat;
   double ccfstart,shift,mean_shift,wsum;
   int first_bin,last_bin,count;
   double *shifts;

   if (argc != 4) get_error("Usage:  pixel_shift  <var> <ref> <out>");

   setbuf(stderr,NULL);

   inform_title("Pixel shift determination");

   read_image(&afits,argv[1]);
   check_image_2d(&afits);

   read_image(&bfits,argv[2]);
   check_image_2d(&bfits);

   get_spectrograph(&afits,&afits.head,spectrograph);

   anumord = afits.npix[2];
   bnumord = bfits.npix[2];

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
   get_shift_info(&sft,spectrograph);

   varflag = ivector(cnumord);
   refflag = ivector(cnumord);

   clear_ivector(varflag,cnumord);
   clear_ivector(refflag,cnumord);

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
      copy_order(&afits,&fitsvar,ord,varflag);
      copy_order(&bfits,&fitsref,ord,refflag);
   }

   for (row=1;row<=cnumord;row++)
   {
      prepare_order(&fitsvar,row,varflag);
      prepare_order(&fitsref,row,refflag);
   }

   vardat = dvector(cnumbin);
   refdat = dvector(cnumbin);
   ccfdat = dvector(cnumbin);
   shifts = dvector(cnumord);

   clear_dvector(vardat,cnumbin);
   clear_dvector(refdat,cnumbin);
   clear_dvector(shifts,cnumord);

   create_sft_table(cnumord,next_tmp_file_name(tmp));
   read_sft_table(&fitstbl,&stbl,tmp);
   remove_tmp_file("%s.fit",tmp);

   inform("Using method: %s",sft.method.txt);

   fprintf(stderr,"Processing row number:   0");
   for (row=1,seq=0,ord=clastord,wsum=mean_shift=0.0,count=0;
        row<=cnumord;row++,seq++,ord--)
   {
      fprintf(stderr,"\b\b\b%3d",row);

      stbl.order[seq] = ord;
      stbl.sel[seq] = sft.sel[ord];
      stbl.status[seq] = 4000 + 16 * varflag[row] + refflag[row];

      if (stbl.status[seq] != 4000) continue;

      stbl.status[seq] = 0;

      collect_pixel_row(&fitsvar,row,&vardat[1]);
      collect_pixel_row(&fitsref,row,&refdat[1]);
      correlate_arrays(vardat,refdat,ccfdat,cnumbin,1.0,1.0,1.0,&ccfstart);
      deposit_pixel_row(&fitsccf,row,&ccfdat[1]);

      stbl.radius[seq] = sft.r[ord];
      stbl.weight[seq] = sft.w[ord];

      switch(sft.method.id)
      {
         case METHOD_GAUSS:
              shift = locate_gauss_maximum(ccfstart,1.0,
                        ccfdat,cnumbin,sft.checkrad,sft.r[ord],NO_INFO,
                        &first_bin,&last_bin,&gauss);
              break;
         case METHOD_PARABOLA:
              shift = locate_parab_maximum(ccfstart,1.0,
                        ccfdat,cnumbin,sft.checkrad,sft.r[ord],NO_INFO,
                        &first_bin,&last_bin,&gauss);
              break;
         case METHOD_SPLINE:
              shift = locate_spline_maximum(ccfstart,1.0,
                        ccfdat,cnumbin,sft.checkrad,1.0e-6,NO_INFO,&gauss);
              break;
         case METHOD_PEAK:
              shift = locate_peak_maximum(ccfstart,1.0,
                                  ccfdat,cnumbin,sft.checkrad,NO_INFO,
                                  &first_bin,&last_bin,&gauss);
              break;
         case METHOD_CCFMAX:
              shift = locate_ccf_maximum(ccfstart,1.0,
                                          ccfdat,cnumbin,&gauss);
              break;
         default: get_error("Bad RV method (%d)!",sft.method);
      }

      if (gauss.status == 0) shifts[row] =  shift;

      stbl.icen[seq] = gauss.icen;
      stbl.xcen[seq] = gauss.xcen;
      stbl.status[seq] = gauss.status;

      switch(sft.method.id)
      {
         case METHOD_GAUSS:
            stbl.xstart[seq] = first_bin;
            stbl.xend[seq] = last_bin;
            stbl.base[seq] = gauss.base;
            stbl.xfwhm[seq] = gauss.xfwhm;
            stbl.chisq[seq] = gauss.chisq;
            stbl.rms[seq] = gauss.rms;
            stbl.niter[seq] = gauss.niter;
            break;
         case METHOD_PARABOLA:
            stbl.xstart[seq] = first_bin;
            stbl.xend[seq] = last_bin;
            stbl.rms[seq] = gauss.rms;
            break;
         case METHOD_PEAK:
            stbl.xstart[seq] = first_bin;
            stbl.xend[seq] = last_bin;
            break;
      }

      if (gauss.status != 0) continue;

      stbl.shift[seq] = shift;

      mean_shift += shift * sft.w[ord];
      wsum += sft.w[ord];
      count++;
   }
   fprintf(stderr,"\n");

   for (seq=0;seq<cnumord;seq++)
   {
      if (stbl.status[seq]<4000 && stbl.status[seq]!=0)
         inform("Order %d: Pixel shift determination failed!",
                 stbl.order[seq]);
   }

   if (count == 0) get_error("No usable orders found!");

   inform("%d orders used.",count);

   mean_shift /= wsum;

   inform("Pixel shift: %10.4f pix",mean_shift);

   write_keyword_double(&fitsccf,&fitsccf.head,"CRVAL1",(double)ccfstart,
                        "Coordinate at reference pixel");
   write_double_array(&fitsccf,&fitsccf.head,"SHIFT",&shifts[1],cnumord);
   write_keyword_double(&fitsccf,&fitsccf.head,"XSHIFT",mean_shift,
                        "Mean pixel shift");
   write_keyword_textual(&fitsccf,&fitsccf.head,"SHIFTREF",bfits.file.name,
                        "Reference image name");

   write_image(&fitsvar,"%s-var",argv[3]);
   inform("Image '%s' created.",fitsvar.file.path);

   write_image(&fitsref,"%s-ref",argv[3]);
   inform("Image '%s' created.",fitsref.file.path);

   write_image(&fitsccf,"%s-ccf",argv[3]);
   inform("Image '%s' created.",fitsccf.file.path);

   write_table(&fitstbl,"%s-tbl",argv[3]);
   inform("Table '%s' created.",fitstbl.file.path);

   free_dvector(vardat);
   free_dvector(refdat);
   free_dvector(ccfdat);
   free_dvector(shifts);

   free_ivector(varflag);
   free_ivector(refflag);

   free_binary_data(&afits);
   free_binary_data(&bfits);

   check_memory();

   return(0);
}
