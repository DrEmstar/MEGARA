/* ------------------------------------------------------------------------

   Program:  austral_flux
   Purpose:  Create a flux data file in Austral format.

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

#define AUSTRAL_ONE 112

spec_type spec;
ccd_type ccd;
fits_file efits,hfits,afits,bfits;
lin_table atbl,btbl;
regre_block areg,breg;

int austral_num(int ord)
{
  return(AUSTRAL_ONE - ord + 1);
}

int confirm_vacuum(lin_table *tbl)
{
  return(tbl->lam[1] == tbl->vaclam[1]);
}

void interpolate_solution(void)
{
   int i,h,m;
   double s;
   char midtime[MAX_TIME+1];
   time_type autc,butc,cutc;

   if (areg.ma != breg.ma)
      get_error("Same number of coefficients expected "
                "in both dispersion solutions!");

   get_keyword_textual(&afits,&afits.extend,"MIDTIME",midtime,MAX_TIME);
   sscanf(midtime,"%d:%d:%lf",&h,&m,&s);
   set_time_hms(&autc,h,m,s);

   get_keyword_textual(&bfits,&bfits.extend,"MIDTIME",midtime,MAX_TIME);
   sscanf(midtime,"%d:%d:%lf",&h,&m,&s);
   set_time_hms(&butc,h,m,s);

   get_keyword_textual(&efits,&efits.head,"MIDTIME",midtime,MAX_TIME);
   sscanf(midtime,"%d:%d:%lf",&h,&m,&s);
   set_time_hms(&cutc,h,m,s);

   if (butc.frac < autc.frac) get_error("Thorium solutions out of sequence!");

   for (i=1;i<=areg.ma;i++)
     areg.a[i] =  linint(autc.frac,areg.a[i],butc.frac,breg.a[i],cutc.frac);
}

double get_wave(int ord,int col)
{
  double **u;
  double wav,ra,rb,rc;
  int p;

  u = dmatrix(1,2);

  u[1][2] = normalized_order(&spec,ord);
  ra = -1.5;
  rb = 1.5;
  while (rb - ra > 1e-12)
  {
    u[1][1] = (ra + rb) * 0.5;
    p = polynomial_val(u,1,2,areg.deg,areg.a,areg.ma);
    if (p < col)
      ra = u[1][1];
    else
      rb = u[1][1];
  }
  rc = (ra + rb) * 0.5;
  wav = denormalized_mlambda(&spec,rc) / (double)ord;

  free_dmatrix(u);

  return(wav);
}

double air_wavelength(double w)
{
  double s,u,n;

  s = 1.0e8 / (w * w);
  u = 8342.13 + (2406030.0 / (130.0 - s)) + (15997.0 / (38.9 - s));
  n = 1.0 + u * 1.0e-8;

  return(w / n);
}

int main(int argc,char **argv)
{
  double **fdat,**sdat,*sn;
  double *snmin,*snmax,*snmean,*snstd,*snmed;
  int nsn;
  int minord;
  char expdate[MAX_DATE+1];
  double utmid,expsec,geojd,baryjd,barycorr;
  int interpolated;
  int lastord;
  int n,m;
  int ibase;
  int vacuum;
  int gainset,gain;
  path_name_type outpath;
  int len;
  char *eptr;
  file_type flux;
  int row,col,ord;
  double firstwave,lastwave,wav,adu,flx,sig,r,totadu,totflx;
  double *pdat,*wdat,*wdat2,**u;
  int keyseq,esize,hsize,slitsize,nout;

  if (argc != 5 && argc != 6)
    get_error("Usage:  austral_flux  <eimg> <himg> <minord> <tbl>\n"
              "        austral_flux  <eimg> <himg> <minord> <tbla> <tblb>");

  interpolated = (argc == 6);

  read_image(&efits,argv[1]);
  check_image_2d(&efits);

  read_image(&hfits,argv[2]);
  check_image_2d(&hfits);

  minord = atoi(argv[3]);

  lastord = (int)efits.crval[2];
  n = efits.npix[1];
  m = efits.npix[2];

  get_specinfo_fits(&efits,&efits.head,&spec);
  get_ccdinfo_fits(&efits,&efits.head,&ccd);

  inform("Generating extracted flux data in Austral format.");

  esize = hsize = 0;

  keyseq = find_fits_keyword(&efits.head,"SLITSIZE");
  if (keyseq > 0) esize = read_keyword_integer(&efits,&efits.head,keyseq);

  keyseq = find_fits_keyword(&hfits.head,"SLITSIZE");
  if (keyseq > 0) hsize = read_keyword_integer(&hfits,&hfits.head,keyseq);

  if (esize == 0 && hsize == 0)
    get_error("FITS keyword SLITSIZE not found!");

  if (esize > 0 && hsize > 0 && esize != hsize)
    get_error("Input images have different SLITSIZE!");

  slitsize = imax(esize,hsize);

  inform("Extraction slit size: %d pix",slitsize);

  get_keyword_textual(&efits,&efits.head,"EXPDATE",expdate,MAX_DATE);
  utmid = get_keyword_double(&efits,&efits.head,"UT_MID");
  expsec = get_keyword_double(&efits,&efits.head,"EXPSEC");
  geojd = get_keyword_double(&efits,&efits.head,"JD_MID");
  baryjd = get_keyword_double(&efits,&efits.head,"JD_BARY");
  barycorr = get_keyword_double(&efits,&efits.head,"RVCORR") * 1000.0;

  read_lin_table(&afits,&atbl,argv[4]);
  read_regression_block(&afits,&afits.extend,"XREGRE",&areg);
  vacuum = confirm_vacuum(&atbl);

  if (interpolated)
  {
    read_lin_table(&bfits,&btbl,argv[5]);
    read_regression_block(&bfits,&bfits.extend,"XREGRE",&breg);
    if (confirm_vacuum(&btbl) != vacuum)
      get_error("Either vacuum or air wavelengths must be used in "
                "both dispersion solutions!");
    interpolate_solution();
  }

  gainset = get_keyword_integer(&efits,&efits.head,"GAINSET");
  for (gain=1;gain<=ccd.ngain && ccd.gain[gain].gainset!=gainset;gain++);
  if (gain > ccd.ngain) get_error("Bad GAINSET value (%d)!",gainset);
  r = ccd.gain[gain].bnoise;

  strcpy(outpath,efits.file.path);
  len = strlen(outpath);
  if (len < 4) get_error("File name too short: '%s'!",efits.file.path);
  eptr = outpath + len - 4;
  if (strcmp(eptr,".fit") != 0)
    get_error("Extension '.fit' expected in '%s'!",efits.file.path);
  strcpy(eptr,".dat");
  set_file_name(&flux,outpath);

  fdat = dmatrix(m,n);
  sdat = dmatrix(m,n);
  sn = dvector(n);
  snmin = dvector(m);
  snmax = dvector(m);
  snmean = dvector(m);
  snstd = dvector(m);
  snmed = dvector(m);
  pdat = dvector(n);
  wdat = dvector(n);
  wdat2 = dvector(n);
  u = dmatrix(1,2);

  for (row=1;row<=m;row++)
  {
    clear_dvector(sn,n);
    for (col=1,nsn=0;col<=n;col++)
    {
      adu = (double)collect_pixel_value(&efits,col,row);
      flx = adu * ccd.gain[gain].invgain;
      totadu = (double)collect_pixel_value(&hfits,col,row);
      totflx = totadu * ccd.gain[gain].invgain;
      if (flx > 0.0)
         sig = sqrt(totflx + slitsize*r*r);
      else
         sig = 1.0e5;
      fdat[row][col] = flx;
      sdat[row][col] = sig;
      if (flx > 0.0) sn[++nsn] = flx / sig;
    }
    if (nsn > 1)
    {
      statistics_double_array
        (sn,nsn,&snmin[row],&snmax[row],&snmean[row],&snstd[row]);
      snmed[row] = get_median(sn,nsn);
    }
    else
    {
      snmin[row] = sn[1];
      snmax[row] = sn[1];
      snmean[row] = sn[1];
      snstd[row] = 0.0;
      snmed[row] = sn[1];
    }
  }

  get_fopen(&flux,"w");

  fprintf(flux.dat,"! Ident: Hercules spectrum\n");
  fprintf(flux.dat,"! Telescope: Mt. John 1m\n");
  fprintf(flux.dat,"! ObsDate (d m y): %s\n",expdate);
  fprintf(flux.dat,"! ObsTime (UT/h) : %8.6f\n",utmid);
  fprintf(flux.dat,"! ExpTime (s)    : %6.2f\n",expsec);
  fprintf(flux.dat,"! geo JD , bary JD , Barycorr(m/s):\n");
  fprintf(flux.dat,"!  %14.6f  %14.6f %10.3f\n",geojd,baryjd,barycorr);
  for (row=1,ord=lastord;row<=m;row++,ord--)
  {
    if (ord > AUSTRAL_ONE) continue;
    if (ord < minord) continue;
    fprintf(flux.dat,"! Echelle order # %d:S/N: Mean, Median, Min, Max:\n",
      austral_num(ord));
    fprintf(flux.dat,"! %3.1f %3.1f %3.1f %3.1f\n",
      snmean[row],snmed[row],snmin[row],snmax[row]);
  }
  fprintf(flux.dat,"! PX#   WAVELENGTH  FLUX   ERROR         MASK\n");

  for (row=1,ord=lastord,nout=0;row<=m;row++,ord--)
  {
    if (ord > AUSTRAL_ONE) continue;
    if (ord < minord) continue;
    ibase = austral_num(ord) * 10000;
    firstwave = get_wave(ord,0);
    lastwave = get_wave(ord,n+1);
    u[1][2] = normalized_order(&spec,ord);
    for (col=1;col<=n;col++)
    {
      wdat[col] = ((n-col) * firstwave + (col-1) * lastwave) / (n-1);
      u[1][1] = normalized_mlambda(&spec,ord*wdat[col]);
      pdat[col] = polynomial_val(u,1,2,areg.deg,areg.a,areg.ma);
    }
    spline(pdat,wdat,n,1e30,1e30,wdat2);
    for (col=1;col<=n;col++)
    {
      splint(pdat,wdat,wdat2,n,(double)col,&wav);
      if (vacuum) wav = air_wavelength(wav);
      adu = (double)collect_pixel_value(&efits,col,row);
      fprintf(flux.dat,"%07d %12.6f %12.3f %10.3f  0\n",
        ibase+col,wav,fdat[row][col],sdat[row][col]);
    }
    nout++;
  }

  inform("%d orders extracted.",nout);

  get_fclose(&flux);

  inform("ASCII file '%s' created.",flux.path);

  free_dvector(pdat);
  free_dvector(wdat);
  free_dvector(wdat2);
  free_dmatrix(u);
  free_dmatrix(fdat);
  free_dmatrix(sdat);
  free_dvector(sn);
  free_dvector(snmin);
  free_dvector(snmax);
  free_dvector(snmean);
  free_dvector(snstd);
  free_dvector(snmed);

  free_binary_data(&efits);
  free_binary_data(&hfits);
  free_binary_data(&afits);
  if (interpolated) free_binary_data(&bfits);

  check_memory();

  return(0);
}
