#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "vector.h"
#include "numeric.h"
#include "astro.h"
#include "config.h"
#include "fits.h"
#include "process.h"

char *thar_format[] =
{
   "Seq",      "",            "1J",    "I6",
   "Order",    "",            "1J",    "I5",
   "Airlam",   "Angstrom",    "1D",    "F10.4",
   "Vaclam",   "Angstrom",    "1D",    "F10.4",
   "Waveno",   "1/cm",        "1D",    "F10.4",
   "Species",  "",            "8A",    "A8",
   "Ustart",   "mm",          "1D",    "F8.3",
   "Ucen",     "mm",          "1D",    "F8.3",
   "Uend",     "mm",          "1D",    "F8.3",
   "Vstart",   "mm",          "1D",    "F8.3",
   "Vcen",     "mm",          "1D",    "F8.3",
   "Vend",     "mm",          "1D",    "F8.3",
   "Unorm",    "",            "1D",    "F8.5",
   "Vnorm",    "",            "1D",    "F8.5",
   "Mnorm",    "",            "1D",    "F8.5",
   "Vfit",     "pixel",       "1D",    "F10.4",
   "Vresid",   "pixel",       "1D",    "F10.4",
   "Mfit",     "pixel",       "1D",    "F10.4",
   "Mresid",   "pixel",       "1D",    "F10.4",
   "Vsel",     "",            "1J",    "I4",
   "Msel",     "",            "1J",    "I4",
   "#"
};

char *trans_format[] =
{
   "Seq",      "",            "1J",    "I6",
   "Order",    "",            "1J",    "I5",
   "Airlam",   "Angstrom",    "1D",    "F10.4",
   "Vaclam",   "Angstrom",    "1D",    "F10.4",
   "Waveno",   "1/cm",        "1D",    "F10.4",
   "Species",  "",            "8A",    "A8",
   "Ustart",   "mm",          "1D",    "F8.3",
   "Ucen",     "mm",          "1D",    "F8.3",
   "Uend",     "mm",          "1D",    "F8.3",
   "Vstart",   "mm",          "1D",    "F8.3",
   "Vcen",     "mm",          "1D",    "F8.3",
   "Vend",     "mm",          "1D",    "F8.3",
   "Xstart",   "pixel",       "1D",    "F7.1",
   "Xcenter",  "pixel",       "1D",    "F7.1",
   "Xend",     "pixel",       "1D",    "F7.1",
   "Ystart",   "pixel",       "1D",    "F7.1",
   "Ycenter",  "pixel",       "1D",    "F7.1",
   "Yend",     "pixel",       "1D",    "F7.1",
   "Base",     "ADU",         "1D",    "F10.3",
   "Icen",     "ADU",         "1D",    "F10.3",
   "Xcen",     "pixel",       "1D",    "F8.3",
   "Ycen",     "pixel",       "1D",    "F8.3",
   "Xfwhm",    "pixel",       "1D",    "F8.3",
   "Yfwhm",    "pixel",       "1D",    "F8.3",
   "XYfactor", "",            "1D",    "F8.5",
   "Chisq",    "",            "1D",    "E12.5",
   "Rms",      "ADU",         "1D",    "E12.5",
   "Niter",    "",            "1J",    "I5",
   "Status",   "",            "1J",    "I6",
   "Xsel",     "",            "1J",    "I4",
   "Ysel",     "",            "1J",    "I4",
   "Usel",     "",            "1J",    "I4",
   "Vsel",     "",            "1J",    "I4",
   "Xnorm",    "",            "1D",    "F8.5",
   "Ynorm",    "",            "1D",    "F8.5",
   "Xfit",     "pixel",       "1D",    "F10.4",
   "Xresid",   "pixel",       "1D",    "F10.4",
   "Yfit",     "pixel",       "1D",    "F10.4",
   "Yresid",   "pixel",       "1D",    "F10.4",
   "Unorm",    "",            "1D",    "F8.5",
   "Vnorm",    "",            "1D",    "F8.5",
   "Ufit",     "pixel",       "1D",    "F10.4",
   "Uresid",   "pixel",       "1D",    "F10.4",
   "Vfit",     "pixel",       "1D",    "F10.4",
   "Vresid",   "pixel",       "1D",    "F10.4",
   "#"
};

char *ord_format[] =
{
   "Seq",      "",            "1J",    "I6",
   "Order",    "",            "1J",    "I5",
   "X",        "pixel",       "1J",    "I4",
   "Ystart",   "pixel",       "1D",    "F7.1",
   "Ycenter",  "pixel",       "1D",    "F7.1",
   "Yend",     "pixel",       "1D",    "F7.1",
   "Base",     "ADU",         "1D",    "F10.3",
   "Icen",     "ADU",         "1D",    "F10.3",
   "Ycen",     "pixel",       "1D",    "F8.3",
   "Yfwhm",    "pixel",       "1D",    "F8.3",
   "Chisq",    "",            "1D",    "E12.5",
   "Rms",      "ADU",         "1D",    "E12.5",
   "Niter",    "",            "1J",    "I5",
   "Status",   "",            "1J",    "I6",
   "Ysel",     "",            "1J",    "I4",
   "Wsel",     "",            "1J",    "I4",
   "Msel",     "",            "1J",    "I4",
   "Xnorm",    "",            "1D",    "F8.5",
   "Ynorm",    "",            "1D",    "F8.5",
   "Mnorm",    "",            "1D",    "F8.5",
   "Yfit",     "pixel",       "1D",    "F10.4",
   "Yresid",   "pixel",       "1D",    "F10.4",
   "Wfit",     "pixel",       "1D",    "F10.4",
   "Wresid",   "pixel",       "1D",    "F10.4",
   "Mfit",     "",            "1D",    "F10.4",
   "Mresid",   "",            "1D",    "F10.4",
   "#"
};

char *lin_format[] =
{
   "Seq",      "",            "1J",    "I6",
   "Order",    "",            "1J",    "I5",
   "Airlam",   "Angstrom",    "1D",    "F10.4",
   "Vaclam",   "Angstrom",    "1D",    "F10.4",
   "Waveno",   "1/cm",        "1D",    "F10.4",
   "Species",  "",            "8A",    "A8",
   "Ustart",   "mm",          "1D",    "F8.3",
   "Ucen",     "mm",          "1D",    "F8.3",
   "Uend",     "mm",          "1D",    "F8.3",
   "Vstart",   "mm",          "1D",    "F8.3",
   "Vcen",     "mm",          "1D",    "F8.3",
   "Vend",     "mm",          "1D",    "F8.3",
   "Xstart",   "pixel",       "1D",    "F7.1",
   "Xcenter",  "pixel",       "1D",    "F7.1",
   "Xend",     "pixel",       "1D",    "F7.1",
   "Base",     "ADU",         "1D",    "F10.3",
   "Icen",     "ADU",         "1D",    "F10.3",
   "Xcen",     "pixel",       "1D",    "F8.3",
   "Xfwhm",    "pixel",       "1D",    "F8.3",
   "Chisq",    "",            "1D",    "E12.5",
   "Rms",      "ADU",         "1D",    "E12.5",
   "Niter",    "",            "1J",    "I5",
   "Status",   "",            "1J",    "I6",
   "Xsel",     "",            "1J",    "I4",
   "Lam",      "Angstrom",    "1D",    "F10.4",
   "Mlam",     "",            "1D",    "F11.4",
   "Mnorm",    "",            "1D",    "F8.5",
   "MLnorm",   "",            "1D",    "F8.5",
   "Xfit",     "pixel",       "1D",    "F10.4",
   "Xresid",   "pixel",       "1D",    "F10.4",
   "#"
};

char *back_format[] =
{
   "Seq",      "",            "1J",    "I6",
   "X",        "pixel",       "1J",    "I4",
   "Y",        "pixel",       "1J",    "I4",
   "Sel",      "",            "1J",    "I3",
   "Order",    "",            "1D",    "F5.1",
   "Pixval",   "ADU",         "1D",    "F8.2",
   "Xnorm",    "",            "1D",    "F8.5",
   "Ynorm",    "",            "1D",    "F8.5",
   "Fit",      "ADU",         "1D",    "F10.4",
   "Resid",    "ADU",         "1D",    "F10.4",
   "#"
};

char *ccf_format[] =
{
   "Seq",      "",            "1J",    "I6",
   "Order",    "",            "1J",    "I5",
   "Sel",      "",            "1J",    "I3",
   "Radius",   "pix",         "1J",    "I6",
   "Weight",   "",            "1D",    "F10.5",
   "Xstart",   "pixel",       "1D",    "F7.1",
   "Xend",     "pixel",       "1D",    "F7.1",
   "Base",     "",            "1D",    "F10.5",
   "Icen",     "",            "1D",    "F10.5",
   "Xcen",     "pixel",       "1D",    "F10.4",
   "Xfwhm",    "pixel",       "1D",    "F10.4",
   "Chisq",    "",            "1D",    "E12.5",
   "Rms",      "",            "1D",    "E12.5",
   "Niter",    "",            "1J",    "I5",
   "Status",   "",            "1J",    "I6",
   "RV_raw",   "km/s",        "1D",    "F10.4",
   "RV",       "km/s",        "1D",    "F10.4",
   "#"
};

void set_thar_columns(thar_table *t,fits_file *f)
{
   t->seq = (int *)f->bin;
   t->order = t->seq + f->nrows;
   t->airlam = (double *)(t->order + f->nrows);
   t->vaclam = t->airlam + f->nrows;
   t->waveno = t->vaclam + f->nrows;
   t->species = (char *)(t->waveno + f->nrows);
   t->ustart = (double *)(t->species + 8*f->nrows);
   t->ucen = t->ustart + f->nrows;
   t->uend = t->ucen + f->nrows;
   t->vstart = t->uend + f->nrows;
   t->vcen = t->vstart + f->nrows;
   t->vend = t->vcen + f->nrows;
   t->unorm = t->vend + f->nrows;
   t->vnorm = t->unorm + f->nrows;
   t->mnorm = t->vnorm + f->nrows;
   t->vfit = t->mnorm + f->nrows;
   t->vresid = t->vfit + f->nrows;
   t->mfit = t->vresid + f->nrows;
   t->mresid = t->mfit + f->nrows;
   t->vsel = (int *)(t->mresid + f->nrows);
   t->msel = t->vsel + f->nrows;
}

void set_trans_columns(trans_table *t,fits_file *f)
{
   t->seq = (int *)f->bin;
   t->order = t->seq + f->nrows;
   t->airlam = (double *)(t->order + f->nrows);
   t->vaclam = t->airlam + f->nrows;
   t->waveno = t->vaclam + f->nrows;
   t->species = (char *)(t->waveno + f->nrows);
   t->ustart = (double *)(t->species + 8*f->nrows);
   t->ucen = t->ustart + f->nrows;
   t->uend = t->ucen + f->nrows;
   t->vstart = t->uend + f->nrows;
   t->vcen = t->vstart + f->nrows;
   t->vend = t->vcen + f->nrows;
   t->xstart = t->vend + f->nrows;
   t->xcenter = t->xstart + f->nrows;
   t->xend = t->xcenter + f->nrows;
   t->ystart = t->xend + f->nrows;
   t->ycenter = t->ystart + f->nrows;
   t->yend = t->ycenter + f->nrows;
   t->base = t->yend + f->nrows;
   t->icen = t->base + f->nrows;
   t->xcen = t->icen + f->nrows;
   t->ycen = t->xcen + f->nrows;
   t->xfwhm = t->ycen + f->nrows;
   t->yfwhm = t->xfwhm + f->nrows;
   t->xyfactor = t->yfwhm + f->nrows;
   t->chisq = t->xyfactor + f->nrows;
   t->rms = t->chisq + f->nrows;
   t->niter = (int *)(t->rms + f->nrows);
   t->status = t->niter + f->nrows;
   t->xsel = t->status + f->nrows;
   t->ysel = t->xsel + f->nrows;
   t->usel = t->ysel + f->nrows;
   t->vsel = t->usel + f->nrows;
   t->xnorm = (double *)(t->vsel + f->nrows);
   t->ynorm = t->xnorm + f->nrows;
   t->xfit = t->ynorm + f->nrows;
   t->xresid = t->xfit + f->nrows;
   t->yfit = t->xresid + f->nrows;
   t->yresid = t->yfit + f->nrows;
   t->unorm = t->yresid + f->nrows;
   t->vnorm = t->unorm + f->nrows;
   t->ufit = t->vnorm + f->nrows;
   t->uresid = t->ufit + f->nrows;
   t->vfit = t->uresid + f->nrows;
   t->vresid = t->vfit + f->nrows;
}

void set_ord_columns(ord_table *t,fits_file *f)
{
   t->seq = (int *)f->bin;
   t->order = t->seq + f->nrows;
   t->x = t->order + f->nrows;
   t->ystart = (double *)(t->x + f->nrows);
   t->ycenter = t->ystart + f->nrows;
   t->yend = t->ycenter + f->nrows;
   t->base = t->yend + f->nrows;
   t->icen = t->base + f->nrows;
   t->ycen = t->icen + f->nrows;
   t->yfwhm = t->ycen + f->nrows;
   t->chisq = t->yfwhm + f->nrows;
   t->rms = t->chisq + f->nrows;
   t->niter = (int *)(t->rms + f->nrows);
   t->status = t->niter + f->nrows;
   t->ysel = t->status + f->nrows;
   t->wsel = t->ysel + f->nrows;
   t->msel = t->wsel + f->nrows;
   t->xnorm = (double *)(t->msel + f->nrows);
   t->ynorm = t->xnorm + f->nrows;
   t->mnorm = t->ynorm + f->nrows;
   t->yfit = t->mnorm + f->nrows;
   t->yresid = t->yfit + f->nrows;
   t->wfit = t->yresid + f->nrows;
   t->wresid = t->wfit + f->nrows;
   t->mfit = t->wresid + f->nrows;
   t->mresid = t->mfit + f->nrows;
}

void set_lin_columns(lin_table *t,fits_file *f)
{
   t->seq = (int *)f->bin;
   t->order = t->seq + f->nrows;
   t->airlam = (double *)(t->order + f->nrows);
   t->vaclam = t->airlam + f->nrows;
   t->waveno = t->vaclam + f->nrows;
   t->species = (char *)(t->waveno + f->nrows);
   t->ustart = (double *)(t->species + 8*f->nrows);
   t->ucen = t->ustart + f->nrows;
   t->uend = t->ucen + f->nrows;
   t->vstart = t->uend + f->nrows;
   t->vcen = t->vstart + f->nrows;
   t->vend = t->vcen + f->nrows;
   t->xstart = t->vend + f->nrows;
   t->xcenter = t->xstart + f->nrows;
   t->xend = t->xcenter + f->nrows;
   t->base = t->xend + f->nrows;
   t->icen = t->base + f->nrows;
   t->xcen = t->icen + f->nrows;
   t->xfwhm = t->xcen + f->nrows;
   t->chisq = t->xfwhm + f->nrows;
   t->rms = t->chisq + f->nrows;
   t->niter = (int *)(t->rms + f->nrows);
   t->status = t->niter + f->nrows;
   t->xsel = t->status + f->nrows;
   t->lam = (double *)(t->xsel + f->nrows);
   t->mlam = t->lam + f->nrows;
   t->mnorm = t->mlam + f->nrows;
   t->mlnorm = t->mnorm + f->nrows;
   t->xfit = t->mlnorm + f->nrows;
   t->xresid = t->xfit + f->nrows;
}

void set_back_columns(back_table *t,fits_file *f)
{
   t->seq = (int *)f->bin;
   t->x = t->seq + f->nrows;
   t->y = t->x + f->nrows;
   t->sel = t->y + f->nrows;
   t->order = (double *)(t->sel + f->nrows);
   t->pixval = t->order + f->nrows;
   t->xnorm = t->pixval + f->nrows;
   t->ynorm = t->xnorm + f->nrows;
   t->fit = t->ynorm + f->nrows;
   t->resid = t->fit + f->nrows;
}

void set_ccf_columns(ccf_table *t,fits_file *f)
{
   t->seq = (int *)f->bin;
   t->order = t->seq + f->nrows;
   t->sel = t->order + f->nrows;
   t->radius = t->sel + f->nrows;
   t->weight = (double *)(t->radius + f->nrows);
   t->xstart = t->weight + f->nrows;
   t->xend = t->xstart + f->nrows;
   t->base = t->xend + f->nrows;
   t->icen = t->base + f->nrows;
   t->xcen = t->icen + f->nrows;
   t->xfwhm = t->xcen + f->nrows;
   t->chisq = t->xfwhm + f->nrows;
   t->rms = t->chisq + f->nrows;
   t->niter = (int *)(t->rms + f->nrows);
   t->status = t->niter + f->nrows;
   t->rv_raw = (double *)(t->status + f->nrows);
   t->rv = t->rv_raw + f->nrows;
}

double normalized_pixel(fits_file *fits,int axis,double pix)
{
   return((pix-1.0)/(double)(fits->npix[axis]-1));
}

double normalized_order(herc_type *herc,double ord)
{
   return(((double)herc->lastord-ord)/(double)(herc->lastord-herc->firstord));
}

double normalized_u(herc_type *herc,double u)
{
   return(u*herc->uscale);
}

double normalized_v(herc_type *herc,double v)
{
   return(v*herc->vscale);
}

double normalized_mlambda(herc_type *herc,double mlam)
{
   return((mlam-herc->mlamstart)*herc->mlamscale);
}

void read_image_parameters_fits(param_type *par,fits_file *fits)
{
   int reg,fib;

   reg = get_keyword_integer(fits,&fits->head,"POSITION");
   fib = get_keyword_integer(fits,&fits->head,"FIBRENO");
   get_parameters(par,reg,fib);
}

void read_image_parameters_file(param_type *par,char *name)
{
   fits_file fits;

   read_image_header(&fits,name);
   read_image_parameters_fits(par,&fits);
}

void load_dmatrix_from_table(double **a,int k,fits_file *t,int col)
{
   double *u;
 
   u = dvector(t->nrows);
   read_table_column_double(t,col,u);
   load_dmatrix_column(a,t->nrows,u,k);
   free_dvector(u);
}

void create_thar_table(int nrows,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vcreate_fits_table(thar_format,nrows,name,ap);
   va_end(ap);
}

void read_thar_table(fits_file *fits,thar_table *table,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vread_table(fits,name,ap);
   va_end(ap);
 
   set_thar_columns(table,fits);
}

void create_trans_table(int nrows,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vcreate_fits_table(trans_format,nrows,name,ap);
   va_end(ap);
}

void read_trans_table(fits_file *fits,trans_table *table,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vread_table(fits,name,ap);
   va_end(ap);
 
   set_trans_columns(table,fits);
}

void create_ord_table(int nrows,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vcreate_fits_table(ord_format,nrows,name,ap);
   va_end(ap);
}

void read_ord_table(fits_file *fits,ord_table *table,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vread_table(fits,name,ap);
   va_end(ap);
 
   set_ord_columns(table,fits);
}

void create_lin_table(int nrows,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vcreate_fits_table(lin_format,nrows,name,ap);
   va_end(ap);
}

void read_lin_table(fits_file *fits,lin_table *table,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vread_table(fits,name,ap);
   va_end(ap);
 
   set_lin_columns(table,fits);
}

void create_back_table(int nrows,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vcreate_fits_table(back_format,nrows,name,ap);
   va_end(ap);
}

void read_back_table(fits_file *fits,back_table *table,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vread_table(fits,name,ap);
   va_end(ap);
 
   set_back_columns(table,fits);
}

void create_ccf_table(int nrows,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vcreate_fits_table(ccf_format,nrows,name,ap);
   va_end(ap);
}

void read_ccf_table(fits_file *fits,ccf_table *table,char *name,...)
{
   va_list ap;
 
   va_start(ap,name);
   vread_table(fits,name,ap);
   va_end(ap);
 
   set_ccf_columns(table,fits);
}

void check_pixel_range(int a,int b,int first,int last)
{
   if (a >= b)
      get_error("check_pixel_range: Bad pixel interval specification!");
   if (a < first || b > last)
      get_error("check_pixel_range: Pixel coordinate out of limits!");
}

void get_ccdinfo_fits(fits_file *xfits,ccd_type *ccd)
{
   char ccdname[MAX_CCDNAME+1];

   get_keyword_textual(xfits,&xfits->head,"CCDNAME",ccdname,MAX_CCDNAME);
   load_ccdinfo(ccdname,ccd);
}

void get_ccdinfo_file(char *xname,ccd_type *ccd)
{
   fits_file xfits;

   read_image_header(&xfits,xname);
   get_ccdinfo_fits(&xfits,ccd);
}

void extract_image_proc(char *xname,char *yname,int ax,int ay,int bx,int by,
                        int info)
{
   char *dest;
   fits_file xfits,yfits;
   int start,y,skip,axis;
   fitskey_type key;

   if (info) inform("Extract:     %s.fit --> %s.fit",xname,yname);

   read_image_header(&xfits,xname);
   check_image_2d(&xfits);

   check_pixel_range(ax,bx,1,xfits.npix[1]);
   check_pixel_range(ay,by,1,xfits.npix[2]);

   yfits = xfits;

   set_fits_name(&yfits,yname);

   for (axis=1;axis<=2;axis++)
   {
      yfits.crpix[axis] = 1.0;
      yfits.crval[axis] = 1.0;
      yfits.cdelt[axis] = 1.0;
      sprintf(key,"CRPIX%d",axis);
      write_keyword_double(&yfits,&yfits.head,key,yfits.crpix[1],
                                                      "Reference pixel");
      sprintf(key,"CRVAL%d",axis);
      write_keyword_double(&yfits,&yfits.head,key,yfits.crval[1],
                                        "Coordinate at reference pixel"); 
      sprintf(key,"CDELT%d",axis);
      write_keyword_double(&yfits,&yfits.head,key,yfits.cdelt[1],
                                       "Coordinate increment per pixel"); 
   }

   yfits.npix[1] = bx - ax + 1;
   yfits.npix[2] = by - ay + 1;
   yfits.totpix = yfits.npix[1]*yfits.npix[2];
   set_img_totbinrec(&yfits);
   allocate_whole_binary(&yfits);

   write_keyword_integer(&yfits,&yfits.head,"NAXIS1",yfits.npix[1],
                                              "Number of columns");
   write_keyword_integer(&yfits,&yfits.head,"NAXIS2",yfits.npix[2],
                                              "Number of rows");

   open_fits_file(&xfits);
   start = xfits.binstart+((ay-1)*xfits.npix[1]+ax-1)*xfits.pixsize;
   seek_fits_file(&xfits,start,SEEK_SET);

   skip = (xfits.npix[1]-yfits.npix[1])*xfits.pixsize;

   for (y=ay,dest=yfits.bin;y<=by;y++,dest+=yfits.npix[1]*yfits.pixsize)
   {
      if (y > ay && skip > 0) seek_fits_file(&xfits,skip,SEEK_CUR);
      read_fits_file(&xfits,xfits.pixsize,yfits.npix[1],dest);
   }

   close_fits_file(&xfits);

   write_raw_image(&yfits,yname);
}

int rotation_argument(int angle)
{
   int alpha;
   div_t q;

   q = div(angle,360);
   alpha = q.rem;
   if (alpha < 0) alpha+=360;
   q = div(alpha,90);
   if (q.rem != 0)
      get_error("Only integer multiples of 90 degrees are allowed!");

   return(q.quot);
}

int flip_argument(char *s)
{
   if (strcasecmp(s,"X") == 0) return(-1);
   if (strcasecmp(s,"Y") == 0) return(-2);

   if (strcasecmp(s,"XY") != 0) get_error("Bad axis specification!");

   return(2); /* same as rotation by 180 degrees */
}

void swap_axes(fits_file *fits)
{
   int n1,n2;

   n1 = fits->npix[2];
   n2 = fits->npix[1];

   write_keyword_integer(fits,&fits->head,"NAXIS1",n1,"Number of columns");
   write_keyword_integer(fits,&fits->head,"NAXIS2",n2,"Number of rows");

   fits->npix[1] = n2;
   fits->npix[2] = n1;
}

void rotate_pixels(fits_file *xfits,fits_file *yfits,
                      int x1,int dx,int y1,int dy,int ry)
{
   char *xptr,*yptr,*xstart,*ystart;
   int xstep,ystep,yjump;
   int row,col;

   if (xfits->pixsize != yfits->pixsize)
      get_error("rotate_pixels: Same pixel size expected!");

   xstart = xfits->bin + (x1-1)*xfits->pixsize;
   ystart = yfits->bin + (y1-1)*yfits->pixsize;

   xstep = dx*xfits->pixsize;
   ystep = dy*yfits->pixsize;
   yjump = ry*yfits->pixsize;

   for (row=1,yptr=ystart;row<=xfits->npix[2];row++,yptr+=yjump)
   {
      xstart = xfits->bin +((row-1)*xfits->npix[1]+x1-1)*xfits->pixsize;
      for(col=1,xptr=xstart;col<=xfits->npix[1];col++,xptr+=xstep,yptr+=ystep)
         memcpy(yptr,xptr,xfits->pixsize);
   }
}

void multiply_imgmat(char *amat,char *bmat,char *cmat)
{
   int a[4],b[4],c[4];

   sscanf(amat,"%d,%d,%d,%d",&a[0],&a[1],&a[2],&a[3]);
   sscanf(bmat,"%d,%d,%d,%d",&b[0],&b[1],&b[2],&b[3]);

   c[0] = a[0]*b[0] + a[1]*b[2];
   c[1] = a[0]*b[1] + a[1]*b[3];
   c[2] = a[2]*b[0] + a[3]*b[2];
   c[3] = a[2]*b[1] + a[3]*b[3];

   sprintf(cmat,"%d,%d,%d,%d",c[0],c[1],c[2],c[3]);
}

void rotate_fits_image(char *xname,char *yname,int q)
{
   fits_file xfits,yfits;
   char rotmat[MAX_IMGMAT+1],xmat[MAX_IMGMAT+1],ymat[MAX_IMGMAT+1];

   read_raw_image(&xfits,xname);
   check_image_2d(&xfits);
   get_keyword_textual(&xfits,&xfits.head,"IMGMAT",xmat,MAX_IMGMAT);

   yfits = xfits;
   set_fits_name(&yfits,yname);
   allocate_whole_binary(&yfits);

   switch (q)
   {
      case -1: rotate_pixels(&xfits,&yfits,xfits.npix[1],-1,1,1,0);
               strcpy(rotmat,"-1,0,0,1");
               break;
      case -2: rotate_pixels(&xfits,&yfits,xfits.npix[1],
                                                -1,xfits.totpix,-1,0);
               strcpy(rotmat,"1,0,0,-1");
               break;
      case 0:  rotate_pixels(&xfits,&yfits,1,1,1,1,0);
               strcpy(rotmat,"1,0,0,1");
               break;
      case 1:  swap_axes(&yfits);
               rotate_pixels(&xfits,&yfits,1,1,xfits.npix[2],
                                       xfits.npix[2],-xfits.totpix-1);
               strcpy(rotmat,"0,-1,1,0");
               break;
      case 2:  rotate_pixels(&xfits,&yfits,1,1,xfits.totpix,-1,0);
               strcpy(rotmat,"-1,0,0,-1");
               break;
      case 3:  swap_axes(&yfits);
               rotate_pixels(&xfits,&yfits,xfits.npix[1],-1,1,
                                       xfits.npix[2],-xfits.totpix+1);
               strcpy(rotmat,"0,1,-1,0");
               break;
      default: get_error("rotate_fits_image: "
                         "Bad rotation argument (%d)!",q);
   }

   multiply_imgmat(rotmat,xmat,ymat);
   write_keyword_textual(&yfits,&yfits.head,"IMGMAT",ymat,
                          "Pixel transformation matrix");

   write_raw_image(&yfits,yname);

   free_binary_data(&xfits);
}

void rotate_image_proc(char *xname,char *yname,int angle,int pmod)
{
   if (pmod) inform("Rotate:      %s.fit --> %s.fit",xname,yname);

   rotate_fits_image(xname,yname,rotation_argument(angle));
}

void flip_image_proc(char *xname,char *yname,char *axis,int pmod)
{
   if (pmod) inform("Flip:        %s.fit --> %s.fit",xname,yname);

   rotate_fits_image(xname,yname,flip_argument(axis));
}

void convert_image_proc(char *xname,char *yname)
{
   fits_file xfits,yfits;
   int pix;
   int bzeroval;
   float *dest;
   char *src;

   logfile("convert_image %s %s",xname,yname);
   inform("Convert:     %s.fit ---> %s.fit",xname,yname);

   read_image_header(&xfits,xname);

   bzeroval = get_keyword_double(&xfits,&xfits.head,"BZERO");

   if (xfits.bitpix == -32)
   {
      get_system("cp %s.fit %s.fit",xname,yname);
      return;
   }

   yfits = xfits;
   set_fits_name(&yfits,yname);

   write_keyword_integer(&yfits,&yfits.head,"BITPIX",-32,
                                            "32-bit real format");
   yfits.bitpix = -32;
   yfits.pixsize = 4;
   set_img_totbinrec(&yfits);

   allocate_whole_binary(&xfits);
   allocate_whole_binary(&yfits);

   read_image_pixels(&xfits);

   for (pix=1,src=xfits.bin,dest=yfits.pix;
        pix<=xfits.totpix;
        pix++,src+=xfits.pixsize,dest++)
   {
      switch(xfits.bitpix)
      {
         case   8: *dest = (float)(*src);
                   break;
         case  16: *dest = (float)(*((short *)src))+bzeroval;
                   break;
         case  32: *dest = (float)(*((long *)src));
                   break;
         case -32: *dest = *((float *)src);
                   break;
         case -64: *dest = (float)(*((double *)src));
                   break;
         default:  get_error("Bad BITPIX value (%d) in '%s'!",
                   xfits.bitpix,xfits.file.path);
      }
   }

   write_fits_card(&yfits,&yfits.head,"BZERO   =                    0 / OFFSET FACTOR (DEFAULT=0)");

   write_image(&yfits,yname);

   free_binary_data(&xfits);
}

void subtract_bias(char *name)
{
   fits_file fits;
   double bias;

   inform("BIAS:        %s.fit",name);

   read_image_header(&fits,name);
   if (get_keyword_logical(&fits,&fits.head,"BIASDONE") == 'T')
   {
      inform("Bias already subtracted in '%s'!",fits.file.path);
   }
   else
   {
      bias = get_keyword_double(&fits,&fits.head,"ADU_BIAS");
      add_constant(name,-bias,name);
      write_image_descriptor_logical(name,"BIASDONE",'T',
                                     "Bias has been subtracted");
   }
}

void correlate_2d_proc(char *xname,char *yname,char *zname)
{
   fits_file x,y;
   int axis,n,m,nm;
   float *sig,*ref,**res;
   double xstart[3],ystart[3];

   read_image(&x,xname);
   check_real_image(&x);
   check_image_2d(&x);
   check_fft_image(&x);

   read_image(&y,yname);
   check_real_image(&y);
   check_image_2d(&y);
   check_fft_image(&y);

   check_same_image_size(&x,&y);
   check_same_image_step(&x,&y);

   n = x.npix[1];
   m = x.npix[2];
   nm = n*m;

   for (axis=1;axis<=2;axis++)
   {
      xstart[axis] = get_axis_start(&x,axis);
      ystart[axis] = get_axis_start(&y,axis);
   }

   sig = (float *)x.bin;
   ref = (float *)y.bin;

   sig--;
   ref--;

   res = fmatrix(m,n);

   ccf2d(sig,ref,n,m,res[1]);

   shift_fmatrix(res,n,m,n/2,m/2);

   memcpy(x.bin,&res[1][1],nm*sizeof(float));

   set_fits_name(&x,zname);

   write_keyword_double(&x,&x.head,"CRPIX1",(double)(n/2+1),
                         "Reference Pixel");
   write_keyword_double(&x,&x.head,"CRPIX2",(double)(m/2+1),
                         "Reference Pixel");
   write_keyword_double(&x,&x.head,"CRVAL1",xstart[1]-ystart[1],
                         "Coordinate at reference pixel");
   write_keyword_double(&x,&x.head,"CRVAL2",xstart[2]-ystart[2],
                         "Coordinate at reference pixel");

   write_image(&x,zname);
   update_image_minmax(zname);

   free_fmatrix(res);

   free_binary_data(&y);
}

void regression_polynomial_proc(char *tabname,char *yname,char *xname,
        char *deglist,char *selname,char *fitname,char *residname,
        char *rkey,double kappa)
{
   fits_file fits;
   double **x,*y,*fit,*resid;
   int *sel,axis;
   regre_block reg;

   read_table(&fits,tabname);

   reg.ycol = get_table_column(&fits,yname);
   reg.scol = get_table_column(&fits,selname);
   reg.fcol = get_table_column(&fits,fitname);
   reg.rcol = get_table_column(&fits,residname);
   reg.kappa = kappa;

   reg.ndim = count_items(xname);
   if (reg.ndim != count_items(deglist)) get_error("Bad number of columns!");
   if (reg.ndim > MAX_REGRE) get_error("Too many axes!");

   get_several_table_columns(&fits,xname,reg.xcol,reg.ndim);
   collect_int(deglist,reg.ndim,reg.deg);
   list_integers(reg.deg,reg.ndim,reg.arg,MAX_REGRESSION_ARG);

   reg.ma = number_of_coefficients(reg.deg,reg.ndim);

   x=dmatrix(fits.nrows,reg.ndim);
   y=dvector(fits.nrows);
   sel=ivector(fits.nrows);
   fit=dvector(fits.nrows);
   resid=dvector(fits.nrows);

   for (axis=1;axis<=reg.ndim;axis++)
      load_dmatrix_from_table(x,axis,&fits,reg.xcol[axis]);
   read_table_column_double(&fits,reg.ycol,y);
   read_table_column_integer(&fits,reg.scol,sel);

   best_regression_polynomial(x,fits.nrows,reg.ndim,y,sel,reg.deg,
                              reg.a,reg.ma,&reg.rms,&reg.nsel,
                              fit,resid,reg.kappa,YES_INFO);

   write_table_column_integer(&fits,reg.scol,sel);
   write_table_column_double(&fits,reg.fcol,fit);
   write_table_column_double(&fits,reg.rcol,resid);

   write_regression_block(&fits,&fits.extend,rkey,&reg);

   write_table(&fits,tabname);

   free_dvector(resid);
   free_dvector(fit);
   free_ivector(sel);
   free_dvector(y);
   free_dmatrix(x);
}

int collect_jd_number(char *s)
{
   int jd;
 
   jd = atoi(s);
   if (jd<0 || jd >9999) get_error("Bad Julian Day!");
 
   return(jd);
}
 
int collect_image_number(char *s)
{
   int img;
 
   img = atoi(s);
   if (img<1 || img>9999) get_error("Bad image number!");
 
   return(img);
}
