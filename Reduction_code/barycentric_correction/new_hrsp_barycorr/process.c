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

static keydef_type keydef_list[] =
{
   /*  0 */ {"OBJECT_NAME",            TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  1 */ {"EXPOSURE_TYPE",          TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  2 */ {"OBSERVATION_DATE",       TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  3 */ {"OBSERVATION_START_TIME", TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  4 */ {"EXPOSURE_SECONDS",       TYPEDBL, 0, "", 0, 0.0, '\x00'},
   /*  5 */ {"FLUX_WEIGHTED_MID_SEC",  TYPEDBL, 0, "", 0, 0.0, '\x00'},
   /*  6 */ {"ADU_BIAS",               TYPEDBL, 0, "", 0, 0.0, '\x00'},
   /*  7 */ {"CCD_CAMERA_NAME",        TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  8 */ {"SPECTROGRAPH_NAME",      TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  9 */ {"CCD_GAIN_SETTING",       TYPEINT, 0, "", 0, 0.0, '\x00'},
   /* 10 */ {"FIBRE_NUMBER",           TYPESTR, 0, "", 0, 0.0, '\x00'},
   /* 11 */ {"",                       0,       0, "", 0, 0.0, '\x00'}
};

char insert_title[] =
   "-------------------------- Inserted keywords --------------------------";

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
   "Lam",      "Angstrom",    "1D",    "F10.4",
   "Mlam",     "",            "1D",    "F11.4",
   "Unorm",    "",            "1D",    "F8.5",
   "Vnorm",    "",            "1D",    "F8.5",
   "Mnorm",    "",            "1D",    "F8.5",
   "MLnorm",   "",            "1D",    "F8.5",
   "Vfit",     "pixel",       "1D",    "F10.4",
   "Vresid",   "pixel",       "1D",    "F10.4",
   "Mfit",     "pixel",       "1D",    "F10.4",
   "Mresid",   "pixel",       "1D",    "F10.4",
   "Ufit",     "pixel",       "1D",    "F10.4",
   "Uresid",   "pixel",       "1D",    "F10.4",
   "Vsel",     "",            "1J",    "I4",
   "Msel",     "",            "1J",    "I4",
   "Usel",     "",            "1J",    "I4",
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
   "RV_raw",   "km/s",        "1D",    "F10.5",
   "RV",       "km/s",        "1D",    "F10.5",
   "#"
};

char *sft_format[] =
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
   "Shift",    "pixel",       "1D",    "F10.4",
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
   t->lam = t->vend + f->nrows;
   t->mlam = t->lam + f->nrows;
   t->unorm = t->mlam + f->nrows;
   t->vnorm = t->unorm + f->nrows;
   t->mnorm = t->vnorm + f->nrows;
   t->mlnorm = t->mnorm + f->nrows;
   t->vfit = t->mlnorm + f->nrows;
   t->vresid = t->vfit + f->nrows;
   t->mfit = t->vresid + f->nrows;
   t->mresid = t->mfit + f->nrows;
   t->ufit = t->mresid + f->nrows;
   t->uresid = t->ufit + f->nrows;
   t->vsel = (int *)(t->uresid + f->nrows);
   t->msel = t->vsel + f->nrows;
   t->usel = t->msel + f->nrows;
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

void set_sft_columns(sft_table *t,fits_file *f)
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
   t->shift = (double *)(t->status + f->nrows);
}

double normalized_pixel(fits_file *fits,int axis,double pix)
{
   return((pix-1.0)/(double)(fits->npix[axis]-1));
}

double normalized_order(spec_type *spec,double ord)
{
   return(((double)spec->lastord-ord)/(double)(spec->lastord-spec->firstord));
}

double normalized_u(spec_type *spec,double u)
{
   return(u*spec->uscale);
}

double normalized_v(spec_type *spec,double v)
{
   return(v*spec->vscale);
}

double normalized_mlambda(spec_type *spec,double mlam)
{
   return((mlam-spec->mlamstart)*spec->mlamscale);
}

double denormalized_mlambda(spec_type *spec,double mlnorm)
{
   return(spec->mlamstart + mlnorm * spec->mlamrange);
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

void create_sft_table(int nrows,char *name,...)
{
   va_list ap;

   va_start(ap,name);
   vcreate_fits_table(sft_format,nrows,name,ap);
   va_end(ap);
}

void read_sft_table(fits_file *fits,sft_table *table,char *name,...)
{
   va_list ap;

   va_start(ap,name);
   vread_table(fits,name,ap);
   va_end(ap);

   set_sft_columns(table,fits);
}

void check_pixel_range(int a,int b,int first,int last)
{
   if (a >= b)
      get_error("check_pixel_range: Bad pixel interval specification!");
   if (a < first || b > last)
      get_error("check_pixel_range: Pixel coordinate out of limits!");
}

void get_spectrograph(fits_file *xfits,fits_head *head,char *spectrograph)
{
   int row;

   row = find_fits_keyword(head,"SPECNAME");
   if (row < 1) row = find_fits_keyword(head,"INSTRUME");
   if (row < 1) row = find_fits_keyword(head,"SPECGRPH");
   if (row < 1) get_error
     ("get_spectrograph: Cannot find the spectrograph name in '%s'!",
     xfits->file.path);

  read_keyword_textual(xfits,&xfits->head,row,spectrograph,MAX_SPECNAME);
  get_lower_case(spectrograph);
}

void get_ccdinfo_fits(fits_file *xfits,fits_head *head,ccd_type *ccd)
{
   char spectrograph[MAX_SPECNAME+1];
   char ccdname[MAX_CCDNAME+1];

   get_spectrograph(xfits,head,spectrograph);
   get_keyword_textual(xfits,head,"CCDNAME",ccdname,MAX_CCDNAME);
   load_ccd_info(ccdname,spectrograph,ccd);
}

void get_ccdinfo_file(char *xname,ccd_type *ccd)
{
   fits_file xfits;
   fits_head *head;

   head = read_fits_header(&xfits,xname);
   get_ccdinfo_fits(&xfits,head,ccd);
}

void get_specinfo_fits(fits_file *xfits,fits_head *head,spec_type *spec)
{
   char specname[MAX_SPECNAME+1];

   get_keyword_textual(xfits,head,"SPECNAME",specname,MAX_SPECNAME);
   load_spec_info(specname,spec);
}

void get_specinfo_file(char *xname,spec_type *spec)
{
   fits_file xfits;
   fits_head *head;

   head = read_fits_header(&xfits,xname);
   get_specinfo_fits(&xfits,head,spec);
}

int standard_name_ok(char *filename)
{
  if (strlen(filename) != 12) return(0);
  if (strcasecmp(filename+8,".fit") != 0) return(0);
  if (!alpha_ok(*filename)) return(0);
  return(numeric_string_ok(filename+1,7));
}

int kiwispec_name_ok(char *filename)
{
  int namelen,i;
  char buffer[29];

  namelen = strlen(filename);
  if (namelen < 28 || namelen > 64) return(0);
  for (i=0;i<namelen;i++)
    if (filename[i] < 0x20 || filename[i] > 0x7E) return(0);
  strcpy(buffer,filename+namelen-28);
  get_lower_case(buffer);

  if (strncmp(buffer,"_kiwispec_",10) != 0) return(0);

  if (!numeric_string_ok(buffer+10,4)) return(0);
  if (month_number(buffer+14) < 1) return(0);
  if (!numeric_string_ok(buffer+17,2)) return(0);
  if (buffer[19] != '_') return(0);

  if (!numeric_string_ok(buffer+20,4)) return(0);

  if (strcmp(buffer+24,".fit") != 0) return(0);

  return(1);
}

void check_standard_name(char *filename)
{
  if (!standard_name_ok(filename))
    get_error("Invalid standard FITS file name '%s'!",filename);
}

void check_kiwispec_name(char *filename)
{
  if (!kiwispec_name_ok(filename))
    get_error("Invalid KiwiSpec file name '%s'!",filename);
}

int check_file_type(char *filename)
{
  if (standard_name_ok(filename)) return(STANDARD_TYPE);
  if (kiwispec_name_ok(filename)) return(KIWISPEC_TYPE);
  return(INVALID_TYPE);
}

void parse_standard_name(char *filename,int *mjd,int *imgno)
{
  char buffer[8];

  check_standard_name(filename);
  memmove(buffer,filename+1,7);
  buffer[7] = 0;
  *imgno = atoi(buffer+4);
  buffer[4] = 0;
  *mjd = atoi(buffer);
}

void parse_kiwispec_name(char *filename,int *mjd,int *imgno)
{
  char buffer[15];
  int day,month,year;

  check_kiwispec_name(filename);
  memmove(buffer,filename+strlen(filename)-18,14);
  buffer[14] = 0;
  *imgno = atoi(buffer+10);
  buffer[9] = 0;
  day = atoi(buffer+7);
  buffer[7] = 0;
  month = month_number(buffer+4);
  buffer[4] = 0;
  year = atoi(buffer);
  *mjd = (int)(date_to_jd(day,month,year) % 10000L);
}

void kiwispec_object(char *filename,char *object)
{
  check_kiwispec_name(filename);
  while(*filename != '_') *(object++) = *(filename++);
  *object = 0;
}

void kiwispec_exposure_type(char *object,char *exptype)
{
  file_name_type buffer;

  if (strlen(object) > MAX_FILE_NAME)
    get_error("kiwispec_exposure_type: Object name too long: '%s'!",object);

  strcpy(buffer,object);
  get_lower_case(buffer);
  strcpy(exptype,"Stellar");
  if (strstr(buffer,"thorium") != NULL) strcpy(exptype,"Thorium Light");
  if (strstr(buffer,"white") != NULL) strcpy(exptype,"White Light");
}

void fix_kiwispec_header(fits_file *fits,char *filename)
{
  int naxis,naxis3;
  int seq,cam,par26,par27,mode,bias,msecs;
  char buffer[CARD_BYTES+1];
  int day,month,year,hour,min;
  double sec,expt,expm;
  double tstart,tmid;
  char object[128],exptype[32];

  seq = get_fits_keyword(fits,&fits->head,"NAXIS");
  naxis = read_keyword_integer(fits,&fits->head,seq);
  if (naxis != 2)
  {
    if (naxis != 3) get_error("fix_kiwispec_header: "
        "Invalid keyword NAXIS = '%d' in '%s'!",naxis,filename);
    replace_value_field(&fits->head,seq,"2");
    seq = get_fits_keyword(fits,&fits->head,"NAXIS3");
    naxis3 = read_keyword_integer(fits,&fits->head,seq);
    if (naxis3 != 1) get_error("fix_kiwispec_header: "
        "Invalid keyword NAXIS3 = '%d' in '%s'!",naxis3,filename);
    comment_out_fits_card(&fits->head,seq);
  }

  ensure_comment_card(fits,&fits->head,insert_title);

  seq = find_fits_keyword(&fits->head,"FILENAME");
  if (seq < 1)
    write_keyword_textual(fits,&fits->head,"FILENAME",filename,"");
  else
  {
    read_keyword_textual(fits,&fits->head,seq,buffer,CARD_BYTES);
    if (strcasecmp(buffer,filename) != 0)
      get_error("fix_kiwispec_header: "
        "Invalid keyword FILENAME = '%s' in '%s'!",buffer,filename);
  }

  kiwispec_object(filename,object);
  kiwispec_exposure_type(object,exptype);

  seq = find_fits_keyword(&fits->head,"OBJECT");
  if (seq < 1) write_keyword_textual
     (fits,&fits->head,"OBJECT",object,"Object name");

  seq = find_fits_keyword(&fits->head,"KIWIEXPT");
  if (seq < 1)
  {
    seq = find_fits_keyword(&fits->head,"IMGTYPE");
    if (seq > 0)
    {
      read_keyword_textual(fits,&fits->head,seq,buffer,CARD_BYTES);
      get_lower_case(buffer);
      strcpy(exptype,"Other");
      if (strstr(buffer,"light") != NULL) strcpy(exptype,"Stellar");
      if (strstr(buffer,"white") != NULL) strcpy(exptype,"White Light");
      if (strstr(buffer,"th-ar") != NULL) strcpy(exptype,"Thorium Light");
    }
    write_keyword_textual
      (fits,&fits->head,"KIWIEXPT",exptype,"Exposure type");
  }

  seq = get_fits_keyword(fits,&fits->head,"DATE-OBS");
  read_keyword_textual(fits,&fits->head,seq,buffer,CARD_BYTES);
  if (!date_string_ok(buffer,&year,&month,&day)) get_error
    ("fix_kiwispec_header: Invalid keyword DATE-OBS = '%s' in '%s'!",
    buffer,filename);
  buffer[10] = 0;
  comment_out_fits_card(&fits->head,seq);
  write_keyword_textual
      (fits,&fits->head,"DATE-OBS",buffer,"UTC Date of observation");

  seq = find_fits_keyword(&fits->head,"REC-STRT");
  if (seq < 1)
  {
    get_keyword_textual(fits,&fits->head,"TIME",buffer,CARD_BYTES);
    if (!time_string_ok(buffer,2,&hour,&min,&sec)) get_error
      ("fix_kiwispec_header: Invalid keyword TIME = '%s' in '%s'!",
      buffer,filename);
    buffer[11] = 0;
    write_keyword_textual
      (fits,&fits->head,"REC-STRT",buffer,"UTC Time of exposure start");
  }
  else
  {
    read_keyword_textual(fits,&fits->head,seq,buffer,CARD_BYTES);
    if (!time_string_ok(buffer,2,&hour,&min,&sec)) get_error
      ("fix_kiwispec_header: Invalid keyword REC-STRT = '%s' in '%s'!",
      buffer,filename);
  }

  tstart = build_decsec_plus(hour,min,sec);

  seq = find_fits_keyword(&fits->head,"EXPTIME");
  if (seq > 0)
    expt = read_keyword_double(fits,&fits->head,seq);
  else
  {
    msecs = get_keyword_integer(fits,&fits->head,"PARAM24");
    expt = (double)msecs / 1000.0;
    write_keyword_double
       (fits,&fits->head,"EXPTIME",expt,"Exposure time (seconds)");
  }

  seq = find_fits_keyword(&fits->head,"MIDEXP");
  if (seq < 1)
  {
    expm = expt / 2.0;
    seq = find_fits_keyword(&fits->head,"EXPMIDPT");
    if (seq > 0)
    {
      read_keyword_textual(fits,&fits->head,seq,buffer,CARD_BYTES);
      if (!date_string_ok(buffer,&year,&month,&day)) get_error
        ("fix_kiwispec_header: Invalid keyword EXPMIDPT = '%s' in '%s'!",
        buffer,filename);
      if (!time_string_ok(buffer+11,-1,&hour,&min,&sec)) get_error
        ("fix_kiwispec_header: Invalid keyword EXPMIDPT = '%s' in '%s'!",
        buffer,filename);
      tmid = build_decsec_plus(hour,min,sec);
      expm = tmid - tstart;
      if (expm < -43200.0) expm += 86400.0;
      if (expm < 0.0) expm = 0.0;
    }
    write_keyword_double(fits,&fits->head,"MIDEXP",expm,
       "Flux-weighted mid-exposure time (seconds)");
  }

  seq = find_fits_keyword(&fits->head,"SPECGRPH");
  if (seq < 1) write_keyword_textual
      (fits,&fits->head,"SPECGRPH","KiwiSpec","Spectrograph name");

  seq = find_fits_keyword(&fits->head,"DETECTOR");
  if (seq < 1)
  {
    cam = get_keyword_integer(fits,&fits->head,"PARAM48");
    sprintf(buffer,"SI%d",cam);
    write_keyword_textual
      (fits,&fits->head,"DETECTOR",buffer,"Detector name");
  }

  seq = find_fits_keyword(&fits->head,"CCDMODE");
  if (seq < 1)
  {
    par26 = get_keyword_integer(fits,&fits->head,"PARAM26");
    par27 = get_keyword_integer(fits,&fits->head,"PARAM27");
    switch (par26 + par27)
    {
      case 90: mode = 0;
               break;
      case 92: mode = 1;
               break;
      case 93: mode = 2;
               break;
      case 40: mode = 3;
               break;
      case 42: mode = 4;
               break;
      case 43: mode = 5;
               break;
      case 10: mode = 6;
               break;
      case 12: mode = 7;
               break;
      case 13: mode = 8;
               break;
      case  2: mode = 9;
               break;
      default: get_error("Invalid CCD camera mode!");
    }
    write_keyword_integer
      (fits,&fits->head,"CCDMODE",mode,"CCD camera mode");
  }

  seq = find_fits_keyword(&fits->head,"FIBRESET");
  if (seq < 1) write_keyword_integer
      (fits,&fits->head,"FIBRESET",1,"Optical fibre setting");

  seq = find_fits_keyword(&fits->head,"OFFSET");
  if (seq < 1)
  {
    bias = get_keyword_integer(fits,&fits->head,"PARAM28");
    write_keyword_integer
      (fits,&fits->head,"OFFSET",bias,"CCD offset (bias) in ADU");
  }
}

void extract_image_proc(char *xname,char *yname,int ax,int ay,int bx,int by,
                        int info)
{
   byte *dest;
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
   if (q.rem != 0) get_error
     ("rotation_argument: Only integer multiples of 90 degrees are allowed!");

   return(q.quot);
}

int flip_argument(int id)
{
   if (id < 4) return(1-id); /* NONE = 0, X = -1, Y = -2 */
   return(2); /* Flip XY is the same as a rotation by 180 degrees */
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
   byte *xptr,*yptr,*xstart,*ystart;
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
         memmove(yptr,xptr,xfits->pixsize);
   }
}

void rotate_fits_image(char *xname,char *yname,int q)
{
   fits_file xfits,yfits;

   read_raw_image(&xfits,xname);
   check_image_2d(&xfits);

   yfits = xfits;
   set_fits_name(&yfits,yname);
   allocate_whole_binary(&yfits);

   switch (q)
   {
      case -1: rotate_pixels(&xfits,&yfits,xfits.npix[1],-1,1,1,0);
               break;
      case -2: rotate_pixels(&xfits,&yfits,xfits.npix[1],
                                                -1,xfits.totpix,-1,0);
               break;
      case 0:  rotate_pixels(&xfits,&yfits,1,1,1,1,0);
               break;
      case 1:  swap_axes(&yfits);
               rotate_pixels(&xfits,&yfits,1,1,xfits.npix[2],
                                       xfits.npix[2],-xfits.totpix-1);
               break;
      case 2:  rotate_pixels(&xfits,&yfits,1,1,xfits.totpix,-1,0);
               break;
      case 3:  swap_axes(&yfits);
               rotate_pixels(&xfits,&yfits,xfits.npix[1],-1,1,
                                       xfits.npix[2],-xfits.totpix+1);
               break;
      default: get_error("rotate_fits_image: "
                         "Bad rotation argument (%d)!",q);
   }

   write_raw_image(&yfits,yname);

   free_binary_data(&xfits);
}

void rotate_image_proc(char *xname,char *yname,int angle,int pmod)
{
   if (pmod) inform("Rotate:      %s.fit --> %s.fit",xname,yname);

   rotate_fits_image(xname,yname,rotation_argument(angle));
}

void flip_image_proc(char *xname,char *yname,int axis,int pmod)
{
   if (pmod) inform("Flip:        %s.fit --> %s.fit",xname,yname);

   rotate_fits_image(xname,yname,flip_argument(axis));
}

void quad_align_proc(char *xname,char *yname)
{
   fits_file xfits,yfits;
   int width,height,bitpix,wrow,hrow,quadw,quadh;
   int col,row;
   byte *aptr,*bptr,*uptr,*vptr;
   double apix,bpix,upix,vpix;

   logfile("quad_align %s %s",xname,yname);
   inform("Align:       %s.fit ---> %s.fit",xname,yname);

   read_image_header(&xfits,xname);
   check_image_2d(&xfits);

   width = xfits.npix[1];
   height = xfits.npix[2];
   bitpix = xfits.bitpix;

   wrow = find_fits_keyword(&xfits.head,"QUADWDTH");
   hrow = find_fits_keyword(&xfits.head,"QUADHGHT");

   if (wrow < 1 || hrow < 1)
   {
      get_system("cp %s.fit %s.fit",xname,yname);
      inform("Quadrant size not specified. File copy performed.");
      return;
   }

   quadw = get_keyword_integer(&xfits,&xfits.head,"QUADWDTH");
   quadh = get_keyword_integer(&xfits,&xfits.head,"QUADHGHT");

   if (quadw < 1 || quadw > width - 3)
      get_error("Bad QUADWDTH (%d)!",quadw);
   if (quadh < 1 || quadh > height - 3)
      get_error("Bad QUADHGHT (%d)!",quadh);

   yfits = xfits;
   set_fits_name(&yfits,yname);

   allocate_whole_binary(&xfits);
   allocate_whole_binary(&yfits);

   read_image_pixels(&xfits);

   copy_pixel_rows(&xfits,1,&yfits,1,quadh);
   copy_pixel_rows(&xfits,quadh+1,&yfits,quadh+3,height-quadh-2);

   aptr = yfits.bin + (quadh - 1) * yfits.rowsize;
   uptr = aptr + yfits.rowsize;
   vptr = uptr + yfits.rowsize;
   bptr = vptr + yfits.rowsize;
   for (col=1;col<=width;col++)
   {
      apix = read_pixel_double(aptr,bitpix);
      bpix = read_pixel_double(bptr,bitpix);
      upix = (2 * apix + bpix) / 3.0;
      vpix = (apix + 2 * bpix) / 3.0;
      write_pixel_double(uptr,bitpix,upix);
      write_pixel_double(vptr,bitpix,vpix);
      aptr += yfits.pixsize;
      uptr += yfits.pixsize;
      vptr += yfits.pixsize;
      bptr += yfits.pixsize;
   }

   copy_pixel_cols(&yfits,1,&xfits,1,quadw);
   copy_pixel_cols(&yfits,quadw+1,&xfits,quadw+3,width-quadw-2);

   memmove(yfits.bin,xfits.bin,height*xfits.rowsize);

   aptr = yfits.bin + (quadw - 1) * yfits.pixsize;
   uptr = aptr + yfits.pixsize;
   vptr = uptr + yfits.pixsize;
   bptr = vptr + yfits.pixsize;
   for (row=1;row<=height;row++)
   {
      apix = read_pixel_double(aptr,bitpix);
      bpix = read_pixel_double(bptr,bitpix);
      upix = (2 * apix + bpix) / 3.0;
      vpix = (apix + 2 * bpix) / 3.0;
      write_pixel_double(uptr,bitpix,upix);
      write_pixel_double(vptr,bitpix,vpix);
      aptr += yfits.rowsize;
      uptr += yfits.rowsize;
      vptr += yfits.rowsize;
      bptr += yfits.rowsize;
   }

   inform("2 pixel rows and 2 pixel columns inserted.");
   write_image(&yfits,yname);

   free_binary_data(&xfits);
}

void convert_image_proc(char *xname,char *yname)
{
   fits_file xfits,yfits;
   int pix;
   double xpix;
   float *dest;
   byte *src;

   logfile("convert_image %s %s",xname,yname);
   inform("Convert:     %s.fit ---> %s.fit",xname,yname);

   read_image_header(&xfits,xname);

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
      xpix = read_pixel_double(src,xfits.bitpix);
      *dest = (float)(xpix * xfits.bscale + xfits.bzero);
   }

   write_keyword_integer(&yfits,&yfits.head,"BSCALE",1,"Pixel value scale");
   write_keyword_integer(&yfits,&yfits.head,"BZERO",0,"Pixel value offset");

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

void fft_fits(fits_file *xfits,fits_file *rfits,fits_file *ifits)
{
   float *data;
   int32_t n[3];
   float *rval,*ival,*re,*im;
   int i,pixcount,ndat;

   check_real_image(xfits);
   check_image_2d(xfits);
   check_fft_image(xfits);

   check_real_image(rfits);
   check_image_2d(rfits);
   check_fft_image(rfits);

   check_real_image(ifits);
   check_image_2d(ifits);
   check_fft_image(ifits);

   check_same_image_size(xfits,rfits);
   check_same_image_size(xfits,ifits);

   n[1] = xfits->npix[2];
   n[2] = xfits->npix[1];
   pixcount = n[1] * n[2];
   ndat = 2 * pixcount;

   data = fvector(ndat);
   clear_fvector(data,ndat);

   rval = xfits->pix;
   re = &data[1];
   for (i=0;i<pixcount;i++,rval++,re+=2) *re = *rval;

   fourn(data,n,2,1,1);

   rval = rfits->pix;
   ival = ifits->pix;
   re = &data[1];
   im = &data[2];
   for (i=0;i<pixcount;i++,rval++,ival++,re+=2,im+=2)
   {
      *rval = *re;
      *ival = *im;
   }

   free_fvector(data);
}

void ifft_fits(fits_file *xfits,fits_file *rfits,fits_file *ifits)
{
   float *data;
   int32_t n[3];
   float *rval,*ival,*re,*im;
   int i,pixcount,ndat;
   float scale;

   check_real_image(xfits);
   check_image_2d(xfits);
   check_fft_image(xfits);

   check_real_image(rfits);
   check_image_2d(rfits);
   check_fft_image(rfits);

   check_real_image(ifits);
   check_image_2d(ifits);
   check_fft_image(ifits);

   check_same_image_size(xfits,rfits);
   check_same_image_size(xfits,ifits);

   n[1] = xfits->npix[2];
   n[2] = xfits->npix[1];
   pixcount = n[1] * n[2];
   ndat = 2 * pixcount;

   data = fvector(ndat);
   clear_fvector(data,ndat);

   rval = rfits->pix;
   ival = ifits->pix;
   re = &data[1];
   im = &data[2];
   for (i=0;i<pixcount;i++,rval++,ival++,re+=2,im+=2)
   {
      *re = *rval;
      *im = *ival;
   }

   fourn(data,n,2,-1,1);

   scale = 1.0 / (float)pixcount;

   rval = xfits->pix;
   re = &data[1];
   for (i=0;i<pixcount;i++,rval++,re+=2) *rval = *re * scale;

   free_fvector(data);
}

void complex_product_fits(fits_file *xr,fits_file *xi,
   fits_file *yr,fits_file *yi,fits_file *zr,fits_file *zi,int conj)
{
   int i;
   float *ar,*ai,*br,*bi,*cr,*ci;
   float arbr,aibi,arbi,aibr;

   check_real_image(xr);
   check_image_2d(xr);
   check_fft_image(xr);

   check_real_image(xi);
   check_image_2d(xi);
   check_fft_image(xi);

   check_same_image_size(xr,xi);

   check_real_image(yr);
   check_image_2d(yr);
   check_fft_image(yr);

   check_real_image(yi);
   check_image_2d(yi);
   check_fft_image(yi);

   check_same_image_size(yr,yi);

   check_real_image(zr);
   check_image_2d(zr);
   check_fft_image(zr);

   check_real_image(zi);
   check_image_2d(zi);
   check_fft_image(zi);

   check_same_image_size(zr,zi);

   check_same_image_size(xr,yr);
   check_same_image_size(xr,zr);

   ar = xr->pix;
   ai = xi->pix;
   br = yr->pix;
   bi = yi->pix;
   cr = zr->pix;
   ci = zi->pix;

   for (i=0;i<xr->totpix;i++,ar++,ai++,br++,bi++,cr++,ci++)
   {
      arbr = (*ar) * (*br);
      aibi = (*ai) * (*bi);
      arbi = (*ar) * (*bi);
      aibr = (*ai) * (*br);
      if (conj)
         *cr = arbr + aibi, *ci = aibr - arbi;
      else
         *cr = arbr - aibi, *ci = aibr + arbi;
   }
}

void correlate_2d_proc(char *xname,char *rname,char *cname)
{
   fits_file x,r;
   int axis,n,m,nm;
   float *sig,*ref,**res;
   double xstart[3],rstart[3];
   float maxpix;
   int maxcol,maxrow;

   read_image(&x,xname);
   check_real_image(&x);
   check_image_2d(&x);
   check_fft_image(&x);

   read_image(&r,rname);
   check_real_image(&r);
   check_image_2d(&r);
   check_fft_image(&r);

   check_same_image_size(&x,&r);
   check_same_image_step(&x,&r);

   normalize_pixel_sum(&x);
   normalize_pixel_sum(&r);

   n = x.npix[1];
   m = x.npix[2];
   nm = n*m;

   for (axis=1;axis<=2;axis++)
   {
      xstart[axis] = get_axis_start(&x,axis);
      rstart[axis] = get_axis_start(&r,axis);
   }

   sig = x.pix - 1;
   ref = r.pix - 1;

   res = fmatrix(m,n);

   ccf2d(sig,ref,n,m,res[1]);

   shift_fmatrix(res,n,m,n/2,m/2);

   memmove(x.bin,&res[1][1],nm*sizeof(float));

   set_fits_name(&x,cname);

   write_keyword_double(&x,&x.head,"CRPIX1",(double)(n/2+1),
                         "Reference Pixel");
   write_keyword_double(&x,&x.head,"CRPIX2",(double)(m/2+1),
                         "Reference Pixel");
   write_keyword_double(&x,&x.head,"CRVAL1",xstart[1]-rstart[1],
                         "Coordinate at reference pixel");
   write_keyword_double(&x,&x.head,"CRVAL2",xstart[2]-rstart[2],
                         "Coordinate at reference pixel");

   update_fits_minmax(&x);
   locate_fits_max(&x,&maxpix,&maxcol,&maxrow);
   write_keyword_integer(&x,&x.head,"MAXCOL",maxcol,
                         "Maximum pixel column number");
   write_keyword_integer(&x,&x.head,"MAXROW",maxrow,
                         "Maximum pixel row number");

   write_image(&x,cname);

   free_fmatrix(res);

   free_binary_data(&r);
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

int keydef_list_count(void)
{
   int n;

   for (n=0;*keydef_list[n].item != 0;n++);

   return(n);
}

void clear_keydef_list(void)
{
   int i;

   for (i=0;*keydef_list[i].item != 0;i++)
   {
      keydef_list[i].found = 0;
      memset(keydef_list[i].sval,0,MAXCFGLINE+1);
      keydef_list[i].ival = 0;
      keydef_list[i].dval = 0.0;
      keydef_list[i].lval = '\x00';
   }
}

void populate_values(keydef_type *k,char *val)
{
   int n;
   cfglinetype theval,str;

   k->found = 0;
   memset(k->sval,0,sizeof(cfglinetype));
   k->ival = 0;
   k->dval = 0.0;
   k->lval = '\x00';

   memset(theval,0,sizeof(cfglinetype));

   n = strlen(val);
   if (n == 0) return;

   if (*val == '\x27')
   {
      if (n > 2) memmove(theval,val+1,n-2);
   }
   else
   {
      memmove(theval,val,n);
   }
   k->found = (*theval != 0);
   if (k->found)
   {
      strcpy(k->sval,theval);
      sscanf(theval,"%s",str);
      k->ival = str_to_int_def(str,0);
      k->dval = str_to_dbl_def(str,0.0);
      k->lval = str_to_log(str);
   }
}

void read_keydef_cfg(fits_file *fits)
{
   char spectrograph[MAX_SPECNAME+1];
   file_type cfg;
   int row;
   keydef_type *k;
   fitskey_type fitskey;
   fitsval_type fitsval;

   clear_keydef_list();

   get_spectrograph(fits,&fits->head,spectrograph);
   open_keydef_cfg(&cfg,spectrograph);

   for (k=keydef_list;*k->item != 0;k++)
   {
      expect_next_string_value(&cfg,k->item,fitskey,MAX_FITS_KEYLEN);
      row = find_fits_keyword(&fits->head,fitskey);
      if (row > 0)
      {
         read_keyword_value(fits,&fits->head,row,fitsval);
         populate_values(k,fitsval);
      }
   }

   get_fclose(&cfg);
}

void get_keydef_item(int seq,keydef_type *key)
{
   int n;

   n = keydef_list_count();
   if (seq >= 0 && seq < n)
      *key = keydef_list[seq];
   else
      *key = keydef_list[n];
}

int get_fibre_number(fits_file *fits)
{
   return(get_keyword_integer(fits,&fits->head,"FIBRENO"));
}

int count_thorium_ref(char *spectrograph)
{
   file_type f;
   int n,k,km;
   char row[256],*cptr;

   get_system("ls -1 %s/spec/%s/img/Thorium??.fit > ls.out",
     hrsp_dir(),spectrograph);
   open_file(&f,"r","ls.out");
   n = km = 0;
   while (fgets(row,sizeof(row),f.dat) != NULL)
   {
     if (strlen(row) < 16)
       get_error("count_thorium_ref: Path name too short!");
     if (strlen(row) > 250)
       get_error("count_thorium_ref: Path name too long!");
     cptr = strrchr(row,(int)'/');
     if (cptr == NULL)
        get_error("count_thorium_ref: Unexpected path name!");
     if (strncmp(cptr+1,"Thorium",7) != 0)
        get_error("count_thorium_ref: Unexpected file name!");
     if (strncmp(cptr+10,".fit",4) != 0)
        get_error("count_thorium_ref: Unexpected file extension!");
     if (!digit_ok(cptr[8]) || !digit_ok(cptr[9]))
        get_error("count_thorium_ref: Sequence number expected!");
     cptr[10] = '\x00';
     k = atoi(cptr+8);
     if (k > km) km = k;
     n++;
   }
   if (km != n)
     get_error("count_thorium_ref: Images out of sequence!");
   get_fclose(&f);
   remove_file("ls.out");

   return(n);
}

void get_dispersion
  (spec_type *spec,regre_block *ureg,int m,double *lamcen,double *disp)
{
   double **p;
   int i;
   double mlnorm[66],ucen[66],mlnorm2[66];
   double ucen_mid,mlnorm_mid,eps=1.0e-6,ua,ub;

   p = dmatrix(1,2);

   for (i=1;i<=65;i++)
   {
      p[1][1] = (double)(i-1) / 64.0;
      p[1][2] = normalized_order(spec,m);
      mlnorm[i] = p[1][1];
      ucen[i] = polynomial_val(p,1,2,ureg->deg,ureg->a,ureg->ma);
   }

   spline(ucen,mlnorm,65,1e30,1e30,mlnorm2);

   ucen_mid = spec->usize / 2.0;

   splint(ucen,mlnorm,mlnorm2,65,ucen_mid,&mlnorm_mid);

   p[1][1] = mlnorm_mid - eps;
   ua = polynomial_val(p,1,2,ureg->deg,ureg->a,ureg->ma);
   p[1][1] = mlnorm_mid + eps;
   ub = polynomial_val(p,1,2,ureg->deg,ureg->a,ureg->ma);

   *lamcen = denormalized_mlambda(spec,mlnorm_mid) / (double)m;
   *disp = 2.0 * eps * spec->mlamrange / (double)m / (ub-ua);

   free_dmatrix(p);
}

void get_dispersion_list(spec_type *spec,double *lamcen,double *disp)
{
   fits_file tfits;
   regre_block ureg;
   int ord;

   read_fits_header(&tfits,"%s/spec/%s/tbl/thar",hrsp_dir(),spec->specname);
   read_regression_block(&tfits,&tfits.extend,"UREGRE",&ureg);
   for (ord=spec->firstord;ord<=spec->lastord;ord++)
   {
      get_dispersion(spec,&ureg,ord,&lamcen[ord],&disp[ord]);
   }
}

double true_dispersion(regre_block *xreg,double mlnorm,double mnorm)
{
  int i,j,k,n,m;
  double disp;

  n = xreg->deg[1];
  m = xreg->deg[2];

  for (j=0,k=1,disp=0.0;j<=m;j++)
  {
    for (i=0;i<=n;i++,k++)
    {
      if (i==0) continue;
      disp += xreg->a[k] * (double)i * intpow(mlnorm,i-1) * intpow(mnorm,j);
    }
  }

  return(disp);
}

/* ---------------------------------------------------------------------------
   Relativistic Doppler effect:

      lam/lam0 = D,

   where lam and lam0 are the observed and emitted wavelengths, and D is
   the Doppler factor:

      D = sqrt( (1+beta) / (1-beta) ),

   where beta is the velocity in terms of the speed of light:

      beta = v/c

   In order to find beta from D, one can write:

      beta = (D^2 -  1) / (D^2 + 1),

   or in a more convenient form:

      beta = U / (U + 2),

   where U is a new quantity:

      U = D^2 - 1

   On a logarithmic scale:

     log(lam) - log(lam0) = log(D)

   Another way of expressing the Doppler shift is:

     (lam - lam0) / lam0 = z,

   where z is the redshift:

      z = D - 1

   *Computational note*
   For small velocities, the Doppler factor D is close to 1 and should not
   be used directly in radial-velocity calculations, as it can lead to a
   signifficant loss of precision. The logarithm of D should be calculated
   instead, using the log1p function in C. Similarly, when calulating beta
   from log(D), one should first obtain the quantity U, using the expm1
   function in C. The same applies to the redshift z, which should be
   obtained from log(D) without any loss in precision, using the expm1
   function.
   ------------------------------------------------------------------------ */

double log_doppler(double vel)
{
  double beta;

  beta = vel / LIGHT_SPEED;
  return((log1p(beta) - log1p(-beta)) * 0.5);
}

double doppler_vel(double logdopp)
{
  double u,beta;

  u = expm1(2 * logdopp);
  beta = u / (u + 2.0);
  return(beta * LIGHT_SPEED);
}

double redshift(double vel)
{
  return(expm1(log_doppler(vel)));
}

int create_file_list(dir_name_type dirname, file_name_type filename)
{
   int m,maxlen;
   file_type list;

   get_system("\\rm -f find.out");
   get_system("find %s -maxdepth 1 -iname \"%s\" | "
      "rev | cut -d/ -f1 | rev > find.out",dirname,filename);

   get_system("\\rm -f sed.a sed.b");
   for (m=1;m<=12;m++)
   {
     get_system("echo s/%s/ %02d/g >> sed.a",short_month_name(m),m);
     get_system("echo s/ %02d/%s/g >> sed.b",m,short_month_name(m));
   }

   get_system("\\rm -f find.kiwispec");
   get_system_err(2,"grep -i _KiwiSpec_ find.out | sed -f sed.a | "
      "sort -t _ -k 3,4 | sed -f sed.b > find.kiwispec");

   get_system("\\rm -f sed.a sed.b");

   get_system("\\rm -f find.other");
   get_system_err(2,"grep -v -i _KiwiSpec_ find.out > find.other");

   get_system("\\rm -f find.all");
   get_system("cat find.other find.kiwispec > find.all");

   get_system("\\rm -f wc.out");
   get_system("wc -L find.all | cut -f 1 -d \" \" > wc.out");

   open_file(&list,"r","wc.out");
   if (fscanf(list.dat," %d",&maxlen) != 1)
      get_error("Unable to determine the longest filename!");
   get_fclose(&list);

   get_system("\\rm -f find.out find.kiwispec find.other wc.out");

   return(maxlen);
}
