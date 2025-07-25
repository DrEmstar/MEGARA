/* ------------------------------------------------------------------------

   Program:  dispersion_data
   Purpose:  Calculate the dispersions for each echelle order.

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
   fits_file xfits;
   regre_block xreg;
   double *wcen,*dcen,*gcen;
   int numord,lastord,row,ord;
   double ra,rb,rc,m,p,pcen;
   double **u;

   if (argc != 2) get_error("Usage:  dispersion_data  <img>");

   read_raw_image(&xfits,argv[1]);

   if (xfits.cdelt[2] != -1.0)
     get_error("Decreasing order numbers expected in '%s'!",xfits.file.path);

   get_specinfo_fits(&xfits,&xfits.head,&spec);

   numord = xfits.npix[2];
   lastord = (int)xfits.crval[2];

   pcen = (xfits.npix[1] + 1) * 0.5;

   wcen = dvector(numord);
   dcen = dvector(numord);
   gcen = dvector(numord);

   u = dmatrix(1,2);

   read_regression_block(&xfits,&xfits.head,"XREGRE",&xreg);

   for (row=1,ord=lastord;row<=numord;row++,ord--)
   {
     u[1][2] = m = normalized_order(&spec,ord);
     ra = -1.5;
     rb = 1.5;
     while (rb - ra > 1e-12)
     {
       u[1][1] = rc = (ra + rb) * 0.5;
       p = polynomial_val(u,1,2,xreg.deg,xreg.a,xreg.ma);
       if (p < pcen)
         ra = u[1][1];
       else
         rb = u[1][1];
     }
     rc = (ra + rb) * 0.5;
     wcen[row] = denormalized_mlambda(&spec,rc) / (double)ord;
     dcen[row] = true_dispersion(&xreg,rc,m) * (double)ord * spec.mlamscale;
     gcen[row] = 1.0 / (dcen[row] * wcen[row]);
   }

   write_double_array(&xfits,&xfits.head,"WAVECENT",&wcen[1],numord);
   write_double_array(&xfits,&xfits.head,"DISPCENT",&dcen[1],numord);
   write_double_array(&xfits,&xfits.head,"LOGDISP",&gcen[1],numord);

   write_raw_image(&xfits,argv[1]);

   free_dvector(wcen);
   free_dvector(dcen);
   free_dvector(gcen);

   free_dmatrix(u);

   check_memory();

   return(0);
}
