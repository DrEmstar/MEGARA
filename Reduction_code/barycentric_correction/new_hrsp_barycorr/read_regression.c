/* ------------------------------------------------------------------------

   Program:  read_regression
   Purpose:  Read a regression block from a FITS header.

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
   fits_file fits;
   fits_head *head;
   regre_block reg;
   int i;

   if (argc != 3) get_error("Usage:  read_regression  <file> <key>");

   head = read_fits_header(&fits,argv[1]);

   read_regression_block(&fits,head,argv[2],&reg);

   printf("%d\n",reg.ndim);
   printf("%d\n",reg.ycol);
   for (i=1;i<=reg.ndim;i++)
   {
      if (i > 1) printf(" ");
      printf("%d",reg.xcol[i]);
   }
   printf("\n");
   for (i=1;i<=reg.ndim;i++)
   {
      if (i > 1) printf(" ");
      printf("%d",reg.deg[i]);
   }
   printf("\n");
   printf("%d\n",reg.scol);
   printf("%d\n",reg.fcol);
   printf("%d\n",reg.rcol);
   printf("%d\n",reg.nsel);
   printf("%22.15e\n",reg.rms);
   printf("%22.15e\n",reg.kappa);
   printf("%d\n",reg.ma);
   for (i=1;i<=reg.ma;i++) printf("%22.15e\n",reg.a[i]);

   check_memory();

   return(0);
}
