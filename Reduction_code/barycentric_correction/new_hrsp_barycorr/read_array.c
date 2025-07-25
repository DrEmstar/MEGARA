/* ------------------------------------------------------------------------

   Program:  read_array
   Purpose:  Read a double-precision array from a FITS header.

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
   double *a;
   int i,n;

   if (argc != 4) get_error("Usage:  read_array  <file> <key> <n>");

   n = atoi(argv[3]);

   a = dvector(n);

   head = read_fits_header(&fits,argv[1]);

   get_double_array(&fits,head,argv[2],&a[1],n);

   for (i=1;i<=n;i++) printf("%22.15e\n",a[i]);

   free_dvector(a);

   check_memory();

   return(0);
}
