/* ------------------------------------------------------------------------

   Program:  copy_keyword
   Purpose:  Copy a FITS keyword from one image to another.

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
   if (argc != 4)
     get_error("Usage:  copy_keyword  <src> <dest> <key>");

   copy_keyword_file(argv[1],argv[2],argv[3],YES_INFO);

   check_memory();

   return(0);
}
