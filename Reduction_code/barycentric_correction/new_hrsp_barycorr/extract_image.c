/* ------------------------------------------------------------------------

   Program:  extract_image
   Purpose:  Extract a rectangular sub-area from a two-dimensional image

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
   if (argc != 7)
      get_error("Usage:  extract_image  <in> <out> <ax> <ay> <bx> <by>");

   extract_image_proc(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]),
                                      atoi(argv[5]),atoi(argv[6]),YES_INFO);

   check_memory();

   return(0);
}
