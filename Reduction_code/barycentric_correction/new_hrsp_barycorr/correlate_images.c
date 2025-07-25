/* ------------------------------------------------------------------------

   Program:  correlate_images
   Purpose:  Cross-correlate two images (two-dimensional FFT).

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
   fits_file cfits;

   if (argc != 4) get_error("Usage:  correlate_images  <img> <ref> <ccf>");

   inform("Image correlation.");

   correlate_2d_proc(argv[1],argv[2],argv[3]);

   read_image_header(&cfits,argv[3]);
   inform("Image '%s' created.",cfits.file.path);

   check_memory();

   return(0);
}
