/* ------------------------------------------------------------------------

   Program:  quad_align
   Purpose:  Insert two missing rows and two missing columns between
             image quadrants. This will fix the quadrant misalignment on
             some images taken with the CCD486 camera.

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
   if (argc != 3) get_error("Usage: quad_align <in> <out>");

   quad_align_proc(argv[1],argv[2]);

   check_memory();

   return(0);
}
