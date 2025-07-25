/* ------------------------------------------------------------------------

   Program:  reduce_white
   Purpose:  Reduce a white lamp image (or a set of images)

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
   if (argc != 3) get_error("Usage:  reduce_white  <JD> <img>");

   logfile("reduce_white %s %s",argv[1],argv[2]);

   inform("");
   inform("Reduction of white lamp image");
   inform("=============================");

   get_system("form_white %s %s",argv[1],argv[2]);
   inform("");
   get_system("define_echelle rwhite");
   inform("");
   get_system("filter_white rwhite cwhite");
   inform("");
   get_system("extract_echelle cwhite ewhite");
   inform("");
   get_system("form_rflat cwhite");
   inform("");
   get_system("form_eflat ewhite");

   check_memory();

   return(0);
}
