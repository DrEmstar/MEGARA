/* ------------------------------------------------------------------------

   Program:  show_flux
   Purpose:  List the spectral intensity values in a given spectral order.

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

#define MAX_ORD 1000

int main(int argc,char **argv)
{
  fits_file fits;
  int first_bin,last_bin;
  int numpix,numord,firstord,lastord,ord,seq,i;
  double *startval,*stepval;
  double x,y;

  if (argc != 5) get_error("Usage:  show_flux <img> <start> <end> <ord>");

  read_image(&fits,argv[1]);
  check_image_2d(&fits);

  numpix = fits.npix[1];
  numord = fits.npix[2];

  first_bin = atoi(argv[2]);
  last_bin = atoi(argv[3]);

  if (first_bin < 1 || first_bin > numpix)
    get_error("Start pixel out of range!");

  if (last_bin < 1 || last_bin > numpix)
    get_error("End pixel out of range!");

  firstord = nint(world_coordinate(&fits,2,(double)numord));
  lastord = nint(world_coordinate(&fits,2,1.0));
  if (firstord < 1 || firstord > MAX_ORD)
    get_error("First order number out of limits in '%s'!",fits.file.path);
  if (lastord < 1 || lastord > MAX_ORD)
    get_error("Last order number out of limits in '%s'!",fits.file.path);
  if (lastord-firstord != numord-1)
    get_error("Unexpected order numbers in '%s'!",fits.file.path);

  startval = dvector(numord);
  stepval = dvector(numord);

  get_double_array(&fits,&fits.head,"STARTVAL",&startval[1],numord);
  get_double_array(&fits,&fits.head,"STEPVAL",&stepval[1],numord);

  ord = atoi(argv[4]);
  if (ord < firstord || ord > lastord)
    get_error("Order number %d not found in '%s'!",ord,fits.file.path);

  seq = pixel_coordinate(&fits,2,(double)ord);

  for (i=1;i<=numpix;i++)
  {
    if (i < first_bin || i > last_bin) continue;
    x = startval[seq]+(i-1)*stepval[seq];
    y = (double)collect_pixel_value(&fits,i,seq);
    printf("%4d %17.10e %17.10e\n",i,x,y);
  }

  free_binary_data(&fits);

  free_dvector(startval);
  free_dvector(stepval);

  check_memory();

  return(0);
}
