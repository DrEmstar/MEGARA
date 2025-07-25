/* ------------------------------------------------------------------------

   Program:  unpack_echelle
   Purpose:  Create a set of one-dimensional images for every order.

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
   fits_file xfits,yfits;
   double *startval,*stepval;
   file_name_type tmp;
   int row,order,firstord,lastord;
   unsigned char *xptr;

   if (argc != 2) get_error("Usage:  unpack_echelle  <img>");

   setbuf(stderr,NULL);

   inform_title("Unpacking echelle orders");

   read_image(&xfits,argv[1]);
   check_image_2d(&xfits);

   startval=dvector(xfits.npix[2]);
   stepval=dvector(xfits.npix[2]);

   get_double_array(&xfits,&xfits.head,"STARTVAL",&startval[1],xfits.npix[2]);
   get_double_array(&xfits,&xfits.head,"STEPVAL",&stepval[1],xfits.npix[2]);

   create_image(next_tmp_file_name(tmp),&xfits.head,xfits.bitpix,
                1,xfits.npix[1],1.0,1.0,1.0);

   inform("%d rows to process.",xfits.npix[2]);

   firstord = lastord = (int)world_coordinate(&xfits,2,1.0);

   fprintf(stderr,"Processing row number:   0");
   for (row=1;row<=xfits.npix[2];row++)
   {
      fprintf(stderr,"\b\b\b%3d",row);
      order = (int)world_coordinate(&xfits,2,(double)row);
      if (order < firstord) firstord = order;
      if (order > lastord) lastord = order;
      xptr = xfits.bin + (row-1)*xfits.npix[1]*xfits.pixsize;
      read_image(&yfits,tmp);
      memmove(yfits.bin,xptr,xfits.npix[1]*xfits.pixsize);
      set_fits_name(&xfits,"%s%03d",argv[1],order);
      write_keyword_double(&yfits,&yfits.head,"CRVAL1",startval[row],
                                            "Coordinate at reference pixel");
      write_keyword_double(&yfits,&yfits.head,"CDELT1",stepval[row],
                                            "Coordinate increment per pixel");
      write_keyword_integer(&yfits,&yfits.head,"ORDER",order,
                                            "Absolute order number");
      write_image(&yfits,"%s%03d",argv[1],order);
   }
   fprintf(stderr,"\n");

   remove_tmp_file("%s.fit",tmp);

   free_dvector(startval);
   free_dvector(stepval);

   inform("%d images '%s%03d' ... '%s%03d' created.",
      xfits.npix[2],argv[1],firstord,argv[1],lastord);

   free_binary_data(&xfits);

   check_memory();

   return(0);
}
