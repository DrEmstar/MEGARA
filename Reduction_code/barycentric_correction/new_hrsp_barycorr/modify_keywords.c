/* ------------------------------------------------------------------------

   Program:  modify_keywords
   Purpose:  Modify the contents of the FITS header in a given FITS file,
             as specified in the 'fitskey.cfg' configuration file.

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

char cfg[] = "fitskey.cfg";
fits_file xfits;
fitsval_type filename;

void read_fitskey_cfg(void)
{
   char **keys,**vals;
   char card[CARD_BYTES+1];
   int n,i;

   keys = allocate_cfglines();
   vals = allocate_cfglines();

   read_matching_cfg(cfg,filename,keys,vals,&n);

   for (i=0;i<n;i++)
   {
      if (!fits_keyword_name_ok(keys[i]))
         get_error("Bad FITS keyword name '%s' in '%s'!",keys[i],cfg);
      if (strlen(vals[i]) > MAX_FITS_VALLEN) get_error
         ("FITS keyword value '%s' too long in '%s'!",vals[i],cfg);
      if (strlen(vals[i]) <= 20 && *vals[i] != '\x27')
         sprintf(card,"%-8s= %20s%50s",keys[i],vals[i]," ");
      else
         sprintf(card,"%-8s= %-70s",keys[i],vals[i]);
      write_fits_card(&xfits,&xfits.head,card);
      logfile("FITS keyword '%s' written.",keys[i]);
   }

   logfile("%d FITS keywords written.",n);

   free_cfglines(keys);
   free_cfglines(vals);
}

int main(int argc,char **argv)
{
   if (argc != 2) get_error("Usage: modify_keywords <img>");

   logfile("modify_keywords %s",argv[1]);

   read_image_header(&xfits,argv[1]);

   if (!file_exists(cfg))
   {
      inform("File '%s' not found. No changes to the FITS header made.",cfg);
      exit(0);
   }

   read_raw_image(&xfits,argv[1]);
   ensure_not_original_file(&xfits);
   get_keyword_filename(&xfits,filename);
   read_fitskey_cfg();
   write_raw_image(&xfits,argv[1]);

   check_memory();

   return(0);
}
