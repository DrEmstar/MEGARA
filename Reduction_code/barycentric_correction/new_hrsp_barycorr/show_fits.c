/* ------------------------------------------------------------------------

   Program:  show_fits
   Purpose:  Generate a list of keywords found in the FITS header of a 
             given FITS file

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
   file_type f;
   char record[RECORD_BYTES];
   char text[CARD_BYTES+1],*card,*recstop;
   int head,i;

   if (argc != 2) get_error("Usage:  show_fits  <file>");

   open_file(&f,"r","%s",argv[1]);

   head=0;
   recstop = record+RECORD_BYTES;
   while (fread(record,1,RECORD_BYTES,f.dat) == RECORD_BYTES)
   {
      if (head == 0)
      {
         if (same_text(record,"SIMPLE  = ")) head = 1;
         if (same_text(record,"XTENSION= ")) head = 2;
         if (head == 2) printf("\n");
      }
      if (head > 0)
      {
         for (card=record;card<recstop && head>0;card+=CARD_BYTES)
         {
            memmove(text,card,CARD_BYTES);
            text[CARD_BYTES] = 0;
            for (i=0;i<CARD_BYTES;i++)
              if (text[i]<'\x20' || text[i]>'\x7E') text[i] = '?';
            printf("%s\n",text);
            if (same_text(card,"END     ")) head=0;
         }
      }
   }

   get_fclose(&f);

   check_memory();

   return(0);
}
