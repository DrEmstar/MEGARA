/* ------------------------------------------------------------------------

   Program:  list_keywords
   Purpose:  Generate a list of basic FITS keywords for every FITS
             image found in a given directory.

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

typedef struct
{
   fitskey_type key;
   int len;
   int flush;
}
fieldtype;

fieldtype fields[] =
{
   {"FILE", 12, -1},
   {"BITPIX", 6, 1},
   {"NAXIS", 5, 1},
   {"NAXIS1", 6, 1},
   {"NAXIS2", 6, 1},
   {"BSCALE", 8, 1},
   {"BZERO", 8, 1},
   {"", 0, 0}
};

#define SED_A "s/Jan/ 01/g\ns/Apr/ 04/g\ns/May/ 05/g"
#define SED_B "s/ 01/Jan/g\ns/ 04/Apr/g\ns/ 05/May/g"

void print_dashed_line(void)
{
   int i,j;

   for (i=0;*fields[i].key != 0;i++)
   {
      if (i > 0) printf("-");
      for (j=0;j<fields[i].len;j++) printf("-");
   }
   printf("\n");
}

void print_captions(void)
{
   int i;

   printf("%*s",fields[0].flush*fields[0].len,"File Name");
   for (i=1;*fields[i].key != 0;i++)
      printf(" %*s",fields[i].flush*fields[i].len,fields[i].key);
   printf("\n");
}

void print_values(fits_file *xfits)
{
   int i,row;
   fitsval_type val;

   for (i=0;*fields[i].key != 0;i++)
   {
      memset(val,0,sizeof(fitsval_type));
      if (i > 0)
      {
         row = find_fits_keyword(&xfits->head,fields[i].key);
         if (row > 0)
            read_keyword_value(xfits,&xfits->head,row,val);
         else
            memset(val,'-',fields[i].len);
      }
      else
        strncpy(val,xfits->file.name,sizeof(fitsval_type)-1);
      if (strlen(val) > fields[i].len)
      {
         val[fields[i].len] = 0;
         val[fields[i].len-1] = '*';
      }
      if (i > 0) printf(" ");
      printf("%*s",fields[i].flush*fields[i].len,val);
   }
   printf("\n");
}

int main(int argc,char **argv)
{
   file_type list;
   fits_file xfits;
   char row[512];
   path_name_type pathname;
   dir_name_type dirname;
   file_name_type filename;

   if (argc != 2)
      get_error("Usage:  list_keywords  <filename>");

   logfile("list_keywords %s",argv[1]);

   if (strlen(argv[1]) > MAX_PATH_NAME) get_error("Path name too long!");
   strcpy(pathname,argv[1]);
   split_path(pathname,dirname,filename);

   fields[0].len = create_file_list(dirname,filename);

   print_dashed_line();
   print_captions();
   print_dashed_line();

   open_file(&list,"r","find.all");

   while (fgets(row,sizeof(row),list.dat) != NULL)
   {
      trim_right(row);
      if (strlen(row) > MAX_PATH_NAME) get_error("File name too long!");
      strcpy(filename,row);
      sprintf(pathname,"%s/%s",dirname,filename);
      read_image_header(&xfits,pathname);
      print_values(&xfits);
   }

   get_fclose(&list);

   get_system("\\rm -f find.all");

   check_memory();

   return(0);
}
