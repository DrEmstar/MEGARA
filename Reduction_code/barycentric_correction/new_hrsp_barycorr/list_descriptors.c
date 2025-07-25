/* ------------------------------------------------------------------------

   Program:  list_descriptors
   Purpose:  Generate a list of fundamental descriptors for every
             FITS image found in a given directory.

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

#define FIELDCOUNT 12

typedef struct
{
   int seq;
   int len;
   int type;
   char caption[80];
}
fieldtype;

fieldtype fields[FIELDCOUNT] =
{
   {-1, 12, TYPESTR, "File Name"},
   { 2, 10, TYPESTR, "Date"},
   { 3, 11, TYPESTR, "Time"},
   { 0, 12, TYPESTR, "Object"},
   { 1, 14, TYPESTR, "Exposure Type"},
   { 4,  6, TYPEDBL, "Exp"},
   { 5,  6, TYPEDBL, "Mid"},
   { 6,  6, TYPEDBL, "Bias"},
   { 7, 10, TYPESTR, "CCD"},
   { 8, 10, TYPESTR, "Spec"},
   { 9,  1, TYPEINT, "G"},
   {10,  1, TYPEINT, "F"}
};

void print_dashed_line(void)
{
   int i,j;

   for (i=0;i<FIELDCOUNT;i++)
   {
      if (i > 0) printf("-");
      for (j=0;j<fields[i].len;j++) printf("-");
   }
   printf("\n");
}

char *sfmt(int i,char *fmt)
{
   if (fields[i].type == TYPESTR)
      sprintf(fmt,"%%%ds",-fields[i].len);
   else
      sprintf(fmt,"%%%ds",fields[i].len);
   return(fmt);
}

void print_captions(void)
{
   int i;
   char fmt[16];

   for (i=0;i<FIELDCOUNT;i++)
   {
      if (i > 0) printf(" ");
      printf(sfmt(i,fmt),fields[i].caption);
   }
   printf("\n");
}

void print_values(char *fname)
{
   int i;
   char fmt[16];
   keydef_type k;
   char val[80],object[80];
   starlabel_type lab;

   strcpy(k.sval,fname);
   k.found = 1;
   for (i=0;i<FIELDCOUNT;i++)
   {
      if (i > 0)
      {
         printf(" ");
         get_keydef_item(fields[i].seq,&k);
      }
      memset(val,0,80);
      if (k.found)
      {
         switch(fields[i].type)
         {
            case TYPESTR: strcpy(val,k.sval);
                          break;
            case TYPEINT: sprintf(val,"%d",k.ival);
                          break;
            case TYPEDBL: sprintf(val,"%3.1f",k.dval);
                          break;
         }
      }
      else
      {
         memset(val,'-',fields[i].len);
      }

      if (i == 3) strcpy(object,val);
      if (i == 4 && strncasecmp(val,"stellar",7) == 0)
        if (parse_star_name(object,&lab) != 0) strcpy(val,"********");

      if (strlen(val) > fields[i].len)
      {
         val[fields[i].len] = 0;
         val[fields[i].len-1] = '*';
      }
      printf(sfmt(i,fmt),val);
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
      get_error("Usage:  list_descriptors  <filename>");

   logfile("list_descriptors %s",argv[1]);

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
      if (strlen(row) > MAX_FILE_NAME) get_error("File name too long!");
      strcpy(filename,row);
      sprintf(pathname,"%s/%s",dirname,filename);
      read_image_header(&xfits,pathname);
      switch (check_file_type(filename))
      {
        case KIWISPEC_TYPE: fix_kiwispec_header(&xfits,filename);
                            break;
      }
      read_keydef_cfg(&xfits);
      print_values(xfits.file.name);
   }

   get_fclose(&list);

   get_system("\\rm -f find.all");

   check_memory();

   return(0);
}
