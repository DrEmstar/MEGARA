/* ------------------------------------------------------------------------

   Program:  prepare_echelle
   Purpose:  Prepare raw echelle spectra for data reduction

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
   int width;
   int height;
   int left;
   int bottom;
}
geometry;

void get_geometry(char *a,geometry *g)
{
   char *s,*b;

   memset(g,0,sizeof(geometry));

   if (strcasecmp(a,"F") == 0) return;

   s = a;

   b = s;
   if (!digit_ok(*s)) get_error("Bad geometry '%s'!",a);
   while(digit_ok(*s)) s++;
   if (*s != 'x') get_error("Bad geometry '%s'!",a);
   g->width = atoi(b);
   if (g->width < 128) get_error("Bad geometry '%s'!",a);
   s++;

   b = s;
   if (!digit_ok(*b)) get_error("Bad geometry '%s'!",a);
   while(digit_ok(*s)) s++;
   if (*s != '+') get_error("Bad geometry '%s'!",a);
   g->height = atoi(b);
   if (g->height < 128) get_error("Bad geometry '%s'!",a);
   s++;

   b = s;
   if (!digit_ok(*b)) get_error("Bad geometry '%s'!",a);
   while(digit_ok(*s)) s++;
   if (*s != '+') get_error("Bad geometry '%s'!",a);
   g->left = atoi(b);
   s++;

   b = s;
   if (!digit_ok(*b)) get_error("Bad geometry '%s'!",a);
   while(digit_ok(*s)) s++;
   if (*s != '\x00') get_error("Bad geometry '%s'!",a);
   g->bottom = atoi(b);
}

int main(int argc,char **argv)
{
   file_type list;
   fits_file yfits;
   char row[MAX_PATH_NAME+8];
   path_name_type pathname;
   dir_name_type dirname;
   file_name_type filename;
   char rname[9];
   int imgtype,mjd,sel[MAX_IMGNO+1],imgno,count,nprep;
   file_name_type acpy,aspl,aext,arot,aflp;
   ccd_type ccd;
   geometry crop;
   int acol,arow,bcol,brow,full_frame;

   if (argc != 4)
      get_error("Usage:  prepare_echelle <crop> <filename> <img>");

   logfile("prepare_echelle %s %s %s",argv[1],argv[2],argv[3]);

   inform_title("Image preparation");

   if (strlen(argv[1]) > 19) get_error("Crop argument too long!");
   get_geometry(argv[1],&crop);
   acol = crop.left + 1;
   arow = crop.bottom + 1;
   bcol = acol + crop.width - 1;
   brow = arow + crop.height - 1;
   full_frame = (crop.width == 0);

   if (strlen(argv[2]) > MAX_PATH_NAME) get_error("Path name too long!");
   strcpy(pathname,argv[2]);
   split_path(pathname,dirname,filename);

   get_system("\\rm -f find.out");
   get_system("find %s -maxdepth 1 -iname \"%s\" > find.out",dirname,filename);

   count = collect_integer_list(argv[3],sel,MAX_IMGNO,"Image number");
   if (count < 1) get_error("No images selected!");

   open_file(&list,"r","find.out");
   nprep = 0;
   while (fgets(row,sizeof(row),list.dat) != NULL)
   {
      trim_right(row);
      if (strlen(row) > MAX_PATH_NAME) get_error("File name too long!");
      strcpy(pathname,row);
      split_path(pathname,dirname,filename);
      mjd = imgno = 0;
      switch (imgtype = check_file_type(filename))
      {
        case STANDARD_TYPE: parse_standard_name(filename,&mjd,&imgno);
                            break;
        case KIWISPEC_TYPE: parse_kiwispec_name(filename,&mjd,&imgno);
                            break;
        default: get_error("Unsupported file name: '%s'!",filename);
      }

      if (!sel[imgno]) continue;

      sprintf(rname,"r%04d%03d",mjd,imgno);

      next_tmp_file_name(acpy);
      inform("\nCopy FITS:   %s --> %s.fit",pathname,acpy);
      get_system("cp %s %s.fit",pathname,acpy);
      get_system("chmod u+w,a-x %s.fit",acpy);

      inform("Checking the required FITS keywords: %s.fit",acpy);
      read_raw_image(&yfits,acpy);
      switch (imgtype)
      {
        case KIWISPEC_TYPE: fix_kiwispec_header(&yfits,filename);
                            break;
      }
      write_raw_image(&yfits,acpy);

      get_system("modify_keywords %s",acpy);
      next_tmp_file_name(aspl);
      get_system("quad_align %s %s",acpy,aspl);
      remove_tmp_file("%s.fit",acpy);

      write_image_descriptor_integer(aspl,"DATAMIN",0,
                                     "Data minimum not defined");
      write_image_descriptor_integer(aspl,"DATAMAX",0,
                                     "Data maximum not defined");

      get_system("hrsp_descriptors %s",aspl);

      get_ccdinfo_file(aspl,&ccd);

      next_tmp_file_name(aext);
      if (full_frame)
      {
        inform("Extract:     %s.fit --> %s.fit",aspl,aext);
        get_system("cp %s.fit %s.fit",aspl,aext);
      }
      else
        extract_image_proc(aspl,aext,acol,arow,bcol,brow,YES_INFO);
      write_image_descriptor_textual
        (aext,"CROPARG",argv[1],"Crop argument");
      remove_tmp_file("%s.fit",aspl);
      next_tmp_file_name(arot);
      rotate_image_proc(aext,arot,ccd.rotate,YES_INFO);
      remove_tmp_file("%s.fit",aext);
      next_tmp_file_name(aflp);
      flip_image_proc(arot,aflp,ccd.flip.id,YES_INFO);
      remove_tmp_file("%s.fit",arot);
      convert_image_proc(aflp,rname);
      remove_tmp_file("%s.fit",aflp);
      subtract_bias(rname);
      update_image_minmax(rname);
      nprep++;
   }
   get_fclose(&list);
   get_system("\\rm -f find.out");

   inform("\n%d image(s) prepared.",nprep);

   check_memory();

   return(0);
}
