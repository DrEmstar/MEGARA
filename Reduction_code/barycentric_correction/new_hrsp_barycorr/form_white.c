/* ------------------------------------------------------------------------

   Program:  form_white
   Purpose:  Create a reference white lamp image from a set of input 
             images.

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
   int jd,img,count,seq,sel[MAX_IMGNO+1];
   char rname[9];
   file_name_type sum,newsum,mean;
   fits_file fits;
   char card[80];

   if (argc != 3) get_error("Usage:  form_white  <JD> <img>");

   logfile("form_white %s %s",argv[1],argv[2]);

   jd = collect_jd_number(argv[1]);
   if (strlen(argv[2]) > 68) get_error("The image list argument too long!");
   count = collect_integer_list(argv[2],sel,MAX_IMGNO,"Image number");
   if (count < 1) get_error("No images selected!");
   if (count == 1)
     inform("1 image selected.");
   else
     inform("%d images selected.",count);

   for (img=1,seq=0;img<=MAX_IMGNO;img++)
   {
     if (!sel[img]) continue;
     sprintf(rname,"r%04d%03d",jd,img);
     seq++;
     if (seq == 1)
     {
        next_tmp_file_name(sum);
        inform("\nCopy:   %s.fit --> %s.fit",rname,sum);
        get_system("cp %s.fit %s.fit",rname,sum);
     }
     else
     {
        next_tmp_file_name(newsum);
        inform("\nAdd:   %s + %s = %s",sum,rname,newsum);
        add_image_file(sum,rname,newsum,YES_INFO);
        remove_tmp_file("%s.fit",sum);
        strcpy(sum,newsum);
     }
   }

   next_tmp_file_name(mean);
   inform("\nMean:   %s / %d = %s",sum,count,mean);
   image_div_scalar_file(sum,mean,(double)count,YES_INFO);
   remove_tmp_file("%s.fit",sum);

   read_image(&fits,mean);

   write_keyword_integer
     (&fits,&fits.head,"WHITEJD",jd,"Incomplete JD for white image");

   sprintf(card,"WHITEIMG= '%s'",argv[2]);
   write_fits_card(&fits,&fits.head,card);

   update_fits_minmax(&fits);
   write_image(&fits,"rwhite");
   inform("\nImage '%s' created.",fits.file.path);

   remove_tmp_file("%s.fit",mean);

   check_memory();

   return(0);
}
