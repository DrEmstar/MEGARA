/* ------------------------------------------------------------------------

   Program:  identify_object
   Purpose:  Identify the object observed in a given image.

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
   fits_file fits;
   int hipnum,usrnum;
   char objtype[MAX_OBJTYPE+1];
   char starname[MAX_STARNAME+1],skyname[MAX_SKYNAME+1],starstr[256];
   starlabel_type lab;

   if (argc != 2) get_error("Usage:  identify_object  <img>");

   inform_title("Object identification");

   load_astro();
   load_cats();

   read_raw_image(&fits,argv[1]);

   get_keyword_textual(&fits,&fits.head,"STARNAME",starname,MAX_STARNAME);
   get_keyword_textual(&fits,&fits.head,"SKYNAME",skyname,MAX_SKYNAME);

   inform("Object name: %s",starname);

   hipnum = usrnum = 0;

   if (strcasecmp(starname,skyname) == 0)
      strcpy(objtype,"SKY");
   else
   {
      strcpy(objtype,"STAR");
      parse_star_name_err(starname,&lab);
   }

   write_keyword_textual(&fits,&fits.head,"OBJTYPE",objtype,"Object type");

   inform("Object type: %s",objtype);

   if (strcmp(objtype,"STAR") == 0)
   {
      print_star_label(&lab,starstr);
      inform("Star recognized: %s",starstr);
      if (lab.kind == USR_STAR)
      {
         usrnum = lab.catnum;
      }
      else
      {
        any_to_hip(&lab);
        hipnum = lab.catnum;
        print_star_label(&lab,starstr);
        inform("Hipparcos catalog number: %s",starstr);
      }

      write_keyword_integer(&fits,&fits.head,"USRNUM",usrnum,
        "User catalogue number");
      write_keyword_integer(&fits,&fits.head,"HIPNUM",hipnum,
        "Hipparcos catalogue number");
      write_keyword_textual(&fits,&fits.head,"HIPCOMP",lab.comp,
        "Hipparcos multiple/double component");
   }

   write_raw_image(&fits,argv[1]);

   free_astro();
   free_cats();

   check_memory();

   return(0);
}
