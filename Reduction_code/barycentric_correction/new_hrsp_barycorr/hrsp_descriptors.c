/* ------------------------------------------------------------------------

   Program:  hrsp_descriptors
   Purpose:  Prepare basic HRSP descriptors using either the 
             corresponding FITS keywords from the raw image 
             (the keyword names are defined in 'keydef.cfg'), or the new 
             values specified in the file 'descriptors.cfg'.

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

char block_head[] =
   "-------------------------- HRSP DESCRIPTORS ---------------------------";

char block_bar[] =
   "-----------------------------------------------------------------------";

cfglinetype starname,source,expdate,expstart,ccdname,specname,skyname;
int gainset,fibreno;
double expsec,expmid,adu_bias;

keydef_type descriptor_list[] =
{
   /*  0 */ {"STARNAME", TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  1 */ {"SOURCE",   TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  2 */ {"EXPDATE",  TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  3 */ {"EXPSTART", TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  4 */ {"EXPSEC",   TYPEDBL, 0, "", 0, 0.0, '\x00'},
   /*  5 */ {"EXPMID",   TYPEDBL, 0, "", 0, 0.0, '\x00'},
   /*  6 */ {"ADU_BIAS", TYPEDBL, 0, "", 0, 0.0, '\x00'},
   /*  7 */ {"CCDNAME",  TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  8 */ {"SPECNAME", TYPESTR, 0, "", 0, 0.0, '\x00'},
   /*  9 */ {"GAINSET",  TYPEINT, 0, "", 0, 0.0, '\x00'},
   /* 10 */ {"FIBRENO",  TYPEINT, 0, "", 0, 0.0, '\x00'},
   /* 11 */ {"SKYNAME",  TYPESTR, 0, "", 0, 0.0, '\x00'},
   /* 12 */ {"",         0,       0, "", 0, 0.0, '\x00'}
};
int sky_descriptor=11;
int void_descriptor=12;

fits_file xfits;
fitsval_type filename;

datetime_type utc;
char midtime[MAX_TIME+1];
int izero=0;
char false='F';

typedef struct
{
   char name[9];
   int type;
   void *val;
   char comm[50];
}
descriptor_type;

descriptor_type descriptor_block[] =
  {{"CCDNAME",  TYPESTR, ccdname,   "CCD camera model name"},
   {"GAINSET",  TYPEINT, &gainset,  "CCD camera gain setting"},
   {"SPECNAME", TYPESTR, specname,  "Spectrograph name"},
   {"FIBRENO",  TYPEINT, &fibreno,  "Fibre number"},
   {"STARNAME", TYPESTR, starname,  "Star name"},
   {"SOURCE",   TYPESTR, source,    "Source name"},
   {"SKYNAME",  TYPESTR, skyname,   "Name used for sky spectra"},
   {"EXPDATE",  TYPESTR, expdate,   "Exposure date"},
   {"EXPSTART", TYPESTR, expstart,  "Exposure start time (UTC)"},
   {"EXPSEC",   TYPEDBL, &expsec,   "Exposure duration (seconds)"},
   {"EXPMID",   TYPEDBL, &expmid,   "Mid-exposure (seconds from EXPSTART)"},
   {"MIDTIME",  TYPESTR, midtime,   "Universal time at mid-exposure"},
   {"UT_MID",   TYPEDBL, &utc.time.frac, "UTC decimal hours at mid-exposure"},
   {"JD_MID",   TYPEDBL, &utc.jd,   "Julian Day at mid-exposure"},
   {"ADU_BIAS", TYPEDBL, &adu_bias, "Constant bias (ADU)"},
   {"BIASDONE", TYPELOG, &false,    "Bias has not been subtracted"},
   {"CROPARG",  TYPESTR, "",        "Crop argument"},
   {"VOFFSET",  TYPEINT, &izero,    "Vertical offset"},
   {"OBJTYPE",  TYPESTR, "",        "Object type"},
   {"USRNUM",   TYPEINT, &izero,    "User catalogue number"},
   {"HIPNUM",   TYPEINT, &izero,    "Hipparcos catalogue number"},
   {"HIPCOMP",  TYPESTR, "",        "Hipparcos multiple/double component"},
   {"RA_OBS",   TYPESTR, "",        "Apparent right ascension (hour,min,sec)"},
   {"DEC_OBS",  TYPESTR, "",        "Apparent declination (deg,min,sec)"},
   {"RVCORR",   TYPEINT, &izero,    "Barycentric RV correction (km/s)"},
   {"JDCORR",   TYPEINT, &izero,    "Barycentric JD correction (days)"},
   {"JD_BARY",  TYPEINT, &izero,   "Corrected barycentric JD at mid-exposure"},
   {"",0,NULL,""}};

void read_descriptors_cfg(void)
{
   char **keys,**vals;
   char cfg[] = "descriptors.cfg";
   int n,i;
   keydef_type *kptr;

   if (!file_exists(cfg)) return;

   keys = allocate_cfglines();
   vals = allocate_cfglines();

   read_matching_cfg(cfg,filename,keys,vals,&n);

   for (i=0;i<n;i++)
   {
      if (strlen(keys[i]) > MAX_FITS_KEYLEN) get_error
         ("FITS keyword '%s' too long in '%s'!",keys[i],cfg);
      if (strlen(vals[i]) > MAX_FITS_VALLEN) get_error
         ("FITS keyword value '%s' too long in '%s'!",vals[i],cfg);
      for (kptr=descriptor_list;
         *kptr->item != 0 && strcasecmp(keys[i],kptr->item) != 0;kptr++);
      if (*kptr->item == (char)0) get_error
         ("Unexpected keyword name '%s' in '%s'!",keys[i],cfg);
      populate_values(kptr,vals[i]);
   }

   free_cfglines(keys);
   free_cfglines(vals);
}

void assign_descriptors(void)
{
   int i;
   keydef_type k;

   for (i=0;i<sky_descriptor;i++)
   {
      if (descriptor_list[i].found)
         k = descriptor_list[i];
      else
         get_keydef_item(i,&k);
      switch (i)
      {
         case 0: strcpy(starname,k.sval);
                 break;
         case 1: strcpy(source,k.sval);
                 break;
         case 2: strcpy(expdate,k.sval);
                 break;
         case 3: strcpy(expstart,k.sval);
                 break;
         case 4: expsec = k.dval;
                 break;
         case 5: expmid = k.dval;
                 break;
         case 6: adu_bias = k.dval;
                 break;
         case 7: strcpy(ccdname,k.sval);
                 get_lower_case(ccdname);
                 break;
         case 8: strcpy(specname,k.sval);
                 get_lower_case(specname);
                 break;
         case 9: gainset = k.ival;
                 break;
         case 10: fibreno = k.ival;
                  break;
      }
   }
   k = descriptor_list[sky_descriptor];
   if (k.found)
      strcpy(skyname,k.sval);
   else
      strcpy(skyname,"sky");
}

void check_descriptors(void)
{
   ccd_type ccd;
   spec_type spec;
   int i,gain_ok;

   if (strlen(specname)<1 || strlen(specname)>MAX_SPECNAME)
      get_error("Bad spectrograph name (%s)!",specname);
   load_spec_info(specname,&spec);
   if (strlen(ccdname)<1 || strlen(ccdname)>MAX_CCDNAME)
      get_error("Bad CCD camera name (%s)!",ccdname);
   load_ccd_info(ccdname,specname,&ccd);
   for (i=1,gain_ok=0;i<=ccd.ngain;i++)
     if (gainset == ccd.gain[i].gainset) gain_ok = -1;
   if (!gain_ok) get_error("Bad CCD gain setting (%d)!",gainset);
   if (fibreno < 1 || fibreno > spec.nfib)
      get_error("Bad fibre number (%d)!",fibreno);
   if (strlen(starname)<1 || strlen(starname)>MAX_STARNAME)
      get_error("Bad star name (%s)!",starname);
   if (strlen(source)<1 || strlen(source)>MAX_SOURCE)
      get_error("Bad source name (%s)!",source);
   if (strlen(skyname)<1 || strlen(skyname)>MAX_SKYNAME)
      get_error("Bad sky name (%s)!",skyname);
   if (strlen(expdate) != MAX_DATE)
      get_error("Bad exposure date (%s)!",expdate);
   if (strlen(expstart) < 8)
      get_error("Bad exposure start time (%s)!",expstart);
   if (expsec <= 0) get_error("Bad exposure time (%f)!",expsec);
   if (expmid <= 0) expmid = expsec / 2;
   if (adu_bias < 0) get_error("Bad ADU bias (%f)!",adu_bias);
}

void calculate_mid_time(void)
{
   int hour,min,day,month,year;
   double sec;
   datetime_type utc_start;

   sscanf(expstart,"%d:%d:%lf",&hour,&min,&sec);
   sscanf(expdate,"%d-%d-%d",&year,&month,&day);
   set_datetime(&utc_start,day,month,year,hour,min,sec);
   shift_datetime(&utc_start,&utc,expmid);

   sprintf(midtime,"%02d:%02d:%04.1f",utc.time.hour,utc.time.min,utc.time.sec);
}

void write_descriptor_block(void)
{
   descriptor_type *h;
   int block_found;

   block_found = (find_fits_comment(&xfits.head,block_head) > 0);
   if (!block_found)
   {
      insert_blank_card(&xfits,&xfits.head);
      insert_comment_card(&xfits,&xfits.head,block_head);
      insert_blank_card(&xfits,&xfits.head);
   }

   for (h=descriptor_block;*h->name != 0;h++)
   {
      switch(h->type)
      {
         case TYPESTR: write_keyword_textual
           (&xfits,&xfits.head,h->name,(char *)h->val,h->comm);
            break;
         case TYPEINT: write_keyword_integer
           (&xfits,&xfits.head,h->name,*((int *)h->val),h->comm);
            break;
         case TYPEDBL: write_keyword_double
           (&xfits,&xfits.head,h->name,*((double *)h->val),h->comm);
            break;
         case TYPELOG: write_keyword_logical
           (&xfits,&xfits.head,h->name,*((char *)h->val),h->comm);
            break;
      }
   }

   if (!block_found)
   {
      insert_blank_card(&xfits,&xfits.head);
      insert_comment_card(&xfits,&xfits.head,block_bar);
      insert_blank_card(&xfits,&xfits.head);
   }
}

int main(int argc,char **argv)
{
   if (argc != 2) get_error("Usage:  hrsp_descriptors  <img>");

   logfile("hrsp_descriptors %s",argv[1]);

   read_raw_image(&xfits,argv[1]);
   ensure_not_original_file(&xfits);
   get_keyword_filename(&xfits,filename);
   read_keydef_cfg(&xfits);
   read_descriptors_cfg();
   assign_descriptors();
   check_descriptors();
   calculate_mid_time();
   write_descriptor_block();
   write_raw_image(&xfits,argv[1]);

   check_memory();

   return(0);
}
