/* ------------------------------------------------------------------------

   Program:  reduce_star
   Purpose:  Reduce a given stellar spectrum.

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
   int jd,tha,star,thb;
   char rnamea[9],rname[9],rnameb[9],cname[9],xname[9],bname[9];
   char pnamea[9],pname[9],pnameb[9];
   char hnames[10],enames[10],unames[10],nnames[10];
   char enamea[10],enameb[10],unamea[10],unameb[10];
   char tnamea[10],tnameb[10],tnames[10];
   char wnames[10],gnames[10];

   if (argc != 5)
      get_error("Usage:  reduce_star  <JD> <tha> <star> <thb>");

   inform("");
   inform("Reduction of stellar image");
   inform("==========================");

   jd = collect_jd_number(argv[1]);
   tha = collect_image_number(argv[2]);
   star = collect_image_number(argv[3]);
   thb = collect_image_number(argv[4]);

   sprintf(rnamea,"r%04d%03d",jd,tha);
   sprintf(rname, "r%04d%03d",jd,star);
   sprintf(rnameb,"r%04d%03d",jd,thb);
 
   sprintf(cname,"c%04d%03d",jd,star);
   sprintf(xname,"x%04d%03d",jd,star);
   sprintf(bname,"b%04d%03d",jd,star);
 
   sprintf(pnamea,"p%04d%03d",jd,tha);
   sprintf(pname, "p%04d%03d",jd,star);
   sprintf(pnameb,"p%04d%03d",jd,thb);
 
   sprintf(hnames,"h%04d%03ds",jd,star);
   sprintf(enames,"e%04d%03ds",jd,star);
   sprintf(unames,"u%04d%03ds",jd,star);
   sprintf(nnames,"n%04d%03ds",jd,star);
 
   sprintf(enamea,"e%04d%03da",jd,star);
   sprintf(enameb,"e%04d%03db",jd,star);
   sprintf(unamea,"u%04d%03da",jd,star);
   sprintf(unameb,"u%04d%03db",jd,star);

   sprintf(tnamea,"t%04d%03da",jd,star);
   sprintf(tnameb,"t%04d%03db",jd,star);
   sprintf(tnames,"t%04d%03ds",jd,star);
 
   sprintf(wnames,"w%04d%03ds",jd,star);
   sprintf(gnames,"g%04d%03ds",jd,star);

   get_system("identify_object %s",rname);
   get_system("compute_barycentric %s",rname);

   get_system("vertical_offset %s",rname);
   get_system("extract_echelle %s %s",rname,hnames);
   get_system("filter_star %s %s %s %s",rname,cname,xname,bname);
   get_system("flat_full %s %s",cname,pname);
   get_system("extract_echelle %s %s",pname,enames);
   get_system("flat_extracted %s %s",enames,unames);
   get_system("normalize_spectrum %s %s",unames,nnames);

   get_system("copy_keyword %s %s %s",rname,rnamea,"VOFFSET");
   get_system("flat_full %s %s",rnamea,pnamea);
   get_system("extract_echelle %s %s",pnamea,enamea);
   get_system("flat_extracted %s %s",enamea,unamea);

   get_system("copy_keyword %s %s %s",rname,rnameb,"VOFFSET");
   get_system("flat_full %s %s",rnameb,pnameb);
   get_system("extract_echelle %s %s",pnameb,enameb);
   get_system("flat_extracted %s %s",enameb,unameb);

   get_system("thorium_lines %s %s",unamea,tnamea);
   get_system("thorium_lines %s %s",unameb,tnameb);
   get_system("dispersion_solution %s",tnamea);
   get_system("dispersion_solution %s",tnameb);
   get_system("interpolate_solution %s %s %s",tnamea,tnameb,nnames);

   get_system("rebin_echelle %s %s %s",nnames,wnames,gnames);

   check_memory();

   return(0);
}
