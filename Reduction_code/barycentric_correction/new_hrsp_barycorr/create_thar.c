/* ------------------------------------------------------------------------

   Program:  create_thar
   Purpose:  Create a table of Th-Ar lines

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
   fits_file tblfits;
   thar_table table;
   file_type file;
   spec_type spec;
   file_name_type tmp;
   int nrows,seq,ndat,found,ref;
   char species[20];
   regre_block vreg,mreg,ureg;

   if (argc != 2) get_error("Usage:   create_thar <specname>");

   load_spec_info(argv[1],&spec);
   nrows = number_of_lines("%s/spec/%s/dat/thar.dat",hrsp_dir(),spec.specname);
   create_thar_table(nrows,next_tmp_file_name(tmp));
   read_thar_table(&tblfits,&table,tmp);
   remove_tmp_file("%s.fit",tmp);

   open_file(&file,"r","%s/spec/%s/dat/thar.dat",hrsp_dir(),spec.specname);
   for (seq=0,found=0;seq<nrows;seq++)
   {
      ndat = fscanf(file.dat," %d %lf %lf %lf %8s %lf %lf %lf %lf %lf %lf",
                   &table.order[seq],&table.airlam[seq],&table.vaclam[seq],
                   &table.waveno[seq],species,
                   &table.ustart[seq],&table.ucen[seq],&table.uend[seq],
                   &table.vstart[seq],&table.vcen[seq],&table.vend[seq]);
      if (ndat != 11)
         get_error("Reading '%s':\nBad data format!",file.path);
      strcat(species,"        "); /* 8 blanks */
      memmove(table.species+8*seq,species,8);
      if (fabs(table.airlam[seq]-spec.refair) < 1e-8 && 
          fabs(table.vaclam[seq]-spec.refvac) < 1e-8 && 
          table.order[seq] == spec.reford) found++,ref=seq;
   }
   get_fclose(&file);

   if (found == 0)
      get_error("Reading '%s':\nReference line %9.4f not found!",
                 file.path,spec.refair);
   if (found > 1)
      get_error("Reading '%s':\nMore than one reference line %9.4f found!",
                 file.path,spec.refair);

   if (fabs(table.ucen[ref]-spec.uref) > 1e-8 || 
       fabs(table.vcen[ref]-spec.vref) > 1e-8)
      get_error("Reading '%s':\n"
                "Bad coordinates for the reference line %9.4f!",
                 file.path,spec.refair);

   for (seq=0;seq<nrows;seq++)
   {
      table.lam[seq] = table.vaclam[seq];
      table.mlam[seq] = table.order[seq] * table.lam[seq];
      table.unorm[seq] = normalized_u(&spec,table.ucen[seq]);
      table.vnorm[seq] = normalized_v(&spec,table.vcen[seq]);
      table.mnorm[seq] = normalized_order(&spec,(double)table.order[seq]);
      table.mlnorm[seq] = normalized_mlambda(&spec,table.mlam[seq]);
      table.vsel[seq] = table.msel[seq] = table.usel[seq] = 1;
   }

   set_file_name(&file,"%s/spec/%s/tbl/thar",hrsp_dir(),spec.specname);
   write_table(&tblfits,file.path);

   inform("Table '%s' created.",tblfits.file.path);

   inform("Two-dimensional fit: vcen = f(unorm,mnorm)");
   regression_polynomial_proc(file.path,"Vcen","Unorm,Mnorm",
                              spec.echdef_vcen.degarg,"Vsel","Vfit",
                             "Vresid","VREGRE",spec.echdef_vcen.kappa);

   inform("Two-dimensional fit: order = f(unorm,vnorm)");
   regression_polynomial_proc(file.path,"Order","Unorm,Vnorm",
                              spec.echdef_ord.degarg,"Msel","Mfit",
                             "Mresid","MREGRE",spec.echdef_ord.kappa);

   inform("Two-dimensional fit: ucen = f(mlnorm,mnorm)");
   regression_polynomial_proc(file.path,"Ucen","MLnorm,Mnorm",
                              spec.echdef_ucen.degarg,"Usel","Ufit",
                             "Uresid","UREGRE",spec.echdef_ucen.kappa);

   read_table_header(&tblfits,"%s/spec/%s/tbl/thar",hrsp_dir(),spec.specname);

   read_regression_block(&tblfits,&tblfits.extend,"VREGRE",&vreg);
   read_regression_block(&tblfits,&tblfits.extend,"MREGRE",&mreg);
   read_regression_block(&tblfits,&tblfits.extend,"UREGRE",&ureg);

   inform("");
   inform("                     SUMMARY");

   inform_table("   Fit        N_tot     N_sel      R.M.S. Error");

   inform("v = f(u,m)    %4d      %4d     %12.3e",nrows,vreg.nsel,vreg.rms);
   inform("m = f(u,v)    %4d      %4d     %12.3e",nrows,mreg.nsel,mreg.rms);
   inform("u = f(ml,m)   %4d      %4d     %12.3e",nrows,ureg.nsel,ureg.rms);
   inform("");

   check_memory();

   return(0);
}
