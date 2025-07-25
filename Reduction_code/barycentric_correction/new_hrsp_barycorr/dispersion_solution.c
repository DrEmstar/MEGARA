/* ------------------------------------------------------------------------

   Program:  dispersion_solution
   Purpose:  Compute the dispersion solution for a given thorium table

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
   spec_type spec;
   fits_file linfits;
   lin_table lintbl;
   dispsol_type dsol;
   int seq;

   if (argc != 2) get_error("Usage:  dispersion_solution  <tbl>");

   inform_title("Dispersion solution");

   read_lin_table(&linfits,&lintbl,argv[1]);
   get_specinfo_fits(&linfits,&linfits.extend,&spec);
   get_dispsol_info(&dsol,spec.specname);

   if (dsol.wave.id == AIR_WAVE)
      copy_double_array(lintbl.airlam,linfits.nrows,lintbl.lam);
   else
      copy_double_array(lintbl.vaclam,linfits.nrows,lintbl.lam);
 
   for (seq=0;seq<linfits.nrows;seq++)
   {
      lintbl.mlam[seq] = lintbl.lam[seq] * lintbl.order[seq];
      lintbl.mnorm[seq] = normalized_order(&spec,(double)lintbl.order[seq]);
      lintbl.mlnorm[seq] = normalized_mlambda(&spec,lintbl.mlam[seq]);
      lintbl.xsel[seq] = !lintbl.status[seq];
      if (lintbl.order[seq]<dsol.firstord || lintbl.order[seq]>dsol.lastord)
         lintbl.xsel[seq] = 0;
      if (lintbl.mlam[seq]<dsol.mlamstart || lintbl.mlam[seq]>dsol.mlamend)
         lintbl.xsel[seq] = 0;
   }

   write_table(&linfits,argv[1]);

   inform("Two-dimensional fit: xcen = f(mlnorm,mnorm)");
   regression_polynomial_proc(argv[1],"Xcen","MLnorm,Mnorm",dsol.reg.degarg,
                    "Xsel","Xfit","Xresid","XREGRE",dsol.reg.kappa);

   check_memory();

   return(0);
}
