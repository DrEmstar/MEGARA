/* ------------------------------------------------------------------------

   Program:  astronomical_almanac
   Purpose:  Generate some tables from "Astronomical Almanac"

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
#include "nutation.h"
#include "astro.h"

#define SEC_C     10000L
#define MIN_C    600000L
#define HOUR_C 36000000L
#define DAY_C 864000000L

void print_astro_date(date_type *date,int year)
{
  switch(date->year - year)
  {
    case -1: printf("Jan %2d",(int)(date->jd - date_to_jd(31,12,year-1)));
             break;
    case  0: printf("%3s %2d",short_month_name(date->month),date->day);
             break;
    case  1: printf("Dec %2d",(int)(date->jd - date_to_jd(30,11,year)));
             break;
  }
}

void print_jul_day(double jd)
{
  printf("%6.1f",fmod(jd,10000.0));
}

void print_rot_matrix(mat r)
{
  int i,j,m;
  double a;

  for (i=0;i<3;i++)
  {
    for (j=0;j<3;j++)
    {
      if (i == j)
        a = r[i][j] - 1.0;
      else
        a = r[i][j];

      m = (int)floor(a * 1e10 + 0.5);
      printf("%+10d",m);
    }
  }
}

void split_hmsf(double decsec,int prec,int *h,int *m,int *s,int *f)
{
  double c,hour_c,min_c,u;

  if (prec < 1) prec = 1;
  if (prec > 6) prec = 6;
  c = pow(10.0,(double)prec);
  hour_c = c * 3600.0;
  min_c = c * 60.0;

  u = floor(decsec * c + 0.5);
  *h = (int)floor(u/hour_c);
  u = fmod(u,hour_c);
  *m = (int)floor(u/min_c);
  u = fmod(u,min_c);
  *s = (int)floor(u/c);
  *f = (int)fmod(u,c);
}

void ut_st(int year)
{
  int32_t jdlast;
  epoch_type epoch,next;
  int mh,mm,ms,mf;
  int ah,am,as,af;
  int h,m,s,f;
  double deltas;
  datetime_type ut1,tdt;
  int i;
  double gsd;

  jdlast = date_to_jd(1,1,year+1);
  set_epoch_ut1(&epoch,31,12,year-1,0,0,0);
  while (epoch.ut1.date.jd <= jdlast)
  {
    print_astro_date(&epoch.ut1.date,year);
    printf("  ");
    print_jul_day(epoch.ut1.date.jd0);
    split_hmsf(epoch.gmst.decsec,4,&mh,&mm,&ms,&mf);
    split_hmsf(epoch.gast.decsec,4,&ah,&am,&as,&af);
    if (fabs(epoch.gmst.frac - epoch.gast.frac) < 0.5)
    {
      if (epoch.gmst.frac < epoch.gast.frac)
        h = mh, m = mm;
      else
        h = ah, m = am;
    }
    else
    {
      if (epoch.gmst.frac > epoch.gast.frac)
        h = mh, m = mm;
      else
        h = ah, m = am;
    }
    if (mm != m) ms += 60;
    if (am != m) as += 60;
    printf("   %2d %02d %02d.%04d   %02d.%04d",h,m,as,af,ms,mf);
    printf("   %+8.4f",rad2tsec(epoch.equeq));
    add_ut1_seconds(&epoch,&next,86400.0);
    deltas = (next.gmst.frac - epoch.gmst.frac) + 1.0;
    if (deltas < 1.0) deltas += 1.0;
    shift_datetime(&epoch.ut1,&ut1,-epoch.gmst.decsec/deltas);
    for (i=1;i<3;i++)
    {
      shift_datetime(&ut1,&ut1,FULL_TSEC/deltas);
      shift_datetime(&ut1,&tdt,epoch.deltat);
      if (ut1.date.jd == epoch.ut1.date.jd)
      {
        if (i > 1) printf("%51s"," ");
        gsd = floor(get_gsd(&ut1,&tdt) + 0.5);
        printf("  ");
        print_jul_day(gsd);
        printf("  ");
        print_astro_date(&ut1.date,year);
        split_hmsf(ut1.time.decsec,4,&h,&m,&s,&f);
        printf("  %2d %02d %02d.%04d",h,m,s,f);
        printf("\n");
      }
    }
    epoch = next;
  }
}

void ut_era(int year)
{
  int32_t jdlast;
  epoch_type epoch;
  double a,q,d,m,s;

  jdlast = date_to_jd(1,1,year+1);
  set_epoch_ut1(&epoch,31,12,year-1,0,0,0);
  while (epoch.ut1.date.jd <= jdlast)
  {
    print_astro_date(&epoch.ut1.date,year);
    printf("  ");
    print_jul_day(epoch.ut1.date.jd0);
    a = floor(rad2asec(epoch.erot) * 1e4 + 0.5);
    s = modf(a/6e5,&q) * 60.0;
    m = modf(q/60.0,&d) * 60;
    printf("  %3.0f %02.0f %07.4f",d,m,s);
    a = floor(rad2asec(fabs(epoch.equorig)) * 1e4 + 0.5);
    s = modf(a/6e5,&m) * 60.0;
    m += 1.0e-6;
    if (epoch.equorig < 0.0) m = -m;
    printf("  %+2.0f %07.4f",m,s);
    printf("\n");
    add_ut1_seconds(&epoch,&epoch,86400.0);
  }
}

void gcrs_equ(int year)
{
  int32_t jdlast;
  epoch_type epoch;

  jdlast = date_to_jd(1,1,year+1);
  set_epoch_tdt(&epoch,31,12,year-1,0,0,0);
  while (epoch.tdt.date.jd <= jdlast)
  {
    print_astro_date(&epoch.tdt.date,year);
    print_rot_matrix(epoch.gcrs_to_true_eq);
    printf("\n");
    add_tdt_seconds(&epoch,&epoch,86400.0);
  }
}

void gcrs_cio(int year)
{
  int32_t jdlast;
  epoch_type epoch;

  jdlast = date_to_jd(1,1,year+1);
  set_epoch_tdt(&epoch,31,12,year-1,0,0,0);
  while (epoch.tdt.date.jd <= jdlast)
  {
    print_jul_day(epoch.tdt.date.jd0);
    print_rot_matrix(epoch.gcrs_to_cio);
    printf("\n");
    add_tdt_seconds(&epoch,&epoch,86400.0);
  }
}

void nut_obliq(int year)
{
  int32_t jdlast;
  epoch_type epoch;

  jdlast = date_to_jd(1,1,year+1);
  set_epoch_tdt(&epoch,31,12,year-1,0,0,0);
  while (epoch.tdt.date.jd <= jdlast)
  {
    print_astro_date(&epoch.tdt.date,year);
    printf("%+10.4f%+10.4f%10.4f    ",
      rad2asec(epoch.delpsi),rad2asec(epoch.deleps),
      fmod(rad2asec(epoch.epsilon+epoch.deleps),60.0));
    print_jul_day(epoch.tdt.date.jd0);
    printf("    %+10.4f%+10.4f%+10.4f",
      rad2asec(epoch.xpol),rad2asec(epoch.ypol),rad2asec(epoch.orig));
    printf("\n");
    add_tdt_seconds(&epoch,&epoch,86400.0);
  }
}

void earth_xyz(int year)
{
  int32_t jdlast;
  epoch_type epoch;

  jdlast = date_to_jd(1,1,year+1);
  set_epoch_tdb(&epoch,31,12,year-1,0,0,0);
  while (epoch.tdb.date.jd <= jdlast)
  {
    print_astro_date(&epoch.tdb.date,year);
    printf("%+15.9f%+15.9f%+15.9f%+12d%+12d%+12d",
      epoch.reb[0],epoch.reb[1],epoch.reb[2],
      (int)floor(epoch.veb[0]*1e9+0.5),
      (int)floor(epoch.veb[1]*1e9+0.5),
      (int)floor(epoch.veb[2]*1e9+0.5));
    printf("\n");
    add_tdb_seconds(&epoch,&epoch,86400.0);
  }
}

int main(int argc,char **argv)
{
  int year,tblno;

  if (argc != 3) get_error("Usage:  astronomical_almanac <year> <tblno>");

  year = atoi(argv[1]);
  tblno = atoi(argv[2]);

  load_astro();

  switch (tblno)
  {
    case 1: ut_st(year);
            break;
    case 2: ut_era(year);
            break;
    case 3: gcrs_equ(year);
            break;
    case 4: gcrs_cio(year);
            break;
    case 5: nut_obliq(year);
            break;
    case 6: earth_xyz(year);
            break;
  }

  free_astro();

  check_memory();

  return(0);
}
