#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "datetime.h"

#define GREG_JD 2299161L /* 15 October 1582 */

monthtype months[] =
   { {"January",   "Jan" },
     {"February",  "Feb" },
     {"March",     "Mar" },
     {"April",     "Apr" },
     {"May",       "May" },
     {"June",      "Jun" },
     {"July",      "Jul" },
     {"August",    "Aug" },
     {"September", "Sep" },
     {"October",   "Oct" },
     {"November",  "Nov" },
     {"December",  "Dec" } };

/* ----------------------------- Month Names  ----------------------------- */

char *full_month_name(int m)
{
  if (m>=1 && m<=12)
    return(months[m-1].full);
  else
    return(NULL);
}

char *short_month_name(int m)
{
  if (m>=1 && m<=12)
    return(months[m-1].abbr);
  else
    return(NULL);
}

int month_number(char *mname)
{
  int m;

  for (m=0;m<12;m++)
    if (strncasecmp(mname,months[m].abbr,3) == 0) return(m+1);

  return(0);
}

/* ----------------------------- Julian Days  ----------------------------- */

int32_t gregorian_date_to_jd(int d,int m,int y)
{
   int32_t mm,jd;

   mm = ((int32_t)m - 14L) / 12L;

   jd =  (1461L * ((int32_t)y + 4800L + mm)) / 4L
       + (367L * ((int32_t)m - 2L - mm * 12L)) / 12L
       - (3L * (((int32_t)y + 4900L + mm) / 100L)) / 4L
       + (int32_t)d - 32075L;

   return(jd);
}

int32_t julian_date_to_jd(int d,int m,int y)
{
   int32_t jd;

   jd =  367L * (int32_t)y
       - (7L * ((int32_t)y + 5001L + ((int32_t)m - 9L) / 7L)) / 4L
       + (275L * (int32_t)m) / 9L
       + (int32_t)d + 1729777L;

   return(jd);
}

int32_t date_to_jd(int d,int m,int y)
{
   int32_t jd;

   jd = gregorian_date_to_jd(d,m,y);
   if (jd >= GREG_JD)
      return(jd);
   else
      return(julian_date_to_jd(d,m,y));
}

void jd_to_gregorian_date(int32_t jd,int *d,int *m,int *y)
{
   int32_t h,n,i,j,b;

   h = jd + 68569L;
   n = (4L * h) / 146097L;
   h = h - (146097L * n + 3L) / 4L;
   i = (4000L * (h + 1L)) / 1461001L;
   h = h - (1461L * i) / 4L + 31L;
   j = (80L * h) / 2447L;
   b = j / 11L;
   *d = (int)(h - (2447L * j) / 80L);
   *m = (int)(j + 2L - 12L * b);
   *y = (int)(100L * (n - 49L) + i + b);
}

void jd_to_julian_date(int32_t jd,int *d,int *m,int *y)
{
   int32_t h,n,i,j,k,b;

   j = jd + 1402L;
   k = (j - 1L) / 1461L;
   h = j - 1461L * k;
   n = (h - 1L) / 365L - h / 1461L;
   i = h - 365L * n + 30L;
   j = (80L * i) / 2447L;
   b = j / 11L;
   *d = (int)(i - (2447L * j) / 80L);
   *m = (int)(j + 2L - 12L * b);
   *y = (int)(4L * k + n + b - 4716L);
}

void jd_to_date(int32_t jd,int *d,int *m,int *y)
{
   if (jd >= GREG_JD)
      jd_to_gregorian_date(jd,d,m,y);
   else
      jd_to_julian_date(jd,d,m,y);
}

/* --------------------------------- Date --------------------------------- */

void set_date_dmy(date_type *u,int d,int m,int y)
{
  int32_t jd;
  int dd,mm,yy;

  jd = date_to_jd(d,m,y);
  jd_to_date(jd,&dd,&mm,&yy);
  if (dd != d || mm != m || yy != y)
    get_error("set_date: Invalid date!");
  u->day = d;
  u->month = m;
  u->year = y;
  u->jd = jd;
  u->jd0 = (double)jd - 0.5;
}

void set_date_jd(date_type *u,int32_t jd)
{
  jd_to_date(jd,&u->day,&u->month,&u->year);
  u->jd = jd;
  u->jd0 = (double)jd - 0.5;
}

void add_days(date_type *a,date_type *b,int d)
{
  set_date_jd(b,a->jd + (int32_t)d);
}

void next_day(date_type *a,date_type *b)
{
  add_days(a,b,1);
}

/* --------------------------------- Time --------------------------------- */

void set_time_decsec_adj(time_type *t,double decsec,int adj)
{
  if (abs(adj) > 1)
    get_error("set_time_decsec_adj: Invalid day length adjustment!");
  t->adj = adj;
  t->full_day = (double)(86400 + adj);

  if (decsec < 0.0 || decsec >= t->full_day)
    get_error("set_time_decsec_adj: Time out of range!");
  t->decsec = decsec;
  t->frac = decsec / t->full_day;
  t->rad = t->frac * TWOPI;
  if (decsec < FULL_TSEC)
    split_decsec_plus(decsec,&t->hour,&t->min,&t->sec);
  else
  {
    t->hour = 23;
    t->min = 59;
    t->sec = (decsec - FULL_TSEC) + 60.0;
  }
}

void set_time_decsec(time_type *t,double decsec)
{
  set_time_decsec_adj(t,decsec,0);
}

void set_time_hms_adj(time_type *t,int hour,int min,double sec,int adj)
{
  double lim;

  if (abs(adj) > 1)
    get_error("set_time_hms_adj: Invalid day length adjustment!");
  if (hour < 0 || hour > 23)
    get_error("set_time_hms_adj: Hours out of range!");
  if (min < 0 || min > 59)
    get_error("set_time_hms_adj: Minutes out of range!");
  if (hour == 23 && min == 59)
    lim = (double)(60 + adj);
  else
    lim = 60.0;
  if (sec < 0.0 || sec >= lim)
    get_error("set_time_hms_adj: Seconds out of range!");
  set_time_decsec_adj(t,build_decsec_plus(hour,min,sec),adj);
}

void set_time_hms(time_type *t,int hour,int min,double sec)
{
  set_time_hms_adj(t,hour,min,sec,0);
}

void set_time_frac_adj(time_type *t,double frac,int adj)
{
  if (frac < 0.0 || frac >= 1.0)
    get_error("set_time_frac_adj: Day fraction out of range!");
  set_time_decsec_adj(t,frac * (double)(86400+adj),adj);
}

void set_time_frac(time_type *t,double frac)
{
  set_time_frac_adj(t,frac,0);
}

/* ------------------------------- DateTime ------------------------------- */

void update_time_args(datetime_type *u)
{
  u->jd0 = u->date.jd0;
  u->jdf = u->time.frac;
  u->jd = u->jd0 + u->jdf;
  u->djd = (u->jd0 - 2451545.0) + u->jdf;
  u->jcen = u->djd / 36525.0;
  u->jyr = 2000.0 + u->djd / 365.25;
}

void set_datetime_adj(datetime_type *u,int day,int month,int year,
  int hour,int min,double dsec,int adj)
{
  set_date_dmy(&u->date,day,month,year);
  set_time_hms_adj(&u->time,hour,min,dsec,adj);
  update_time_args(u);
}

void set_datetime(datetime_type *u,int day,int month,int year,
  int hour,int min,double dsec)
{
  set_datetime_adj(u,day,month,year,hour,min,dsec,0);
}

void split_jd(double jda,double jdb,double *jd0,double *jdf)
{
  double af,ai,bf,bi,cf,ci;

  af = modf(jda,&ai);
  bf = modf(jdb,&bi);

  cf = modf(af+bf+2.5,&ci);
  *jd0 = ai + bi + ci - 2.5;
  *jdf = cf;
}

void set_datetime_jd(datetime_type *u,double jda,double jdb)
{
  double jd0,jdf;

  split_jd(jda,jdb,&jd0,&jdf);
  set_date_jd(&u->date,(int32_t)(jd0+0.5));
  set_time_frac(&u->time,jdf);
  update_time_args(u);
}

void set_julian_year(datetime_type *u,double year)
{
  double jda,jdb;

  jda = 2451545.0;
  jdb = (year - 2000.0) * 365.25;
  set_datetime_jd(u,jda,jdb);
}

void set_besselian_year(datetime_type *u,double year)
{
  double jda,jdb;

  jda = 2415020.0;
  jdb = 0.31352 + (year - 1900.0) * 365.242198781;
  set_datetime_jd(u,jda,jdb);
}

void shift_datetime(datetime_type *a,datetime_type *b,double sec)
{
  double decsec,days;

  decsec = a->time.decsec + sec;

  days = floor(decsec/FULL_TSEC);
  decsec = full_range(decsec,FULL_TSEC);

  set_date_jd(&b->date,a->date.jd + (int32_t)days);
  set_time_decsec(&b->time,decsec);
  update_time_args(b);
}

int compare_datetime(datetime_type *a,datetime_type *b)
{
  if (a->date.jd == b->date.jd)
  {
    if (a->time.decsec == b->time.decsec) return(0);
    if (a->time.decsec < b->time.decsec)
      return(-1);
    else
      return(1);
  }
  else
  {
    if (a->date.jd < b->date.jd)
      return(-1);
    else
      return(1);
  }
}

/* ------------------------------- Interval ------------------------------- */

void set_interval_years(interval_type *t,double years)
{
  t->julian_cents = years / 100.0;
  t->julian_years = years;
  t->days = years * 365.25;
  t->seconds = years * 31557600.0;
}

void set_interval_days(interval_type *t,double days)
{
  t->julian_cents = days / 36525.0;
  t->julian_years = days / 365.25;
  t->days = days;
  t->seconds = days * 86400.0;
}
