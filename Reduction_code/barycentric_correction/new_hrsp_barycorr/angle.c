#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"

double full_range(double x,double r)
{
  x = fmod(x,r);
  if (x < 0.0) x += r;
  if (x == r) x = 0.0;
  return(x);
}

double half_range(double x,double h)
{
  double r;

  r = 2.0 * h;
  x = fmod(x,r);
  if (x < -h) x += r;
  if (x >= h) x -= r;
  return(x);
}

double full_frac(double frac)
{
  return(full_range(frac,1.0));
}

double half_frac(double frac)
{
  return(half_range(frac,0.5));
}

double build_decsec_plus(int deg,int min,double sec)
{
  return(((double)deg * 3600.0 + (double)min * 60.0) + sec);
}

double build_decsec(char sgn,int deg,int min,double sec)
{
  if (sgn == '-')
    return(-build_decsec_plus(deg,min,sec));
  else
    return(build_decsec_plus(deg,min,sec));
}

void split_decsec(double decsec,char *sgn,int *deg,int *min,double *sec)
{
  double fsec,isec;

  fsec = modf(fabs(decsec),&isec);

  *sgn = sign_char(decsec);
  *deg = (int)floor(isec/3600.0);
  isec = fmod(isec,3600.0);
  *min = (int)floor(isec/60.0);;
  *sec = fmod(isec,60.0) + fsec;
}

void split_decsec_plus(double decsec,int *deg,int *min,double *sec)
{
  char sgn;

  split_decsec(decsec,&sgn,deg,min,sec);
}

/* ------------------------------- Radians  ------------------------------- */

double deg2rad(double x)
{
  return(x * DEG2R);
}

double hr2rad(double x)
{
  return(x * HOUR2R);
}

double asec2rad(double x)
{
  return(x * ASEC2R);
}

double tsec2rad(double x)
{
  return(x * TSEC2R);
}

double rad2deg(double x)
{
  return(x * R2DEG);
}

double rad2hr(double x)
{
  return(x * R2HOUR);
}

double rad2asec(double x)
{
  return(x * R2ASEC);
}

double rad2tsec(double x)
{
  return(x * R2TSEC);
}

/* -------------------------------- Theta  -------------------------------- */

void set_theta_frac(theta_type *the,double frac)
{
  the->frac = full_frac(frac);
  the->rad = the->frac * TWOPI;
  the->decsec = the->frac * FULL_ASEC;
  split_decsec_plus(the->decsec,&the->deg,&the->min,&the->sec);
}

void set_theta_rad(theta_type *the,double rad)
{
  set_theta_frac(the,rad/TWOPI);
}

void set_theta_decsec(theta_type *the,double decsec)
{
  set_theta_frac(the,decsec/FULL_ASEC);
}

void set_theta_dms(theta_type *the,int deg,int min,double sec)
{
  if (deg < 0 || deg >= 360)
    get_error("set_theta_dms: Degrees out of range!");
  if (min < 0 || min >= 60)
    get_error("set_theta_dms: Minutes out of range!");
  if (sec < 0.0 || sec >= 60.0)
    get_error("set_theta_dms: Seconds out of range!");
  set_theta_decsec(the,build_decsec_plus(deg,min,sec));
}

/* -------------------------------- Alpha  -------------------------------- */

void set_alpha_frac(alpha_type *alp,double frac)
{
  alp->frac = full_frac(frac);
  alp->rad = alp->frac * TWOPI;
  alp->decsec = alp->frac * FULL_TSEC;
  split_decsec_plus(alp->decsec,&alp->hour,&alp->min,&alp->sec);
}

void set_alpha_rad(alpha_type *alp,double rad)
{
  set_alpha_frac(alp,rad/TWOPI);
}

void set_alpha_decsec(alpha_type *alp,double decsec)
{
  set_alpha_frac(alp,decsec/FULL_TSEC);
}

void set_alpha_hms(alpha_type *alp,int hour,int min,double sec)
{
  if (hour < 0 || hour >= 24)
    get_error("set_alpha_hms: Hours out of range!");
  if (min < 0 || min >= 60)
    get_error("set_alpha_hms: Minutes out of range!");
  if (sec < 0.0 || sec >= 60.0)
    get_error("set_alpha_hms: Seconds out of range!");
  set_alpha_decsec(alp,build_decsec_plus(hour,min,sec));
}

char *print_alpha(alpha_type *alp,char *str,int prec)
{
  sprintf(str,"%2d %02d %0*.*f",alp->hour,alp->min,prec+3,prec,alp->sec);
  return(str);
}

/* -------------------------------- Delta  -------------------------------- */

void set_delta_frac(delta_type *del,double frac)
{
  del->frac = half_frac(frac);
  del->rad = del->frac * TWOPI;
  del->decsec =del->frac * FULL_ASEC;
  split_decsec(del->decsec,&del->sgn,&del->deg,&del->min,&del->sec);
}

void set_delta_rad(delta_type *del,double rad)
{
  set_delta_frac(del,rad/TWOPI);
}

void set_delta_decsec(delta_type *del,double decsec)
{
  set_delta_frac(del,decsec/FULL_ASEC);
}

void set_delta_dms(delta_type *del,char sgn,int deg,int min,double sec)
{
  if (sgn != '+' && sgn != '-' && sgn != ' ')
    get_error("set_delta_dms: Invalid sign character!");
  if (deg < 0 || deg > 180)
    get_error("set_delta_dms: Degrees out of range!");
  if (min < 0 || min >= 60)
    get_error("set_delta_dms: Minutes out of range!");
  if (sec < 0.0 || sec >= 60.0)
    get_error("set_delta_dms: Seconds out of range!");
  if (deg == 180 && (min != 0 || sec != 0.0 || sgn != '-'))
    get_error("set_delta_dms: Angle out of range!");
  set_delta_decsec(del,build_decsec(sgn,deg,min,sec));
}

char *print_delta(delta_type *del,char *str,int prec)
{
  char deg[8];

  sprintf(deg,"%c%d",del->sgn,del->deg);
  sprintf(str,"%3s %02d %0*.*f",deg,del->min,prec+3,prec,del->sec);
  return(str);
}

/* -------------------------------- Lambda -------------------------------- */

void set_lambda_frac(lambda_type *lam,double frac)
{
  lam->frac = half_frac(frac);
  lam->rad = lam->frac * TWOPI;
  lam->decsec = lam->frac * FULL_TSEC;
  split_decsec(lam->decsec,&lam->sgn,&lam->hour,&lam->min,&lam->sec);
}

void set_lambda_rad(lambda_type *lam,double rad)
{
  set_lambda_frac(lam,rad/TWOPI);
}

void set_lambda_decsec(lambda_type *lam,double decsec)
{
  set_lambda_frac(lam,decsec/FULL_TSEC);
}

void set_lambda_hms(lambda_type *lam,char sgn,int hour,int min,double sec)
{
  if (sgn != '+' && sgn != '-' && sgn != ' ')
    get_error("set_alpha_hms: Invalid sign character!");
  if (hour < 0 || hour > 12)
    get_error("set_alpha_hms: Hours out of range!");
  if (min < 0 || min >= 60)
    get_error("set_alpha_hms: Minutes out of range!");
  if (sec < 0.0 || sec >= 60.0)
    get_error("set_alpha_hms: Seconds out of range!");
  if (hour == 12 && (min != 0 || sec != 0.0 || sgn != '-'))
    get_error("set_delta_dms: Angle out of range!");
  set_lambda_decsec(lam,build_decsec(sgn,hour,min,sec));
}

/* -------------------------------- Angle  -------------------------------- */

void set_angle_frac(angle_type *a,double frac)
{
  set_theta_frac(&a->theta,frac);
  set_alpha_frac(&a->alpha,frac);
  set_delta_frac(&a->delta,frac);
  set_lambda_frac(&a->lambda,frac);
}

void set_angle_rad(angle_type *a,double rad)
{
  set_angle_frac(a,rad/TWOPI);
}

void set_angle_arcsec(angle_type *a,double arcsec)
{
  set_angle_frac(a,arcsec/FULL_ASEC);
}

void set_angle_timesec(angle_type *a,double timesec)
{
  set_angle_frac(a,timesec/FULL_TSEC);
}

void set_angle_deg(angle_type *a,double deg)
{
  set_angle_frac(a,deg/360.0);
}

void set_angle_theta(angle_type *a,int deg,int min,double sec)
{
  theta_type the;

  set_theta_dms(&the,deg,min,sec);
  set_angle_frac(a,the.frac);
}

void set_angle_alpha(angle_type *a,int hour,int min,double sec)
{
  alpha_type alp;

  set_alpha_hms(&alp,hour,min,sec);
  set_angle_frac(a,alp.frac);
}

void set_angle_delta(angle_type *a,char sgn,int deg,int min,double sec)
{
  delta_type del;

  set_delta_dms(&del,sgn,deg,min,sec);
  set_angle_frac(a,del.frac);
}

void set_angle_lambda(angle_type *a,char sgn,int hour,int min,double sec)
{
  lambda_type lam;

  set_lambda_hms(&lam,sgn,hour,min,sec);
  set_angle_frac(a,lam.frac);
}
