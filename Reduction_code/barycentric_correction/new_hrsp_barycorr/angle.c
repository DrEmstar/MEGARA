#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"

double range(double x,double r)
{
   double u;

   u = fmod(x,r);
   if (u < 0) u+=r;
   return(u);
}

double getrad(double x)
{
   return(x*PI/180.0);
}

void break_angle(double frac,char *z,double *h,double *m,double *s)
{
   *z=sign_character(frac);
   frac=fabs(frac);
   *h=floor(frac);
   frac=(frac-*h)*60;
   *m=floor(frac);
   *s=(frac-*m)*60;
   *s = floor((*s)*10000+0.5)/10000;
   if (*s == 60) *s=0,*m+=1.0;
   if (*m == 60) *m=0,*h+=1.0;
}

double glue_angle(char z,double h,double m,double s)
{
   if (z == '-')
      return(-(h+m/60+s/3600));
   else
      return(h+m/60+s/3600);
}

angletype form_angle(double rad)
{
   angletype w;
   char z;

   w.theta.rad = range(rad,TWOPI);
   w.theta.frac = w.theta.rad*180/PI;
   break_angle(w.theta.frac,&z,&w.theta.deg,&w.theta.min,&w.theta.sec);

   w.alpha.rad = w.theta.rad;
   w.alpha.frac = w.alpha.rad*12/PI;
   break_angle(w.alpha.frac,&z,&w.alpha.hour,&w.alpha.min,&w.alpha.sec);

   if (w.theta.rad > PI)
      w.delta.rad = w.theta.rad - TWOPI;
   else
      w.delta.rad = w.theta.rad;
   w.delta.frac = w.delta.rad*180/PI;
   break_angle(w.delta.frac,
               &w.delta.sgn,&w.delta.deg,&w.delta.min,&w.delta.sec);

   w.lambda.rad = w.delta.rad;
   w.lambda.frac = w.lambda.rad*12/PI;
   break_angle(w.lambda.frac,
               &w.lambda.sgn,&w.lambda.hour,&w.lambda.min,&w.lambda.sec);

   return(w);
}

angletype form_delta(char sgn,double deg,double min,double sec)
{
   return(form_angle(glue_angle(sgn,deg,min,sec)*PI/180));
}

angletype form_lambda(char sgn,double hour,double min,double sec)
{
   return(form_angle(glue_angle(sgn,hour,min,sec)*PI/12));
}

angletype form_alpha(double hour,double min,double sec)
{
   return(form_angle(glue_angle(' ',hour,min,sec)*PI/12));
}
