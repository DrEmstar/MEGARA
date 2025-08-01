#ifndef ANGLE_DEFINED

   #define ANGLE_DEFINED

   #define TWOPI   6.28318530717958648
   #define PI      3.14159265358979324
   #define HALFPI  1.57079632679489662

   typedef struct
   {
      double hour,min,sec;
      double frac,rad;
   }
   alphatype;

   typedef struct
   {
      double deg,min,sec;
      double frac,rad;
   }
   thetatype;

   typedef struct
   {
      char sgn;
      double hour,min,sec;
      double frac,rad;
   }
   lambdatype;

   typedef struct
   {
      char sgn;
      double deg,min,sec;
      double frac,rad;
   }
   deltatype;

   typedef struct
   {
      thetatype theta;
      alphatype alpha;
      deltatype delta;
      lambdatype lambda;
   }
   angletype;

#endif

/* ------------------------ Function Declaration ------------------------ */
double range(double x,double r);
double getrad(double x);
void break_angle(double frac,char *z,double *h,double *m,double *s);
double glue_angle(char z,double h,double m,double s);
angletype form_angle(double rad);
angletype form_delta(char sgn,double deg,double min,double sec);
angletype form_lambda(char sgn,double hour,double min,double sec);
angletype form_alpha(double hour,double min,double sec);
