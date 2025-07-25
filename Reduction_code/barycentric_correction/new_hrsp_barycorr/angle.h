#ifndef ANGLE_DEFINED

  #define ANGLE_DEFINED

  #define TWOPI  6.28318530717958647692528
  #define PI     3.14159265358979323846264
  #define HALFPI 1.57079632679489661923132

  #define R2DEG  57.2957795130823208767981
  #define R2HOUR 3.81971863420548805845321
  #define R2ASEC 206264.806247096355156473
  #define R2TSEC 13750.9870831397570104315

  #define DEG2R  1.74532925199432957692369e-2
  #define HOUR2R 2.61799387799149436538553e-1
  #define ASEC2R 4.84813681109535993589914e-6
  #define TSEC2R 7.27220521664303990384871e-5

  #define FULL_ASEC 1296000.0
  #define FULL_TSEC 86400.0

  typedef struct
  {
    int deg,min;
    double sec,frac,rad,decsec;
  }
  theta_type;

  typedef struct
  {
    int hour,min;
    double sec,frac,rad,decsec;
  }
  alpha_type;

  typedef struct
  {
    char sgn;
    int deg,min;
    double sec,frac,rad,decsec;
  }
  delta_type;

  typedef struct
  {
    char sgn;
    int hour,min;
    double sec,frac,rad,decsec;
  }
  lambda_type;

  typedef struct
  {
    theta_type theta;
    alpha_type alpha;
    delta_type delta;
    lambda_type lambda;
  }
  angle_type;

#endif

double full_range(double x,double r);
double half_range(double x,double h);
double full_frac(double frac);
double half_frac(double frac);
double build_decsec_plus(int deg,int min,double sec);
double build_decsec(char sgn,int deg,int min,double sec);
void split_decsec(double decsec,char *sgn,int *deg,int *min,double *sec);
void split_decsec_plus(double decsec,int *deg,int *min,double *sec);
double deg2rad(double x);
double hr2rad(double x);
double asec2rad(double x);
double tsec2rad(double x);
double rad2deg(double x);
double rad2hr(double x);
double rad2asec(double x);
double rad2tsec(double x);
void set_theta_frac(theta_type *the,double frac);
void set_theta_rad(theta_type *the,double rad);
void set_theta_decsec(theta_type *the,double decsec);
void set_theta_dms(theta_type *the,int deg,int min,double sec);
void set_alpha_frac(alpha_type *alp,double frac);
void set_alpha_rad(alpha_type *alp,double rad);
void set_alpha_decsec(alpha_type *alp,double decsec);
void set_alpha_hms(alpha_type *alp,int hour,int min,double sec);
char *print_alpha(alpha_type *alp,char *str,int prec);
void set_delta_frac(delta_type *del,double frac);
void set_delta_rad(delta_type *del,double rad);
void set_delta_decsec(delta_type *del,double decsec);
void set_delta_dms(delta_type *del,char sgn,int deg,int min,double sec);
char *print_delta(delta_type *del,char *str,int prec);
void set_lambda_frac(lambda_type *lam,double frac);
void set_lambda_rad(lambda_type *lam,double rad);
void set_lambda_decsec(lambda_type *lam,double decsec);
void set_lambda_hms(lambda_type *lam,char sgn,int hour,int min,double sec);
void set_angle_frac(angle_type *a,double frac);
void set_angle_rad(angle_type *a,double rad);
void set_angle_arcsec(angle_type *a,double arcsec);
void set_angle_timesec(angle_type *a,double timesec);
void set_angle_deg(angle_type *a,double deg);
void set_angle_theta(angle_type *a,int deg,int min,double sec);
void set_angle_alpha(angle_type *a,int hour,int min,double sec);
void set_angle_delta(angle_type *a,char sgn,int deg,int min,double sec);
void set_angle_lambda(angle_type *a,char sgn,int hour,int min,double sec);
