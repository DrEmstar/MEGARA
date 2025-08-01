#ifndef ASTRO_DEFINED

   #define ASTRO_DEFINED

   #define LIGHT_SPEED 299792.458

   typedef struct
   {
      double day,month,year;
      double jd;
   }
   datetype;

   typedef struct
   {
      datetype date;
      angletype time;
      double JD,T,yr;
   }
   momentype;

   typedef struct
   {
      momentype ut,tdt;
      double epsilon,delpsi,deleps;
      mat prec,nut,precnut;
      mat invprec,invnut,invprecnut;
      mat mean_to_gal,true_to_gal;
      mat gal_to_mean,gal_to_true;
      angletype gmst,gast;
   }
   epochtype;

   typedef struct
   {
      angletype longitude,latitude;
      double height;
      vec r,v;
      angletype lmst,last;
   }
   topotype;

   typedef struct
   {
      long hd,hip,hr;
      char id[11];
      epochtype equinox,epoch;
      int position;
      angletype ra,dec;
      double sra,sdec;
      double mura,mudec,smura,smudec;
      double vmag,bvmag;
      char sp[13];
      double par,rv,spar,srv;
      char quality[3];
      double rho[10];
      double absmag;
      angletype l,b;
      double sl,sb;
      vec r,sr,v,sv;
      mat loc_to_eq,eq_to_loc;
      mat eq_to_gal,gal_to_eq;
      mat loc_to_gal,gal_to_loc;
      double jval[37];
      double *jacob[7];
      double sigma[7];
      int sflag;
   }
   startype;

#endif

/* ------------------------ Function Declaration ------------------------ */
int leap_year(int year);
void prep_feb(int year);
int year_length(int year);
long date_order(long d,long m,long y);
double julian_century(double jd);
double get_dt(double year);
long gregorian_date_to_jd(long d,long m,long y);
long julian_date_to_jd(long d,long m,long y);
long date_to_jd(long d,long m,long y);
void jd_to_gregorian_date(long jd,long *d,long *m,long *y);
void jd_to_julian_date(long jd,long *d,long *m,long *y);
void jd_to_date(long jd,long *d,long *m,long *y);
momentype make_moment(double day,double month,double year,
                      double hour,double min,double sec);
void get_precarg(double precarg[3],double t);
void prerot(double precarg[3],mat prec);
double obliquity(double t);
double fundamental_argument(int arg,double t);
double nutarg(int row,double t);
void nutate(double *delpsi,double *deleps,double t);
void nutrot(mat nut,double epsilon,double delpsi,double deleps);
void compute_precnut(epochtype *w);
void compute_gst(epochtype *w);
void compute_lst(topotype *t,epochtype *w);
void galactic_matrices(void);
void compute_galactic(epochtype *w);
epochtype set_ut(double day,double month,double year,
                 double hour,double min,double sec);
epochtype set_tdt(double day,double month,double year,
                  double hour,double min,double sec);
epochtype set_jd(double jd);
epochtype julian_epoch(double year);
double jd_to_julian_year(double jd);
topotype set_topocentre(char l_sgn,double l_hr,double l_min,double l_sec,
                        char f_sgn,double f_dg,double f_min,double f_sec,
                        double height,epochtype *w);
void astro_initialize(void);
void equatorial_matrices(startype *s);
void local_matrices(startype *s);
void fix_star(startype *s);
void unfix_star(startype *s);
startype new_position(startype s,epochtype equinox,int pos,epochtype epoch);
void load_earth(void);
void free_earth(void);
void interpolate_earth(double jd,vec rh,vec vh,vec rb,vec vb);
int hd_to_hip(int hd);
int hr_to_hd(int hr);
int find_bright_star(char *star,char *mask);
int get_hipparcos_number(char *starname);
double rounded(float x);
startype find_star(int hip);
int get_object_info(char *starname,int *hipnum,int *usrnum,startype *star);
void get_barycentric_correction(epochtype epoch,startype star1,topotype topo,
      double *drv,double *djd);
void get_sun_correction(epochtype epoch,topotype topo,double *drv,double *djd);
