#ifndef ASTRO_DEFINED

  #define ASTRO_DEFINED

  #define ICRS    1
  #define MEAN_EQ 2
  #define TRUE_EQ 3

  #define BARY_CENT 1
  #define GEO_CENT  2
  #define TOPO_CENT 3

  #define GEOMETRIC 1
  #define APPARENT  2

  #define PRIMARY_STAR   1
  #define SECONDARY_STAR 2
  #define CENTER_OF_MASS 3

  #define FULL_SPACE_MOTION   1
  #define PROPER_MOTION_ONLY  2

  #define MAX_CONST_COUNT 400
  #define EPHEM_ID_LENGTH 84
  #define CNAMLEN 6
  #define NTARGS 13

  #define FW_GAM 0
  #define FW_PHI 1
  #define FW_PSI 2
  #define FW_EPS 3

  #define MERCURY     0
  #define VENUS       1
  #define EMBARY      2
  #define MARS        3
  #define JUPITER     4
  #define SATURN      5
  #define URANUS      6
  #define NEPTUNE     7
  #define PLUTO       8
  #define MOON        9
  #define SUN        10
  #define NUTATION   11
  #define LIBRATION  12

  #define HIP_STAR      1
  #define HD_STAR       2
  #define HR_STAR       3
  #define FLM_STAR      4
  #define BYR_STAR      5
  #define USR_STAR      6

  #define POSONLY 1
  #define POSVEL  2

  #define LIGHT_SPEED    299792.458

  #define GREEK_LETTER_EXPECTED     1
  #define BAD_GREEK_LETTER          2
  #define CONSTELL_LETTER_EXPECTED  3
  #define BAD_CONSTELL_NAME         4
  #define BAD_FLAMSTEED_NUMBER      5
  #define BAD_BAYER_SEQ             6
  #define BAD_CATALOG_NUMBER        7
  #define BAD_MULTI_COMP            8
  #define MULTI_COMP_END_STR        9
  #define STAR_NAME_ALPHANUM       10

  typedef struct
  {
    char ephem_id[EPHEM_ID_LENGTH+1];
    char ephem_start_epoch[EPHEM_ID_LENGTH+1];
    char ephem_final_epoch[EPHEM_ID_LENGTH+1];
    char cnam[MAX_CONST_COUNT][CNAMLEN+1];
    double cval[MAX_CONST_COUNT];
    int32_t ccount;
    double first_epoch;
    double last_epoch;
    double days_per_record;
    double au,invau;
    double emrat,gamma;
    int32_t pt[NTARGS];
    int32_t ncf[NTARGS];
    int32_t na[NTARGS];
    int32_t numde;
    int32_t ksize;
    int32_t recsize;
  }
  jpl_header;

  typedef struct
  {
    date_type utc;
    datetime_type tai;
    int leap_old,adj,leap_new;
  }
  leap_type;

  typedef struct
  {
    date_type utc;
    double deltat;
  }
  deltat_type;

  typedef struct
  {
    datetime_type tai,tdt,tdb,utc,ut1;
    int leap;
    double deltat;
    double epsilon,delpsi,deleps;
    double xpol,ypol,orig;
    double equorig,erot,equeq;
    alpha_type gmst,gast;
    mat gcrs_to_mean_eq,gcrs_to_true_eq;
    mat mean_eq_to_gcrs,true_eq_to_gcrs;
    mat gcrs_to_earth,earth_to_gcrs;
    mat gcrs_to_cio;
    vec reh,veh,reb,veb;
  }
  epoch_type;

  typedef struct
  {
    double a,i,ome,Ome,e,P,T;
    double b;
    mat r;
  }
  orbit_t;

  typedef struct
  {
    angle_type lon,lat;
    double height;
    epoch_type epoch;
    alpha_type lmst,last;
    vec r,v,re,ve,rb,vb,rh,vh;
  }
  topo_type;

  typedef struct
  {
    int kind;
    int catnum;
    int flam;
    char bayer[4];
    int seq;
    char cons[4];
    char comp[4];
    int hip_ok;
  }
  starlabel_type;

  typedef struct
  {
    int ref,cent,pos;
    epoch_type equinox,epoch;
    topo_type topo;
    angle_type ra,dec;
    double mura,mudec;
    double par,rv;
    double mass_ratio;
    orbit_t orb;
    int component;
    mat xyz_to_loc,loc_to_xyz;
    vec r,v;
    int model;
  }
  star_type;

  typedef struct
  {
    char nom[32];
    char gen[32];
    char abbr[4];
  }
  constellationtype;

  typedef struct
  {
    char full[16];
    char abbr[4];
  }
  greektype;

  typedef struct
  {
    int hr;
    int flam;
    char bayer[4];
    int seq;
    char cons[4];
    int hd;
    char comp[4];
    int valid;
  }
  brighttype;

  typedef struct
  {
    char ccdm[12];
    char comp[4];
    int hd;
    int dblhd;
  }
  ccdmtype;

  typedef struct
  {
    int seq;
    char ccdm[12];
    char comp[4];
    int hip;
    double rv;
  }
  barbiertype;

  typedef struct
  {
    int hip;
    double ra,dec,par,mura,mudec;
    char aref[4];
    char ccdm[12];
    int ncomp;
    char annex[4];
    int hd;
    char comp[4];
    int valid;
  }
  hiptype;

  typedef struct
  {
    starlabel_type prim,seco;
    double epoch;
    double ra,dec,par,mura,mudec,rv;
    double mass_ratio;
    double a,i,ome,Ome,e,P,T;
  }
  binary_type;

#endif

double get_fw_arg(int argseq,double t);
void get_fw_all(double *fw,double t);
void fw_to_matrix(double *fw,mat r);
double cio_sum(int k,double *fa);
double cio_locator(double t,double x,double y);
double earth_rot_angle(datetime_type *ut1);
double gmstp(datetime_type *tdt);
double get_gsd(datetime_type *ut1,datetime_type *tdt);
void deflection_get_apparent(vec p,vec r_hel,vec app);
void deflection_get_true(vec app,vec r_hel,vec p);
void aberration_get_apparent(vec p,vec v_obs,vec app);
void aberration_get_true(vec app,vec v_obs,vec p);
void load_leap(void);
void free_leap(void);
int get_leap_count(void);
void get_leap_item(leap_type *a,int seq);
int get_leap_total_utc(int32_t utc_jd);
int get_leap_adj_utc(int32_t utc_jd);
int get_leap_total_tai(datetime_type *tai);
int leap_seq(datetime_type *tai);
void load_deltat(void);
void free_deltat(void);
int get_deltat_count(void);
void get_deltat_item(deltat_type *a,int seq);
double get_deltat_utc(datetime_type *utc);
void load_jpl(void);
void free_jpl(void);
void jpl_interpolate
  (double *buf,double *t,int ncf,int ncm,int na,int ifl,double *pv);
void jpl_target(datetime_type *tdb,int targ,int fl,double *pv);
void jpl_posvel(datetime_type *tdb,int targ,vec r,vec v);
void jpl_earth(datetime_type *tdb,vec rh,vec vh,vec rb,vec vb);
void load_astro(void);
void free_astro(void);
void clear_epoch(epoch_type *w);
void make_epoch_data(epoch_type *w);
double get_tdb_correction(datetime_type *tdt);
void tai_to_tdt(datetime_type *tai,datetime_type *tdt);
void tdt_to_tai(datetime_type *tdt,datetime_type *tai);
void tdt_to_tdb(datetime_type *tdt,datetime_type *tdb);
void tdb_to_tdt(datetime_type *tdb,datetime_type *tdt);
void utc_to_tai(datetime_type *utc,datetime_type *tai);
void tai_to_utc(datetime_type *tai,datetime_type *utc);
void set_datetime_utc(datetime_type *utc,
  int day,int month,int year,int hour,int min,double sec);
void use_epoch_utc(epoch_type *w);
void use_epoch_tai(epoch_type *w);
void use_epoch_tdt(epoch_type *w);
void use_epoch_tdb(epoch_type *w);
void derive_ut1(epoch_type *w);
void finalize_epoch(epoch_type *w);
void set_epoch_utc(epoch_type *w,
  int day,int month,int year,int hour,int min,double sec);
void set_epoch_tdt(epoch_type *w,
  int day,int month,int year,int hour,int min,double sec);
void set_epoch_tdb(epoch_type *w,
  int day,int month,int year,int hour,int min,double sec);
void new_epoch_ut1(epoch_type *w);
void set_epoch_ut1(epoch_type *w,
  int day,int month,int year,int hour,int min,double sec);
void set_julian_year_tdt(epoch_type *w,double year);
void add_tdt_seconds(epoch_type *a,epoch_type *b,double secs);
void add_tdb_seconds(epoch_type *a,epoch_type *b,double secs);
void add_ut1_seconds(epoch_type *a,epoch_type *b,double secs);
void set_orbit(orbit_t *orb,double a,double i,double ome,double Ome,
  double e,interval_type *P,datetime_type *T);
void get_apparent_position(orbit_t *orb,datetime_type *t,vec r);
void clear_topo(topo_type *t);
void set_topo(topo_type *topo,angle_type *lon,angle_type *lat,
  double height,epoch_type *epoch);
void clear_star(star_type *s);
void local_matrices(star_type *s);
void space_vectors(star_type *s);
void sky_angles(star_type *s);
void cent_bary(star_type *s);
void cent_obs(star_type *s);
void apply_deflection(star_type *s);
void remove_deflection(star_type *s);
void apply_aberration(star_type *s);
void remove_aberration(star_type *s);
void rotate_to_gcrs(star_type *s);
void rotate_to_obs(star_type *s);
void space_motion(star_type *s,double dt);
void locate_components(star_type *s,vec ab,vec ac,vec cb);
void select_primary_star(star_type *s);
void select_secondary_star(star_type *s);
void select_center_of_mass(star_type *s);
void select_component(star_type *s,int comp);
void move_star(star_type *a,star_type *b,int ref,int cent,int pos,
  epoch_type *eq,epoch_type *ep,topo_type *topo);
void full_corr(star_type *astar,star_type *bstar,epoch_type *epoch,
  topo_type *topo,double *vorb,double *vrot,double *dv,double *dt);
void star_corr(star_type *astar,star_type *bstar,epoch_type *epoch,
  topo_type *topo,double *dv,double *dt);
void sun_corr(epoch_type *epoch,topo_type *topo,double *dv,double *dt);
void clear_label(starlabel_type *lab);
void clear_lab_bayer(starlabel_type *lab);
void clear_lab_cons(starlabel_type *lab);
void clear_lab_comp(starlabel_type *lab);
void set_hip_label(starlabel_type *lab,int hip,char *comp);
void set_hr_label(starlabel_type *lab,int hr,char *comp);
void set_usr_label(starlabel_type *lab,int usr);
void compact_star_name(char *instr,char *outstr,int maxlen);
int get_starname_kind(char *starname);
int parse_greek_letter(char *greekname,char **endptr,starlabel_type *lab);
int parse_constellation(char *consname,char **endptr,starlabel_type *lab);
int parse_flamsteed_number(char *flmno,char **endptr,starlabel_type *lab);
int parse_bayer_seq(char *bayerseq,char **endptr,starlabel_type *lab);
int parse_catalog_number(char *catno,char **endptr,starlabel_type *lab);
int parse_bayer_star(char *starname,char **endptr,starlabel_type *lab);
int parse_flamsteed_star(char *starname,char **endptr,starlabel_type *lab);
int parse_catalog_star(char *starname,char **endptr,starlabel_type *lab);
int parse_multiple_component(char *comptr,starlabel_type *lab);
int parse_base_name(char *starname,char **endptr,starlabel_type *lab);
int parse_star_name(char *starname,starlabel_type *lab);
void bad_star_name(char *starname,int parse_code);
void parse_star_name_err(char *starname,starlabel_type *lab);
void load_bright(void);
void free_bright(void);
int bright_loaded(void);
void load_ccdm(void);
void free_ccdm(void);
int ccdm_loaded(void);
void load_barbier(void);
void free_barbier(void);
int barbier_loaded(void);
void load_hip_main(void);
void free_hip_main(void);
void load_hip_com(void);
void free_hip_com(void);
void make_comp_list(void);
void free_comp_list(void);
void load_hip(void);
void free_hip(void);
int hip_loaded(void);
void load_binary(void);
void free_binary(void);
int binary_loaded(void);
void load_cats(void);
void free_cats(void);
char *print_star_label(starlabel_type *lab,char *buf);
char *print_comp_choice(char *comp,char *buf);
void check_bsc_entry(brighttype *fptr,int n,starlabel_type *lab);
void check_hip_entry(hiptype *fptr,int n,starlabel_type *lab);
void bayer_to_hd(starlabel_type *lab);
void flamsteed_to_hd(starlabel_type *lab);
void hr_to_hd(starlabel_type *lab);
void hip_to_hip(starlabel_type *lab);
void hd_to_hip(starlabel_type *lab);
void any_to_hip(starlabel_type *lab);
int same_star(starlabel_type *a,starlabel_type *b);
void load_hip_star(starlabel_type *lab,star_type *star);
void load_usr_star(starlabel_type *lab,star_type *star);
void load_star(starlabel_type *lab,star_type *star);
