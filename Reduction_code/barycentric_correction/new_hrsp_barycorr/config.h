#ifndef CONFIG_DEFINED

   #define CONFIG_DEFINED

   #define MAXCFGLINE 120
   #define MAXCFGROWS 1000

   #define MAX_NVALS 10
   #define MAX_VALSIZE 32
   #define MAX_KEYSIZE 32

   #define MAX_FIBRES  8

   #define MAX_CCDNAME 16
   #define MAX_SPECNAME 16
   #define MAX_OBSNAME 16
   #define MAXGAINS 8

   #define MAX_DEGARG  8
   #define MAX_REGRE   2

   #define MAXORD 200
   #define MAXGR 16

   #define TYPESTR 1
   #define TYPEINT 2
   #define TYPEDBL 3
   #define TYPELOG 4

   #define AIR_WAVE 1
   #define VAC_WAVE 2

   #define FLAT_NONE      1
   #define FLAT_EXTRACTED 2
   #define FLAT_FULL      3

   #define BACK_CUBIC     1
   #define BACK_SPLINE    2
   #define BACK_POLYFIT   3

   #define PREP_NONE   1
   #define PREP_FLIP   2
   #define PREP_MEAN   3

   #define METHOD_GAUSS    1
   #define METHOD_PARABOLA 2
   #define METHOD_SPLINE   3
   #define METHOD_PEAK     4
   #define METHOD_CCFMAX   5

   typedef char cfglinetype[MAXCFGLINE+1];
   typedef char keytype[MAX_KEYSIZE+1];
   typedef char valtype[MAX_NVALS+1][MAX_VALSIZE+1];

   typedef struct
   {
     int id;
     char txt[MAX_VALSIZE+1];
   }
   txtpar;

   typedef struct
   {
      keytype name;
      int nval;
   }
   key_info;

   typedef struct
   {
      double ax,ay,bx,by;
   }
   fbox_type;

   typedef struct
   {
      double dx,dy;
   }
   fibre_type;

   typedef struct
   {
      int naxis;
      int deg[MAX_REGRE+1];
      char degarg[MAX_DEGARG+1];
      double kappa;
   }
   regre_type;

   typedef struct
   {
      char specname[MAX_SPECNAME+1];
      char fullname[MAX_SPECNAME+1];
      double usize,vsize,uscale,vscale;
      int firstord,lastord,nord;
      int bincount;
      double binsize;
      double mscale;
      double mlamstart,mlamend,mlamrange,mlamscale;
      double refair,refvac;
      int reford;
      double uref,vref;
      int nfib;
      fibre_type fib[MAX_FIBRES+1];
      regre_type echdef_vcen,echdef_ord,echdef_ucen;
   }
   spec_type;

   typedef struct
   {
      int gainset;
      double invgain,bnoise;
   }
   gain_type;

   typedef struct
   {
      char ccdname[MAX_CCDNAME+1];
      char fullname[MAX_CCDNAME+1];
      int chipsize[3];
      double pixsize[3];
      double maxpix;
      int ngain;
      gain_type gain[MAXGAINS+1];
      int rotate;
      txtpar flip;
   }
   ccd_type;

   typedef struct
   {
      int trans_rad[3];
      regre_type trans;
      int echdef_ncols;
      regre_type echdef_ycen,echdef_fwhm,echdef_ord;
      txtpar back_method;
      int echback_ncols;
      regre_type echback;
      int medfilt_rx,medfilt_ry;
      int crfilt_nboxes,crfilt_rslit;
      double crfilt_nsigma;
      int crfilt_npasses;
      regre_type crfilt;
      int flat_rfilt;
      regre_type flat;
      int flat_smooth;
      double flat_thres_full,flat_thres_ext;
      txtpar flat_method;
      double extract_slit;
   }
   param_type;

   typedef struct
   {
      int firstord,lastord;
      double mlamstart,mlamend;
      txtpar wave;
      regre_type reg;
   }
   dispsol_type;

   typedef struct
   {
      int first,last;
   }
   block_type;

   typedef struct
   {
      int sel[MAXORD+1];
      double mlamstart,mlamend;
      txtpar ordprep;
      int cosbell;
      txtpar method;
      int checkrad;
      double w[MAXORD+1];
      int r[MAXORD+1];
   }
   radvel_type;

   typedef struct
   {
      cfglinetype item;
      int type;
      int found;
      cfglinetype sval;
      int ival;
      double dval;
      char lval;
   }
   keydef_type;

   typedef struct
   {
      char name[MAX_OBSNAME+1];
      angle_type lon,lat;
      double height;
   }
   obs_type;

   typedef struct
   {
     int left_margin;
     int right_margin;
     int top_margin;
     int bottom_margin;
     int numrows,numcols;
     int rowgap,colgap;
     int linewidth,fontsize,ticksize;
     int minsep[3];
     int global_limits,window_limits;
     int plots_per_page;
     int boxw,boxh;
     int fontcent;
   }
   plt_type;

   typedef struct
   {
      int sel[MAXORD+1];
      int xleft,xright;
      txtpar ordprep;
      int cosbell;
      txtpar method;
      int checkrad;
      double w[MAXORD+1];
      int r[MAXORD+1];
   }
   shift_type;

#endif

int cfg_white_space_ok(char c);
int cfg_quot_ok(char c);
int cfg_end_of_line_ok(char c);
int keyname_char_ok(char c);
int keyval_char_ok(char c);
char get_char(void);
char next_char(void);
void scan_config_line(char *s,keytype key,int *n,valtype val);
void get_next_keyword(file_type *f,keytype key,int *n,valtype val);
void expect_next_keyword(file_type *f,char *s,int m,keytype key,valtype val);
void get_string_value(char *s,int maxlen,valtype val,int k);
int get_integer_value(valtype val,int k);
void get_integer_value_array(int *a,valtype val,int k,int n);
double get_double_value(valtype val,int k);
void get_double_value_array(double *a,valtype val,int k,int n);
void expect_next_string_value(file_type *f,char *s,char *x,int maxlen);
int expect_next_integer_value(file_type *f,char *s);
double expect_next_double_value(file_type *f,char *s);
void search_text_arg(char *arg,char list[][MAX_VALSIZE+1],int *seq,char *name);
void expect_next_text_param
  (file_type *f,char *s,txtpar *p,char list[][MAX_VALSIZE+1]);
int expect_yes_or_no(file_type *file,char *key);
void open_config_file(file_type *file,char *spectrograph,char *name,...);
void open_keydef_cfg(file_type *file,char *spectrograph);
void get_regre_value(regre_type *reg,valtype val,int k,int n);
void load_spec_info(char *specname,spec_type *spec);
void load_ccd_info(char *ccdname,char *spectrograph,ccd_type *ccd);
void get_observatory(obs_type *obs,char *spectrograph);
int micro(double mm);
void get_plot_info(plt_type *plt,char *spectrograph);
void get_parameters(param_type *par,char *spectrograph);
void get_dispsol_info(dispsol_type *dsol,char *spectrograph);
void get_radvel_info(radvel_type *rvcfg,char *spectrograph);
void get_shift_info(shift_type *sft,char *spectrograph);
void read_lsout(char * str);
int matched_ok(char *f,char *t);
void read_matching_cfg
  (char *cfgname,char *fname,char **keys,char **vals,int *n);
char **allocate_cfglines(void);
void free_cfglines(char **list);
