#ifndef CONFIG_DEFINED

   #define CONFIG_DEFINED

   #define MAXCFGLINE 120

   #define MAX_NVALS 10
   #define MAX_VALSIZE 32
   #define MAX_KEYSIZE 32

   #define MAX_REGIONS 8
   #define MAX_FIBRES  8

   #define MAX_CCDNAME 16
   #define MAXBADLINES 20
   #define MAXGAINS 8

   #define MAX_DEGARG  8
   #define MAX_REGRE   2

   #define MAXBLOCKS 16

   #define MAXORD 200
   #define MAXGR 16

   typedef char keytype[MAX_KEYSIZE+1];
   typedef char valtype[MAX_NVALS+1][MAX_VALSIZE+1];

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
      double usize,vsize,uscale,vscale;
      int firstord,lastord,nord;
      double mscale;
      double mlamstart,mlamend,mlamrange,mlamscale;
      double refair,refvac;
      int reford;
      double uref,vref;
      int nreg;
      fbox_type reg[MAX_REGIONS+1];
      int nfib;
      fibre_type fib[MAX_FIBRES+1];
      regre_type echdef_vcen,echdef_ord;
   }
   herc_type;

   typedef struct
   {
      int ax,ay,bx,by;
   }
   box_type;

   typedef struct
   {
      int gainset;
      double invgain,bnoise;
   }
   gain_type;

   typedef struct
   {
      char name[MAX_CCDNAME+1];
      int chipsize[3];
      double pixsize[3];
      double maxpix;
      int nbadlin;
      box_type badlin[MAXBADLINES+1];
      int ngain;
      gain_type gain[MAXGAINS+1];
      int rotate;
      char flip[3];
   }
   ccd_type;

   typedef struct
   {
      int reg,fib;
      int trans_rad[3];
      regre_type trans;
      int echdef_ncols;
      regre_type echdef_ycen,echdef_fwhm,echdef_ord;
      int echback_ncols;
      regre_type echback;
      int medfilt_rx,medfilt_ry;
      int crfilt_nboxes,crfilt_rslit;
      double crfilt_nsigma;
      int crfilt_npasses;
      regre_type crfilt;
      int flat_rfilt;
      regre_type flat;
      double extract_slit;
      double mlamstart,mlamend,mlamrange;
      int nbins;
   }
   param_type;

   typedef struct
   {
      int firstord,lastord;
      double mlamstart,mlamend;
      char wave[8];
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
      int nblock;
      block_type block[MAXBLOCKS+1];
      int nomit;
      block_type omit[MAXBLOCKS+1];
      int rvok[MAXORD];
      double mlamstart,mlamend;
      int cosbell;
      char method[16];
      char weights[80];
      double w[MAXORD];
      char radius[80];
      int r[MAXORD];
   }
   radvel_type;   

#endif

/* ------------------------ Function Declaration ------------------------ */
int white_space_ok(char c);
int quot_ok(char c);
int end_of_line_ok(char c);
int keyname_char_ok(char c);
int keyval_char_ok(char c);
char get_char(void);
char next_char(void);
void scan_config_line(char *s,keytype key,int *n,valtype val);
void get_next_keyword(file_type *f,keytype key,int *n,valtype val);
void expect_next_keyword(file_type *f,char *s,int m,keytype key,valtype val);
int expect_next_keyword_list(file_type *f,key_info *s,keytype key,valtype val);
void get_string_value(char *s,int maxlen,valtype val,int k);
int get_integer_value(valtype val,int k);
void get_integer_value_array(int *a,valtype val,int k,int n);
void get_box_type_value(box_type *box,valtype val);
double get_double_value(valtype val,int k);
void get_double_value_array(double *a,valtype val,int k,int n);
void expect_next_string_value(file_type *f,char *s,char *x,int maxlen);
int expect_next_integer_value(file_type *f,char *s);
double expect_next_double_value(file_type *f,char *s);
void open_config_file(file_type *file,char *dir,char *name,...);
void open_keydef_cfg(file_type *file);
void get_regre_value(regre_type *reg,valtype val,int k,int n);
void get_herc_info(herc_type *herc);
void load_ccdinfo(char *ccdname,ccd_type *ccd);
void get_parameters(param_type *par,int reg,int fib);
void get_dispsol_info(dispsol_type *dsol,int id);
void get_radvel_info(radvel_type *rvel,int id);
