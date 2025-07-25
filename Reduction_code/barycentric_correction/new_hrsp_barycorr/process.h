#ifndef PROCESS_DEFINED

   #define PROCESS_DEFINED

   #define MAX_STARNAME 16
   #define MAX_SOURCE   16
   #define MAX_SKYNAME  16
   #define MAX_OBJTYPE  16
   #define MAX_DATE     10
   #define MAX_TIME     11

   #define NO_INFO    0
   #define YES_INFO   1

   typedef struct
   {
      int *seq;
      int *order;
      double *airlam,*vaclam,*waveno;
      char *species;
      double *ustart,*ucen,*uend,*vstart,*vcen,*vend;
      double *lam,*mlam;
      double *unorm,*vnorm,*mnorm,*mlnorm;
      double *vfit,*vresid,*mfit,*mresid,*ufit,*uresid;
      int *vsel,*msel,*usel;
   }
   thar_table;

   typedef struct
   {
      int *seq;
      int *order;
      double *airlam,*vaclam,*waveno;
      char *species;
      double *ustart,*ucen,*uend,*vstart,*vcen,*vend;
      double *xstart,*xcenter,*xend,*ystart,*ycenter,*yend;
      double *base,*icen,*xcen,*ycen,*xfwhm,*yfwhm,*xyfactor;
      double *chisq,*rms;
      int *niter,*status;
      int *xsel,*ysel,*usel,*vsel;
      double *xnorm,*ynorm,*xfit,*xresid,*yfit,*yresid;
      double *unorm,*vnorm,*ufit,*uresid,*vfit,*vresid;
   }
   trans_table;

   typedef struct
   {
      int *seq;
      int *order,*x;
      double *ystart,*ycenter,*yend;
      double *base,*icen,*ycen,*yfwhm;
      double *chisq,*rms;
      int *niter,*status;
      int *ysel,*wsel,*msel;
      double *xnorm,*ynorm,*mnorm;
      double *yfit,*yresid,*wfit,*wresid,*mfit,*mresid;
   }
   ord_table;

   typedef struct
   {
      int *seq;
      int *order;
      double *airlam,*vaclam,*waveno;
      char *species;
      double *ustart,*ucen,*uend,*vstart,*vcen,*vend;
      double *xstart,*xcenter,*xend;
      double *base,*icen,*xcen,*xfwhm;
      double *chisq,*rms;
      int *niter,*status;
      int *xsel;
      double *lam,*mlam,*mnorm,*mlnorm,*xfit,*xresid;
   }
   lin_table;

   typedef struct
   {
      int *seq;
      int *x,*y;
      int *sel;
      double *order,*pixval;
      double *xnorm,*ynorm,*fit,*resid;
   }
   back_table;

   typedef struct
   {
      int *seq;
      int *order,*sel,*radius;
      double *weight,*xstart,*xend;
      double *base,*icen,*xcen,*xfwhm;
      double *chisq,*rms;
      int *niter,*status;
      double *rv_raw,*rv;
   }
   ccf_table;

   typedef struct
   {
      int *seq;
      int *order,*sel,*radius;
      double *weight,*xstart,*xend;
      double *base,*icen,*xcen,*xfwhm;
      double *chisq,*rms;
      int *niter,*status;
      double *shift;
   }
   sft_table;

#endif

void set_thar_columns(thar_table *t,fits_file *f);
void set_trans_columns(trans_table *t,fits_file *f);
void set_ord_columns(ord_table *t,fits_file *f);
void set_lin_columns(lin_table *t,fits_file *f);
void set_back_columns(back_table *t,fits_file *f);
void set_ccf_columns(ccf_table *t,fits_file *f);
void set_sft_columns(sft_table *t,fits_file *f);
double normalized_pixel(fits_file *fits,int axis,double pix);
double normalized_order(spec_type *spec,double ord);
double normalized_u(spec_type *spec,double u);
double normalized_v(spec_type *spec,double v);
double normalized_mlambda(spec_type *spec,double mlam);
double denormalized_mlambda(spec_type *spec,double mlnorm);
void load_dmatrix_from_table(double **a,int k,fits_file *t,int col);
void create_thar_table(int nrows,char *name,...);
void read_thar_table(fits_file *fits,thar_table *table,char *name,...);
void create_trans_table(int nrows,char *name,...);
void read_trans_table(fits_file *fits,trans_table *table,char *name,...);
void create_ord_table(int nrows,char *name,...);
void read_ord_table(fits_file *fits,ord_table *table,char *name,...);
void create_lin_table(int nrows,char *name,...);
void read_lin_table(fits_file *fits,lin_table *table,char *name,...);
void create_back_table(int nrows,char *name,...);
void read_back_table(fits_file *fits,back_table *table,char *name,...);
void create_ccf_table(int nrows,char *name,...);
void read_ccf_table(fits_file *fits,ccf_table *table,char *name,...);
void create_sft_table(int nrows,char *name,...);
void read_sft_table(fits_file *fits,sft_table *table,char *name,...);
void check_pixel_range(int a,int b,int first,int last);
void get_spectrograph(fits_file *xfits,fits_head *head,char *spectrograph);
void get_ccdinfo_fits(fits_file *xfits,fits_head *head,ccd_type *ccd);
void get_ccdinfo_file(char *xname,ccd_type *ccd);
void get_specinfo_fits(fits_file *xfits,fits_head *head,spec_type *spec);
void get_specinfo_file(char *xname,spec_type *spec);
int standard_name_ok(char *filename);
int kiwispec_name_ok(char *filename);
void check_standard_name(char *filename);
void check_kiwispec_name(char *filename);
int check_file_type(char *filename);
void parse_standard_name(char *filename,int *mjd,int *imgno);
void parse_kiwispec_name(char *filename,int *mjd,int *imgno);
void kiwispec_object(char *filename,char *object);
void kiwispec_exposure_type(char *object,char *exptype);
void fix_kiwispec_header(fits_file *fits,char *filename);
void extract_image_proc(char *xname,char *yname,int ax,int ay,int bx,int by,
                        int info);
int rotation_argument(int angle);
int flip_argument(int id);
void swap_axes(fits_file *fits);
void rotate_pixels(fits_file *xfits,fits_file *yfits,
                      int x1,int dx,int y1,int dy,int ry);
void rotate_fits_image(char *xname,char *yname,int q);
void rotate_image_proc(char *xname,char *yname,int angle,int pmod);
void flip_image_proc(char *xname,char *yname,int axis,int pmod);
void quad_align_proc(char *xname,char *yname);
void convert_image_proc(char *xname,char *yname);
void subtract_bias(char *name);
void fft_fits(fits_file *xfits,fits_file *rfits,fits_file *ifits);
void ifft_fits(fits_file *xfits,fits_file *rfits,fits_file *ifits);
void complex_product_fits(fits_file *xr,fits_file *xi,
   fits_file *yr,fits_file *yi,fits_file *zr,fits_file *zi,int conj);
void correlate_2d_proc(char *xname,char *rname,char *cname);
void regression_polynomial_proc(char *tabname,char *yname,char *xname,
        char *deglist,char *selname,char *fitname,char *residname,
        char *rkey,double kappa);
int collect_jd_number(char *s);
int collect_image_number(char *s);
int keydef_list_count(void);
void clear_keydef_list(void);
void populate_values(keydef_type *k,char *val);
void read_keydef_cfg(fits_file *fits);
void get_keydef_item(int seq,keydef_type *key);
int get_fibre_number(fits_file *fits);
int count_thorium_ref(char *spectrograph);
void get_dispersion
  (spec_type *spec,regre_block *ureg,int m,double *lamcen,double *disp);
void get_dispersion_list(spec_type *spec,double *lamcen,double *disp);
double true_dispersion(regre_block *xreg,double mlnorm,double mnorm);
double log_doppler(double vel);
double doppler_vel(double logdopp);
double redshift(double vel);
int create_file_list(dir_name_type dirname, file_name_type filename);
