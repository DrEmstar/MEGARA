#ifndef PROCESS_DEFINED

   #define PROCESS_DEFINED

   #define MAX_STARNAME 16
   #define MAX_SOURCE   16
   #define MAX_SKYNAME  16
   #define MAX_DATE     10
   #define MAX_TIME     10
   #define MAX_HERCFIB  64
   #define MAX_IMGMAT   20

   #define NO_INFO    0
   #define YES_INFO   1

   typedef struct
   {
      int *seq;
      int *order;
      double *airlam,*vaclam,*waveno;
      char *species;
      double *ustart,*ucen,*uend,*vstart,*vcen,*vend;
      double *unorm,*vnorm,*mnorm,*vfit,*vresid,*mfit,*mresid;
      int *vsel,*msel;
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

#endif

/* ------------------------ Function Declaration ------------------------ */
void set_thar_columns(thar_table *t,fits_file *f);
void set_trans_columns(trans_table *t,fits_file *f);
void set_ord_columns(ord_table *t,fits_file *f);
void set_lin_columns(lin_table *t,fits_file *f);
void set_back_columns(back_table *t,fits_file *f);
void set_ccf_columns(ccf_table *t,fits_file *f);
double normalized_pixel(fits_file *fits,int axis,double pix);
double normalized_order(herc_type *herc,double ord);
double normalized_u(herc_type *herc,double u);
double normalized_v(herc_type *herc,double v);
double normalized_mlambda(herc_type *herc,double mlam);
void read_image_parameters_fits(param_type *par,fits_file *fits);
void read_image_parameters_file(param_type *par,char *name);
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
void check_pixel_range(int a,int b,int first,int last);
void get_ccdinfo_fits(fits_file *xfits,ccd_type *ccd);
void get_ccdinfo_file(char *xname,ccd_type *ccd);
void extract_image_proc(char *xname,char *yname,int ax,int ay,int bx,int by,
                        int info);
int rotation_argument(int angle);
int flip_argument(char *s);
void swap_axes(fits_file *fits);
void rotate_pixels(fits_file *xfits,fits_file *yfits,
                      int x1,int dx,int y1,int dy,int ry);
void multiply_imgmat(char *amat,char *bmat,char *cmat);
void rotate_fits_image(char *xname,char *yname,int q);
void rotate_image_proc(char *xname,char *yname,int angle,int pmod);
void flip_image_proc(char *xname,char *yname,char *axis,int pmod);
void convert_image_proc(char *xname,char *yname);
void subtract_bias(char *name);
void correlate_2d_proc(char *xname,char *yname,char *zname);
void regression_polynomial_proc(char *tabname,char *yname,char *xname,
        char *deglist,char *selname,char *fitname,char *residname,
        char *rkey,double kappa);
int collect_jd_number(char *s);
int collect_image_number(char *s);
