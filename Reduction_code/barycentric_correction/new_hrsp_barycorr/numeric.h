#ifndef NUMERIC_DEFINED

   #define NUMERIC_DEFINED

   typedef double (*functype)(double);

   typedef struct
   {
      double base,icen,xcen,ycen,xfwhm,yfwhm,xyfactor;
      double chisq,rms;
      int niter,status;
   }
   gaussian;

#endif

/* ------------------------ Function Declaration ------------------------ */
double ran1(long *idum);
double gasdev(long *idum);
int *ivector(int n);
float *fvector(int n);
double *dvector(int n);
char **cmatrix(int m,int n);
int **imatrix(int m,int n);
float **fmatrix(int m,int n);
double **dmatrix(int m,int n);
void free_ivector(int *v);
void free_fvector(float *v);
void free_dvector(double *v);
void free_fmatrix(float **a);
void free_dmatrix(double **a);
void free_imatrix(int **a);
void free_cmatrix(char **a);
void clear_imatrix(int **a,int m,int n);
void clear_cmatrix(char **a,int m,int n);
void clear_dmatrix(double **a,int m,int n);
void clear_ivector(int *a,int n);
void clear_fvector(float *a,int n);
void clear_dvector(double *a,int n);
void copy_dvector(double *a,int n,double *b);
void dvector_difference(double *a,double *b,double *c,int n);
void load_dmatrix_column(double **a,int m,double *x,int k);
float fselect(int k,int n,float *arr);
double dselect(int k,int n,double *arr);
void shift_fmatrix(float **imag,int nx,int ny,int dx,int dy);
void spline(double *x,double *y,int n,double yp1,double ypn,double *y2);
void spline_value(double *xa,double *ya,double *y2a,int n,
                  int klo,int khi,double x,double *y);
int locate_node(double *x,int n,double u);
int nearest_node(double *x,int n,double u);
void splint(double *xa,double *ya,double *y2a,int n,double x,double *y);
double splineint(double *x,double *y,double *y2,int n,int k,double u,double v);
void integrate_spline(double *x,double *y,double *y2,int n,
                      double *u,int m,double *v);
void spline_derivative(double *x,double *y,double *y2,int n,
                       double *u,int m,double *v);
double golden(double ax,double bx,double cx,double (*f)(double),
              double tol,double *xmin);
double linint(double xa,double ya,double xb,double yb,double x);
double interpolation(double *x,double *y,int n,double u);
double pythag(double a,double b);
void svdcmp(double **a,int m,int n,double *w,double **v,int pmod);
void svedit(double *w,int n);
void svbksb(double **u,double *w,double **v,int m,int n,double *b,double *x);
double intpow(double x,int n);
double compute_polynomial(double *a,int n,double x);
double compute_sub_polynomial(double *a,int n,int g,int s,int m,double x);
void accumulate_coefficients(double *a,int ma,double *b,int mb,
                             double u,int axis);
int number_of_coefficients(int *deg,int naxis);
double polynomial_func(double **x,int dat,int naxis,int *deg,int seq);
int polynomial_matrix(double **x,int ndat,int naxis,double *y,int *sel,
                       int *deg,double **a,int ncols,double *b);
double polynomial_val(double **x,int dat,int naxis,int *deg,double *a,int ma);
void polynomial_fit(double **x,int ndat,int naxis,int *deg,
                    double *a,int ma,double *fit);
double get_chisq(double *del,int *sel,int ndat);
void regression_polynomial(double **x,int ndat,int naxis,double *y,int *sel,
              int *deg,double *a,int ma,double *rms,int *nsel,
              double *fit,double *del,int pmod);
void best_regression_polynomial(double **x,int ndat,int naxis,double *y,
              int *sel,int *deg,double *a,int ma,double *rms, int *nsel,
              double *fit,double *del,double kappa,int pmod);
int polynomial_matrix_1d(double *x,double *y,int ndat,int *sel,
                       int deg,double **a,double *b);
void polynomial_fit_1d(double *x,int ndat,int deg,
                    double *a,double *fit);
void regression_polynomial_1d(double *x,double *y,int ndat,int *sel,
              int deg,double *a,double *rms,int *nsel,
              double *fit,double *del,int pmod);
void best_regression_polynomial_1d(double *x,double *y,int ndat,
              int *sel,int deg,double *a,double *rms, int *nsel,
              double *fit,double *del,double kappa,int pmod);
int gaussj(double **a, int n, double **b, int m);
void covsrt(double **covar, int ma, int *ia, int mfit);
void mrqcof(double *x, double *y, double *sigma, int ndata,
   double *a, int *ia, int ma,
   double **alpha, double *beta, double *chisq,
   void (*funcs)(double, double *, double *, double *, int));
int mrqmin(double *x, double *y, double *sigma, int ndata,
   double *a, int *ia, int ma, 
   double **covar, double **alpha, double *chisq,
   void (*funcs)(double, double *, double *, double *, int), double *alambda);
int nonlinfit(double *x,double *y,int ndata,double *a,int *ia,double *s,
   double *eps,int ma,int nmax,double *chisq,double *rms,int *niter,
   void (*funcs)(double,double *,double *,double *,int));
void fgauss(double x,double *a,double *f,double *dfda,int na);
void fgauss2d(double x,double *a,double *f,double *dfda,int na);
void clear_gaussian(gaussian *g);
void fit_gaussian_1d(double *x,double *f,int n,gaussian *g);
void center_gaussian_1d(float *pix,int npix,double start,double step,
                        int a,int b,double w,gaussian *g);
void locate_gaussian_1d(float *pix,int npix,double start,double step,
                        double a,double c,double b,double w,gaussian *g);
void fit_gaussian_2d(double **xy,double *f,int n,gaussian *g);
void center_gaussian_2d(float *pix,int npix1,int npix2,
                        double start1,double start2,double step1,double step2,
                        int a1,int a2,int b1,int b2,
                        double wx,double wy,gaussian *g);
void locate_gaussian_2d(float *pix,int npix1,int npix2,
                        double start1,double start2,double step1,double step2,
                        double a1,double a2,double c1,double c2,
                        double b1,double b2,double wx,double wy,gaussian *g);
int peak_ok(double *x,int n);
void init_parab(double *y,int n,int *imax,double *xparab,
                double *yparab,int pmod);
void define_core_peak(int imax,double xparab,int rad,int *kg,int *ng,
                      int *b1,int *b2);
void extract_core_peak(double *u,double *v,double *y,int imax,int kg,int ng);
void init_gauss(gaussian *g,double *v,int ng,double xparab,double yparab);
double get_shift(gaussian *g,int imax,double x1,double step);
double locate_gauss_maximum(double x1,double step,double *y,int n,
                            int rad,int pmod,int *b1,int *b2,gaussian *g);
double locate_parab_maximum(double x1,double step,double *y,int n,
                            int rad,int pmod,int *b1,int *b2,gaussian *g);
double locate_spline_maximum(double x1,double step,double *y,int n,
                             double tol,int pmod,gaussian *g);
double locate_peak_maximum(double x1,double step,double *y,int n,
                           int pmod,int *b1,int *b2,gaussian *g);
void four1(double *data,long nn,int isign);
void twofft(double *data1,double *data2,double *fft1,double *fft2,long n);
void realft(double *data,long n,int isign);
void convlv(double *data,long n,double *respns,long m,int isign,double *ans);
void correl(double *data1,double *data2,int n,double *ans);
void fourn(float *data,long *nn,int ndim,int isign,int pmod);
void ccf2d(float *img,float *ref,int n1,int n2,float *ccf);
void rebinning(double *g,int n,double xstart,double xstep,
               double *h,int m,double ustart,double ustep,
               functype func);
void rebin_continuum(double xstep,
                     double *h,int m,double ustart,double ustep,
                     functype func);
void rebin_normal(double *g,int n,double xstart,double xstep,
                  double *h,int m,double ustart,double ustep,
                  functype func);
void correlate_arrays(double *xdat,double *ydat,double *zdat,int n,
                      double xstart,double ystart,double step,double *zstart);
void statistics_double_array(double *a,int n,double *min,double *max,
                             double *mean,double *sig);
void best_statistics_double_array(double *a,int n,double kappa,
                                  double *min,double *max,
                                  double *mean,double *sig,int *nsel);
