#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "vector.h"
#include "numeric.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define R 0.61803399
#define C (1.0-R)
#define SHFT2(a,b,c) {(a)=(b);(b)=(c);}
#define SHFT3(a,b,c,d) {(a)=(b);(b)=(c);(c)=(d);}

#define FWHM_TOLLERANCE 4.0

#define FWHM_FACTOR 1.6651092223153955  /* k = 2*sqrt(log(2)) = FWHM / S */

#define GFIT_SCALE_MIN   1.0e-6

#define GFIT_RMS_MAX     0.2

#define G1DFIT_EPS       1.0e-8
#define G2DFIT_EPS       1.0e-8

#define G1DFIT_MAXPASS   50
#define G2DFIT_MAXPASS   50

#define GAUSSJ_SINGMAT_1     1
#define GAUSSJ_SINGMAT_2     2

#define MRQMIN_FLAG_OFFSET  10

#define NONLINFIT_NO_CONVERGENCE  21

#define G1DFIT_SMALL_NDAT    1111
#define G1DFIT_SMALL_XSCALE  1121
#define G1DFIT_SMALL_FSCALE  1131
#define G1DFIT_LARGE_RMS     1141
#define G1DFIT_BAD_XFWHM     1151
#define G1DFIT_BAD_XCEN      1161

#define G1DCENT_BAD_XREGION  1211
#define G1DCENT_XREGION_OUT  1221
#define G1DCENT_SMALL_NDAT   1231
#define G1DCENT_ZERO_FRANGE  1241
#define G1DCENT_BAD_XFWHM    1251

#define G1DLOC_BAD_XREGION   1311
#define G1DLOC_XREGION_OUT   1321
#define G1DLOC_BAD_XMAX      1331
#define G1DLOC_BAD_MAX       1341

#define G2DFIT_SMALL_NDAT    2111
#define G2DFIT_SMALL_XSCALE  2121
#define G2DFIT_SMALL_YSCALE  2122
#define G2DFIT_SMALL_FSCALE  2131
#define G2DFIT_LARGE_RMS     2141
#define G2DFIT_BAD_XFWHM     2151
#define G2DFIT_BAD_YFWHM     2152
#define G2DFIT_BAD_XCEN      2161
#define G2DFIT_BAD_YCEN      2162

#define G2DCENT_BAD_XREGION  2211
#define G2DCENT_BAD_YREGION  2212
#define G2DCENT_XREGION_OUT  2221
#define G2DCENT_YREGION_OUT  2222
#define G2DCENT_SMALL_NXDAT  2231
#define G2DCENT_SMALL_NYDAT  2232
#define G2DCENT_ZERO_FRANGE  2241
#define G2DCENT_BAD_XFWHM    2251
#define G2DCENT_BAD_YFWHM    2252

#define G2DLOC_BAD_XREGION   2311
#define G2DLOC_BAD_YREGION   2312
#define G2DLOC_XREGION_OUT   2321
#define G2DLOC_YREGION_OUT   2322
#define G2DLOC_BAD_XMAX      2331
#define G2DLOC_BAD_YMAX      2332
#define G2DLOC_BAD_MAX       2341

#define CLEAN_GAUSSJ {\
   free_ivector(ipiv);\
   free_ivector(indxr);\
   free_ivector(indxc);}

#define CLEAN_MRQMIN {\
   free_dmatrix(oneda);\
   free_dvector(da);\
   free_dvector(beta);\
   free_dvector(atry);}

#define CLEAN_NONLINFIT {\
   free_dmatrix(covar);\
   free_dmatrix(alpha);\
   free_dvector(b);\
   free_dvector(sigma);}

#define CLEAN_CENTGAUSS1D {\
   free_dvector(x);\
   free_dvector(f);}

#define CLEAN_CENTGAUSS2D {\
   free_dmatrix(xy);\
   free_dvector(f);}

static float fswaptemp;
#define  FSWAP(a,b) {fswaptemp=(a);(a)=(b);(b)=fswaptemp;}

static double dswaptemp;
#define  DSWAP(a,b) {dswaptemp=(a);(a)=(b);(b)=dswaptemp;}

static double dsqrarg;
#define DSQR(x) ((dsqrarg=x) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static int iminarg1,iminarg2;
#define IMIN(a,b) ((iminarg2=b) < (iminarg1=a) ? iminarg2 : iminarg1)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) ((dmaxarg2=b) > (dmaxarg1=a) ? dmaxarg2 : dmaxarg1)

static double dminarg1,dminarg2;
#define DMIN(a,b) ((dminarg2=b) < (dminarg1=a) ? dminarg2 : dminarg1)

#define SIGN(a,b) (b < 0 ? -fabs(a) : fabs(a))

#define ODD(k) (k & 1)

static double **xypoint;

static double *x_vec = NULL;
static double *y_vec = NULL;
static double *y2_vec = NULL;
static int y_count = 0;

double ran1(int32_t *idum)
{
   int j;
   int32_t k;
   static int32_t iy=0L;
   static int32_t iv[NTAB];
   double random;

   if (*idum<=0 || !iy)
   {
      if (-(*idum) < 1)
         *idum = 1;
      else
         *idum = -(*idum);
      for (j=NTAB+7;j>=0;j--)
      {
         k=(*idum)/IQ;
         *idum=IA*(*idum-k*IQ)-IR*k;
         if (*idum < 0) *idum += IM;
         if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if (*idum < 0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j]=*idum;
   random=AM*iy;
   if (random > RNMX)
      return(RNMX);
   else
      return(random);
}

double gasdev(int32_t *idum)
{
   static int iset=0;
   static double gset;
   double fac,rsq,v1,v2;

   if (iset == 0)
   {
      do
      {
         v1=2.0*ran1(idum)-1.0;
         v2=2.0*ran1(idum)-1.0;
         rsq=v1*v1+v2*v2;
      }
      while(rsq>=1.0 || rsq==0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return(v2*fac);
   }
   else
   {
      iset=0;
      return(gset);
   }
}

int *ivector(int n)
{
   int *v;

   v = (int *)get_space((n+1)*sizeof(int));
   return(v);
}

float *fvector(int n)
{
   float *v;

   v = (float *)get_space((n+1)*sizeof(float));
   return(v);
}

double *dvector(int n)
{
   double *v;

   v = (double *)get_space((n+1)*sizeof(double));
   return(v);
}

char **cmatrix(int m,int n)
{
   int i;
   char **a;

   a = (char **)get_space((m+1)*sizeof(char *));
   a[1] = (char *)get_space((m*n+1)*sizeof(char));
   for (i=2;i<=m;i++) a[i]=a[i-1]+n;
   return(a);
}

int **imatrix(int m,int n)
{
   int i;
   int **a;

   a = (int **)get_space((m+1)*sizeof(int *));
   a[1] = (int *)get_space((m*n+1)*sizeof(int));
   for (i=2;i<=m;i++) a[i]=a[i-1]+n;
   return(a);
}

float **fmatrix(int m,int n)
{
   int i;
   float **a;

   a = (float **)get_space((m+1)*sizeof(float *));
   a[1] = (float *)get_space((m*n+1)*sizeof(float));
   for (i=2;i<=m;i++) a[i]=a[i-1]+n;
   return(a);
}

double **dmatrix(int m,int n)
{
   int i;
   double **a;

   a = (double **)get_space((m+1)*sizeof(double *));
   a[1] = (double *)get_space((m*n+1)*sizeof(double));
   for (i=2;i<=m;i++) a[i]=a[i-1]+n;
   return(a);
}

void free_ivector(int *v)
{
   get_free((void *)v);
}

void free_fvector(float *v)
{
   get_free((void *)v);
}

void free_dvector(double *v)
{
   get_free((void *)v);
}

void free_fmatrix(float **a)
{
   get_free((void *)a[1]);
   get_free((void *)a);
}

void free_dmatrix(double **a)
{
   get_free((void *)a[1]);
   get_free((void *)a);
}

void free_imatrix(int **a)
{
   get_free((void *)a[1]);
   get_free((void *)a);
}

void free_cmatrix(char **a)
{
   get_free((void *)a[1]);
   get_free((void *)a);
}

void clear_imatrix(int **a,int m,int n)
{
   memset(&a[1][1],0,m*n*sizeof(int));
}

void clear_cmatrix(char **a,int m,int n)
{
   memset(&a[1][1],0,m*n*sizeof(char));
}

void clear_dmatrix(double **a,int m,int n)
{
   memset(&a[1][1],0,m*n*sizeof(double));
}

void clear_ivector(int *a,int n)
{
   memset(&a[1],0,n*sizeof(int));
}

void clear_fvector(float *a,int n)
{
   memset(&a[1],0,n*sizeof(float));
}

void clear_dvector(double *a,int n)
{
   memset(&a[1],0,n*sizeof(double));
}

void copy_dvector(double *a,int n,double *b)
{
   memmove(&b[1],&a[1],n*sizeof(double));
}

void dvector_difference(double *a,double *b,double *c,int n)
{
   int i;

   for (i=1;i<=n;i++) c[i] = a[i] - b[i];
}

void load_dmatrix_column(double **a,int m,double *x,int k)
{
   int i;

   for (i=1;i<=m;i++) a[i][k] = x[i];
}

float fselect(int k,int n,float *arr)
{
   int i,ir,j,l,mid;
   float a;

   l = 1;
   ir = n;
   for (;;)
   {
      if (ir <= l+1)
      {
         if (ir == l+1 && arr[ir] < arr[l]) FSWAP(arr[l],arr[ir]);
         return(arr[k]);
      }
      else
      {
         mid = (l+ir) >> 1;
         FSWAP(arr[mid],arr[l+1]);
         if (arr[l+1] > arr[ir]) FSWAP(arr[l+1],arr[ir]);
         if (arr[l] > arr[ir]) FSWAP(arr[l],arr[ir]);
         if (arr[l+1] > arr[l]) FSWAP(arr[l+1],arr[l]);
         i = l+1;
         j = ir;
         a = arr[l];
         for (;;)
         {
            do i++; while (arr[i] < a);
            do j--; while (arr[j] > a);
            if (j < i) break;
            FSWAP(arr[i],arr[j]);
         }
         arr[l] = arr[j];
         arr[j] = a;
         if (j >= k) ir = j-1;
         if (j <= k) l=i;
      }
   }
}

double dselect(int k,int n,double *arr)
{
   int i,ir,j,l,mid;
   double a;

   l = 1;
   ir = n;
   for (;;)
   {
      if (ir <= l+1)
      {
         if (ir == l+1 && arr[ir] < arr[l]) DSWAP(arr[l],arr[ir]);
         return(arr[k]);
      }
      else
      {
         mid = (l+ir) >> 1;
         DSWAP(arr[mid],arr[l+1]);
         if (arr[l+1] > arr[ir]) DSWAP(arr[l+1],arr[ir]);
         if (arr[l] > arr[ir]) DSWAP(arr[l],arr[ir]);
         if (arr[l+1] > arr[l]) DSWAP(arr[l+1],arr[l]);
         i = l+1;
         j = ir;
         a = arr[l];
         for (;;)
         {
            do i++; while (arr[i] < a);
            do j--; while (arr[j] > a);
            if (j < i) break;
            DSWAP(arr[i],arr[j]);
         }
         arr[l] = arr[j];
         arr[j] = a;
         if (j >= k) ir = j-1;
         if (j <= k) l=i;
      }
   }
}

double get_median(double *a,int n)
{
  if ((n/2)*2 == n)
    return((dselect(n/2,n,a) + dselect(n/2+1,n,a)) / 2.0);
  else
    return(dselect((n+1)/2,n,a));
}

void shift_fmatrix(float **imag,int nx,int ny,int dx,int dy)
{
   float **buf;
   div_t q;
   int m,sx,sy,ax,ay,px,py,gx,gy,iy,ky;

   q = div(dx,nx);
   sx = q.rem;
   if (sx < 0) sx += nx;
   ax = nx-sx;
   px = sx + 1;
   gx = ax + 1;

   q = div(dy,ny);
   sy = q.rem;
   if (sy < 0) sy += ny;
   ay = ny-sy;
   py = sy + 1;
   gy = ay + 1;

   if (sx==0 && sy==0) return;

   buf = fmatrix(ny,nx);

   m = sizeof(float);

   for (iy=1,ky=py;ky<=ny;iy++,ky++) memmove(&buf[ky][px],&imag[iy][1],ax*m);
   for (iy=gy,ky=1;iy<=ny;iy++,ky++) memmove(&buf[ky][px],&imag[iy][1],ax*m);
   if (sx > 0)
   {
     for (iy=1,ky=py;ky<=ny;iy++,ky++) memmove(&buf[ky][1],&imag[iy][gx],sx*m);
     for (iy=gy,ky=1;iy<=ny;iy++,ky++) memmove(&buf[ky][1],&imag[iy][gx],sx*m);
   }

   memmove(&imag[1][1],&buf[1][1],nx*ny*m);

   free_fmatrix(buf);
}

void spline(double *x,double *y,int n,double yp1,double ypn,double *y2)
{
   int i,k;
   double p,qn,sig,un,*u;

   for (i=1;i<n;i++)
      if (x[i]>=x[i+1]) get_error("spline: X values not increasing!");

   u = dvector(n-1);
   if (yp1 > 0.99e30)
      y2[1]=u[1]=0.0;
   else
      y2[1]=-0.5, u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
   for (i=2;i<=n-1;i++)
   {
      sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p = sig*y2[i-1]+2.0;
      y2[i] = (sig-1.0)/p;
      u[i] = (y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
   }
   if (ypn > 0.99e30)
      qn=un=0.0;
   else
      qn=0.5, un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
   y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0);
   for (k=n-1;k>=1;k--) y2[k] = y2[k]*y2[k+1]+u[k];
   free_dvector(u);
}

void spline_value(double *xa,double *ya,double *y2a,int n,
                  int klo,int khi,double x,double *y)
{
   double h,a,b;

   h = xa[khi]-xa[klo];
   if (h == 0.0) get_error("spline_value: Bad xa input!");
   a = (xa[khi]-x)/h;
   b = (x-xa[klo])/h;
   if (a>=0.0 && b>=0) /* Interpolation */
   {
     *y = a*ya[klo] + b*ya[khi]
        + ((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
   }
   else /* Extrapolation (linear)*/
   {
     *y = ya[klo] + (ya[khi] - ya[klo]) * b;
   }
}

int locate_node(double *x,int n,double u)
{
   int klo,khi,k;

   klo=1,khi=n;
   while (khi-klo > 1)
   {
      k = (khi+klo)>>1;
      if (x[k] > u)
         khi=k;
      else
         klo=k;
   }

   return(klo);
}

int nearest_node(double *x,int n,double u)
{
   int k;

   if (n == 1) return(1);
   k = locate_node(x,n,u);
   if (u-x[k] < x[k+1]-u)
      return(k);
   else
      return(k+1);
}

void splint(double *xa,double *ya,double *y2a,int n,double x,double *y)
{
   int k;

   k = locate_node(xa,n,x);
   spline_value(xa,ya,y2a,n,k,k+1,x,y);
}

double f_norm(double t,double b,double h)
{
   return((b - t) / h);
}

double g_norm(double t,double a,double h)
{
   return((t - a) / h);
}

double splineint(double *x,double *y,double *y2,int n,int k,double u,double v)
{
   double a,b,h,fu,fv,fu2,fv2,fu4,fv4,gu,gv,gu2,gv2,gu4,gv4;
   double fint,gint,fint3,gint3,s;

   a = x[k];
   b = x[k+1];
   h = b - a;

   if (v <= a || u >= b) return(0.0);
   if (u < a) u = a;
   if (v > b) v = b;

   fu = f_norm(u,b,h);
   fv = f_norm(v,b,h);

   fu2 = fu*fu;
   fv2 = fv*fv;

   fu4 = fu2*fu2;
   fv4 = fv2*fv2;

   gu = g_norm(u,a,h);
   gv = g_norm(v,a,h);

   gu2 = gu*gu;
   gv2 = gv*gv;

   gu4 = gu2*gu2;
   gv4 = gv2*gv2;

   fint = (fu2-fv2)*h/2;
   gint = (gv2-gu2)*h/2;

   fint3 = (fu4-fv4)*h/4;
   gint3 = (gv4-gu4)*h/4;

   s = fint*y[k] + gint*y[k+1] +
       ((fint3-fint)*y2[k] + (gint3-gint)*y2[k+1])*h*h/6;

   return(s);
}

void integrate_spline(double *x,double *y,double *y2,int n,
                      double *u,int m,double *v)
{
   int g,k,j,i;

   g = locate_node(x,n,u[1]);

   for (j=1;j<m;j++)
   {
      if (u[j+1] < u[j]) get_error("integrate_spline: Negative interval!");
      k = g;
      g =  locate_node(x,n,u[j+1]);
      if (k == g)
      {
         v[j] = splineint(x,y,y2,n,k,u[j],u[j+1]);
      }
      else
      {
         v[j] = splineint(x,y,y2,n,k,u[j],x[k+1]);
         for (i=k+1;i<g;i++) v[j] += splineint(x,y,y2,n,i,x[i],x[i+1]);
         v[j] += splineint(x,y,y2,n,g,x[g],u[j+1]);
      }
   }
}

void spline_derivative(double *x,double *y,double *y2,int n,
                       double *u,int m,double *v)
{
   int i,k;
   double a,b,h,ff,gg;

   for (i=1;i<=m;i++)
   {
      k =  locate_node(x,n,u[i]);
      a = x[k];
      b = x[k+1];
      h = b-a;
      ff = DSQR(f_norm(u[i],b,h));
      gg = DSQR(g_norm(u[i],a,h));
      v[i] = (y[k+1]-y[k])/h + ((3*gg-1)*y2[k+1]-(3*ff-1)*y2[k])*h/6;
   }
}

double golden(double ax,double bx,double cx,double (*func)(double),
              double tol,double *xmin)
{
   double f1,f2,x0,x1,x2,x3;

   x0 = ax;
   x3 = cx;
   if (fabs(cx-bx) > fabs(bx-ax))
   {
      x1 = bx;
      x2 = bx+C*(cx-bx);
   }
   else
   {
      x2 = bx;
      x1 = bx-C*(bx-ax);
   }
   f1 = func(x1);
   f2 = func(x2);
   while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2)))
   {
      if (f2 < f1)
      {
         SHFT3(x0,x1,x2,R*x1+C*x3);
         SHFT2(f1,f2,func(x2));
      }
      else
      {
         SHFT3(x3,x2,x1,R*x2+C*x0);
         SHFT2(f2,f1,func(x1));
      }
   }
   if (f1 < f2)
   {
      *xmin = x1;
      return(f1);
   }
   else
   {
      *xmin = x2;
      return(f2);
   }
}

double linint(double xa,double ya,double xb,double yb,double x)
{
   if (xa == xb)
      return((ya+yb)/2);
   else
      return(ya + (yb-ya)/(xb-xa)*(x-xa));
}

double interpolation(double *x,double *y,int n,double u)
{
   int i,j;
   double v,r;

   for (i=0,v=0.0;i<n;i++)
   {
      for (j=0,r=1.0;j<n;j++) if (j != i) r *= (u-x[j])/(x[i]-x[j]);
      v += r*y[i];
   }

   return(v);
}

double pythag(double a,double b)
{
   double absa,absb;

   absa = fabs(a);
   absb = fabs(b);

   if (a == 0) return(absb);
   if (b == 0) return(absa);

   if (absa > absb)
      return(absa*sqrt(1.0+DSQR(absb/absa)));
   else
      return(absb*sqrt(1.0+DSQR(absa/absb)));
}

void svdcmp(double **a,int m,int n,double *w,double **v,int pmod)
{
   int flag,i,its,j,jj,k,l,nm;
   double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
   int pnew,pold,start;

   if (pmod) setbuf(stderr,NULL);

   rv1 = dvector(n);
   g = scale = anorm = 0.0;
   l = nm = 0;

   if (pmod) fprintf(stderr,"Singular value decomposition:   0%%");
   for (i=1,pold=0;i<=n;i++)
   {
      if (pmod)
      {
         pnew = nint((i-1)*34.0/(double)(n-1));
         if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
      }
      l = i+1;
      rv1[i] = scale*g;
      g = s = scale = 0.0;
      if (i <= m)
      {
         for (k=i;k<=m;k++) scale += fabs(a[k][i]);
         if (scale)
         {
            for (k=i;k<=m;k++)
            {
               a[k][i] /= scale;
               s += a[k][i]*a[k][i];
            }
            f = a[i][i];
            g = -SIGN(sqrt(s),f);
            h = f*g - s;
            a[i][i] = f-g;
            for (j=l;j<=n;j++)
            {
               for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
               f = s/h;
               for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
            }
            for (k=i;k<=m;k++) a[k][i] *= scale;
         }
      }
      w[i] = scale*g;
      g = s = scale = 0.0;
      if (i <= m && i != n)
      {
         for (k=l;k<=n;k++) scale += fabs(a[i][k]);
         if (scale)
         {
            for (k=l;k<=n;k++)
            {
               a[i][k] /= scale;
               s += a[i][k]*a[i][k];
            }
            f = a[i][l];
            g = -SIGN(sqrt(s),f);
            h = f*g-s;
            a[i][l] = f-g;
            for (k=l;k<=n;k++) rv1[k] = a[i][k]/h;
            for (j=l;j<=m;j++)
            {
               for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
               for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
            }
            for (k=l;k<=n;k++) a[i][k] *= scale;
         }
      }
      anorm = DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
   }

   for (i=n;i>=1;i--)
   {
      if (i < n)
      {
         if (g)
         {
            for (j=l;j<=n;j++) v[j][i] = (a[i][j]/a[i][l])/g;
            for (j=l;j<=n;j++)
            {
               for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
               for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
            }
         }
         for (j=l;j<=n;j++) v[i][j] = v[j][i] = 0.0;
      }
      v[i][i] = 1.0;
      g = rv1[i];
      l = i;
   }

   for (i=(start=IMIN(m,n)),pold=34;i>=1;i--)
   {
      if (pmod)
      {
         pnew = 34 + nint((start-i)*33.0/(double)(start-1));
         if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
      }
      l = i+1;
      g = w[i];
      for (j=l;j<=n;j++) a[i][j] = 0.0;
      if (g)
      {
         g = 1.0/g;
         for (j=l;j<=n;j++)
         {
            for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
            f = (s/a[i][i])*g;
            for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
         }
         for (j=i;j<=m;j++) a[j][i] *= g;
      }
      else
      {
         for (j=i;j<=m;j++) a[j][i] = 0.0;
      }
      ++a[i][i];
   }

   for (k=n,pold=60;k>=1;k--)
   {
      if (pmod)
      {
         pnew = 67 + nint((n-k)*33.0/(double)(n-1));
         if (pnew > pold) fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
      }
      for (its=1;its<=30;its++)
      {
         flag = 1;
         for (l=k;l>=1;l--)
         {
            nm = l-1;
            if (fabs(rv1[l])+anorm == anorm)
            {
               flag = 0;
               break;
            }
            if (fabs(w[nm])+anorm == anorm) break;
         }
         if (flag)
         {
            c = 0.0;
            s = 1.0;
            for (i=l;i<=k;i++)
            {
               f = s*rv1[i];
               rv1[i] = c*rv1[i];
               if (fabs(f)+anorm == anorm) break;
               g = w[i];
               h = pythag(f,g);
               w[i] = h;
               h = 1.0/h;
               c= g*h;
               s = -f*h;
               for (j=1;j<=m;j++)
               {
                  y = a[j][nm];
                  z = a[j][i];
                  a[j][nm] = y*c + z*s;
                  a[j][i]= z*c - y*s;
               }
            }
         }
         z = w[k];
         if (l == k)
         {
            if (z < 0.0)
            {
               w[k] = -z;
               for (j=1;j<=n;j++) v[j][k] = -v[j][k];
            }
            break;
         }
         if (its == 30) get_error("svdcmp: No convergence achieved!");
         x = w[l];
         nm = k-1;
         y = w[nm];
         g = rv1[nm];
         h = rv1[k];
         f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
         g = pythag(f,1.0);
         f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
         c = s = 1.0;
         for (j=l;j<=nm;j++)
         {
            i = j+1;
            g = rv1[i];
            y = w[i];
            h = s*g;
            g = c*g;
            z = pythag(f,h);
            rv1[j] = z;
            c = f/z;
            s = h/z;
            f = x*c + g*s;
            g = g*c - x*s;
            h = y*s;
            y *= c;
            for (jj=1;jj<=n;jj++)
            {
               x = v[jj][j];
               z = v[jj][i];
               v[jj][j] = x*c + z*s;
               v[jj][i] = z*c - x*s;
            }
            z = pythag(f,h);
            w[j] = z;
            if (z)
            {
               z = 1.0/z;
               c = f*z;
               s = h*z;
            }
            f = c*g + s*y;
            x = c*y - s*g;
            for (jj=1;jj<=m;jj++)
            {
               y = a[jj][j];
               z = a[jj][i];
               a[jj][j] = y*c + z*s;
               a[jj][i] = z*c - y*s;
            }
         }
         rv1[l] = 0.0;
         rv1[k] = f;
         w[k] = x;
      }
   }
   if (pmod) fprintf(stderr,"\n");

   free_dvector(rv1);
}

void svedit(double *w,int n)
{
   double wmax,thresh;
   int i;

   for (i=1,wmax=0.0;i<=n;i++)
      if (w[i] > wmax) wmax = w[i];
   thresh = wmax * 1.0e-10;
   for (i=1;i<=n;i++)
      if (w[i] < thresh) w[i] = 0.0;
}

void svbksb(double **u,double *w,double **v,int m,int n,double *b,double *x)
{
   int jj,j,i;
   double s,*tmp;

   tmp = dvector(n);
   for (j=1;j<=n;j++)
   {
      s = 0.0;
      if (w[j])
      {
         for (i=1;i<=m;i++) s += u[i][j]*b[i];
         s /= w[j];
      }
      tmp[j] = s;
   }
   for (j=1;j<=n;j++)
   {
      s = 0.0;
      for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
      x[j] = s;
   }
   free_dvector(tmp);
}

double intpow(double x,int n)
{
   double y,x2,x4;

   switch(n)
   {
      case 0:  y = 1.0;
               break;
      case 1:  y = x;
               break;
      case 2:  y = x*x;
               break;
      case 3:  y = x*x*x;
               break;
      case 4:  x2 = x*x;
               y = x2*x2;
               break;
      case 5:  x2 = x*x;
               y = x2*x2*x;
               break;
      case 6:  x2 = x*x;
               y = x2*x2*x2;
               break;
      case 7:  x2 = x*x;
               y = x2*x2*x2*x;
               break;
      case 8:  x2 = x*x;
               x4 = x2*x2;
               y = x4*x4;
               break;
      case 9:  x2 = x*x;
               x4 = x2*x2;
               y = x4*x4*x;
               break;
      default: get_error("intpow: Power too high!");
               y = 0;
   }

   return(y);
}

double compute_polynomial(double *a,int n,double x)
{
   double p;
   int i;

   for (i=n-1,p=a[n];i>=1;i--) p = p*x + a[i];

   return(p);
}

double compute_sub_polynomial(double *a,int n,int g,int s,int m,double x)
{
   double *b;
   double p;
   int i;

   b = dvector(m);
   for (i=1;i<=m;i++) b[i] = a[g+(i-1)*s];
   p = compute_polynomial(b,m,x);
   free_dvector(b);

   return(p);
}

void accumulate_coefficients(double *a,int ma,double *b,int mb,
                             double u,int axis)
{
   int i,m;

   m = ma/mb;
   if (mb*m != ma) get_error("accumulate_coefficients: Bad degree!");
   switch(axis)
   {
      case 1: for (i=1;i<=mb;i++)
                 b[i]=compute_sub_polynomial(a,ma,(i-1)*m+1,1,m,u);
              break;
      case 2: for (i=1;i<=mb;i++)
                 b[i]=compute_sub_polynomial(a,ma,i,mb,m,u);
              break;
      default: get_error("accumulate_coefficients: Bad axis!");
   }
}

int number_of_coefficients(int *deg,int naxis)
{
   int axis,m;

   for (axis=m=1;axis<=naxis;axis++) m *= (deg[axis]+1);

   return(m);
}

double polynomial_func(double **x,int dat,int naxis,int *deg,int seq)
{
   double fval;
   div_t q;
   int axis,count;

   fval = 1.0;
   count = seq-1;

   for (axis=1;axis<=naxis;axis++)
   {
      q = div(count,deg[axis]+1);
      fval *= intpow(x[dat][axis],q.rem);
      count = q.quot;
   }

   return(fval);
}

int polynomial_matrix(double **x,int ndat,int naxis,double *y,int *sel,
                       int *deg,double **a,int ncols,double *b)
{
   int dat,col,nsel;

   clear_dmatrix(a,ndat,ncols);
   clear_dvector(b,ndat);

   for (dat=1,nsel=0;dat<=ndat;dat++)
   {
      if (!sel[dat]) continue;
      for (col=1;col<=ncols;col++)
      {
         a[dat][col] = polynomial_func(x,dat,naxis,deg,col);
      }
      b[dat] = y[dat];
      nsel++;
   }

   return(nsel);
}

double polynomial_val(double **x,int dat,int naxis,int *deg,double *a,int ma)
{
   double fval;
   int seq;

   for (seq=1,fval=0.0;seq<=ma;seq++)
      fval += a[seq]*polynomial_func(x,dat,naxis,deg,seq);

   return(fval);
}

void polynomial_fit(double **x,int ndat,int naxis,int *deg,
                    double *a,int ma,double *fit)
{
   int dat;

   for (dat=1;dat<=ndat;dat++)
   {
      fit[dat] = polynomial_val(x,dat,naxis,deg,a,ma);
   }
}

double get_chisq(double *del,int *sel,int ndat)
{
   int dat;
   double chisq;

   for (dat=1,chisq=0.0;dat<=ndat;dat++)
      if (sel[dat]) chisq += del[dat]*del[dat];

   return(chisq);
}

void regression_polynomial(double **x,int ndat,int naxis,double *y,int *sel,
              int *deg,double *a,int ma,double *rms,int *nsel,
              double *fit,double *del,int pmod)
{
   double **xi,*eta;
   int *sig;
   double **u,*b,*w,**v;
   int ntot,i,j,xsize;

   if (pmod) setbuf(stderr,NULL);
   if (pmod) fprintf(stderr,"Polynomial regression -- ");

   for (i=1,ntot=0;i<=ndat;i++) ntot += sel[i];

   xi = dmatrix(ntot,naxis);
   eta = dvector(ntot);
   sig = ivector(ntot);

   xsize = naxis*sizeof(double);

   for (i=j=1;i<=ndat;i++)
   {
      if (!sel[i]) continue;
      memmove(&xi[j][1],&x[i][1],xsize);
      eta[j] = y[i];
      sig[j++] = 1;
   }

   u = dmatrix(ntot,ma);
   b = dvector(ntot);
   v = dmatrix(ma,ma);
   w = dvector(ma);

   *nsel = polynomial_matrix(xi,ntot,naxis,eta,sig,deg,u,ma,b);
   if (*nsel <= ma) get_error("Too few data points!");

   svdcmp(u,ntot,ma,w,v,pmod);
   svedit(w,ma);
   svbksb(u,w,v,ntot,ma,b,a);

   polynomial_fit(x,ndat,naxis,deg,a,ma,fit);
   dvector_difference(y,fit,del,ndat);
   *rms = sqrt(get_chisq(del,sel,ndat)/(*nsel-ma));

   if (pmod) inform("   %d points selected, r.m.s. error: %10.3e",*nsel,*rms);

   free_dmatrix(u);
   free_dvector(b);
   free_dmatrix(v);
   free_dvector(w);
   free_dmatrix(xi);
   free_dvector(eta);
   free_ivector(sig);
}

void best_regression_polynomial(double **x,int ndat,int naxis,double *y,
              int *sel,int *deg,double *a,int ma,double *rms, int *nsel,
              double *fit,double *del,double kappa,int pmod)
{
   double dmax,level;
   int k,kmax;

   do
   {
      regression_polynomial(x,ndat,naxis,y,sel,deg,a,ma,rms,nsel,fit,del,pmod);
      if (kappa == 0.0) break;
      level = kappa*(*rms);
      for (k=1,kmax=0,dmax=-1;k<=ndat;k++)
      {
         if (!sel[k]) continue;
         if (fabs(del[k]) > dmax) dmax=fabs(del[k]), kmax = k;
      }
      if (dmax > level)
      {
         sel[kmax] = 0;
         inform("   A data point rejected with a residual of %10.3e "
                "(%4.1f times sigma)",del[kmax],fabs(del[kmax])/(*rms));
      }
   }
   while(dmax > level);
}

int polynomial_matrix_1d(double *x,double *y,int ndat,int *sel,
                       int deg,double **a,double *b)
{
   int dat,col,nsel;
   int ncols;

   ncols = deg + 1;

   clear_dmatrix(a,ndat,ncols);
   clear_dvector(b,ndat);

   for (dat=1,nsel=0;dat<=ndat;dat++)
   {
      if (!sel[dat]) continue;
      for (col=2,a[dat][1]=1.0;col<=ncols;col++)
      {
         a[dat][col] = a[dat][col-1]*x[dat];
      }
      b[dat] = y[dat];
      nsel++;
   }

   return(nsel);
}

void polynomial_fit_1d(double *x,int ndat,int deg,
                    double *a,double *fit)
{
   int dat;

   for (dat=1;dat<=ndat;dat++)
   {
      fit[dat] = compute_polynomial(a,deg+1,x[dat]);
   }
}

void regression_polynomial_1d(double *x,double *y,int ndat,int *sel,
              int deg,double *a,double *rms,int *nsel,
              double *fit,double *del,int pmod)
{
   double **u,*b,*w,**v;
   int ma;

   if (pmod) inform("1-D polynomial regression:");

   ma = deg+1;

   u = dmatrix(ndat,ma);
   b = dvector(ndat);
   v = dmatrix(ma,ma);
   w = dvector(ma);

   *nsel = polynomial_matrix_1d(x,y,ndat,sel,deg,u,b);
   if (*nsel < ma) get_error("Too few data points!");

   svdcmp(u,ndat,ma,w,v,pmod);
   svedit(w,ma);
   svbksb(u,w,v,ndat,ma,b,a);

   polynomial_fit_1d(x,ndat,deg,a,fit);
   dvector_difference(y,fit,del,ndat);
   if (*nsel > ma)
      *rms = sqrt(get_chisq(del,sel,ndat)/(*nsel-ma));
   else
      *rms = 0.0;

   if (pmod) inform("%d points selected, r.m.s. error: %10.3e",*nsel,*rms);

   free_dmatrix(u);
   free_dvector(b);
   free_dmatrix(v);
   free_dvector(w);
}

void best_regression_polynomial_1d(double *x,double *y,int ndat,
              int *sel,int deg,double *a,double *rms, int *nsel,
              double *fit,double *del,double kappa,int pmod)
{
   double dmax,level;
   int k,kmax;

   do
   {
      regression_polynomial_1d(x,y,ndat,sel,deg,a,rms,nsel,fit,del,pmod);
      if (kappa == 0.0) break;
      level = kappa*(*rms);
      for (k=1,kmax=0,dmax=-1;k<=ndat;k++)
      {
         if (!sel[k]) continue;
         if (fabs(del[k]) > dmax) dmax=fabs(del[k]), kmax = k;
      }
      if (dmax > level) sel[kmax] = 0;
   }
   while(dmax > level);
}

int gaussj(double **a, int n, double **b, int m)
{
   int *indxc,*indxr,*ipiv;
   int i,icol,irow,j,k,l,ll;
   double big,dum,pivinv;

   indxc=ivector(n);
   indxr=ivector(n);
   ipiv=ivector(n);

   for (j=1;j<=n;j++) ipiv[j]=0;
   icol=irow=1;
   for (i=1;i<=n;i++)
   {
      big=0.0;
      for (j=1;j<=n;j++)
      {
         if (ipiv[j] != 1)
         {
            for (k=1;k<=n;k++)
            {
               if (ipiv[k] == 0)
               {
                  if (fabs(a[j][k]) >= big) big=fabs(a[j][k]),irow=j,icol=k;
               }
               else
               {
                  if (ipiv[k] > 1)
                  {
                     CLEAN_GAUSSJ;
                     return(GAUSSJ_SINGMAT_1);
                  }
               }
            }
         }
      }
      ipiv[icol]++;
      if (irow != icol)
      {
         for (l=1;l<=n;l++) DSWAP(a[irow][l],a[icol][l]);
         for (l=1;l<=m;l++) DSWAP(b[irow][l],b[icol][l]);
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if (a[icol][icol] == 0.0)
      {
         CLEAN_GAUSSJ;
         return(GAUSSJ_SINGMAT_2);
      }
      pivinv=1.0/a[icol][icol];
      a[icol][icol]=1.0;
      for (l=1;l<=n;l++) a[icol][l] *= pivinv;
      for (l=1;l<=m;l++) b[icol][l] *= pivinv;
      for (ll=1;ll<=n;ll++)
      {
         if (ll != icol)
         {
            dum=a[ll][icol];
            a[ll][icol]=0.0;
            for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
            for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
         }
      }
   }
   for (l=n;l>=1;l--)
      if (indxr[l] != indxc[l])
         for (k=1;k<=n;k++) DSWAP(a[k][indxr[l]],a[k][indxc[l]]);

   CLEAN_GAUSSJ;

   return(0);
}

void covsrt(double **covar, int ma, int *ia, int mfit)
{
   int i,j,k;

   for (i=mfit+1;i<=ma;i++) for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
   k=mfit;
   for (j=ma;j>=1;j--)
   {
      if (ia[j])
      {
         for (i=1;i<=ma;i++) DSWAP(covar[i][k],covar[i][j]);
         for (i=1;i<=ma;i++) DSWAP(covar[k][i],covar[j][i]);
         k--;
      }
   }
}

void mrqcof(double *x, double *y, double *sigma, int ndata,
   double *a, int *ia, int ma,
   double **alpha, double *beta, double *chisq,
   void (*funcs)(double, double *, double *, double *, int))
{
   int i,j,k,l,m,mfit=0;
   double ymod,wt,sigma2i,dy,*dyda;

   dyda=dvector(ma);
   for (j=1;j<=ma;j++) if (ia[j]) mfit++;
   for (j=1;j<=mfit;j++)
   {
      for (k=1;k<=j;k++) alpha[j][k]=0.0;
      beta[j]=0.0;
   }
   *chisq=0.0;
   for (i=1;i<=ndata;i++)
   {
      (*funcs)(x[i],a,&ymod,dyda,ma);
      sigma2i=1.0/(sigma[i]*sigma[i]);
      dy=y[i]-ymod;
      for (j=0,l=1;l<=ma;l++)
      {
         if (ia[l])
         {
            wt=dyda[l]*sigma2i;
            for (j++,k=0,m=1;m<=l;m++) if (ia[m]) alpha[j][++k] += wt*dyda[m];
            beta[j] += dy*wt;
         }
      }
      *chisq += dy*dy*sigma2i;
   }
   for (j=2;j<=mfit;j++) for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
   free_dvector(dyda);
}

int mrqmin(double *x, double *y, double *sigma, int ndata,
   double *a, int *ia, int ma,
   double **covar, double **alpha, double *chisq,
   void (*funcs)(double, double *, double *, double *, int), double *alambda)
{
   int j,k,l,flag;
   static int mfit;
   static double ochisq,*atry,*beta,*da,**oneda;

   if (alambda == NULL)
   {
      CLEAN_MRQMIN;
      return(0);
   }

   if (*alambda < 0.0)
   {
      atry=dvector(ma);
      beta=dvector(ma);
      da=dvector(ma);
      for (mfit=0,j=1;j<=ma;j++) if (ia[j]) mfit++;
      oneda=dmatrix(mfit,1);
      *alambda=0.001;
      mrqcof(x,y,sigma,ndata,a,ia,ma,alpha,beta,chisq,funcs);
      ochisq=(*chisq);
      for (j=1;j<=ma;j++) atry[j]=a[j];
   }

   for (j=1;j<=mfit;j++)
   {
      for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
      covar[j][j]=alpha[j][j]*(1.0+(*alambda));
      oneda[j][1]=beta[j];
   }

   flag = gaussj(covar,mfit,oneda,1);
   if (flag != 0)
   {
      CLEAN_MRQMIN;
      return(flag);
   }

   for (j=1;j<=mfit;j++) da[j]=oneda[j][1];

   if (*alambda == 0.0)
   {
      covsrt(covar,ma,ia,mfit);
      CLEAN_MRQMIN;
      return(0);
   }

   for (j=0,l=1;l<=ma;l++) if (ia[l]) atry[l]=a[l]+da[++j];
   mrqcof(x,y,sigma,ndata,atry,ia,ma,covar,da,chisq,funcs);

   if (*chisq < ochisq)
   {
      *alambda *= 0.1;
      ochisq = (*chisq);
      for (j=1;j<=mfit;j++)
      {
         for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
         beta[j]=da[j];
      }
      for (l=1;l<=ma;l++) a[l]=atry[l];
   }
   else
   {
      *alambda *= 10.0;
      *chisq = ochisq;
   }

   return(0);
}

int nonlinfit(double *x,double *y,int ndata,double *a,int *ia,double *s,
   double *eps,int ma,int nmax,double *chisq,double *rms,int *niter,
   void (*funcs)(double,double *,double *,double *,int))
{
   double *sigma,*b,**covar,**alpha;
   double alam;
   int done,i,flag;

   sigma=dvector(ndata);
   b=dvector(ma);
   alpha=dmatrix(ma,ma);
   covar=dmatrix(ma,ma);

   for (i=1;i<=ndata;i++) sigma[i]=1.0;

   for (*niter=1,alam=-1.0;*niter<=nmax;(*niter)++)
   {
      copy_dvector(a,ma,b);
      flag = mrqmin(x,y,sigma,ndata,a,ia,ma,covar,alpha,chisq,funcs,&alam);
      if (flag != 0)
      {
         CLEAN_NONLINFIT;
         return(flag);
      }
      for (i=1,done=1;i<=ma;i++) if (fabs(a[i]-b[i]) > eps[i]) done=0;
      if (done) break;
   }

   if (!done)
   {
      mrqmin(x,y,sigma,ndata,a,ia,ma,covar,alpha,chisq,funcs,NULL);
      CLEAN_NONLINFIT;
      return(NONLINFIT_NO_CONVERGENCE);
   }

   alam=0.0;
   flag = mrqmin(x,y,sigma,ndata,a,ia,ma,covar,alpha,chisq,funcs,&alam);
   if (flag != 0)
   {
      CLEAN_NONLINFIT;
      return(flag+MRQMIN_FLAG_OFFSET);
   }

   for (i=1;i<=ma;i++) s[i]=sqrt(covar[i][i])*sqrt((*chisq)/(ndata-ma));

   *rms = sqrt((*chisq)/(ndata-ma));

   CLEAN_NONLINFIT;

   return(0);
}

void fgauss(double x,double *a,double *f,double *dfda,int na)
{
   int i;
   double fac,ex,arg;

   *f = a[1];
   dfda[1] = 1.0;
   for (i=2;i<na;i+=3)
   {
      arg = (x-a[i+1]) / a[i+2];
      ex = exp(-arg*arg);
      fac = 2.0*a[i]*ex*arg/a[i+2];
      *f += a[i]*ex;
      dfda[i] = ex;
      dfda[i+1] = fac;
      dfda[i+2] = fac*arg;
   }
}

void fgauss2d(double x,double *a,double *f,double *dfda,int na)
{
   int i,k;
   double fac,fac1,fac2,ex,arg1,arg2;

   k = (int)x;
   *f = a[1];
   dfda[1] = 1.0;
   for (i=2;i<na;i+=6)
   {
      arg1 = (xypoint[k][1]-a[i+1]) / a[i+2];
      arg2 = (xypoint[k][2]-a[i+3]) / a[i+4];
      ex = exp(-arg1*arg1-arg2*arg2+2*a[i+5]*arg1*arg2);
      fac = 2.0*a[i]*ex;
      fac1 = fac*(arg1-a[i+5]*arg2)/a[i+2];
      fac2 = fac*(arg2-a[i+5]*arg1)/a[i+4];
      *f += a[i]*ex;
      dfda[i] = ex;
      dfda[i+1] = fac1;
      dfda[i+2] = fac1*arg1;
      dfda[i+3] = fac2;
      dfda[i+4] = fac2*arg2;
      dfda[i+5] = fac*arg1*arg2;
   }
}

void clear_gaussian(gaussian *gauss)
{
   memset(gauss,0,sizeof(gaussian));
}

void fit_gaussian_1d(double *x,double *f,int n,gaussian *gauss)
{
   int i;
   double xref,xscale,fref,fscale,xfwhm_min,xfwhm_max;
   double a[5],s[5],eps[5];
   int ia[5];
   double *xi,*eta;

   xref = gauss->xcen;
   xscale = gauss->xfwhm / FWHM_FACTOR;
   fref = gauss->base;
   fscale = gauss->icen;

   xfwhm_min = gauss->xfwhm/FWHM_TOLLERANCE;
   xfwhm_max = gauss->xfwhm*FWHM_TOLLERANCE;

   clear_gaussian(gauss);

   if (n < 5)
   {
      gauss->status = G1DFIT_SMALL_NDAT;
      return;
   }
   if (xscale < GFIT_SCALE_MIN)
   {
      gauss->status = G1DFIT_SMALL_XSCALE;
      return;
   }
   if (fscale < GFIT_SCALE_MIN)
   {
      gauss->status = G1DFIT_SMALL_FSCALE;
      return;
   }

   xi = dvector(n);
   eta = dvector(n);

   for (i=1;i<=n;i++)
   {
      xi[i] = (x[i] - xref) / xscale;
      eta[i] = (f[i] - fref) / fscale;
   }

   a[1] = 0.0;
   a[2] = 1.0;
   a[3] = 0.0;
   a[4] = 1.0;

   for (i=1;i<=4;i++) ia[i]=1,eps[i]=G1DFIT_EPS;

   gauss->status = nonlinfit(xi,eta,n,a,ia,s,eps,4,G1DFIT_MAXPASS,
                         &gauss->chisq,&gauss->rms,&gauss->niter,fgauss);

   if (gauss->status != 0) clear_dvector(a,4);

   gauss->base = a[1] * fscale + fref;
   gauss->icen = a[2] * fscale;
   gauss->xcen = a[3] * xscale + xref;
   gauss->xfwhm = a[4] * xscale * FWHM_FACTOR;
   gauss->chisq *= fscale*fscale;
   gauss->rms *= fscale;

   free_dvector(xi);
   free_dvector(eta);

   if (gauss->status == 0 && gauss->rms > GFIT_RMS_MAX*fscale)
      gauss->status = G1DFIT_LARGE_RMS;

   if (gauss->status == 0 && (gauss->xfwhm<xfwhm_min ||
                              gauss->xfwhm>xfwhm_max))
      gauss->status = G1DFIT_BAD_XFWHM;

   if (gauss->status == 0 && fabs(gauss->xcen-xref)/gauss->xfwhm > 0.5)
      gauss->status = G1DFIT_BAD_XCEN;
}

void center_gaussian_1d(float *pix,int npix,double start,double step,
                        int a,int b,double w,gaussian *gauss)
{
   double *x,*f;
   double fmin,fmax;
   int i,k,imax,n;

   clear_gaussian(gauss);

   if (b <= a)
   {
      gauss->status = G1DCENT_BAD_XREGION;
      return;
   }

   if (a < 1 || b > npix)
   {
      gauss->status = G1DCENT_XREGION_OUT;
      return;
   }

   n = b - a + 1;

   if (n < 5)
   {
      gauss->status = G1DCENT_SMALL_NDAT;
      return;
   }

   x = dvector(n);
   f = dvector(n);

   for (k=a,i=1;k<=b;k++,i++)
   {
      x[i] = start + (k-1)*step;
      f[i] = (double)pix[k-1];
   }

   for (i=imax=1,fmin=fmax=f[1];i<=n;i++)
   {
      if (f[i] < fmin) fmin = f[i];
      if (f[i] > fmax) fmax = f[i], imax = i;
   }

   if (fmax == fmin)
   {
      gauss->status = G1DCENT_ZERO_FRANGE;
      CLEAN_CENTGAUSS1D;
      return;
   }

   gauss->base = fmin;
   gauss->icen = fmax - fmin;
   gauss->xcen = x[imax];
   gauss->xfwhm = w;

   if (x[n]-x[1] < gauss->xfwhm)
   {
      gauss->status = G1DCENT_BAD_XFWHM;
      CLEAN_CENTGAUSS1D;
      return;
   }

   fit_gaussian_1d(x,f,n,gauss);

   CLEAN_CENTGAUSS1D;
}

void locate_gaussian_1d(float *pix,int npix,double start,double step,
                        double a,double c,double b,double w,gaussian *gauss)
{
   int ka,kb,k,kmax;
   float maxpix;

   clear_gaussian(gauss);

   if (a >= c || b <= c)
   {
      gauss->status = G1DLOC_BAD_XREGION;
      return;
   }

   ka = nint((a-start)/step) + 1;
   kb = nint((b-start)/step) + 1;

   if (ka < 1 || kb > npix)
   {
      gauss->status = G1DLOC_XREGION_OUT;
      return;
   }

   for (k=kmax=ka,maxpix=pix[ka-1];k<=kb;k++)
      if (pix[k-1] > maxpix) maxpix=pix[k-1], kmax=k;

   if (kmax==ka || kmax==kb)
   {
      gauss->status = G1DLOC_BAD_XMAX;
      return;
   }

   ka = kmax - nint((c-a)/step);
   kb = kmax + nint((b-c)/step);

   if (ka < 1 || kb > npix)
   {
      gauss->status = G1DLOC_XREGION_OUT;
      return;
   }

   for (k=ka;k<=kb;k++)
   {
      if (k == kmax) continue;
      if (pix[k-1] >= maxpix)
      {
         gauss->status = G1DLOC_BAD_MAX;
         return;
      }
   }

   center_gaussian_1d(pix,npix,start,step,ka,kb,w,gauss);
}

void fit_gaussian_2d(double **xy,double *f,int n,gaussian *gauss)
{
   int i;
   double xref,xscale,yref,yscale,fref,fscale;
   double xfwhm_min,xfwhm_max,yfwhm_min,yfwhm_max;
   double a[8],s[8],eps[8];
   int ia[8];
   double **xi,*eta,*u;

   xref = gauss->xcen;
   xscale = gauss->xfwhm / FWHM_FACTOR;
   yref = gauss->ycen;
   yscale = gauss->yfwhm / FWHM_FACTOR;
   fref = gauss->base;
   fscale = gauss->icen;

   xfwhm_min = gauss->xfwhm/FWHM_TOLLERANCE;
   xfwhm_max = gauss->xfwhm*FWHM_TOLLERANCE;
   yfwhm_min = gauss->yfwhm/FWHM_TOLLERANCE;
   yfwhm_max = gauss->yfwhm*FWHM_TOLLERANCE;

   clear_gaussian(gauss);

   if (n < 25)
   {
      gauss->status = G2DFIT_SMALL_NDAT;
      return;
   }
   if (xscale < GFIT_SCALE_MIN)
   {
      gauss->status = G2DFIT_SMALL_XSCALE;
      return;
   }
   if (yscale < GFIT_SCALE_MIN)
   {
      gauss->status = G2DFIT_SMALL_YSCALE;
      return;
   }
   if (fscale < GFIT_SCALE_MIN)
   {
      gauss->status = G2DFIT_SMALL_FSCALE;
      return;
   }

   xi = dmatrix(n,2);
   eta = dvector(n);

   u = dvector(n);

   xypoint = xi;

   for (i=1;i<=n;i++)
   {
      xi[i][1] = (xy[i][1] - xref) / xscale;
      xi[i][2] = (xy[i][2] - yref) / yscale;
      eta[i] = (f[i] - fref) / fscale;
      u[i] = (double)i;
   }

   a[1] = 0.0;
   a[2] = 1.0;
   a[3] = 0.0;
   a[4] = 1.0;
   a[5] = 0.0;
   a[6] = 1.0;
   a[7] = gauss->xyfactor;

   for (i=1;i<=7;i++) ia[i]=1,eps[i]=G2DFIT_EPS;

   gauss->status = nonlinfit(u,eta,n,a,ia,s,eps,7,G2DFIT_MAXPASS,
                         &gauss->chisq,&gauss->rms,&gauss->niter,fgauss2d);

   if (gauss->status != 0) clear_dvector(a,7);

   gauss->base = a[1] * fscale + fref;
   gauss->icen = a[2] * fscale;
   gauss->xcen = a[3] * xscale + xref;
   gauss->xfwhm = a[4] * xscale * FWHM_FACTOR;
   gauss->ycen = a[5] * yscale + yref;
   gauss->yfwhm = a[6] * yscale * FWHM_FACTOR;
   gauss->xyfactor = a[7];
   gauss->chisq *= fscale*fscale;
   gauss->rms *= fscale;

   free_dmatrix(xi);
   free_dvector(eta);
   free_dvector(u);

   if (gauss->status == 0 && gauss->rms > GFIT_RMS_MAX*fscale)
      gauss->status = G2DFIT_LARGE_RMS;

   if (gauss->status == 0 && (gauss->xfwhm<xfwhm_min ||
                              gauss->xfwhm>xfwhm_max))
      gauss->status = G2DFIT_BAD_XFWHM;
   if (gauss->status == 0 && (gauss->yfwhm<yfwhm_min ||
                              gauss->yfwhm>yfwhm_max))
      gauss->status = G2DFIT_BAD_YFWHM;

   if (gauss->status == 0 && fabs(gauss->xcen-xref)/gauss->xfwhm > 0.5)
      gauss->status = G2DFIT_BAD_XCEN;
   if (gauss->status == 0 && fabs(gauss->ycen-yref)/gauss->yfwhm > 0.5)
      gauss->status = G2DFIT_BAD_YCEN;
}

void center_gaussian_2d(float *pix,int npix1,int npix2,
                        double start1,double start2,double step1,double step2,
                        int a1,int a2,int b1,int b2,
                        double wx,double wy,gaussian *gauss)
{
   double **xy,*f;
   double fmin,fmax;
   int i,k1,k2,imax,n1,n2,n;
   float *ptr,*src;

   clear_gaussian(gauss);

   if (b1 <= a1)
   {
      gauss->status = G2DCENT_BAD_XREGION;
      return;
   }

   if (b2 <= a2)
   {
      gauss->status = G2DCENT_BAD_YREGION;
      return;
   }

   if (a1 < 1 || b1 > npix1)
   {
      gauss->status = G2DCENT_XREGION_OUT;
      return;
   }

   if (a2 < 1 || b2 > npix2)
   {
      gauss->status = G2DCENT_YREGION_OUT;
      return;
   }

   n1 = b1 - a1 + 1;
   n2 = b2 - a2 + 1;

   if (n1 < 5)
   {
      gauss->status = G2DCENT_SMALL_NXDAT;
      return;
   }

   if (n2 < 5)
   {
      gauss->status = G2DCENT_SMALL_NYDAT;
      return;
   }

   n = n1*n2;

   xy = dmatrix(n,2);
   f = dvector(n);

   ptr = pix + (a2-1)*npix1 + a1 - 1;
   for (k2=a2,i=1;k2<=b2;k2++,ptr+=npix1)
   {
      for (k1=a1,src=ptr;k1<=b1;k1++,i++,src++)
      {
         xy[i][1] = start1 + (k1-1)*step1;
         xy[i][2] = start2 + (k2-1)*step2;
         f[i] = (double)(*src);
      }
   }

   for (i=imax=1,fmin=fmax=f[1];i<=n;i++)
   {
      if (f[i] < fmin) fmin = f[i];
      if (f[i] > fmax) fmax = f[i], imax = i;
   }

   if (fmax == fmin)
   {
      gauss->status = G2DCENT_ZERO_FRANGE;
      CLEAN_CENTGAUSS2D;
      return;
   }

   gauss->base = fmin;
   gauss->icen = fmax - fmin;
   gauss->xcen = xy[imax][1];
   gauss->xfwhm = wx;
   gauss->ycen = xy[imax][2];
   gauss->yfwhm = wy;
   gauss->xyfactor = 0.0;

   if (xy[1][n]-xy[1][1] < gauss->xfwhm)
   {
      gauss->status = G2DCENT_BAD_XFWHM;
      CLEAN_CENTGAUSS2D;
      return;
   }

   if (xy[2][n]-xy[2][1] < gauss->yfwhm)
   {
      gauss->status = G2DCENT_BAD_YFWHM;
      CLEAN_CENTGAUSS2D;
      return;
   }

   fit_gaussian_2d(xy,f,n,gauss);

   CLEAN_CENTGAUSS2D;
}

void locate_gaussian_2d(float *pix,int npix1,int npix2,
                     double start1,double start2,double step1,double step2,
                     double a1,double a2,double c1,double c2,
                     double b1,double b2,double wx,double wy,gaussian *gauss)
{
   int ka1,ka2,kb1,kb2,kmax1,kmax2,k1,k2,ka,ks,k;
   float maxpix;

   clear_gaussian(gauss);

   if (a1 >= c1 || b1 <= c1)
   {
      gauss->status = G2DLOC_BAD_XREGION;
      return;
   }

   if (a2 >= c2 || b2 <= c2)
   {
      gauss->status = G2DLOC_BAD_YREGION;
      return;
   }

   ka1 = nint((a1-start1)/step1) + 1;
   ka2 = nint((a2-start2)/step2) + 1;
   kb1 = nint((b1-start1)/step1) + 1;
   kb2 = nint((b2-start2)/step2) + 1;

   if (ka1 < 1 || kb1 > npix1)
   {
      gauss->status = G2DLOC_XREGION_OUT;
      return;
   }

   if (ka2 < 1 || kb2 > npix2)
   {
      gauss->status = G2DLOC_YREGION_OUT;
      return;
   }

   ka = (ka2-1)*npix1 + ka1;
   maxpix = pix[ka-1];
   for (k2=ka2,ks=ka;k2<=kb2;k2++,ks+=npix1)
   {
      for (k1=ka1,k=ks;k1<=kb1;k1++,k++)
      {
         if (pix[k-1] > maxpix) maxpix=pix[k-1], kmax1=k1, kmax2=k2;
      }
   }

   if (kmax1==ka1 || kmax1==kb1)
   {
      gauss->status = G2DLOC_BAD_XMAX;
      return;
   }

   if (kmax2==ka2 || kmax2==kb2)
   {
      gauss->status = G2DLOC_BAD_YMAX;
      return;
   }

   ka1 = kmax1 - nint((c1-a1)/step1);
   ka2 = kmax2 - nint((c2-a2)/step2);
   kb1 = kmax1 + nint((b1-c1)/step1);
   kb2 = kmax2 + nint((b2-c2)/step2);

   if (ka1 < 1 || kb1 > npix1)
   {
      gauss->status = G2DLOC_XREGION_OUT;
      return;
   }

   if (ka2 < 1 || kb2 > npix2)
   {
      gauss->status = G2DLOC_YREGION_OUT;
      return;
   }

   ka = (ka2-1)*npix1 + ka1;
   for (k2=ka2,ks=ka;k2<=kb2;k2++,ks+=npix1)
   {
      for (k1=ka1,k=ks;k1<=kb1;k1++,k++)
      {
         if (k1 == kmax1 && k2 == kmax2) continue;
         if (pix[k-1] >= maxpix)
         {
            gauss->status = G2DLOC_BAD_MAX;
            return;
         }
      }
   }

   center_gaussian_2d(pix,npix1,npix2,start1,start2,step1,step2,
                      ka1,ka2,kb1,kb2,wx,wy,gauss);
}

int peak_ok(double *x,int n)
{
   int i;

   for (i=2;i<=n;i++) if (x[i] <= x[i-1]) break;
   if (i >= n) return(0);
   for (i=i+1;i<=n;i++) if (x[i] >= x[i-1]) return(0);

   return(1);
}

void init_parab(double *y,int ny,int *imax,double *xparab,
                double *yparab,int crad,int pmod)
{
   int i,j,n;
   double ymax;
   double *u,*v,*fit,*del;
   int *sel,nsel;
   double a[4],rms;

   n = 2 * crad + 1;

   for (i=1,*imax=1,ymax=y[1];i<=ny;i++) if (y[i] > ymax) ymax=y[i],*imax=i;

   if (*imax <= crad || *imax > ny-crad)
      get_error("CCF maximum too close to the image edges!");

   *xparab = 0.0;
   *yparab = ymax;

   if (crad < 1) return;

   u = dvector(n);
   v = dvector(n);
   fit = dvector(n);
   del = dvector(n);
   sel = ivector(n);

   for (i=-crad,j=1;i<=crad;i++,j++)
   {
      u[j] = i;
      v[j] = y[*imax+i];
      sel[j] = 1;
   }

   if (!peak_ok(v,n))
      get_error("The core of the CCF peak not well defined!");

   regression_polynomial_1d(u,v,n,sel,2,a,&rms,&nsel,fit,del,pmod);

   *xparab = -a[2]/(2*a[3]);
   *yparab = a[1] - a[2]*a[2]/(4*a[3]);

   free_dvector(u);
   free_dvector(v);
   free_dvector(fit);
   free_dvector(del);
   free_ivector(sel);
}

void define_core_peak(int imax,double xparab,int rad,int *kg,int *ng,
                      int *b1,int *b2)
{
   double h;

   h = floor(2*xparab+0.5)/2;
   *kg = imax + (int)floor(h) - rad;
   *ng = 2*rad + 1;
   if (floor(h) != h) (*ng)++;
   *b1 = (*kg);
   *b2 = (*kg) + (*ng) - 1;
}

void extract_core_peak(double *u,double *v,double *y,int imax,int kg,int ng)
{
   int i,j;

   for (i=1,j=kg;i<=ng;i++,j++)
   {
      u[i] = j - imax;
      v[i] = y[j];
   }
}

void init_gauss(gaussian *gauss,double *v,int ng,double xparab,double yparab)
{
   int i,ka,kb;
   double half;

   gauss->base = DMIN(v[1],v[ng]);
   gauss->icen = yparab - gauss->base;
   gauss->xcen = xparab;

   half = gauss->base + gauss->icen/2;

   for (i=1;i<=ng;i++) if (v[i] > half) break;
   ka = i;
   for (i=ka;i<=ng;i++) if (v[i] < half) break;
   kb = i;

   gauss->xfwhm = (double)(kb-ka+1);
}

double get_shift(gaussian *gauss,int imax,double x1,double step)
{
   gauss->xcen += (double)imax;
   return(x1 + (gauss->xcen-1)*step);
}

double locate_gauss_maximum(double x1,double step,double *y,int n,int crad,
                            int rad,int pmod,int *b1,int *b2,gaussian *gauss)
{
   int imax;
   double *u,*v;
   double xparab,yparab;
   int kg,ng;

   init_parab(y,n,&imax,&xparab,&yparab,crad,pmod);

   define_core_peak(imax,xparab,rad,&kg,&ng,b1,b2);

   u = dvector(ng);
   v = dvector(ng);

   extract_core_peak(u,v,y,imax,kg,ng);

   clear_gaussian(gauss);

   init_gauss(gauss,v,ng,xparab,yparab);

   fit_gaussian_1d(u,v,ng,gauss);

   free_dvector(u);
   free_dvector(v);

   return(get_shift(gauss,imax,x1,step));
}

double locate_parab_maximum(double x1,double step,double *y,int n,int crad,
                            int rad,int pmod,int *b1,int *b2,gaussian *gauss)
{
   int i,imax;
   double *u,*v,*fit,*del;
   int *sel,nsel;
   double a[4],rms;
   double xparab,yparab;
   int kg,ng;

   init_parab(y,n,&imax,&xparab,&yparab,crad,pmod);

   define_core_peak(imax,xparab,rad,&kg,&ng,b1,b2);

   u = dvector(ng);
   v = dvector(ng);
   fit = dvector(ng);
   del = dvector(ng);
   sel = ivector(ng);

   extract_core_peak(u,v,y,imax,kg,ng);
   for (i=1;i<=ng;i++) sel[i] = 1;

   clear_gaussian(gauss);

   regression_polynomial_1d(u,v,ng,sel,2,a,&rms,&nsel,fit,del,pmod);

   gauss->icen = a[1] - a[2]*a[2]/(4*a[3]);
   gauss->xcen = -a[2]/(2*a[3]);
   gauss->rms = rms;
   gauss->status = 0;

   free_dvector(u);
   free_dvector(v);
   free_dvector(fit);
   free_dvector(del);
   free_ivector(sel);

   return(get_shift(gauss,imax,x1,step));
}

double inter_func(double u)
{
   double v;

   splint(x_vec,y_vec,y2_vec,y_count,u,&v);
   return(-v);
}

double locate_spline_maximum(double x1,double step,double *y,int n,int crad,
                             double tol,int pmod,gaussian *gauss)
{
   int i,imax;
   double xparab,yparab;

   y_vec = y;
   y_count = n;

   init_parab(y,n,&imax,&xparab,&yparab,crad,pmod);

   x_vec = dvector(n);
   y2_vec = dvector(n);

   for (i=1;i<=n;i++) x_vec[i] = (double)(i-imax);

   spline(x_vec,y,n,1e30,1e30,y2_vec);

   gauss->icen = -golden(-2,0,2,inter_func,tol,&gauss->xcen);
   gauss->status = 0;

   free_dvector(y2_vec);
   free_dvector(x_vec);

   return(get_shift(gauss,imax,x1,step));
}

double locate_peak_maximum(double x1,double step,double *y,int n,int crad,
                           int pmod,int *b1,int *b2,gaussian *gauss)
{
   int imax,i;
   double xparab,yparab;
   double a,b,c;

   init_parab(y,n,&imax,&xparab,&yparab,crad,pmod);

   for (i=imax-2;i<=imax+2;i++)
      if (y[i]>=y[i-1] && y[i]>y[i+1]) break;
   if (abs(i-imax) > 1)
      get_error("locate_peak_maximum: CCF maximum not found!");

   *b1 = i-1;
   *b2 = i+1;

   a = (y[i-1]+y[i+1]-2*y[i])/2;
   b = (y[i+1]-y[i-1])/2;
   c = y[i];

   gauss->xcen = -b/(2*a);
   gauss->icen = c - b*b/(4*a);
   gauss->status = 0;

   return(get_shift(gauss,i,x1,step));
}

double locate_ccf_maximum
   (double x1,double step,double *y,int n,gaussian *gauss)
{
   int i,imax;
   double ymax;

   for (i=imax=1,ymax=y[1];i<=n;i++) if (y[i] > ymax) ymax=y[i],imax=i;

   gauss->xcen = x1 + (imax-1)*step;
   gauss->icen = ymax;
   gauss->status = 0;

   return(get_shift(gauss,imax,x1,step));
}

void four1(double *data,int32_t nn,int isign)
{
   int32_t n,mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

   n = nn<<1;
   j=1;
   for (i=1;i<n;i+=2)
   {
      if (j > i)
      {
         DSWAP(data[j],data[i]);
         DSWAP(data[j+1],data[i+1]);
      }
      m = n>>1;
      while (m>=2 && j>m)
      {
         j -= m;
         m >>= 1;
      }
      j += m;
   }
   mmax=2;
   while (n > mmax)
   {
      istep = mmax<<1;
      theta = isign*(TWOPI/mmax);
      wpr=-2.0*DSQR(sin(0.5*theta));
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2)
      {
         for (i=m;i<=n;i+=istep)
         {
            j=i+mmax;
            tempr=wr*data[j]-wi*data[j+1];
            tempi=wr*data[j+1]+wi*data[j];
            data[j]=data[i]-tempr;
            data[j+1]=data[i+1]-tempi;
            data[i] += tempr;
            data[i+1] += tempi;
         }
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
      }
      mmax=istep;
   }
}

void twofft(double *data1,double *data2,double *fft1,double *fft2,int32_t n)
{
   int32_t nn3,nn2,jj,j;
   double rep,rem,aip,aim;

   nn3 = 1+(nn2=2+n+n);
   for (j=1,jj=2;j<=n;j++,jj+=2) fft1[jj-1]=data1[j],fft1[jj]=data2[j];
   four1(fft1,n,1);
   fft2[1]=fft1[2];
   fft1[2]=fft2[2]=0.0;
   for (j=3;j<=n+1;j+=2)
   {
      rep = 0.5*(fft1[j]+fft1[nn2-j]);
      rem = 0.5*(fft1[j]-fft1[nn2-j]);
      aip = 0.5*(fft1[j+1]+fft1[nn3-j]);
      aim = 0.5*(fft1[j+1]-fft1[nn3-j]);
      fft1[j] = rep;
      fft1[j+1] = aim;
      fft1[nn2-j] = rep;
      fft1[nn3-j] = -aim;
      fft2[j] = aip;
      fft2[j+1] = -rem;
      fft2[nn2-j] = aip;
      fft2[nn3-j] = rem;
   }
}

void realft(double *data,int32_t n,int isign)
{
   int32_t i,i1,i2,i3,i4,np3;
   double c1=0.5,c2,h1r,h1i,h2r,h2i;
   double wr,wi,wpr,wpi,wtemp,theta;

   theta = PI/(double)(n>>1);
   if (isign == 1)
      c2=-0.5, four1(data,n>>1,1);
   else
      c2=0.5, theta=-theta;
   wtemp = sin(0.5*theta);
   wpr = -2.0*wtemp*wtemp;
   wpi = sin(theta);
   wr = 1.0+wpr;
   wi = wpi;
   np3 = n+3;
   for (i=2;i<=(n>>2);i++)
   {
      i4 = 1+(i3=np3-(i2=1+(i1=i+i-1)));
      h1r = c1*(data[i1]+data[i3]);
      h1i = c1*(data[i2]-data[i4]);
      h2r = -c2*(data[i2]+data[i4]);
      h2i = c2*(data[i1]-data[i3]);
      data[i1] = h1r+wr*h2r-wi*h2i;
      data[i2] = h1i+wr*h2i+wi*h2r;
      data[i3] = h1r-wr*h2r+wi*h2i;
      data[i4] = -h1i+wr*h2i+wi*h2r;
      wr = (wtemp=wr)*wpr-wi*wpi+wr;
      wi = wi*wpr+wtemp*wpi+wi;
   }
   if (isign == 1)
   {
      data[1] = (h1r=data[1])+data[2];
      data[2] = h1r-data[2];
   }
   else
   {
      data[1] = c1*((h1r=data[1])+data[2]);
      data[2] = c1*(h1r-data[2]);
      four1(data,n>>1,-1);
   }
}

void convlv(double *data,int32_t n,double *respns,int32_t m,int isign,double *ans)
{
   int32_t i,no2;
   double dum,mag2,*fft;

   fft = dvector(n<<1);
   for (i=1;i<=(m-1)/2;i++) respns[n+1-i]=respns[m+1-i];
   for (i=(m+3)/2;i<=n-(m-1)/2;i++) respns[i]=0.0;

   twofft(data,respns,fft,ans,n);

   no2 = n>>1;
   for (i=2;i<=n+2;i+=2)
   {
      if (isign == 1)
      {
         ans[i-1] = (fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
         ans[i] = (fft[i]*dum+fft[i-1]*ans[i])/no2;
      }
      else
      {
         if ((mag2=ans[i-1]*ans[i-1]+ans[i]*ans[i]) == 0.0)
            get_error("convlv: Deconvolving at reponse zero!");
         ans[i-1] = (fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
         ans[i] = (fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
      }
   }
   ans[2] = ans[n+1];
   realft(ans,n,-1);
   free_dvector(fft);
}

void correl(double *data1,double *data2,int n,double *ans)
{
   int nhalf,ndouble,i;
   double dum,*fft;

   if (!power_of_two(n)) get_error("correl: %d is not a power of two!",n);

   nhalf = n>>1;
   ndouble = n<<1;

   fft = dvector(ndouble);
   twofft(data1,data2,fft,ans,n);
   for (i=2;i<=n+2;i+=2)
   {
      ans[i-1] = (fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/nhalf;
      ans[i] = (fft[i]*dum-fft[i-1]*ans[i])/nhalf;
   }
   ans[2] = ans[n+1];
   realft(ans,n,-1);
   free_dvector(fft);
}

void fourn(float *data,int32_t *nn,int ndim,int isign,int pmod)
{
   int idim;
   int32_t i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
   int32_t ibit,k1,k2,n,nprev,nrem,ntot;
   float tempi,tempr;
   double theta,wi,wpi,wpr,wr,wtemp;
   int npass,pass,pold,pnew;
   double percent;

   if (pmod) setbuf(stderr,NULL);

   for (ntot=1,idim=1;idim<=ndim;idim++) ntot *= nn[idim];

   if (pmod)
   {
      nprev = 1;
      for (idim=ndim,npass=0;idim>=1;idim--)
      {
         n=nn[idim];
         nrem=ntot/(n*nprev);
         ip1=nprev<<1;
         ip2=ip1*n;
         ip3=ip2*nrem;
         ifp1=ip1;
         while (ifp1 < ip2)
         {
            ifp2=ifp1<<1;
            for (i3=1;i3<=ifp1;i3+=ip1)
               for (i1=i3;i1<=i3+ip1-2;i1+=2)
                  for (i2=i1;i2<=ip3;i2+=ifp2) npass++;
            ifp1=ifp2;
         }
         nprev *= n;
      }
      percent = 100.0/(double)npass;
      if (isign > 0)
         fprintf(stderr,"Fast Fourier Transform (%d-D):   0%%",ndim);
      else
         fprintf(stderr,"Inverse Fourier Transform (%d-D):   0%%",ndim);
   }

   nprev = 1;
   for (idim=ndim,pass=0,pold=-1;idim>=1;idim--)
   {
      n=nn[idim];
      nrem=ntot/(n*nprev);
      ip1=nprev<<1;
      ip2=ip1*n;
      ip3=ip2*nrem;
      i2rev=1;
      for (i2=1;i2<=ip2;i2+=ip1)
      {
         if (i2 < i2rev)
         {
            for (i1=i2;i1<=i2+ip1-2;i1+=2)
            {
               for (i3=i1;i3<=ip3;i3+=ip2)
               {
                  i3rev=i2rev+i3-i2;
                  FSWAP(data[i3],data[i3rev]);
                  FSWAP(data[i3+1],data[i3rev+1]);
               }
            }
         }
         ibit=ip2>>1;
         while (ibit>=ip1 && i2rev>ibit)
         {
            i2rev -= ibit;
            ibit = ibit>>1;
          }
         i2rev += ibit;
      }
      ifp1=ip1;
      while (ifp1 < ip2)
      {
         ifp2=ifp1<<1;
         theta = isign*TWOPI/(ifp2/ip1);
         wtemp=sin(0.5*theta);
         wpr=-2.0*wtemp*wtemp;
         wpi=sin(theta);
         wr=1.0;
         wi=0.0;
         for (i3=1;i3<=ifp1;i3+=ip1)
         {
            for (i1=i3;i1<=i3+ip1-2;i1+=2)
            {
               for (i2=i1;i2<=ip3;i2+=ifp2)
               {
                  k1=i2,k2=i2+ifp1;
                  tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
                  tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
                  data[k2]=data[k1]-tempr;
                  data[k2+1]=data[k1+1]-tempi;
                  data[k1] += tempr;
                  data[k1+1] += tempi;
                  if (pmod)
                  {
                     pass++;
                     pnew = nint((double)pass*percent);
                     if (pnew > pold)
                        fprintf(stderr,"\b\b\b\b%3d%%",pnew),pold=pnew;
                  }
               }
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
         }
         ifp1=ifp2;
      }
      nprev *= n;
   }
   if (pmod) fprintf(stderr,"\n");
}

void ccf2d(float *img,float *ref,int n1,int n2,float *ccf)
{
   int32_t nn[3];
   float *cimg,*cref;
   float cimg_re;
   int i,re,im;

   inform("Two dimensional cross-correlation start.");

   nn[1]=n1,nn[2]=n2;

   cimg = fvector(n1*n2*2);
   cref = fvector(n1*n2*2);

   for (i=1,re=1,im=2;i<=n1*n2;i++,re+=2,im+=2)
      cimg[re]=img[i],cref[re]=ref[i],cimg[im]=cref[im]=0.0;
   fourn(cimg,nn,2,1,1);
   fourn(cref,nn,2,1,1);

   for (i=1,re=1,im=2;i<=n1*n2;i++,re+=2,im+=2)
   {
      cimg_re = cimg[re];
      cimg[re] = cimg[re]*cref[re] + cimg[im]*cref[im];
      cimg[im] = cimg[im]*cref[re] - cimg_re*cref[im];
   }

   fourn(cimg,nn,2,-1,1);

   for (i=1,re=1;i<=n1*n2;i++,re+=2) ccf[i]=cimg[re]/(float)(n1*n2);

   free_fvector(cref);
   free_fvector(cimg);

   inform("Cross-correlation done.");
}

void linear_rebinning(double *g,int n,double xstart,double xstep,
  double *h,int m,double ustart,double ustep,double (*func)(double))
{
  int i,j,ia,ib;
  double ua,ub,xa,xb,ra,rb,s;

  clear_dvector(h,m);

  for (j=1;j<=m;j++)
  {
    ua = ustart + (j-1) * ustep;
    ub = ua + ustep;
    xa = (*func)(ua);
    xb = func(ub);
    ra = (xa - xstart) / xstep;
    rb = (xb - xstart) / xstep;
    if (rb <= ra) get_error("linear_rebinning: Argument out of order!");
    ia = (int)ceil(ra);
    ib = (int)ceil(rb);
    if (ia < 1 || ib > n) continue;
    if (ib == ia)
    {
      h[j] = (rb - ra) * g[ia];
    }
    else
    {
      s = (ia-ra) * g[ia] + (rb-ib+1) * g[ib];
      for (i=ia+1;i<=ib-1;i++) s += g[i];
      h[j] = s;
    }
  }
}

void rebin_spectrum(double *g,int n,double xstart,double xstep,
  double *h,int m,double ustart,double ustep,double (*func)(double))
{
  double *c,*q;
  int i,j;

  c = dvector(n+1);
  q = dvector(m+1);

  for (i=1;i<=n;i++) c[i] = 1.0;

  linear_rebinning(g,n,xstart,xstep,h,m,ustart,ustep,func);
  linear_rebinning(c,n,xstart,xstep,q,m,ustart,ustep,func);

  for (j=1;j<=m;j++) if (q[j] != 0.0) h[j] /= q[j];

  free_dvector(q);
  free_dvector(c);
}

void spline_rebinning(double *g,int n,double xstart,double xstep,
               double *h,int m,double ustart,double ustep,
               double (*func)(double))
{
   double *x,*g2,*g1,*f,*f2,*gc,*p;
   double c;
   int i;

   x = dvector(n+1);
   g2 = dvector(n);
   g1 = dvector(n);
   f = dvector(n+1);
   f2 = dvector(n+1);
   gc = dvector(n+1);
   p = dvector(m+1);

   for (i=1;i<=n;i++) x[i] = xstart +(i-1)*xstep;

   spline(x,g,n,1e30,1e30,g2);
   spline_derivative(x,g,g2,n,x,n,g1);

   for (i=1;i<=n+1;i++) x[i] = xstart +(i-1.5)*xstep;
   for (i=2,f[1]=0;i<=n+1;i++) f[i] = g1[i-1] + f[i-1];

   spline(x,f,n+1,1e30,1e30,f2);
   integrate_spline(x,f,f2,n+1,x,n+1,gc);
   for (i=1;i<=n;i++) gc[i] = (g[i] - gc[i]) / xstep;
   c = dselect((n+1)/2,n,gc);
   for (i=1;i<=n+1;i++) f[i] += c;

   for (i=1;i<=m+1;i++) p[i] = func(ustart+(i-1.5)*ustep);
   integrate_spline(x,f,f2,n+1,p,m+1,h);

   free_dvector(x);
   free_dvector(g2);
   free_dvector(g1);
   free_dvector(f);
   free_dvector(f2);
   free_dvector(gc);
   free_dvector(p);
}

void rebin_continuum(double xstep,
                     double *h,int m,double ustart,double ustep,
                     double (*func)(double))
{
   double *p,f;
   int i;

   p = dvector(m+1);

   f = 1.0/xstep;

   for (i=1;i<=m+1;i++) p[i] = func(ustart+(i-1.5)*ustep);
   for (i=1;i<=m;i++) h[i] = (p[i+1]-p[i])*f;

   free_dvector(p);
}

void rebin_normal(double *g,int n,double xstart,double xstep,
                  double *h,int m,double ustart,double ustep,
                  double (*func)(double))
{
   double *s;
   int i;

   s = dvector(m);

   spline_rebinning(g,n,xstart,xstep,h,m,ustart,ustep,func);
   rebin_continuum(xstep,s,m,ustart,ustep,func);

   for (i=1;i<=m;i++) h[i] /= s[i];

   free_dvector(s);
}

void correlate_arrays(double *xdat,double *ydat,double *zdat,int n,
                      double xstart,double ystart,double step,double *zstart)
{
   double *sig,*ref,*res;
   int n1,n2,n4,icen;

   for (n1=1;n1<n;n1=n1<<1);
   n2 = n1<<1;
   n4 = n2<<1;

   sig = dvector(n2);
   ref = dvector(n2);
   res = dvector(n4);

   clear_dvector(sig,n2);
   clear_dvector(ref,n2);

   copy_dvector(xdat,n,sig);
   copy_dvector(ydat,n,ref);

   correl(sig,ref,n2,res);

   if (ODD(n))
      icen = (1 + n) >> 1;
   else
      icen = n>>1;

   copy_dvector(res,n-icen+1,zdat+icen-1);
   copy_dvector(res+n2-icen+1,icen-1,zdat);

   *zstart = xstart - ystart - (icen-1)*step;

   free_dvector(sig);
   free_dvector(ref);
   free_dvector(res);
}

void statistics_double_array(double *a,int n,double *min,double *max,
                             double *mean,double *sig)
{
   int i;
   double sum;

   if (n < 2) get_error("statistics_double_array: Too few data values!");
   for (i=1,*min=*max=a[1],sum=0.0;i<=n;i++)
   {
      sum += a[i];
      if (a[i] < *min)
         *min = a[i];
      else
         if (a[i] > *max) *max = a[i];
   }
   *mean = sum/n;
   for (i=1,sum=0.0;i<=n;i++)
      sum += DSQR(a[i]-(*mean));
   *sig = sqrt(sum/(n-1));
}
