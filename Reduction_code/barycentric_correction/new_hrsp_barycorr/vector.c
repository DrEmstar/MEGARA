#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "vector.h"

int next_index(int i)
{
   if (i < 2)
      return(i+1);
   else
      return(0);
}

void set_vector(vec a,double x,double y,double z)
{
   a[0]=x, a[1]=y, a[2]=z;
}

void copy_vector(vec a,vec b)
{
   memcpy((void *)b,(void *)a,3*sizeof(double));
}

double scalar_product(vec a,vec b)
{
   int i;
   double w;

   for (i=0,w=0.0;i<3;i++) w += a[i]*b[i];
   return(w);
}

double vector_length(vec a)
{
   return(sqrt(scalar_product(a,a)));
}

void vector_times_scalar(vec a,double u,vec b)
{
   int i;

   for (i=0;i<3;i++) b[i]=u*a[i];
}

void unit_vector(vec a,vec b)
{
   double u;

   u = vector_length(a);
   if (a == 0) get_error("unit_vector: Zero-vector detected!");
   vector_times_scalar(a,1.0/u,b);
}

void vector_sum(vec a,vec b,vec c)
{
   int i;

   for (i=0;i<3;i++) c[i]=a[i]+b[i];
}

void vector_product(vec a,vec b,vec c)
{
   int i,j,k;

   for (i=0;i<3;i++)
   {
      j = next_index(i);
      k = next_index(j);
      c[i] = a[j]*b[k] - a[k]*b[j];
   }
}

void copy_matrix(mat a,mat b)
{
   memcpy((void *)b,(void *)a,9*sizeof(double));
}

void matrix_times_scalar(mat a,double u,mat b)
{
   int i,j;

   for (i=0;i<3;i++) for(j=0;j<3;j++) b[i][j] = u*a[i][j];
}

void matrix_times_vector(mat g,vec a,vec b)
{
   int i,j;

   for (i=0;i<3;i++)
      for (j=0,b[i]=0;j<3;j++) b[i] += g[i][j]*a[j];
}

void matrix_product(mat a,mat b,mat c)
{
   int i,j,k;

   for (i=0;i<3;i++)
      for (j=0;j<3;j++)
         for (k=c[i][j]=0;k<3;k++) c[i][j] += a[i][k]*b[k][j];
}

void transpose_matrix(mat a,mat b)
{
   int i,j;

   for (i=0;i<3;i++) for(j=0;j<3;j++) b[i][j] = a[j][i];
}

double cofactor(mat a,int i,int j)
{
   int m,n,p,q;

   m = next_index(i);
   n = next_index(m);
   p = next_index(j);
   q = next_index(p);

   return(a[m][p]*a[n][q]-a[n][p]*a[m][q]);
}

double determinant(mat a)
{
   int i;
   double det;

   for (i=0,det=0;i<3;i++) det += a[i][0]*cofactor(a,i,0);

   return(det);
}

void inverse_matrix(mat a,mat b)
{
   int i,j;
   mat c;
   double det;

   det = determinant(a);
   if (det == 0.0) get_error("inverse_matrix: Singular matrix!");
   for (i=0;i<3;i++) for(j=0;j<3;j++) b[i][j] = cofactor(a,i,j);
   transpose_matrix(b,c);
   matrix_times_scalar(c,1/det,b);
}

void rotate(int axis,double angle,mat rot)
{
   int i,j,k;

   i = axis;
   j = next_index(i);
   k = next_index(j);

   rot[i][i] = 1.0;
   rot[i][j] = rot[i][k] = rot[j][i] = rot[k][i] = 0.0;
   rot[j][j] = rot[k][k] = cos(angle);
   rot[j][k] = sin(angle);
   rot[k][j] = -rot[j][k];
}

void sphere_to_xyz(double a,double h,double dist,vec r)
{
   r[0] = dist*cos(h)*cos(a);
   r[1] = dist*cos(h)*sin(a);
   r[2] = dist*sin(h);
}

void xyz_to_sphere(vec r,double *a,double *h,double *dist)
{
   double rho;

   if (r[0]==0 && r[1]==0)
      *a = 0.0;
   else
      *a = atan2(r[1],r[0]);
   if (*a < 0) *a += TWOPI;

   rho = sqrt(r[0]*r[0]+r[1]*r[1]);
   if (rho==0 && r[2]==0)
      *h = 0.0;
   else
      *h = atan2(r[2],rho);

   *dist = vector_length(r);
}
