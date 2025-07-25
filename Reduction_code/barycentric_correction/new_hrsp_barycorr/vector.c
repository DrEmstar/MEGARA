#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "vector.h"

void set_vector(vec a,double x,double y,double z)
{
/* ------------------------------------------------------------------------
   Assign the three components (X, Y and Z) of a given vector.
   ------------------------------------------------------------------------ */

   a[0]=x;
   a[1]=y;
   a[2]=z;
}

void zero_vector(vec a)
{
/* ------------------------------------------------------------------------
   Define a zero-vector by setting all three components (X, Y and Z) to 0.
   ------------------------------------------------------------------------ */

  set_vector(a,0,0,0);
}

void set_ort(int axis,vec e)
{
/* ------------------------------------------------------------------------
   Create a unit vector (ort) along a given axis (X, Y or Z).
   ------------------------------------------------------------------------ */

   switch(axis)
   {
      case XAXIS: set_vector(e,1.0,0.0,0.0);
                  break;
      case YAXIS: set_vector(e,0.0,1.0,0.0);
                  break;
      case ZAXIS: set_vector(e,0.0,0.0,1.0);
                  break;
      default: get_error("set_ort: Invalid axis!");
   }
}

void copy_vector(vec a,vec b)
{
/* ------------------------------------------------------------------------
   Make a copy of a given vector, so that b = a.
   ------------------------------------------------------------------------ */

   memmove((void *)b,(void *)a,3*sizeof(double));
}

void opposite_vector(vec a,vec b)
{
/* ------------------------------------------------------------------------
   Invert the direction of a given vector, so that b = -a.
   ------------------------------------------------------------------------ */

   int i;

   for (i=0;i<3;i++) b[i] = -a[i];
}

double scalar_product(vec a,vec b)
{
/* ------------------------------------------------------------------------
   Calculate the vector dot product (scalar product) between two vectors.
   ------------------------------------------------------------------------ */

   int i;
   double w;

   for (i=0,w=0.0;i<3;i++) w += a[i]*b[i];
   return(w);
}

double vector_length(vec a)
{
/* ------------------------------------------------------------------------
   Return the length (magnitude) of a given vector.
   ------------------------------------------------------------------------ */

   return(sqrt(scalar_product(a,a)));
}

void vector_times_scalar(vec a,double u,vec b)
{
/* ------------------------------------------------------------------------
   Multiply a given vector by a scalar.
   ------------------------------------------------------------------------ */

   int i;

   for (i=0;i<3;i++) b[i]=u*a[i];
}

void vector_div_scalar(vec a,double u,vec b)
{
/* ------------------------------------------------------------------------
   Divide a given vector by a scalar.
   ------------------------------------------------------------------------ */

   int i;

   for (i=0;i<3;i++) b[i]=a[i]/u;
}

void unit_vector(vec a,vec b)
{
/* ------------------------------------------------------------------------
   Return the unit vector for a given vector. A zero vector on input
   will generate an error.
   ------------------------------------------------------------------------ */

   double u;

   u = vector_length(a);
   if (u == 0) get_error("unit_vector: Zero-vector detected!");
   vector_div_scalar(a,u,b);
}

void vector_sum(vec a,vec b,vec c)
{
/* ------------------------------------------------------------------------
   Add two vectors together.
   ------------------------------------------------------------------------ */

   int i;

   for (i=0;i<3;i++) c[i]=a[i]+b[i];
}

void vector_difference(vec a,vec b,vec c)
{
/* ------------------------------------------------------------------------
   Subtract two vectors.
   ------------------------------------------------------------------------ */

   int i;

   for (i=0;i<3;i++) c[i]=a[i]-b[i];
}

void vector_product(vec a,vec b,vec c)
{
/* ------------------------------------------------------------------------
   Calculate the vector cross product between two vectors.
   ------------------------------------------------------------------------ */

   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];
}

void copy_matrix(mat a,mat b)
{
/* ------------------------------------------------------------------------
   Make a copy of a given 3x3 matrix, so that b = a.
   ------------------------------------------------------------------------ */

   memmove((void *)b,(void *)a,9*sizeof(double));
}

void matrix_times_scalar(mat a,double u,mat b)
{
/* ------------------------------------------------------------------------
   Multiply a given 3x3 matrix by a scalar.
   ------------------------------------------------------------------------ */

   int i,j;

   for (i=0;i<3;i++) for(j=0;j<3;j++) b[i][j] = u*a[i][j];
}

void matrix_times_vector(mat g,vec a,vec b)
{
/* ------------------------------------------------------------------------
   Multiply a gven 3-D vector by a 3x3 matrix. If the matrix describes
   a rotation of coordinates, then the product will transform the input
   vector from the old coordinates to the new ones.
   ------------------------------------------------------------------------ */

   int i,j;

   for (i=0;i<3;i++)
      for (j=0,b[i]=0;j<3;j++) b[i] += g[i][j]*a[j];
}

void matrix_product(mat a,mat b,mat c)
{
/* ------------------------------------------------------------------------
   Multiply two matrices.
   ------------------------------------------------------------------------ */

   int i,j,k;

   for (i=0;i<3;i++)
      for (j=0;j<3;j++)
         for (k=c[i][j]=0;k<3;k++) c[i][j] += a[i][k]*b[k][j];
}

void matrix_triple_product(mat a,mat b,mat c,mat x)
{
/* ------------------------------------------------------------------------
   Multiply three matrices (matrix triple product).
   ------------------------------------------------------------------------ */

  mat u;

  matrix_product(b,c,u);
  matrix_product(a,u,x);
}

void matrix_quad_product(mat a,mat b,mat c,mat d,mat x)
{
/* ------------------------------------------------------------------------
   Multiply four matrices (matrix quadruple product).
   ------------------------------------------------------------------------ */

  mat u,v;

  matrix_product(c,d,u);
  matrix_product(b,u,v);
  matrix_product(a,v,x);
}

void matrix_quint_product(mat a,mat b,mat c,mat d,mat e,mat x)
{
/* ------------------------------------------------------------------------
   Multiply five matrices (matrix quintuple product).
   ------------------------------------------------------------------------ */

  mat u,v,w;

  matrix_product(d,e,u);
  matrix_product(c,u,v);
  matrix_product(b,v,w);
  matrix_product(a,w,x);
}

void matrix_row(mat a,int i,vec b)
{
/* ------------------------------------------------------------------------
   Extract a row-vector from a given matrix.
   ------------------------------------------------------------------------ */

  int j;

  for (j=0;j<3;j++) b[j] = a[i][j];
}

void matrix_col(mat a,int j,vec b)
{
/* ------------------------------------------------------------------------
   Extract a column-vector from a given matrix.
   ------------------------------------------------------------------------ */

  int i;

  for (i=0;i<3;i++) b[i] = a[i][j];
}

void transpose_matrix(mat a,mat b)
{
/* ------------------------------------------------------------------------
   Transpose a given matrix 3x3 by swapping its columns and rows.
   ------------------------------------------------------------------------ */

   int i,j;

   for (i=0;i<3;i++) for(j=0;j<3;j++) b[i][j] = a[j][i];
}

void rotate(int axis,double angle,mat rot)
{
/* ------------------------------------------------------------------------
   Rotate the XYZ reference frame around a coordinate axis (X, Y or Z) by
   a given angle and return the corresponding rotation matrix.
   The rotation matrix can be used to transform any given vector from the
   old coordinates to the new ones:

      [new_vector] = [rotation_matrix] x [old_vector]
   ------------------------------------------------------------------------ */

   int i,j,k;

   i = axis;
   j = (i + 1) % 3;
   k = (j + 1) % 3;

   rot[i][i] = 1.0;
   rot[i][j] = rot[i][k] = rot[j][i] = rot[k][i] = 0.0;
   rot[j][j] = rot[k][k] = cos(angle);
   rot[j][k] = sin(angle);
   rot[k][j] = -rot[j][k];
}

void sphere_to_xyz(double a,double h,double dist,vec r)
{
/* ------------------------------------------------------------------------
   Convert a given set of spherical coordinates (azimuth, elevation and
   distance) to a vector in XYZ space.
   ------------------------------------------------------------------------ */

   double rho;

   rho = dist * cos(h);

   r[0] = rho * cos(a);
   r[1] = rho * sin(a);
   r[2] = dist * sin(h);
}

void xyz_to_sphere(vec r,double *a,double *h,double *dist)
{
/* ------------------------------------------------------------------------
   Convert a given set of X, Y and Z coordinates (a 3-D vector) into
   spherical coordinates: azimuth, elevation and distance. The transform
   works for all input vectors. A zero input vector will be converted to
   zero azimuth, zero elevation and zero distance. An input vector along
   the Z-axis will have a zero azimuth.
   ------------------------------------------------------------------------ */

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


void arbitrary_rotation(vec axis,double psi,mat rot)
{
/* ------------------------------------------------------------------------
   Rotate the XYZ reference frame around an arbitrary axis (vector in 3-D
   space) by a given angle and return the corresponding rotation matrix.
   The rotation matrix can be used to transform any given vector from the
   old coordinates to the new ones:

      [new_vector] = [rotation_matrix] x [old_vector]
   ------------------------------------------------------------------------ */

   double alp,del,rho;
   mat arot,ainv,drot,dinv,rpsi;

   xyz_to_sphere(axis,&alp,&del,&rho);

   rotate(ZAXIS,alp,arot);
   rotate(ZAXIS,-alp,ainv);

   rotate(YAXIS,HALFPI-del,drot);
   rotate(YAXIS,del-HALFPI,dinv);

   rotate(ZAXIS,psi,rpsi);

   matrix_quint_product(ainv,dinv,rpsi,drot,arot,rot);
}
