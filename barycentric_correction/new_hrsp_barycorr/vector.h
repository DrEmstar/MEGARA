#ifndef VECTOR_DEFINED

   #define VECTOR_DEFINED

   #define XAXIS 0
   #define YAXIS 1
   #define ZAXIS 2

   typedef double vec[3];
   typedef double mat[3][3];

#endif

/* ------------------------ Function Declaration ------------------------ */
int next_index(int i);
void set_vector(vec a,double x,double y,double z);
void copy_vector(vec a,vec b);
double scalar_product(vec a,vec b);
double vector_length(vec a);
void vector_times_scalar(vec a,double u,vec b);
void unit_vector(vec a,vec b);
void vector_sum(vec a,vec b,vec c);
void vector_product(vec a,vec b,vec c);
void copy_matrix(mat a,mat b);
void matrix_times_scalar(mat a,double u,mat b);
void matrix_times_vector(mat g,vec a,vec b);
void matrix_product(mat a,mat b,mat c);
void transpose_matrix(mat a,mat b);
double cofactor(mat a,int i,int j);
double determinant(mat a);
void inverse_matrix(mat a,mat b);
void rotate(int axis,double angle,mat rot);
void sphere_to_xyz(double a,double h,double dist,vec r);
void xyz_to_sphere(vec r,double *a,double *h,double *dist);
