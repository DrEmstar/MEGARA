#ifndef VECTOR_DEFINED

   #define VECTOR_DEFINED

   #define XAXIS 0
   #define YAXIS 1
   #define ZAXIS 2

   typedef double vec[3];
   typedef double mat[3][3];

#endif

void set_vector(vec a,double x,double y,double z);
void zero_vector(vec a);
void set_ort(int axis,vec e);
void copy_vector(vec a,vec b);
void opposite_vector(vec a,vec b);
double scalar_product(vec a,vec b);
double vector_length(vec a);
void vector_times_scalar(vec a,double u,vec b);
void vector_div_scalar(vec a,double u,vec b);
void unit_vector(vec a,vec b);
void vector_sum(vec a,vec b,vec c);
void vector_difference(vec a,vec b,vec c);
void vector_product(vec a,vec b,vec c);
void copy_matrix(mat a,mat b);
void matrix_times_scalar(mat a,double u,mat b);
void matrix_times_vector(mat g,vec a,vec b);
void matrix_product(mat a,mat b,mat c);
void matrix_triple_product(mat a,mat b,mat c,mat x);
void matrix_quad_product(mat a,mat b,mat c,mat d,mat x);
void matrix_quint_product(mat a,mat b,mat c,mat d,mat e,mat x);
void matrix_row(mat a,int i,vec b);
void matrix_col(mat a,int j,vec b);
void transpose_matrix(mat a,mat b);
void rotate(int axis,double angle,mat rot);
void sphere_to_xyz(double a,double h,double dist,vec r);
void xyz_to_sphere(vec r,double *a,double *h,double *dist);
void arbitrary_rotation(vec axis,double psi,mat rot);
