#ifndef MATRIXUTILH
#define MATRIXUTILH

void double_mat_dot(double*, double*, double*);

void double_vec_cross(VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*);

double double_vec_dot(VECTOR_DOUBLE*, VECTOR_DOUBLE*);

void double_vec_diff(VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*);

int check_inverse1(double*, double*);

void map_to_wigner(CRYSTAL*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*);

void map_to_brillouin(CRYSTAL*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_INT*, FILES);

void getcart(VECTOR_INT*, VECTOR_DOUBLE*, int*, CRYSTAL*);


/*

void double_vec_sum(VECTOR_DOUBLE *vec1, VECTOR_DOUBLE *vec2, VECTOR_DOUBLE *vec3);

void rotate_vector2(int *symmoperator, VECTOR_DOUBLE *vector_in, VECTOR_DOUBLE *vectorR);
*/

void rotate_vector3(double*, VECTOR_DOUBLE*, VECTOR_DOUBLE*);

void rotate_vector_int(int*, VECTOR_INT*, VECTOR_INT*);

void rotate_vector_int_latt(int*, int[3], int, VECTOR_INT*, VECTOR_INT*, CRYSTAL*, JOB_PARAM*, FILES);

/*

void rotate_vector_int_transpose(int *, VECTOR_INT *, VECTOR_INT *);
*/

int check_vec(VECTOR_DOUBLE*, VECTOR_DOUBLE*);

int check_mat(double*, double*) ;

#endif
