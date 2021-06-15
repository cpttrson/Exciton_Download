/*
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <fstream>
#include <mpi.h>
#include "mycomplex.h"
#include "mylogical.h"
#include "conversion_factors.h"
*/
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"

using namespace std;

void double_mat_dot(double *mat1, double *mat2, double *mat3)

{
  *(mat3 + 0) = *(mat1 + 0) * *(mat2 + 0) + *(mat1 + 1) * *(mat2 + 3) + *(mat1 + 2) * *(mat2 + 6);
  *(mat3 + 1) = *(mat1 + 0) * *(mat2 + 1) + *(mat1 + 1) * *(mat2 + 4) + *(mat1 + 2) * *(mat2 + 7);
  *(mat3 + 2) = *(mat1 + 0) * *(mat2 + 2) + *(mat1 + 1) * *(mat2 + 5) + *(mat1 + 2) * *(mat2 + 8);
  *(mat3 + 3) = *(mat1 + 3) * *(mat2 + 0) + *(mat1 + 4) * *(mat2 + 3) + *(mat1 + 5) * *(mat2 + 6);
  *(mat3 + 4) = *(mat1 + 3) * *(mat2 + 1) + *(mat1 + 4) * *(mat2 + 4) + *(mat1 + 5) * *(mat2 + 7);
  *(mat3 + 5) = *(mat1 + 3) * *(mat2 + 2) + *(mat1 + 4) * *(mat2 + 5) + *(mat1 + 5) * *(mat2 + 8);
  *(mat3 + 6) = *(mat1 + 6) * *(mat2 + 0) + *(mat1 + 7) * *(mat2 + 3) + *(mat1 + 8) * *(mat2 + 6);
  *(mat3 + 7) = *(mat1 + 6) * *(mat2 + 1) + *(mat1 + 7) * *(mat2 + 4) + *(mat1 + 8) * *(mat2 + 7);
  *(mat3 + 8) = *(mat1 + 6) * *(mat2 + 2) + *(mat1 + 7) * *(mat2 + 5) + *(mat1 + 8) * *(mat2 + 8);
}

void double_vec_cross(VECTOR_DOUBLE *a, VECTOR_DOUBLE *b, VECTOR_DOUBLE *n)

{

  n->comp1 = a->comp2 * b->comp3 - b->comp2 * a->comp3;
  n->comp2 = a->comp3 * b->comp1 - b->comp3 * a->comp1;
  n->comp3 = a->comp1 * b->comp2 - b->comp1 * a->comp2;

}

double double_vec_dot(VECTOR_DOUBLE *vec1, VECTOR_DOUBLE *vec2)

{

  int i;
  double product = k_zero;

  product = vec1->comp1 * vec2->comp1 + vec1->comp2 * vec2->comp2 + vec1->comp3 * vec2->comp3;

  return product;

}

void double_vec_diff(VECTOR_DOUBLE *vec1, VECTOR_DOUBLE *vec2, VECTOR_DOUBLE *vec3)

{

  vec3->comp1 = vec1->comp1 - vec2->comp1;
  vec3->comp2 = vec1->comp2 - vec2->comp2;
  vec3->comp3 = vec1->comp3 - vec2->comp3;

}

int check_mat(double *mat1, double *mat2)

{
  double mag = 0.001;
  if (fabs(*(mat1 + 0) - *(mat2 + 0)) < mag && fabs(*(mat1 + 1) - *(mat2 + 1)) < mag && fabs(*(mat1 + 2) - *(mat2 + 2))
      < mag && fabs(*(mat1 + 3) - *(mat2 + 3)) < mag && fabs(*(mat1 + 4) - *(mat2 + 4)) < mag && fabs(*(mat1 + 5)
      - *(mat2 + 5)) < mag && fabs(*(mat1 + 6) - *(mat2 + 6)) < mag && fabs(*(mat1 + 7) - *(mat2 + 7)) < mag && fabs(
      *(mat1 + 8) - *(mat2 + 8)) < mag)
    return 1;
  else
    return 0;

}

int check_inverse1(double *p_irr, double *q_irr)

{

  double A[3][3], eps = 0.0000001;
  int n = 0;

  A[0][0] = *(p_irr) * *(q_irr) + *(p_irr + 1) * *(q_irr + 3) + *(p_irr + 2) * *(q_irr + 6);
  A[0][1] = *(p_irr + 3) * *(q_irr) + *(p_irr + 4) * *(q_irr + 3) + *(p_irr + 5) * *(q_irr + 6);
  A[0][2] = *(p_irr + 6) * *(q_irr) + *(p_irr + 7) * *(q_irr + 3) + *(p_irr + 8) * *(q_irr + 6);
  A[1][0] = *(p_irr) * *(q_irr + 1) + *(p_irr + 1) * *(q_irr + 4) + *(p_irr + 2) * *(q_irr + 7);
  A[1][1] = *(p_irr + 3) * *(q_irr + 1) + *(p_irr + 4) * *(q_irr + 4) + *(p_irr + 5) * *(q_irr + 7);
  A[1][2] = *(p_irr + 6) * *(q_irr + 1) + *(p_irr + 7) * *(q_irr + 4) + *(p_irr + 8) * *(q_irr + 7);
  A[2][0] = *(p_irr) * *(q_irr + 2) + *(p_irr + 1) * *(q_irr + 5) + *(p_irr + 2) * *(q_irr + 8);
  A[2][1] = *(p_irr + 3) * *(q_irr + 2) + *(p_irr + 4) * *(q_irr + 5) + *(p_irr + 5) * *(q_irr + 8);
  A[2][2] = *(p_irr + 6) * *(q_irr + 2) + *(p_irr + 7) * *(q_irr + 5) + *(p_irr + 8) * *(q_irr + 8);

  if (fabs(A[0][0] - 1.0) < eps && fabs(A[0][1]) < eps && fabs(A[0][2]) < eps && fabs(A[1][0]) < eps && fabs(A[1][1]
      - 1.0) < eps && fabs(A[1][2]) < eps && fabs(A[2][0]) < eps && fabs(A[2][1]) < eps && fabs(A[2][2] - 1.0) < eps) {
    n = 1;
  }

  return n;

}

void map_to_wigner(CRYSTAL *crystal, VECTOR_DOUBLE *vec, VECTOR_DOUBLE *mapped_vec, VECTOR_DOUBLE *map_vec)

{

  double mag1, mag2, mag3;
  double m, n, p;

  mag1 = (vec->comp1 * crystal->reciprocal_cell[0].comp1 + vec->comp2 * crystal->reciprocal_cell[0].comp2 + vec->comp3 * crystal->reciprocal_cell[0].comp3) / 2.0 / pi;
  mag2 = (vec->comp1 * crystal->reciprocal_cell[1].comp1 + vec->comp2 * crystal->reciprocal_cell[1].comp2 + vec->comp3 * crystal->reciprocal_cell[1].comp3) / 2.0 / pi;
  mag3 = (vec->comp1 * crystal->reciprocal_cell[2].comp1 + vec->comp2 * crystal->reciprocal_cell[2].comp2 + vec->comp3 * crystal->reciprocal_cell[2].comp3) / 2.0 / pi;
  m = floor(mag1 + double(0.49999));
  n = floor(mag2 + double(0.49999));
  p = floor(mag3 + double(0.49999));
  if (mag1 + double(0.49999) - m > 0.999999)
    m += 1.0;
  if (mag2 + double(0.49999) - n > 0.999999)
    n += 1.0;
  if (mag3 + double(0.49999) - p > 0.999999)
    p += 1.0;
  map_vec->comp1 = -(m * crystal->primitive_cell[0].comp1 + n * crystal->primitive_cell[1].comp1 + p * crystal->primitive_cell[2].comp1);
  map_vec->comp2 = -(m * crystal->primitive_cell[0].comp2 + n * crystal->primitive_cell[1].comp2 + p * crystal->primitive_cell[2].comp2);
  map_vec->comp3 = -(m * crystal->primitive_cell[0].comp3 + n * crystal->primitive_cell[1].comp3 + p * crystal->primitive_cell[2].comp3);
  mapped_vec->comp1 = vec->comp1 + map_vec->comp1;
  mapped_vec->comp2 = vec->comp2 + map_vec->comp2;
  mapped_vec->comp3 = vec->comp3 + map_vec->comp3;

}

void map_to_brillouin(CRYSTAL *crystal, VECTOR_DOUBLE *vec, VECTOR_DOUBLE *mapped_vec, VECTOR_DOUBLE *map_vec, VECTOR_INT *map, FILES file)

{

  int i, j, k;
  double x, y, z, mag1, mag2, mag3, magk, magG, magtmp;
  double m, n, p;
  VECTOR_DOUBLE tmp_vec;

  mag1 = (vec->comp1 * crystal->primitive_cell[0].comp1 + vec->comp2 * crystal->primitive_cell[0].comp2 + vec->comp3 * crystal->primitive_cell[0].comp3) / 2.0 / pi;
  mag2 = (vec->comp1 * crystal->primitive_cell[1].comp1 + vec->comp2 * crystal->primitive_cell[1].comp2 + vec->comp3 * crystal->primitive_cell[1].comp3) / 2.0 / pi;
  mag3 = (vec->comp1 * crystal->primitive_cell[2].comp1 + vec->comp2 * crystal->primitive_cell[2].comp2 + vec->comp3 * crystal->primitive_cell[2].comp3) / 2.0 / pi;
  m = floor(mag1 + double(0.49999));
  n = floor(mag2 + double(0.49999));
  p = floor(mag3 + double(0.49999));
  if (mag1 + double(0.49999) - m > 0.999999)
    m += 1.0;
  if (mag2 + double(0.49999) - n > 0.999999)
    n += 1.0;
  if (mag3 + double(0.49999) - p > 0.999999)
    p += 1.0;
  map_vec->comp1 = -(m * crystal->reciprocal_cell[0].comp1 + n * crystal->reciprocal_cell[1].comp1 + p * crystal->reciprocal_cell[2].comp1);
  map_vec->comp2 = -(m * crystal->reciprocal_cell[0].comp2 + n * crystal->reciprocal_cell[1].comp2 + p * crystal->reciprocal_cell[2].comp2);
  map_vec->comp3 = -(m * crystal->reciprocal_cell[0].comp3 + n * crystal->reciprocal_cell[1].comp3 + p * crystal->reciprocal_cell[2].comp3);
  mapped_vec->comp1 = vec->comp1 + map_vec->comp1;
  mapped_vec->comp2 = vec->comp2 + map_vec->comp2;
  mapped_vec->comp3 = vec->comp3 + map_vec->comp3;
  
  magk = sqrt(mapped_vec->comp1 * mapped_vec->comp1 + mapped_vec->comp2 * mapped_vec->comp2 + mapped_vec->comp3 * mapped_vec->comp3);

  map->comp1 = 0;
  map->comp2 = 0;
  map->comp3 = 0;
  for (i = -3; i < 4; i++) {
    for (j = -3; j < 4; j++) {
      for (k = -3; k < 4; k++) {
        if (i == 0 && j == 0 && k == 0) continue;
        x = i * crystal->reciprocal_cell[0].comp1 + j * crystal->reciprocal_cell[1].comp1 + k * crystal->reciprocal_cell[2].comp1;
        y = i * crystal->reciprocal_cell[0].comp2 + j * crystal->reciprocal_cell[1].comp2 + k * crystal->reciprocal_cell[2].comp2;
        z = i * crystal->reciprocal_cell[0].comp3 + j * crystal->reciprocal_cell[1].comp3 + k * crystal->reciprocal_cell[2].comp3;
        tmp_vec.comp1 = mapped_vec->comp1 - x;
        tmp_vec.comp2 = mapped_vec->comp2 - y;
        tmp_vec.comp3 = mapped_vec->comp3 - z;
        magG = sqrt(x * x + y * y + z * z);
        magtmp = sqrt(tmp_vec.comp1 * tmp_vec.comp1 + tmp_vec.comp2 * tmp_vec.comp2 + tmp_vec.comp3 * tmp_vec.comp3);
        if (magtmp  < magk) {
          magk = magtmp;
          mapped_vec->comp1 = tmp_vec.comp1;
          mapped_vec->comp2 = tmp_vec.comp2;
          mapped_vec->comp3 = tmp_vec.comp3;
          map->comp1 = i;
          map->comp2 = j;
          map->comp3 = k;
          //fprintf(file.out,"%3d %3d %3d not mapped %f %f   %f %f %f   %f %f %f\n",i,j,k,magG,magk,\
          mapped_vec->comp1, mapped_vec->comp2, mapped_vec->comp3,x,y,z);
         }
       }
      }
     }

}


/*
void double_vec_sum(VECTOR_DOUBLE *vec1, VECTOR_DOUBLE *vec2, VECTOR_DOUBLE *vec3)

{

  vec3->comp1 = vec1->comp1 + vec2->comp1;
  vec3->comp2 = vec1->comp2 + vec2->comp2;
  vec3->comp3 = vec1->comp3 + vec2->comp3;

}

void rotate_vector2(int *symmoperator, VECTOR_DOUBLE *vector_in, VECTOR_DOUBLE *vectorR)

{

  int *p_irr;

  p_irr = symmoperator;
  vectorR->comp1 = (*p_irr) * vector_in->comp1 + (*(p_irr + 1)) * vector_in->comp2 +\
 (*(p_irr + 2))
      * vector_in->comp3;
  vectorR->comp2 = (*(p_irr + 3)) * vector_in->comp1 + (*(p_irr + 4)) * vector_in->comp2 +\
 (*(p_irr + 5))
      * vector_in->comp3;
  vectorR->comp3 = (*(p_irr + 6)) * vector_in->comp1 + (*(p_irr + 7)) * vector_in->comp2 +\
 (*(p_irr + 8))
      * vector_in->comp3;

  return;

}

*/

void getcart(VECTOR_INT *vec, VECTOR_DOUBLE *vec2, int *is, CRYSTAL *crystal)

{
     vec2->comp1 = vec->comp1 * crystal->reciprocal_cell[0].comp1 / (double)is[0] +
                   vec->comp2 * crystal->reciprocal_cell[1].comp1 / (double)is[1] +
                   vec->comp3 * crystal->reciprocal_cell[2].comp1 / (double)is[2] ;
     vec2->comp2 = vec->comp1 * crystal->reciprocal_cell[0].comp2 / (double)is[0] +
                   vec->comp2 * crystal->reciprocal_cell[1].comp2 / (double)is[1] +
                   vec->comp3 * crystal->reciprocal_cell[2].comp2 / (double)is[2] ;
     vec2->comp3 = vec->comp1 * crystal->reciprocal_cell[0].comp3 / (double)is[0] +
                   vec->comp2 * crystal->reciprocal_cell[1].comp3 / (double)is[1] +
                   vec->comp3 * crystal->reciprocal_cell[2].comp3 / (double)is[2] ;

   return ;

}

void rotate_vector3(double *symmoperator, VECTOR_DOUBLE *vector_in, VECTOR_DOUBLE *vectorR)

{

  double *p_irr;

  p_irr = symmoperator;
  vectorR->comp1 = (*p_irr) * vector_in->comp1 + (*(p_irr + 1)) * vector_in->comp2 +\
 (*(p_irr + 2))
      * vector_in->comp3;
  vectorR->comp2 = (*(p_irr + 3)) * vector_in->comp1 + (*(p_irr + 4)) * vector_in->comp2 +\
 (*(p_irr + 5))
      * vector_in->comp3;
  vectorR->comp3 = (*(p_irr + 6)) * vector_in->comp1 + (*(p_irr + 7)) * vector_in->comp2 +\
 (*(p_irr + 8))
      * vector_in->comp3;

  return;
}

void rotate_vector_int(int *symmoperator, VECTOR_INT *vector_in, VECTOR_INT *vectorR)

{

  int *p_irr;

  p_irr = symmoperator;

  vectorR->comp1 = (*p_irr)       * vector_in->comp1 + (*(p_irr + 1)) * vector_in->comp2 + (*(p_irr + 2)) * vector_in->comp3;
  vectorR->comp2 = (*(p_irr + 3)) * vector_in->comp1 + (*(p_irr + 4)) * vector_in->comp2 + (*(p_irr + 5)) * vector_in->comp3;
  vectorR->comp3 = (*(p_irr + 6)) * vector_in->comp1 + (*(p_irr + 7)) * vector_in->comp2 + (*(p_irr + 8)) * vector_in->comp3;

  return;

}

void rotate_vector_int_latt(int *symmoperator, int is[3], int trs, VECTOR_INT *vector_in, VECTOR_INT *vectorR, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int *p_irr;

  p_irr = symmoperator;

  vectorR->comp1 = (*p_irr)       * vector_in->comp1 + (*(p_irr + 1)) * vector_in->comp2 + (*(p_irr + 2)) * vector_in->comp3;
  vectorR->comp2 = (*(p_irr + 3)) * vector_in->comp1 + (*(p_irr + 4)) * vector_in->comp2 + (*(p_irr + 5)) * vector_in->comp3;
  vectorR->comp3 = (*(p_irr + 6)) * vector_in->comp1 + (*(p_irr + 7)) * vector_in->comp2 + (*(p_irr + 8)) * vector_in->comp3;

  if (trs == 1) { vectorR->comp1 *= -1; vectorR->comp2 *= -1; vectorR->comp3 *= -1; }
              
  switch (crystal->type[0]) {

  case 'C':

  vectorR->comp1 -= (vectorR->comp1 / is[0]) * is[0];
  vectorR->comp2 -= (vectorR->comp2 / is[1]) * is[1];
  vectorR->comp3 -= (vectorR->comp3 / is[2]) * is[2];
  if (vectorR->comp1 < 0) vectorR->comp1 += is[0];
  if (vectorR->comp2 < 0) vectorR->comp2 += is[1];
  if (vectorR->comp3 < 0) vectorR->comp3 += is[2];
  if (vectorR->comp1 == is[0]) vectorR->comp1 = 0;
  if (vectorR->comp2 == is[1]) vectorR->comp2 = 0;
  if (vectorR->comp3 == is[2]) vectorR->comp3 = 0;

  break;

  case 'S':

  vectorR->comp1 -= (vectorR->comp1 / is[0]) * is[0];
  vectorR->comp2 -= (vectorR->comp2 / is[1]) * is[1];
  if (vectorR->comp1 < 0) vectorR->comp1 += is[0];
  if (vectorR->comp2 < 0) vectorR->comp2 += is[1];
  if (vectorR->comp1 == is[0]) vectorR->comp1 = 0;
  if (vectorR->comp2 == is[1]) vectorR->comp2 = 0;

  break;

  case 'P':

  vectorR->comp3 -= (vectorR->comp3 / is[2]) * is[2];
  if (vectorR->comp3 < 0) vectorR->comp3 += is[2];
  if (vectorR->comp3 == is[2]) vectorR->comp3 = 0;

  break;

  default:

  if (job->taskid == 0) fprintf(file.out,"Routine rotate_vector_int_latt called with incorrect Crystal->type\n"); 
  MPI_Finalize();
  exit(1);

  break;

 } // close switch

  return;

}

void rotate_vector_int_latt1(int *symmoperator, int is[3], int trs, VECTOR_INT *vector_in, VECTOR_INT *vectorR, VECTOR_INT *vectorM, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int *p_irr;

  p_irr = symmoperator;

  vectorR->comp1 = (*p_irr)       * vector_in->comp1 + (*(p_irr + 1)) * vector_in->comp2 + (*(p_irr + 2)) * vector_in->comp3;
  vectorR->comp2 = (*(p_irr + 3)) * vector_in->comp1 + (*(p_irr + 4)) * vector_in->comp2 + (*(p_irr + 5)) * vector_in->comp3;
  vectorR->comp3 = (*(p_irr + 6)) * vector_in->comp1 + (*(p_irr + 7)) * vector_in->comp2 + (*(p_irr + 8)) * vector_in->comp3;

  if (trs == 1) { vectorR->comp1 *= -1; vectorR->comp2 *= -1; vectorR->comp3 *= -1; }
              
  vectorM->comp1 = 0;
  vectorM->comp2 = 0;
  vectorM->comp3 = 0;

  switch (crystal->type[0]) {

  case 'C':

  vectorR->comp1 -= (vectorR->comp1 / is[0]) * is[0];
  vectorR->comp2 -= (vectorR->comp2 / is[1]) * is[1];
  vectorR->comp3 -= (vectorR->comp3 / is[2]) * is[2];
  if (vectorR->comp1 < 0) { vectorR->comp1 += is[0]; vectorM->comp1 =  is[0]; }
  if (vectorR->comp2 < 0) { vectorR->comp2 += is[1]; vectorM->comp2 =  is[1]; }
  if (vectorR->comp3 < 0) { vectorR->comp3 += is[2]; vectorM->comp3 =  is[2]; }
  if (vectorR->comp1 == is[0]) { vectorR->comp1 = 0; vectorM->comp1 = -is[0]; }
  if (vectorR->comp2 == is[1]) { vectorR->comp2 = 0; vectorM->comp2 = -is[1]; }
  if (vectorR->comp3 == is[2]) { vectorR->comp3 = 0; vectorM->comp3 = -is[2]; }

  break;

  case 'S':

  vectorR->comp1 -= (vectorR->comp1 / is[0]) * is[0];
  vectorR->comp2 -= (vectorR->comp2 / is[1]) * is[1];
  if (vectorR->comp1 < 0) { vectorR->comp1 += is[0]; vectorM->comp1 =  is[0]; }
  if (vectorR->comp2 < 0) { vectorR->comp2 += is[1]; vectorM->comp2 =  is[1]; }
  if (vectorR->comp1 == is[0]) { vectorR->comp1 = 0; vectorM->comp1 = -is[0]; }
  if (vectorR->comp2 == is[1]) { vectorR->comp2 = 0; vectorM->comp2 = -is[1]; }

  break;

  case 'P':

  vectorR->comp3 -= (vectorR->comp3 / is[2]) * is[2];
  if (vectorR->comp3 < 0)      { vectorR->comp3 += is[2]; vectorM->comp3 =  is[2]; }
  if (vectorR->comp3 == is[2]) { vectorR->comp3 = 0;      vectorM->comp3 = -is[2]; }

  break;

  default:

  if (job->taskid == 0) fprintf(file.out,"Routine rotate_vector_int_latt called with incorrect Crystal->type\n"); 
  MPI_Finalize();
  exit(1);

  break;

 } // close switch

  return;

}

void rotate_vector_int_latt_C(int *symmoperator, int is[3], int trs, VECTOR_INT *vector_in, VECTOR_INT *vectorR)

{

  int *p_irr;

  p_irr = symmoperator;

  vectorR->comp1 = (*p_irr)       * vector_in->comp1 + (*(p_irr + 1)) * vector_in->comp2 + (*(p_irr + 2)) * vector_in->comp3;
  vectorR->comp2 = (*(p_irr + 3)) * vector_in->comp1 + (*(p_irr + 4)) * vector_in->comp2 + (*(p_irr + 5)) * vector_in->comp3;
  vectorR->comp3 = (*(p_irr + 6)) * vector_in->comp1 + (*(p_irr + 7)) * vector_in->comp2 + (*(p_irr + 8)) * vector_in->comp3;

  if (trs == 1) { vectorR->comp1 *= -1; vectorR->comp2 *= -1; vectorR->comp3 *= -1; }
              
  vectorR->comp1 -= (vectorR->comp1 / is[0]) * is[0];
  vectorR->comp2 -= (vectorR->comp2 / is[1]) * is[1];
  vectorR->comp3 -= (vectorR->comp3 / is[2]) * is[2];
  if (vectorR->comp1 < 0) vectorR->comp1 += is[0];
  if (vectorR->comp2 < 0) vectorR->comp2 += is[1];
  if (vectorR->comp3 < 0) vectorR->comp3 += is[2];
  if (vectorR->comp1 == is[0]) vectorR->comp1 = 0;
  if (vectorR->comp2 == is[1]) vectorR->comp2 = 0;
  if (vectorR->comp3 == is[2]) vectorR->comp3 = 0;

  return;

}

/*
void rotate_vector_int_transpose(int *symmoperator, VECTOR_INT *vector_in, VECTOR_INT *vectorR)

{

  int *p_irr;

  p_irr = symmoperator;
  vectorR->comp1 = (*p_irr) * vector_in->comp1 + (*(p_irr + 3)) * vector_in->comp2 +\
 (*(p_irr + 6))
      * vector_in->comp3;
  vectorR->comp2 = (*(p_irr + 1)) * vector_in->comp1 + (*(p_irr + 4)) * vector_in->comp2 +\
 (*(p_irr + 7))
      * vector_in->comp3;
  vectorR->comp3 = (*(p_irr + 2)) * vector_in->comp1 + (*(p_irr + 5)) * vector_in->comp2 +\
 (*(p_irr + 8))
      * vector_in->comp3;

  return;
}

*/

int check_vec(VECTOR_DOUBLE *Rvec_tmp, VECTOR_DOUBLE *Rvec_ai)

{

  double mag = 0.00001;

  if (fabs(Rvec_tmp->comp1 - Rvec_ai->comp1) < mag && fabs(Rvec_tmp->comp2 - Rvec_ai->comp2) < mag && fabs(
      Rvec_tmp->comp3 - Rvec_ai->comp3) < mag)

    return 1;

  else

    return 0;

}

/*
int check_kvec(VECTOR_DOUBLE *Rvec_tmp, VECTOR_DOUBLE *Rvec_ai)

{
  double mag = 0.001;

  if (fabs((*Rvec_tmp).comp1 - (*Rvec_ai).comp1) < mag && fabs((*Rvec_tmp).comp2 - (*Rvec_ai).comp2) < mag && fabs(
      (*Rvec_tmp).comp3 - (*Rvec_ai).comp3) < mag) 

    return 1;

  else

    return 0;
}

int check_product(int *p_irr, int *q_irr, int *r_irr)

{

  double A[3][3], eps = 0.0000001;
  int n = 0;

  A[0][0] = *(p_irr) * *(q_irr) + *(p_irr + 1) * *(q_irr + 3) + *(p_irr + 2) * *(q_irr + 6);
  A[0][1] = *(p_irr + 3) * *(q_irr) + *(p_irr + 4) * *(q_irr + 3) + *(p_irr + 5) * *(q_irr + 6);
  A[0][2] = *(p_irr + 6) * *(q_irr) + *(p_irr + 7) * *(q_irr + 3) + *(p_irr + 8) * *(q_irr + 6);
  A[1][0] = *(p_irr) * *(q_irr + 1) + *(p_irr + 1) * *(q_irr + 4) + *(p_irr + 2) * *(q_irr + 7);
  A[1][1] = *(p_irr + 3) * *(q_irr + 1) + *(p_irr + 4) * *(q_irr + 4) + *(p_irr + 5) * *(q_irr + 7);
  A[1][2] = *(p_irr + 6) * *(q_irr + 1) + *(p_irr + 7) * *(q_irr + 4) + *(p_irr + 8) * *(q_irr + 7);
  A[2][0] = *(p_irr) * *(q_irr + 2) + *(p_irr + 1) * *(q_irr + 5) + *(p_irr + 2) * *(q_irr + 8);
  A[2][1] = *(p_irr + 3) * *(q_irr + 2) + *(p_irr + 4) * *(q_irr + 5) + *(p_irr + 5) * *(q_irr + 8);
  A[2][2] = *(p_irr + 6) * *(q_irr + 2) + *(p_irr + 7) * *(q_irr + 5) + *(p_irr + 8) * *(q_irr + 8);

  if (fabs(A[0][0] - *(r_irr)) < eps && fabs(A[0][1] - *(r_irr + 1)) < eps && fabs(A[0][2] - *(r_irr + 2)) < eps
      && fabs(A[1][0] - *(r_irr + 3)) < eps && fabs(A[1][1] - *(r_irr + 4)) < eps && fabs(A[1][2] - *(r_irr + 5)) < eps
      && fabs(A[2][0] - *(r_irr + 6)) < eps && fabs(A[2][1] - *(r_irr + 7)) < eps && fabs(A[2][2] - *(r_irr + 8)) < eps) 

      n = 1;

  return n;

}
*/

