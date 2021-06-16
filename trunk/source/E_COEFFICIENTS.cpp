#include "mycomplex.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "RECURSION.h"
#include "E_COEFFICIENTS.h"

using namespace std;

void E_coefficients(int ip, int jp, int gi, int gj, int index_i, int index_j, int gausposi, int gausposj, double *C1x, double *C1y, double *C1z, double *C1_max, double *sab, double *pab_inv, double *R_AB_1esqrd, VECTOR_DOUBLE *R_AB, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  // **************************************************************************************************
  // * Compute various coefficients for shells index_i and index_j                                    *
  // **************************************************************************************************

int i4, j4;
int t, m, n;
int imax, jmax, count;
int off1, off2, off3;
double ab, pab, p32, expnta, expntb;
double SAB, KAB;
double PAx, PAy, PAz, PBx, PBy, PBz;
double C1x_max, C1y_max, C1z_max;

  C1x_max = k_zero;
  C1y_max = k_zero;
  C1z_max = k_zero;
 *C1_max  = k_zero;
  imax    = shells->imax_sh[index_i];
  jmax    = shells->imax_sh[index_j];
  off3    = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2    = (imax + jmax + 1) * off3;
  off1    = (jmax + 1) * off2;
  count = 0;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
    for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
      pab = gaussians->expo_sh[gausposi + i4] + gaussians->expo_sh[gausposj + j4];
      ab  = gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4];
      p32 = pab * sqrt(pab);
      KAB = ab * *R_AB_1esqrd / pab ;
      SAB = pi32 * exp(-KAB) / p32 ;
      pab_inv[count] = k_one / pab;
      sab[count]     = gaussians->c_sh[gausposi + i4] * gaussians->c_sh[gausposj + j4] * SAB;
      expnta = gaussians->expo_sh[gausposi + i4];
      expntb = gaussians->expo_sh[gausposj + j4];
      R_AB[count].comp1 = (expnta * (atoms->cell_vector[ip].comp1 + R->vec_ai[gi].comp1) + \
      expntb * (atoms->cell_vector[jp].comp1 + R->vec_ai[gj].comp1)) / pab;
      R_AB[count].comp2 = (expnta * (atoms->cell_vector[ip].comp2 + R->vec_ai[gi].comp2) + \
      expntb * (atoms->cell_vector[jp].comp2 + R->vec_ai[gj].comp2)) / pab;
      R_AB[count].comp3 = (expnta * (atoms->cell_vector[ip].comp3 + R->vec_ai[gi].comp3) + \
      expntb * (atoms->cell_vector[jp].comp3 + R->vec_ai[gj].comp3)) / pab;
      PAx = R_AB[count].comp1 - atoms->cell_vector[ip].comp1 - R->vec_ai[gi].comp1;
      PAy = R_AB[count].comp2 - atoms->cell_vector[ip].comp2 - R->vec_ai[gi].comp2;
      PAz = R_AB[count].comp3 - atoms->cell_vector[ip].comp3 - R->vec_ai[gi].comp3;
      PBx = R_AB[count].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
      PBy = R_AB[count].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
      PBz = R_AB[count].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
      for (m = 0; m <= imax; m++) {
        for (n = 0; n <= jmax; n++) {
          for (t = 0; t <= imax + jmax; t++) {
            C1x[m * off1 + n * off2 + t * off3 + count] = e(m, n, t, pab, PAx, PBx);
            C1y[m * off1 + n * off2 + t * off3 + count] = e(m, n, t, pab, PAy, PBy);
            C1z[m * off1 + n * off2 + t * off3 + count] = e(m, n, t, pab, PAz, PBz);
            C1x_max = (C1x_max > fabs(C1x[m * off1 + n * off2 + t * off3 + count])) ? \
            C1x_max : fabs(C1x[m * off1 + n * off2 + t * off3 + count]);
            C1y_max = (C1y_max > fabs(C1y[m * off1 + n * off2 + t * off3 + count])) ? \
            C1y_max : fabs(C1y[m * off1 + n * off2 + t * off3 + count]);
            C1z_max = (C1z_max > fabs(C1z[m * off1 + n * off2 + t * off3 + count])) ? \
            C1z_max : fabs(C1z[m * off1 + n * off2 + t * off3 + count]);
          }
         }
        } // end loop to set up C1 factors
       *C1_max = (*C1_max > C1x_max * C1y_max * C1z_max) ? *C1_max : C1x_max * C1y_max * C1z_max;
        count++;
      }
     }

}

void E_coefficients1(int ip, int jp, int gi, int gj, int index_i, int index_j, int gausposi, int gausposj, int im, int jm, double *C1x, double *C1y, double *C1z, double *C1_max, double *sab, double *pab_inv, double *R_AB_1esqrd, VECTOR_DOUBLE *R_AB, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

int i4, j4;
int t, m, n;
int count;
int off1, off2, off3;
int imax, jmax;
double ab, pab, p32, expnta, expntb;
double SAB, KAB;
double PAx, PAy, PAz, PBx, PBy, PBz;
double C1x_max, C1y_max, C1z_max;

  // **************************************************************************************************
  // * Compute various coefficients for shells index_i and index_j for 2C integrals                   *
  // **************************************************************************************************

  C1x_max = k_zero;
  C1y_max = k_zero;
  C1z_max = k_zero;
 *C1_max  = k_zero;
  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3    = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2    = (imax + im + jmax + jm + 1) * off3;
  off1    = (jmax + jm + 1) * off2;
  count = 0;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
    for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
      pab = gaussians->expo_sh[gausposi + i4] + gaussians->expo_sh[gausposj + j4];
      ab  = gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4];
      p32 = pab * sqrt(pab);
      KAB = ab * *R_AB_1esqrd / pab ;
      SAB = pi32 * exp(-KAB) / p32 ;
      pab_inv[count] = k_one / pab;
      sab[count]     = gaussians->c_sh[gausposi + i4] * gaussians->c_sh[gausposj + j4] * SAB;
      expnta = gaussians->expo_sh[gausposi + i4];
      expntb = gaussians->expo_sh[gausposj + j4];
      R_AB[count].comp1 = (expnta * (atoms->cell_vector[ip].comp1 + R->vec_ai[gi].comp1) + \
      expntb * (atoms->cell_vector[jp].comp1 + R->vec_ai[gj].comp1)) / pab;
      R_AB[count].comp2 = (expnta * (atoms->cell_vector[ip].comp2 + R->vec_ai[gi].comp2) + \
      expntb * (atoms->cell_vector[jp].comp2 + R->vec_ai[gj].comp2)) / pab;
      R_AB[count].comp3 = (expnta * (atoms->cell_vector[ip].comp3 + R->vec_ai[gi].comp3) + \
      expntb * (atoms->cell_vector[jp].comp3 + R->vec_ai[gj].comp3)) / pab;
      PAx = R_AB[count].comp1 - atoms->cell_vector[ip].comp1 - R->vec_ai[gi].comp1;
      PAy = R_AB[count].comp2 - atoms->cell_vector[ip].comp2 - R->vec_ai[gi].comp2;
      PAz = R_AB[count].comp3 - atoms->cell_vector[ip].comp3 - R->vec_ai[gi].comp3;
      PBx = R_AB[count].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
      PBy = R_AB[count].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
      PBz = R_AB[count].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
      for (m = 0; m <= imax + im; m++) {
        for (n = 0; n <= jmax + jm; n++) {
          for (t = 0; t <= imax + im + jmax + jm; t++) {
            C1x[m * off1 + n * off2 + t * off3 + count] = e(m, n, t, pab, PAx, PBx);
            C1y[m * off1 + n * off2 + t * off3 + count] = e(m, n, t, pab, PAy, PBy);
            C1z[m * off1 + n * off2 + t * off3 + count] = e(m, n, t, pab, PAz, PBz);
            C1x_max = (C1x_max > fabs(C1x[m * off1 + n * off2 + t * off3 + count])) ? \
            C1x_max : fabs(C1x[m * off1 + n * off2 + t * off3 + count]);
            C1y_max = (C1y_max > fabs(C1y[m * off1 + n * off2 + t * off3 + count])) ? \
            C1y_max : fabs(C1y[m * off1 + n * off2 + t * off3 + count]);
            C1z_max = (C1z_max > fabs(C1z[m * off1 + n * off2 + t * off3 + count])) ? \
            C1z_max : fabs(C1z[m * off1 + n * off2 + t * off3 + count]);
          }
         }
        } // end loop to set up C1 factors
       *C1_max = (*C1_max > C1x_max * C1y_max * C1z_max) ? *C1_max : C1x_max * C1y_max * C1z_max;
        count++;
      }
     }

}

void E_coefficients_1c(int index_i, int gausposi, double *C1x, double *C1y, double *C1z, double *C1_max, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

int i4, j4;
int m, t;
int count;
int off1, off2;
int imax, jmax;
double pa, C1x_max, C1y_max, C1z_max;

  C1x_max = k_zero;
  C1y_max = k_zero;
  C1z_max = k_zero;
 *C1_max  = k_zero;
  imax = shells->imax_sh[index_i];
  off2    = shells->ng_sh[index_i];
  off1    = (imax + 1) * off2;

  count = 0;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
    pa = gaussians->expo_sh[gausposi + i4];
    for (m = 0; m <= imax; m++) {
      for (t = 0; t <= imax; t++) {
        C1x[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
        C1y[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
        C1z[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
       *C1_max = (*C1_max > C1x_max * C1y_max * C1z_max) ? *C1_max : C1x_max * C1y_max * C1z_max;
       }
      }
     count++;
     }

}

void E_coefficients_2c(int ip, int jp, int gi, int gj, int index_i, int index_j, int gausposi, int gausposj, double *C1x, double *C1y, double *C1z, double *C1_max, double *C2x, double *C2y, double *C2z, double *C2_max, double *sab_fac, double *ab_inv, double *R_AB_1esqrd, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

int i4, j4;
int m, t;
int count;
int off1, off2, off3, off4;
int imax, jmax;
double ab, pab, p32, pa, pb;
double SAB, KAB;
double C1x_max, C1y_max, C1z_max;
double C2x_max, C2y_max, C2z_max;

  C1x_max = k_zero;
  C1y_max = k_zero;
  C1z_max = k_zero;
 *C1_max  = k_zero;
  C2x_max = k_zero;
  C2y_max = k_zero;
  C2z_max = k_zero;
 *C2_max  = k_zero;
  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off4    = shells->ng_sh[index_j];
  off3    = (jmax + 1) * off4;
  off2    = shells->ng_sh[index_i];
  off1    = (imax + 1) * off2;

  count = 0;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
    pa = gaussians->expo_sh[gausposi + i4];
    for (m = 0; m <= imax; m++) {
      for (t = 0; t <= imax; t++) {
        C1x[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
        C1y[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
        C1z[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
       *C1_max = (*C1_max > C1x_max * C1y_max * C1z_max) ? *C1_max : C1x_max * C1y_max * C1z_max;
       }
      }
     count++;
     }

  count = 0;
  for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
    pb = gaussians->expo_sh[gausposj + j4];
    for (m = 0; m <= jmax; m++) {
      for (t = 0; t <= jmax; t++) {
        C2x[m * off3 + t * off4 + count] = e(m, 0, t, pb, k_zero, k_zero);
        C2y[m * off3 + t * off4 + count] = e(m, 0, t, pb, k_zero, k_zero);
        C2z[m * off3 + t * off4 + count] = e(m, 0, t, pb, k_zero, k_zero);
       *C2_max = (*C2_max > C2x_max * C2y_max * C2z_max) ? *C2_max : C2x_max * C2y_max * C2z_max;
       }
      }
     count++;
     }

}
