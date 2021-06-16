#include <cstring>
#include <mpi.h>
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "LIMITS.h"
#include "PARALLEL.h" 
#include "RECURSION.h"
#include "DFT.h"
#include "ROTATIONS_MOLECULE.h"
#include "CARTESIAN_TO_SH.h"
#include "MCMURCHIE_DAVIDSON.h"
#include "E_COEFFICIENTS.h"
#include "INCOMPLETE_GAMMA.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "INTEGRALS_2C_MOLECULE.h"

using namespace std;

void integrals_crystal_coulomb_ijkl(INTEGRAL_LIST *integral_list, QUAD_TRAN *quad, int *start_index, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Compute four centre coulomb integrals (i|O)(j|gj)(k|O)(l|gl) over one quad             *
  // ******************************************************************************************
  
int index_i, index_j, index_k, index_l, index_G, index_R, index_S, i4, j4, k4, l4;
int bfposi, bfposj, bfposk, bfposl;
int bfposi1, bfposj1, bfposk1, bfposl1;
int gausposi, gausposj, gausposk, gausposl;
int int_count, count;
int i, j, k, l, m, n, p, pm, q, r;
int mm, mm0, mm1, mm2, mm3, mm4;
int ip, jp, kp, lp, gi, gj, gk, gl;
int i1, ijk;
int imax, jmax, kmax, lmax;
int shelposi, shelposj, shelposk, shelposl;
int sheli, shelj, shelk, shell;
int sheli1, shelj1, shelk1, shell1;
int shell_index;
int offset1, offset2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int dimFsh;
double pi_vol = pi / crystal->primitive_cell_volume;
double gamma_0, root_gamma_0, gamma_1, gamma_1_inv, gamma_inv, small;
double p_inv_ex, fac1, fac2;
double dot_product;
double shell_sum;
double R_AB_1esqrd, R_CD_1esqrd, Rsqrd;
double C1_max, C2_max, C_max;
double fgtuv_max, fgtuv_temp;
double time1, time2;
double *F_cart, *p_F_cart, *F_sh, *p_F_sh, *fgtuv, *p_fgtuv;
double en[R->last_ewald_vector][55], em[R->last_ewald_vector][55], in[55], hn[55];
VECTOR_DOUBLE R_AB_1e, R_CD_1e, r_12, s_12, t_12, Rvec_tmp;

  time1 = MPI_Wtime();

  gamma_1_inv = G->gamma_0_inv;
  gamma_1 = k_one / G->gamma_0_inv;
  gamma_0 = k_one / G->gamma_0_inv;
  root_gamma_0 = sqrt(gamma_0);

  ip = quad->cell1[0];
  gi = 0;
  jp = quad->cell2[0];
  gj = quad->latt2[0];

  kp = quad->cell3[0];
  gk = quad->latt3[0];
  lp = quad->cell4[0];
  gl = quad->latt4[0];

  R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
  R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
  R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
  R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

  R_CD_1e.comp1 = atoms->cell_vector[kp].comp1 + R->vec_ai[gk].comp1 - atoms->cell_vector[lp].comp1 - R->vec_ai[gl].comp1;
  R_CD_1e.comp2 = atoms->cell_vector[kp].comp2 + R->vec_ai[gk].comp2 - atoms->cell_vector[lp].comp2 - R->vec_ai[gl].comp2;
  R_CD_1e.comp3 = atoms->cell_vector[kp].comp3 + R->vec_ai[gk].comp3 - atoms->cell_vector[lp].comp3 - R->vec_ai[gl].comp3;
  R_CD_1esqrd = double_vec_dot(&R_CD_1e, &R_CD_1e);

  int_count = 0;
  shell_index = 0;

  bfposi   = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip];
  gausposi = atoms->gausposn_sh[ip];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
    sheli  = shells->type1_sh[index_i];
    sheli1 = shells->type_sh[index_i];
    imax   = shells->imax_sh[index_i];
    bfposj   = 0;
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp];
    gausposj = atoms->gausposn_sh[jp];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
      shelj  = shells->type1_sh[index_j];
      shelj1 = shells->type_sh[index_j];
      jmax   = shells->imax_sh[index_j];
      double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1x_max[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1y_max[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1z_max[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
        shells,gaussians,job,file);
        bfposk   = 0;
        bfposk1  = 0;
        shelposk = atoms->shelposn_sh[kp];
        gausposk = atoms->gausposn_sh[kp];
        for (index_k = shelposk; index_k < shelposk + atoms->nshel_sh[kp]; index_k++) {
          shelk  = shells->type1_sh[index_k];
          shelk1 = shells->type_sh[index_k];
          kmax   = shells->imax_sh[index_k];
          bfposl   = 0;
          bfposl1  = 0;
          shelposl = atoms->shelposn_sh[lp];
          gausposl = atoms->gausposn_sh[lp];
          for (index_l = shelposl; index_l < shelposl + atoms->nshel_sh[lp]; index_l++) {
            shell  = shells->type1_sh[index_l];
            shell1 = shells->type_sh[index_l];
            lmax   = shells->imax_sh[index_l];
            if (start_index[shell_index] == 0) {
              bfposl   += shells->type1_sh[index_l];
              bfposl1  += shells->type_sh[index_l];
              gausposl += shells->ng_sh[index_l];
              shell_index++;
              continue;
             }
              shell_index++;
              double scd[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
              double C2x[shells->ng_sh[index_k] * shells->ng_sh[index_l] * (kmax+lmax+1) * (kmax+1) * (lmax+1)];
              double C2y[shells->ng_sh[index_k] * shells->ng_sh[index_l] * (kmax+lmax+1) * (kmax+1) * (lmax+1)];
              double C2z[shells->ng_sh[index_k] * shells->ng_sh[index_l] * (kmax+lmax+1) * (kmax+1) * (lmax+1)];
              double C2x_max[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
              double C2y_max[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
              double C2z_max[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
              double pcd_inv[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
              VECTOR_DOUBLE R_CD[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
              E_coefficients(kp,lp,gk,gl,index_k,index_l,gausposk,gausposl,C2x,C2y,C2z,&C2_max,scd,pcd_inv,&R_CD_1esqrd,R_CD,R,atoms,\
              shells,gaussians,job,file);
              C_max = C1_max * C2_max;

              mm  = imax + jmax + kmax + lmax;
              mm4 = shells->ng_sh[index_k] * shells->ng_sh[index_l];
              mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
              mm2 = (mm + 1) * mm3;
              mm1 = (mm + 1) * mm2;
              mm0 = (mm + 1) * mm1;
             
              fgtuv = (double *) calloc(mm0, sizeof(double));
              if (fgtuv == NULL) {
              fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
              exit(1);
             }
             
              F_cart = (double *) malloc(sheli * shelj * shelk * shell * sizeof(double));
              if (F_cart == NULL) {
              fprintf(stderr, "ERROR: not enough memory for double F_cart! \n");
              exit(1);
             }
         
             F_sh = (double *) malloc(sheli1 * shelj1 * shelk1 * shell1 * sizeof(double));
             if (F_sh == NULL) {
             fprintf(stderr, "ERROR: not enough memory for double F_sh! \n");
             exit(1);
            }
         
             double sum5[16 + 1], sum2[16 + 1], x;
             double fac1, fac2, Rzsqrd, expfac, x1, x2, sign;
             double D_erf[16 + 1], D_exp[16 + 1], derivative_1[16 + 1], derivative_2[16 + 1];

  // ******************************************************************************************
  // * Reciprocal space part for C or S periodicity                                           *
  // ******************************************************************************************
  
             switch (crystal->type[0]) {

               case 'C':

               fgtuv_max = k_zero;
               for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
                 C_max = C1x_max[i4] * C1y_max[i4] * C1z_max[i4];
                 for (j4 = 0; j4 < shells->ng_sh[index_k] * shells->ng_sh[index_l]; j4++) {
                   C_max *= C2x_max[j4] * C2y_max[j4] * C2z_max[j4];
                   s_12.comp1 = R_AB[i4].comp1 - R_CD[j4].comp1;
                   s_12.comp2 = R_AB[i4].comp2 - R_CD[j4].comp2;
                   s_12.comp3 = R_AB[i4].comp3 - R_CD[j4].comp3;
                   fac1 = sab[i4] * scd[j4];
                   p_inv_ex = pab_inv[i4] + pcd_inv[j4];
                   p_fgtuv = fgtuv + i4 * mm4 + j4;
                  *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pcd_inv[j4]);
                   count = 1;
                   for (index_S = 1; index_S < G->number_of_shells; index_S++) {
                     fac2 = fac1 * G->EXPFAC[index_S];
                     for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                       dot_product = G->vec_b2[count].comp1 * s_12.comp1 + G->vec_b2[count].comp2 * s_12.comp2 + \
                       G->vec_b2[count].comp3 * s_12.comp3;
                       for (i = 0; i <= mm; i++) {
                         for (j = 0; j <= mm; j++) {
                           for (k = 0; k <= mm; k++) {
                             ijk = i + j + k;
                             if (ijk > mm) break;
                             fgtuv_temp = fac2 * G->x[count * 9 + i] * G->y[count * 9 + j] * G->z[count * 9 + k] * \
                             cosfactor(ijk, dot_product);
                             p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                            *p_fgtuv += fgtuv_temp;
                             fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                             fgtuv_max = k_one;
                            }
                           }
                          }
                         count++;
                        }
                       }
                      }
                     }
 
               break;
 
               case 'S':
 
               fgtuv_max = k_zero;
               for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
                 C_max = C1x_max[i4] * C1y_max[i4] * C1z_max[i4];
                 for (j4 = 0; j4 < shells->ng_sh[index_k] * shells->ng_sh[index_l]; j4++) {
                   C_max *= C2x_max[j4] * C2y_max[j4] * C2z_max[j4];
                   s_12.comp1 = R_AB[i4].comp1 - R_CD[j4].comp1;
                   s_12.comp2 = R_AB[i4].comp2 - R_CD[j4].comp2;
                   s_12.comp3 = R_AB[i4].comp3 - R_CD[j4].comp3;
                   Rsqrd = double_vec_dot(&s_12, &s_12);
                   Rzsqrd = s_12.comp3 * s_12.comp3;
                   fac1 = sab[i4] * scd[j4];
                   p_inv_ex = pab_inv[i4] + pcd_inv[j4];
                   p_fgtuv = fgtuv + i4 * mm4 + j4;
                   x = root_gamma_0 * s_12.comp3;
                   erf_exp_derivative(D_erf, D_exp, mm, root_gamma_0, x, G);
                   sum5[0] = D_exp[0] / rtpi / root_gamma_0 + s_12.comp3 * D_erf[0];
                   for (n = 1; n <= mm; n++) {
                     sum5[n] = D_exp[n] / rtpi / root_gamma_0 + s_12.comp3 * D_erf[n] + (double) n * D_erf[n - 1];
                    }
                   for (k = 0; k <= mm; k++) {
                     p_fgtuv = fgtuv + k * mm3  + i4 * mm4 + j4;
                    *p_fgtuv -= fac1 * two * pi / crystal->primitive_cell_volume * sum5[k]; 
                   }
                   count = 1;
                   for (index_S = 1; index_S < G->number_of_shells; index_S++) {
                     x1 = G->shell_mag[9 * index_S + 1] / two / root_gamma_0 + root_gamma_0 * s_12.comp3;
                     erfc_derivative(derivative_1, mm, root_gamma_0, x1, G);
                     x2 = G->shell_mag[9 * index_S + 1] / two / root_gamma_0 - root_gamma_0 * s_12.comp3;
                     erfc_derivative(derivative_2, mm, -root_gamma_0, x2, G);
                     expfac = exp(G->shell_mag[9 * index_S + 1] * s_12.comp3);
                     for (n = 0; n <= mm; n++) {
                       sign = k_one;
                       sum2[n] = k_zero;
                       for (k = 0; k <= n; k++) {
                         sum2[n] += G->A->a[n][k] * G->shell_mag[9 * index_S + k] * (expfac * derivative_1[n - k] + sign / expfac * derivative_2[n - k]) / \
                         G->shell_mag[9 * index_S + 1];
                         sign *= -k_one;
                        }
                       }
                     for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                       dot_product = G->vec_b2[count].comp1 * s_12.comp1 + G->vec_b2[count].comp2 * s_12.comp2;
                       for (i = 0; i <= mm; i++) {
                         for (j = 0; j <= mm; j++) {
                           for (k = 0; k <= mm; k++) {
                             ijk = i + j + k;
                             if (ijk > mm) break;
                             fac2 = fac1 * two * pi / crystal->primitive_cell_volume * sum2[k];
                             fgtuv_temp = fac2 * G->x[count * 9 + i] * G->y[count * 9 + j] * cosfactor(i + j, dot_product);
                             p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                            *p_fgtuv += fgtuv_temp;
                             fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                             fgtuv_max = k_one;
                            }
                           }
                          }
                         count++;
                        }
                       }
                      }
                     }
               break;
 
               case 'P':
 
               if (job->taskid == 0)
               fprintf(file.out,"Electron-electron repulsion not coded for 1-D systems\n");
               MPI_Finalize();
               exit(1);
          
               break;
          
               case 'M':
          
               if (job->taskid == 0)
               fprintf(file.out,"integrals_crystal_coulomb_ijkl routine is for periodic systems only\n");
               MPI_Finalize();
               exit(1);
          
               break;
          
             } // close switch

  // ******************************************************************************************
  // * Real space part for any periodicity                                                    *
  // ******************************************************************************************
  
  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
    C_max = C1x_max[i4] * C1y_max[i4] * C1z_max[i4];
    for (j4 = 0; j4 < shells->ng_sh[index_k] * shells->ng_sh[index_l]; j4++) {
      C_max *= C2x_max[j4] * C2y_max[j4] * C2z_max[j4];
      p_inv_ex = pab_inv[i4] + pcd_inv[j4];
      s_12.comp1 = R_AB[i4].comp1 - R_CD[j4].comp1;
      s_12.comp2 = R_AB[i4].comp2 - R_CD[j4].comp2;
      s_12.comp3 = R_AB[i4].comp3 - R_CD[j4].comp3;
      map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      fac1 = sab[i4] * scd[j4];
      count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
          r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
          r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
          r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
          Rsqrd = r_12.comp1 * r_12.comp1 + r_12.comp2 * r_12.comp2 + r_12.comp3 * r_12.comp3;
          f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
          f000m(&en[count][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
          for (i1 = 0; i1 <= mm; i1++) {
            en[count][i1] -= em[count][i1];
           }
            shell_sum += fabs(en[count][0]);
            count++;
           } // end loop on index_R
          } // end loop on index_S
           for (index_R = 0; index_R < count; index_R++) {
             r_12.comp1 = t_12.comp1 + R->vec_ai[index_R].comp1;
             r_12.comp2 = t_12.comp2 + R->vec_ai[index_R].comp2;
             r_12.comp3 = t_12.comp3 + R->vec_ai[index_R].comp3;
             for (i = 0; i <= mm; i++) {
               for (j = 0; j <= mm; j++) {
                 for (k = 0; k <= mm; k++) {
                   ijk = i + j + k;
                   if (ijk > mm) break;
                   fgtuv_temp = fac1 * ftuvn(i,j,k,0,&en[index_R][0],r_12);
                   p_fgtuv  = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                  *p_fgtuv += fgtuv_temp;
                  }
                 }
                }
               } 
              } 
             } 

  // ******************************************************************************************
  // * McMurchie-Davidson loop and transform to SH                                            *
  // ******************************************************************************************
  
      mcmurchie_davidson_ijkl(F_cart, index_i, index_j, index_k, index_l, C1x, C1y, C1z, C2x, C2y, C2z, fgtuv, shells, job, file);
      free(fgtuv);
      dimFsh = sheli1 * shelj1 * shelk1 * shell1;
      ResetDoubleArray(F_sh,&dimFsh);
      cartesian_to_sh_ijkl(F_cart, F_sh, index_i, index_j, index_k, index_l, shells, job, file);
      //four_centre_cartesian_to_sh_ijkl(F_cart, F_sh, index_i, index_j, index_k, index_l, shells, job, file);
  
      p_F_sh = F_sh;
      for (i = 0; i < sheli1; i++) {
        for (j = 0; j < shelj1; j++) {
          for (k = 0; k < shelk1; k++) {
            for (l = 0; l < shell1; l++) {
              if (fabs(*p_F_sh) > integral_rejection_threshold) {
              integral_list->value[int_count] = *p_F_sh ;
              integral_list->i[int_count] = bfposi1 + i;
              integral_list->j[int_count] = bfposj1 + j;
              integral_list->k[int_count] = bfposk1 + k;
              integral_list->l[int_count] = bfposl1 + l;
              int_count++;
              }
              p_F_sh++;
             }
            }
           }
          }
         free(F_cart);
         free(F_sh);
         bfposl   += shell;
         bfposl1  += shell1;
         gausposl += shells->ng_sh[index_l];
        } // close loop over index_l
       bfposk   += shelk;
       bfposk1  += shelk1;
       gausposk += shells->ng_sh[index_k];
      } // close loop over index_k
     bfposj   += shelj;
     bfposj1  += shelj1;
     gausposj += shells->ng_sh[index_j];
    } // close loop over index_j
   bfposi   += sheli;
   bfposi1  += sheli1;
   gausposi += shells->ng_sh[index_i];
  } // close loop over index_i

  integral_list->num = int_count;

  time2 = MPI_Wtime();

  if (job->taskid == 0 && job->verbosity > 1) fprintf(file.out, "Time for %d integrals no sym %10.4e\n", int_count,time2-time1);

}

void integrals_crystal_exchange_ijkl(INTEGRAL_LIST *integral_list, QUAD_TRAN *quad, int *start_index, int *count_index, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Compute four centre exchange integrals (i|O)(j|gj)(k|gk)(l|gl) over one quad           *
  // ******************************************************************************************
  
int index_i, index_j, index_k, index_l, i4, j4, k4, l4;
int bfposi, bfposj, bfposk, bfposl;
int bfposi1, bfposj1, bfposk1, bfposl1;
int gausposi, gausposj, gausposk, gausposl;
int int_count;
int i, j, k, l, m, n, p, pm, q;
int mm, mm0, mm1, mm2, mm3, mm4;
int ip, jp, kp, lp, gi, gj, gk, gl;
int imax, jmax, kmax, lmax;
int shelposi, shelposj, shelposk, shelposl;
int sheli, shelj, shelk, shell;
int sheli1, shelj1, shelk1, shell1;
int shell_index;
int offset1, offset2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int dimFsh;
double p_inv_ex, fac1;
double R_AB_1esqrd, R_CD_1esqrd, Rsqrd;
double C1_max, C2_max, C_max;
double fgtuv_max, fgtuv_temp;
double time1, time2;
double em[1][55];
double *F_cart, *p_F_cart, *F_sh, *p_F_sh, *fgtuv, *p_fgtuv;
VECTOR_DOUBLE R_AB_1e, R_CD_1e, s_12;

  time1 = MPI_Wtime();

  ip = quad->cell1[0];
  gi = 0;
  jp = quad->cell2[0];
  gj = quad->latt2[0];

  kp = quad->cell3[0];
  gk = quad->latt3[0];
  lp = quad->cell4[0];
  gl = quad->latt4[0];

  R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
  R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
  R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
  R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

  R_CD_1e.comp1 = atoms->cell_vector[kp].comp1 + R->vec_ai[gk].comp1 - atoms->cell_vector[lp].comp1 - R->vec_ai[gl].comp1;
  R_CD_1e.comp2 = atoms->cell_vector[kp].comp2 + R->vec_ai[gk].comp2 - atoms->cell_vector[lp].comp2 - R->vec_ai[gl].comp2;
  R_CD_1e.comp3 = atoms->cell_vector[kp].comp3 + R->vec_ai[gk].comp3 - atoms->cell_vector[lp].comp3 - R->vec_ai[gl].comp3;
  R_CD_1esqrd = double_vec_dot(&R_CD_1e, &R_CD_1e);

  int_count = 0;
  shell_index = 0;

  bfposi   = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip];
  gausposi = atoms->gausposn_sh[ip];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
    sheli  = shells->type1_sh[index_i];
    sheli1 = shells->type_sh[index_i];
    imax   = shells->imax_sh[index_i];
    bfposj   = 0;
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp];
    gausposj = atoms->gausposn_sh[jp];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
      shelj  = shells->type1_sh[index_j];
      shelj1 = shells->type_sh[index_j];
      jmax   = shells->imax_sh[index_j];
      double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
      shells,gaussians,job,file);
      bfposk   = 0;
      bfposk1  = 0;
      shelposk = atoms->shelposn_sh[kp];
      gausposk = atoms->gausposn_sh[kp];
      for (index_k = shelposk; index_k < shelposk + atoms->nshel_sh[kp]; index_k++) {
        shelk  = shells->type1_sh[index_k];
        shelk1 = shells->type_sh[index_k];
        kmax   = shells->imax_sh[index_k];
        bfposl   = 0;
        bfposl1  = 0;
        shelposl = atoms->shelposn_sh[lp];
        gausposl = atoms->gausposn_sh[lp];
        for (index_l = shelposl; index_l < shelposl + atoms->nshel_sh[lp]; index_l++) {
          shell  = shells->type1_sh[index_l];
          shell1 = shells->type_sh[index_l];
          lmax   = shells->imax_sh[index_l];
          if (start_index[shell_index] == 0) {
            bfposl   += shells->type1_sh[index_l];
            bfposl1  += shells->type_sh[index_l];
            gausposl += shells->ng_sh[index_l];
            shell_index++;
            continue;
           }
            shell_index++;
            double scd[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
            double C2x[shells->ng_sh[index_k] * shells->ng_sh[index_l] * (kmax+lmax+1) * (kmax+1) * (lmax+1)];
            double C2y[shells->ng_sh[index_k] * shells->ng_sh[index_l] * (kmax+lmax+1) * (kmax+1) * (lmax+1)];
            double C2z[shells->ng_sh[index_k] * shells->ng_sh[index_l] * (kmax+lmax+1) * (kmax+1) * (lmax+1)];
            double pcd_inv[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
            VECTOR_DOUBLE R_CD[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
            E_coefficients(kp,lp,gk,gl,index_k,index_l,gausposk,gausposl,C2x,C2y,C2z,&C2_max,scd,pcd_inv,&R_CD_1esqrd,R_CD,R,atoms,\
            shells,gaussians,job,file);
            C_max = C1_max * C2_max;

            mm  = imax + jmax + kmax + lmax;
            mm4 = shells->ng_sh[index_k] * shells->ng_sh[index_l];
            mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
            mm2 = (mm + 1) * mm3;
            mm1 = (mm + 1) * mm2;
            mm0 = (mm + 1) * mm1;
   
            fgtuv = (double *) calloc(mm0, sizeof(double));
            if (fgtuv == NULL) {
            fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
            exit(1);
           }

            F_cart = (double *) malloc(sheli * shelj * shelk * shell * sizeof(double));
            if (F_cart == NULL) {
            fprintf(stderr, "ERROR: not enough memory for double F_cart! \n");
            exit(1);
           }
         
           F_sh = (double *) malloc(sheli1 * shelj1 * shelk1 * shell1 * sizeof(double));
           if (F_sh == NULL) {
           fprintf(stderr, "ERROR: not enough memory for double F_sh! \n");
           exit(1);
          }
         
  // ******************************************************************************************
  // * Real space part for any periodicity                                                    *
  // ******************************************************************************************
  
           fgtuv_max = k_zero;
           for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
             for (j4 = 0; j4 < shells->ng_sh[index_k] * shells->ng_sh[index_l]; j4++) {
               s_12.comp1 = R_AB[i4].comp1 - R_CD[j4].comp1;
               s_12.comp2 = R_AB[i4].comp2 - R_CD[j4].comp2;
               s_12.comp3 = R_AB[i4].comp3 - R_CD[j4].comp3;
               fac1 = sab[i4] * scd[j4];
               p_inv_ex = pab_inv[i4] + pcd_inv[j4];
               Rsqrd = double_vec_dot(&s_12, &s_12);
               f000m(&em[0][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
               for (m = 0; m <= mm; m++) {
                 for (n = 0; n <= mm; n++) {
                   for (pm = 0; pm <= mm; pm++) {
                     if (m + n + pm > mm) break;
                       fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&em[0][0],s_12);
                       p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
                      *p_fgtuv += fgtuv_temp;
                       fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                       fgtuv_max = k_one;
                      }
                     }
                    }
                   }
                  }

  // ******************************************************************************************
  // * McMurchie-Davidson loop and transform to SH                                            *
  // ******************************************************************************************
  
         mcmurchie_davidson_ijkl(F_cart, index_i, index_j, index_k, index_l, C1x, C1y, C1z, C2x, C2y, C2z, fgtuv, shells, job, file);
         free(fgtuv);
         dimFsh = sheli1 * shelj1 * shelk1 * shell1;
         ResetDoubleArray(F_sh,&dimFsh);
         //four_centre_cartesian_to_sh_ijkl(F_cart, F_sh, index_i, index_j, index_k, index_l, shells, job, file);
         cartesian_to_sh_ijkl(F_cart, F_sh, index_i, index_j, index_k, index_l, shells, job, file);

         p_F_sh = F_sh;
         for (i = 0; i < sheli1; i++) {
           for (j = 0; j < shelj1; j++) {
             for (k = 0; k < shelk1; k++) {
               for (l = 0; l < shell1; l++) {
                 if (fabs(*p_F_sh) > integral_rejection_threshold) {
                 integral_list->value[int_count] = *p_F_sh ;
                 integral_list->i[int_count] = bfposi1 + i;
                 integral_list->j[int_count] = bfposj1 + j;
                 integral_list->k[int_count] = bfposk1 + k;
                 integral_list->l[int_count] = bfposl1 + l;
                 int_count++;
                 }
                 p_F_sh++;
                }
               }
              }
             }
         free(F_cart);
         free(F_sh);
         bfposl   += shell;
         bfposl1  += shell1;
         gausposl += shells->ng_sh[index_l];
        } // close loop over index_l
       bfposk   += shelk;
       bfposk1  += shelk1;
       gausposk += shells->ng_sh[index_k];
      } // close loop over index_k
     bfposj   += shelj;
     bfposj1  += shelj1;
     gausposj += shells->ng_sh[index_j];
    } // close loop over index_j
   bfposi   += sheli;
   bfposi1  += sheli1;
   gausposi += shells->ng_sh[index_i];
  } // close loop over index_i

  integral_list->num = int_count;

  time2 = MPI_Wtime();

  if (job->taskid == 0 && job->verbosity > 1) fprintf(file.out, "Time for %d integrals no sym %10.4e\n", int_count,time2-time1);

}

void integrals_crystal_coulomb_screen_ijij(double *F_sh, int ip, int jp, int gj, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Compute four centre coulomb screening integrals (i|O)(j|gj)(i|O)(j|gj) over one quad   *
  // ******************************************************************************************

int index_i, index_j, index_G, index_R, index_S, i4, j4;
int bfposi, bfposj;
int bfposi1, bfposj1;
int gausposi, gausposj;
int count;
int i, j, k, l, m, n, p, pm, q;
int mm, mm0, mm1, mm2, mm3, mm4;
int i1, ijk;
int gi;
int imax, jmax;
int shelposi, shelposj;
int sheli, shelj;
int sheli1, shelj1;
int offset1, offset2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int nd1 = atoms->bfnnumb[ip];
int nd2 = atoms->bfnnumb[jp];
int nd3 = atoms->bfnnumb_sh[ip];
int nd4 = atoms->bfnnumb_sh[jp];
double p_inv_ex, fac1, fac2, dot_product, shell_sum;
double gamma_1;
double R_AB_1esqrd, Rsqrd;
double C1_max, C_max;
double fgtuv_max, fgtuv_temp;
double time1, time2;
double em[R->last_ewald_vector][55];
double en[R->last_ewald_vector][55];
double pi_vol = pi / crystal->primitive_cell_volume;
double gamma_1_inv = G->gamma_0_inv;
double *F_cart, *fgtuv, *p_fgtuv;
VECTOR_DOUBLE R_AB_1e, r_12, s_12, t_12, Rvec_tmp;

  time1 = MPI_Wtime();

  gamma_1 = k_one / G->gamma_0_inv;

  gi = 0;
  
  R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
  R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
  R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
  R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

  bfposi   = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip];
  gausposi = atoms->gausposn_sh[ip];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
    sheli  = shells->type1_sh[index_i];
    sheli1 = shells->type_sh[index_i];
    imax   = shells->imax_sh[index_i];
    bfposj   = 0;
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp];
    gausposj = atoms->gausposn_sh[jp];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
      shelj  = shells->type1_sh[index_j];
      shelj1 = shells->type_sh[index_j];
      jmax   = shells->imax_sh[index_j];
      double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
      shells,gaussians,job,file);

      mm  = imax + jmax + imax + jmax;
      mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
      mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
      mm2 = (mm + 1) * mm3;
      mm1 = (mm + 1) * mm2;
      mm0 = (mm + 1) * mm1;

      fgtuv = (double *) calloc(mm0, sizeof(double));
      if (fgtuv == NULL) {
      fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
      exit(1);
     }

      F_cart = (double *) malloc(sheli * shelj * sheli * shelj * sizeof(double));
      if (F_cart == NULL) {
      fprintf(stderr, "ERROR: not enough memory for double F_cart! \n");
      exit(1);
     }
       
  // ******************************************************************************************
  // * Reciprocal space part for C or S periodicity                                           *
  // ******************************************************************************************
  
      switch (crystal->type[0]) {

        case 'C':

        fgtuv_max = k_zero;
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
            s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
            s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
            s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
            fac1 = sab[i4] * sab[j4];
            p_inv_ex = pab_inv[i4] + pab_inv[j4];
            p_fgtuv = fgtuv + i4 * mm4 + j4;
           *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pab_inv[j4]);
            count = 1;
            for (index_S = 1; index_S < G->number_of_shells; index_S++) {
              fac2 = fac1 * G->EXPFAC[index_S];
              for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                dot_product = G->vec_b2[count].comp1 * s_12.comp1 + G->vec_b2[count].comp2 * s_12.comp2 + \
                G->vec_b2[count].comp3 * s_12.comp3;
                for (i = 0; i <= mm; i++) {
                  for (j = 0; j <= mm; j++) {
                    for (k = 0; k <= mm; k++) {
                      ijk = i + j + k;
                      if (ijk > mm) break;
                      fgtuv_temp = fac2 * G->x[count * 9 + i] * G->y[count * 9 + j] * G->z[count * 9 + k] * \
                      cosfactor(ijk, dot_product);
                      p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                     *p_fgtuv += fgtuv_temp;
                      fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                      fgtuv_max = k_one;
                     }
                    }
                   }
                  count++;
                 }
                }
               }
              }

        break;
     
        case 'S':
        case 'P':
     
        break;
     
        case 'M':
     
        if (job->taskid == 0)
        fprintf(file.out,"integrals_crystal_coulomb_screen_ijij routine is for periodic systems only\n");
        MPI_Finalize();
        exit(1);
     
        break;
     
      } // close switch

  // ******************************************************************************************
  // * Real space part for any periodicity                                                    *
  // ******************************************************************************************
  
  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
    for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
      s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
      s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
      s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
      map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      fac1 = sab[i4] * sab[j4];
      p_inv_ex = pab_inv[i4] + pab_inv[j4];
      count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
          r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
          r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
          r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
          Rsqrd = double_vec_dot(&r_12, &r_12);
          f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
          f000m(&en[count][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
          for (i1 = 0; i1 <= mm; i1++) {
            en[count][i1] -= em[count][i1];
           }
            shell_sum += fabs(en[count][0]);
            count++;
           } // end loop on index_R
          } // end loop on index_S

            for (index_R = 0; index_R < count; index_R++) {
              r_12.comp1 = t_12.comp1 + R->vec_ai[index_R].comp1;
              r_12.comp2 = t_12.comp2 + R->vec_ai[index_R].comp2;
              r_12.comp3 = t_12.comp3 + R->vec_ai[index_R].comp3;
              for (m = 0; m <= mm; m++) {
                for (n = 0; n <= mm; n++) {
                  for (pm = 0; pm <= mm; pm++) {
                    if (m + n + pm > mm) break;
                    fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&en[index_R][0],r_12);
                    p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
                   *p_fgtuv += fgtuv_temp;
                    fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                    fgtuv_max = k_one;
                   }
                  }
                 }
                }
               }
              }

  // ******************************************************************************************
  // * McMurchie-Davidson loop and transform to SH                                            *
  // ******************************************************************************************
  
        mcmurchie_davidson_ijij(F_cart, index_i, index_j, bfposi, bfposj, nd2, C1x, C1y, C1z, fgtuv, shells, job, file);
        free(fgtuv);
        //four_centre_cartesian_to_sh_ijij(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        cartesian_to_sh_ijij(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        free(F_cart);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i

}

void integrals_crystal_exchange_screen_ijij(double *F_sh, int ip, int jp, int gj, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Compute four centre coulomb screening integrals (i|O)(j|gj)(i|O)(j|gj) over one quad   *
  // ******************************************************************************************

  int index_i, index_j, i4, j4;
  int bfposi, bfposj;
  int bfposi1, bfposj1;
  int gausposi, gausposj;
  int count;
  int i, j, k, l, m, n, p, pm, q;
  int mm, mm0, mm1, mm2, mm3, mm4;
  int gi;
  int imax, jmax;
  int shelposi, shelposj;
  int sheli, shelj;
  int sheli1, shelj1;
  int offset1, offset2;
  int dim1 = atoms->number_of_atoms_in_unit_cell;
  int dim2 = dim1 * dim1;
  int nd1 = atoms->bfnnumb[ip];
  int nd2 = atoms->bfnnumb[jp];
  int nd3 = atoms->bfnnumb_sh[ip];
  int nd4 = atoms->bfnnumb_sh[jp];
  double p_inv_ex, fac1;
  double R_AB_1esqrd, Rsqrd;
  double C1_max, C_max;
  double fgtuv_max, fgtuv_temp;
  double time1, time2;
  double em[1][55];
  double *F_cart, *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, s_12;

  time1 = MPI_Wtime();

  gi = 0;

  R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
  R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
  R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
  R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

  bfposi   = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip];
  gausposi = atoms->gausposn_sh[ip];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
    sheli  = shells->type1_sh[index_i];
    sheli1 = shells->type_sh[index_i];
    imax   = shells->imax_sh[index_i];
    bfposj   = 0;
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp];
    gausposj = atoms->gausposn_sh[jp];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
      shelj  = shells->type1_sh[index_j];
      shelj1 = shells->type_sh[index_j];
      jmax   = shells->imax_sh[index_j];
      double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
      double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
      shells,gaussians,job,file);

      mm  = imax + jmax + imax + jmax;
      mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
      mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
      mm2 = (mm + 1) * mm3;
      mm1 = (mm + 1) * mm2;
      mm0 = (mm + 1) * mm1;

      fgtuv = (double *) calloc(mm0, sizeof(double));
      if (fgtuv == NULL) {
      fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
      exit(1);
     }

      F_cart = (double *) malloc(sheli * shelj * sheli * shelj * sizeof(double));
      if (F_cart == NULL) {
      fprintf(stderr, "ERROR: not enough memory for double F_cart! \n");
      exit(1);
     }

  // ******************************************************************************************
  // * Real space part for any periodicity                                                    *
  // ******************************************************************************************
  
      fgtuv_max = k_zero;
      for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
        for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
          s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
          s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
          s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
          fac1 = sab[i4] * sab[j4];
          p_inv_ex = pab_inv[i4] + pab_inv[j4];
          Rsqrd = double_vec_dot(&s_12, &s_12);
          f000m(&em[0][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
          for (m = 0; m <= mm; m++) {
            for (n = 0; n <= mm; n++) {
              for (pm = 0; pm <= mm; pm++) {
                if (m + n + pm > mm) break;
                fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&em[0][0],s_12);
                p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
               *p_fgtuv += fgtuv_temp;
                fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                fgtuv_max = k_one;
               }
              }
             }
            }
           }

  // ******************************************************************************************
  // * McMurchie-Davidson loop and transform to SH                                            *
  // ******************************************************************************************
  
      mcmurchie_davidson_ijij(F_cart, index_i, index_j, bfposi, bfposj, nd2, C1x, C1y, C1z, fgtuv, shells, job, file);
      free(fgtuv);
      cartesian_to_sh_ijij(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
      //four_centre_cartesian_to_sh_ijij(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
      free(F_cart);
      bfposj   += shelj;
      bfposj1  += shelj1;
      gausposj += shells->ng_sh[index_j];
     } // close loop over index_j
    bfposi   += sheli;
    bfposi1  += sheli1;
    gausposi += shells->ng_sh[index_i];
   } // close loop over index_i

}

void integrals_crystal_screen_complex(Complex *F_sh, int ip, int jp, int gj, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, Q_LATTICE *q_G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Compute four centre coulomb screening integrals (i|O)(j|gj)(i|O)(j|gj) over one quad   *
  // ******************************************************************************************

  // compute 4-centre integrals for a quad of atoms (ij|ij) using pair_p as indices
  // starting point is pair_p->posn[*p1]

  int index_i, index_j, index_G, index_R, index_S, i4, j4;
  int bfposi, bfposj;
  int bfposi1, bfposj1;
  int gausposi, gausposj;
  int count;
  int i, j, k, l, m, n, p, pm, q;
  int mm, mm0, mm1, mm2, mm3, mm4;
  int i1, ijk;
  int gi;
  int imax, jmax;
  int shelposi, shelposj;
  int sheli, shelj;
  int sheli1, shelj1;
  int offset1, offset2;
  int dim1 = atoms->number_of_atoms_in_unit_cell;
  int dim2 = dim1 * dim1;
  int nd1 = atoms->bfnnumb[ip];
  int nd2 = atoms->bfnnumb[jp];
  int nd3 = atoms->bfnnumb_sh[ip];
  int nd4 = atoms->bfnnumb_sh[jp];
  double p_inv_ex, fac1, fac2, expfac, dot_product, shell_sum;
  double gamma_1;
  double R_AB_1esqrd, Rsqrd;
  double C1_max, C_max;
  //double fgtuv_max, fgtuv_temp;
  Complex fac, fgtuv_max, fgtuv_temp;
  double time1, time2;
  double em[R->last_ewald_vector][55];
  double en[R->last_ewald_vector][55];
  double pi_vol = pi / crystal->primitive_cell_volume;
  double gamma_1_inv = G->gamma_0_inv;
  Complex *F_cart, *p_F_cart, *p_F_sh, *fgtuv, *p_fgtuv;
  //double *F_cart, *p_F_cart, *p_F_sh, *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, r_12, s_12, t_12, Rvec_tmp;

  int size;
  size = 4 * job->l_max + 1;
  gamma_1 = k_one / G->gamma_0_inv;

  time1 = MPI_Wtime();

    gi = 0;

    F_cart = (Complex *) calloc(nd1 * nd2 * nd1 * nd2, sizeof(Complex));
    if (F_cart == NULL) {
    fprintf(stderr, "ERROR: not enough memory for Complex F_cart! \n");
    exit(1);
   }
    
    R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
    R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
    R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
    R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

    //fprintf(file.out,"screen %3d %3d  %3d %10.4lf\n",ip,jp,gj,sqrt(R_AB_1esqrd));

    bfposi   = 0;
    bfposi1  = 0;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type1_sh[index_i];
      sheli1 = shells->type_sh[index_i];
      imax   = shells->imax_sh[index_i];
      bfposj   = 0;
      bfposj1  = 0;
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        shelj  = shells->type1_sh[index_j];
        shelj1 = shells->type_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
        shells,gaussians,job,file);

        mm  = imax + jmax + imax + jmax;
        mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;

        fgtuv = (Complex *) calloc(mm0, sizeof(Complex));
        if (fgtuv == NULL) {
        fprintf(stderr, "ERROR: not enough memory for Complex fgtuv! \n");
        exit(1);
       }

        fgtuv_max = Complex(k_zero, k_zero);

         switch (crystal->type[0]) {

           case 'C':

          expfac = four * pi / crystal->primitive_cell_volume;
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
            s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
            s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
            s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
            fac1 = sab[i4] * sab[j4];
           //for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
             //s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
             //s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
             //s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
             //p_inv_ex = pab_inv[i4] + pc_inv[j4];
             //fac1 = sab[i4] * sc[j4];
             p_fgtuv = fgtuv + i4 * mm4 + j4;
             if (fabs(q_G->sqr[0]) < 0.00001)
            *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pab_inv[j4]);
             //fprintf(file.out,"pi_vol 3C %f\n",\
             -pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pab_inv[j4]));
             int start = 0;
             if (fabs(q_G->sqr[0]) < 0.00001) start = 1;
             for (index_G = start; index_G < q_G->last_vector; index_G++) {
               fac2 = fac1 * expfac * exp(-q_G->sqr[index_G] / four * q_G->gamma_0_inv) / q_G->sqr[index_G];
               //uncomment for G only 
               //fac2 = fac1 * expfac * exp(-q_G->sqr[index_G] / four * (pab_inv[i4] + pc_inv[j4]) ) / q_G->sqr[index_G];
               dot_product = q_G->vec[index_G].comp1 * s_12.comp1 + q_G->vec[index_G].comp2 * s_12.comp2 + \
               q_G->vec[index_G].comp3 * s_12.comp3;
               for (i = 0; i <= mm; i++) {
                 for (j = 0; j <= mm; j++) {
                   for (k = 0; k <= mm; k++) {
                     ijk = i + j + k;
                     if (ijk > mm) break;
                     fgtuv_temp = fac2 * q_G->x[index_G * size + i] * q_G->y[index_G * size + j] * q_G->z[index_G * size + k] * \
                     conj(cosfactor_complex(ijk, dot_product)); // use complex conjugate for e^-iq.r
             //fgtuv_temp = Complex(k_zero,k_zero);
                     //fprintf(file.out,"fgtuv %f\n", fgtuv_temp.real());
                     p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                    *p_fgtuv += fgtuv_temp;
                    }
                   }
                  }
                 }
                }
               }
                 //for (i=0;i<mm0;i++) fprintf(file.out,"2C complex FGTUV0 %3d  %3d %3d %3d %10.4f %10.4f\n",\
                 i,index_i,index_j,index_k,fgtuv[i].real(),fgtuv[i].imag());
                 //fprintf(file.out,"\n");
                 //for (i=0;i<mm0;i++) fgtuv[i] = Complex(k_zero, k_zero);


 //           fgtuv_max = k_zero;
 //       for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
 //         for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
 //           s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
 //           s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
 //           s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
 //           fac1 = sab[i4] * sab[j4];
 //           //fprintf(file.out,"i4 %3d j4 %3d  %10.2e %10.2e %10.2e\n",i4,j4,sab[i4],sab[j4],fac1);
 //           p_inv_ex = pab_inv[i4] + pab_inv[j4];
 //           p_fgtuv = fgtuv + i4 * mm4 + j4;
 //          *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pab_inv[j4]);
 //           count = 1;
 //           for (index_S = 1; index_S < G->number_of_shells; index_S++) {
 //             fac2 = fac1 * G->EXPFAC[index_S];
 //             //double test = G->EXPFAC[index_S] * fac1 * \
 //             (k_one > G->shell_mag[index_S * 9 + mm] ? k_one : G->shell_mag[index_S * 9 + mm]);
 //             //if (fabs(test) * C_max < 1.0e-18) break;
 //             //fac2 = fac1 * expfac * exp(-G->sqr[index_S] / four * p_inv_ex) / G->sqr[index_S];
 //             for (index_G = 0; index_G < G->num[index_S]; index_G++) {
 //               dot_product = G->vec_b2[count].comp1 * s_12.comp1 + G->vec_b2[count].comp2 * s_12.comp2 + \
 //               G->vec_b2[count].comp3 * s_12.comp3;
 //               for (i = 0; i <= mm; i++) {
 //                 for (j = 0; j <= mm; j++) {
 //                   for (k = 0; k <= mm; k++) {
 //                     ijk = i + j + k;
 //                     if (ijk > mm) break;
 //                     fgtuv_temp = fac2 * G->x[count * 9 + i] * G->y[count * 9 + j] * G->z[count * 9 + k] * \
 //                     cosfactor(ijk, dot_product);
 //                     p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
 //                    *p_fgtuv += fgtuv_temp;
 //                     fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
 //                     fgtuv_max = k_one;
 //                    }
 //                   }
 //                  }
 //                 count++;
 //                }
 //               }
 //              }
 //             }
 //
 
       break;
 
       case 'S':
       case 'P':
 
       break;

       case 'M':

       if (job->taskid == 0)
       fprintf(file.out,"integrals_crystal_screen_complex routine is for periodic systems only\n");
       MPI_Finalize();
       exit(1);

       break;

      } // close switch

  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
    for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
      s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
      s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
      s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
      map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      fac1 = sab[i4] * sab[j4];
      //fprintf(file.out,"i4 %3d j4 %3d  %10.2e %10.2e %10.2e\n",i4,j4,sab[i4],sab[j4],fac1);
      p_inv_ex = pab_inv[i4] + pab_inv[j4];
      count = 0;

    //for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
      //s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
      //s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
      //s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
      //p_inv_ex = pab_inv[i4] + pc_inv[j4];
      //map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      //fac1 = sab[i4] * sc[j4];
      //count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
          //r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
          //r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
          //r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
          r_12.comp1 = s_12.comp1 + R->vec_ai[count].comp1;
          r_12.comp2 = s_12.comp2 + R->vec_ai[count].comp2;
          r_12.comp3 = s_12.comp3 + R->vec_ai[count].comp3;
          Rsqrd = r_12.comp1 * r_12.comp1 + r_12.comp2 * r_12.comp2 + r_12.comp3 * r_12.comp3;
          f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
          f000m(&en[count][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
          for (i1 = 0; i1 <= mm; i1++) {
            en[count][i1] -= em[count][i1];
           }
            shell_sum += fabs(en[count][0]);
            count++;
           } // end loop on index_R
             } // end loop on index_S
              for (index_R = 0; index_R < count; index_R++) {
                //r_12.comp1 = t_12.comp1 + R->vec_ai[index_R].comp1;
                //r_12.comp2 = t_12.comp2 + R->vec_ai[index_R].comp2;
                //r_12.comp3 = t_12.comp3 + R->vec_ai[index_R].comp3;
                r_12.comp1 = s_12.comp1 + R->vec_ai[index_R].comp1;
                r_12.comp2 = s_12.comp2 + R->vec_ai[index_R].comp2;
                r_12.comp3 = s_12.comp3 + R->vec_ai[index_R].comp3;
                dot_product = q_G->vec[0].comp1 * R->vec_ai[index_R].comp1 + q_G->vec[0].comp2 * R->vec_ai[index_R].comp2 + \
                q_G->vec[0].comp3 * R->vec_ai[index_R].comp3;
                //fac = fac1 * Complex(cos(dot_product), -sin(dot_product));
                fac = fac1 * Complex(cos(dot_product), sin(dot_product));  // sign of sin(dot ??
                //non_recursive_ftuvn(mm, index_R, f, en, &r_12);
                //replace for G functions fgtuv_temp = sab[i4] * sc[j4] * ftuvn(m,n,pm,0,&em[0][0],s_12);
                for (i = 0; i <= mm; i++) {
                  for (j = 0; j <= mm; j++) {
                    for (k = 0; k <= mm; k++) {
                      ijk = i + j + k;
                      if (ijk > mm) break;
                      fgtuv_temp = fac * ftuvn(i,j,k,0,&en[index_R][0],r_12);
           //fgtuv_temp = Complex(k_zero,k_zero);
                      //fgtuv_temp = fac * f[i][j][k][0];
                      //fgtuv_temp = fac1 * f[i][j][k][0];
                      p_fgtuv  = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                     *p_fgtuv += fgtuv_temp;
                     }
                    }
                   }
                  } 
                 } 
                } 
                 //for (i=0;i<mm0;i++) fprintf(file.out,"3C complex FGTUV1 %3d  %3d %3d %3d %10.4f %10.4f\n",\
                 i,index_i,index_j,index_k,fgtuv[i].real(),fgtuv[i].imag());
                 //for (i=0;i<mm0;i++) fprintf(file.out,"%3d  %3d %3d %3d %18.12f %18.12f\n",\
                 i,index_i,index_j,index_k,fgtuv[i].real(),fgtuv[i].imag());
                 //for (i=0;i<mm0;i++) fprintf(file.out,"%3d  %3d %3d %3d %10.4f %10.4f\n",\
                 i,index_i,index_j,index_k,fgtuv[i].real(),fgtuv[i].imag());
                 //fprintf(file.out,"\n");

 // for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
 //   for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
 //     s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
 //     s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
 //     s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
 //     map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
 //     fac1 = sab[i4] * sab[j4];
 //     //fprintf(file.out,"i4 %3d j4 %3d  %10.2e %10.2e %10.2e\n",i4,j4,sab[i4],sab[j4],fac1);
 //     p_inv_ex = pab_inv[i4] + pab_inv[j4];
 //     count = 0;
 //     for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
 //       shell_sum = k_zero;
 //       for (index_R = 0; index_R < R->num[index_S]; index_R++) {
 //         r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
 //         r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
 //         r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
 //         Rsqrd = double_vec_dot(&r_12, &r_12);
 //         f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
 //         f000m(&en[count][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
 //         for (i1 = 0; i1 <= mm; i1++) {
 //           en[count][i1] -= em[count][i1];
 //          }
 //           shell_sum += fabs(en[count][0]);
 //           count++;
 //          } // end loop on index_R
 //            } // end loop on index_S
 //             for (index_R = 0; index_R < count; index_R++) {
 //               r_12.comp1 = t_12.comp1 + R->vec_ai[index_R].comp1;
 //               r_12.comp2 = t_12.comp2 + R->vec_ai[index_R].comp2;
 //               r_12.comp3 = t_12.comp3 + R->vec_ai[index_R].comp3;
 //           for (m = 0; m <= mm; m++) {
 //             for (n = 0; n <= mm; n++) {
 //               for (pm = 0; pm <= mm; pm++) {
 //                 if (m + n + pm > mm) break;
 //                 fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&en[index_R][0],r_12);
 //                 p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
 //                *p_fgtuv += fgtuv_temp;
 //                 fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
 //                 fgtuv_max = k_one;
 //                }
 //               }
 //              }
 //             }
 //            }
 //           }
 //      //mcmurchie_davidson_3c_reversed_complex(Coulomb_cart,fgtuv,index_i,index_j,index_k,bfposi,bfposj,bfposk,nd2,nd3,\
 //      C1x,C1y,C1z,C2x,C2y,C2z,shells,shells_ax,job,file);
 //      mcmurchie_davidson_ija_complex(Coulomb_cart,fgtuv,index_i,index_j,index_k,bfposi,bfposj,bfposk,nd2,nd3,\
 //      C1x,C1y,C1z,C2x,C2y,C2z,shells,shells_ax,job,file);
 //      //for (int ggg = 0; ggg < dimtp; ggg++) fprintf(file.out,"ggg %3d %10.4lf\n",ggg,Coulomb_cart[ggg]);
 //      three_center_cartesian_to_sh_shell_ax_reversed_complex(Coulomb_cart,Coulomb,index_i,index_j,index_k,bfposi1,bfposj1,bfposk1,\
 //      bfposi,bfposj,bfposk,nd2,nd3,nd5,nd6,shells,shells_ax,job,file);
 //      //for (int ggg = 0; ggg < dimtp; ggg++) fprintf(file.out,"ggg %3d %10.4lf\n",ggg,Coulomb[ggg]);
 //      free(fgtuv);
        mcmurchie_davidson_ijij_complex(F_cart, index_i, index_j, bfposi, bfposj, nd2, C1x, C1y, C1z, fgtuv, shells, job, file);
        //mcmurchie_davidson_screen_complex(F_cart, index_i, index_j, bfposi, bfposj, nd2, C1x, C1y, C1z, fgtuv, shells, job, file);

//int num = sheli1 * shelj1 * sheli1 * shelj1;

   //for (int ggg = 0; ggg < num; ggg++) fprintf(file.out,"ggg %3d %16.10f %16.10f\n",ggg,(F_cart[ggg]).real(),(F_cart[ggg]).imag());
        free(fgtuv);
//void four_centre_cartesian_to_sh_ijij_complex(Complex*, Complex*, int, int, int, int, int, int, int, int, SHELL*, JOB_PARAM*,FILES);
        //cartesian_to_sh_ijij_complex(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        four_centre_cartesian_to_sh_ijij_complex(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        //two_centre_cartesian_to_sh_shell_ij_complex(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        //two_center_cartesian_to_sh_shell1_complex(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i
//for(i=0;i<nd1;i++){for(j=0;j<nd2;j++) { for(k=0;k<nd1;k++) { for(l=0;l<nd2;l++) { fprintf(file.out,"F %3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_sh[i*nd2*nd1*nd2+j*nd1*nd2+k*nd2+l]);}}}}
//for(i=0;i<nd3;i++) { for(j=0;j<nd4;j++) { for(k=0;k<nd3;k++) { for(l=0;l<nd4;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_sh[i*nd4*nd3*nd4+j*nd3*nd4+k*nd4+l]);}}}}

      free(F_cart);

}

void expand_screening_integral_matrix_complex(Complex *P, Complex *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int dim1, dim2, count, i, j, p, q, r, s;
  
  dim1 = 0;
  dim2 = 0;
  for (p = 0; p < pair_p->nump; p++) {
    q = pair_p->posn[p];
    rotate_permute_expand_pair_complex(p, pair_p, &P[dim1], &F[dim2], atoms, shells, symmetry, job, file);
    dim1 += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
    dim2 += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]] * pair_p->numb[p];
   }

  if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"full screening integral matrix\n");
    count = 0;
    for (p = 0; p < pair_p->nump; p++) {
      q = pair_p->posn[p];
      for (r = 0; r < pair_p->numb[p]; r++) {
        fprintf(file.out,"pair %d [%3d %3d] gj %d \n",p,pair_p->cell1[q + r],pair_p->cell2[q + r],pair_p->latt2[q + r]);
        for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
          for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
          fprintf(file.out,"%5.2e %5.2e",(F[count]).real(),(F[count]).imag());
          count++;
         }
        fprintf(file.out,"\n");
       }
      fprintf(file.out,"\n");
     }
    }
   fprintf(file.out,"\n");
  }

}

void pack_write_integrals_crystal_coulomb_ijkl(INTEGRAL_LIST *integral_list, QUAD_TRAN *quad, FILE *integrals, JOB_PARAM *job, FILES file)

{

int i, *packed_integral_indices, packed_quads[8 * quad->tot + 2];

  packed_integral_indices = (int *) malloc(2 * integral_list->num * sizeof(int));
  if (packed_integral_indices == NULL) { if (job->taskid == 0)
  fprintf(stderr, "ERROR: not enough memory for packed_integral_indices in pack_write_coulomb_2e_integrals\n");
  MPI_Finalize(); exit(1); }

  packed_quads[0] = quad->tot;
  packed_quads[1] = integral_list->num;
  for (i = 0; i < quad->tot; i++) {
  packed_quads[8 * i + 2] = quad->cell1[i];
  packed_quads[8 * i + 3] = quad->cell2[i];
  packed_quads[8 * i + 4] = quad->cell3[i];
  packed_quads[8 * i + 5] = quad->cell4[i];
  packed_quads[8 * i + 6] = quad->latt2[i];
  packed_quads[8 * i + 7] = quad->latt4[i];
  packed_quads[8 * i + 8] = quad->p[i];
  packed_quads[8 * i + 9] = quad->k[i];
 }

  for (i = 0; i < integral_list->num; i++) {
  packed_integral_indices[2 * i]     = integral_list->i[i] * 256 + integral_list->j[i];
  packed_integral_indices[2 * i + 1] = integral_list->k[i] * 256 + integral_list->l[i];
 }

  fwrite(packed_quads,sizeof(int), 8 * quad->tot + 2, integrals);
  fwrite(packed_integral_indices,sizeof(int), 2 * integral_list->num,integrals);
  fwrite(&integral_list->value[0],sizeof(double),integral_list->num,integrals);

  free(packed_integral_indices);

}

void read_unpack_integrals_crystal_coulomb_ijkl(INTEGRAL_LIST *integral_list, QUAD_TRAN *quad, FILE *integrals, JOB_PARAM *job, FILES file)

{

int i, *packed_integral_indices;
int packed_quads[8 * quad->tot];
size_t result;

  result = fread(packed_quads,sizeof(int),8 * quad->tot,integrals);
   
  for (i = 0; i < quad->tot; i++) {
  quad->cell1[i] = packed_quads[8 * i];
  quad->cell2[i] = packed_quads[8 * i + 1];
  quad->cell3[i] = packed_quads[8 * i + 2];
  quad->cell4[i] = packed_quads[8 * i + 3];
  quad->latt2[i] = packed_quads[8 * i + 4];
  quad->latt4[i] = packed_quads[8 * i + 5];
  quad->p[i] = packed_quads[8 * i + 6];
  quad->k[i] = packed_quads[8 * i + 7];
 }

  packed_integral_indices = (int *) malloc(2 * integral_list->num * sizeof(int));
  if (packed_integral_indices == NULL) { if (job->taskid == 0)
  fprintf(stderr, "ERROR: not enough memory for packed_integral_indices in pack_write_coulomb_2e_integrals\n");
  MPI_Finalize(); exit(1); }

  result = fread(packed_integral_indices,sizeof(int),2 * integral_list->num,integrals);
  result = fread(&integral_list->value[0],sizeof(double),integral_list->num,integrals);

  for (i = 0; i < integral_list->num; i++) {
  integral_list->i[i] = packed_integral_indices[2 * i] / 256;
  integral_list->j[i] = packed_integral_indices[2 * i]  - integral_list->i[i] * 256;
  integral_list->k[i] = packed_integral_indices[2 * i + 1] / 256;
  integral_list->l[i] = packed_integral_indices[2 * i + 1]  - integral_list->k[i] * 256;
 }

  free(packed_integral_indices);

}

void pack_write_integrals_crystal_exchange_ijkl(INTEGRAL_LIST *integral_list, QUAD_TRAN *quad, FILE *integrals, JOB_PARAM *job, FILES file)

{

int i, *packed_integral_indices, packed_quads[9 * quad->tot + 2];

  packed_integral_indices = (int *) malloc(2 * integral_list->num * sizeof(int));
  if (packed_integral_indices == NULL) { if (job->taskid == 0)
  fprintf(stderr, "ERROR: not enough memory for packed_integral_indices in pack_write_exchange_2e_integrals\n");
  MPI_Finalize(); exit(1); }

  packed_quads[0] = quad->tot;
  packed_quads[1] = integral_list->num;
  for (i = 0; i < quad->tot; i++) {
  packed_quads[9 * i + 2] = quad->cell1[i];
  packed_quads[9 * i + 3] = quad->cell2[i];
  packed_quads[9 * i + 4] = quad->cell3[i];
  packed_quads[9 * i + 5] = quad->cell4[i];
  packed_quads[9 * i + 6] = quad->latt2[i];
  packed_quads[9 * i + 7] = quad->latt3[i];
  packed_quads[9 * i + 8] = quad->latt4[i];
  packed_quads[9 * i + 9] = quad->p[i];
  packed_quads[9 * i + 10] = quad->k[i];
 }

  for (i = 0; i < integral_list->num; i++) {
  packed_integral_indices[2 * i]     = integral_list->i[i] * 256 + integral_list->j[i];
  packed_integral_indices[2 * i + 1] = integral_list->k[i] * 256 + integral_list->l[i];
 }

  fwrite(packed_quads,sizeof(int), 9 * quad->tot + 2, integrals);
  fwrite(packed_integral_indices,sizeof(int), 2 * integral_list->num,integrals);
  fwrite(&integral_list->value[0],sizeof(double),integral_list->num,integrals);

  free(packed_integral_indices);

}

void read_unpack_integrals_crystal_exchange_ijkl(INTEGRAL_LIST *integral_list, QUAD_TRAN *quad, FILE *integrals, JOB_PARAM *job, FILES file)

{

int i, *packed_integral_indices, packed_quads[9 * quad->tot];
size_t result;

  result = fread(packed_quads,sizeof(int),9 * quad->tot,integrals);
  
  for (i = 0; i < quad->tot; i++) {
  quad->cell1[i] = packed_quads[9 * i];
  quad->cell2[i] = packed_quads[9 * i + 1];
  quad->cell3[i] = packed_quads[9 * i + 2];
  quad->cell4[i] = packed_quads[9 * i + 3];
  quad->latt2[i] = packed_quads[9 * i + 4];
  quad->latt3[i] = packed_quads[9 * i + 5];
  quad->latt4[i] = packed_quads[9 * i + 6];
  quad->p[i] = packed_quads[9 * i + 7];
  quad->k[i] = packed_quads[9 * i + 8];
 }

  packed_integral_indices = (int *) malloc(2 * integral_list->num * sizeof(int));
  if (packed_integral_indices == NULL) { if (job->taskid == 0)
  fprintf(stderr, "ERROR: not enough memory for packed_integral_indices in pack_write_exchange_2e_integrals\n");
  MPI_Finalize(); exit(1); }

  result = fread(packed_integral_indices,sizeof(int),2 * integral_list->num,integrals);
  result = fread(&integral_list->value[0],sizeof(double),integral_list->num,integrals);

  for (i = 0; i < integral_list->num; i++) {
  integral_list->i[i] = packed_integral_indices[2 * i] / 256;
  integral_list->j[i] = packed_integral_indices[2 * i]  - integral_list->i[i] * 256;
  integral_list->k[i] = packed_integral_indices[2 * i + 1] / 256;
  integral_list->l[i] = packed_integral_indices[2 * i + 1]  - integral_list->k[i] * 256;
 }

  free(packed_integral_indices);

}

