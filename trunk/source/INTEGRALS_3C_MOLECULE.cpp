#include <mpi.h>
#include "myconstants.h"
#include "USER_DATA.h"
#include "CARTESIAN_TO_SH.h"
#include "MATRIX_UTIL.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "INCOMPLETE_GAMMA.h"
#include "MCMURCHIE_DAVIDSON.h"
#include "RECURSION.h"
#include "E_COEFFICIENTS.h"
#include "INTEGRALS_3C_MOLECULE.h"

using namespace std;

void integrals_molecule_ija(int p, TRIPLE_TRAN *triple, double *Coulomb, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int i, j, k, m, n, t, u, v, tp, up, vp;
int i1, i4, j4;
int n1, n2, n3, n4, n5, n6;
int nd1, nd2, nd3, nd4, nd5, nd6;
int dim1, dim2, dimtp, dimtf;
int tmax, umax, vmax, tpmax, upmax, vpmax;
int off1, off2, off3;
int imax, jmax, kmax;
int count;
int mm, mm0, mm1, mm2, mm3, mm4;
int shelposi, shelposj, shelposk;
int sheli, shelj, shelk;
int sheli1, shelj1, shelk1;
int index_i, index_j, index_k;
int bfposi, bfposj, bfposk;
int bfposi1, bfposj1, bfposk1;
int gausposi, gausposj, gausposk;
int tpupvp2, tpupvpsign;
double Rsqrd, p_inv_ex, fn[55];
double gamma_0, root_gamma_0;
double fgtuv_max, fgtuv_temp;
double fac1;

int ip, jp, kp, gi, gj, gk, q;
int pm;
double R_AB_1esqrd;
double C1_max, C2_max;
VECTOR_DOUBLE R_AB_1e;
double *Coulomb_cart, *p_Coulomb_cart;
VECTOR_DOUBLE R_AB, s_12;
double time1, time2;

int ijk;
int index_R, index_S, index_G;
double en[R->last_ewald_vector][55], em[R->last_ewald_vector][55], in[55], hn[55];
double sum5[16 + 1], sum2[16 + 1], x;
double fac2, fac3, Rzsqrd, expfac, x1, x2, sign, dot_product;
double D_erf[16 + 1], D_exp[16 + 1], derivative_1[16 + 1], derivative_2[16 + 1];
double shell_sum;
double *p_fgtuv;
double gamma_1;
double pi_vol = pi / crystal->primitive_cell_volume;
double gamma_1_inv = G->gamma_0_inv;
VECTOR_DOUBLE r_12, t_12, Rvec_tmp;

  gamma_0 = k_one / G->gamma_0_inv;
  root_gamma_0 = sqrt(gamma_0);
  gamma_1 = k_one / G->gamma_0_inv;

  q = triple->posn[p];
  dimtp = atoms->bfnnumb[triple->cell1[q]] * atoms->bfnnumb[triple->cell2[q]] * atoms_ax->bfnnumb[triple->cell3[q]];
  AllocateDoubleArray(&Coulomb_cart,&dimtp,job);
  ResetDoubleArray(Coulomb_cart,&dimtp);

  ip = triple->cell1[q];
  jp = triple->cell2[q];
  kp = triple->cell3[q];
  gi = 0;
  gj = 0;
  gk = 0;
  gi =  triple->latt1[q];
  gj =  triple->latt2[q];
  nd1 = atoms->bfnnumb[ip];
  nd2 = atoms->bfnnumb[jp];
  nd3 = atoms_ax->bfnnumb[kp];
  nd4 = atoms->bfnnumb_sh[ip];
  nd5 = atoms->bfnnumb_sh[jp];
  nd6 = atoms_ax->bfnnumb_sh[kp];
  R_AB_1e.comp1 = atoms->cell_vector[jp].comp1 + R->vec_ai[gj].comp1 - atoms->cell_vector[ip].comp1 - R->vec_ai[gi].comp1;
  R_AB_1e.comp2 = atoms->cell_vector[jp].comp2 + R->vec_ai[gj].comp2 - atoms->cell_vector[ip].comp2 - R->vec_ai[gi].comp2;
  R_AB_1e.comp3 = atoms->cell_vector[jp].comp3 + R->vec_ai[gj].comp3 - atoms->cell_vector[ip].comp3 - R->vec_ai[gi].comp3;
  R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
  bfposi   = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip];
  gausposi = atoms->gausposn_sh[ip];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
    sheli  = shells->cart[index_i];
    sheli1 = shells->shar[index_i];
    imax   = shells->imax_sh[index_i];
    bfposj   = 0;
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp];
    gausposj = atoms->gausposn_sh[jp];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
      shelj  = shells->cart[index_j];
      shelj1 = shells->shar[index_j];
      jmax   = shells->imax_sh[index_j];
      double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax + jmax + 1) * (imax + 1) * (jmax + 1)];
      double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax + jmax + 1) * (imax + 1) * (jmax + 1)];
      double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax + jmax + 1) * (imax + 1) * (jmax + 1)];
      double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
      E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
      shells,gaussians,job,file);
      bfposk   = 0;
      bfposk1  = 0;
      shelposk = atoms_ax->shelposn_sh[kp];
      gausposk = atoms_ax->gausposn_sh[kp];
      for (index_k = shelposk; index_k < shelposk + atoms_ax->nshel_sh[kp]; index_k++) {
        shelk  = shells_ax->cart[index_k];
        shelk1 = shells_ax->shar[index_k];
        kmax   = shells_ax->imax_sh[index_k];
        double pc_inv[shells_ax->ng_sh[index_k]];
        double sc[shells_ax->ng_sh[index_k]];
        double C2x[shells_ax->ng_sh[index_k] * (kmax + 1) * (kmax + 1)];
        double C2y[shells_ax->ng_sh[index_k] * (kmax + 1) * (kmax + 1)];
        double C2z[shells_ax->ng_sh[index_k] * (kmax + 1) * (kmax + 1)];
        E_coefficients_1c(index_k,gausposk,C2x,C2y,C2z,&C2_max,shells_ax,gaussians_ax,job,file);
        for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
        sc[j4] = pi32 / sqrt(gaussians_ax->expo_sh[gausposk + j4]) / gaussians_ax->expo_sh[gausposk + j4] * \
        gaussians_ax->c_sh[gausposk +j4];
        pc_inv[j4] = k_one / gaussians_ax->expo_sh[gausposk + j4];
       }
        mm  = imax + jmax + kmax;
        mm4 = shells_ax->ng_sh[index_k];
        mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;
        double *fgtuv;
        fgtuv = (double *) calloc(mm0, sizeof(double));
        if (fgtuv == NULL) {
        if (job->taskid == 0)
        fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
        MPI_Finalize();
        exit(1);
       }

 switch (crystal->type[0]) {

   case 'C':

  fgtuv_max = k_zero;
  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
   for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
     s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
     s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
     s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
     p_inv_ex = pab_inv[i4] + pc_inv[j4];
     fac1 = sab[i4] * sc[j4];
     p_fgtuv = fgtuv + i4 * mm4 + j4;
    *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pc_inv[j4]);
     count = 1;
     for (index_S = 1; index_S < G->number_of_shells; index_S++) {
       fac2 = fac1 * G->EXPFAC[index_S];
       for (index_G = 0; index_G < G->num[index_S]; index_G++) {
         dot_product = G->vec_b2[count].comp1 * s_12.comp1 + G->vec_b2[count].comp2 * s_12.comp2 + G->vec_b2[count].comp3 * \
         s_12.comp3;
         for (i = 0; i <= mm; i++) {
           for (j = 0; j <= mm; j++) {
             for (k = 0; k <= mm; k++) {
               ijk = i + j + k;
               if (ijk > mm) break;
               fgtuv_temp = fac2 * G->x[count * 9 + i] * G->y[count * 9 + j] * G->z[count * 9 + k] * cosfactor(ijk, dot_product);
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

  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
    for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
      s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
      s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
      s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
      p_inv_ex = pab_inv[i4] + pc_inv[j4];
      map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      fac1 = sab[i4] * sc[j4];
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
                      fgtuv_temp = fac1 * ftuvn(i,j,k,0,&em[index_R][0],r_12);
                      fgtuv[i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4] += fgtuv_temp;
                     }
                    }
                   }
                  } 
                 } 
                } 

           break;

           case 'S':

          //need to check reversal of ax and non ax atoms here to agree with triples reversed order
          fgtuv_max = k_zero;
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
           for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
             s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
             s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
             s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
             p_inv_ex = pab_inv[i4] + pc_inv[j4];
             Rsqrd = double_vec_dot(&s_12, &s_12);
             f000m(&em[0][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
                Rzsqrd = s_12.comp3 * s_12.comp3;
                fac1 = sab[i4] * sc[j4];
                //fac1 = sa[i4] * sbc[j4];
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
                      sum2[n] += G->A->a[n][k] * G->shell_mag[9 * index_S + k] * \
                      (expfac * derivative_1[n - k] + sign / expfac * derivative_2[n - k]) / \
                      G->shell_mag[9 * index_S + 1];
                      sign *= -k_one;
                     }
                    }
                  //double test = G->EXPFAC[index_S] * fac1 * \
                  (k_one > G->shell_mag[index_S * 9 + mm] ? k_one : G->shell_mag[index_S * 9 + mm]);
                  //if (fabs(test) * C_max < 1.0e-18) break;
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

          fgtuv_max = k_zero;
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
           for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
             s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
             s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
             s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
             p_inv_ex = pab_inv[i4] + pc_inv[j4];
             Rsqrd = double_vec_dot(&s_12, &s_12);
             f000m(&em[0][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
             Rzsqrd = s_12.comp3 * s_12.comp3;
             fac1 = sab[i4] * sc[j4];
             map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
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
                //non_recursive_ftuvn(mm, index_R, f, en, &r_12);
                for (i = 0; i <= mm; i++) {
                  for (j = 0; j <= mm; j++) {
                    for (k = 0; k <= mm; k++) {
                      ijk = i + j + k;
                      if (ijk > mm) break;
                      fgtuv_temp = fac1 * ftuvn(i,j,k,0,&en[index_R][0],r_12);
                      fgtuv[i * mm1 + j * mm2 + k * mm3  + i4] += fgtuv_temp;
                      //fgtuv_temp = fac1 * f[i][j][k][0];
                      //p_fgtuv  = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4;
                     //*p_fgtuv += fgtuv_temp;
                     }
                    }
                   }
                  } 
                 } 
                } 

       break;

           case 'M':

          fgtuv_max = k_zero;
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
           for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
             s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
             s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
             s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
             p_inv_ex = pab_inv[i4] + pc_inv[j4];
             Rsqrd = double_vec_dot(&s_12, &s_12);
             f000m(&em[0][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
             for (m = 0; m <= mm; m++) {
               for (n = 0; n <= mm; n++) {
                 for (pm = 0; pm <= mm; pm++) {
                   if (m + n + pm > mm) break;
                   fgtuv_temp = sab[i4] * sc[j4] * ftuvn(m,n,pm,0,&em[0][0],s_12);
                   fgtuv[m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4] += fgtuv_temp;
                   fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                   fgtuv_max = k_one;
                  }
                 }
                }
               }
              }

       break;

       } // close switch

       mcmurchie_davidson_ija(Coulomb_cart,index_i,index_j,index_k,bfposi,bfposj,bfposk,nd2,nd3,C1x,C1y,C1z,C2x,C2y,C2z, \
       fgtuv,shells,shells_ax,job,file);
       cartesian_to_sh_ija(Coulomb_cart,Coulomb,index_i,index_j,index_k,bfposi1,bfposj1,bfposk1,\
       bfposi,bfposj,bfposk,nd2,nd3,nd5,nd6,shells,shells_ax,job,file);
       free(fgtuv);
       bfposk   += shelk;
       bfposk1  += shelk1;
       gausposk += shells_ax->ng_sh[index_k];
      } // close loop over index_k
     bfposj   += shelj;
     bfposj1  += shelj1;
     gausposj += shells->ng_sh[index_j];
    } // close loop over index_j
   bfposi   += sheli;
   bfposi1  += sheli1;
   gausposi += shells->ng_sh[index_i];
  } // close loop over index_i

  free(Coulomb_cart);

  time2 = MPI_Wtime();

}

