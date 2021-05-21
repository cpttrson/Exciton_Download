//#include <cstdlib>
//#include "conversion_factors.h"
#include <mpi.h>
#include "myconstants.h"
#include "USER_DATA.h"
#include "CARTESIAN_TO_SH.h"
#include "MATRIX_UTIL.h"
#include "ALLOCATE_MEMORY.h"
#include "INCOMPLETE_GAMMA.h"
#include "MCMURCHIE_DAVIDSON.h"
#include "RECURSION.h"
#include "E_COEFFICIENTS.h"
#include "INTEGRALS_3C_CRYSTAL.h"

using namespace std;


void three_centre_coulomb1_reversed2_crystal_test1(int p, TRIPLE_TRAN *triple, int *start_index, Complex *Coulomb, REAL_LATTICE *R, Q_LATTICE *q_G, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int i, j, k, m, n, t, u, v, tp, up, vp;
int i1, i4, j4;
int ip, jp, kp, gi, gj, gk, q;
int nd1, nd2, nd3, nd4, nd5, nd6;
int ijk;
int size = 4 * job->l_max + 1;
int index_R, index_S, index_G;
int dimtp;
int imax, jmax, kmax;
int count, start;
int mm, mm0, mm1, mm2, mm3, mm4;
int shelposi, shelposj, shelposk;
int sheli, shelj, shelk;
int sheli1, shelj1, shelk1;
int index_i, index_j, index_k;
int bfposi, bfposj, bfposk;
int bfposi1, bfposj1, bfposk1;
int gausposi, gausposj, gausposk;
int shell_index;
double Rsqrd, p_inv_ex;
double time1, time2;
double f[13][13][13][13];
double en[R->last_ewald_vector][55], em[R->last_ewald_vector][55], in[55], hn[55];
double sum5[16 + 1], sum2[16 + 1], x;
double fac2, fac3, Rzsqrd, expfac, x1, x2, sign, dot_product;
double D_erf[16 + 1], D_exp[16 + 1], derivative_1[16 + 1], derivative_2[16 + 1];
double shell_sum;
double *p_fgtuv;
double gamma_1 = k_one / q_G->gamma_0_inv;
double pi_vol = pi / crystal->primitive_cell_volume;
double gamma_1_inv = q_G->gamma_0_inv;
double fac1;
double R_AB_1esqrd;
double C1_max, C2_max;
VECTOR_DOUBLE R_AB_1e;
VECTOR_DOUBLE R_AB, s_12;
VECTOR_DOUBLE r_12, t_12, Rvec_tmp;
Complex fgtuv_max, fgtuv_temp;
Complex *Coulomb_cart, *p_Coulomb_cart;
Complex fac;

  ip = triple->cell1[0];
  jp = triple->cell2[0];
  kp = triple->cell3[0];

  gi =  triple->latt1[0];
  gj =  triple->latt2[0];
  gk = 0;

  //printf("task %3d triple %3d %3d %3d   %3d %3d %3d\n",job->taskid,ip, jp, kp, gi,gj,gk); fflush(stdout);

  dimtp = atoms->bfnnumb[ip] * atoms->bfnnumb[jp] * atoms_ax->bfnnumb[kp];
  AllocateComplexArray(&Coulomb_cart,&dimtp,job);
  ResetComplexArray(Coulomb_cart,&dimtp);

  shell_index = 0;

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
        if (start_index[shell_index] == 0) {
          bfposk   += shells_ax->type1_sh[index_k];
          bfposk1  += shells_ax->type_sh[index_k];
          gausposk += shells_ax->ng_sh[index_k];
          shell_index++;
          continue;
         }
          shell_index++;
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
        Complex *fgtuv, *p_fgtuv;
        fgtuv = (Complex *) calloc(mm0, sizeof(Complex));
        if (fgtuv == NULL) {
        if (job->taskid == 0)
        fprintf(stderr, "ERROR: not enough memory for Complex fgtuv! \n");
        MPI_Finalize();
        exit(1);
       }

       switch (crystal->type[0]) {

         case 'C':

        fgtuv_max = k_zero;
        expfac = four * pi / crystal->primitive_cell_volume;
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
         for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
           s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
           s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
           s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
           fac1 = sab[i4] * sc[j4];
           p_fgtuv = fgtuv + i4 * mm4 + j4;
           if (fabs(q_G->sqr[0]) < 0.00001)
          *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pc_inv[j4]);
           //printf("fgtuv %3d %3d %14.8f  %14.8f %14.8f\n",\
           i4,j4,-pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pc_inv[j4]),gamma_1_inv - pab_inv[i4] - pc_inv[j4],q_G->gamma_0_inv);
           start = 0;
           if (fabs(q_G->sqr[0]) < 0.00001) start = 1;
             for (index_G = start; index_G < q_G->last_vector; index_G++) {
               fac2 = fac1 * expfac * exp(-q_G->sqr[index_G] / four / gamma_1) / q_G->sqr[index_G];
               dot_product = q_G->vec[index_G].comp1 * s_12.comp1 + q_G->vec[index_G].comp2 * s_12.comp2 + \
               q_G->vec[index_G].comp3 * s_12.comp3;
               for (i = 0; i <= mm; i++) {
                 for (j = 0; j <= mm; j++) {
                   for (k = 0; k <= mm; k++) {
                     ijk = i + j + k;
                     if (ijk > mm) break;
                     fgtuv_temp = fac2 * q_G->x[index_G * size + i] * q_G->y[index_G * size + j] * q_G->z[index_G * size + k] * \
                     conj(cosfactor_complex(ijk, dot_product)); // use complex conjugate for e^-iq.r
                     //cosfactor_complex(ijk, dot_product);
                     p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                    *p_fgtuv += fgtuv_temp;
                     //fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                     //fgtuv_max = k_one;
                    }
                   }
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
                r_12.comp1 = s_12.comp1 + R->vec_ai[index_R].comp1;
                r_12.comp2 = s_12.comp2 + R->vec_ai[index_R].comp2;
                r_12.comp3 = s_12.comp3 + R->vec_ai[index_R].comp3;
                dot_product = q_G->vec[0].comp1 * R->vec_ai[index_R].comp1 + q_G->vec[0].comp2 * R->vec_ai[index_R].comp2 + \
                q_G->vec[0].comp3 * R->vec_ai[index_R].comp3;
                fac = fac1 * Complex(cos(dot_product), sin(dot_product));
                //non_recursive_ftuvn(mm, index_R, f, en, &r_12);
                for (i = 0; i <= mm; i++) {
                  for (j = 0; j <= mm; j++) {
                    for (k = 0; k <= mm; k++) {
                      ijk = i + j + k;
                      if (ijk > mm) break;
                        fgtuv_temp = fac * ftuvn(i,j,k,0,&en[index_R][0],r_12);
                        //fgtuv_temp = fac * f[i][j][k][0];
                        p_fgtuv  = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                       *p_fgtuv += fgtuv_temp;
                      }
                     }
                    }
                   } 
                  } 
                 } 

           break;

           case 'S':

       break;

           case 'M':

       break;

       } // close switch

       mcmurchie_davidson_3c_reversed_complex(Coulomb_cart,fgtuv,index_i,index_j,index_k,bfposi,bfposj,bfposk,nd2,nd3,\
       C1x,C1y,C1z,C2x,C2y,C2z,shells,shells_ax,job,file);
       three_center_cartesian_to_sh_shell_ax_reversed_complex(Coulomb_cart,Coulomb,index_i,index_j,index_k,bfposi1,bfposj1,bfposk1,\
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

