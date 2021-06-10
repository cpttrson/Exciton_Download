/*
#include "PARALLEL.h"
#include "PAIRS_QUADS.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "ROTATION_OPERATORS.h"
#include "ROTATIONS_MOLECULE.h"
#include "FOURIER_TRANSFORM.h"
#include "INTEGRALS_3C_MOLECULE.h"
#include "DENSITY_FITTING_MOLECULE.h"
*/
#include <cstring>
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "INCOMPLETE_GAMMA.h"
#include "E_COEFFICIENTS.h"
#include "RECURSION.h"
#include "MCMURCHIE_DAVIDSON.h"
#include "CARTESIAN_TO_SH.h"
#include "INTEGRALS_2C_CRYSTAL.h"

using namespace std;

void two_centre_coulomb1_crystal(ComplexMatrix *V_q, REAL_LATTICE *R, Q_LATTICE *q_G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int i, j, k, m, n, p, t, u, v;
  int ip, jp, gi, gj;
  int i1, i4, j4;
  int index_i, index_j, index_G, index_R, index_S;
  int bfposi, bfposj, bfposi1, bfposj1;
  int bfposip, bfposjp;
  int imax, jmax, im, jm;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4;
  int count, size;
  int mm, mm0, mm1, mm2, mm3, mm4, ijk;
  double Rsqrd, fn[55];
  double expfac, fac2, gamma_1, dot_product;
  double en[R->last_ewald_vector][55], em[R->last_ewald_vector][55];
  double f[13][13][13][13];
  double R_AB_1esqrd, C1_max, C2_max;
  double shell_sum, pi_vol = pi / crystal->primitive_cell_volume;
  Complex fgtuv_max, fgtuv_temp, temp;
  Complex *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, r_12, s_12;

  ResetComplexMatrix(V_q);
  
  size = 4 * job->l_max + 1;
  gamma_1 = k_one / q_G->gamma_0_inv;
  //fprintf(file.out,"gamma_1 %10.4f pi_vol %10.4f\n",gamma_1,pi_vol);

  switch (crystal->type[0]) {

  case 'C':

  for (ip = 0; ip < atoms->number_of_atoms_in_unit_cell; ip++) {
  //for (jp = 0; jp < atoms->number_of_atoms_in_unit_cell; jp++) {
    for (jp = 0; jp <= ip; jp++) {

    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    Complex coulomb_cart[nd1 * nd2];
    Complex coulomb_sh[nd3 * nd4];
    for (i = 0; i < nd1 * nd2; i++) coulomb_cart[i] = Complex(k_zero, k_zero);
    for (i = 0; i < nd3 * nd4; i++) coulomb_sh[i]   = Complex(k_zero, k_zero);
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->shar[index_i];
      sheli1 = shells->cart[index_i];
      imax   = shells->imax_sh[index_i];
      //E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->shar[index_j];
        shelj1 = shells->cart[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab_fac[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double ab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1y[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1z[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C2x[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2y[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2z[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        count = 0;
        for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
            ab_inv[count] = k_one / gaussians->expo_sh[gausposi + i4] + k_one / gaussians->expo_sh[gausposj + j4];
            sab_fac[count] = pi32 * pi32 / sqrt(gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4]) / \
            gaussians->expo_sh[gausposi + i4] / gaussians->expo_sh[gausposj + j4] * gaussians->c_sh[gausposi + i4] * \
            gaussians->c_sh[gausposj + j4];
            count++;
           }
          }
        E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
        E_coefficients_1c(index_j,gausposj,C2x,C2y,C2z,&C2_max,shells,gaussians,job,file);
        mm  = imax + jmax;
        mm4 = shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;
        fgtuv = (Complex *) calloc(mm0, sizeof(Complex));
        //fgtuv = (double *) calloc(mm0, sizeof(double));
        if (fgtuv == NULL) {
        if (job->taskid == 0)
        fprintf(stderr, "ERROR: not enough memory for Complex fgtuv! \n");
        MPI_Finalize();
        exit(1);
       }

          fgtuv_max = k_zero;
          expfac = four * pi / crystal->primitive_cell_volume;
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
             s_12.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
             s_12.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
             s_12.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
             //p_inv_ex = ab_inv[i4] + pc_inv[j4];
             //fac1 = sab[i4] * sc[j4];
             // *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pc_inv[j4]);
             //p_fgtuv = fgtuv + i4 * mm4 + j4;
             if (fabs(q_G->sqr[0]) < 0.00001)
             fgtuv[i4] -= pi_vol * sab_fac[i4] * (q_G->gamma_0_inv - ab_inv[i4]);
   //fprintf(file.out,"pi vol fac %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",pi_vol, sab_fac[i4], q_G->gamma_0_inv, - ab_inv[i4],\
             - pi_vol * sab_fac[i4] * (q_G->gamma_0_inv - ab_inv[i4]),fgtuv[i4]);
int start = 0;
if (fabs(q_G->sqr[0]) < 0.00001) start = 1;

             //for (index_G = 1; index_G < q_G->last_vector; index_G++) {
             for (index_G = start; index_G < q_G->last_vector; index_G++) {
               fac2 = sab_fac[i4] * expfac * exp(-q_G->sqr[index_G] / four * q_G->gamma_0_inv) / q_G->sqr[index_G];
               // uncomment for G only fac2 = sab_fac[i4] * expfac * exp(-q_G->sqr[index_G] / four * ab_inv[i4]) / q_G->sqr[index_G];
               dot_product = q_G->vec[index_G].comp1 * s_12.comp1 + q_G->vec[index_G].comp2 * s_12.comp2 + \
               q_G->vec[index_G].comp3 * s_12.comp3;
               //fprintf(file.out,"fac2 dot %f %f\n", fac2, dot_product);
               for (i = 0; i <= mm; i++) {
                 for (j = 0; j <= mm; j++) {
                   for (k = 0; k <= mm; k++) {
                     ijk = i + j + k;
                     if (ijk > mm) break;
                     fgtuv_temp = fac2 * q_G->x[index_G * size + i] * q_G->y[index_G * size + j] * q_G->z[index_G * size + k] * \
                     conj(cosfactor_complex(ijk, dot_product));
                     //check conj(cosfactor_complex(ijk, dot_product));
                     ////cosfactor_complex(ijk, dot_product);
//fgtuv_temp = Complex(k_zero,k_zero);
                     //fprintf(file.out,"fgtuv %10.4f %10.4lf\n", (fgtuv_temp).real(),(fgtuv_temp).imag());
                     p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4;
                    *p_fgtuv += fgtuv_temp;
                     //fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                     //fgtuv_max = k_one;
                    }
                   }
                  }
                 }
                }
                 //for (i=0;i<mm0;i++) fprintf(file.out,"2C FGTUV0 %3d %10.4f %10.4f\n",i,fgtuv[i].real(),fgtuv[i].imag());
                 //for (i=0;i<mm0;i++) fprintf(file.out,"2C FGTUV0 %3d %14.8f %14.8f\n",i,(fgtuv[i]).real(),(fgtuv[i]).imag());

  Complex fac;
                //printf("q_G %10.4lf %10.4lf %10.4lf\n",q_G->vec[0].comp1,q_G->vec[0].comp2,q_G->vec[0].comp3);

  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
      s_12.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
      s_12.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
      s_12.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
      count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
//fprintf(file.out,"%3d %3d %3d %3d\n",count,index_S,index_R,R->num[index_S]); fflush(file.out);
          r_12.comp1 = s_12.comp1 + R->vec_ai[count].comp1;
          r_12.comp2 = s_12.comp2 + R->vec_ai[count].comp2;
          r_12.comp3 = s_12.comp3 + R->vec_ai[count].comp3;
          Rsqrd = r_12.comp1 * r_12.comp1 + r_12.comp2 * r_12.comp2 + r_12.comp3 * r_12.comp3;
          f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
          f000m(&en[count][0], Rsqrd / ab_inv[i4], k_one / ab_inv[i4], mm);
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
                //dot_product = q_G->vec[0].comp1 * r_12.comp1 + q_G->vec[0].comp2 * r_12.comp2 + \
                q_G->vec[0].comp3 * r_12.comp3;
                fac = sab_fac[i4] * Complex(cos(dot_product), sin(dot_product));
                //non_recursive_ftuvn(mm, index_R, f, en, &r_12);
                //replace for G functions fgtuv_temp = sab[i4] * sc[j4] * ftuvn(m,n,pm,0,&em[0][0],s_12);
                for (i = 0; i <= mm; i++) {
                  for (j = 0; j <= mm; j++) {
                    for (k = 0; k <= mm; k++) {
                      ijk = i + j + k;
                      if (ijk > mm) break;
                  fgtuv_temp = fac * ftuvn(i,j,k,0,&en[index_R][0],r_12);
                      //fgtuv_temp = fac * f[i][j][k][0];
                     //fprintf(file.out,"FGTUV %14.8lf %14.8lf %14.8lf  %14.8f %14.8lf\n", \
                     dot_product,(fac).real(),(fac).imag(),(fgtuv_temp).real(),(fgtuv_temp).imag());
                      //fgtuv_temp = sab_fac[i4] * f[i][j][k][0];
//fgtuv_temp = Complex(k_zero,k_zero);
                      p_fgtuv  = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4;
                     *p_fgtuv += fgtuv_temp;
                      //fprintf(file.out,"fgtuv1 %3d %3d %20.10lf %20.10lf\n",k,index_R,fgtuv_temp.real(),fgtuv_temp.imag());
                     }
                    }
                   }
                  } 
                 } 

                 //for (i=0;i<mm0;i++) fprintf(file.out,"%3d  %18.12f %18.12f\n",\
                 i,fgtuv[i].real(),fgtuv[i].imag());

              //for (i=0;i<mm0;i++) fprintf(file.out,"2C FGTUV1 %3d %14.8f %14.8f\n",i,(fgtuv[i]).real(),(fgtuv[i]).imag());
            //mcmurchie_davidson_2c(coulomb_cart,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,fgtuv,shells,job,file);
              //mcmurchie_davidson_2c_complex(coulomb_cart,fgtuv,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,\
              shells,job,file);
              mcmurchie_davidson_ij_complex(coulomb_cart,fgtuv,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,\
              shells,job,file);
              free(fgtuv);
              //two_center_cartesian_to_sh_shell_complex(coulomb_cart,&V_q->a[0][0],index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
	      //two_centre_cartesian_to_sh_shell_ij_complex(coulomb_cart,coulomb_sh,index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
	      four_centre_cartesian_to_sh_ijij_complex(coulomb_cart,coulomb_sh,index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
              //two_center_cartesian_to_sh_shell_complex(coulomb_cart,coulomb_sh,index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
              //two_center_cartesian_to_sh_shell_complex(coulomb_cart,&V_q->a[0][0],index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
              bfposj   += shelj;
              bfposj1  += shelj1;
              gausposj += shells->ng_sh[index_j];
             } // close loop over index_j
            bfposi   += sheli;
            bfposi1  += sheli1;
            gausposi += shells->ng_sh[index_i];
           } // close loop over index_i

              bfposip = atoms->bfnposn_sh[ip];
              bfposjp = atoms->bfnposn_sh[jp];
              count = 0;
              for (i = 0; i < nd3;i++) {
                for (j = 0; j < nd4; j++) {
                  V_q->a[bfposip + i][bfposjp + j] = coulomb_sh[count];
                  count++;
                  //fprintf(file.out,"CCART %3d %3d  %3d %3d %10.4f %10.4lf\n",\
                  bfposip,bfposjp,bfposip+i,bfposjp+j,(V_q->a[bfposip+i][bfposjp+j]).real(),(V_q->a[bfposip+i][bfposjp+j]).imag());
                 }
                }
              count = 0;
              if (ip > jp) {
              //for (i = 0; i < nd4;i++) {
                //for (j = 0; j < nd3; j++) {
              for (i = 0; i < nd3;i++) {
                for (j = 0; j < nd4; j++) {
                  V_q->a[bfposjp + j][bfposip + i] = conj(coulomb_sh[count]);
                  count++;
                 }
                }
               }

          } // close loop on jp
         } // close loop on ip

  break;

  case 'S':

  break;

  case 'P':

  break;

  case 'M':

  break;

      } // close switch

}

