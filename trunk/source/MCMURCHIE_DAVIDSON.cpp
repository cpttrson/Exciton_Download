/*
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <sys/types.h>
#include "mylogical.h"
#include "mycomplex.h"
#include "conversion_factors.h"
#include "PAIRS_QUADS1.h"
#include "TOOLS.h"
CHANGES2015#include "INTEGRALS1.h"
#include "MATRIX_UTIL.h"
#include "SYMMETRY.h"
#include "PAIRS_QUADS.h"
#include "LIMITS.h"
#include "INCOMPLETE_GAMMA.h" 
#include "PARALLEL.h"
#include "RECURSION.h"
#include "INCOMPLETE_GAMMA.h"
*/
#include "myconstants.h"
#include "USER_DATA.h"
#include "MCMURCHIE_DAVIDSON.h"

using namespace std;

void mcmurchie_davidson(double *F_cart, int index_i, int index_j, int index_k, int index_l, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, double *fgtuv, SHELL *shells, JOB_PARAM *job, FILES file)
//void mcmurchie_davidson(double *F_cart, double *Fc, double *Fe, int index_i, int index_j, int index_k, int index_l, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, double *fgtuv, SHELL *shells, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int n7, n8, n9, n10, n11, n12;
  int t, u, v, tp, up, vp, tpupvp2;
  int tmax, umax, vmax, tpmax, upmax, vpmax;
  int sheli, shelj, shelk, shell;
  int dim1, dim2, dim3, dim4;
  int off1, off2, off3, off4, off5, off6;
  int imax, jmax, kmax, lmax;
  int mm, mm0, mm1, mm2, mm3, mm4;
  double c1fac[9][9][9];
  //double c2fac[9][9][9][shells->ng_sh[index_i] * shells->ng_sh[index_j]];
  double tpupvpsign, p_inv_ex, fac1;

  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  kmax   = shells->imax_sh[index_k];
  lmax   = shells->imax_sh[index_l];

  off6 = shells->ng_sh[index_k] * shells->ng_sh[index_l];
  off5 = (kmax + lmax + 1) * off6;
  off4 = (lmax + 1) * off5;
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + jmax + 1) * off3;
  off1 = (jmax + 1) * off2;

  mm  = imax + jmax + kmax + lmax;
  mm4 = shells->ng_sh[index_k] * shells->ng_sh[index_l];
  mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  mm0 = (mm + 1) * mm1;

  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  shelk = shells->type1_sh[index_k];
  shell = shells->type1_sh[index_l];

  for (i = 0; i < sheli * shelj * shelk * shell; i++) {
    F_cart[i] = k_zero;
   }

  dim4 = shell;
  dim3 = shelk * dim4;
  dim2 = shelj * dim3;
  dim1 = sheli * dim2;

  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      //CHANGES2015tmax = shells->tuv[sheli][i][0] + shells->tuv[shelj][j][0];
      //umax = shells->tuv[sheli][i][1] + shells->tuv[shelj][j][1];
      //vmax = shells->tuv[sheli][i][2] + shells->tuv[shelj][j][2];
      //n1 = shells->tuv[sheli][i][0];
      //n2 = shells->tuv[sheli][i][1];
      //n3 = shells->tuv[sheli][i][2];
      //n4 = shells->tuv[shelj][j][0];
      //n5 = shells->tuv[shelj][j][1];
      //n6 = shells->tuv[shelj][j][2];
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (k = 0; k < shelk; k++) {
        for (l = 0; l < shell; l++) {
          //CHANGES2015tpmax = shells->tuv[shelk][k][0] + shells->tuv[shell][l][0];
          //upmax = shells->tuv[shelk][k][1] + shells->tuv[shell][l][1];
          //vpmax = shells->tuv[shelk][k][2] + shells->tuv[shell][l][2];
          //n7 = shells->tuv[shelk][k][0];
          //n8 = shells->tuv[shelk][k][1];
          //n9 = shells->tuv[shelk][k][2];
          //n10 = shells->tuv[shell][l][0];
          //n11 = shells->tuv[shell][l][1];
          //n12 = shells->tuv[shell][l][2];
          tpmax = shells->tuv[kmax][k][0] + shells->tuv[lmax][l][0];
          upmax = shells->tuv[kmax][k][1] + shells->tuv[lmax][l][1];
          vpmax = shells->tuv[kmax][k][2] + shells->tuv[lmax][l][2];
          n7 = shells->tuv[kmax][k][0];
          n8 = shells->tuv[kmax][k][1];
          n9 = shells->tuv[kmax][k][2];
          n10 = shells->tuv[lmax][l][0];
          n11 = shells->tuv[lmax][l][1];
          n12 = shells->tuv[lmax][l][2];
          
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
            for (t = 0; t <= tmax; t++) {
              for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                  c1fac[t][u][v] = k_zero ; 
                 }
                }
               }
                for (tp = 0; tp <= tpmax; tp++) {
                  for (up = 0; up <= upmax; up++) {
                    for (vp = 0; vp <= vpmax; vp++) {
                      tpupvpsign = -k_one;
                      tpupvp2 = tp + up + vp + 2;
                      if ((tpupvp2 / 2) * 2 == tpupvp2)
                        tpupvpsign = k_one;
                        for (j4 = 0; j4 < shells->ng_sh[index_k] * shells->ng_sh[index_l]; j4++) {
                          fac1 = C2x[n7 * off4 + n10 * off5 + tp * off6 + j4] * \
                                 C2y[n8 * off4 + n11 * off5 + up * off6 + j4] * \
                                 C2z[n9 * off4 + n12 * off5 + vp * off6 + j4] * tpupvpsign;
                          if (fabs(fac1) < 1e-20) continue;
                          for (t = 0; t <= tmax; t++) {
                            for (u = 0; u <= umax; u++) {
                              for (v = 0; v <= vmax; v++) {
                                c1fac[t][u][v] += fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4] * fac1;
                               }
                              }
                             } // end t u v loop
                            } // End j4 loop
                           }
                          }
                         } // end tp up vp loop
                for (t = 0; t <= tmax; t++) {
                  for (u = 0; u <= umax; u++) {
                    for (v = 0; v <= vmax; v++) {
                      F_cart[dim2 * i + dim3 * j + dim4 * k + l] += c1fac[t][u][v] * \
                      C1x[n1 * off1 + n4 * off2 + t * off3 + i4] * \
                      C1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                      C1z[n3 * off1 + n6 * off2 + v * off3 + i4];
                     }
                    }
                   }
                  } // End i4 loop
                 }
                }
               }
              } // end ijkl loop

}

void mcmurchie_davidson_2c_complex(Complex *F_cart, Complex *fgtuv, int index_i, int index_j, int bfposi, int bfposj, int nd2, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, SHELL *shells, JOB_PARAM *job, FILES file)

{

int i, j;
int i4, j4;
int n1, n2, n3, n4, n5, n6;
int off1, off2, off3, off4;
int mm, mm1, mm2, mm3, mm4;
int imax, jmax;
int t, u, v, tmax, umax, vmax;
int tp, up, vp, tpmax, upmax, vpmax;
int tpupvpsign, tpupvp2;
int sheli, shelj, sheli1, shelj1;
int nsheli, nshelj;
int count;
Complex c1fac[15][15][15];
//Complex *p_fgtuv;

//fprintf(file.out,"%3d %3d\n",index_i,index_j);

  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  off2 = shells->ng_sh[index_i];
  off1 = (imax + 1) * off2;
  off4 = shells->ng_sh[index_j];
  off3 = (jmax + 1) * off4;
  sheli = shells->cart[index_i];
  shelj = shells->cart[index_j];
  mm  = imax + jmax;
  mm4 = shells->ng_sh[index_j];
  mm3 = shells->ng_sh[index_i] * mm4;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      tmax =  shells->tuv[imax][i][0];
      umax =  shells->tuv[imax][i][1];
      vmax =  shells->tuv[imax][i][2];
      tpmax = shells->tuv[jmax][j][0];
      upmax = shells->tuv[jmax][j][1];
      vpmax = shells->tuv[jmax][j][2];
      //printf("%3d %3d %3d %3d    %3d %3d %3d\n",index_i,index_j,i,j,tmax,umax,vmax);
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      count = (bfposi + i) * nd2 + bfposj + j;
      for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
        for (t = 0; t <= tmax; t++) {
          for (u = 0; u <= umax; u++) {
            for (v = 0; v <= vmax; v++) {
              c1fac[t][u][v] = Complex(k_zero, k_zero) ;
             }
            }
           }
            for (t = 0; t <= tmax; t++) {
              for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                  for (tp = 0; tp <= tpmax; tp++) {
                    for (up = 0; up <= upmax; up++) {
                      for (vp = 0; vp <= vpmax; vp++) {
                        tpupvpsign = -k_one;
                        tpupvp2 = tp + up + vp + 2;
                        if ((tpupvp2 / 2) * 2 == tpupvp2)
                          tpupvpsign = k_one;
                          for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
                            //p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                            c1fac[t][u][v] += C2x[n4 * off3 + tp * off4 + j4] * C2y[n5 * off3 + up * off4 + j4] * \
                            C2z[n6 * off3 + vp * off4 + j4] * tpupvpsign * \
                            fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4];
                            //fprintf(file.out,"MM %14.8lf %14.8lf\n",\
                           (fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4]).real(),\
                           (fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4]).imag());
                            //c1fac[t][u][v] += C2x[n4 * off3 + tp * off4 + j4] * C2y[n5 * off3 + up * off4 + j4] * \
                            C2z[n6 * off3 + vp * off4 + j4] * *p_fgtuv * tpupvpsign;
                           } // end j4 loop
                          }
                         }
                        } // end tp up vp loop
                       }
                      }
                     } // end t u v loop
                      for (t = 0; t <= tmax; t++) {
                        for (u = 0; u <= umax; u++) {
                          for (v = 0; v <= vmax; v++) {
                            F_cart[count] += c1fac[t][u][v] * C1x[n1 * off1 + t * off2 + i4] * C1y[n2 * off1 + u * off2 + i4] * \
                            C1z[n3 * off1 + v * off2 + i4];
//printf("c1fac %3d %3d %3d %10.4lf %10.4lf\n",t,u,v,c1fac[t][u][v],F_cart[i * shelj + j]);
//j4 = 0;
//printf("%3d %3d %3d %10.4lf %10.4lf %10.4lf      %10.4lf\n",t,u,v,\
C1x[n1 * off1 + t * off2 + i4], C1y[n2 * off1 + u * off2 + i4],C1z[n3 * off1 + v * off2 + i4], c1fac[t][u][v]);
                           }
                          }
                         } // end t u v loop
                        } // end i4 loop
//fprintf(file.out,"%10.4lf ",F_cart[count]);
//printf("c2 %3d %3d %3d %10.4lf\n",bfposi,bfposj,count,F_cart[count]);
                       }
//fprintf(file.out,"\n");
                      }

}

void mcmurchie_davidson_screen(double *F_cart, int index_i, int index_j, int bfposi, int bfposj, int nd2, double *C1x, double *C1y, double *C1z, double *fgtuv, SHELL *shells, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int n7, n8, n9, n10, n11, n12;
  int t, u, v, tp, up, vp, tpupvp2;
  int tmax, umax, vmax, tpmax, upmax, vpmax;
  int sheli, shelj, shelk, shell;
  int dim1, dim2, dim3, dim4;
  int off1, off2, off3, off4, off5, off6;
  int imax, jmax, kmax, lmax;
  int mm, mm0, mm1, mm2, mm3, mm4;
  double c1fac[9][9][9];
  double tpupvpsign, p_inv_ex, fac1;
  double *p_F_cart, *p_fgtuv;

  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  kmax   = shells->imax_sh[index_i];
  lmax   = shells->imax_sh[index_j];
  //kmax   = shells->imax_sh[index_k];
  //lmax   = shells->imax_sh[index_l];

  //off6 = shells->ng_sh[index_k] * shells->ng_sh[index_l];
  //off5 = (kmax + lmax + 1) * off6;
  //off4 = (lmax + 1) * off5;
  off6 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off5 = (imax + jmax + 1) * off6;
  off4 = (jmax + 1) * off5;
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + jmax + 1) * off3;
  off1 = (jmax + 1) * off2;

  //mm  = imax + jmax + kmax + lmax;
  //mm4 = shells->ng_sh[index_k] * shells->ng_sh[index_l];
  mm  = imax + jmax + imax + jmax;
  mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  mm0 = (mm + 1) * mm1;

  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  shelk = shells->type1_sh[index_i];
  shell = shells->type1_sh[index_j];
  //shelk = shells->type1_sh[index_k];
  //shell = shells->type1_sh[index_l];

  p_F_cart = F_cart;
  for (i = 0; i < sheli * shelj * shelk * shell; i++) {
   *p_F_cart = k_zero;
    p_F_cart++;
   }

  dim4 = shell;
  dim3 = shelk * dim4;
  dim2 = shelj * dim3;
  dim1 = sheli * dim2;
  //dim2 = shelj;

  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (k = 0; k < shelk; k++) {
        for (l = 0; l < shell; l++) {
          //k = i; l = j;
          //kmax = imax; lmax = jmax;
          tpmax = shells->tuv[kmax][k][0] + shells->tuv[lmax][l][0];
          upmax = shells->tuv[kmax][k][1] + shells->tuv[lmax][l][1];
          vpmax = shells->tuv[kmax][k][2] + shells->tuv[lmax][l][2];
          n7 = shells->tuv[kmax][k][0];
          n8 = shells->tuv[kmax][k][1];
          n9 = shells->tuv[kmax][k][2];
          n10 = shells->tuv[lmax][l][0];
          n11 = shells->tuv[lmax][l][1];
          n12 = shells->tuv[lmax][l][2];
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
            for (t = 0; t <= tmax; t++) {
              for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                  c1fac[t][u][v] = k_zero ; 
                 }
                }
               }
             for (t = 0; t <= tmax; t++) {
               for (u = 0; u <= umax; u++) {
                 for (v = 0; v <= vmax; v++) {
                   for (tp = 0; tp <= tpmax; tp++) {
                     for (up = 0; up <= upmax; up++) {
                       for (vp = 0; vp <= vpmax; vp++) {
                         tpupvpsign = -k_one;
                         tpupvp2 = tp + up + vp + 2;
                         if ((tpupvp2 / 2) * 2 == tpupvp2)
                           tpupvpsign = k_one;
                           for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
                           //for (j4 = 0; j4 < shells->ng_sh[index_k] * shells->ng_sh[index_l]; j4++) {
                            p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                            c1fac[t][u][v] += *p_fgtuv * C1x[n7 * off4 + n10 * off5 + tp * off6 + j4] * \
                            C1y[n8 * off4 + n11 * off5 + up * off6 + j4] * C1z[n9 * off4 + n12 * off5 + vp * off6 + j4] * tpupvpsign;
                            } // End j4 loop
                           }
                          }
                         } // end tp up vp loop
                        }
                       }
                      } // end t u v loop
                     for (t = 0; t <= tmax; t++) {
                       for (u = 0; u <= umax; u++) {
                         for (v = 0; v <= vmax; v++) {
                           //p_F_cart = F_cart + (bfposi + i) * nd2 + bfposj + j;
                           p_F_cart = F_cart + dim2 * i + dim3 * j + dim4 * k + l ;
                          *p_F_cart += c1fac[t][u][v] * C1x[n1 * off1 + n4 * off2 + t * off3 + i4] * \
                           C1y[n2 * off1 + n5 * off2 + u * off3 + i4] * C1z[n3 * off1 + n6 * off2 + v * off3 + i4];
                          }
                         }
                        }
                       } // End i4 loop
                      }
                     }
                    }
                   } // end ij loop

}

void mcmurchie_davidson_3c_reversed(double *F_cart, int index_i, int index_j, int index_k, int bfposi, int bfposj, int bfposk, int nd2, int nd3, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, double *fgtuv, SHELL *shells, SHELL *shells_ax, JOB_PARAM *job, FILES file)

{

int i, j, k;
int i4, j4;
int n1, n2, n3, n4, n5, n6, n7, n8, n9;
int off1, off2, off3, off5, off6;
int mm, mm1, mm2, mm3, mm4;
int imax, jmax, kmax;
int t, u, v, tmax, umax, vmax;
int tp, up, vp, tpmax, upmax, vpmax;
int tpupvpsign, tpupvp2;
int sheli, shelj, shelk;
int count;
//double c1fac[3][3][3];
double c1fac[15][15][15];
double *p_fgtuv;

  //imax   = shells_ax->imax_sh[index_i];
  //jmax   = shells->imax_sh[index_j];
  //kmax   = shells->imax_sh[index_k];
  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  kmax   = shells_ax->imax_sh[index_k];
  //off6 = shells->ng_sh[index_j] * shells->ng_sh[index_k];
  //off5 = (jmax + kmax + 1) * off6;
  //off4 = (kmax + 1) * off5;
  //off2 = shells_ax->ng_sh[index_i];
  //off1 = (imax + 1) * off2;
  off6 = shells_ax->ng_sh[index_k];
  off5 = (kmax + 1) * off6;
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + jmax + 1) * off3;
  off1 = (jmax + 1) * off2;
//check!

  //sheli = shells_ax->cart[index_i];
  //shelj = shells->cart[index_j];
  //shelk = shells->cart[index_k];
  sheli = shells->cart[index_i];
  shelj = shells->cart[index_j];
  shelk = shells_ax->cart[index_k];
  mm  = imax + jmax + kmax;
  //mm4 = shells->ng_sh[index_j] * shells->ng_sh[index_k];
  //mm3 = shells_ax->ng_sh[index_i] * mm4;
  mm4 = shells_ax->ng_sh[index_k];
  mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  //for (i = 0; i < sheli; i++) {
    //tmax =  shells_ax->tuv[imax][i][0];
    //umax =  shells_ax->tuv[imax][i][1];
    //vmax =  shells_ax->tuv[imax][i][2];
    //n1 = shells_ax->tuv[imax][i][0];
    //n2 = shells_ax->tuv[imax][i][1];
    //n3 = shells_ax->tuv[imax][i][2];
    //for (j = 0; j < shelj; j++) {
      //for (k = 0; k < shelk; k++) {
      //tpmax = shells->tuv[jmax][j][0] + shells->tuv[kmax][k][0];
      //upmax = shells->tuv[jmax][j][1] + shells->tuv[kmax][k][1];
      //vpmax = shells->tuv[jmax][j][2] + shells->tuv[kmax][k][2];
      //n4 = shells->tuv[jmax][j][0];
      //n5 = shells->tuv[jmax][j][1];
      //n6 = shells->tuv[jmax][j][2];
      //n7 = shells->tuv[kmax][k][0];
      //n8 = shells->tuv[kmax][k][1];
      //n9 = shells->tuv[kmax][k][2];

  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
    tmax =  shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
    umax =  shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
    vmax =  shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
    n1 = shells->tuv[imax][i][0];
    n2 = shells->tuv[imax][i][1];
    n3 = shells->tuv[imax][i][2];
    n4 = shells->tuv[jmax][j][0];
    n5 = shells->tuv[jmax][j][1];
    n6 = shells->tuv[jmax][j][2];
    for (k = 0; k < shelk; k++) {
      tpmax = shells_ax->tuv[kmax][k][0];
      upmax = shells_ax->tuv[kmax][k][1];
      vpmax = shells_ax->tuv[kmax][k][2];
      n7 = shells_ax->tuv[kmax][k][0];
      n8 = shells_ax->tuv[kmax][k][1];
      n9 = shells_ax->tuv[kmax][k][2];
      count = (bfposi + i) * nd2 * nd3 + (bfposj + j) * nd3 + bfposk + k;
//fprintf(file.out,"count %3d   %3d %3d %3d  %3d %3d   %3d %3d %3d\n",count, bfposi, bfposj, bfposk,nd2,nd3,i,j,k);
//fflush(file.out);
      //for (i4 = 0; i4 < shells_ax->ng_sh[index_i]; i4++) {
      for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
        for (t = 0; t <= tmax; t++) {
          for (u = 0; u <= umax; u++) {
            for (v = 0; v <= vmax; v++) {
              c1fac[t][u][v] = k_zero ;
             }
            }
           }
            for (t = 0; t <= tmax; t++) {
              for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                  for (tp = 0; tp <= tpmax; tp++) {
                    for (up = 0; up <= upmax; up++) {
                      for (vp = 0; vp <= vpmax; vp++) {
                        tpupvpsign = -k_one;
                        tpupvp2 = tp + up + vp + 2;
                        if ((tpupvp2 / 2) * 2 == tpupvp2)
                          tpupvpsign = k_one;
                          for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
                          //for (j4 = 0; j4 < shells->ng_sh[index_j] * shells->ng_sh[index_k]; j4++) {
                            p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                            //fprintf(file.out,"%3d \n",(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4);
                            //fprintf(file.out,"%3d %3d %3d\n",n4 *off4 + n7 * off5 + tp * off6 + j4, n5 * off4+n8*off5+up*off6+j4, \
                            n6 * off4 + n9 * off5 + vp * off6 + j4);
                            //c1fac[t][u][v] += C2x[n4 * off4 + n7 * off5 + tp * off6 + j4] * \
                            C2y[n5 * off4 + n8 * off5 + up * off6 + j4] * \
                            C2z[n6 * off4 + n9 * off5 + vp * off6 + j4] * *p_fgtuv * tpupvpsign;
                            c1fac[t][u][v] += C2x[n7 * off5 + tp * off6 + j4] * C2y[n8 * off5 + up * off6 + j4] * \
                            C2z[n9 * off5 + vp * off6 + j4] * *p_fgtuv * tpupvpsign;
                           } // end j4 loop
                          }
                         }
                        } // end tp up vp loop
                       }
                      }
                     } // end t u v loop
                      for (t = 0; t <= tmax; t++) {
                        for (u = 0; u <= umax; u++) {
                          for (v = 0; v <= vmax; v++) {
                            //F_cart[count] += c1fac[t][u][v] * C1x[n1 * off1 + t * off2 + i4] * C1y[n2 * off1 + u * off2 + i4] * \
                            C1z[n3 * off1 + v * off2 + i4];
                            F_cart[count] += c1fac[t][u][v] * C1x[n1 * off1 + n4 * off2 + t * off3 + i4] * \
                            C1y[n2 * off1 + n5 * off2 + u * off3 + i4] * C1z[n3 * off1 + n6 * off2 + v * off3 + i4];
                            //fprintf(file.out,"%3d  %3d %3d %3d %10.4lf\n",count,n1 * off1 + t * off2 + i4, n2 * off1 + u*off2+i4,\
                            n3 * off1 + v * off2 + i4,F_cart[count]);
                            //fflush(file.out);
                           }
                          }
                         } // end t u v loop
                        } // end i4 loop
                       }
                      }
                     }

}

void mcmurchie_davidson_3c_reversed_complex(Complex *F_cart, Complex *fgtuv, int index_i, int index_j, int index_k, int bfposi, int bfposj, int bfposk, int nd2, int nd3, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, SHELL *shells, SHELL *shells_ax, JOB_PARAM *job, FILES file)

{

int i, j, k;
int i4, j4;
int n1, n2, n3, n4, n5, n6, n7, n8, n9;
int off1, off2, off3, off5, off6;
int mm, mm1, mm2, mm3, mm4;
int imax, jmax, kmax;
int t, u, v, tmax, umax, vmax;
int tp, up, vp, tpmax, upmax, vpmax;
int tpupvpsign, tpupvp2;
int sheli, shelj, shelk;
int count;
//double c1fac[3][3][3];
//double c1fac[15][15][15];
//double *p_fgtuv;
Complex *p_fgtuv;
Complex c1fac[15][15][15];

  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  kmax   = shells_ax->imax_sh[index_k];
  off6 = shells_ax->ng_sh[index_k];
  off5 = (kmax + 1) * off6;
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + jmax + 1) * off3;
  off1 = (jmax + 1) * off2;
//check!

  sheli = shells->cart[index_i];
  shelj = shells->cart[index_j];
  shelk = shells_ax->cart[index_k];
  mm  = imax + jmax + kmax;
  mm4 = shells_ax->ng_sh[index_k];
  mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;

  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
    tmax =  shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
    umax =  shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
    vmax =  shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
    n1 = shells->tuv[imax][i][0];
    n2 = shells->tuv[imax][i][1];
    n3 = shells->tuv[imax][i][2];
    n4 = shells->tuv[jmax][j][0];
    n5 = shells->tuv[jmax][j][1];
    n6 = shells->tuv[jmax][j][2];
    for (k = 0; k < shelk; k++) {
      tpmax = shells_ax->tuv[kmax][k][0];
      upmax = shells_ax->tuv[kmax][k][1];
      vpmax = shells_ax->tuv[kmax][k][2];
      n7 = shells_ax->tuv[kmax][k][0];
      n8 = shells_ax->tuv[kmax][k][1];
      n9 = shells_ax->tuv[kmax][k][2];
      count = (bfposi + i) * nd2 * nd3 + (bfposj + j) * nd3 + bfposk + k;
      //fprintf(file.out,"count %3d   %3d %3d %3d  %3d %3d   %3d %3d %3d\n",count, bfposi, bfposj, bfposk,nd2,nd3,i,j,k);
      //for (i4 = 0; i4 < shells_ax->ng_sh[index_i]; i4++) {
      for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
        for (t = 0; t <= tmax; t++) {
          for (u = 0; u <= umax; u++) {
            for (v = 0; v <= vmax; v++) {
              c1fac[t][u][v] = Complex(k_zero, k_zero) ;
             }
            }
           }
            for (t = 0; t <= tmax; t++) {
              for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                  for (tp = 0; tp <= tpmax; tp++) {
                    for (up = 0; up <= upmax; up++) {
                      for (vp = 0; vp <= vpmax; vp++) {
                        tpupvpsign = -k_one;
                        tpupvp2 = tp + up + vp + 2;
                        if ((tpupvp2 / 2) * 2 == tpupvp2)
                          tpupvpsign = k_one;
                          for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
                          //for (j4 = 0; j4 < shells->ng_sh[index_j] * shells->ng_sh[index_k]; j4++) {
                            //p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                            //fprintf(file.out,"%3d \n",(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4);
                            //fprintf(file.out,"%3d %3d %3d\n",n4 *off4 + n7 * off5 + tp * off6 + j4, n5 * off4+n8*off5+up*off6+j4, \
                            n6 * off4 + n9 * off5 + vp * off6 + j4);
                            //c1fac[t][u][v] += C2x[n4 * off4 + n7 * off5 + tp * off6 + j4] * \
                            C2y[n5 * off4 + n8 * off5 + up * off6 + j4] * \
                            C2z[n6 * off4 + n9 * off5 + vp * off6 + j4] * *p_fgtuv * tpupvpsign;
                            c1fac[t][u][v] += C2x[n7 * off5 + tp * off6 + j4] * C2y[n8 * off5 + up * off6 + j4] * \
                            C2z[n9 * off5 + vp * off6 + j4] * tpupvpsign * \
                            fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4];
                           } // end j4 loop
                          }
                         }
                        } // end tp up vp loop
                       }
                      }
                     } // end t u v loop
                      for (t = 0; t <= tmax; t++) {
                        for (u = 0; u <= umax; u++) {
                          for (v = 0; v <= vmax; v++) {
                            //F_cart[count] += c1fac[t][u][v] * C1x[n1 * off1 + t * off2 + i4] * C1y[n2 * off1 + u * off2 + i4] * \
                            C1z[n3 * off1 + v * off2 + i4];
                            F_cart[count] += c1fac[t][u][v] * C1x[n1 * off1 + n4 * off2 + t * off3 + i4] * \
                            C1y[n2 * off1 + n5 * off2 + u * off3 + i4] * C1z[n3 * off1 + n6 * off2 + v * off3 + i4];
                            //fprintf(file.out,"%3d  %3d %3d %3d %10.4lf\n",count,n1 * off1 + t * off2 + i4, n2 * off1 + u*off2+i4,\
                            n3 * off1 + v * off2 + i4,F_cart[count]);
                            //fflush(file.out);
                           }
                          }
                         } // end t u v loop
                        } // end i4 loop
                       }
                      }
                     }

}

void mcmurchie_davidson_screen_complex(Complex *F_cart, int index_i, int index_j, int bfposi, int bfposj, int nd2, double *C1x, double *C1y, double *C1z, Complex *fgtuv, SHELL *shells, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int n7, n8, n9, n10, n11, n12;
  int t, u, v, tp, up, vp, tpupvp2;
  int tmax, umax, vmax, tpmax, upmax, vpmax;
  int sheli, shelj, shelk, shell;
  int dim1, dim2, dim3, dim4;
  int off1, off2, off3, off4, off5, off6;
  int imax, jmax, kmax, lmax;
  int mm, mm0, mm1, mm2, mm3, mm4;
  Complex c1fac[9][9][9];
  double tpupvpsign, p_inv_ex, fac1;
  Complex *p_F_cart, *p_fgtuv;

  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  kmax   = shells->imax_sh[index_i];
  lmax   = shells->imax_sh[index_j];
  //kmax   = shells->imax_sh[index_k];
  //lmax   = shells->imax_sh[index_l];

  //off6 = shells->ng_sh[index_k] * shells->ng_sh[index_l];
  //off5 = (kmax + lmax + 1) * off6;
  //off4 = (lmax + 1) * off5;
  off6 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off5 = (imax + jmax + 1) * off6;
  off4 = (jmax + 1) * off5;
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + jmax + 1) * off3;
  off1 = (jmax + 1) * off2;

  //mm  = imax + jmax + kmax + lmax;
  //mm4 = shells->ng_sh[index_k] * shells->ng_sh[index_l];
  mm  = imax + jmax + imax + jmax;
  mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  mm0 = (mm + 1) * mm1;

  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  shelk = shells->type1_sh[index_i];
  shell = shells->type1_sh[index_j];
  //shelk = shells->type1_sh[index_k];
  //shell = shells->type1_sh[index_l];

  p_F_cart = F_cart;
  for (i = 0; i < sheli * shelj * shelk * shell; i++) {
   *p_F_cart = k_zero;
    p_F_cart++;
   }

  dim4 = shell;
  dim3 = shelk * dim4;
  dim2 = shelj * dim3;
  dim1 = sheli * dim2;
  //dim2 = shelj;

  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (k = 0; k < shelk; k++) {
        for (l = 0; l < shell; l++) {
          //k = i; l = j;
          //kmax = imax; lmax = jmax;
          tpmax = shells->tuv[kmax][k][0] + shells->tuv[lmax][l][0];
          upmax = shells->tuv[kmax][k][1] + shells->tuv[lmax][l][1];
          vpmax = shells->tuv[kmax][k][2] + shells->tuv[lmax][l][2];
          n7 = shells->tuv[kmax][k][0];
          n8 = shells->tuv[kmax][k][1];
          n9 = shells->tuv[kmax][k][2];
          n10 = shells->tuv[lmax][l][0];
          n11 = shells->tuv[lmax][l][1];
          n12 = shells->tuv[lmax][l][2];
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
            for (t = 0; t <= tmax; t++) {
              for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                  c1fac[t][u][v] = k_zero ; 
                 }
                }
               }
             for (t = 0; t <= tmax; t++) {
               for (u = 0; u <= umax; u++) {
                 for (v = 0; v <= vmax; v++) {
                   for (tp = 0; tp <= tpmax; tp++) {
                     for (up = 0; up <= upmax; up++) {
                       for (vp = 0; vp <= vpmax; vp++) {
                         tpupvpsign = -k_one;
                         tpupvp2 = tp + up + vp + 2;
                         if ((tpupvp2 / 2) * 2 == tpupvp2)
                           tpupvpsign = k_one;
                           for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
                           //for (j4 = 0; j4 < shells->ng_sh[index_k] * shells->ng_sh[index_l]; j4++) {
                             p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                             c1fac[t][u][v] += *p_fgtuv * C1x[n7 * off4 + n10 * off5 + tp * off6 + j4] * \
                             C1y[n8 * off4 + n11 * off5 + up * off6 + j4] * C1z[n9 * off4 + n12 * off5 + vp * off6 + j4] * tpupvpsign ;
                            } // End j4 loop
                           }
                          }
                         } // end tp up vp loop
                        }
                       }
                      } // end t u v loop
                     for (t = 0; t <= tmax; t++) {
                       for (u = 0; u <= umax; u++) {
                         for (v = 0; v <= vmax; v++) {
                           //p_F_cart = F_cart + (bfposi + i) * nd2 + bfposj + j;
                           p_F_cart = F_cart + dim2 * i + dim3 * j + dim4 * k + l ;
                          *p_F_cart += c1fac[t][u][v] * C1x[n1 * off1 + n4 * off2 + t * off3 + i4] * \
                           C1y[n2 * off1 + n5 * off2 + u * off3 + i4] * C1z[n3 * off1 + n6 * off2 + v * off3 + i4];
                          }
                         }
                        }
                       } // End i4 loop
                      }
                     }
                    }
                   } // end ij loop

}

