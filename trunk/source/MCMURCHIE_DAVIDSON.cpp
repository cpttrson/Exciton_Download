
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

#include <cstdlib>
#include "myconstants.h"
#include "USER_DATA.h"
#include "MCMURCHIE_DAVIDSON.h"

using namespace std;

void mcmurchie_davidson_ij(double *F_cart, int index_i, int index_j, int bfposi, int bfposj, int nd2, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, double *fgtuv, SHELL *shells, JOB_PARAM *job, FILES file)

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
double c1fac[15][15][15];

  sheli = shells->cart[index_i];
  shelj = shells->cart[index_j];
  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  off2 = shells->ng_sh[index_i];
  off1 = (imax + 1) * off2;
  off4 = shells->ng_sh[index_j];
  off3 = (jmax + 1) * off4;
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
                          for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
                            //p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                            //c1fac[t][u][v] += C2x[n4 * off3 + tp * off4 + j4] * C2y[n5 * off3 + up * off4 + j4] * C2z[n6 * off3 + vp * off4 + j4] * \
                           *p_fgtuv * tpupvpsign;
                            c1fac[t][u][v] += fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4] * \
			    C2x[n4 * off3 + tp * off4 + j4] * \
			    C2y[n5 * off3 + up * off4 + j4] * \
			    C2z[n6 * off3 + vp * off4 + j4] * (double) tpupvpsign;
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
                            F_cart[count] += c1fac[t][u][v] * \
		            C1x[n1 * off1 + t * off2 + i4] * \
			    C1y[n2 * off1 + u * off2 + i4] * \
                            C1z[n3 * off1 + v * off2 + i4];
                           }
                          }
                         } // end t u v loop
                        } // end i4 loop
                       }
                      }

}

void mcmurchie_davidson_ij_complex(Complex *F_cart, Complex *fgtuv, int index_i, int index_j, int bfposi, int bfposj, int nd2, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, SHELL *shells, JOB_PARAM *job, FILES file)

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

  sheli = shells->cart[index_i];
  shelj = shells->cart[index_j];
  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  off2 = shells->ng_sh[index_i];
  off1 = (imax + 1) * off2;
  off4 = shells->ng_sh[index_j];
  off3 = (jmax + 1) * off4;
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
                            //c1fac[t][u][v] += C2x[n4 * off3 + tp * off4 + j4] * C2y[n5 * off3 + up * off4 + j4] * \
                            C2z[n6 * off3 + vp * off4 + j4] * tpupvpsign * \
                            fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4];
                            c1fac[t][u][v] += fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4] * \
		            C2x[n4 * off3 + tp * off4 + j4] * \
			    C2y[n5 * off3 + up * off4 + j4] * \
                            C2z[n6 * off3 + vp * off4 + j4] * (double) tpupvpsign;
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
                            F_cart[count] += c1fac[t][u][v] * \
		            C1x[n1 * off1 + t * off2 + i4] * \
			    C1y[n2 * off1 + u * off2 + i4] * \
                            C1z[n3 * off1 + v * off2 + i4];
                           }
                          }
                         } // end t u v loop
                        } // end i4 loop
                       }
                      }

}

void mcmurchie_davidson_ija(double *F_cart, int index_i, int index_j, int index_k, int bfposi, int bfposj, int bfposk, int nd2, int nd3, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, double *fgtuv, SHELL *shells, SHELL *shells_ax, JOB_PARAM *job, FILES file)

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
double c1fac[15][15][15];

  sheli = shells->cart[index_i];
  shelj = shells->cart[index_j];
  shelk = shells_ax->cart[index_k];
  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  kmax   = shells_ax->imax_sh[index_k];
  off6 = shells_ax->ng_sh[index_k];
  off5 = (kmax + 1) * off6;
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + jmax + 1) * off3;
  off1 = (jmax + 1) * off2;
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
                            //p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                            //c1fac[t][u][v] += C2x[n7 * off5 + tp * off6 + j4] * C2y[n8 * off5 + up * off6 + j4] * \
                            C2z[n9 * off5 + vp * off6 + j4] * *p_fgtuv * tpupvpsign;
                            c1fac[t][u][v] += fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4] * \
			    C2x[n7 * off5 + tp * off6 + j4] * \
			    C2y[n8 * off5 + up * off6 + j4] * \
                            C2z[n9 * off5 + vp * off6 + j4] * (double) tpupvpsign;
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
                            F_cart[count] += c1fac[t][u][v] * \
		            C1x[n1 * off1 + n4 * off2 + t * off3 + i4] * \
                            C1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
			    C1z[n3 * off1 + n6 * off2 + v * off3 + i4];
                           }
                          }
                         } // end t u v loop
                        } // end i4 loop
                       }
                      }
                     }

}

void mcmurchie_davidson_ija_complex(Complex *F_cart, Complex *fgtuv, int index_i, int index_j, int index_k, int bfposi, int bfposj, int bfposk, int nd2, int nd3, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, SHELL *shells, SHELL *shells_ax, JOB_PARAM *job, FILES file)

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
Complex c1fac[15][15][15];

  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  kmax   = shells_ax->imax_sh[index_k];
  off6 = shells_ax->ng_sh[index_k];
  off5 = (kmax + 1) * off6;
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + jmax + 1) * off3;
  off1 = (jmax + 1) * off2;
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
                            //c1fac[t][u][v] += C2x[n7 * off5 + tp * off6 + j4] * C2y[n8 * off5 + up * off6 + j4] * \
                            C2z[n9 * off5 + vp * off6 + j4] * tpupvpsign * \
                            fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4];
                            c1fac[t][u][v] += fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4] * \
		            C2x[n7 * off5 + tp * off6 + j4] * \
			    C2y[n8 * off5 + up * off6 + j4] * \
                            C2z[n9 * off5 + vp * off6 + j4] * (double) tpupvpsign;
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
                            F_cart[count] += c1fac[t][u][v] * \
			    C1x[n1 * off1 + n4 * off2 + t * off3 + i4] * \
                            C1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
			    C1z[n3 * off1 + n6 * off2 + v * off3 + i4];
                           }
                          }
                         } // end t u v loop
                        } // end i4 loop
                       }
                      }
                     }

}

void mcmurchie_davidson_ijkl(double *F_cart, int index_i, int index_j, int index_k, int index_l, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, double *fgtuv, SHELL *shells, JOB_PARAM *job, FILES file)

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
  double c1fac[15][15][15];
  double tpupvpsign, fac1;

  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  shelk = shells->type1_sh[index_k];
  shell = shells->type1_sh[index_l];
  for (i = 0; i < sheli * shelj * shelk * shell; i++) {
    F_cart[i] = k_zero;
   }
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
  dim4 = shell;
  dim3 = shelk * dim4;
  dim2 = shelj * dim3;
  dim1 = sheli * dim2;
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

void mcmurchie_davidson_ijij(double *F_cart, int index_i, int index_j, int bfposi, int bfposj, int nd2, double *C1x, double *C1y, double *C1z, double *fgtuv, SHELL *shells, JOB_PARAM *job, FILES file)

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
  double c1fac[15][15][15];
  double tpupvpsign;

  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  shelk = shells->type1_sh[index_i];
  shell = shells->type1_sh[index_j];
  for (i = 0; i < sheli * shelj * shelk * shell; i++) {
    F_cart[i] = k_zero;
   }
  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  kmax   = shells->imax_sh[index_i];
  lmax   = shells->imax_sh[index_j];
  off6 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off5 = (imax + jmax + 1) * off6;
  off4 = (jmax + 1) * off5;
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + jmax + 1) * off3;
  off1 = (jmax + 1) * off2;
  mm  = imax + jmax + imax + jmax;
  mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  mm0 = (mm + 1) * mm1;
  dim4 = shell;
  dim3 = shelk * dim4;
  dim2 = shelj * dim3;
  dim1 = sheli * dim2;
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
                            //p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                            //c1fac[t][u][v] += *p_fgtuv * C1x[n7 * off4 + n10 * off5 + tp * off6 + j4] * \
                            C1y[n8 * off4 + n11 * off5 + up * off6 + j4] * C1z[n9 * off4 + n12 * off5 + vp * off6 + j4] * tpupvpsign;
                            //p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                            c1fac[t][u][v] += fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4] * \
			    C1x[n7 * off4 + n10 * off5 + tp * off6 + j4] * \
                            C1y[n8 * off4 + n11 * off5 + up * off6 + j4] * \
			    C1z[n9 * off4 + n12 * off5 + vp * off6 + j4] * tpupvpsign;
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
                   } // end ij loop

}

void mcmurchie_davidson_ijij_complex(Complex *F_cart, int index_i, int index_j, int bfposi, int bfposj, int nd2, double *C1x, double *C1y, double *C1z, Complex *fgtuv, SHELL *shells, JOB_PARAM *job, FILES file)

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
  double tpupvpsign;
  Complex c1fac[15][15][15];

  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  shelk = shells->type1_sh[index_i];
  shell = shells->type1_sh[index_j];
  for (i = 0; i < sheli * shelj * shelk * shell; i++) {
    F_cart[i] = k_zero;
   }
  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  kmax   = shells->imax_sh[index_i];
  lmax   = shells->imax_sh[index_j];
  off6 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off5 = (imax + jmax + 1) * off6;
  off4 = (jmax + 1) * off5;
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + jmax + 1) * off3;
  off1 = (jmax + 1) * off2;
  mm  = imax + jmax + imax + jmax;
  mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  mm0 = (mm + 1) * mm1;
  dim4 = shell;
  dim3 = shelk * dim4;
  dim2 = shelj * dim3;
  dim1 = sheli * dim2;
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
                             //p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                             //c1fac[t][u][v] += *p_fgtuv * C1x[n7 * off4 + n10 * off5 + tp * off6 + j4] * \
                             C1y[n8 * off4 + n11 * off5 + up * off6 + j4] * C1z[n9 * off4 + n12 * off5 + vp * off6 + j4] * tpupvpsign ;
                             //p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                             c1fac[t][u][v] += fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4] * \
		             C1x[n7 * off4 + n10 * off5 + tp * off6 + j4] * \
                             C1y[n8 * off4 + n11 * off5 + up * off6 + j4] * \
			     C1z[n9 * off4 + n12 * off5 + vp * off6 + j4] * tpupvpsign ;
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
                   } // end ij loop

}

