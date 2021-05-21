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
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "myconstants.h"
#include "SYMMETRY.h"
#include "PAIRS_QUADS.h"
//#include "PAIRS_QUADS1.h"
#include "TOOLS.h"
#include "LIMITS.h"
//#include "ATOM_SCF.h" // for f000m
#include "PARALLEL.h" // for f000m
//#include "INTEGRALS.h"
#include "INCOMPLETE_GAMMA.h"
#include "RECURSION.h"
#include "INTEGRALS1.h"

using namespace std;

void fock_element_1e(INT_1E *one_ints, int dim, PAIR_TRAN *pair_p, int Function[6], REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int ip, jp, gj, p, index_R, index_G, index_S;
  int i, j, k, nn, atm0, count;
  //CHANGES2015int i, j, k, mm, nn, atm0, count;
  int index_i, index_j, i4, j4;
  int dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, dim9;
  int bfposi, bfposj, bfposi1, bfposj1;
  int n1, n2, n3, n4, n5, n6;
  int t, u, v;
  //int t, u, v, tp, up, vp, tuv;
  int tmax, umax, vmax;
  //int tmax, umax, vmax, tpmax, upmax, vpmax;
  int imax, jmax;
  //int imax, jmax, kmax, lmax;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4;
  //int sheli_lim, shelj_lim;
  int nsheli, nshelj;
  int *p_i, *p_j, *p_i1, *p_j1;
  //int number_of_R_vec;
  double *p_rot1, *p_rot2;
  double r1, r2, r3, r4, r5, r6;
  double gamma_0;
  //double gamma_0, root_gamma_0;
  double pab, ab, p32;
  double SAB, KAB, R_AB_1esqrd, Rsqrd, safac, sbfac, sabfac, sabfac1;
  //double tpupvpsign;
  double fn[55], gn[55];
  //double ITOL1, ITOL2, ITOL3, ITOL4, ITOL5;
  double E1x[11][7][7], E1y[11][7][7], E1z[11][7][7];
  //CHANGES2015double E1x[8][8][8], E1y[8][8][8], E1z[8][8][8];
  double PAx, PAy, PAz, PBx, PBy, PBz;
  double expnta, expntb, GdotR;
  double shell_sum;
  double *p_ElecNuc, *p_Overlap, *p_Kinetic, *p_Momentum, *p_Grad_Grad, *p_Dipole;
  //double *p_ElecNuc, *p_Overlap, *p_Kinetic, *p_Momentum;
  //double *Kinetic, *ElecNuc, *Momentum, *Overlap;
  double *Kinetic, *ElecNuc, *Momentum, *Overlap, *Grad_Grad, *Dipole;
  VECTOR_DOUBLE Rvec_tmp, t_12;
  VECTOR_DOUBLE R_AB, R_AB_1e, r_12, s_12;

  int begin_p[job->numtasks], end_p[job->numtasks];

     mpi_begin_end(begin_p,end_p,pair_p->nump,job->numtasks,job,file);
     //printf("process int_1e %d begin %d end %d\n", job->taskid, begin_p[job->taskid],end_p[job->taskid]);

    dim3 = 0;
    for (p = 0; p < begin_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];
    dim3 += nd1 * nd2;
   }

    dim4 = dim3;
    dim5 = dim3;
    dim6 = dim3 * 3;
    dim7 = dim3;
    dim8 = dim3 * 6;

    for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    dim += atoms->bfnnumb[pair_p->cell1[pair_p->posn[p]]] * atoms->bfnnumb[pair_p->cell2[pair_p->posn[p]]] ;
   }

    Kinetic = (double *) malloc(dim * sizeof(double));
    if (Kinetic == NULL) { fprintf(stderr, "ERROR: not enough memory for overlap in\n"); exit(1); }
    p_Kinetic = Kinetic;
    for (i=0;i<dim;i++) { *p_Kinetic = k_zero; p_Kinetic++;}

    ElecNuc = (double *) malloc(dim * sizeof(double));
    if (ElecNuc == NULL) { fprintf(stderr, "ERROR: not enough memory for overlap in\n"); exit(1); }
    p_ElecNuc = ElecNuc;
    for (i=0;i<dim;i++) { *p_ElecNuc = k_zero; p_ElecNuc++;}

    Momentum = (double *) malloc(3 * dim * sizeof(double));
    if (Momentum == NULL) { fprintf(stderr, "ERROR: not enough memory for overlap in\n"); exit(1); }
    p_Momentum = Momentum;
    for (i=0;i<3 * dim;i++) { *p_Momentum = k_zero; p_Momentum++;}

    Overlap = (double *) malloc(dim * sizeof(double));
    if (Overlap == NULL) { fprintf(stderr, "ERROR: not enough memory for overlap in\n"); exit(1); }
    p_Overlap = Overlap;
    for (i=0;i<dim;i++) { *p_Overlap = k_zero; p_Overlap++;}

    Grad_Grad = (double *) malloc(6 * dim * sizeof(double));
    if (Grad_Grad == NULL) { fprintf(stderr, "ERROR: not enough memory for Grad_Grad in\n"); exit(1); }
    p_Grad_Grad = Grad_Grad;
    for (i=0;i<6 * dim;i++) { *p_Grad_Grad = k_zero; p_Grad_Grad++;}

    Dipole = (double *) malloc(3 * dim * sizeof(double));
    if (Dipole == NULL) { fprintf(stderr, "ERROR: not enough memory for Dipole in fock_element_1e\n"); exit(1); }
    p_Dipole = Dipole;
    for (i=0;i<3 * dim;i++) { *p_Dipole = k_zero; p_Dipole++;}

  // Function[0] Fock Operator
  // Function[1] Kinetic Energy Operator
  // Function[2] Electron-Nuclear Potential Energy = 1 (finite systems) = 2 (extended systems) = 3 (atoms)
  // Function[3] Momentum Operator Matrix Elements
  // Function[4] Overlap
  // Function[5] Grad_Grad xx, xy, xz, yy, yz, zz components
  // Function[6] Dipole Operator Matrix Elements

  gamma_0 = k_one / G->gamma_0_inv;
  //root_gamma_0 = sqrt(gamma_0);

  //ITOL1 = 0.0000001;
  //ITOL2 = 0.0000001;
  //ITOL3 = 0.0000001;
  //ITOL4 = 0.0000001;
  //ITOL5 = 0.0000001;

  double time1 = MPI_Wtime();

  int iii = 0;

  dim1 = 0;
  dim2 = 0;
  dim9 = 0;

  if (job->verbosity > 1)
  fprintf(file.out,"Function %d %d %d %d %d %d\n",Function[0],Function[1],Function[2],Function[3],Function[4], Function[5]);

  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {

    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    if (job->verbosity > 1)
    fprintf(file.out,"ip jp gj numb %3d %3d %3d %3d %3d\n",p,ip,jp,gj,pair_p->numb[p]) ;
    shelposi = atoms->shelposn[ip];
    gausposi = atoms->gausposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli = shells->type[index_i];
    imax  = shells->imax[index_i];
      shelposj = atoms->shelposn[jp];
      gausposj = atoms->gausposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
         shelj = shells->type[index_j];
         jmax  = shells->imax[index_j];
           //CHANGES2015 mm = imax + jmax;
        for (i4 = 0; i4 < shells->ng[index_i]; i4++) {
          for (j4 = 0; j4 < shells->ng[index_j]; j4++) {
            pab = gaussians->expo[gausposi + i4] + gaussians->expo[gausposj + j4];
            ab = gaussians->expo[gausposi + i4] * gaussians->expo[gausposj + j4];
            KAB = ab * R_AB_1esqrd / pab;
            p32 = pab * sqrt(pab);
            SAB = pi32 * exp(-KAB) / p32;
            expnta = gaussians->expo[gausposi + i4];
            expntb = gaussians->expo[gausposj + j4];
  R_AB.comp1 = (expnta * atoms->cell_vector[ip].comp1 + expntb * (atoms->cell_vector[jp].comp1 + R->vec_ai[gj].comp1)) / pab;
  R_AB.comp2 = (expnta * atoms->cell_vector[ip].comp2 + expntb * (atoms->cell_vector[jp].comp2 + R->vec_ai[gj].comp2)) / pab;
  R_AB.comp3 = (expnta * atoms->cell_vector[ip].comp3 + expntb * (atoms->cell_vector[jp].comp3 + R->vec_ai[gj].comp3)) / pab;
            PAx = R_AB.comp1 - atoms->cell_vector[ip].comp1;
            PAy = R_AB.comp2 - atoms->cell_vector[ip].comp2;
            PAz = R_AB.comp3 - atoms->cell_vector[ip].comp3;
            PBx = R_AB.comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
            PBy = R_AB.comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
            PBz = R_AB.comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
            for (i = 0; i <= imax + 2; i++) {
              for (j = 0; j <= jmax + 2; j++) {
                for (t = 0; t <= imax + jmax + 2; t++) {
                  E1x[t][i][j] = e(i, j, t, pab, PAx, PBx);
                  E1y[t][i][j] = e(i, j, t, pab, PAy, PBy);
                  E1z[t][i][j] = e(i, j, t, pab, PAz, PBz);
                }
              }
            } // end loop to set up E factors

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
                int i1 = 0, j1 = 0;
                if (sheli == 4 && i >  0) i1 =  1;
                if (sheli == 4 && i == 0) i1 = -4;
                if (shelj == 4 && j >  0) j1 =  1;
                if (shelj == 4 && j == 0) j1 = -4;
                tmax = shells->tuv[imax][i-i1][0] + shells->tuv[jmax][j-j1][0];
                umax = shells->tuv[imax][i-i1][1] + shells->tuv[jmax][j-j1][1];
                vmax = shells->tuv[imax][i-i1][2] + shells->tuv[jmax][j-j1][2];
                n1 = shells->tuv[imax][i-i1][0];
                n2 = shells->tuv[imax][i-i1][1];
                n3 = shells->tuv[imax][i-i1][2];
                n4 = shells->tuv[jmax][j-j1][0];
                n5 = shells->tuv[jmax][j-j1][1];
                n6 = shells->tuv[jmax][j-j1][2];
                //CHANGES2015fprintf(file.out,"%3d %3d   %3d %3d %3d    %3d %3d %3d\n",i,j,n1,n2,n3,n4,n5,n6);

                if (sheli == 1)
                 safac = gaussians->sc[gausposi + i4];
                if (sheli == 3)
                 safac = gaussians->pc[gausposi + i4];
                if (sheli == 4)
                 safac = gaussians->sc[gausposi + i4];
                if (sheli == 4 && i > 0)
                 safac = gaussians->pc[gausposi + i4];
                if (sheli == 6)
                 safac = gaussians->dc[gausposi + i4];
                if (sheli == 10)
                 safac = gaussians->fc[gausposi + i4];
                if (sheli == 15)
                 safac = gaussians->gc[gausposi + i4];

                if (shelj == 1)
                 sbfac = gaussians->sc[gausposj + j4];
                if (shelj == 3)
                 sbfac = gaussians->pc[gausposj + j4];
                if (shelj == 4)
                 sbfac = gaussians->sc[gausposj + j4];
                if (shelj == 4 && j > 0)
                 sbfac = gaussians->pc[gausposj + j4];
                if (shelj == 6)
                 sbfac = gaussians->dc[gausposj + j4];
                if (shelj == 10)
                 sbfac = gaussians->fc[gausposj + j4];
                if (shelj == 15)
                 sbfac = gaussians->gc[gausposj + j4];

                 sabfac = safac * sbfac * SAB;
                 sabfac1 = - sabfac / two;

        if (Function[1] == 1) {
              p_Kinetic = Kinetic + dim1 + (bfposi + i) * nd2 + bfposj + j;
              if (n4 > 1)
             *p_Kinetic += (double)(n4 * (n4 - 1)) * E1x[0][n1][n4 - 2] * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac1;
             // *p_Kinetic += (double)(n4 * (n4 + 1)) * E1x[0][n1][n4 - 2] * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac1;
              if (n5 > 1)
             *p_Kinetic += (double)(n5 * (n5 - 1)) * E1x[0][n1][n4] * E1y[0][n2][n5 - 2] * E1z[0][n3][n6] * sabfac1;
             // *p_Kinetic += (double)(n5 * (n5 + 1)) * E1x[0][n1][n4] * E1y[0][n2][n5 - 2] * E1z[0][n3][n6] * sabfac1;
              if (n6 > 1)
             *p_Kinetic += (double)(n6 * (n6 - 1)) * E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 - 2] * sabfac1;
             // *p_Kinetic += (double)(n6 * (n6 + 1)) * E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 - 2] * sabfac1;

             *p_Kinetic += (-two * gaussians->expo[gausposj + j4] * (double)(2 * (n4 + n5 + n6) + 3) * \
                            E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6] +  \
                           four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * \
                           (E1x[0][n1][n4 + 2] * E1y[0][n2][n5] * E1z[0][n3][n6] + E1x[0][n1][n4] * E1y[0][n2][n5 + 2] * \
                            E1z[0][n3][n6] + E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 + 2])) * sabfac1;
                            }

         if (Function[2] == 1) {

         switch (crystal->type[0]) {

           case 'C':
           case 'S':
           case 'P':

/*
           number_of_R_vec = R->last_vector; // need proper limit here
           r3 = k_zero;
           r5 = k_zero;
           r6 = k_zero;
            for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
              s_12.comp1 = R_AB.comp1 - atoms->cell_vector[atm0].comp1;
              s_12.comp2 = R_AB.comp2 - atoms->cell_vector[atm0].comp2;
              s_12.comp3 = R_AB.comp3 - atoms->cell_vector[atm0].comp3;
              for (t = 0; t <= tmax; t++) {
               for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                 nn = t + u + v;
                  count = 1;
                   r1 = E1x[t][n1][n4] * E1y[u][n2][n5] * E1z[v][n3][n6];
                    for (index_S = 1; index_S < G->number_of_shells; index_S++) {
                     r2 = r1 * G->EXPFAC[index_S] * (double)-atoms->atomic_number[atm0];
                      for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                       GdotR = double_vec_dot(&G->vec_b2[count], &s_12);
                        r3 += r2 * G->x[count * 9 + t] * G->y[count * 9 + u] * G->z[count * 9 + v] * cosfactor(nn, GdotR);
                        count++;
                       } // end loop over index_G
                      } // end loop over index_S
                      r4 = E1x[t][n1][n4] * E1y[u][n2][n5] * E1z[v][n3][n6];
                       for (index_R = 0; index_R < number_of_R_vec; index_R++) {
                        r_12.comp1 = s_12.comp1 + R->vec_ai[index_R].comp1;
                        r_12.comp2 = s_12.comp2 + R->vec_ai[index_R].comp2;
                        r_12.comp3 = s_12.comp3 + R->vec_ai[index_R].comp3;
                         Rsqrd = double_vec_dot(&r_12, &r_12);
                          f000m(fn, pab * Rsqrd, pab, nn);
                          f000m(gn, gamma_0 * Rsqrd, gamma_0, nn);
                          r5 += r4 * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * (double)-atoms->atomic_number[atm0];
                          fprintf(file.out,"r5 %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n",r5,ftuvn(t,u,v,0,fn,r_12),ftuvn(t, u, v, 0, gn, r_12),\
                          pab,pab*Rsqrd,sqrt(Rsqrd),r4 * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * (double)-atoms->atomic_number[atm0]);
                          } // end loop over index_R
                        }
                      }
                    } // end t u v loop
                    r6 += E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6] * -pi * (G->gamma_0_inv - k_one / pab) /
                        crystal->primitive_cell_volume * double(-atoms->atomic_number[atm0]);
                  } // end loop over atm0
                   p_ElecNuc = ElecNuc + dim1 + (bfposi + i) * nd2 + bfposj + j;
                  *p_ElecNuc += (r3 + r5 + r6) * sabfac;
                  //printf("r5 %e %e %e %e\n",r3,r5,r6,*p_ElecNuc);
*/

           //number_of_R_vec = R->max_vector; // need proper limit here
           //May2013number_of_R_vec = R->last_vector; // need proper limit here
           r3 = k_zero;
           r5 = k_zero;
           r6 = k_zero;
            for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
              s_12.comp1 = R_AB.comp1 - atoms->cell_vector[atm0].comp1;
              s_12.comp2 = R_AB.comp2 - atoms->cell_vector[atm0].comp2;
              s_12.comp3 = R_AB.comp3 - atoms->cell_vector[atm0].comp3;
              for (t = 0; t <= tmax; t++) {
               for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                 nn = t + u + v;
                  count = 1;
                   r1 = E1x[t][n1][n4] * E1y[u][n2][n5] * E1z[v][n3][n6];
                    for (index_S = 1; index_S < G->number_of_shells; index_S++) {
                     r2 = r1 * G->EXPFAC[index_S] * (double)-atoms->atomic_number[atm0];
                      for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                       GdotR = double_vec_dot(&G->vec_b2[count], &s_12);
                        r3 += r2 * G->x[count * 9 + t] * G->y[count * 9 + u] * G->z[count * 9 + v] * cosfactor(nn, GdotR);
                        count++;
                       } // end loop over index_G
                      } // end loop over index_S
                        }
                      }
                    } // end t u v loop
                      //r4 = E1x[t][n1][n4] * E1y[u][n2][n5] * E1z[v][n3][n6];

                       //for (index_R = 0; index_R < number_of_R_vec; index_R++) {
                        map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
        count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {

                        //r_12.comp1 = s_12.comp1 + R->vec_ai[index_R].comp1;
                        //r_12.comp2 = s_12.comp2 + R->vec_ai[index_R].comp2;
                        //r_12.comp3 = s_12.comp3 + R->vec_ai[index_R].comp3;
                        //r_12.comp1 = s_12.comp1 + R->vec_ai[count].comp1;
                        //r_12.comp2 = s_12.comp2 + R->vec_ai[count].comp2;
                        //r_12.comp3 = s_12.comp3 + R->vec_ai[count].comp3;
                        r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
                        r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
                        r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
                         Rsqrd = double_vec_dot(&r_12, &r_12);
              for (t = 0; t <= tmax; t++) {
               for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                  nn = t + u + v;
                  r4 = E1x[t][n1][n4] * E1y[u][n2][n5] * E1z[v][n3][n6];
                  f000m(fn, pab * Rsqrd, pab, nn);
                  f000m(gn, gamma_0 * Rsqrd, gamma_0, nn);
                  r5 += r4 * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * (double)-atoms->atomic_number[atm0];
                  shell_sum += r4 * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * (double)-atoms->atomic_number[atm0];
                  //fprintf(file.out,"r5 %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e   %3d   %3d %3d %3d\n", \
                  r5,ftuvn(t,u,v,0,fn,r_12),ftuvn(t,u,v,0,gn,r_12),\
                  pab,pab*Rsqrd,sqrt(Rsqrd),r4*(ftuvn(t, u, v, 0, fn, r_12)-ftuvn(t, u, v, 0, gn, r_12))*(double)-atoms->atomic_number[atm0],index_R,t,u,v);
                  //} // end loop over index_R
                 }
                }
               } // end t u v loop
            count++;
           } // end loop on index_R
           //fprintf(file.out,"shell sum %3d %3d %15.5e %15.5e\n",index_R,index_S,shell_sum,r5);
              } // end loop over index_S
                    r6 += E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6] * -pi * (G->gamma_0_inv - k_one / pab) /
                        crystal->primitive_cell_volume * double(-atoms->atomic_number[atm0]);
                  } // end loop over atm0
                   p_ElecNuc = ElecNuc + dim1 + (bfposi + i) * nd2 + bfposj + j;
                  //*p_ElecNuc += r5 * sabfac;
                  *p_ElecNuc += (r3 + r5 + r6) * sabfac;
                  //printf("r5 %e %e %e %e\n",r3,r5,r6,*p_ElecNuc);

           break;

           case 'M':

               for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
                 s_12.comp1 = R_AB.comp1 - atoms->cell_vector[atm0].comp1;
                 s_12.comp2 = R_AB.comp2 - atoms->cell_vector[atm0].comp2;
                 s_12.comp3 = R_AB.comp3 - atoms->cell_vector[atm0].comp3;
                 Rsqrd = double_vec_dot(&s_12, &s_12);
                 p_ElecNuc = ElecNuc + dim1 + (bfposi + i) * nd2 + bfposj + j;
                  for (t = 0; t <= tmax; t++) {
                    for (u = 0; u <= umax; u++) {
                      for (v = 0; v <= vmax; v++) {
                         nn = t + u + v;
                         f000m(fn, pab * Rsqrd, pab, nn);
                        *p_ElecNuc -= E1x[t][n1][n4] * E1y[u][n2][n5] * E1z[v][n3][n6] * sabfac *
                        ftuvn(t,u,v,0,fn,s_12) * atoms->atomic_number[atm0];
                      }
                    }
                  } // end t u v loop
                 } // end loop over atm0

           break;

      } // close switch

         } // end if(Function[2]

     if (Function[3] == 1) {

double screen = k_one;
//screen = ((tanh((R_AB.comp3 + 6.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 - 2.0) / 1.0)) / two   )  ; 
//screen = ((tanh((R_AB.comp3 + 35.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 + 12.7) / 1.0)) / two   )  ; 
//screen = ((tanh((R_AB.comp3 + 30.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 - 2.0) / 1.0)) / two   )  ; 
//screen = ((tanh((R_AB.comp3 +  6.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 - 30.0) / 1.0)) / two   )  ; 
//screen = ((tanh((R_AB.comp3 -  0.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 - 30.0) / 1.0)) / two   )  ; 
//screen = ((tanh((R_AB.comp3 -  6.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 - 30.0) / 1.0)) / two   )  ; 
//screen = ((tanh((R_AB.comp3 - 12.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 - 30.0) / 1.0)) / two   )  ; 
//screen = ((tanh((R_AB.comp3 + 25.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 + 5.0) / 1.0)) / two   )  ; 
//screen = ((tanh((R_AB.comp3 + 20.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 + 10.0) / 1.0)) / two   )  ; 
//screen = ((tanh((R_AB.comp3 + 30.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 + 0.0) / 1.0)) / two   )  ; 
//double screen = ((tanh((R_AB.comp3 + 23.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 + 9.0) / 1.0)) / two   )  ; 
//double screen = ((tanh((R_AB.comp3 + 15.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 + 1.0) / 1.0)) / two   )  ; 
//double screen = ((tanh((R_AB.comp3 + 7.0) / 1.0)) / two)  - ((  tanh((R_AB.comp3 -7.0) / 1.0)) / two   )  ; 
//double screen = ((tanh((R_AB.comp3 + 6.5) / 1.0)) / two)  - ((  tanh((R_AB.comp3 -7.5) / 1.0)) / two   )  ; 
//double screen = ((tanh((R_AB.comp3 + 3.8) / 1.0)) / two)  - ((  tanh((R_AB.comp3 - 10.2) / 1.0)) / two   )  ; 

         p_Momentum = Momentum + dim2 + (bfposi + i) * nd2 + bfposj + j;
       if (n4 >= 1)
        //*p_Momentum += (double(n4) * E1x[0][n1][n4 - 1] * E1y[0][n2][n5] * E1z[0][n3][n6]) * sabfac;
        //*p_Momentum -= two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 + 1] * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac;
        *p_Momentum += (double(n4) * E1x[0][n1][n4 - 1] * E1y[0][n2][n5] * E1z[0][n3][n6]) * sabfac * screen;
        *p_Momentum -= two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 + 1] * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac * screen;
         p_Momentum = Momentum + dim2 + nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
       if (n5 >= 1)
        //*p_Momentum += (double(n5) * E1x[0][n1][n4] * E1y[0][n2][n5 - 1] * E1z[0][n3][n6]) * sabfac;
        //*p_Momentum -= two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5 + 1] * E1z[0][n3][n6] * sabfac;
        *p_Momentum += (double(n5) * E1x[0][n1][n4] * E1y[0][n2][n5 - 1] * E1z[0][n3][n6]) * sabfac * screen;
        *p_Momentum -= two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5 + 1] * E1z[0][n3][n6] * sabfac * screen;
         p_Momentum = Momentum + dim2 + 2 * nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
       if (n6 >= 1)
        //*p_Momentum += (double(n6) * E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 - 1]) * sabfac;
        //*p_Momentum -= two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 + 1] * sabfac;
        *p_Momentum += (double(n6) * E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 - 1]) * sabfac * screen;
        *p_Momentum -= two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 + 1] * sabfac * screen;

                } // end if (Function[3]

        if (Function[4] == 1) {

                 p_Overlap = Overlap + dim1 + (bfposi + i) * nd2 + bfposj + j;
                *p_Overlap += E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac;

               }

        if (Function[5] == 1) {

            p_Grad_Grad = Grad_Grad + dim9 + (bfposi + i) * nd2 + bfposj + j;
           if (n4 > 1)
            *p_Grad_Grad += (double)(n4 * (n4 - 1)) * E1x[0][n1][n4 - 2] * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac;
            *p_Grad_Grad += (double)(n4 - (n4 + 1) * two * gaussians->expo[gausposj + j4]) * E1x[0][n1][n4] * E1y[0][n2][n5] * \
                            E1z[0][n3][n6] * sabfac;
            *p_Grad_Grad += four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 + 2] * E1y[0][n2][n5] * \
                            E1z[0][n3][n6] * sabfac;

            p_Grad_Grad = Grad_Grad + dim9 + nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
           if (n4 > 0 && n5 > 0)
            *p_Grad_Grad += (double)(n4 * n5) * E1x[0][n1][n4 - 1] * E1y[0][n2][n5 - 1] * E1z[0][n3][n6] * sabfac;
           if (n4 > 0)
            *p_Grad_Grad += (double)(n4) * E1x[0][n1][n4 - 1] * E1y[0][n2][n5 + 1] * E1z[0][n3][n6] * sabfac;
           if (n5 > 0)
            *p_Grad_Grad -= (double)(n5) * two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 + 1] * E1y[0][n2][n5 - 1] * \
                            E1z[0][n3][n6] * sabfac;
            *p_Grad_Grad += four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 + 1] * \
                            E1y[0][n2][n5 + 1] * E1z[0][n3][n6] * sabfac;

            p_Grad_Grad = Grad_Grad + dim9 + 2 * nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
           if (n4 > 0 && n6 > 0)
            *p_Grad_Grad += (double)(n4 * n6) * E1x[0][n1][n4 - 1] * E1y[0][n2][n5] * E1z[0][n3][n6 - 1] * sabfac;
           if (n4 > 0)
            *p_Grad_Grad += (double)(n4) * E1x[0][n1][n4 - 1] * E1y[0][n2][n5] * E1z[0][n3][n6 + 1] * sabfac;
           if (n6 > 0)
            *p_Grad_Grad -= (double)(n6) * two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 + 1] * E1y[0][n2][n5] * \
                            E1z[0][n3][n6 + 1] * sabfac;
            *p_Grad_Grad += four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 + 1] * E1y[0][n2][n5] * \
                            E1z[0][n3][n6 + 1] * sabfac;

/*
            p_Grad_Grad = Grad_Grad + dim9 + 3 * nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
           if (n4 > 0 && n5 > 0)
            *p_Grad_Grad += (double)(n4 * n5) * E1x[0][n1][n4 - 1] * E1y[0][n2][n5 - 1] * E1z[0][n3][n6] * sabfac;
           if (n5 > 0)
            *p_Grad_Grad += (double)(n5) * E1x[0][n1][n4 + 1] * E1y[0][n2][n5 - 1] * E1z[0][n3][n6] * sabfac;
           if (n4 > 0)
            *p_Grad_Grad -= (double)(n4) * two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 - 1] * E1y[0][n2][n5 + 1] * \
                            E1z[0][n3][n6] * sabfac;
            *p_Grad_Grad += four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 + 1] * \
                            E1y[0][n2][n5 + 1] * E1z[0][n3][n6] * sabfac;
*/

            p_Grad_Grad = Grad_Grad + dim9 + 3 * nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
           if (n5 > 1)
            *p_Grad_Grad += (double)(n5 * (n5 - 1)) * E1x[0][n1][n4] * E1y[0][n2][n5 - 2] * E1z[0][n3][n6] * sabfac;
            *p_Grad_Grad += (double)(n5 - (n5 + 1) * two * gaussians->expo[gausposj + j4]) * E1x[0][n1][n4] * E1y[0][n2][n5] * \
                            E1z[0][n3][n6] * sabfac;
            *p_Grad_Grad += four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5 + 2] * \
                            E1z[0][n3][n6] * sabfac;

            p_Grad_Grad = Grad_Grad + dim9 + 4 * nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
           if (n5 > 0 && n6 > 0)
            *p_Grad_Grad += (double)(n5 * n6) * E1x[0][n1][n4] * E1y[0][n2][n5 - 1] * E1z[0][n3][n6 - 1] * sabfac;
           if (n5 > 0)
            *p_Grad_Grad += (double)(n5) * E1x[0][n1][n4] * E1y[0][n2][n5 - 1] * E1z[0][n3][n6 + 1] * sabfac;
           if (n6 > 0)
            *p_Grad_Grad -= (double)(n6) * two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5 + 1] * \
                            E1z[0][n3][n6 - 1] * sabfac;
            *p_Grad_Grad += four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5 + 1] * \
                            E1z[0][n3][n6 + 1] * sabfac;

/*
            p_Grad_Grad = Grad_Grad + dim9 + 6 * nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
           if (n6 > 0 && n4 > 0)
            *p_Grad_Grad += (double)(n6 * n4) * E1x[0][n1][n4 - 1] * E1y[0][n2][n5] * E1z[0][n3][n6 - 1] * sabfac;
           if (n6 > 0)
            *p_Grad_Grad += (double)(n4) * E1x[0][n1][n4 + 1] * E1y[0][n2][n5] * E1z[0][n3][n6 - 1] * sabfac;
           if (n4 > 0)
            *p_Grad_Grad -= (double)(n4) * two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 - 1] * E1y[0][n2][n5] * \
                            E1z[0][n3][n6 + 1] * sabfac;
            *p_Grad_Grad += four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * E1x[0][n1][n4 + 1] * E1y[0][n2][n5] * \
                            E1z[0][n3][n6 + 1] * sabfac;

            p_Grad_Grad = Grad_Grad + dim9 + 7 * nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
           if (n6 > 0 && n5 > 0)
            *p_Grad_Grad += (double)(n6 * n5) * E1x[0][n1][n4] * E1y[0][n2][n5 - 1] * E1z[0][n3][n6 - 1] * sabfac;
           if (n6 > 0)
            *p_Grad_Grad += (double)(n5) * E1x[0][n1][n4] * E1y[0][n2][n5 + 1] * E1z[0][n3][n6 - 1] * sabfac;
           if (n5 > 0)
            *p_Grad_Grad -= (double)(n5) * two * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5 - 1] * \
                            E1z[0][n3][n6 + 1] * sabfac;
            *p_Grad_Grad += four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5 + 1] * \
                            E1z[0][n3][n6 + 1] * sabfac;
*/

            p_Grad_Grad = Grad_Grad + dim9 + 5 * nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
           if (n6 > 1)
            *p_Grad_Grad += (double)(n6 * (n6 - 1)) * E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 - 2] * sabfac;
            *p_Grad_Grad += (double)(n6 - (n6 + 1) * two * gaussians->expo[gausposj + j4]) * E1x[0][n1][n4] * E1y[0][n2][n5] * \
                            E1z[0][n3][n6] * sabfac;
            *p_Grad_Grad += four * gaussians->expo[gausposj + j4] * gaussians->expo[gausposj + j4] * E1x[0][n1][n4] * E1y[0][n2][n5] * \
                            E1z[0][n3][n6 + 2] * sabfac;

          } // end if (Function[5]

     if (Function[6] == 1) {

         p_Dipole = Dipole + dim2 + (bfposi + i) * nd2 + bfposj + j;
        *p_Dipole += (E1x[1][n1][n4] + R_AB.comp1 * E1x[0][n1][n4]) * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac;
        //*p_Dipole += (E1x[0][n1 + 1][n4] + PAx * E1x[0][n1][n4]) * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac;

         p_Dipole = Dipole + dim2 + nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
        *p_Dipole += E1x[0][n1][n4] * (E1y[1][n2][n5] + R_AB.comp2 * E1y[0][n2][n5]) * E1z[0][n3][n6] * sabfac;
        //*p_Dipole += E1x[0][n1][n4] * (E1y[0][n2 + 1][n5] + PAy * E1y[0][n2][n5]) * E1z[0][n3][n6] * sabfac;

         p_Dipole = Dipole + dim2 + 2 * nd1 * nd2 + (bfposi + i) * nd2 + bfposj + j;
        *p_Dipole += E1x[0][n1][n4] * E1y[0][n2][n5] * (E1z[1][n3][n6] + R_AB.comp3 * E1z[0][n3][n6]) * sabfac;
        //*p_Dipole += E1x[0][n1][n4] * E1y[0][n2][n5] * (E1z[0][n3 + 1][n6] + PAz * E1z[0][n3][n6]) * sabfac;

                } // end if (Function[6]

              } // close loop over j
            } // close loop over i
          } // End j4 loop
        } // End i4 loop
        bfposj += shelj;
       gausposj += shells->ng[index_j];
      } // close loop over index_j
      bfposi += sheli;
     gausposi += shells->ng[index_i];
    } // close loop over index_i

if (job->verbosity > 1) {
p_Kinetic = Dipole + dim2;
fprintf(file.out,"Dipole\n");
for(i=0;i<atoms->bfnnumb[0];i++) {
for(j=0;j<atoms->bfnnumb[0];j++) {
  fprintf(file.out,"%6.2lf",*p_Kinetic);
   printf("Dipole proc %d    %d %d   %d %d %16.8lf\n",job->taskid,pair_p->cell1[p],pair_p->cell2[p], i,j,*p_Kinetic);
  p_Kinetic++;
}
fprintf(file.out,"\n");}
}

    dim1 +=     nd1 * nd2;
    dim2 += 3 * nd1 * nd2;
    dim9 += 6 * nd1 * nd2;

  } // close loop on p

//if (job->verbosity >= 1) {
//p_Kinetic = Dipole;
//fprintf(file.out,"Dipole\n");
//for(i=0;i<atoms->bfnnumb[0];i++) {
//for(j=0;j<atoms->bfnnumb[0];j++) {
  //fprintf(file.out,"%11.5lf",*p_Kinetic);
   //printf("Dipole proc %d %d %d %16.8lf\n",job->taskid, i,j,*p_Kinetic);
  //p_Kinetic++;
//}
//fprintf(file.out,"\n");}
//}

 if (Function[0] == 1) {
  dim2 = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];
    nd3 = atoms->bfnnumb[ip];
    nd4 = atoms->bfnnumb[jp];
    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli1 = shells->type[index_i];
    sheli  = shells->type1[index_i];
    nsheli = *(shells->num_ij + shells->ord[index_i]);
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
      shelj1 = shells->type[index_j];
      shelj  = shells->type1[index_j];
      nshelj = *(shells->num_ij + shells->ord[index_j]);
        p_i1 = shells->ind_i + shells->opp[index_i];
        p_i  = shells->ind_j + shells->opp[index_i];
        p_rot1 = shells->rot + shells->opp[index_i];
        for (i = 0; i < nsheli; i++) {
          p_j1 = shells->ind_i + shells->opp[index_j];
          p_j  = shells->ind_j + shells->opp[index_j];
          p_rot2 = shells->rot + shells->opp[index_j];
          for (j = 0; j < nshelj; j++) {
            p_Kinetic = Kinetic + dim2 + (bfposi1 + *p_i1) * nd4 + bfposj1 + *p_j1;
            p_ElecNuc = ElecNuc + dim2 + (bfposi1 + *p_i1) * nd4 + bfposj1 + *p_j1;
//May2013
            //one_ints->Fock[dim3 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_Kinetic ;
            one_ints->Fock[dim3 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * (*p_Kinetic + *p_ElecNuc);
            //one_ints->Fock[dim3 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_ElecNuc;
//May2013
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj   += shelj;
        bfposj1  += shelj1;
      } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
    } // close loop over index_i

    if (job->verbosity > 1) {
      count = 0;
      fprintf(file.out,"Fock 1e %3d %3d %3d \n",ip,jp,gj);
      //fprintf(file.out,"Fock 1e\n");
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%16.10lf",one_ints->Fock[dim3 + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }

    dim3 += nd1 * nd2;
    dim2 += nd3 * nd4;
  } // close loop on p
 }

 if (Function[1] == 1) {
  dim2 = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];
    nd3 = atoms->bfnnumb[ip];
    nd4 = atoms->bfnnumb[jp];
    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli1 = shells->type[index_i];
    sheli  = shells->type1[index_i];
    nsheli = *(shells->num_ij + shells->ord[index_i]);
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
      shelj1 = shells->type[index_j];
      shelj  = shells->type1[index_j];
      nshelj = *(shells->num_ij + shells->ord[index_j]);
        p_i1 = shells->ind_i + shells->opp[index_i];
        p_i  = shells->ind_j + shells->opp[index_i];
        p_rot1 = shells->rot + shells->opp[index_i];
        for (i = 0; i < nsheli; i++) {
          p_j1 = shells->ind_i + shells->opp[index_j];
          p_j  = shells->ind_j + shells->opp[index_j];
          p_rot2 = shells->rot + shells->opp[index_j];
          for (j = 0; j < nshelj; j++) {
            p_Kinetic = Kinetic + dim2 + (bfposi1 + *p_i1) * nd4 + bfposj1 + *p_j1;
//May2013
            one_ints->Kinetic[dim4 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_Kinetic;
//May2013
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj   += shelj;
        bfposj1  += shelj1;
      } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
    } // close loop over index_i

    if (job->verbosity > 1) {
      count = 0;
      fprintf(file.out,"Kinetic\n");
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.3lf",one_ints->Kinetic[dim4 + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }

    dim4+= nd1 * nd2;
    dim2 += nd3 * nd4;
  } // close loop on p
 }

 if (Function[2] == 1) {
  dim2 = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];
    nd3 = atoms->bfnnumb[ip];
    nd4 = atoms->bfnnumb[jp];
    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli1 = shells->type[index_i];
    sheli  = shells->type1[index_i];
    nsheli = *(shells->num_ij + shells->ord[index_i]);
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
      shelj1 = shells->type[index_j];
      shelj  = shells->type1[index_j];
      nshelj = *(shells->num_ij + shells->ord[index_j]);
        p_i1 = shells->ind_i + shells->opp[index_i];
        p_i  = shells->ind_j + shells->opp[index_i];
        p_rot1 = shells->rot + shells->opp[index_i];
        for (i = 0; i < nsheli; i++) {
          p_j1 = shells->ind_i + shells->opp[index_j];
          p_j  = shells->ind_j + shells->opp[index_j];
          p_rot2 = shells->rot + shells->opp[index_j];
          for (j = 0; j < nshelj; j++) {
            p_ElecNuc = ElecNuc + dim2 + (bfposi1 + *p_i1) * nd4 + bfposj1 + *p_j1;
//May2013
            //one_ints->ElecNuc[dim5 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * 0.0 * *p_ElecNuc;
            one_ints->ElecNuc[dim5 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_ElecNuc;
      //fprintf(file.out,"%3d %3d   %3d %3d   %12.4e %12.4e \n",(bfposi1 + *p_i1) * nd4 , bfposj1 + *p_j1,\
      (bfposi + *p_i) * nd2 , bfposj + *p_j, *p_rot1 * *p_rot2, *p_ElecNuc);
//May2013
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj   += shelj;
        bfposj1  += shelj1;
      } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
    } // close loop over index_i

    if (job->verbosity > 1) {
      count = 0;
      fprintf(file.out,"ElecNuc\n");
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.3lf",one_ints->ElecNuc[dim5 + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }

    dim5+= nd1 * nd2;
    dim2 += nd3 * nd4;
  } // close loop on p
 }

 if (Function[3] == 1) {
  dim2 = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];
    nd3 = atoms->bfnnumb[ip];
    nd4 = atoms->bfnnumb[jp];
    for (k = 0; k < 3; k++) {
    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli1 = shells->type[index_i];
    sheli  = shells->type1[index_i];
    nsheli = *(shells->num_ij + shells->ord[index_i]);
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
//fprintf(file.out,"%3d %3d %3d   %3d \n",p,index_i,index_j,dim6);
      shelj1 = shells->type[index_j];
      shelj  = shells->type1[index_j];
      nshelj = *(shells->num_ij + shells->ord[index_j]);
        p_i1 = shells->ind_i + shells->opp[index_i];
        p_i  = shells->ind_j + shells->opp[index_i];
        p_rot1 = shells->rot + shells->opp[index_i];
        for (i = 0; i < nsheli; i++) {
          p_j1 = shells->ind_i + shells->opp[index_j];
          p_j  = shells->ind_j + shells->opp[index_j];
          p_rot2 = shells->rot + shells->opp[index_j];
          for (j = 0; j < nshelj; j++) {
         p_Momentum = Momentum + dim2 + k * nd3 * nd4 + (bfposi1 + *p_i1) * nd4 + bfposj1 + *p_j1;
         one_ints->Momentum[dim6 + k * nd1 * nd2 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_Momentum;
//fprintf(file.out,"%3d %3d %3d   %3d %3d %3d\n",p,index_i,index_j,dim6,k * nd1 * nd2 + (bfposi + *p_i) * nd2 + bfposj + *p_j,dim6 + k * nd1 * nd2 + (bfposi + *p_i) * nd2 + bfposj + *p_j);
     //fprintf(file.out,"%3d %3d   %3d %3d   %3d %3d   %12.4e %12.4e \n",bfposj1, *p_j1, (bfposi1 + *p_i1) * nd2 , bfposj1 + *p_j1,\
      (bfposi + *p_i) * nd4 , bfposj + *p_j, *p_rot1 * *p_rot2, *p_Momentum);
     //fprintf(file.out,"%3d   %3d %3d     %3d   %3d %3d   %12.4e %12.4e \n", k*nd3*nd4,\
      (bfposi1 + *p_i1) * nd4 , bfposj1 + *p_j1, k*nd1*nd2,(bfposi + *p_i) * nd2,bfposj + *p_j,*p_rot1 * *p_rot2, *p_Momentum);
           p_j++;
           p_j1++;
           p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj   += shelj;
        bfposj1  += shelj1;
      } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
    } // close loop over index_i
   }  // close loop on k

 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"momentum\n");
    for(k = 0; k < 3; k++) {
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.3lf",one_ints->Momentum[dim6 + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }
      }

   dim6 += 3 * nd1 * nd2;
   dim2 += 3 * nd3 * nd4;
  } // close loop on p
 }

 if (Function[4] == 1) {
  dim2 = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];
    nd3 = atoms->bfnnumb[ip];
    nd4 = atoms->bfnnumb[jp];
    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli1 = shells->type[index_i];
    sheli  = shells->type1[index_i];
    nsheli = *(shells->num_ij + shells->ord[index_i]);
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
      shelj1 = shells->type[index_j];
      shelj  = shells->type1[index_j];
      nshelj = *(shells->num_ij + shells->ord[index_j]);
        p_i1 = shells->ind_i + shells->opp[index_i];
        p_i  = shells->ind_j + shells->opp[index_i];
        p_rot1 = shells->rot + shells->opp[index_i];
        for (i = 0; i < nsheli; i++) {
          p_j1 = shells->ind_i + shells->opp[index_j];
          p_j  = shells->ind_j + shells->opp[index_j];
          p_rot2 = shells->rot + shells->opp[index_j];
          for (j = 0; j < nshelj; j++) {
            p_Overlap = Overlap + dim2 + (bfposi1 + *p_i1) * nd4 + bfposj1 + *p_j1;
            one_ints->Overlap[dim7 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_Overlap;
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj   += shelj;
        bfposj1  += shelj1;
      } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
    } // close loop over index_i

    if (job->verbosity > 1) {
      count = 0;
      fprintf(file.out,"Overlap   Pair[%2d][%2d]  [%2d]\n",ip,jp,gj);
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%10.4lf",one_ints->Overlap[dim7 + count]);
          //fprintf(file.out,"%18.12lf",one_ints->Overlap[dim7 + count]);
          //printf("%20.15e ",one_ints->Overlap[dim7 + count]);
           count++;
          }
         fprintf(file.out,"\n");
         //printf("\n");
        }
        fprintf(file.out,"\n");
        //printf("\n");
       }

    dim7 += nd1 * nd2;
    dim2 += nd3 * nd4;
  } // close loop on p
 }

 if (Function[5] == 1) {
  dim2 = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];
    nd3 = atoms->bfnnumb[ip];
    nd4 = atoms->bfnnumb[jp];
    for (k = 0; k < 6; k++) {
    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli1 = shells->type[index_i];
    sheli  = shells->type1[index_i];
    nsheli = *(shells->num_ij + shells->ord[index_i]);
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
      shelj1 = shells->type[index_j];
      shelj  = shells->type1[index_j];
      nshelj = *(shells->num_ij + shells->ord[index_j]);
        p_i1 = shells->ind_i + shells->opp[index_i];
        p_i  = shells->ind_j + shells->opp[index_i];
        p_rot1 = shells->rot + shells->opp[index_i];
        for (i = 0; i < nsheli; i++) {
          p_j1 = shells->ind_i + shells->opp[index_j];
          p_j  = shells->ind_j + shells->opp[index_j];
          p_rot2 = shells->rot + shells->opp[index_j];
          for (j = 0; j < nshelj; j++) {
         p_Grad_Grad = Grad_Grad + dim2 + k * nd3 * nd4 + (bfposi1 + *p_i1) * nd4 + bfposj1 + *p_j1;
         one_ints->Grad_Grad[dim8 + k * nd1 * nd2 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_Grad_Grad;
           p_j++;
           p_j1++;
           p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj   += shelj;
        bfposj1  += shelj1;
      } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
    } // close loop over index_i
   }  // close loop on k

 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"grad_grad\n");
    for(k = 0; k < 6; k++) {
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.2lf",one_ints->Grad_Grad[dim8 + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }
      }

   dim8 += 6 * nd1 * nd2;
   dim2 += 6 * nd3 * nd4;
  } // close loop on p
 }

 if (Function[6] == 1) {
  dim2 = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];
    nd3 = atoms->bfnnumb[ip];
    nd4 = atoms->bfnnumb[jp];
    for (k = 0; k < 3; k++) {
    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli1 = shells->type[index_i];
    sheli  = shells->type1[index_i];
    nsheli = *(shells->num_ij + shells->ord[index_i]);
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
      shelj1 = shells->type[index_j];
      shelj  = shells->type1[index_j];
      nshelj = *(shells->num_ij + shells->ord[index_j]);
        p_i1 = shells->ind_i + shells->opp[index_i];
        p_i  = shells->ind_j + shells->opp[index_i];
        p_rot1 = shells->rot + shells->opp[index_i];
        for (i = 0; i < nsheli; i++) {
          p_j1 = shells->ind_i + shells->opp[index_j];
          p_j  = shells->ind_j + shells->opp[index_j];
          p_rot2 = shells->rot + shells->opp[index_j];
          for (j = 0; j < nshelj; j++) {
         p_Dipole = Dipole + dim2 + k * nd3 * nd4 + (bfposi1 + *p_i1) * nd4 + bfposj1 + *p_j1;
         one_ints->Dipole[dim6 + k * nd1 * nd2 + (bfposi + *p_i) * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_Dipole;
           p_j++;
           p_j1++;
           p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj   += shelj;
        bfposj1  += shelj1;
      } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
    } // close loop over index_i
   }  // close loop on k


 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"dipole\n");
    for(k = 0; k < 3; k++) {
      fprintf(file.out,"xyz %3d ip %3d jp %3d\n",k,ip,jp);
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%8.4lf",one_ints->Dipole[dim6 + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }
      }

   dim6 += 3 * nd1 * nd2;
   dim2 += 3 * nd3 * nd4;
  } // close loop on p
 }

   free(ElecNuc);
   free(Kinetic);
   free(Momentum);
   free(Overlap);
   free(Grad_Grad);
   free(Dipole);

   double time2 = MPI_Wtime();

  if (job->verbosity > 1)
  fprintf(file.out, "Time for 1e integrals %10.4e %d\n", (double) time2 - time1, iii);

}
