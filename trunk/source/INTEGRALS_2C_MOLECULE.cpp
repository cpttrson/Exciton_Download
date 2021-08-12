

  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include "myconstants.h"
#include "USER_DATA.h"
#include "PARALLEL.h"
#include "MATRIX_UTIL.h"
#include "RECURSION.h"
#include "CARTESIAN_TO_SH.h"
#include "E_COEFFICIENTS.h"
#include "INCOMPLETE_GAMMA.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "MCMURCHIE_DAVIDSON.h"
#include "INTEGRALS_2C_MOLECULE.h"

using namespace std;

void integrals_molecule_ij(double *Coulomb, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, REAL_LATTICE *R, JOB_PARAM *job, FILES file)

{

int ip, jp, gi, gj, p;
int i, j, i4, j4, count;
int index_i, index_j;
int dim, dim3;
int bfposi, bfposj, bfposi1, bfposj1;
int imax, jmax, im, jm;
int gausposi, gausposj, shelposi, shelposj;
int sheli, shelj, sheli1, shelj1;
int nd1, nd2, nd3, nd4, nd12, nd34;
int begin_p[job->numtasks], end_p[job->numtasks];
double R_AB_1esqrd;
double E1_max;
double time1, time2;
VECTOR_DOUBLE R_AB_1e;

  time1 = MPI_Wtime();

  mpi_begin_end(begin_p,end_p,pair_p->nump,job->numtasks,job,file);
  //printf("process int_1e %d begin %d end %d\n", job->taskid, begin_p[job->taskid],end_p[job->taskid]);

  im = 0; jm = 0;
  dim3 = 0;
  for (p = 0; p < begin_p[job->taskid]; p++) {
  ip = pair_p->cell1[pair_p->posn[p]];
  jp = pair_p->cell2[pair_p->posn[p]];
  nd3 = atoms->bfnnumb_sh[ip];
  nd4 = atoms->bfnnumb_sh[jp];
  dim3 += nd3 * nd4;
 }

  dim = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    double coulomb_cart[nd12]; 
    for (i = 0; i < nd12; i++) coulomb_cart[i] = k_zero;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type_sh[index_i];
      sheli1 = shells->type1_sh[index_i];
      imax   = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->type_sh[index_j];
        shelj1 = shells->type1_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+im+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,\
        R_AB,R,atoms,shells,gaussians,job,file);
        two_centre_coulomb(coulomb_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd2,im,jm,E1x,E1y,E1z,&R_AB_1e,\
        atoms,shells,gaussians,job,file);
        cartesian_to_sh_ij(coulomb_cart,&Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i

 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"Coulomb\n");
      fprintf(file.out,"p %3d ip %3d jp %3d gj %3d\n",p,ip,jp,gj);
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.2lf",Coulomb[dim3 + dim + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }

    dim += nd34;

  } // close loop on p

  time2 = MPI_Wtime();

}

void nuclear_repulsion_energy(double *nuc_nuc, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  /******************************************************************************************
   * Nuclear-nuclear contribution to total energy                                           *
   ******************************************************************************************/

  int count;
  int index_R, index_S, index_G;
  int atm_no, atm_n1;
  double result1, result2, result3;
  double gamma_0_inv = pow(crystal->primitive_cell_volume , two_thirds) / eight;
  double gamma_0 = k_one / G->gamma_0_inv;
  double root_gamma_0 = sqrt(gamma_0);
  double GdotR;
  double Rsqrd, root_Rsqrd;
  double shell_sum;
  double h, z;
  VECTOR_DOUBLE s_12;

  *nuc_nuc = k_zero;

  for (atm_no = 0; atm_no < atoms->number_of_atoms_in_unit_cell; atm_no++) {
    for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {

      result1 = k_zero;
      result2 = k_zero;

      s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1;
      s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2;
      s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3;

         switch (crystal->type[0]) {

           case 'C':

        s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1;
        s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2;
        s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3;

        count = 1;
        for (index_S = 1; index_S < G->number_of_shells; index_S++) {
        for (index_G = 0; index_G < G->num[index_S]; index_G++) {
          GdotR = double_vec_dot(&G->vec_b2[count], &s_12);
          count++;
          if (crystal->type[0] == 'C')
            result1 += G->EXPFAC[index_S] * cos(GdotR);
       }
      }
        count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
        s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1 + R->vec_ai[count].comp1;
        s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2 + R->vec_ai[count].comp2;
        s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3 + R->vec_ai[count].comp3;
        Rsqrd = double_vec_dot(&s_12, &s_12);
        root_Rsqrd = sqrt(Rsqrd);
        if (atm_n1 == atm_no && Rsqrd < 0.00001) {
          result2 -= two * root_gamma_0 / rtpi;
        } else {
          result2 += erfc(root_gamma_0 * root_Rsqrd) / root_Rsqrd;
          shell_sum += erfc(root_gamma_0 * root_Rsqrd) / root_Rsqrd;
        } // close else { result2
            count++;
           } // end loop on index_R
              } // end loop over index_S

        result3 = -pi * G->gamma_0_inv / crystal->primitive_cell_volume;

       break;

       case 'S':

        result1 = k_zero;
        result2 = k_zero;

        s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1;
        s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2;
        s_12.comp3 = k_zero;
        z = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3;

        count = 1;
        for (index_S = 1; index_S < G->number_of_shells; index_S++) {
          h = sqrt(G->sqr[index_S]);
          for (index_G = 0; index_G < G->num[index_S]; index_G++) {
            GdotR = double_vec_dot(&G->vec_b2[count], &s_12);
            result1 += two * pi / crystal->primitive_cell_volume * cos(GdotR) * (exp(h * z) * \
            erfc(z * root_gamma_0 + h / two / root_gamma_0) + \
            exp(-h * z) * erfc(-z * root_gamma_0 + h / two / root_gamma_0)) / h;
            count++;
          }
         }

        count = 0;
        for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
          shell_sum = k_zero;
          for (index_R = 0; index_R < R->num[index_S]; index_R++) {
            s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1 + R->vec_ai[count].comp1;
            s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2 + R->vec_ai[count].comp2;
            s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3 + R->vec_ai[count].comp3;
            Rsqrd = double_vec_dot(&s_12, &s_12);
            root_Rsqrd = sqrt(Rsqrd);
            if (atm_n1 == atm_no && Rsqrd < 0.00001) {
            result2 -= two * root_gamma_0 / rtpi;
           } else {
            result2 += erfc(root_gamma_0 * root_Rsqrd) / root_Rsqrd;
            shell_sum += erfc(root_gamma_0 * root_Rsqrd) / root_Rsqrd;
           } // close else { result2
            count++;
           } // end loop on index_R
          } // end loop over index_S

           result3 = -two * pi / crystal->primitive_cell_volume * \
           (z * erf(root_gamma_0 * z) + exp(-gamma_0 * z * z) / rtpi / root_gamma_0);

       break;

       case 'P':

       result1 =  k_zero;
       result2 =  k_zero;
       result3 =  k_zero;

       if (job->taskid == 0)
       fprintf(file.out,"Nuclear-nuclear repulsion not coded for 1-D systems\n");

       break;

       case 'M':

         result1 =  k_zero;
         result2 =  k_zero;
         result3 =  k_zero;

         s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1;
         s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2;
         s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3;

         if (atm_no != atm_n1) {
         Rsqrd = double_vec_dot(&s_12, &s_12);
         result2 =  k_one / sqrt(Rsqrd);
        }

       break;

      } // close switch

      *nuc_nuc += double(atoms->atomic_number[atm_no] * atoms->atomic_number[atm_n1]) * (result1 + result2 + result3) / two;

    }
  }

  job->nuc_nuc = *nuc_nuc;

}

void fock_element_1e1(INT_1E *one_ints, int dim, PAIR_TRAN *pair_p, int num_p, int Function[8], REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int ip, jp, gi, gj, p;
  int i, j, i4, j4, k, l, mm, nn, atm0, count;
  int index_i, index_j;
  int dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, dim9, atm1, atm2;
  int bfposi, bfposj, bfposi1, bfposj1;
  int imax, jmax, im, jm;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4, nd12, nd34;
  int nsheli, nshelj;
  int begin_p[job->numtasks], end_p[job->numtasks];
  double R_AB_1esqrd;
  double E1_max;
  double time1, time2;
  VECTOR_DOUBLE R_AB_1e;

  // Function[0] Fock Operator
  // Function[1] Kinetic Energy Operator
  // Function[2] Electron-Nuclear Potential Energy = 1 (finite systems) = 2 (extended systems) = 3 (atoms)
  // Function[3] Momentum Operator Matrix Elements
  // Function[4] Overlap
  // Function[5] 
  // Function[6] Dipole Operator Matrix Elements
  // Function[7] Coulomb Operator Matrix Elements

  time1 = MPI_Wtime();

  //mpi_begin_end(begin_p,end_p,num_p,job->numtasks,job,file);
  //printf("process int_1e %d begin %d end %d\n", job->taskid, begin_p[job->taskid],end_p[job->taskid]);

  //dim3 = 0;
  //for (p = 0; p < begin_p[job->taskid]; p++) {
  //ip = pair_p->cell1[pair_p->posn[p]];
  //jp = pair_p->cell2[pair_p->posn[p]];
  //nd3 = atoms->bfnnumb_sh[ip];
  //nd4 = atoms->bfnnumb_sh[jp];
  //dim3 += nd3 * nd4;
 //}
    //dim = 0;
    //dim4 = 0;
  
  //for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {


  int *counter, winRank, myCounter, increment;
  const int izero = 0;
  MPI_Win win;
  winRank = 0;
  CreateCounter(job->taskid, winRank, &counter, izero, &win, MPI_COMM_WORLD);
  myCounter = 0;

  if (job->taskid > 0 || job->numtasks == 1)
  while (myCounter < pair_p->nump) {

    increment = (pair_p->nump - myCounter) / job->numtasks / 2;
    if (increment <  1) increment =  1;
    if (increment > 64) increment = 64;
    myCounter = GetCounter(winRank, increment , &win);
    //printf("proc %d myCounter %d\n",job->taskid,myCounter);
    int begin_p1 = (myCounter             < pair_p->nump) ? myCounter : pair_p->nump;
    int end_p1   = (myCounter + increment < pair_p->nump) ? myCounter + increment : pair_p->nump;
    //printf("proc %3d Fock begin %6d end %6d out of %6d\n",job->taskid,begin_p1,end_p1,pair_p->nump);
    dim = 0;
    dim3 = 0;
    dim4 = 0;
    for (p = 0; p < begin_p1; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    dim3 += nd3 * nd4;
   }

  for (p = begin_p1; p < end_p1; p++) {

    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    double kinetic_cart[nd12]; 
    double elecnuc_cart[nd12]; 
    double momentum_cart[3 * nd12]; 
    double overlap_cart[nd12]; 
    double dipole_cart[3 * nd12]; 
    double coulomb_cart[nd12]; 
    // order matters for these loops to set im and jm
    im = 0;
    jm = 0;
    if (Function[0] || Function[2]) {
    for (i = 0; i < nd12; i++) elecnuc_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[4]) {
    for (i = 0; i < nd12; i++) overlap_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[7]) {
    for (i = 0; i < nd12; i++) coulomb_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[3]) {
    for (i = 0; i < 3 * nd12; i++) momentum_cart[i] = k_zero;
    im = 1; jm = 1;
   }
    if (Function[6]) {
    for (i = 0; i < 3 * nd12; i++) dipole_cart[i] = k_zero;
    im = 1; jm = 1;
   }
    if (Function[0] || Function[1]) {
    for (i = 0; i < nd12; i++) kinetic_cart[i] = k_zero;
    im = 2; jm = 2;
   }
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type_sh[index_i];
      sheli1 = shells->type1_sh[index_i];
      imax   = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->type_sh[index_j];
        shelj1 = shells->type1_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+im+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,\
        R_AB,R,atoms,shells,gaussians,job,file);
    if (Function[0] || Function[1]) {
        Kinetic(kinetic_cart, index_i, index_j, bfposi1, bfposj1, gausposi, gausposj, nd2, im, jm, E1x, E1y, E1z, sab, \
        shells, gaussians, job, file);
        cartesian_to_sh_ij(kinetic_cart,&one_ints->Kinetic[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[0] || Function[2]) {
        ElecNuc(elecnuc_cart, index_i, index_j, bfposi1, bfposj1, nd2,  im, jm, E1x, E1y, E1z, R_AB, pab_inv, sab, R, G, \
        crystal, atoms,shells,job,file);
        cartesian_to_sh_ij(elecnuc_cart,&one_ints->ElecNuc[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[3]) {
        Momentum(momentum_cart, index_i, index_j, bfposi1, bfposj1, gausposi, gausposj, nd1, nd2, im, jm, E1x, E1y, E1z, sab, \
        shells, gaussians, job,file);
        cartesian_to_sh_ij_vector(momentum_cart,&one_ints->Momentum[3 * (dim3 + dim)],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
       }
    if (Function[4]) {
        Overlap(overlap_cart, index_i, index_j, bfposi1, bfposj1, nd2, im, jm, E1x, E1y, E1z, sab, shells, job, file);
        cartesian_to_sh_ij(overlap_cart,&one_ints->Overlap[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[6]) {
        Dipole(dipole_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd1,nd2,ip,im,jm,E1x,E1y,E1z,sab,R_AB,atoms,shells,\
        gaussians,job,file);
        cartesian_to_sh_ij_vector(dipole_cart,&one_ints->Dipole[3 * (dim3 + dim)],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
       }
    if (Function[7]) {
        two_centre_coulomb(coulomb_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd2,im,jm,E1x,E1y,E1z,&R_AB_1e,\
        atoms,shells,gaussians,job,file);
        cartesian_to_sh_ij(coulomb_cart,&one_ints->Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
       }
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i

 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"ElecNuc\n");
      fprintf(file.out,"XYZ %3d ip %3d jp %3d gj %3d\n",p,ip,jp,gj);
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.2lf",one_ints->Overlap[dim3 + dim + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }

    dim += nd34;
    dim4 += 3 * nd34;

  } // close loop on p

    if (Function[0]) {
      for (i = dim3; i < dim3 + dim; i++) {
        one_ints->Fock[i] = one_ints->Kinetic[i] + one_ints->ElecNuc[i];
       }
      }

  } // close while (
        
   DestroyCounter(job->taskid, winRank, &win, counter);

   time2 = MPI_Wtime();

}

void fock_element_1e2(INT_1E *one_ints, PAIR_TRAN *pair_p, int Function[8], REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int ip, jp, gi, gj, p;
int i, j, i4, j4, count;
int index_i, index_j;
int dim, dim3, dim4, dimall, dimall3, dimg;
int bfposi, bfposj, bfposi1, bfposj1;
int imax, jmax, im, jm;
int gausposi, gausposj, shelposi, shelposj;
int sheli, shelj, sheli1, shelj1;
int nd1, nd2, nd3, nd4, nd12, nd34;
int nsheli, nshelj;
int begin_p[job->numtasks], end_p[job->numtasks];
double R_AB_1esqrd;
double E1_max;
double time1, time2;
VECTOR_DOUBLE R_AB_1e;
INT_1E one_ints_buffer;

  // Function[0] Fock Operator
  // Function[1] Kinetic Energy Operator
  // Function[2] Electron-Nuclear Potential Energy = 1 (finite systems) = 2 (extended systems) = 3 (atoms)
  // Function[3] Momentum Operator Matrix Elements
  // Function[4] Overlap
  // Function[5] 
  // Function[6] Dipole Operator Matrix Elements
  // Function[7] Coulomb Operator Matrix Elements

  time1 = MPI_Wtime();

  mpi_begin_end(begin_p,end_p,pair_p->nump,job->numtasks,job,file);
  //printf("process int_1e %d begin %d end %d\n", job->taskid, begin_p[job->taskid],end_p[job->taskid]);
  array_dimensions(&dimall, &dimg, pair_p, atoms, job, file); 
  allocate_INT_1E(&one_ints_buffer, dimall, Function, job, file);
  dimall3 = 3 * dimall;

  dim3 = 0;
  for (p = 0; p < begin_p[job->taskid]; p++) {
  ip = pair_p->cell1[pair_p->posn[p]];
  jp = pair_p->cell2[pair_p->posn[p]];
  nd3 = atoms->bfnnumb_sh[ip];
  nd4 = atoms->bfnnumb_sh[jp];
  dim3 += nd3 * nd4;
 }

  dim = 0;
  dim4 = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    double kinetic_cart[nd12]; 
    double elecnuc_cart[nd12]; 
    double momentum_cart[3 * nd12]; 
    double overlap_cart[nd12]; 
    double dipole_cart[3 * nd12]; 
    double coulomb_cart[nd12]; 
    // order matters for these loops to set im and jm
    im = 0;
    jm = 0;
    if (Function[0] || Function[2]) {
    for (i = 0; i < nd12; i++) elecnuc_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[4]) {
    for (i = 0; i < nd12; i++) overlap_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[7]) {
    for (i = 0; i < nd12; i++) coulomb_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[3]) {
    for (i = 0; i < 3 * nd12; i++) momentum_cart[i] = k_zero;
    im = 1; jm = 1;
   }
    if (Function[6]) {
    for (i = 0; i < 3 * nd12; i++) dipole_cart[i] = k_zero;
    im = 1; jm = 1;
   }
    if (Function[0] || Function[1]) {
    for (i = 0; i < nd12; i++) kinetic_cart[i] = k_zero;
    im = 2; jm = 2;
   }
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type_sh[index_i];
      sheli1 = shells->type1_sh[index_i];
      imax   = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->type_sh[index_j];
        shelj1 = shells->type1_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+im+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,\
        R_AB,R,atoms,shells,gaussians,job,file);
    if (Function[0] || Function[1]) {
        Kinetic(kinetic_cart, index_i, index_j, bfposi1, bfposj1, gausposi, gausposj, nd2, im, jm, E1x, E1y, E1z, sab, \
        shells, gaussians, job, file);
        cartesian_to_sh_ij(kinetic_cart,&one_ints_buffer.Kinetic[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[0] || Function[2]) {
        ElecNuc(elecnuc_cart, index_i, index_j, bfposi1, bfposj1, nd2,  im, jm, E1x, E1y, E1z, R_AB, pab_inv, sab, R, G, \
        crystal, atoms,shells,job,file);
        cartesian_to_sh_ij(elecnuc_cart,&one_ints_buffer.ElecNuc[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[3]) {
        Momentum(momentum_cart, index_i, index_j, bfposi1, bfposj1, gausposi, gausposj, nd1, nd2, im, jm, E1x, E1y, E1z, sab, \
        shells, gaussians, job,file);
        cartesian_to_sh_ij_vector(momentum_cart,&one_ints_buffer.Momentum[3 * dim3 + dim4],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
       }
    if (Function[4]) {
        Overlap(overlap_cart, index_i, index_j, bfposi1, bfposj1, nd2, im, jm, E1x, E1y, E1z, sab, shells, job, file);
        cartesian_to_sh_ij(overlap_cart,&one_ints_buffer.Overlap[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[6]) {
        Dipole(dipole_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd1,nd2,ip,im,jm,E1x,E1y,E1z,sab,R_AB,atoms,shells,\
        gaussians,job,file);
        cartesian_to_sh_ij_vector(dipole_cart,&one_ints_buffer.Dipole[3 * dim3 + dim4],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
       }
    if (Function[7]) {
        two_centre_coulomb(coulomb_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd2,im,jm,E1x,E1y,E1z,&R_AB_1e,\
        atoms,shells,gaussians,job,file);
        cartesian_to_sh_ij(coulomb_cart,&one_ints_buffer.Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
       }
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i

 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"Coulomb\n");
      fprintf(file.out,"p %3d ip %3d jp %3d gj %3d\n",p,ip,jp,gj);
      for(int k=0;k<3;k++) {
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.2lf %6.2f",one_ints_buffer.Dipole[3 * dim3 + dim4 + count],dipole_cart[count]);
          //fprintf(file.out,"%6.2lf",one_ints_buffer.Coulomb[dim3 + dim + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }
        fprintf(file.out,"\n");
       }

    dim += nd34;
    dim4 += 3 * nd34;

  } // close loop on p

  // Function[0] Fock Operator
  // Function[1] Kinetic Energy Operator
  // Function[2] Electron-Nuclear Potential Energy = 1 (finite systems) = 2 (extended systems) = 3 (atoms)
  // Function[3] Momentum Operator Matrix Elements
  // Function[4] Overlap
  // Function[5] 
  // Function[6] Dipole Operator Matrix Elements
  // Function[7] Coulomb Operator Matrix Elements

  if (Function[1]) MPI_Allreduce(&one_ints_buffer.Kinetic[0],  &one_ints->Kinetic[0],  dimall,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[2]) MPI_Allreduce(&one_ints_buffer.ElecNuc[0],  &one_ints->ElecNuc[0],  dimall,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[3]) MPI_Allreduce(&one_ints_buffer.Momentum[0], &one_ints->Momentum[0], dimall3, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[4]) MPI_Allreduce(&one_ints_buffer.Overlap[0],  &one_ints->Overlap[0],  dimall,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[6]) MPI_Allreduce(&one_ints_buffer.Dipole[0],   &one_ints->Dipole[0],   dimall3,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[7]) MPI_Allreduce(&one_ints_buffer.Coulomb[0],  &one_ints->Coulomb[0],  dimall  ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer, Function, job, file);

  if (Function[0]) {
    for (i = 0; i < dimall; i++) {
      one_ints->Fock[i] = one_ints->Kinetic[i] + one_ints->ElecNuc[i];
     }
    }

   time2 = MPI_Wtime();

}

void Kinetic(double *Kinetic, int index_i, int index_j, int bfposi, int bfposj, int gausposi, int gausposj, int nd2,  int im, int jm, double *E1x, double *E1y, double *E1z, double *sab, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, i1, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int sheli, shelj;
  int imax, jmax;
  int off1, off2, off3;
  int count;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      i1 = 0;
      count = (bfposi + i) * nd2 + bfposj + j;
      for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
      for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
      if (n4 > 1)
      Kinetic[count] -= (double)(n4 * (n4 - 1)) * E1x[n1 * off1 + (n4 - 2) * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * \
                         E1z[n3 * off1 + n6 * off2 + i1] * sab[i1] / two;
      if (n5 > 1)
      Kinetic[count] -= (double)(n5 * (n5 - 1)) * E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + (n5 - 2) * off2 + i1] * \
                         E1z[n3 * off1 + n6 * off2 + i1] * sab[i1] / two;
      if (n6 > 1)
      Kinetic[count] -= (double)(n6 * (n6 - 1)) * E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * \
                         E1z[n3 * off1 + (n6 - 2) * off2 + i1] * sab[i1] / two;
      Kinetic[count] -= (-two * gaussians->expo_sh[gausposj + j4] * (double)(2 * (n4 + n5 + n6) + 3) * \
                    E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] + \
                    four * gaussians->expo_sh[gausposj + j4] * gaussians->expo_sh[gausposj + j4] * \
                   (E1x[n1 * off1 + (n4 + 2) * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] + \
                    E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + (n5 + 2) * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] + \
                    E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + (n6 + 2) * off2 + i1])) * sab[i1] / two;
       i1++;
      } // close loop over j4
     } // close loop over i4
    } // close loop over i
   } // close loop over i

}

void ElecNuc(double *ElecNuc, int index_i, int index_j, int bfposi, int bfposj, int nd2,  int im, int jm, double *E1x, double *E1y, double *E1z, VECTOR_DOUBLE *R_AB, double *pab_inv, double *sab, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, CRYSTAL *crystal, ATOM *atoms, SHELL *shells, JOB_PARAM *job, FILES file)

{

  int i, j, k, i4, nn, t, u, v;
  int n1, n2, n3, n4, n5, n6;
  int tmax, umax, vmax;
  int sheli, shelj;
  int off1, off2, off3;
  int imax, jmax;
  int count;
  int atm0;
  int count1, index_G, index_R, index_S;
  double Rsqrd, fn[55], gn[55];
  double gamma_0, root_gamma_0;
  double r1, r2, r3, r4, GdotR;
  VECTOR_DOUBLE r_12, s_12, t_12, Rvec_tmp;

  double em[1][55];
  double f[13][13][13][13];
  double temp, fgtuv_max, fgtuv_temp;
  double *fgtuv;
  int mm, mm0, mm1, mm2, mm3, mm4;

  double Az;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];

  mm  = imax + jmax;
  mm3 = off3;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  mm0 = (mm + 1) * mm1;

  fgtuv = (double *) malloc(mm0 * sizeof(double));
  if (fgtuv == NULL) {
  if (job->taskid == 0)
  fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
  MPI_Finalize();
  exit(1);
 }

int n;
double r5, Rzsqrd, sum5[mm + 1], sum2[mm + 1], x;
double fac1, fac2, expfac, x1, x2, sign;
double D_erf[16 + 1], D_exp[16 + 1], derivative_1[16 + 1], derivative_2[16 + 1];

  switch (crystal->type[0]) {

  case 'M':

  for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
    for (k = 0; k < mm0; k++) fgtuv[k] = k_zero;
    for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
      s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
      s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
      s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
      Rsqrd = double_vec_dot(&s_12, &s_12);
      f000m(&em[0][0], Rsqrd / pab_inv[i4], k_one / pab_inv[i4], mm);
      non_recursive_ftuvn(mm, 0, f, em, &s_12);
      temp = sab[i4] * atoms->atomic_number[atm0];
      for (t = 0; t <= mm; t++) {
        for (u = 0; u <= mm; u++) {
          for (v = 0; v <= mm; v++) {
            if (t + u + v > mm) break;
            fgtuv_temp = f[t][u][v][0] * temp;
            fgtuv[t * mm1 + u * mm2 + v * mm3  + i4] += fgtuv_temp;
            fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
            fgtuv_max = k_one;
           }
          }
         }
        }

  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (t = 0; t <= tmax; t++) {
        for (u = 0; u <= umax; u++) {
          for (v = 0; v <= vmax; v++) {
            for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
               ElecNuc[count] -= E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                                 E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * fgtuv[t * mm1 + u * mm2 + v * mm3  + i4];
              }
             }
            }
           } // end t u v loop
          } // close loop over j
         } // close loop over i
        } // end loop over atm0

        free(fgtuv);

  break;

  case 'C':

  gamma_0 = k_one / G->gamma_0_inv;
  root_gamma_0 = sqrt(gamma_0);

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      r3 = k_zero;
      for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
          s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
          s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
          Rsqrd = double_vec_dot(&s_12, &s_12);
          r4 = atoms->atomic_number[atm0] * sab[i4];
          for (t = 0; t <= tmax; t++) {
            for (u = 0; u <= umax; u++) {
              for (v = 0; v <= vmax; v++) {
                nn = t + u + v;
                count1 = 1;
                r1 = r4 * E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                E1z[n3 * off1 + n6 * off2 + v * off3+i4];
                for (index_S = 1; index_S < G->number_of_shells; index_S++) {
                  r2 = r1 * G->EXPFAC[index_S]; 
                  for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                    GdotR = double_vec_dot(&G->vec_b2[count1], &s_12);
                    r3 += r2 * G->x[count1 * 9 + t] * G->y[count1 * 9 + u] * G->z[count1 * 9 + v] * cosfactor(nn, GdotR);
                    count1++;
                   } // end loop over index_G
                  } // end loop over index_S
                 }
                }
               } // end t u v loop
              } // end loop over i4
             } // end loop over atm0
            ElecNuc[count] -= r3;
           } // close loop over j
          } // close loop over i


  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          r4 = sab[i4] * (double)-atoms->atomic_number[atm0];
          r1 = r4 * -pi * (G->gamma_0_inv - pab_inv[i4]) / crystal->primitive_cell_volume;
          s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
          s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
          s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
          map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
          count1 = 0;
          for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
            double shell_sum = k_zero;
            for (index_R = 0; index_R < R->num[index_S]; index_R++) {
              r_12.comp1 = t_12.comp1 + R->vec_ai[count1].comp1;
              r_12.comp2 = t_12.comp2 + R->vec_ai[count1].comp2;
              r_12.comp3 = t_12.comp3 + R->vec_ai[count1].comp3;
              Rsqrd = double_vec_dot(&r_12, &r_12);
              for (t = 0; t <= tmax; t++) {
                for (u = 0; u <= umax; u++) {
                  for (v = 0; v <= vmax; v++) {
                    nn = t + u + v;
                    f000m(fn, Rsqrd / pab_inv[i4], k_one / pab_inv[i4], nn);
                    f000m(gn, gamma_0 * Rsqrd, gamma_0, nn);
		    shell_sum += E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                    E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * r4;
                   }
                  }
                 } // end t u v loop
                count1++;
               } // end index_R loop
               ElecNuc[count] += shell_sum;
               if (fabs(shell_sum) < 1e-13 && sab[i4] < 1e-13) break;
              } // end index_S loop
          ElecNuc[count] += E1x[n1 * off1 + n4 * off2 + i4] * E1y[n2 * off1 + n5 * off2 + i4] * E1z[n3 * off1 + n6 * off2 + i4] * r1;
            } // end loop over i4
           } // end loop over atm0
          } // close loop over j
         } // close loop over i

  break;

  case 'S':

  gamma_0 = k_one / G->gamma_0_inv;
  root_gamma_0 = sqrt(gamma_0);
  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      r3 = k_zero;
      for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
          s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
          s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
          Rzsqrd = s_12.comp3 * s_12.comp3;
          Rsqrd = double_vec_dot(&s_12, &s_12);
          r4 = -atoms->atomic_number[atm0] * sab[i4];
          //r5 = r4 * -two * pi / crystal->primitive_cell_volume * \
         (exp(-gamma_0 * Rzsqrd) / root_gamma_0 / rtpi + s_12.comp3 * erf(root_gamma_0 * s_12.comp3));
          x = root_gamma_0 * s_12.comp3;
          erf_exp_derivative(D_erf, D_exp, mm, root_gamma_0, x, G);
          sum5[0] = D_exp[0] / rtpi / root_gamma_0 + s_12.comp3 * D_erf[0];
          for (n = 1; n <= mm; n++) {
          sum5[n] = D_exp[n] / rtpi / root_gamma_0 + s_12.comp3 * D_erf[n] + (double) n * D_erf[n - 1];
         }
          r5 = k_zero;
          for (v = 0; v <= vmax; v++) {
            r5 += -two * pi / crystal->primitive_cell_volume * sum5[v] * r4 * \
            E1x[n1 * off1 + n4 * off2 + i4] * E1y[n2 * off1 + n5 * off2 + i4] * E1z[n3 * off1 + n6 * off2 + v * off3 + i4];
           }
	    ElecNuc[count] += r5;
          for (t = 0; t <= tmax; t++) {
            for (u = 0; u <= umax; u++) {
              for (v = 0; v <= vmax; v++) {
                nn = t + u;
                count1 = 1;
              r1 = r4 * E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
              E1z[n3 * off1 + n6 * off2 + v * off3 + i4];
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
                      sum2[n] += G->A->a[n][k] * G->shell_mag[9 * index_S + k] * (expfac * derivative_1[n - k] + sign / expfac * \
                      derivative_2[n - k]) / G->shell_mag[9 * index_S + 1];
                      sign *= -k_one;
                     }
                    }
                  //fac2 = fac1 * G->EXPFAC[index_S];
                  //double shell_sum = k_zero;
                  r2 = r1 *  two * pi / crystal->primitive_cell_volume * sum2[v];
                  for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                    GdotR = double_vec_dot(&G->vec_b2[count1], &s_12);
                    r3 += r2 * G->x[count1 * 9 + t] * G->y[count1 * 9 + u] * cosfactor(nn, GdotR);
                    count1++;
                   } // end loop over index_G
                  //r3 += shell_sum;
                  //if (fabs(shell_sum) < 1e-24 && sab[i4] < 1e-24) continue;
                  } // end loop over index_S
                 }
                }
               } // end t u v loop
              } // end loop over i4
             } // end loop over atm0
            ElecNuc[count] += r3;
           } // close loop over j
          } // close loop over i


  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          r4 = sab[i4] * (double)-atoms->atomic_number[atm0];
          s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
          s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
          s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
          map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
          count1 = 0;
          for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
            double shell_sum = k_zero;
            for (index_R = 0; index_R < R->num[index_S]; index_R++) {
              r_12.comp1 = t_12.comp1 + R->vec_ai[count1].comp1;
              r_12.comp2 = t_12.comp2 + R->vec_ai[count1].comp2;
              r_12.comp3 = t_12.comp3 + R->vec_ai[count1].comp3;
              Rsqrd = double_vec_dot(&r_12, &r_12);
              for (t = 0; t <= tmax; t++) {
                for (u = 0; u <= umax; u++) {
                  for (v = 0; v <= vmax; v++) {
                    nn = t + u + v;
                    f000m(fn, Rsqrd / pab_inv[i4], k_one / pab_inv[i4], nn);
                    f000m(gn, gamma_0 * Rsqrd, gamma_0, nn);
		    shell_sum += E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                    E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * r4;
                    //fprintf(file.out,"shellsum %3d %10.3e\n",index_S,r4*(ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)));
		    //ElecNuc[count] += E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                    E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * r4;
                   }
                  }
                 } // end t u v loop
                count1++;
               } // end index_R loop
               ElecNuc[count] += shell_sum;
               //fprintf(file.out," %20.14lf\n", shell_sum);
               //if (fabs(shell_sum) < 1e-13 && sab[i4] < 1e-13) break;
              } // end index_S loop
	     //ElecNuc[count] += E1x[n1 * off1 + n4 * off2 + i4] * E1y[n2 * off1 + n5 * off2 + i4] * E1z[n3 * off1 + n6 * off2 + i4] * r1;
             //fprintf(file.out,"Recip real other all %10.3e  %10.3e %10.3e %16.8e\n",r3,ElecNuc[count]-r3,r1,ElecNuc[count]);
            } // end loop over i4
           } // end loop over atm0
          } // close loop over j
         } // close loop over i

  break;

  case 'P':

  if (job->taskid == 0)

  fprintf(file.out,"Electron-nuclear attraction for 1-D not coded\n");

  break;

      } // close switch

}

void Momentum(double *Momentum, int index_i, int index_j, int bfposi, int bfposj, int gausposi, int gausposj, int nd1, int nd2, int im, int jm, double *E1x, double *E1y, double *E1z, double *sab, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, i1, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int sheli, shelj;
  int off1, off2, off3;
  int imax, jmax;
  int countx, county, countz, nd12;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  nd12 = nd1 * nd2;
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      countx = (bfposi + i) * nd2 + bfposj + j;
      county = nd12 + countx;
      countz = nd12 + county;
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      i1 = 0;
      for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
      for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
        if (n4 > 0)
        Momentum[countx] += double(n4) * E1x[n1 * off1 + (n4 - 1) * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        Momentum[countx] -= two * gaussians->expo_sh[gausposj + j4] * \
                            E1x[n1 * off1 + (n4 + 1) * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        if (n5 > 0)
        Momentum[county] += double(n5) * E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + (n5 - 1) * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        Momentum[county] -= two * gaussians->expo_sh[gausposj + j4] * \
                            E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + (n5 + 1) * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        if (n6 > 0)
        Momentum[countz] += double(n6) * E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + (n6 - 1) * off2 + i1] * sab[i1];
        Momentum[countz] -= two * gaussians->expo_sh[gausposj + j4] * \
                            E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + (n6 + 1) * off2 + i1] * sab[i1];
        i1++;
       } // close loop over j4
      } // close loop over i4
     } // close loop over i
    } // close loop over i

}

void Overlap(double *Overlap, int index_i, int index_j, int bfposi, int bfposj, int nd2,  int im, int jm, double *E1x, double *E1y, double *E1z, double *sab, SHELL *shells, JOB_PARAM *job, FILES file)

{

  int i, j, i1, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int tmax, umax, vmax;
  int sheli, shelj;
  int off1, off2, off3;
  int imax, jmax;
  int count;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
      Overlap[count] += E1x[n1 * off1 + n4 * off2 + i4] * E1y[n2 * off1 + n5 * off2 + i4] * E1z[n3 * off1 + n6 * off2 + i4] * sab[i4];
     } // close loop over i4
    } // close loop over j
   } // close loop over i

}

void Dipole(double *Dipole, int index_i, int index_j, int bfposi, int bfposj, int gausposi, int gausposj, int nd1, int nd2, int ip, int im, int jm, double *E1x, double *E1y, double *E1z, double *sab, VECTOR_DOUBLE *R_AB, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, i1, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int sheli, shelj;
  int off1, off2, off3;
  int imax, jmax;
  int countx, county, countz, nd12;
  double PAx, PAy, PAz;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  nd12 = nd1 * nd2;
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      countx = (bfposi + i) * nd2 + bfposj + j;
      county = nd12 + countx;
      countz = nd12 + county;
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      i1 = 0;
      for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
      for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
        Dipole[countx] += (E1x[n1 * off1 + n4 * off2 + off3 + i1] + R_AB[i1].comp1 * E1x[n1 * off1 + n4 * off2 + i1]) * \
                           E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        Dipole[county] += (E1y[n2 * off1 + n5 * off2 + off3 + i1] + R_AB[i1].comp2 * E1y[n2 * off1 + n5 * off2 + i1]) * \
                           E1x[n1 * off1 + n4 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        Dipole[countz] += (E1z[n3 * off1 + n6 * off2 + off3 + i1] + R_AB[i1].comp3 * E1z[n3 * off1 + n6 * off2 + i1]) * \
                           E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * sab[i1];
        i1++;
       } // close loop over j4
      } // close loop over i4
     } // close loop over i
    } // close loop over i

}

void two_centre_coulomb(double *Coulomb, int index_i, int index_j, int bfposi1, int bfposj1, int gausposi, int gausposj, int nd2,  int im, int jm, double *E1x, double *E1y, double *E1z, VECTOR_DOUBLE *R_AB_1e, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

int t, u, v, i4, j4;
int imax, jmax;
int count;
int mm, mm0, mm1, mm2, mm3;
double Rsqrd;
double em[1][55];
double fgtuv_max, fgtuv_temp;
double *fgtuv;

  double prefac[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
  double ab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
  count = 0;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
  for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
  prefac[count] = gaussians->c_sh[gausposi + i4] * gaussians->c_sh[gausposj + j4] * \
  pi32 * pi32 / gaussians->expo_sh[gausposi + i4] / 
  gaussians->expo_sh[gausposj + j4] / sqrt(gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4]); 
  ab[count] = k_one / gaussians->expo_sh[gausposi + i4] + k_one /  gaussians->expo_sh[gausposj + j4];
  count++; 
 }
}

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  mm  = imax + jmax;
  mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  mm0 = (mm + 1) * mm1;

  fgtuv = (double *) calloc(mm0, sizeof(double));
  if (fgtuv == NULL) {
  if (job->taskid == 0)
  fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
  MPI_Finalize();
  exit(1);
 }

  Rsqrd = double_vec_dot(R_AB_1e, R_AB_1e);
  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
  f000m(&em[0][0], Rsqrd / ab[i4], k_one / ab[i4], mm);
  for (t = 0; t <= mm; t++) {
    for (u = 0; u <= mm; u++) {
      for (v = 0; v <= mm; v++) {
        if (t + u + v > mm) break;
          fgtuv_temp = ftuvn(t,u,v,0,&em[0][0],R_AB_1e[0]) * prefac[i4];
          fgtuv[t * mm1 + u * mm2 + v * mm3  + i4] += fgtuv_temp;
          fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
          fgtuv_max = k_one;
         }
        }
       }
      }

  double C1_max, C2_max;
  double C1x[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
  double C1y[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
  double C1z[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
  double C2x[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
  double C2y[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
  double C2z[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
  E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
  E_coefficients_1c(index_j,gausposj,C2x,C2y,C2z,&C2_max,shells,gaussians,job,file);
  mcmurchie_davidson_ij(Coulomb,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,fgtuv,shells,job,file);
  free(fgtuv);

}
