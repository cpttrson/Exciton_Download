/*
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include "mycomplex.h"
#include "conversion_factors.h"
#include "mylogical.h"
#include "SYMMETRY.h"
#include "PAIRS_QUADS1.h"
#include "TOOLS.h"
*/
#include <cstring>
#include <mpi.h>
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "PAIRS_QUADS.h"
#include "LIMITS.h"
#include "INCOMPLETE_GAMMA.h" 
#include "PARALLEL.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "RECURSION.h"
#include "MCMURCHIE_DAVIDSON.h"
#include "E_COEFFICIENTS.h"
#include "CARTESIAN_TO_SH.h"
#include "INTEGRALS_4C_MOLECULE.h"

using namespace std;

void integrals_molecule_ijkl(INTEGRAL_LIST *integral_list, int ip, int jp, int kp, int lp, int *start_index, int *count_index, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // compute 4-centre integrals for a quad of atoms using pair_p and pair_q as indices
  // starting point is pair_p->posn[*p1] and pair_q->posn[*q1]

  int index_i, index_j, index_k, index_l, i4, j4, k4, l4;
  int bfposi, bfposj, bfposk, bfposl;
  int bfposi1, bfposj1, bfposk1, bfposl1;
  int gausposi, gausposj, gausposk, gausposl;
  int count;
  int i, j, k, l, m, n, p, pm, q;
  int mm, mm0, mm1, mm2, mm3, mm4;
  int gi, gj, gk, gl;
  int imax, jmax, kmax, lmax;
  int shelposi, shelposj, shelposk, shelposl;
  int sheli, shelj, shelk, shell;
  int sheli1, shelj1, shelk1, shell1;
  int offset1, offset2;
  int dim1 = atoms->number_of_atoms_in_unit_cell;
  int dim2 = dim1 * dim1;
  double p_inv_ex, fac1;
  double R_AB_1esqrd, R_CD_1esqrd, Rsqrd;
  double C1_max, C2_max, C_max;
  double fgtuv_max, fgtuv_temp;
  double time1, time2;
  double time3, time4,time5,time6,time7,time8;
  double em[1][55];
  double f[13][13][13][13];
  double *F_cart, *p_F_cart, *F_sh, *p_F_sh, *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, R_CD_1e, s_12;

int shell_index;

  time1 = MPI_Wtime();
  time3 = k_zero;;
  time5 = k_zero;;
  time7 = k_zero;;

    count = 0;

    gi = 0;
    gj = 0;
    gk = 0;
    gl = 0;


    F_cart = (double *) malloc(50625 * sizeof(double)); // 50625 = 15^4 for g functions
    if (F_cart == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double F_cart! \n");
    exit(1);
   }

    F_sh = (double *) malloc(6561 * sizeof(double));    // 6561 = 9^4 for g functions
    if (F_sh == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double F_sh! \n");
    exit(1);
   }

    //fprintf(file.out,"p q   ip jp kp lp               %3d %3d    %3d %3d %3d %3d\n",*p1,*q1,ip,jp,kp,lp); fflush(file.out) ;
    //fprintf(file.out,"count_index %3d %3d\n",*count_index,start_index[*count_index]);
    //printf("ip jp kp lp               %3d %3d %3d %3d\n",ip,jp,kp,lp) ;
    //fprintf(file.out,"ip jp kp lp               %3d %3d %3d %3d\n",ip,jp,kp,lp) ;

    R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
    R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
    R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
    R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

    R_CD_1e.comp1 = atoms->cell_vector[kp].comp1 - atoms->cell_vector[lp].comp1;
    R_CD_1e.comp2 = atoms->cell_vector[kp].comp2 - atoms->cell_vector[lp].comp2;
    R_CD_1e.comp3 = atoms->cell_vector[kp].comp3 - atoms->cell_vector[lp].comp3;
    R_CD_1esqrd = double_vec_dot(&R_CD_1e, &R_CD_1e);

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
        //time1 = MPI_Wtime();
        E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
        shells,gaussians,job,file);
        //time2 = MPI_Wtime();
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
            //if (start_index[*count_index] == 0) {
            //fprintf(file.out,"START_index %3d %3d  %3d %3d %3d %3d\n",shell_index,start_index[shell_index],index_i,index_j,index_k,index_l);
            if (start_index[shell_index] == 0) {
              bfposl   += shells->type1_sh[index_l];
              bfposl1  += shells->type_sh[index_l];
              gausposl += shells->ng_sh[index_l];
            //fprintf(file.out,"start_index %3d %3d\n",shell_index,start_index[shell_index]);
                shell_index++;
              //(*count_index)++;
              continue;
             }
                shell_index++;
              //(*count_index)++;
                double scd[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
                double C2x[shells->ng_sh[index_k] * shells->ng_sh[index_l] * (kmax+lmax+1) * (kmax+1) * (lmax+1)];
                double C2y[shells->ng_sh[index_k] * shells->ng_sh[index_l] * (kmax+lmax+1) * (kmax+1) * (lmax+1)];
                double C2z[shells->ng_sh[index_k] * shells->ng_sh[index_l] * (kmax+lmax+1) * (kmax+1) * (lmax+1)];
                double pcd_inv[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
                VECTOR_DOUBLE R_CD[shells->ng_sh[index_k] * shells->ng_sh[index_l]];
                E_coefficients(kp,lp,gk,gl,index_k,index_l,gausposk,gausposl,C2x,C2y,C2z,&C2_max,scd,pcd_inv,&R_CD_1esqrd,R_CD,\
                R,atoms,shells,gaussians,job,file);
                C_max = C1_max * C2_max;

                //if (C_max < 1.0e-16) {
                  //bfposl   += shells->type1_sh[index_l];
                  //bfposl1  += shells->type_sh[index_l];
                  //gausposl += shells->ng_sh[index_l];
                  //continue;
                 //}

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

  time4 = MPI_Wtime();
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
      //non_recursive_ftuvn(mm, 0, f, em, &s_12);
      for (m = 0; m <= mm; m++) {
        for (n = 0; n <= mm; n++) {
          for (pm = 0; pm <= mm; pm++) {
            if (m + n + pm > mm) break;
              //CHANGES2015 CHANGEMEBACK
              fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&em[0][0],s_12);
              //fgtuv_temp = fac1 * f[m][n][pm][0];
               fgtuv[m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4] += fgtuv_temp;
              //p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
             //*p_fgtuv += fgtuv_temp;
              fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
              fgtuv_max = k_one;
             }
            }
           }
          }
         }
  time3 += MPI_Wtime() - time4;

              //if (fgtuv_max * C_max < 1.0e-11) {
              //bfposl   += shells->type1_sh[index_l];
              //bfposl1  += shells->type_sh[index_l];
              //gausposl += shells->ng_sh[index_l];
              //continue;
             //}

  //int ptrc = ptr2 + bfposk1 * nd4 + bfposl1;
  //int ptre = ptr3 + bfposj1 * nd4 + bfposl1;

  time6 = MPI_Wtime();
  mcmurchie_davidson(F_cart, index_i, index_j, index_k, index_l, C1x, C1y, C1z, C2x, C2y, C2z, fgtuv, shells, job, file);
  //mcmurchie_davidson(F_cart, &F[ptrc], &F[ptre], index_i, index_j, index_k, index_l, C1x, C1y, C1z, C2x, C2y, C2z, fgtuv, shells, job, file);
  time5 += MPI_Wtime() - time6;
  free(fgtuv);
  time8 = MPI_Wtime();
  int dimFsh = sheli1 * shelj1 * shelk1 * shell1;
  ResetDoubleArray(F_sh,&dimFsh);
  four_centre_cartesian_to_sh_atom_ijkl(F_cart, F_sh, index_i, index_j, index_k, index_l, shells, job, file);
  //four_centre_cartesian_to_sh(F_cart, F_sh, index_i, index_j, index_k, index_l, shells, job, file);
  time7 += MPI_Wtime() - time8;

//fprintf(file.out,"INDEX %3d %3d %3d %3d    %3d %3d %3d %3d\n",index_i,index_j,index_k,index_l,bfposi1,bfposj1,bfposk1,bfposl1);

  //if (pair_t->lattd[r + *qt] < pair_p->latt_max) {
  //if (pair_p->ptr[pair_t->lattd[r + *qt] * dim2 + pair_t->cell1[q] * dim1 + pair_t->cell2[r + *qt]] > -1) {
  //offset1 = pair_p->Off[*p1];
  //offset2 = pair_p->off[pair_p->ptr[pair_t->lattd[r + *qt] * dim2 + pair_t->cell1[q] * dim1 + pair_t->cell2[r + *qt]]];
  p_F_sh = F_sh;
  for (i = 0; i < sheli1; i++) {
    for (j = 0; j < shelj1; j++) {
      for (k = 0; k < shelk1; k++) {
        for (l = 0; l < shell1; l++) {
          if (fabs(*p_F_sh) > integral_rejection_threshold) {
//if (fabs(*p_F_sh) > 100000000.0) fprintf(file.out,"%3d %3d %3d %3d   %3d %3d %3d %3d  %10.4lf\n",index_i,index_j,index_k,index_l,i,j,k,l,*p_F_sh);
//ADDED
//*count_index = 0;
          integral_list->value[*count_index + count] = *p_F_sh ;
          integral_list->i[*count_index + count] = bfposi1 + i;
          integral_list->j[*count_index + count] = bfposj1 + j;
          integral_list->k[*count_index + count] = bfposk1 + k;
          integral_list->l[*count_index + count] = bfposl1 + l;
//fprintf(file.out,"COUNT %3d %3d %3d  %3d   %3d %3d %3d %3d %12.6lf\n",ip,jp,kp,lp,bfposi1+i,bfposj1+j,bfposk1+k,bfposl1+l,*p_F_sh);
//ADDED
          ////integral_list->value[count] = *p_F_sh ;
          ////integral_list->i[count] = bfposi1 + i;
          ////integral_list->j[count] = bfposj1 + j;
          ////integral_list->k[count] = bfposk1 + k;
          ////integral_list->l[count] = bfposl1 + l;
          //integral_list->gj[count] = pair_p->latt2[p];
          //integral_list->gk[count] = pair_t->latt1[q];
          //integral_list->gl[count] = pair_t->latt2[r + *qt];
          //integral_list->offset1[count] = offset1;
          //integral_list->offset2[count] = offset2;
          count++;
         }
          p_F_sh++;
         }
        }
       }
      }
     //}
   // }

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

  integral_list->num = count;
//for(i=0;i<6;i++) { for(j=0;j<6;j++) { for(k=0;k<6;k++) { for(l=0;l<6;l++) { 
//fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",i,j,k,l,F_cart[i*216+j*36+k*6+l]);}}}}
//for(i=0;i<5;i++) { for(j=0;j<5;j++) { for(k=0;k<5;k++) { for(l=0;l<5;l++) { 
//fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",i,j,k,l,F_sh[i*125+j*25+k*5+l]);}}}}

  free(F_cart);
  free(F_sh);

  time2 = MPI_Wtime() - time1;

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out, "Time for %d integrals no sym %10.4e\n", count,time2-time1);
  //if (job->taskid == 1) printf("Time for %d integrals no sym %10.4e %10.4e %10.4e %10.4e\n", count,time2,time3,time5,time7);
  //printf("Time for %d integrals no sym %10.4e\n", count,time2-time1);

}

void integrals_molecule_screen(double *F_sh, int ip, int jp, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // compute 4-centre integrals for a quad of atoms (ij|ij) using pair_p as indices
  // starting point is pair_p->posn[*p1]

  int index_i, index_j, i4, j4;
  int bfposi, bfposj;
  int bfposi1, bfposj1;
  int gausposi, gausposj;
  int count;
  int i, j, k, l, m, n, p, pm, q;
  int mm, mm0, mm1, mm2, mm3, mm4;
  int gi, gj;
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
  double f[13][13][13][13];
  double *F_cart, *p_F_cart, *p_F_sh, *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, s_12;

  time1 = MPI_Wtime();

    gi = 0;
    gj = 0;

    F_cart = (double *) calloc(nd1 * nd2 * nd1 * nd2, sizeof(double));
    if (F_cart == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double F_cart! \n");
    exit(1);
   }
    
    R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
    R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
    R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
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
            //non_recursive_ftuvn(mm, 0, f, em, &s_12);
            for (m = 0; m <= mm; m++) {
              for (n = 0; n <= mm; n++) {
                for (pm = 0; pm <= mm; pm++) {
                  if (m + n + pm > mm) break;
                  //CHANGES2015 CHANGEMEBACK
                  fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&em[0][0],s_12);
                  //fgtuv_temp = fac1 * f[m][n][pm][0];
                  p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
                 *p_fgtuv += fgtuv_temp;
                  fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                  fgtuv_max = k_one;
                 }
                }
               }
              }
             }
        mcmurchie_davidson_screen(F_cart, index_i, index_j, bfposi, bfposj, nd2, C1x, C1y, C1z, fgtuv, shells, job, file);
        free(fgtuv);
        four_centre_cartesian_to_sh_atom_ijkl_screen(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        //two_centre_cartesian_to_sh(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i
//for(i=0;i<nd1;i++) { for(j=0;j<nd2;j++) { for(k=0;k<nd1;k++) { for(l=0;l<nd2;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_cart[i*nd2*nd1*nd2+j*nd1*nd2+k*nd2+l]);}}}}
//for(i=0;i<nd3;i++) { for(j=0;j<nd4;j++) { for(k=0;k<nd3;k++) { for(l=0;l<nd4;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_sh[i*nd4*nd3*nd4+j*nd3*nd4+k*nd4+l]);}}}}

      free(F_cart);

}

void pack_write_molecule_2e_integrals(INTEGRAL_LIST *integral_list, QUAD_TRAN *quad, FILE *integrals, JOB_PARAM *job, FILES file)

{

int i, *packed_integral_indices, packed_quads[6 * quad->tot + 2];

          packed_integral_indices = (int *) malloc(2 * integral_list->num * sizeof(int));
          if (packed_integral_indices == NULL) { if (job->taskid == 0)
          fprintf(stderr, "ERROR: not enough memory for packed_integral_indices in pack_write_molecule_2e_integrals\n");
          MPI_Finalize(); exit(1); }

          packed_quads[0] = quad->tot;
          packed_quads[1] = integral_list->num;
        for (i = 0; i < quad->tot; i++) {
          packed_quads[6 * i + 2] = quad->cell1[i];
          packed_quads[6 * i + 3] = quad->cell2[i];
          packed_quads[6 * i + 4] = quad->cell3[i];
          packed_quads[6 * i + 5] = quad->cell4[i];
          packed_quads[6 * i + 6] = quad->p[i];
          packed_quads[6 * i + 7] = quad->k[i];
          //printf("PACKED %3d %3d %3d %3d   %3d %3d\n",\
          quad->cell1[i],quad->cell2[i],quad->cell3[i],quad->cell4[i],quad->p[i],quad->k[i]);
          //fprintf(file.out,"PACKED %3d %3d %3d %3d   %3d %3d\n",\
          quad->cell1[i],quad->cell2[i],quad->cell3[i],quad->cell4[i],quad->p[i],quad->k[i]);
         }

        for (i = 0; i < integral_list->num; i++) {
          packed_integral_indices[2 * i]     = integral_list->i[i] * 256 + integral_list->j[i];
          packed_integral_indices[2 * i + 1] = integral_list->k[i] * 256 + integral_list->l[i];
          //fprintf(file.out,"INT %3d %10.4f\n",i,integral_list->value[i]);
          //printf("INT %3d %10.4f\n",i,integral_list->value[i]);
         }
          fwrite(packed_quads,sizeof(int), 6 * quad->tot + 2, integrals);
          fwrite(packed_integral_indices,sizeof(int), 2 * integral_list->num,integrals);
          fwrite(&integral_list->value[0],sizeof(double),integral_list->num,integrals);

          free(packed_integral_indices);

}

void read_unpack_molecule_2e_integrals(INTEGRAL_LIST *integral_list, QUAD_TRAN *quad, FILE *integrals, JOB_PARAM *job, FILES file)

{

int i, *packed_integral_indices, packed_quads[6 * quad->tot];
size_t result;

          result = fread(packed_quads,sizeof(int),6 * quad->tot,integrals);
          
        for (i = 0; i < quad->tot; i++) {
          quad->cell1[i] = packed_quads[6 * i];
          quad->cell2[i] = packed_quads[6 * i + 1];
          quad->cell3[i] = packed_quads[6 * i + 2];
          quad->cell4[i] = packed_quads[6 * i + 3];
          quad->p[i] = packed_quads[6 * i + 4];
          quad->k[i] = packed_quads[6 * i + 5];
        }

          packed_integral_indices = (int *) malloc(2 * integral_list->num * sizeof(int));
          if (packed_integral_indices == NULL) { if (job->taskid == 0)
          fprintf(stderr, "ERROR: not enough memory for packed_integral_indices in pack_write_molecule_2e_integrals\n");
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

