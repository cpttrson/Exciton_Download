/*
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <fstream>
#include <mpi.h>
#include "mycomplex.h"
#include "mylogical.h"
#include "conversion_factors.h"
#include "LIMITS.h"
#include "PRINT_UTIL.h"
#include "CRYSTAL1.h"
#include "SYMMETRY_ADAPTATION.h"
#include "ROTATION_OPERATORS.h"
#include "TOOLS.h"
#include "ALLOCATE_MEMORY.h"
*/
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "ROTATIONS_MOLECULE.h"

using namespace std;

void rotate_single(int op, int ip, Complex *F, Complex *F_rotate, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Rotates a contracted triple by operator op with no permutation                         *
  // * Called with F index as auxillary basis                                                 *
  // ******************************************************************************************

int i, nd1, index_i, shelposi, bfposi, oppshift1, sheli1;
int *p_i, *p_i1;
double *p_rot1;
Complex *p_F, *p_F_rotate;

  nd1 = atoms->bfnnumb_sh[ip];
  for (i = 0; i < nd1; i++) {
  F_rotate[i] = Complex(k_zero, k_zero);
 }

  //op    = symmetry->inverse[knet->opr[k]];
  //opinv = knet->opr[k];
  //O  = atom_p->O[i * symmetry->number_of_operators + op];
  //temp  = double_vec_dot(&knet->cart[k],&R->vec_ai[O]) ;
  //if (knet->trs[k] == 1) temp *=-k_one;
  //temp1 = Complex(cos(temp),sin(temp)) ;

  shelposi = atoms->shelposn[ip];
  bfposi = 0;
  for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
    sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
    p_i1   = symmetry->ind_i + oppshift1;
    p_i    = symmetry->ind_j + oppshift1;
    p_rot1 = symmetry->rot   + oppshift1;
    for (i = 0; i < sheli1; i++) {
      F_rotate[bfposi + *p_i1] += *p_rot1 * F[bfposi + *p_i];
      //fprintf(file.out,"index_i %3d i %3d bfposi %3d *p_i1 %3d *p_i %3d F.real %14.8f F.imag %14.8f F_rot %14.8f %14.8f\n",\
      index_i,i,bfposi,*p_i1,*p_i,(F[bfposi + *p_i]).real(),(F[bfposi+*p_i]).imag(),(F_rotate[bfposi + *p_i1]).real(),\
      (F_rotate[bfposi + *p_i1]).imag());
      p_i++;
      p_i1++;
      p_rot1++;
    }
   bfposi += shells->type1[index_i];
  } // close loop over index_i

}

void rotate_permute_expand_pair(int p, PAIR_TRAN *pair_p, double *F_reduced, double *F, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // generates a full operator from a reduced operator using permutation and rotation operations

  int i, j, k, l, nd, nd1, nd2;
  int ip, jp, gj, op, pm;
  int index_i, index_j, shelposi, shelposj, bfposi, bfposj;
  int shift1, shift2, oppshift1, oppshift2;
  int sheli, shelj;
  int sheli1, shelj1;
  int shelpos1, shelpos2, limit1, limit2;
  int *p_i, *p_i1, *p_j, *p_j1;
  double *p_rot1, *p_rot2, *p_F, *p_F_reduced, *F_rotate, *p_F_rotate;

  l = pair_p->posn[p];
  ip = pair_p->cell1[l];
  jp = pair_p->cell2[l];
  nd1 = atoms->bfnnumb_sh[ip];
  nd2 = atoms->bfnnumb_sh[jp];

  F_rotate = (double *) malloc(nd1 * nd2 * sizeof(double));
  if (F_rotate == NULL) {  fprintf(stderr, "ERROR: not enough memory for double F_rotate! \n"); exit(1); }

  for (k = 0; k < pair_p->numb[p]; k++) {

    ip = pair_p->cell1[l + k];
    jp = pair_p->cell2[l + k];
    gj = pair_p->latt2[l + k];
    op = pair_p->k[l + k];
    pm = pair_p->p[l + k];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];

    //fprintf(file.out,"  ip jp %3d %3d %3d   %3d %3d  %3d %3d\n",ip,jp,gj,op,pm,nd1,nd2); fflush(file.out);

    p_F_rotate = F_rotate;
    for (i = 0; i < nd1 * nd2; i++) {
      *p_F_rotate = k_zero;
        p_F_rotate++;
       }
 
      switch (pm) {

        case 0:

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
        //CHANGES2015oppshift1 = *(symmetry->op_shift + op * 5 + shells->ord[index_i]);
        //oppshift2 = *(symmetry->op_shift + op * 5 + shells->ord[index_j]);
        //sheli1 = *(symmetry->num_ij + op * 5 + shells->ord[index_i]);
        //shelj1 = *(symmetry->num_ij + op * 5 + shells->ord[index_j]);
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
        p_i1   = symmetry->ind_i + oppshift1;
        p_i    = symmetry->ind_j + oppshift1;
        p_rot1 = symmetry->rot   + oppshift1;
        for (i = 0; i < sheli1; i++) {
          p_j1   = symmetry->ind_i + oppshift2;
          p_j    = symmetry->ind_j + oppshift2;
          p_rot2 = symmetry->rot   + oppshift2;
          for (j = 0; j < shelj1; j++) {
            p_F_reduced = F_reduced + (bfposi + *p_i)  * nd2 + bfposj + *p_j;
            p_F_rotate  = F_rotate  + (bfposi + *p_i1) * nd2 + bfposj + *p_j1;   // no transpose
           *p_F_rotate += *p_rot1 * *p_rot2 * *p_F_reduced;
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj += shells->type1[index_j];
      } // close loop over index_j
      bfposi += shells->type1[index_i];
    } // close loop over index_i

       break;

        case 1:

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
        //CHANGES2015oppshift1 = *(symmetry->op_shift + op * 5 + shells->ord[index_i]);
        //oppshift2 = *(symmetry->op_shift + op * 5 + shells->ord[index_j]);
        //sheli1 = *(symmetry->num_ij + op * 5 + shells->ord[index_i]);
        //shelj1 = *(symmetry->num_ij + op * 5 + shells->ord[index_j]);
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
        p_i1   = symmetry->ind_i + oppshift1;
        p_i    = symmetry->ind_j + oppshift1;
        p_rot1 = symmetry->rot   + oppshift1;
        for (i = 0; i < sheli1; i++) {
          p_j1   = symmetry->ind_i + oppshift2;
          p_j    = symmetry->ind_j + oppshift2;
          p_rot2 = symmetry->rot   + oppshift2;
          for (j = 0; j < shelj1; j++) {
            p_F_reduced = F_reduced + (bfposj + *p_j)  * nd1 + bfposi + *p_i;     // transpose
            p_F_rotate  = F_rotate  + (bfposi + *p_i1) * nd2 + bfposj + *p_j1;
           *p_F_rotate += *p_rot1 * *p_rot2 * *p_F_reduced;
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj += shells->type1[index_j];
      } // close loop over index_j
      bfposi += shells->type1[index_i];
    } // close loop over index_i

      break;

     } // end switch

    p_F_rotate = F_rotate;
    p_F = F + nd1 * nd2 * k; // Put rotated pairs in consecutive order. Place in proper order only in FT'd arrays

    //fprintf(file.out,"p %3d ip %3d jp %3d gj %3d op %3d pm %3d\n",p,ip,jp,gj,op,pm) ;
    for (i = 0; i < nd1; i++) {
      for (j = 0; j < nd2; j++) {
        *p_F = *p_F_rotate;
        //fprintf(file.out,"%10.4lf ",*p_F) ;
        p_F++;
        p_F_rotate++;
      }
     //fprintf(file.out,"\n");
     }

  } // close loop on k

  free(F_rotate);

}

void rotate_permute_expand_pair_complex(int p, PAIR_TRAN *pair_p, Complex *F_reduced, Complex *F, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // generates a full operator from a reduced operator using permutation and rotation operations

  int i, j, k, l, nd, nd1, nd2;
  int ip, jp, gj, op, pm;
  int index_i, index_j, shelposi, shelposj, bfposi, bfposj;
  int shift1, shift2, oppshift1, oppshift2;
  int sheli, shelj;
  int sheli1, shelj1;
  int shelpos1, shelpos2, limit1, limit2;
  int *p_i, *p_i1, *p_j, *p_j1;
  double *p_rot1, *p_rot2;
  Complex *p_F, *p_F_reduced, *F_rotate, *p_F_rotate;

int lmax;

  lmax = 26;
  for (i = 2; i <= job->l_max; i++) 
  lmax += (2 * i + 1) * (2 * i + 1);
//printf("lmax line 2084 rotation_operators %4d %4d\n",job->lmax,lmax);

  l = pair_p->posn[p];
  ip = pair_p->cell1[l];
  jp = pair_p->cell2[l];
  nd1 = atoms->bfnnumb_sh[ip];
  nd2 = atoms->bfnnumb_sh[jp];

  F_rotate = (Complex *) malloc(nd1 * nd2 * sizeof(Complex));
  if (F_rotate == NULL) {  fprintf(stderr, "ERROR: not enough memory for Complex F_rotate! \n"); exit(1); }

  for (k = 0; k < pair_p->numb[p]; k++) {

    ip = pair_p->cell1[l + k];
    jp = pair_p->cell2[l + k];
    gj = pair_p->latt2[l + k];
    op = pair_p->k[l + k];
    pm = pair_p->p[l + k];
    nd1 = atoms->bfnnumb_sh[ip];
    nd2 = atoms->bfnnumb_sh[jp];

    //fprintf(file.out,"  ip jp %3d %3d %3d   %3d %3d  %3d %3d\n",ip,jp,gj,op,pm,nd1,nd2); fflush(file.out);

    p_F_rotate = F_rotate;
    for (i = 0; i < nd1 * nd2; i++) {
      *p_F_rotate = k_zero;
        p_F_rotate++;
       }
 
      switch (pm) {

        case 0:

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
        //CHANGES2015oppshift1 = *(symmetry->op_shift + op * 5 + shells->ord[index_i]);
        //oppshift2 = *(symmetry->op_shift + op * 5 + shells->ord[index_j]);
        //sheli1 = *(symmetry->num_ij + op * 5 + shells->ord[index_i]);
        //shelj1 = *(symmetry->num_ij + op * 5 + shells->ord[index_j]);
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
        p_i1   = symmetry->ind_i + oppshift1;
        p_i    = symmetry->ind_j + oppshift1;
        p_rot1 = symmetry->rot   + oppshift1;
//if (oppshift1 + sheli1 >= lmax * symmetry->number_of_operators || oppshift2 + shelj1 >= lmax * symmetry->number_of_operators) \
fprintf(file.out,"oppshift1,2 %5d %5d lmax %5d\n",oppshift1,oppshift2,lmax*symmetry->number_of_operators);
        for (i = 0; i < sheli1; i++) {
          p_j1   = symmetry->ind_i + oppshift2;
          p_j    = symmetry->ind_j + oppshift2;
          p_rot2 = symmetry->rot   + oppshift2;
          for (j = 0; j < shelj1; j++) {
            p_F_reduced = F_reduced + (bfposi + *p_i)  * nd2 + bfposj + *p_j;
            //if( (bfposi + *p_i)  * nd2 + bfposj + *p_j >= nd1 * nd2 || (bfposi + *p_i1) * nd2 + bfposj + *p_j1 >= nd1 * nd2) \
            fprintf(file.out,"nd1*nd2 %5d (bfposi %3d  *p_i)  %3d nd2 %3d bfposj %3d *p_j %3d %5d\n",\
            nd1 * nd2,bfposi,*p_i, nd2, bfposj, *p_j,(bfposi + *p_i)  * nd2 + bfposj + *p_j);
            p_F_rotate  = F_rotate  + (bfposi + *p_i1) * nd2 + bfposj + *p_j1;   // no transpose
           *p_F_rotate += *p_rot1 * *p_rot2 * *p_F_reduced;
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj += shells->type1[index_j];
      } // close loop over index_j
      bfposi += shells->type1[index_i];
    } // close loop over index_i

       break;

        case 1:

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
        //CHANGES2015oppshift1 = *(symmetry->op_shift + op * 5 + shells->ord[index_i]);
        //oppshift2 = *(symmetry->op_shift + op * 5 + shells->ord[index_j]);
        //sheli1 = *(symmetry->num_ij + op * 5 + shells->ord[index_i]);
        //shelj1 = *(symmetry->num_ij + op * 5 + shells->ord[index_j]);
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
        p_i1   = symmetry->ind_i + oppshift1;
        p_i    = symmetry->ind_j + oppshift1;
        p_rot1 = symmetry->rot   + oppshift1;
        for (i = 0; i < sheli1; i++) {
          p_j1   = symmetry->ind_i + oppshift2;
          p_j    = symmetry->ind_j + oppshift2;
          p_rot2 = symmetry->rot   + oppshift2;
          for (j = 0; j < shelj1; j++) {
            p_F_reduced = F_reduced + (bfposj + *p_j)  * nd1 + bfposi + *p_i;     // transpose
            p_F_rotate  = F_rotate  + (bfposi + *p_i1) * nd2 + bfposj + *p_j1;
           *p_F_rotate += *p_rot1 * *p_rot2 * *p_F_reduced;
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj += shells->type1[index_j];
      } // close loop over index_j
      bfposi += shells->type1[index_i];
    } // close loop over index_i

      break;

     } // end switch

    p_F_rotate = F_rotate;
    p_F = F + nd1 * nd2 * k; // Put rotated pairs in consecutive order. Place in proper order only in FT'd arrays

    //fprintf(file.out,"p %3d ip %3d jp %3d gj %3d op %3d pm %3d\n",p,ip,jp,gj,op,pm) ;
    for (i = 0; i < nd1; i++) {
      for (j = 0; j < nd2; j++) {
        *p_F = *p_F_rotate;
        //fprintf(file.out,"%10.4lf ",*p_F) ;
        p_F++;
        p_F_rotate++;
      }
     //fprintf(file.out,"\n");
     }

  } // close loop on k

  free(F_rotate);

}

void rotate_permute_expand_tensor1(int p, PAIR_TRAN *pair_p, double *F_reduced, double *F, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{
int i, j, k, t, l, nd1, nd2 ;
int ip, jp, gj, op, pm ;
int index_i, index_j, shelposi, shelposj, bfposi, bfposj ;
int shift1, shift2, oppshift1, oppshift2, oppshift3 ;
int sheli, shelj ;
int sheli1, shelj1 ;
int *p_i, *p_i1, *p_j, *p_j1, *p_k, *p_k1 ;
double *p_rot1, *p_rot2, *p_rot3, *p_F, *p_F_reduced, *F_rotate, *p_F_rotate ;

   l = pair_p->posn[p];
  ip = pair_p->cell1[l];
  jp = pair_p->cell2[l];
  nd1 = atoms->bfnnumb_sh[ip];
  nd2 = atoms->bfnnumb_sh[jp];

      F_rotate = (double *) malloc(3 * nd1 * nd2 * sizeof(double)) ;
      if (F_rotate==NULL) { fprintf(stderr, "ERROR: not enough memory for double F_rotate! \n") ; exit(1) ; }

      for (k = 0; k < pair_p->numb[p]; k++) {
      ip = pair_p->cell1[l + k];
      jp = pair_p->cell2[l + k];
      gj = pair_p->latt2[l + k];
      op = pair_p->k[l + k];
      pm = pair_p->p[l + k];
      nd1 = atoms->bfnnumb_sh[ip];
      nd2 = atoms->bfnnumb_sh[jp];

      p_F_rotate = F_rotate;
      for (i = 0; i < 3 * nd1 * nd2; i++) {
        *p_F_rotate = k_zero;
         p_F_rotate++;
        }

      switch (pm) {

        case 0:

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
        //CHANGES2015oppshift1 = *(symmetry->op_shift + op * 5 + shells->ord[index_i]);
        //oppshift2 = *(symmetry->op_shift + op * 5 + shells->ord[index_j]);
        //oppshift3 = *(symmetry->op_shift + op * 5 + 2);
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
        oppshift3 = *(symmetry->op_shift + op * (job->l_max + 2) + 2);
//oppshift3 = *(symmetry->op_shift + opinv * 5 + 2);
        //CHANGES2015sheli1 = *(symmetry->num_ij + op * 5 + shells->ord[index_i]);
        //shelj1 = *(symmetry->num_ij + op * 5 + shells->ord[index_j]);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
        p_k1   = symmetry->ind_i + oppshift3 ;
        p_k    = symmetry->ind_j + oppshift3 ;
        p_rot3 = symmetry->rot   + oppshift3 ;
        for (t = 0; t < 3; t++) {
          p_i1   = symmetry->ind_i + oppshift1;
          p_i    = symmetry->ind_j + oppshift1;
          p_rot1 = symmetry->rot   + oppshift1;
          for (i = 0; i < sheli1; i++) {
            p_j1   = symmetry->ind_i + oppshift2;
            p_j    = symmetry->ind_j + oppshift2;
            p_rot2 = symmetry->rot   + oppshift2;
            for (j = 0; j < shelj1; j++) {
              p_F_reduced = F_reduced + *p_k  * nd1 * nd2 + (bfposi + *p_i)  * nd2 + bfposj + *p_j;
              p_F_rotate  = F_rotate  + *p_k1 * nd1 * nd2 + (bfposi + *p_i1) * nd2 + bfposj + *p_j1;   // no transpose
              *p_F_rotate += *p_rot1 * *p_rot2 * *p_rot3 * *p_F_reduced ;
              p_j++;
              p_j1++;
              p_rot2++;
            }
            p_i++;
            p_i1++;
            p_rot1++;
          }
          p_k++;
          p_k1++ ;
          p_rot3++ ;
         }
        bfposj += shells->type1[index_j];
       } // close loop over index_j
      bfposi += shells->type1[index_i];
     } // close loop over index_i

        break;

        case 1:

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
        //CHANGES2015oppshift1 = *(symmetry->op_shift + op * 5 + shells->ord[index_i]);
        //oppshift2 = *(symmetry->op_shift + op * 5 + shells->ord[index_j]);
        //oppshift3 = *(symmetry->op_shift + op * 5 + 2);
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
        oppshift3 = *(symmetry->op_shift + op * (job->l_max + 2) + 2);
//oppshift3 = *(symmetry->op_shift + opinv * 5 + 2);
//CHP2014
        //CHANGES2015sheli1 = *(symmetry->num_ij + op * 5 + shells->ord[index_i]);
        //shelj1 = *(symmetry->num_ij + op * 5 + shells->ord[index_j]);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
        p_k1   = symmetry->ind_i + oppshift3 ;
        p_k    = symmetry->ind_j + oppshift3 ;
        p_rot3 = symmetry->rot   + oppshift3 ;
        for (t = 0; t < 3; t++) {
          p_i1   = symmetry->ind_i + oppshift1;
          p_i    = symmetry->ind_j + oppshift1;
          p_rot1 = symmetry->rot   + oppshift1;
          for (i = 0; i < sheli1; i++) {
            p_j1   = symmetry->ind_i + oppshift2;
            p_j    = symmetry->ind_j + oppshift2;
            p_rot2 = symmetry->rot   + oppshift2;
            for (j = 0; j < shelj1; j++) {
              p_F_reduced  = F_reduced + *p_k  * nd1 * nd2 + (bfposj + *p_j)  * nd1 + bfposi + *p_i;     // transpose
              p_F_rotate   = F_rotate  + *p_k1 * nd1 * nd2 + (bfposi + *p_i1) * nd2 + bfposj + *p_j1;
              *p_F_rotate += *p_rot1 * *p_rot2 * *p_rot3 * *p_F_reduced ;
              p_j++;
              p_j1++;
              p_rot2++;
            }
            p_i++;
            p_i1++;
            p_rot1++;
          }
          p_k++;
          p_k1++ ;
          p_rot3++ ;
        }
        bfposj += shells->type1[index_j];
      } // close loop over index_j
      bfposi += shells->type1[index_i];
    } // close loop over index_i

        break;

    } // close switch

    p_F_rotate = F_rotate;
    p_F = F + 3 * nd1 * nd2 * k; // Put rotated pairs in consecutive order. Place in proper order only in FT'd arrays

   for (t = 0; t < 3; t++) {
    for (i = 0; i < nd1; i++) {
      for (j = 0; j < nd2; j++) {
        *p_F = *p_F_rotate;
        //fprintf(file.out,"p latt2 i j *p_rotate%3d %3d %3d %3d %12.6lf\n",p,pair_p->latt2[k],i,j,*p_F_rotate) ;
        p_F++;
        p_F_rotate++;
      }
      //fprintf(file.out,"\n") ;
    } // close loop on t

    }
    //fprintf(file.out,"\n") ;
    //fprintf(file.out,"\n") ;
  } // close loop on k

  free(F_rotate);

}

void rotate_permute_expand_tensor_3(double *matrix_reduced, double *matrix_k, PAIR_TRAN *pair_p, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{
      int i, j, l, m, n, p, dim2, dim3, dim4 = atoms->number_of_sh_bfns_in_unit_cell * atoms->number_of_sh_bfns_in_unit_cell;
      int atm1, atm2, nd1, nd2, pm;
      double *p_matrix_reduced, *p_matrix, *matrix, *p_matrix_k ;
      double time1 = MPI_Wtime();

  // Routine similar to fourier_transform_3 which expands, permutes and rotates pairs for finite systems

      dim2 = 0 ;

      for (i = 0; i < pair_p->nump; i++) {
      atm1 = pair_p->cell1[pair_p->posn[i]];
      atm2 = pair_p->cell2[pair_p->posn[i]];
      nd1 = atoms->bfnnumb_sh[atm1];
      nd2 = atoms->bfnnumb_sh[atm2];
      dim2 += nd1 * nd2 * pair_p->numb[i];
     }

      matrix = (double *) malloc(dim2 * 3 * sizeof(double)) ;
      if (matrix==NULL) { fprintf(stderr, "ERROR: not enough memory for double matrix! \n") ; exit(1) ; }

      p_matrix = matrix ;
      for (i = 0; i < dim2 * 3; i++) {
      *p_matrix = k_zero ; p_matrix++ ; }

      p_matrix_k = matrix_k ;
      for ( i = 0; i < dim4 * 3; i++) {
      *p_matrix_k = k_zero ;
       p_matrix_k++ ; }

      dim3 = 0 ;

      for ( p = 0; p < pair_p->nump; p++) {
      i = pair_p->posn[p];
      nd1 = atoms->bfnnumb_sh[pair_p->cell1[i]];
      nd2 = atoms->bfnnumb_sh[pair_p->cell2[i]];
      rotate_permute_expand_tensor1(p, pair_p, &matrix_reduced[3 * dim3], matrix, atoms, shells, symmetry, job, file);

      dim3 += nd1 * nd2 ;

      //fprintf(file.out,"pair_p %d   [%3d][%3d]  [%3d]\n",i,pair_p->cell1[i],pair_p->cell2[i],pair_p->latt2[i]);
      for (j = 0; j < pair_p->numb[p]; j++) {
      //fprintf(file.out,"pair_p %d   [%3d][%3d]  [%3d]\n",i,pair_p->cell1[i + j],pair_p->cell2[i + j],pair_p->latt2[i + j]);
      nd1 = atoms->bfnnumb_sh[pair_p->cell1[i + j]];
      nd2 = atoms->bfnnumb_sh[pair_p->cell2[i + j]];
      for (n = 0; n < 3; n++) {
      p_matrix      = matrix      + nd1 * nd2 * (3 * j + n) ;
      p_matrix_k = matrix_k + atoms->bfnposn_sh[pair_p->cell1[i + j]] * \
      atoms->number_of_sh_bfns_in_unit_cell + atoms->bfnposn_sh[pair_p->cell2[i + j]] + dim4 * n;
      //fprintf(file.out,"%d\n",n);
      //fprintf(file.out,"M_k address %4d %4d %4d %4d\n",atoms->bfnposn_sh[pair_p->cell1[i + j]] * \
      atoms->number_of_sh_bfns_in_unit_cell , atoms->bfnposn_sh[pair_p->cell2[i + j]] , dim4 * n,\
      atoms->number_of_sh_bfns_in_unit_cell - nd2);
      for (l = 0; l < nd1; l++) {
      for (m = 0; m < nd2; m++) {
      *p_matrix_k += *p_matrix; 
       p_matrix++ ;
       p_matrix_k++ ;
      }
       p_matrix_k += atoms->number_of_sh_bfns_in_unit_cell - nd2;
      } // close loop on l
      } // close loop on n
      } // close loop on j
      } // close loop on p

      free(matrix) ;

}

void rotate_permute_triple_ij_alpha(int *ip, int *jp, int *alpha, int *op, int *pm, double *F_reduced, double *F, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, ATOM *atoms_ax, SHELL *shells_ax, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // generates a full operator <ij|alpha> from a reduced operator using permutation and rotation operations

  int i, j, k, l, nd1, nd2, nd3;
  int index_i, index_j, index_k, shelposi, shelposj, shelposk, bfposi, bfposj, bfposk;
  int shift1, shift2, shift3, oppshift1, oppshift2, oppshift3;
  int sheli, shelj, shelk;
  int sheli1, shelj1, shelk1;
  int shelpos1, shelpos2, shelpos3, limit1, limit2, limit3;
  int *p_i, *p_i1, *p_j, *p_j1, *p_k, *p_k1;
  double *p_rot1, *p_rot2, *p_rot3, *p_F, *p_F_reduced, *F_rotate, *p_F_rotate;

    nd1 = atoms->bfnnumb_sh[*ip];
    nd2 = atoms->bfnnumb_sh[*jp];
    nd3 = atoms_ax->bfnnumb_sh[*alpha];

    for (i = 0; i < nd1 * nd2 * nd3; i++) F[i] = k_zero;
 
      switch (*pm) {

        case 0: // no permutation transpose

    shelposi = atoms->shelposn[*ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[*ip]; index_i++) {
      shelposj = atoms->shelposn[*jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[*jp]; index_j++) {
        shelposk = atoms_ax->shelposn[*alpha];
        bfposk = 0;
        for (index_k = shelposk; index_k < shelposk + atoms_ax->nshel[*alpha]; index_k++) {
          oppshift1 = *(symmetry->op_shift + *op * (job->l_max + 2) + shells->ord[index_i]);
          oppshift2 = *(symmetry->op_shift + *op * (job->l_max + 2) + shells->ord[index_j]);
          oppshift3 = *(symmetry->op_shift + *op * (job->l_max + 2) + shells_ax->ord[index_k]);
          sheli1 = *(symmetry->num_ij + *op * (job->l_max + 2) + shells->ord[index_i]);
          shelj1 = *(symmetry->num_ij + *op * (job->l_max + 2) + shells->ord[index_j]);
          shelk1 = *(symmetry->num_ij + *op * (job->l_max + 2) + shells_ax->ord[index_k]);
          p_i1   = symmetry->ind_i + oppshift1;
          p_i    = symmetry->ind_j + oppshift1;
          p_rot1 = symmetry->rot   + oppshift1;
          for (i = 0; i < sheli1; i++) {
            p_j1   = symmetry->ind_i + oppshift2;
            p_j    = symmetry->ind_j + oppshift2;
            p_rot2 = symmetry->rot   + oppshift2;
            for (j = 0; j < shelj1; j++) {
              p_k1   = symmetry->ind_i + oppshift3;
              p_k    = symmetry->ind_j + oppshift3;
              p_rot3 = symmetry->rot   + oppshift3;
              for (k = 0; k < shelk1; k++) {
                p_F_reduced = F_reduced + (bfposi + *p_i)  * nd2 * nd3 + (bfposj + *p_j)  * nd3 + bfposk + *p_k;
                F[(bfposi + *p_i1) * nd2 * nd3 + (bfposj + *p_j1) * nd3 + bfposk + *p_k1] += *p_rot1 * *p_rot2 * *p_rot3 * *p_F_reduced;
                p_k++;
                p_k1++;
                p_rot3++;
               }
               p_j++;
               p_j1++;
               p_rot2++;
              }
              p_i++;
              p_i1++;
              p_rot1++;
             }
            bfposk += shells_ax->type1[index_k];
           } // close loop over index_k
          bfposj += shells->type1[index_j];
         } // close loop over index_j
        bfposi += shells->type1[index_i];
       } // close loop over index_i

       break;

        case 1: // permutation transpose

    shelposi = atoms->shelposn[*ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[*ip]; index_i++) {
      shelposj = atoms->shelposn[*jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[*jp]; index_j++) {
        shelposk = atoms_ax->shelposn[*alpha];
        bfposk = 0;
        for (index_k = shelposk; index_k < shelposk + atoms_ax->nshel[*alpha]; index_k++) {
          oppshift1 = *(symmetry->op_shift + *op * (job->l_max + 2) + shells->ord[index_i]);
          oppshift2 = *(symmetry->op_shift + *op * (job->l_max + 2) + shells->ord[index_j]);
          oppshift3 = *(symmetry->op_shift + *op * (job->l_max + 2) + shells_ax->ord[index_k]);
          sheli1 = *(symmetry->num_ij + *op * (job->l_max + 2) + shells->ord[index_i]);
          shelj1 = *(symmetry->num_ij + *op * (job->l_max + 2) + shells->ord[index_j]);
          shelk1 = *(symmetry->num_ij + *op * (job->l_max + 2) + shells_ax->ord[index_k]);
          p_i1   = symmetry->ind_i + oppshift1;
          p_i    = symmetry->ind_j + oppshift1;
          p_rot1 = symmetry->rot   + oppshift1;
          for (i = 0; i < sheli1; i++) {
            p_j1   = symmetry->ind_i + oppshift2;
            p_j    = symmetry->ind_j + oppshift2;
            p_rot2 = symmetry->rot   + oppshift2;
            for (j = 0; j < shelj1; j++) {
              p_k1   = symmetry->ind_i + oppshift3;
              p_k    = symmetry->ind_j + oppshift3;
              p_rot3 = symmetry->rot   + oppshift3;
              for (k = 0; k < shelk1; k++) {
                p_F_reduced = F_reduced + (bfposj + *p_j)  * nd1 * nd3 + (bfposi + *p_i)  * nd3 + bfposk + *p_k; 
                F[(bfposi + *p_i1) * nd2 * nd3 + (bfposj + *p_j1) * nd3 + bfposk + *p_k1] += *p_rot1 * *p_rot2 * *p_rot3 * *p_F_reduced;
                p_k++;
                p_k1++;
                p_rot3++;
               }
               p_j++;
               p_j1++;
               p_rot2++;
              }
              p_i++;
              p_i1++;
              p_rot1++;
             }
            bfposk += shells_ax->type1[index_k];
           } // close loop over index_k
          bfposj += shells->type1[index_j];
         } // close loop over index_j
        bfposi += shells->type1[index_i];
       } // close loop over index_i

      break;

     } // end switch

}

void rotate_sum_block(double *F, double *F_rotated, int ip, int jp, int op, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, nd1, nd2;
  int index_i, index_j, shelposi, shelposj, bfposi, bfposj;
  int shift1, shift2, oppshift1, oppshift2;
  int sheli, shelj;
  int sheli1, shelj1;
  int *p_i, *p_i1, *p_j, *p_j1;
  double *p_rot1, *p_rot2, *p_F, *p_F_rotated;
  double *p_F_rotated1;
  nd1 = atoms->bfnnumb_sh[ip];
  nd2 = atoms->bfnnumb_sh[jp];

    //fprintf(file.out,"before rotate_sum_block\n") ;
    //double F_rotated1[nd1*nd2];
    //for (i=0;i<nd1 * nd2;i++) F_rotated1[i] = k_zero;
    //p_F = F ; for (i=0;i<nd1;i++) { for(j=0;j<nd2;j++) {  fprintf(file.out,"%14.8lf",*p_F); p_F++; } 
    //fprintf(file.out,"\n"); }
    //fprintf(file.out,"\n");


    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
        //CHANGES2015oppshift1 = *(symmetry->op_shift + op * 5 + shells->ord[index_i]);
        //oppshift2 = *(symmetry->op_shift + op * 5 + shells->ord[index_j]);
        //sheli1 = *(symmetry->num_ij + op * 5 + shells->ord[index_i]);
        //shelj1 = *(symmetry->num_ij + op * 5 + shells->ord[index_j]);
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
        //fprintf(file.out,"oppshift %3d %3d %3d %3d %3d\n",op,oppshift1,oppshift2,sheli1,shelj1) ; fflush(file.out);
        p_i1 = symmetry->ind_i + oppshift1;
        p_i = symmetry->ind_j + oppshift1;
        p_rot1 = symmetry->rot + oppshift1;
        for (i = 0; i < sheli1; i++) {
          p_j1 = symmetry->ind_i + oppshift2;
          p_j = symmetry->ind_j + oppshift2;
          p_rot2 = symmetry->rot + oppshift2;
          for (j = 0; j < shelj1; j++) {

            //p_F_reduced = F_reduced + (bfposi + *p_i)  * nd2 + bfposj + *p_j;
            //p_F_rotate  = F_rotate  + (bfposi + *p_i1) * nd2 + bfposj + *p_j1;   // no transpose

            p_F         = F         + (bfposi + *p_i)  * nd2 + bfposj + *p_j;
            p_F_rotated = F_rotated + (bfposi + *p_i1) * nd2 + bfposj + *p_j1;
            //p_F_rotated1= F_rotated1+ (bfposi + *p_i1) * nd2 + bfposj + *p_j1;
            //fprintf(file.out,"offset %3d  i j  %3d %3d *p_i,j  %3d %3d bfposi   %3d %3d  %10.4lf %10.4lf %10.4lf\n", \
            (bfposi+*p_i )*nd2  + bfposj + *p_j,i,j,*p_i,*p_j,bfposi,bfposj,*p_rot1,*p_rot2,*p_F) ;
            //fprintf(file.out,"offset1%3d \n", (bfposi+*p_i1)*nd2  + bfposj + *p_j1) ;
            *p_F_rotated += *p_rot1 * *p_rot2 * *p_F;
            //*p_F_rotated1+= *p_rot1 * *p_rot2 * *p_F;
             //fprintf(file.out,"%18.6e\n",*p_rot1 * *p_rot2 * *p_F);
            //F_rotated1[(bfposi + *p_i)  * nd2 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_F;
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj += shells->type1[index_j];
      } // close loop over index_j
      bfposi += shells->type1[index_i];
    } // close loop over index_i

    //fprintf(file.out,"ip %3d jp %3d op %3d\n",ip,jp,op);
    //for (i=0;i<nd1;i++) { for(j=0;j<nd2;j++) {  fprintf(file.out,"%18.10lf",F_rotated[i*nd2+j]); } 
    //fprintf(file.out,"\n"); }
    //fprintf(file.out,"\n");

    //p_F_rotated = F_rotated;
    //p_F = F; // + nd1 * nd2 * k; // Put rotated pairs in consecutive order. Place in proper order only in FT'd arrays

    //fprintf(file.out,"SYM 2E ROTATED ip %3d jp %3d op %3d \n",ip,jp,op) ;
    //p_F = F_rotated;
    //for (i = 0; i < nd1; i++) {
      //for (j = 0; j < nd2; j++) {
        //fprintf(file.out,"%6.3lf ",*p_F) ;
        //p_F++;
       //}
      //fprintf(file.out,"\n");
     //}

    //for (i = 0; i < nd1; i++) {
      //for (j = 0; j < nd2; j++) {
        //*p_F = *p_F_rotate;
        //fprintf(file.out,"p latt2 i j *p_rotate%3d %3d %3d %3d %10.4lf\n",p,pair_p->latt2[k],i,j,*p_F_rotate) ;
        //p_F++;
        //p_F_rotated++;
      //}
    //}

}

void rotate_psi(Complex *F, Complex *F_rotate, int nbands, int k, KPOINT_TRAN *knet, ATOM_TRAN *atom_p, ATOM *atoms, REAL_LATTICE *R, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)
//CHANGE2014void rotate_psi(Complex *F, Complex *F_rotate, int bands[2], int k, KPOINT_TRAN *knet, ATOM_TRAN *atom_p, ATOM *atoms, REAL_LATTICE *R, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int i, ip, j, n, dim;
  int op, opinv, O;
//CHANGE2014  int nbands = bands[1] - bands[0] + 1;
  int index_i, shelposi, bfposi;
  int bfposip;
  int shift1, oppshift1;
  int sheli;
  int sheli1;
  int *p_i, *p_i1;
  double temp;
  double *p_rot1;
  Complex temp1;
  Complex *p_F, *p_F_rotate;

  // rotate and complex conjugate (if necessary) a block of wavefunctions at one k point
 
  dim = nbands * atoms->number_of_sh_bfns_in_unit_cell;

    p_F_rotate = F_rotate;
    for (i = 0; i < dim; i++) {
      *p_F_rotate = k_zero;
        p_F_rotate++;
       }

    op    = symmetry->inverse[knet->opr[k]];
    opinv = knet->opr[k];
    //fprintf(file.out,"rotate_psi k %3d op %3d opinv %3d\n",k,op,opinv);
    for (n = 0; n < nbands; n++) {
    for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
    ip = atom_p->K[i * symmetry->number_of_operators + op];
    O  = atom_p->O[i * symmetry->number_of_operators + op];
    temp  = double_vec_dot(&knet->cart[k],&R->vec_ai[O]) ;
    if (knet->trs[k] == 1) temp *=-k_one;
    temp1 = Complex(cos(temp),sin(temp)) ;
    shelposi = atoms->shelposn[i];
    bfposi  = atoms->bfnposn_sh[i];
    bfposip = atoms->bfnposn_sh[ip];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[i]; index_i++) {
      //CHANGES2015oppshift1 = *(symmetry->op_shift + opinv * 5 + shells->ord[index_i]);
      //sheli1    = *(symmetry->num_ij   + opinv * 5 + shells->ord[index_i]);
      oppshift1 = *(symmetry->op_shift + opinv * (job->l_max + 2) + shells->ord[index_i]);
      sheli1    = *(symmetry->num_ij   + opinv * (job->l_max + 2) + shells->ord[index_i]);
      p_i1 = symmetry->ind_i + oppshift1;
      p_i  = symmetry->ind_j + oppshift1;
      p_rot1 = symmetry->rot + oppshift1;
      for (j = 0; j < sheli1; j++) {
        p_F        = F        + n * atoms->number_of_sh_bfns_in_unit_cell + bfposi + *p_i1;
        p_F_rotate = F_rotate + n * atoms->number_of_sh_bfns_in_unit_cell + bfposip + *p_i;
       *p_F_rotate += *p_rot1 * *p_F * temp1;
        p_i++;
        p_i1++;
        p_rot1++;
       }
      bfposi  += shells->type1[index_i];
      bfposip += shells->type1[index_i];
     } // close loop over index_i
    }
   }

    if (knet->trs[k] == 1) {
    p_F_rotate = F_rotate;
    for (i = 0; i < nbands * atoms->number_of_sh_bfns_in_unit_cell; i++) {
    *p_F_rotate *= Complex (-k_one, k_zero);
    //(*p_F_rotate).imag() *= -k_one;
    p_F_rotate++;
   }
  }
//printf("psi %d %d %d %10.4lf\n",dim,nbands,k,(*F_rotate).real());

}

