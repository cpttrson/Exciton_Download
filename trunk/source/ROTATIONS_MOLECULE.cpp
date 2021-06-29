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
  for (i = 0; i < nd1; i++) F_rotate[i] = Complex(k_zero, k_zero);

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
      p_i++;
      p_i1++;
      p_rot1++;
     }
    bfposi += shells->type1[index_i];
   } // close loop over index_i

}

void rotate_permute_expand_pair(int p, PAIR_TRAN *pair_p, double *F_reduced, double *F, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ********************************************************************************************************
  // Generates a full real operator <i|j> from a reduced operator using permutation and rotation operations *
  // ********************************************************************************************************

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

    for (i = 0; i < nd1 * nd2; i++) F_rotate[i] = k_zero;
 
      switch (pm) {

        case 0:

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
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

    // Put rotated pairs in consecutive order. Place in proper order only in FT'd arrays
    for (i = 0; i < nd1 *nd2; i++) F[nd1 * nd2 * k + i] = F_rotate[i];

  } // close loop on k

  free(F_rotate);

}

void rotate_permute_expand_pair_complex(int p, PAIR_TRAN *pair_p, Complex *F_reduced, Complex *F, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ***********************************************************************************************************
  // Generates a full complex operator <i|j> from a reduced operator using permutation and rotation operations *
  // ***********************************************************************************************************

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

    for (i = 0; i < nd1 * nd2; i++) F_rotate[i] = Complex(k_zero, k_zero);
    //p_F_rotate = F_rotate;
    //for (i = 0; i < nd1 * nd2; i++) {
      //*p_F_rotate = k_zero;
        //p_F_rotate++;
       //}
 
      switch (pm) {

        case 0:

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
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

    // Put rotated pairs in consecutive order. Place in proper order only in FT'd arrays
    for (i = 0; i < nd1 *nd2; i++) F[nd1 * nd2 * k + i] = F_rotate[i];

    //p_F_rotate = F_rotate;
    //p_F = F + nd1 * nd2 * k; // Put rotated pairs in consecutive order. Place in proper order only in FT'd arrays
    //fprintf(file.out,"p %3d ip %3d jp %3d gj %3d op %3d pm %3d\n",p,ip,jp,gj,op,pm) ;
    //for (i = 0; i < nd1; i++) {
      //for (j = 0; j < nd2; j++) {
        //*p_F = *p_F_rotate;
        //fprintf(file.out,"%10.4lf ",*p_F) ;
        //p_F++;
        //p_F_rotate++;
      //}
     //fprintf(file.out,"\n");
     //}

  } // close loop on k

  free(F_rotate);

}

void rotate_permute_expand_tensor1(int p, PAIR_TRAN *pair_p, double *F_reduced, double *F, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ********************************************************************************************************
  // Generates a full operator real <i|j> from a reduced operator using permutation and rotation operations *
  // ********************************************************************************************************

int i, j, k, t, l, nd1, nd2 ;
int ip, jp, gj, op, pm ;
int index_i, index_j, shelposi, shelposj, bfposi, bfposj ;
int shift1, shift2, oppshift1, oppshift2, oppshift3 ;
int sheli, shelj ;
int sheli1, shelj1, shelk1 ;
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

  for (i = 0; i < 3 * nd1 * nd2; i++) F_rotate[i] = k_zero;
      //p_F_rotate = F_rotate;
      //for (i = 0; i < 3 * nd1 * nd2; i++) {
        //*p_F_rotate = k_zero;
         //p_F_rotate++;
        //}

    switch (pm) {

        case 0:

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
      shelposj = atoms->shelposn[jp];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
        oppshift3 = *(symmetry->op_shift + op * (job->l_max + 2) + 2);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
        shelk1 = *(symmetry->num_ij + op * (job->l_max + 2) + 2);
        p_k1   = symmetry->ind_i + oppshift3 ;
        p_k    = symmetry->ind_j + oppshift3 ;
        p_rot3 = symmetry->rot   + oppshift3 ;
        for (t = 0; t < shelk1; t++) {
        //for (t = 0; t < 3; t++) {
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
	      //fprintf(file.out,"op %3d %3d %3d off %3d %3d shift %3d %6.2f %6.2f %6.2f %6.2f %6.2f\n",\
	      t,i,j,*p_k,*p_k1,oppshift3,*p_rot1,*p_rot2,*p_rot3,F_reduced[*p_k  * nd1 * nd2 + (bfposi + *p_i)  * nd2 + bfposj + *p_j],\
	      F_rotate[*p_k1 * nd1 * nd2 + (bfposi + *p_i1) * nd2 + bfposj + *p_j1]);
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
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
        oppshift3 = *(symmetry->op_shift + op * (job->l_max + 2) + 2);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
        shelk1 = *(symmetry->num_ij + op * (job->l_max + 2) + 2);
        p_k1   = symmetry->ind_i + oppshift3 ;
        p_k    = symmetry->ind_j + oppshift3 ;
        p_rot3 = symmetry->rot   + oppshift3 ;
        for (t = 0; t < shelk1; t++) {
        //for (t = 0; t < 3; t++) {
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

    // Put rotated pairs in consecutive order. Place in proper order only in FT'd arrays
    for (i = 0; i < 3 * nd1 *nd2; i++) F[3 * nd1 * nd2 * k + i] = F_rotate[i];

   //p_F_rotate = F_rotate;
   //p_F = F + 3 * nd1 * nd2 * k; // Put rotated pairs in consecutive order. Place in proper order only in FT'd arrays
   //for (t = 0; t < 3; t++) {
    //for (i = 0; i < nd1; i++) {
      //for (j = 0; j < nd2; j++) {
        //*p_F = *p_F_rotate;
        //fprintf(file.out,"p latt2 i j *p_rotate%3d %3d %3d %3d %12.6lf\n",p,pair_p->latt2[k],i,j,*p_F_rotate) ;
        //p_F++;
        //p_F_rotate++;
      //}
      //fprintf(file.out,"\n") ;
    //} // close loop on t
    //}

  } // close loop on k

  free(F_rotate);

}

void rotate_permute_triple_ija(int *ip, int *jp, int *kp, int *op, int *pm, double *F_reduced, double *F, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, ATOM *atoms_ax, SHELL *shells_ax, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)
//void rotate_permute_triple_ij_alpha(int *ip, int *jp, int *alpha, int *op, int *pm, double *F_reduced, double *F, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, ATOM *atoms_ax, SHELL *shells_ax, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ****************************************************************************************************
  // Generates a full operator <ij|a> from a reduced operator using permutation and rotation operations *
  // ****************************************************************************************************

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
  nd3 = atoms_ax->bfnnumb_sh[*kp];

  for (i = 0; i < nd1 * nd2 * nd3; i++) F[i] = k_zero;

    switch (*pm) {

      case 0: // no permutation transpose

  shelposi = atoms->shelposn[*ip];
  bfposi = 0;
  for (index_i = shelposi; index_i < shelposi + atoms->nshel[*ip]; index_i++) {
    shelposj = atoms->shelposn[*jp];
    bfposj = 0;
    for (index_j = shelposj; index_j < shelposj + atoms->nshel[*jp]; index_j++) {
      shelposk = atoms_ax->shelposn[*kp];
      bfposk = 0;
      for (index_k = shelposk; index_k < shelposk + atoms_ax->nshel[*kp]; index_k++) {
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
      shelposk = atoms_ax->shelposn[*kp];
      bfposk = 0;
      for (index_k = shelposk; index_k < shelposk + atoms_ax->nshel[*kp]; index_k++) {
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

  // ************************************************************************************************************************
  // Generates a Fock operator <ij|kl>P_kl or <ij|kl>P_jl from a reduced operator using permutation and rotation operations *
  // ************************************************************************************************************************

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
  shelposi = atoms->shelposn[ip];
  bfposi = 0;
  for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    shelposj = atoms->shelposn[jp];
    bfposj = 0;
    for (index_j = shelposj; index_j < shelposj + atoms->nshel[jp]; index_j++) {
      oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_i]);
      oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord[index_j]);
      sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_i]);
      shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + shells->ord[index_j]);
      p_i1 = symmetry->ind_i + oppshift1;
      p_i = symmetry->ind_j + oppshift1;
      p_rot1 = symmetry->rot + oppshift1;
      for (i = 0; i < sheli1; i++) {
        p_j1 = symmetry->ind_i + oppshift2;
        p_j = symmetry->ind_j + oppshift2;
        p_rot2 = symmetry->rot + oppshift2;
        for (j = 0; j < shelj1; j++) {
          p_F         = F         + (bfposi + *p_i)  * nd2 + bfposj + *p_j;
          p_F_rotated = F_rotated + (bfposi + *p_i1) * nd2 + bfposj + *p_j1;
          *p_F_rotated += *p_rot1 * *p_rot2 * *p_F;
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

}

void rotate_psi(Complex *F, Complex *F_rotate, int nbands, int k, KPOINT_TRAN *knet, ATOM_TRAN *atom_p, ATOM *atoms, REAL_LATTICE *R, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, ip, j, n, dim;
int op, opinv, O;
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

  // ************************************************************************************************************************
  // Rotate and complex conjugate (if necessary) a block of wavefunctions at one k point                                    *
  // ************************************************************************************************************************

  dim = nbands * atoms->number_of_sh_bfns_in_unit_cell;
  for (i = 0; i < dim; i++) F_rotate[i] = Complex(k_zero, k_zero);
  op    = symmetry->inverse[knet->opr[k]];
  opinv = knet->opr[k];
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

      if (knet->trs[k] == 1) { for (i = 0; i < dim; i++) F_rotate[i] *= Complex(-k_one, k_one);  }

}

