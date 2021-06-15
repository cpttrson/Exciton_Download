#include "myconstants.h"
#include "USER_DATA.h"
#include "CARTESIAN_TO_SH.h"

using namespace std;

void cartesian_to_sh_ij(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd2, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)
//void two_centre_cartesian_to_sh_ij(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd2, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)
//void two_centre_cartesian_to_sh_shell_ij(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd2, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Transform one-body operator over atom pair from Cartesian to solid harmonic basis      *
  // ******************************************************************************************

int i, j;
int nsheli, nshelj;
int *p_i, *p_j, *p_i1, *p_j1;
double *p_rot1, *p_rot2;
double *p_Operator_cart;

  nsheli = *(shells->num_ij + shells->ord_sh[index_i]);
  nshelj = *(shells->num_ij + shells->ord_sh[index_j]);
  p_i1 = shells->ind_i + shells->opp_sh[index_i];
  p_i  = shells->ind_j + shells->opp_sh[index_i];
  p_rot1 = shells->rot + shells->opp_sh[index_i];
  for (i = 0; i < nsheli; i++) {
    p_j1 = shells->ind_i + shells->opp_sh[index_j];
    p_j  = shells->ind_j + shells->opp_sh[index_j];
    p_rot2 = shells->rot + shells->opp_sh[index_j];
    for (j = 0; j < nshelj; j++) {
      p_Operator_cart = Operator_cart + (bfposi1 + *p_i1) * nd2 + bfposj1 + *p_j1;
      Operator_sh[(bfposi + *p_i) * nd4 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_Operator_cart;
      p_j++;
      p_j1++;
      p_rot2++;
     }
    p_i++;
    p_i1++;
    p_rot1++;
   }

}

void cartesian_to_sh_ij_vector(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd1, int nd2, int nd3, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)
//void two_centre_cartesian_to_sh_ij_vector(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd1, int nd2, int nd3, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)
//void two_centre_vector_cartesian_to_sh_shell_ij(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd1, int nd2, int nd3, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)

{

  // **************************************************************************************************
  // * Transform three-vector one-body operator over atom pair from Cartesian to solid harmonic basis *
  // **************************************************************************************************

int i, j, k;
int nsheli, nshelj;
int *p_i, *p_j, *p_i1, *p_j1;
double *p_rot1, *p_rot2;
double *p_Operator_cart;

  for (k = 0; k < 3; k++) {
    nsheli = *(shells->num_ij + shells->ord_sh[index_i]);
    nshelj = *(shells->num_ij + shells->ord_sh[index_j]);
    p_i1 = shells->ind_i + shells->opp_sh[index_i];
    p_i  = shells->ind_j + shells->opp_sh[index_i];
    p_rot1 = shells->rot + shells->opp_sh[index_i];
    for (i = 0; i < nsheli; i++) {
      p_j1 = shells->ind_i + shells->opp_sh[index_j];
      p_j  = shells->ind_j + shells->opp_sh[index_j];
      p_rot2 = shells->rot + shells->opp_sh[index_j];
      for (j = 0; j < nshelj; j++) {
        p_Operator_cart = Operator_cart + k * nd1 * nd2 + (bfposi1 + *p_i1) * nd2 + bfposj1 + *p_j1;
          Operator_sh[k * nd3 * nd4 + (bfposi + *p_i) * nd4 + bfposj + *p_j] += *p_rot1 * *p_rot2 * *p_Operator_cart;
        p_j++;
        p_j1++;
        p_rot2++;
       }
      p_i++;
      p_i1++;
      p_rot1++;
     }
    }

}

void cartesian_to_sh_ija(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int index_k, int bfposi, int bfposj, int bfposk, int bfposi1, int bfposj1, int bfposk1, int nd2, int nd3, int nd5, int nd6, SHELL *shells, SHELL *shells_ax, JOB_PARAM *job, FILES file)
//void three_centre_cartesian_to_sh_ija(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int index_k, int bfposi, int bfposj, int bfposk, int bfposi1, int bfposj1, int bfposk1, int nd2, int nd3, int nd5, int nd6, SHELL *shells, SHELL *shells_ax, JOB_PARAM *job, FILES file)
//void three_centre_cartesian_to_sh_shell_ax_reversed(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int index_k, int bfposi, int bfposj, int bfposk, int bfposi1, int bfposj1, int bfposk1, int nd2, int nd3, int nd5, int nd6, SHELL *shells, SHELL *shells_ax, JOB_PARAM *job, FILES file)

{

  // **************************************************************************************************
  // * Transform three-centre two-body operator (ij|alpha) from Cartesian to solid harmonic basis     *
  // **************************************************************************************************

int i, j, k;
int nsheli, nshelj, nshelk;
int *p_i, *p_j, *p_i1, *p_j1, *p_k, *p_k1;
double *p_rot1, *p_rot2, *p_rot3;
double *p_Operator_cart;

  nsheli = *(shells->num_ij + shells->ord_sh[index_i]);
  nshelj = *(shells->num_ij + shells->ord_sh[index_j]);
  nshelk = *(shells->num_ij + shells_ax->ord_sh[index_k]);
  p_i1 = shells->ind_i + shells->opp_sh[index_i];
  p_i  = shells->ind_j + shells->opp_sh[index_i];
  p_rot1 = shells->rot + shells->opp_sh[index_i];
  for (i = 0; i < nsheli; i++) {
    p_j1 = shells->ind_i + shells->opp_sh[index_j];
    p_j  = shells->ind_j + shells->opp_sh[index_j];
    p_rot2 = shells->rot + shells->opp_sh[index_j];
    for (j = 0; j < nshelj; j++) {
      p_k1 = shells->ind_i + shells_ax->opp_sh[index_k];
      p_k  = shells->ind_j + shells_ax->opp_sh[index_k];
      p_rot3 = shells->rot + shells_ax->opp_sh[index_k];
      for (k = 0; k < nshelk; k++) {
      p_Operator_cart = Operator_cart + (bfposi1 + *p_i1) * nd2 * nd3 + (bfposj1 + *p_j1) * nd3 + bfposk1 + *p_k1;
      Operator_sh[(bfposi + *p_i) * nd5 * nd6 + (bfposj + *p_j) * nd6 +  bfposk + *p_k] += *p_rot1 * *p_rot2 * *p_rot3 * *p_Operator_cart;
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

}

void three_centre_cartesian_to_sh_ija_complex(Complex *Operator_cart, Complex *Operator_sh, int index_i, int index_j, int index_k, int bfposi, int bfposj, int bfposk, int bfposi1, int bfposj1, int bfposk1, int nd2, int nd3, int nd5, int nd6, SHELL *shells, SHELL *shells_ax, JOB_PARAM *job, FILES file)
//void three_center_cartesian_to_sh_shell_ax_reversed_complex(Complex *Operator_cart, Complex *Operator_sh, int index_i, int index_j, int index_k, int bfposi, int bfposj, int bfposk, int bfposi1, int bfposj1, int bfposk1, int nd2, int nd3, int nd5, int nd6, SHELL *shells, SHELL *shells_ax, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************************
  // * Transform complex three-centre two-body operator (ij|alpha) from Cartesian to solid harmonic basis *
  // ******************************************************************************************************

int i, j, k;
int nsheli, nshelj, nshelk;
int *p_i, *p_j, *p_i1, *p_j1, *p_k, *p_k1;
double *p_rot1, *p_rot2, *p_rot3;
Complex *p_Operator_cart;

  nsheli = *(shells->num_ij + shells->ord_sh[index_i]);
  nshelj = *(shells->num_ij + shells->ord_sh[index_j]);
  nshelk = *(shells->num_ij + shells_ax->ord_sh[index_k]);
  p_i1 = shells->ind_i + shells->opp_sh[index_i];
  p_i  = shells->ind_j + shells->opp_sh[index_i];
  p_rot1 = shells->rot + shells->opp_sh[index_i];
  for (i = 0; i < nsheli; i++) {
    p_j1 = shells->ind_i + shells->opp_sh[index_j];
    p_j  = shells->ind_j + shells->opp_sh[index_j];
    p_rot2 = shells->rot + shells->opp_sh[index_j];
    for (j = 0; j < nshelj; j++) {
      p_k1 = shells->ind_i + shells_ax->opp_sh[index_k];
      p_k  = shells->ind_j + shells_ax->opp_sh[index_k];
      p_rot3 = shells->rot + shells_ax->opp_sh[index_k];
      for (k = 0; k < nshelk; k++) {
      p_Operator_cart = Operator_cart + (bfposi1 + *p_i1) * nd2 * nd3 + (bfposj1 + *p_j1) * nd3 + bfposk1 + *p_k1;
      Operator_sh[(bfposi + *p_i) * nd5 * nd6 + (bfposj + *p_j) * nd6 +  bfposk + *p_k] += *p_rot1 * *p_rot2 * *p_rot3 * *p_Operator_cart;
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

}

void cartesian_to_sh_ijkl(double *F_cart, double *F_sh, int index_i, int index_j, int index_k, int index_l, SHELL *shells, JOB_PARAM *job, FILES file)
//void four_centre_cartesian_to_sh_ijkl(double *F_cart, double *F_sh, int index_i, int index_j, int index_k, int index_l, SHELL *shells, JOB_PARAM *job, FILES file)
//void four_centre_cartesian_to_sh_atom_ijkl(double *F_cart, double *F_sh, int index_i, int index_j, int index_k, int index_l, SHELL *shells, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Transform four-centre two-body operator (ij|kl) from Cartesian to solid harmonic basis *
  // ******************************************************************************************

int off1, off2, off3, off4, off5, off6;
int nsheli, nshelj, nshelk, nshell;
int *p_i, *p_j, *p_k, *p_l, *p_i1, *p_j1, *p_k1, *p_l1;
int i, j, k, l;
double *p_rot1, *p_rot2, *p_rot3, *p_rot4;
double *p_F_cart, *p_F_sh;

  off3  = shells->type1_sh[index_l];
  off2  = off3 * shells->type1_sh[index_k];
  off1  = off2 * shells->type1_sh[index_j];
  off6  = shells->type_sh[index_l];
  off5  = off6 * shells->type_sh[index_k];
  off4  = off5 * shells->type_sh[index_j];
  nsheli = *(shells->num_ij + shells->ord_sh[index_i]);
  nshelj = *(shells->num_ij + shells->ord_sh[index_j]);
  nshelk = *(shells->num_ij + shells->ord_sh[index_k]);
  nshell = *(shells->num_ij + shells->ord_sh[index_l]);
  p_i1 = shells->ind_i + shells->opp_sh[index_i];
  p_i  = shells->ind_j + shells->opp_sh[index_i];
  p_rot1 = shells->rot + shells->opp_sh[index_i];
  for (i = 0; i < nsheli; i++) {
    p_j1 = shells->ind_i + shells->opp_sh[index_j];
    p_j  = shells->ind_j + shells->opp_sh[index_j];
    p_rot2 = shells->rot + shells->opp_sh[index_j];
    for (j = 0; j < nshelj; j++) {
      p_k1 = shells->ind_i + shells->opp_sh[index_k];
      p_k  = shells->ind_j + shells->opp_sh[index_k];
      p_rot3 = shells->rot + shells->opp_sh[index_k];
      for (k = 0; k < nshelk; k++) {
        p_l1 = shells->ind_i + shells->opp_sh[index_l];
        p_l  = shells->ind_j + shells->opp_sh[index_l];
        p_rot4 = shells->rot + shells->opp_sh[index_l];
        for (l = 0; l < nshell; l++) {
          p_F_cart  = F_cart + *p_i1 * off1 + *p_j1 * off2 + *p_k1 * off3 + *p_l1;
          p_F_sh    = F_sh   + *p_i  * off4 + *p_j  * off5 + *p_k  * off6 + *p_l;
         *p_F_sh   += *p_F_cart * *p_rot1 * *p_rot2 * *p_rot3 * *p_rot4;
          p_l++;
          p_l1++;
          p_rot4++;
         }
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

}

void cartesian_to_sh_ijij(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd2, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)
//void four_centre_cartesian_to_sh_ijij(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd2, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)
//void four_centre_cartesian_to_sh_atom_ijkl_screen(double *Operator_cart, double *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd2, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Transform four-centre two-body operator (ij|ij) from Cartesian to solid harmonic basis *
  // ******************************************************************************************

int i, j, k, l;
int nsheli, nshelj, nshelk, nshell;
int *p_i, *p_j, *p_k, *p_l, *p_i1, *p_j1, *p_k1, *p_l1;
int off1, off2, off3;
double *p_rot1, *p_rot2, *p_rot3, *p_rot4;
double *p_Operator_cart;

  off3  = shells->type1_sh[index_j];
  off2  = off3 * shells->type1_sh[index_i];
  off1  = off2 * shells->type1_sh[index_j];
  nsheli = *(shells->num_ij + shells->ord_sh[index_i]);
  nshelj = *(shells->num_ij + shells->ord_sh[index_j]);
  p_i1 = shells->ind_i + shells->opp_sh[index_i];
  p_i  = shells->ind_j + shells->opp_sh[index_i];
  p_rot1 = shells->rot + shells->opp_sh[index_i];
  for (i = 0; i < nsheli; i++) {
    p_j1 = shells->ind_i + shells->opp_sh[index_j];
    p_j  = shells->ind_j + shells->opp_sh[index_j];
    p_rot2 = shells->rot + shells->opp_sh[index_j];
    for (j = 0; j < nshelj; j++) {
      p_k1 = shells->ind_i + shells->opp_sh[index_i];
      p_k  = shells->ind_j + shells->opp_sh[index_i];
      p_rot3 = shells->rot + shells->opp_sh[index_i];
      for (k = 0; k < nsheli; k++) {
        p_l1 = shells->ind_i + shells->opp_sh[index_j];
        p_l  = shells->ind_j + shells->opp_sh[index_j];
        p_rot4 = shells->rot + shells->opp_sh[index_j];
        for (l = 0; l < nshelj; l++) {
          if (*p_k == *p_i && *p_l == *p_j) {
          p_Operator_cart  = Operator_cart + *p_i1 * off1 + *p_j1 * off2 + *p_k1 * off3 + *p_l1;
          Operator_sh[(bfposi1 + *p_i) * nd4 + bfposj1 + *p_j] += *p_rot1 * *p_rot2 * *p_rot3 * *p_rot4 * *p_Operator_cart;
         }
        p_l++;
        p_l1++;
        p_rot4++;
       }
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

}

void four_centre_cartesian_to_sh_ijij_complex(Complex *Operator_cart, Complex *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd2, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)
//void two_centre_cartesian_to_sh_shell_ij_complex(Complex *Operator_cart, Complex *Operator_sh, int index_i, int index_j, int bfposi, int bfposj, int bfposi1, int bfposj1, int nd2, int nd4, SHELL *shells, JOB_PARAM *job, FILES file)

{

  // **************************************************************************************************
  // * Transform complex four-centre two-body operator (ij|ij) from Cartesian to solid harmonic basis *
  // **************************************************************************************************

int i, j, k, l;
int nsheli, nshelj, nshelk, nshell;
int *p_i, *p_j, *p_k, *p_l, *p_i1, *p_j1, *p_k1, *p_l1;
int off1, off2, off3;
double *p_rot1, *p_rot2, *p_rot3, *p_rot4;
Complex *p_Operator_cart;

  off3  = shells->type1_sh[index_j];
  off2  = off3 * shells->type1_sh[index_i];
  off1  = off2 * shells->type1_sh[index_j];
  nsheli = *(shells->num_ij + shells->ord_sh[index_i]);
  nshelj = *(shells->num_ij + shells->ord_sh[index_j]);
  p_i1 = shells->ind_i + shells->opp_sh[index_i];
  p_i  = shells->ind_j + shells->opp_sh[index_i];
  p_rot1 = shells->rot + shells->opp_sh[index_i];
  for (i = 0; i < nsheli; i++) {
    p_j1 = shells->ind_i + shells->opp_sh[index_j];
    p_j  = shells->ind_j + shells->opp_sh[index_j];
    p_rot2 = shells->rot + shells->opp_sh[index_j];
    for (j = 0; j < nshelj; j++) {
      p_k1 = shells->ind_i + shells->opp_sh[index_i];
      p_k  = shells->ind_j + shells->opp_sh[index_i];
      p_rot3 = shells->rot + shells->opp_sh[index_i];
      for (k = 0; k < nsheli; k++) {
        p_l1 = shells->ind_i + shells->opp_sh[index_j];
        p_l  = shells->ind_j + shells->opp_sh[index_j];
        p_rot4 = shells->rot + shells->opp_sh[index_j];
        for (l = 0; l < nshelj; l++) {
          if (*p_k == *p_i && *p_l == *p_j) {
          p_Operator_cart  = Operator_cart + *p_i1 * off1 + *p_j1 * off2 + *p_k1 * off3 + *p_l1;
          Operator_sh[(bfposi1 + *p_i) * nd4 + bfposj1 + *p_j] += *p_rot1 * *p_rot2 * *p_rot3 * *p_rot4 * *p_Operator_cart;
         }
        p_l++;
        p_l1++;
        p_rot4++;
       }
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

}

