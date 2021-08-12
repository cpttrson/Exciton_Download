

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
#include "mycomplex.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "PRINT_UTIL.h"
#include "ROTATIONS_MOLECULE.h"
#include "SYMMETRY_ADAPTATION.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"

using namespace std;

void count_atom_irrep(IntMatrix *num_irrep_in_atom, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, s, op;
int i1, i2, i3, i4;
int oppshift1;
int atm, atm_rot, index_i, shelposi, ip, ipp;
int *p_ind_i, *p_ind_j, *p_num;
double *p_rot, dot_product;
double character[symmetry->number_of_classes][atoms->number_of_unique_atoms];

  ResetIntMatrix(num_irrep_in_atom);

  for (s = 0; s < symmetry->number_of_classes; s++) {
    for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
      character[s][atm] = k_zero;
     }
    }

  for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
    ip = atom_p->posn[atm];
    shelposi = atoms->shelposn_sh[ip];
    for (ipp = 0; ipp < atom_p->numb[atm]; ipp++) {
      for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
        for (i1 = 0; i1 < shells->shar[index_i]; i1++) {
          op = 0;
          for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
            for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
              oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
              p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
              atm_rot = atom_p->K[(ip + ipp) * symmetry->number_of_operators + op];
              p_ind_i = symmetry->ind_i  + oppshift1;
              p_ind_j = symmetry->ind_j  + oppshift1;
              p_rot   = symmetry->rot    + oppshift1;
              for (i4 = 0; i4 < *p_num; i4++) {
                if (i1 == *p_ind_i && i1 == *p_ind_j && atm_rot == ip + ipp) character[i2][atm] += *p_rot;
                p_ind_i++;
                p_ind_j++;
                p_rot++;
               }
              op++;
             } // close loop on i3
            } // close loop on i2
           } // close loop on i1
          } // close index_i
         } // close loop on ipp
        } // close loop on atm
 
      for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
        for (i = 0; i < symmetry->number_of_classes; i++) {
          dot_product = k_zero;
          for (j = 0; j < symmetry->number_of_classes; j++) {
            dot_product += symmetry->character_table[i * symmetry->number_of_classes + j] * character[j][atm];
           }
            num_irrep_in_atom->a[i][atm] = (int) ((dot_product + 0.000001) / (double) symmetry->grp_dim);
            //fprintf(file.out,"ATOM  IRREP %3d %3d %3d %10.4lf\n",\
            i,atm,num_irrep_in_atom->a[i][atm],dot_product / (double) symmetry->grp_dim);
          }
         }
         
}

/*
void count_atom_irrep_crystal(IntMatrix *num_irrep_in_atom, ATOM *atoms, ATOM_TRAN *atom_p, REAL_LATTICE_TABLES *R_tables, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, s, op;
int i1, i2, i3, i4;
int oppshift1;
int atm, atm_rot, index_i, shelposi, ip, ipp;
int *p_ind_i, *p_ind_j, *p_num;
double *p_rot, dot_product;
double character[symmetry->number_of_classes][atoms->number_of_unique_atoms];

int latt = 1, latt_rot, O_temp, latt_temp;

  ResetIntMatrix(num_irrep_in_atom);

  for (s = 0; s < symmetry->number_of_classes; s++) {
    for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
      character[s][atm] = k_zero;
     }
    }

for (latt = 1; latt < 13; latt++) {

  for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
    ip = atom_p->posn[atm];
    shelposi = atoms->shelposn_sh[ip];
    for (ipp = 0; ipp < atom_p->numb[atm]; ipp++) {
      for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
        for (i1 = 0; i1 < shells->shar[index_i]; i1++) {
          op = 0;
          for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
            for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
              oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
              p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
              atm_rot = atom_p->K[(ip + ipp) * symmetry->number_of_operators + op];
              O_temp = atom_p->O[(ip + ipp) * symmetry->number_of_operators + op];
              latt_temp = R_tables->diffvec[R_tables->lattvec[latt * symmetry->number_of_operators + op] * \
              R_tables->margin_vector + O_temp];
              latt_rot = latt_temp;
              p_ind_i = symmetry->ind_i  + oppshift1;
              p_ind_j = symmetry->ind_j  + oppshift1;
              p_rot   = symmetry->rot    + oppshift1;
              for (i4 = 0; i4 < *p_num; i4++) {
                //if (i1 == *p_ind_i && i1 == *p_ind_j && atm_rot == ip + ipp) character[i2][atm] += *p_rot;
                if (i1 == *p_ind_i && i1 == *p_ind_j && atm_rot == ip + ipp && latt_rot == latt) character[i2][atm] += *p_rot;
                p_ind_i++;
                p_ind_j++;
                p_rot++;
               }
              op++;
             } // close loop on i3
            } // close loop on i2
           } // close loop on i1
          } // close index_i
         } // close loop on ipp
        } // close loop on atm

       } // close loop on latt
 
      for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
        for (i = 0; i < symmetry->number_of_classes; i++) {
          dot_product = k_zero;
          for (j = 0; j < symmetry->number_of_classes; j++) {
            dot_product += symmetry->character_table[i * symmetry->number_of_classes + j] * character[j][atm];
           }
            num_irrep_in_atom->a[i][atm] = (int) ((dot_product + 0.000001) / (double) symmetry->grp_dim);
            fprintf(file.out,"ATOM  IRREP %3d %3d %3d %10.4lf\n",\
            i,atm,num_irrep_in_atom->a[i][atm],dot_product / (double) symmetry->grp_dim);
          }
         }
            if (job->taskid > 0)
            fprintf(file.out,"\n");
         
}
*/

void count_basis_irrep(int *num_irrep_in_basis, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, s, op;
int i1, i2, i3, i4;
int oppshift1;
int atm, atm_rot, index_i, shelposi, ip, ipp;
int *p_ind_i, *p_ind_j, *p_num;
double *p_rot, dot_product;
double character[symmetry->number_of_classes];

  for (s = 0; s < symmetry->number_of_classes; s++)
  character[s] = k_zero;

  for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
    ip = atom_p->posn[atm];
    shelposi = atoms->shelposn_sh[ip];
    for (ipp = 0; ipp < atom_p->numb[atm]; ipp++) {
      for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
        for (i1 = 0; i1 < shells->shar[index_i]; i1++) {
          op = 0;
          for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
            for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
              oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
              p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
              atm_rot = atom_p->K[(ip + ipp) * symmetry->number_of_operators + op];
              p_ind_i = symmetry->ind_i  + oppshift1;
              p_ind_j = symmetry->ind_j  + oppshift1;
              p_rot   = symmetry->rot    + oppshift1;
              for (i4 = 0; i4 < *p_num; i4++) {
              //fprintf(file.out,"CHARACTER %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %10.4lf %10.4lf\n", \
              ipp,index_i,i1,i2,i3,i4,*p_ind_i,*p_ind_j,*symmetry->op_shift,oppshift1,*p_rot,character[i2]);
                if (i1 == *p_ind_i && i1 == *p_ind_j && atm_rot == ip + ipp) character[i2] += *p_rot;
                p_ind_i++;
                p_ind_j++;
                p_rot++;
               }
              op++;
             } // close loop on i3
            } // close loop on i2
           } // close loop on i1
          } // close index_i
         } // close loop on ipp
        } // close loop on atm
 
        for (i = 0; i < symmetry->number_of_classes; i++) {
          //fprintf(file.out,"CHARACTER %3d %10.4lf\n",i,character[i]);
          dot_product = k_zero;
          for (j = 0; j < symmetry->number_of_classes; j++) {
            dot_product += symmetry->character_table[i * symmetry->number_of_classes + j] * character[j];
           }
            num_irrep_in_basis[i] = (int) ((dot_product + 0.000001) / (double) symmetry->grp_dim);
            if (job->taskid > 0)
            fprintf(file.out,"BASIS IRREP class %3d num irrep in basis %3d     %10.4lf\n",\
            i,num_irrep_in_basis[i],dot_product / (double) symmetry->grp_dim);
           }
            //if (job->taskid > 0)
            //fprintf(file.out,"\n");
         
}

/*
void count_basis_irrep_crystal(int *num_irrep_in_basis, ATOM *atoms, ATOM_TRAN *atom_p, REAL_LATTICE_TABLES *R_tables, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, s, op;
int i1, i2, i3, i4;
int oppshift1;
int atm, atm_rot, index_i, shelposi, ip, ipp;
int *p_ind_i, *p_ind_j, *p_num;
double *p_rot, dot_product;
double character[symmetry->number_of_classes];

  for (s = 0; s < symmetry->number_of_classes; s++)
  character[s] = k_zero;

int latt = 1, latt_rot, O_temp, latt_temp;

for (latt = 1; latt < 13; latt++) {

  for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
    ip = atom_p->posn[atm];
    shelposi = atoms->shelposn_sh[ip];
    for (ipp = 0; ipp < atom_p->numb[atm]; ipp++) {
      for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
        for (i1 = 0; i1 < shells->shar[index_i]; i1++) {
          op = 0;
          for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
            for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
              oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
              p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
              atm_rot = atom_p->K[(ip + ipp) * symmetry->number_of_operators + op];
              O_temp = atom_p->O[(ip + ipp) * symmetry->number_of_operators + op];
              latt_temp = R_tables->diffvec[R_tables->lattvec[latt * symmetry->number_of_operators + op] * \
              R_tables->margin_vector + O_temp];
              latt_rot = latt_temp;
              //latt_rot = R_tables->diffvec[latt_temp * R_tables->margin_vector + latt_temp];
              //fprintf(file.out,"atm %3d op %3d latt %3d latt_rot %3d\n",atm,op,latt,latt_temp);
              p_ind_i = symmetry->ind_i  + oppshift1;
              p_ind_j = symmetry->ind_j  + oppshift1;
              p_rot   = symmetry->rot    + oppshift1;
              for (i4 = 0; i4 < *p_num; i4++) {
              //fprintf(file.out,"CHARACTER %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %10.4lf %10.4lf\n", \
              ipp,index_i,i1,i2,i3,i4,*p_ind_i,*p_ind_j,*symmetry->op_shift,oppshift1,*p_rot,character[i2]);
                //if (i1 == *p_ind_i && i1 == *p_ind_j && atm_rot == ip + ipp) character[i2] += *p_rot;
                if (i1 == *p_ind_i && i1 == *p_ind_j && atm_rot == ip + ipp && latt_rot == latt) character[i2] += *p_rot;
                p_ind_i++;
                p_ind_j++;
                p_rot++;
               }
              op++;
             } // close loop on i3
            } // close loop on i2
           } // close loop on i1
          } // close index_i
         } // close loop on ipp
        } // close loop on atm
} 
        for (i = 0; i < symmetry->number_of_classes; i++) {
          //fprintf(file.out,"CHARACTER %3d %10.4lf\n",i,character[i]);
          dot_product = k_zero;
          for (j = 0; j < symmetry->number_of_classes; j++) {
            dot_product += symmetry->character_table[i * symmetry->number_of_classes + j] * character[j];
           }
            num_irrep_in_basis[i] = (int) ((dot_product + 0.000001) / (double) symmetry->grp_dim);
            if (job->taskid > 0)
            fprintf(file.out,"BASIS IRREP %3d %3d     %10.4lf\n",i,num_irrep_in_basis[i], dot_product / (double) symmetry->grp_dim);
           }
            if (job->taskid > 0)
            fprintf(file.out,"\n");
         
}
*/

void count_atom_salc(int s, int ip, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4;
int num_coef, total_atom_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int count;
int offset, atm, atm_rot, index_i, shelposi, bfposi;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
DoubleMatrix *projection_operator1, *vector_new1;
double *vector_rot, *vector;
int dimp1, dimp2;

  total_atom_coef = 0;
  num_salc = 0;
  atm = atoms->uniq[ip];
  int ip1;
  ip1 = atom_p->posn[atm];
  ip = atm;
  salc->num_atom = atom_p->numb[atm];
  shelposi = atoms->shelposn_sh[ip1];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip1]; index_i++) {
    dimp1 = atom_p->numb[atm] * shells->shar[index_i];
    dimp2 = dimp1 * (symmetry->cls_num_k[s] + 1);
    AllocateDoubleArray(&vector,&dimp1,job);
    AllocateDoubleArray(&vector_rot,&dimp1,job);
    AllocateDoubleMatrix(&projection_operator1,&dimp1,&dimp1,job);
    AllocateDoubleMatrix(&vector_new1,&dimp2,&dimp1,job);
    ResetDoubleMatrix(vector_new1);
    num_irp = 0;
    for (i = 0; i < shells->shar[index_i]; i++) {
      ResetDoubleArray(vector,&dimp1);
      op = 0;
      for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
        for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          atm_rot = atom_p->K[ip1 * symmetry->number_of_operators + op];
          offset = (atm_rot - ip1) * shells->shar[index_i];
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i4 = 0; i4 < *p_num; i4++) {
            if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * \
            (double) symmetry->character_table[s * symmetry->number_of_classes + i2];
//if (i == *p_ind_j) fprintf(file.out,"ip1 %3d op %3d i2 %3d i3 %3d i4 %3d offset %3d p_rot %10.4f character %10.4f\n",ip1,op,i2,i3,i4,offset,*p_rot,symmetry->character_table[s * symmetry->number_of_classes + i2]);
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           }
          op++;
         } // close loop on i3
        } // close loop on i2
//for (i2 = 0; i2 < dimp1; i2++) { fprintf(file.out,"%3d %3d %6.3lf",s,index_i,vector[i2]); } fprintf(file.out,"\n");
          norm = k_zero;
          for (i1 = 0; i1 < dimp1; i1++) norm += vector[i1] * vector[i1];
            if (fabs(norm) > 0.00001) {
              for (i1 = 0; i1 < dimp1; i1++) vector[i1] /= sqrt(norm); 
             }
      for (op = 0; op < symmetry->number_of_operators; op++) {
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
        p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
        ResetDoubleArray(vector_rot,&dimp1);
        for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
          atm_rot = atom_p->K[(ip1 + i1) * symmetry->number_of_operators + op];
          offset = (atm_rot - ip1) * shells->shar[index_i];
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i2 = 0; i2 < *p_num; i2++) {
            vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           } // close loop on i2
          } // close loop on i1
          norm = k_zero;
          // Gram-Schmidt orthogonalise
          for (i2 = 0; i2 < dimp1; i2++) vector_new1->a[num_irp][i2] = k_zero;
          ResetDoubleMatrix(projection_operator1);
          for (i1 = 0; i1 < dimp1; i1++) 
          projection_operator1->a[i1][i1] = k_one;
          for (i1 = 0; i1 < num_irp; i1++) {
            for (i2 = 0; i2 < dimp1; i2++) {
               for (i3 = 0; i3 < dimp1; i3++) {
                 projection_operator1->a[i2][i3] -= vector_new1->a[i1][i2] * vector_new1->a[i1][i3];
                }
               }
              }
          for (i1 = 0; i1 < dimp1; i1++) {
            for (i2 = 0; i2 < dimp1; i2++) {
              vector_new1->a[num_irp][i1] += projection_operator1->a[i1][i2] * vector_rot[i2];
             }
            }
          norm = k_zero;
          for (i1 = 0; i1 < dimp1; i1++) 
          norm += vector_new1->a[num_irp][i1] * vector_new1->a[num_irp][i1];
          if (norm > 0.0000001) {
          for (i1 = 0; i1 < dimp1; i1++) vector_new1->a[num_irp][i1] /= sqrt(norm); 
        //fprintf(file.out,"%3d %3d %3d\n",i,s,index_i);
        //for (i2 = 0; i2 < dimp1; i2++) { fprintf(file.out,"%7.3lf",vector_new1->a[num_irp][i2]); } fprintf(file.out,"\n");
          total_atom_coef += shells->shar[index_i];
          //total_coef += shells->shar[index_i];
          //total_coef += dimp1;
          num_irp++;
          num_salc++;
         } // close if (norm
        } // close loop on op
          //fprintf(file.out,"\n");
         } // close loop on i
        DestroyDoubleArray(&vector,&dimp1,job);
        DestroyDoubleArray(&vector_rot,&dimp1,job);
        DestroyDoubleMatrix(&vector_new1,job);
        DestroyDoubleMatrix(&projection_operator1,job);
       } // close index_i
       salc->num_salc = num_salc;
       salc->total_coef = total_atom_coef;
       //fprintf(file.out,"counted %3d %3d %3d\n",total_coef,salc->num_atom,total_coef * salc->num_atom);

}

/*
void count_atom_salc_crystal(int s, IntMatrix *lattice_tran, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, REAL_LATTICE_TABLES *R_tables, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4;
int num_coef, total_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int count;
int offset, atm, atm_rot, index_i, shelposi, bfposi;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
DoubleMatrix *projection_operator1, *vector_new1;
double *vector_rot, *vector;
int dimp1, dimp2;

  //fprintf(file.out,"s %3d pair %3d posn %3d natoms %3d\n",s,pair,posn,count_atoms);
  //for (i2 = 0; i2 < count_atoms; i2++) 
  //fprintf(file.out,"%3d %3d %3d %3d\n",posn,count_atoms,array[i2],offset1[i2]); fflush(file.out);

  total_coef = 0;
  num_salc = 0;
int ip = 0; //FIX
  atm = atoms->uniq[ip];
  int ip1;
  ip1 = ip;
  shelposi = atoms->shelposn_sh[ip1];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip1]; index_i++) {
    dimp1 = salc->num_atom * shells->shar[index_i];
    dimp2 = dimp1 * (symmetry->cls_num_k[s] + 1);
    AllocateDoubleArray(&vector,&dimp1,job);
    AllocateDoubleArray(&vector_rot,&dimp1,job);
    AllocateDoubleMatrix(&projection_operator1,&dimp1,&dimp1,job);
    AllocateDoubleMatrix(&vector_new1,&dimp2,&dimp1,job);
    ResetDoubleMatrix(vector_new1);
    num_irp = 0;
    for (i = 0; i < shells->shar[index_i]; i++) {
      ResetDoubleArray(vector,&dimp1);
      op = 0;
      for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
        for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          atm_rot = atom_p->K[ip1 * symmetry->number_of_operators + op];
          offset = lattice_tran->a[0][op] * shells->shar[index_i];
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i4 = 0; i4 < *p_num; i4++) {
            if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * symmetry->character_table[s * symmetry->number_of_classes + i2];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           }
          op++;
         } // close loop on i3
        } // close loop on i2
          norm = k_zero;
          for (i1 = 0; i1 < dimp1; i1++) norm += vector[i1] * vector[i1];
            if (fabs(norm) > 0.00001) {
              for (i1 = 0; i1 < dimp1; i1++) vector[i1] /= sqrt(norm); 
             }
              //for (i2 = 0; i2 < dimp1; i2++) {
              //fprintf(file.out,"z%6.3lf",vector[i2]); } fprintf(file.out,"\n"); 
      for (op = 0; op < symmetry->number_of_operators; op++) {
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
        p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
        ResetDoubleArray(vector_rot,&dimp1);
        for (i1 = 0; i1 < salc->num_atom; i1++) {
          //atm_rot = atom_p->K[(ip1 + i1) * symmetry->number_of_operators + op];
          atm_rot = atom_p->K[ip1 * symmetry->number_of_operators + op];
          offset = lattice_tran->a[i1][op] * shells->shar[index_i];
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i2 = 0; i2 < *p_num; i2++) {
            vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           } // close loop on i2
          } // close loop on i1
           //for (i2 = 0; i2 < dimp1; i2++) {
           //fprintf(file.out,"Z%6.3lf",vector_rot[i2]); } fprintf(file.out,"\n"); 
          norm = k_zero;
          // Gram-Schmidt orthogonalise
          for (i2 = 0; i2 < dimp1; i2++) vector_new1->a[num_irp][i2] = k_zero;
          ResetDoubleMatrix(projection_operator1);
          for (i1 = 0; i1 < dimp1; i1++) 
          projection_operator1->a[i1][i1] = k_one;
          for (i1 = 0; i1 < num_irp; i1++) {
            for (i2 = 0; i2 < dimp1; i2++) {
               for (i3 = 0; i3 < dimp1; i3++) {
                 projection_operator1->a[i2][i3] -= vector_new1->a[i1][i2] * vector_new1->a[i1][i3];
                }
               }
              }
          for (i1 = 0; i1 < dimp1; i1++) {
            for (i2 = 0; i2 < dimp1; i2++) {
              vector_new1->a[num_irp][i1] += projection_operator1->a[i1][i2] * vector_rot[i2];
             }
            }
          norm = k_zero;
          for (i1 = 0; i1 < dimp1; i1++) 
          norm += vector_new1->a[num_irp][i1] * vector_new1->a[num_irp][i1];
          if (norm > 0.0000001) {
          for (i1 = 0; i1 < dimp1; i1++) vector_new1->a[num_irp][i1] /= sqrt(norm); 
          //fprintf(file.out,"%3d %3d %3d\n",i,s,index_i);
          for (i2 = 0; i2 < dimp1; i2++) {
          fprintf(file.out,"%7.3lf",vector_new1->a[num_irp][i2]); } fprintf(file.out,"\n"); 
          total_coef += dimp1;
          num_irp++;
          num_salc++;
         } // close if (norm
        } // close loop on op
          //fprintf(file.out,"\n");
         } // close loop on i
        DestroyDoubleArray(&vector,&dimp1,job);
        DestroyDoubleArray(&vector_rot,&dimp1,job);
        DestroyDoubleMatrix(&vector_new1,job);
        DestroyDoubleMatrix(&projection_operator1,job);
       } // close index_i
       salc->num_salc = num_salc;
       salc->total_coef = total_coef;
       fprintf(file.out,"s %3d num_salc %3d total_coef %3d\n",s,num_salc,total_coef);

}
*/

void generate_atom_salc(int s, int ip, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4, ip1;
int num_coef, total_atom_coef;
//int num_coef, total_coef, total_atom_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int offset, atm, atm_rot, index_i, shelposi, bfposi;
int dimp1, dimp2;
int count;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
double *vector_rot, *vector;
DoubleMatrix *projection_operator, *vector_new;

  for (i = 0; i < salc->num_atom; i++) { for (j = 0; j < salc->total_coef; j++) { salc->coeff->a[i][j] = k_zero; }}
  //total_coef = 0;
  total_atom_coef = 0;
  num_salc = 0;
  atm = atoms->uniq[ip];
  ip = atm;
  ip1 = atom_p->posn[atm];
  salc->num_atom = atom_p->numb[atm];
  shelposi = atoms->shelposn_sh[ip1];
  bfposi = 0;
  salc->num_irp[s] = 0;
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip1]; index_i++) {
    dimp1 = atom_p->numb[atm] * shells->shar[index_i];
    dimp2 = dimp1 * (symmetry->cls_num_k[s] + 1);
    AllocateDoubleArray(&vector,&dimp1,job);
    AllocateDoubleArray(&vector_rot,&dimp1,job);
    AllocateDoubleMatrix(&projection_operator,&dimp1,&dimp1,job);
    AllocateDoubleMatrix(&vector_new,&dimp2,&dimp1,job);
    num_irp = 0;
    for (i = 0; i < shells->shar[index_i]; i++) {
      ResetDoubleArray(vector,&dimp1);
      op = 0;
      for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
        for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          atm_rot = atom_p->K[ip1 * symmetry->number_of_operators + op];
          offset = (atm_rot - ip1) * shells->shar[index_i]; 
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i4 = 0; i4 < *p_num; i4++) {
            if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * (double) symmetry->character_table[s * symmetry->number_of_classes + i2];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           }
          op++;
         } // close loop on i3
        } // close loop on i2
      norm = k_zero;
      for (i1 = 0; i1 < dimp1; i1++) norm += vector[i1] * vector[i1];
        if (fabs(norm) > 0.00001) {
          for (i1 = 0; i1 < dimp1; i1++) vector[i1] /= sqrt(norm); 
         }
      for (op = 0; op < symmetry->number_of_operators; op++) {
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
        p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
        ResetDoubleArray(vector_rot,&dimp1);
        for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
          atm_rot = atom_p->K[(ip1 + i1) * symmetry->number_of_operators + op];
          offset = (atm_rot - ip1) * shells->shar[index_i];
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i2 = 0; i2 < *p_num; i2++) {
            vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           } // close loop on i2
          } // close loop on i1
          norm = k_zero;
          // Gram-Schmidt orthogonalise
          for (i2 = 0; i2 < dimp1; i2++) vector_new->a[num_irp][i2] = k_zero;
          ResetDoubleMatrix(projection_operator);
          for (i1 = 0; i1 < dimp1; i1++) {
              projection_operator->a[i1][i1] = k_one;
             }
          for (i1 = 0; i1 < num_irp; i1++) {
            for (i2 = 0; i2 < dimp1; i2++) {
              for (i3 = 0; i3 < dimp1; i3++) {
                projection_operator->a[i2][i3] -= vector_new->a[i1][i2] * vector_new->a[i1][i3];
               }
              }
             }
          for (i1 = 0; i1 < dimp1; i1++) {
            for (i2 = 0; i2 < dimp1; i2++) {
              vector_new->a[num_irp][i1] += projection_operator->a[i1][i2] * vector_rot[i2];
             }
            }
          norm = k_zero;
          num_coef = 0;
          for (i1 = 0; i1 < dimp1; i1++) norm += vector_new->a[num_irp][i1] * vector_new->a[num_irp][i1];
	  //fprintf(file.out,"s %3d index_i %3d xyz %3d op %3d norm %18.12f\n",s,index_i,i,op,norm);
          if (norm > 0.0000001) {
          for (i1 = 0; i1 < dimp1; i1++) vector_new->a[num_irp][i1] /= sqrt(norm); 
          for (i2 = 0; i2 < shells->shar[index_i]; i2++) {
            for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
              salc->coeff->a[i1][total_atom_coef] = vector_new->a[num_irp][i1 * shells->shar[index_i] + i2];
              salc->bfn_posn->a[i1][total_atom_coef] = bfposi + i2;
              //total_coef++;
             }
              total_atom_coef++;
              num_coef++;
             }
              salc->num_coef[num_salc] = num_coef;
             (salc->num_irp[s])++;
              salc->irp[num_salc] = s;
	      //fprintf(file.out,"s %3d op %3d num_irp %3d norm %14.8f\n",s,op,salc->num_irp[s],norm);
              num_irp++;
              num_salc++;
             } // close if (norm
            } // close loop on op
           } // close loop on i
        bfposi += shells->shar[index_i];
        DestroyDoubleArray(&vector,&dimp1,job);
        DestroyDoubleArray(&vector_rot,&dimp1,job);
        DestroyDoubleMatrix(&vector_new,job);
        DestroyDoubleMatrix(&projection_operator,job);
        } // close index_i
        salc->num_salc = num_salc;
        //salc->total_coef = total_coef;
        salc->total_coef = total_atom_coef;
        //fprintf(file.out,"needed %3d\n",total_atom_coef);

    //   bfposi = 0;
    //   total_atom_coef = 0;
    //   shelposi = atoms->shelposn_sh[ip1];
    //   if (job->verbosity > 1)
    //   for (j = 0; j < salc->num_irp[s]; j++) {
    //     for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip1]; index_i++) {
    //       for (i2 = 0; i2 < shells->shar[index_i]; i2++) {
    //         for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
    //           fprintf(file.out,"class %3d atm %3d shell %3d xyz %3d atmp %3d salc %3d coeff %14.8f %3d\n",\
    //  	       s,atm,index_i,i2,i1,j,salc->coeff->a[i1][total_atom_coef],total_atom_coef);
    //          }
    //           total_atom_coef++;
    //          }
    //         fprintf(file.out,"\n");
    //         bfposi += shells->shar[index_i];
    //        }
    //         fprintf(file.out,"\n");
    //        }

       int count_j = 0;
       if (job->verbosity > 1)
       for (int r = 0; r < salc->num_irp[s]; r++){
         for (int t = 0; t < salc->num_coef[r]; t++){
           for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
             fprintf(file.out,"irp %3d i1 %3d t %3d count %3d %3d %14.8f\n",\
	     r,i1,t,count_j+t,salc->bfn_posn->a[i1][count_j+t],salc->coeff->a[i1][count_j+t]);
            }
           }
	  fprintf(file.out,"\n");
          count_j += salc->num_coef[r];
         }
}

/*
void generate_atom_salc_crystal(int s, IntMatrix *lattice_tran, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, REAL_LATTICE_TABLES *R_tables, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4;
int num_coef, total_coef, total_atom_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int count;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
int offset, atm, atm_rot, index_i, shelposi, bfposi;
DoubleMatrix *projection_operator, *vector_new;
double *vector_rot, *vector;
int dimp1, dimp2;

int ip, latt, latt_rot, O_temp, latt_temp;
//
if (pair == 0) {
ip = pair_p->cell1[posn];
latt = pair_p->latt1[posn];
}
else if (pair == 1) {
ip = pair_p->cell2[posn];
latt = pair_p->latt2[posn];
}

int array[symmetry->number_of_operators];
int oper[symmetry->number_of_operators];
int count_atoms = 0;
int lat_rot[symmetry->number_of_operators];
for (i1 = 0; i1 < symmetry->number_of_operators; i1++) array[i1] = -1;
for (op = 0; op < symmetry->number_of_operators; op++) {
  atm_rot  = atom_p->K[ip * symmetry->number_of_operators + op];
  O_temp   = atom_p->O[ip * symmetry->number_of_operators + op];
  latt_rot = R_tables->diffvec[R_tables->lattvec[latt * symmetry->number_of_operators + op] * R_tables->margin_vector + O_temp];
  lat_rot[op] = R_tables->diffvec[R_tables->lattvec[latt * symmetry->number_of_operators + op] * R_tables->margin_vector + O_temp];
  for (i2 = 0; i2 <= count_atoms; i2++) {
    //fprintf(file.out,"%3d %3d %3d %3d\n",latt,latt_rot,count_atoms,array[i2]); fflush(file.out);
    if (latt_rot == array[i2]) break;
    if (latt_rot != array[i2] && i2 == count_atoms) { array[i2] = latt_rot; oper[i2] = op; count_atoms++; break; }
   }
  }
int latt_tran[count_atoms][symmetry->number_of_operators];
int latt_rotd[count_atoms][symmetry->number_of_operators];
for (i2 = 0; i2 < count_atoms; i2++) {
for (op = 0; op < symmetry->number_of_operators; op++) {
latt_tran[i2][op] = lat_rot[symmetry->grp_k[oper[i2] * symmetry->number_of_operators + op]];
  //fprintf(file.out,"s %3d i2 %3d op %3d  %3d \n",s,i2,op,latt_tran[i2][op]);
}
}
for (i1 = 0; i1 < count_atoms; i1++) {
for (i2 = 0; i2 < count_atoms; i2++) {
for (op = 0; op < symmetry->number_of_operators; op++) {
if (latt_tran[i2][op] == array[i1]) { latt_rotd[i2][op] = i1; } //fprintf(file.out,"%3d %3d %3d\n",i2,op,latt_rotd[i2][op]); }
}
}
}
//

  //fprintf(file.out,"s %3d pair %3d posn %3d natoms %3d\n",s,pair,posn,count_atoms);
  //for (i2 = 0; i2 < count_atoms; i2++) 
  //fprintf(file.out,"%3d %3d %3d %3d\n",count_atoms,array[i2],posn,pair_p->numb[posn]); fflush(file.out);
  total_coef = 0;
  total_atom_coef = 0;
  num_salc = 0;
  atm = atoms->uniq[ip];
  int ip1;
  ip = ip1;
  shelposi = atoms->shelposn_sh[ip1];
  bfposi = 0;
  salc->num_irp[s] = 0;
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip1]; index_i++) {
    dimp1 = salc->num_atom * shells->shar[index_i];
    dimp2 = dimp1 * (symmetry->cls_num_k[s] + 1);
    AllocateDoubleArray(&vector,&dimp1,job);
    AllocateDoubleArray(&vector_rot,&dimp1,job);
    AllocateDoubleMatrix(&projection_operator,&dimp1,&dimp1,job);
    AllocateDoubleMatrix(&vector_new,&dimp2,&dimp1,job);
    num_irp = 0;
    for (i = 0; i < shells->shar[index_i]; i++) {
      ResetDoubleArray(vector,&dimp1);
        op = 0;
        for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
          for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          atm_rot = atom_p->K[ip1 * symmetry->number_of_operators + op];
          offset = lattice_tran->a[0][op] * shells->shar[index_i];
            p_ind_i = symmetry->ind_i  + oppshift1;
            p_ind_j = symmetry->ind_j  + oppshift1;
            p_rot   = symmetry->rot    + oppshift1;
            for (i4 = 0; i4 < *p_num; i4++) {
              if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * \
              (double) symmetry->character_table[s * symmetry->number_of_classes + i2];
              p_ind_i++;
              p_ind_j++;
              p_rot++;
             }
            op++;
           } // close loop on i3
          } // close loop on i2
            norm = k_zero;
            for (i1 = 0; i1 < dimp1; i1++) norm += vector[i1] * vector[i1];
              if (fabs(norm) > 0.00001) {
                for (i1 = 0; i1 < dimp1; i1++) vector[i1] /= sqrt(norm); 
               }
              //for (i2 = 0; i2 < dimp1; i2++) {
              //fprintf(file.out,"z%6.3lf",vector[i2]); } fprintf(file.out,"\n"); 
        for (op = 0; op < symmetry->number_of_operators; op++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          ResetDoubleArray(vector_rot,&dimp1);
          for (i1 = 0; i1 < salc->num_atom; i1++) {
          atm_rot = atom_p->K[ip1 * symmetry->number_of_operators + op];
          offset = lattice_tran->a[i1][op] * shells->shar[index_i];
            p_ind_i = symmetry->ind_i  + oppshift1;
            p_ind_j = symmetry->ind_j  + oppshift1;
            p_rot   = symmetry->rot    + oppshift1;
            for (i2 = 0; i2 < *p_num; i2++) {
              vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
              p_ind_i++;
              p_ind_j++;
              p_rot++;
             } // close loop on i2
            } // close loop on i1
            norm = k_zero;
            // Gram-Schmidt orthogonalise
            for (i2 = 0; i2 < dimp1; i2++) vector_new->a[num_irp][i2] = k_zero;
            ResetDoubleMatrix(projection_operator);
            for (i1 = 0; i1 < dimp1; i1++) 
            projection_operator->a[i1][i1] = k_one;
            for (i1 = 0; i1 < num_irp; i1++) {
              for (i2 = 0; i2 < dimp1; i2++) {
                 for (i3 = 0; i3 < dimp1; i3++) {
                   projection_operator->a[i2][i3] -= vector_new->a[i1][i2] * vector_new->a[i1][i3];
                  }
                 }
                }
            for (i1 = 0; i1 < dimp1; i1++) {
              for (i2 = 0; i2 < dimp1; i2++) {
                vector_new->a[num_irp][i1] += projection_operator->a[i1][i2] * vector_rot[i2];
               }
              }
            norm = k_zero;
            num_coef = 0;
            for (i1 = 0; i1 < dimp1; i1++) 
            norm += vector_new->a[num_irp][i1] * vector_new->a[num_irp][i1];
            if (norm > 0.0000001) {
            for (i1 = 0; i1 < dimp1; i1++) vector_new->a[num_irp][i1] /= sqrt(norm); 
            for (i2 = 0; i2 < dimp1; i2++) {
            fprintf(file.out,"%7.3lf",vector_new->a[num_irp][i2]); } fprintf(file.out,"\n"); 

            for (i2 = 0; i2 < shells->shar[index_i]; i2++) {
              for (i1 = 0; i1 < salc->num_atom; i1++) {
                salc->coeff->a[i1][total_atom_coef]  = vector_new->a[num_irp][i1 * shells->shar[index_i] + i2];  // FIX num_salc
                salc->bfn_posn->a[i1][total_atom_coef]  = bfposi + i2;
                total_coef++;
               }
                total_atom_coef++;
                num_coef++;
               }

                salc->num_coef[num_salc] = num_coef;
                //fprintf(file.out,"s %3d num_salc %3d num_coef %3d\n",s,num_salc,num_coef);
               (salc->num_irp[s])++;
                salc->irp[num_salc] = s;
                num_irp++;
                num_salc++;
               } // close if (norm
              } // close loop on op
             } // close loop on i
          //fprintf(file.out,"\n");
          bfposi += shells->shar[index_i];
          DestroyDoubleArray(&vector,&dimp1,job);
          DestroyDoubleArray(&vector_rot,&dimp1,job);
          DestroyDoubleMatrix(&vector_new,job);
          DestroyDoubleMatrix(&projection_operator,job);
          } // close index_i
          salc->num_salc = num_salc;
          salc->total_coef = total_coef;

}
*/

void count_shell_salc(int index_i, int s, int atm, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4;
int num_coef, total_coef, total_atom_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int count;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
int offset, atm_rot, shelposi, bfposi;
int atm1;
DoubleMatrix *projection_operator, *vector_new;
double *vector_rot, *vector;
int dimp1, dimp2;
atm1 = atom_p->posn[atm];

  total_coef = 0;
  total_atom_coef = 0;
  num_salc = 0;
  salc->num_atom = atom_p->numb[atm];
  shelposi = atoms->shelposn[atm];
  bfposi = 0;
  //double vector[atom_p->numb[atm] * shells->shar[index_i]];
  //double vector_rot[atom_p->numb[atm] * shells->shar[index_i]];
  //double vector_new[atom_p->numb[atm] * shells->shar[index_i] * symmetry->cls_num_k[s]][atom_p->numb[atm]*shells->shar[index_i]];
  //double projection_operator[atom_p->numb[atm] * shells->shar[index_i]][atom_p->numb[atm] * shells->shar[index_i]];
  dimp1 = atom_p->numb[atm] * shells->shar[index_i];
  dimp2 = dimp1 * (symmetry->cls_num_k[s] + 1);
  AllocateDoubleArray(&vector,&dimp1,job);
  AllocateDoubleArray(&vector_rot,&dimp1,job);
  AllocateDoubleMatrix(&projection_operator,&dimp1,&dimp1,job);
  AllocateDoubleMatrix(&vector_new,&dimp2,&dimp1,job);
  num_irp = 0;
  for (i = 0; i < shells->shar[index_i]; i++) {
    ResetDoubleArray(vector,&dimp1);
    //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector[i1] = 0;
      op = 0;
      for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
        for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          atm_rot = atom_p->K[atm1 * symmetry->number_of_operators + op];
          offset = (atm_rot - atm1) * shells->shar[index_i]; // check atm_rot or ip??
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i4 = 0; i4 < *p_num; i4++) {
            if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * \
            (double) symmetry->character_table[s * symmetry->number_of_classes + i2];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           }
          op++;
         } // close loop on i3
        } // close loop on i2
  norm = k_zero;
  //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) norm += vector[i1] * vector[i1];
  for (i1 = 0; i1 < dimp1; i1++) norm += vector[i1] * vector[i1];
    if (fabs(norm) > 0.00001) {
      //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector[i1] /= sqrt(norm); 
      for (i1 = 0; i1 < dimp1; i1++) vector[i1] /= sqrt(norm); 
     }
      for (op = 0; op < symmetry->number_of_operators; op++) {
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
        p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
        ResetDoubleArray(vector_rot,&dimp1);
        //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector_rot[i1] = k_zero;
        for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
          atm_rot = atom_p->K[(atm1 + i1) * symmetry->number_of_operators + op];
          offset = (atm_rot - atm1) * shells->shar[index_i];
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i2 = 0; i2 < *p_num; i2++) {
            vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           } // close loop on i2
          } // close loop on i1
          norm = k_zero;
          // Gram-Schmidt orthogonalise
          //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) vector_new[num_irp][i2] = k_zero;
          for (i2 = 0; i2 < dimp1; i2++) vector_new->a[num_irp][i2] = k_zero;
          ResetDoubleMatrix(projection_operator);
          for (i1 = 0; i1 < dimp1; i1++)
          projection_operator->a[i1][i1] = k_one;
          //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
            //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
              //projection_operator[i1][i2] = k_zero;
             //}
              //projection_operator[i1][i1] = k_one;
             //}
           for (i1 = 0; i1 < num_irp; i1++) {
             for (i2 = 0; i2 < dimp1; i2++) {
               for (i3 = 0; i3 < dimp1; i3++) {
             //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
               //for (i3 = 0; i3 < atom_p->numb[atm] * shells->shar[index_i]; i3++) {
                 //projection_operator[i2][i3] -= vector_new[i1][i2] * vector_new[i1][i3];
                 projection_operator->a[i2][i3] -= vector_new->a[i1][i2] * vector_new->a[i1][i3];
                }
               }
              }
           for (i1 = 0; i1 < dimp1; i1++) {
             for (i2 = 0; i2 < dimp1; i2++) {
           //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
             //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
               //vector_new[num_irp][i1] += projection_operator[i1][i2] * vector_rot[i2];
               vector_new->a[num_irp][i1] += projection_operator->a[i1][i2] * vector_rot[i2];
             }
            }
          norm = k_zero;
          num_coef = 0;
          //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) 
          //norm += vector_new[num_irp][i1] * vector_new[num_irp][i1];
          for (i1 = 0; i1 < dimp1; i1++) 
          norm += vector_new->a[num_irp][i1] * vector_new->a[num_irp][i1];
          if (norm > 0.0000001) {
          //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector_new[num_irp][i1] /= sqrt(norm); 
          for (i1 = 0; i1 < dimp1; i1++) vector_new->a[num_irp][i1] /= sqrt(norm); 
          //fprintf(file.out,"%3d %3d %3d\n",i,s,index_i);
          //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
          //fprintf(file.out,"%7.3lf",vector_new[num_irp][i2]); } fprintf(file.out,"\n"); 
          for (i2 = 0; i2 < shells->shar[index_i]; i2++) {
          for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
            total_coef++;
           }
            total_atom_coef++;
            num_coef++;
           }
            num_irp++;
            num_salc++;
          } // close if (norm
         } // close loop on op
        } // close loop on i
       DestroyDoubleArray(&vector,&dimp1,job);
       DestroyDoubleArray(&vector_rot,&dimp1,job);
       DestroyDoubleMatrix(&vector_new,job);
       DestroyDoubleMatrix(&projection_operator,job);
       salc->num_salc = num_salc;
       salc->total_coef = total_coef;

}

void generate_shell_salc(int index_i, int s, int atm, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4;
int num_coef, total_coef, total_atom_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int count;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
int offset, atm_rot, shelposi, bfposi;
int atm1;
DoubleMatrix *projection_operator, *vector_new;
double *vector_rot, *vector;
int dimp1, dimp2;
atm1 = atom_p->posn[atm];

  total_coef = 0;
  total_atom_coef = 0;
  num_salc = 0;
  salc->num_atom = atom_p->numb[atm];
  shelposi = atoms->shelposn_sh[atm];
  bfposi = 0;
  salc->num_irp[s] = 0;
  //double vector[atom_p->numb[atm] * shells->shar[index_i]];
  //double vector_rot[atom_p->numb[atm] * shells->shar[index_i]];
  //double vector_new[atom_p->numb[atm] * shells->shar[index_i] * symmetry->cls_num_k[s]][atom_p->numb[atm]*shells->shar[index_i]];
  //double projection_operator[atom_p->numb[atm] * shells->shar[index_i]][atom_p->numb[atm] * shells->shar[index_i]];
  dimp1 = atom_p->numb[atm] * shells->shar[index_i];
  dimp2 = dimp1 * (symmetry->cls_num_k[s] + 1);
  AllocateDoubleArray(&vector,&dimp1,job);
  AllocateDoubleArray(&vector_rot,&dimp1,job);
  AllocateDoubleMatrix(&projection_operator,&dimp1,&dimp1,job);
  AllocateDoubleMatrix(&vector_new,&dimp2,&dimp1,job);
  num_irp = 0;
  for (i = 0; i < shells->shar[index_i]; i++) {
    ResetDoubleArray(vector,&dimp1);
    //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector[i1] = 0;
      op = 0;
      for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
        for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          atm_rot = atom_p->K[atm1 * symmetry->number_of_operators + op];
          offset = (atm_rot - atm1) * shells->shar[index_i]; 
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i4 = 0; i4 < *p_num; i4++) {
            if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * (double) symmetry->character_table[s * \
            symmetry->number_of_classes + i2];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           }
          op++;
         } // close loop on i3
        } // close loop on i2
  norm = k_zero;
  //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) norm += vector[i1] * vector[i1];
  for (i1 = 0; i1 < dimp1; i1++) norm += vector[i1] * vector[i1];
    if (fabs(norm) > 0.00001) {
      //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector[i1] /= sqrt(norm); 
      for (i1 = 0; i1 < dimp1; i1++) vector[i1] /= sqrt(norm); 
     }
      for (op = 0; op < symmetry->number_of_operators; op++) {
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
        p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
        ResetDoubleArray(vector_rot,&dimp1);
        //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector_rot[i1] = k_zero;
        for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
          atm_rot = atom_p->K[(atm1 + i1) * symmetry->number_of_operators + op];
          offset = (atm_rot - atm1) * shells->shar[index_i];
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i2 = 0; i2 < *p_num; i2++) {
            vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           } // close loop on i2
          } // close loop on i1
          norm = k_zero;
          // Gram-Schmidt orthogonalise
          //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) vector_new[num_irp][i2] = k_zero;
          for (i2 = 0; i2 < dimp1; i2++) vector_new->a[num_irp][i2] = k_zero;
          ResetDoubleMatrix(projection_operator);
          for (i1 = 0; i1 < dimp1; i1++)
          projection_operator->a[i1][i1] = k_one;
          //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
            //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
              //projection_operator->a[i1][i2] = k_zero;
             //}
              //projection_operator->a[i1][i1] = k_one;
             //}
           for (i1 = 0; i1 < num_irp; i1++) {
             for (i2 = 0; i2 < dimp1; i2++) {
               for (i3 = 0; i3 < dimp1; i3++) {
             //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
               //for (i3 = 0; i3 < atom_p->numb[atm] * shells->shar[index_i]; i3++) {
                 //projection_operator[i2][i3] -= vector_new[i1][i2] * vector_new[i1][i3];
                 projection_operator->a[i2][i3] -= vector_new->a[i1][i2] * vector_new->a[i1][i3];
                }
               }
              }
           for (i1 = 0; i1 < dimp1; i1++) {
             for (i2 = 0; i2 < dimp1; i2++) {
           //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
             //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
               //vector_new[num_irp][i1] += projection_operator[i1][i2] * vector_rot[i2];
               vector_new->a[num_irp][i1] += projection_operator->a[i1][i2] * vector_rot[i2];
             }
            }
          norm = k_zero;
          num_coef = 0;
          //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) 
          //norm += vector_new[num_irp][i1] * vector_new[num_irp][i1];
          for (i1 = 0; i1 < dimp1; i1++) 
          norm += vector_new->a[num_irp][i1] * vector_new->a[num_irp][i1];
          if (norm > 0.0000001) {
          //for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector_new[num_irp][i1] /= sqrt(norm); 
          for (i1 = 0; i1 < dimp1; i1++) vector_new->a[num_irp][i1] /= sqrt(norm); 
          //fprintf(file.out,"class %3d sph. harm. basis %3d shell%3d\n",s,i,index_i);
          //for (i2 = 0; i2 < dimp1; i2++) {
          //fprintf(file.out,"%7.3lf",vector_new->a[num_irp][i2]); } fprintf(file.out,"\n"); 
          for (i2 = 0; i2 < shells->shar[index_i]; i2++) {
          for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
            salc->atm[total_coef] = i1; 
            //salc->coeff->a[i1][total_atom_coef]  = vector_new[num_irp][i1 * shells->shar[index_i] + i2];  // FIX num_salc
            salc->coeff->a[i1][total_atom_coef]  = vector_new->a[num_irp][i1 * shells->shar[index_i] + i2];  // FIX num_salc
            salc->bfn_posn->a[i1][total_atom_coef]  = bfposi + i2;
            total_coef++;
           }
            total_atom_coef++;
            num_coef++;
           }
            salc->num_coef[num_salc] = num_coef;
           (salc->num_irp[s])++;
            salc->irp[num_salc] = s;
            num_irp++;
            num_salc++;
          } // close if (norm
         } // close loop on op
        } // close loop on i
       //fprintf(file.out,"\n");
       DestroyDoubleArray(&vector,&dimp1,job);
       DestroyDoubleArray(&vector_rot,&dimp1,job);
       DestroyDoubleMatrix(&vector_new,job);
       DestroyDoubleMatrix(&projection_operator,job);
       salc->num_salc = num_salc;
       salc->total_coef = total_coef;

}

void transform_salc_to_atomic_basis(ComplexMatrix *eigenvectors2, double *eigval, int *irrep, ComplexMatrix **eigenvectors_salc, KPOINT_TRAN *knet, int nk[2], ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)
{

int r, s, t, u, v;
int atm, atm1, atm_temp, index_i, shelposi, bfposi, bfn_temp, count_salc, count_irrep, count_irrep1;
int temp, eigval_ord[atoms->number_of_sh_bfns_in_unit_cell], eigval_ord_inverse[atoms->number_of_sh_bfns_in_unit_cell];
double temp_val;
SALC salc1;

  int num_irrep_in_basis[symmetry->number_of_classes];
  count_basis_irrep(num_irrep_in_basis,atoms,atom_p,shells,symmetry,job,file);

  for (r = 0; r < atoms->number_of_sh_bfns_in_unit_cell; r++)
  eigval_ord[r] = r;
  for (r = 1; r < atoms->number_of_sh_bfns_in_unit_cell; r++) {
    for (s = atoms->number_of_sh_bfns_in_unit_cell - 1; s >= r; --s) {
      if ((eigval[s - 1]) > (eigval[s])) {

        temp_val = eigval[s - 1];
        eigval[s - 1] = eigval[s];
        eigval[s] = temp_val;

        temp = eigval_ord[s - 1];
        eigval_ord[s - 1] = eigval_ord[s];
        eigval_ord[s] = temp;

        temp = irrep[s - 1];
        irrep[s - 1] = irrep[s];
        irrep[s] = temp;

      }
    }
  }

  for (r = 0; r < atoms->number_of_sh_bfns_in_unit_cell; r++)
  eigval_ord_inverse[eigval_ord[r]] = r;

  //for (int i = 0; i < atoms->number_of_sh_bfns_in_unit_cell; i++) 
  //fprintf(file.out,"salc eigenvalues in order %3d %3d %10.4f\n",i,irrep[i],eigval[i]);

  count_irrep1 = 0;
  for (v = 0; v < symmetry->number_of_classes; v++) {
    bfposi = 0;
    atm_temp = 0;
    count_irrep = 0;
    for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
      atm1 = atom_p->posn[atm];
      shelposi = atoms->shelposn_sh[atm1];
      bfn_temp = 0;
      for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[atm1]; index_i++) {
        count_shell_salc(index_i,v,atm,&salc1,atoms,atom_p,pair_p,shells,symmetry,job,file);
        allocate_SALC(&salc1,symmetry,job,file);
        generate_shell_salc(index_i,v,atm,&salc1,atoms,atom_p,pair_p,shells,symmetry,job,file);
        count_salc = 0;
        for (r = 0; r < salc1.num_irp[v]; r++) {
          for (s = 0; s < shells->shar[index_i]; s++) {
            for (t = 0; t < atom_p->numb[atm]; t++) {
              bfposi = atoms->bfnposn_sh[atm_temp + t] + bfn_temp + s;
              for (u = 0; u < num_irrep_in_basis[v] * symmetry->irp_dim_k[v]; u++) {
                eigenvectors2->a[eigval_ord_inverse[count_irrep1 + u]][bfposi] += \
                eigenvectors_salc[v]->a[u][count_irrep] * salc1.coeff->a[t][count_salc];
                //eigenvectors2->a[count_irrep1 + u][bfposi] += \
                eigenvectors_salc[v]->a[u][count_irrep] * salc1.coeff->a[t][count_salc];
               } // close loop on u
              } // close loop on t
             count_salc++;
            } // close loop on s
           count_irrep++;
          } // close loop on r
         free_SALC(&salc1,job);
         bfn_temp += shells->shar[index_i];
        } // close loop on index_i
        atm_temp += atom_p->numb[atm];
       } // close loop on atm
      count_irrep1 += num_irrep_in_basis[v] * symmetry->irp_dim_k[v];
     } // close loop on v

}

void transform_degenerate_salc(ComplexMatrix **eigvec_salc, double **eigenvalues_salc, int *iii, int *num_irrep_in_basis, KPOINT_TRAN *knet, int nk[2], SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  double x1, x2, norm_fac, value, largest;
  int ii, jj, kk, ll, mm, n_largest, dim_salc;
  ComplexMatrix *eigvec_tmp;
  DoubleMatrix *rot_mat;
  dim_salc = symmetry->irp_dim_k[*iii] * num_irrep_in_basis[*iii];
  AllocateDoubleMatrix(&rot_mat,&symmetry->irp_dim_k[*iii],&symmetry->irp_dim_k[*iii],job);
  AllocateComplexMatrix(&eigvec_tmp,&dim_salc,&dim_salc,job);
  ResetComplexMatrix(eigvec_tmp);
  //fprintf(file.out,"salc eigenvectors before transformation to (x,y)\n");
  //print_complex_eigenvector_matrix2(eigvec_salc[*iii],eigenvalues_salc[*iii],6,12,1.0,file);

  for (ii = 0; ii < num_irrep_in_basis[*iii]; ii++) { 
    n_largest = 0;
    largest = k_zero;
    for (jj = 0; jj < dim_salc; jj+=2) { 
      value = fabs(eigvec_salc[*iii]->a[2 * ii][jj].real());
      n_largest = (value > largest) ? jj : n_largest;
      largest =   (value > largest) ? value : largest;
     }
      x1 = (eigvec_salc[*iii]->a[2 * ii][n_largest]).real();
      x2 = (eigvec_salc[*iii]->a[2 * ii][n_largest + 1]).real();
      //fprintf(file.out,"largest %3d %3d %10.4lf %10.4lf ",ii,n_largest,x1,x2);
      norm_fac = -sqrt(x1 * x1 + x2 * x2);
      if (fabs(norm_fac) > k_zero) {
      x1 /= norm_fac;
      x2 /= norm_fac;
     }
      //fprintf(file.out,"%10.4lf %10.4lf\n",x1,x2);
      rot_mat->a[0][0] =  x1;
      rot_mat->a[0][1] =  x2;
      rot_mat->a[1][0] = -x2;
      rot_mat->a[1][1] =  x1;
      for (jj = 0; jj < 2; jj++) {
        for (kk = 0; kk < num_irrep_in_basis[*iii]; kk++) { 
          for (ll = 0; ll < 2; ll++) {
            for (mm = 0; mm < 2; mm++) {
              eigvec_tmp->a[2 * ii + jj][2 * kk + ll] += rot_mat->a[ll][mm] * eigvec_salc[*iii]->a[2 * ii + jj][2 * kk + mm];
             }
            }
           }
          }
         } // close loop on ii
      for (jj = 0; jj < dim_salc; jj++) {
        for (kk = 0; kk < dim_salc; kk++) { 
          eigvec_salc[*iii]->a[jj][kk] = eigvec_tmp->a[jj][kk];
       }
      }
        //fprintf(file.out,"salc eigenvectors after transformation to (x,y)\n");
        //print_complex_eigenvector_matrix2(eigvec_salc[*iii],eigenvalues_salc[*iii],6,12,1.0,file);
        DestroyDoubleMatrix(&rot_mat,job);
        DestroyComplexMatrix(&eigvec_tmp,job);

}

/*
void count_atom_salc(int s, int ip, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4;
int num_coef, total_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int count;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
int offset, atm, atm_rot, index_i, shelposi, bfposi;

  total_coef = 0;
  num_salc = 0;
  atm = atoms->uniq[ip];
int ip1;
ip1 = atom_p->posn[atm];
ip = atm;
  salc->num_atom = atom_p->numb[atm];
  shelposi = atoms->shelposn_sh[ip1];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip1]; index_i++) {
    double vector[atom_p->numb[atm] * shells->shar[index_i]];
    double vector_rot[atom_p->numb[atm] * shells->shar[index_i]];
    double vector_new[atom_p->numb[atm] * shells->shar[index_i] * symmetry->cls_num_k[s]][atom_p->numb[atm]*shells->shar[index_i]];
    double projection_operator[atom_p->numb[atm] * shells->shar[index_i]][atom_p->numb[atm] * shells->shar[index_i]];
    num_irp = 0;
    for (i = 0; i < shells->shar[index_i]; i++) {
      for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector[i1] = k_zero;
        op = 0;
        for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
          for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
            oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
            p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
            atm_rot = atom_p->K[ip1 * symmetry->number_of_operators + op];
            offset = (atm_rot - ip1) * shells->shar[index_i];
            p_ind_i = symmetry->ind_i  + oppshift1;
            p_ind_j = symmetry->ind_j  + oppshift1;
            p_rot   = symmetry->rot    + oppshift1;
            for (i4 = 0; i4 < *p_num; i4++) {
              if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * \
              (double) symmetry->character_table[s * symmetry->number_of_classes + i2];
              p_ind_i++;
              p_ind_j++;
              p_rot++;
             }
            op++;
           } // close loop on i3
          } // close loop on i2
            norm = k_zero;
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) norm += vector[i1] * vector[i1];
              if (fabs(norm) > 0.00001) {
                for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector[i1] /= sqrt(norm); 
               }
        for (op = 0; op < symmetry->number_of_operators; op++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector_rot[i1] = k_zero;
          for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
            atm_rot = atom_p->K[(ip1 + i1) * symmetry->number_of_operators + op];
            offset = (atm_rot - ip1) * shells->shar[index_i];
            p_ind_i = symmetry->ind_i  + oppshift1;
            p_ind_j = symmetry->ind_j  + oppshift1;
            p_rot   = symmetry->rot    + oppshift1;
            for (i2 = 0; i2 < *p_num; i2++) {
              vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
              p_ind_i++;
              p_ind_j++;
              p_rot++;
             } // close loop on i2
            } // close loop on i1
            norm = k_zero;
            // Gram-Schmidt orthogonalise
            for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) vector_new[num_irp][i2] = k_zero;
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
              for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
                projection_operator[i1][i2] = k_zero;
               }
                projection_operator[i1][i1] = k_one;
               }
            for (i1 = 0; i1 < num_irp; i1++) {
              for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
                 for (i3 = 0; i3 < atom_p->numb[atm] * shells->shar[index_i]; i3++) {
                   projection_operator[i2][i3] -= vector_new[i1][i2] * vector_new[i1][i3];
                  }
                 }
                }
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
              for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
                vector_new[num_irp][i1] += projection_operator[i1][i2] * vector_rot[i2];
               }
              }
            norm = k_zero;
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) 
            norm += vector_new[num_irp][i1] * vector_new[num_irp][i1];
            if (norm > 0.0000001) {
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector_new[num_irp][i1] /= sqrt(norm); 
            //fprintf(file.out,"%3d %3d %3d\n",i,s,index_i);
            //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
            //fprintf(file.out,"%7.3lf",vector_new[num_irp][i2]); } fprintf(file.out,"\n"); 
            fprintf(file.out,"numirp %3d\n",num_irp);
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) total_coef++;
            num_irp++;
            num_salc++;
           } // close if (norm
          } // close loop on op
            //fprintf(file.out,"\n");
           } // close loop on i
         } // close index_i
         //if (num_salc == 0) num_salc = 1;
         //if (total_coef == 0) total_coef = 1;
         salc->num_salc = num_salc;
         salc->total_coef = total_coef;

}
*/

/*
void generate_atom_salc(int s, int ip, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4;
int num_coef, total_coef, total_atom_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int count;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
int offset, atm, atm_rot, index_i, shelposi, bfposi;

  total_coef = 0;
  total_atom_coef = 0;
  num_salc = 0;
  atm = atoms->uniq[ip];
int ip1;
ip = atm;
ip1 = atom_p->posn[atm];
  salc->num_atom = atom_p->numb[atm];
  shelposi = atoms->shelposn_sh[ip1];
  bfposi = 0;
  salc->num_irp[s] = 0;
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip1]; index_i++) {
    double vector[atom_p->numb[atm] * shells->shar[index_i]];
    double vector_rot[atom_p->numb[atm] * shells->shar[index_i]];
    double vector_new[atom_p->numb[atm] * shells->shar[index_i] * symmetry->cls_num_k[s]][atom_p->numb[atm]*shells->shar[index_i]];
    double projection_operator[atom_p->numb[atm] * shells->shar[index_i]][atom_p->numb[atm] * shells->shar[index_i]];
    num_irp = 0;
    for (i = 0; i < shells->shar[index_i]; i++) {
      for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector[i1] = k_zero;
        op = 0;
        for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
          for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
            oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
            p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
            atm_rot = atom_p->K[ip1 * symmetry->number_of_operators + op];
            offset = (atm_rot - ip1) * shells->shar[index_i]; 
            p_ind_i = symmetry->ind_i  + oppshift1;
            p_ind_j = symmetry->ind_j  + oppshift1;
            p_rot   = symmetry->rot    + oppshift1;
            for (i4 = 0; i4 < *p_num; i4++) {
              if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * \
              (double) symmetry->character_table[s * symmetry->number_of_classes + i2];
              p_ind_i++;
              p_ind_j++;
              p_rot++;
             }
            op++;
           } // close loop on i3
          } // close loop on i2
            norm = k_zero;
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) norm += vector[i1] * vector[i1];
              if (fabs(norm) > 0.00001) {
                for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector[i1] /= sqrt(norm); 
               }
        for (op = 0; op < symmetry->number_of_operators; op++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector_rot[i1] = k_zero;
          for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
            atm_rot = atom_p->K[(ip1 + i1) * symmetry->number_of_operators + op];
            offset = (atm_rot - ip1) * shells->shar[index_i];
            p_ind_i = symmetry->ind_i  + oppshift1;
            p_ind_j = symmetry->ind_j  + oppshift1;
            p_rot   = symmetry->rot    + oppshift1;
            for (i2 = 0; i2 < *p_num; i2++) {
              vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
              p_ind_i++;
              p_ind_j++;
              p_rot++;
             } // close loop on i2
            } // close loop on i1
            norm = k_zero;
            // Gram-Schmidt orthogonalise
            for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) vector_new[num_irp][i2] = k_zero;
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
              for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
                projection_operator[i1][i2] = k_zero;
               }
                projection_operator[i1][i1] = k_one;
               }
            for (i1 = 0; i1 < num_irp; i1++) {
              for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
                 for (i3 = 0; i3 < atom_p->numb[atm] * shells->shar[index_i]; i3++) {
                   projection_operator[i2][i3] -= vector_new[i1][i2] * vector_new[i1][i3];
                  }
                 }
                }
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
              for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
                vector_new[num_irp][i1] += projection_operator[i1][i2] * vector_rot[i2];
               }
              }
            norm = k_zero;
            num_coef = 0;
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) 
            norm += vector_new[num_irp][i1] * vector_new[num_irp][i1];
            if (norm > 0.0000001) {
            for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector_new[num_irp][i1] /= sqrt(norm); 
            //fprintf(file.out,"%3d %3d %3d\n",i,s,index_i);
            //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
            //fprintf(file.out,"%7.3lf",vector_new[num_irp][i2]); } fprintf(file.out,"\n"); 
            for (i2 = 0; i2 < shells->shar[index_i]; i2++) {
              for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
                salc->coeff->a[i1][total_atom_coef]  = vector_new[num_irp][i1 * shells->shar[index_i] + i2];  // FIX num_salc
                salc->bfn_posn->a[i1][total_atom_coef]  = bfposi + i2;
                total_coef++;
               }
                total_atom_coef++;
                num_coef++;
               }
                salc->num_coef[num_salc] = num_coef;
               (salc->num_irp[s])++;
                salc->irp[num_salc] = s;
                num_irp++;
                num_salc++;
               } // close if (norm
              } // close loop on op
             } // close loop on i
              //fprintf(file.out,"\n");
           bfposi += shells->shar[index_i];
          } // close index_i
         salc->num_salc = num_salc;
         salc->total_coef = total_coef;

}
*/
/*
void count_shell_salc(int index_i, int s, int atm, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4;
int num_coef, total_coef, total_atom_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int count;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
int offset, atm_rot, shelposi, bfposi;
int atm1;
atm1 = atom_p->posn[atm];

  total_coef = 0;
  total_atom_coef = 0;
  num_salc = 0;
  salc->num_atom = atom_p->numb[atm];
  shelposi = atoms->shelposn[atm];
  bfposi = 0;
  double vector[atom_p->numb[atm] * shells->shar[index_i]];
  double vector_rot[atom_p->numb[atm] * shells->shar[index_i]];
  double vector_new[atom_p->numb[atm] * shells->shar[index_i] * symmetry->cls_num_k[s]][atom_p->numb[atm]*shells->shar[index_i]];
  double projection_operator[atom_p->numb[atm] * shells->shar[index_i]][atom_p->numb[atm] * shells->shar[index_i]];
  num_irp = 0;
  for (i = 0; i < shells->shar[index_i]; i++) {
    for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector[i1] = 0;
      op = 0;
      for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
        for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          atm_rot = atom_p->K[atm1 * symmetry->number_of_operators + op];
          offset = (atm_rot - atm1) * shells->shar[index_i]; // check atm_rot or ip??
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i4 = 0; i4 < *p_num; i4++) {
            if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * \
            (double) symmetry->character_table[s * symmetry->number_of_classes + i2];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           }
          op++;
         } // close loop on i3
        } // close loop on i2
  norm = k_zero;
  for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) norm += vector[i1] * vector[i1];
    if (fabs(norm) > 0.00001) {
      for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector[i1] /= sqrt(norm); 
     }
      for (op = 0; op < symmetry->number_of_operators; op++) {
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
        p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
        for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector_rot[i1] = k_zero;
        for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
          atm_rot = atom_p->K[(atm1 + i1) * symmetry->number_of_operators + op];
          offset = (atm_rot - atm1) * shells->shar[index_i];
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i2 = 0; i2 < *p_num; i2++) {
            vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           } // close loop on i2
          } // close loop on i1
          norm = k_zero;
          // Gram-Schmidt orthogonalise
          for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) vector_new[num_irp][i2] = k_zero;
          for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
            for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
              projection_operator[i1][i2] = k_zero;
             }
              projection_operator[i1][i1] = k_one;
             }
           for (i1 = 0; i1 < num_irp; i1++) {
             for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
               for (i3 = 0; i3 < atom_p->numb[atm] * shells->shar[index_i]; i3++) {
                 projection_operator[i2][i3] -= vector_new[i1][i2] * vector_new[i1][i3];
                }
               }
              }
           for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
             for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
               vector_new[num_irp][i1] += projection_operator[i1][i2] * vector_rot[i2];
             }
            }
          norm = k_zero;
          num_coef = 0;
          for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) 
          norm += vector_new[num_irp][i1] * vector_new[num_irp][i1];
          if (norm > 0.0000001) {
          for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector_new[num_irp][i1] /= sqrt(norm); 
          //fprintf(file.out,"%3d %3d %3d\n",i,s,index_i);
          //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
          //fprintf(file.out,"%7.3lf",vector_new[num_irp][i2]); } fprintf(file.out,"\n"); 
          for (i2 = 0; i2 < shells->shar[index_i]; i2++) {
          for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
            total_coef++;
           }
            total_atom_coef++;
            num_coef++;
           }
            num_irp++;
            num_salc++;
          } // close if (norm
         } // close loop on op
        } // close loop on i
       salc->num_salc = num_salc;
       salc->total_coef = total_coef;

}

void generate_shell_salc(int index_i, int s, int atm, SALC *salc, ATOM *atoms, ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int num_salc;
int num_irp;
int i1, i2, i3, i4;
int num_coef, total_coef, total_atom_coef;
int i, j, k, l, m, op;
int oppshift1, sheli1;
int count;
int *p_ind_i, *p_ind_j, *num_i, *p_num, *p_ost;
double *p_rot, *p_rot1;
double norm;
int offset, atm_rot, shelposi, bfposi;
int atm1;
atm1 = atom_p->posn[atm];

  total_coef = 0;
  total_atom_coef = 0;
  num_salc = 0;
  salc->num_atom = atom_p->numb[atm];
  shelposi = atoms->shelposn_sh[atm];
  bfposi = 0;
  salc->num_irp[s] = 0;
  double vector[atom_p->numb[atm] * shells->shar[index_i]];
  double vector_rot[atom_p->numb[atm] * shells->shar[index_i]];
  double vector_new[atom_p->numb[atm] * shells->shar[index_i] * symmetry->cls_num_k[s]][atom_p->numb[atm]*shells->shar[index_i]];
  double projection_operator[atom_p->numb[atm] * shells->shar[index_i]][atom_p->numb[atm] * shells->shar[index_i]];
  num_irp = 0;
  for (i = 0; i < shells->shar[index_i]; i++) {
    for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector[i1] = 0;
      op = 0;
      for (i2 = 0; i2 < symmetry->number_of_classes; i2++) {
        for (i3 = 0; i3 < symmetry->cls_num_k[i2]; i3++) {
          oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
          p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
          atm_rot = atom_p->K[atm1 * symmetry->number_of_operators + op];
          offset = (atm_rot - atm1) * shells->shar[index_i]; 
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i4 = 0; i4 < *p_num; i4++) {
            if (i == *p_ind_j) vector[offset + *p_ind_i] += *p_rot * (double) symmetry->character_table[s * \
            symmetry->number_of_classes + i2];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           }
          op++;
         } // close loop on i3
        } // close loop on i2
  norm = k_zero;
  for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) norm += vector[i1] * vector[i1];
    if (fabs(norm) > 0.00001) {
      for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector[i1] /= sqrt(norm); 
     }
      for (op = 0; op < symmetry->number_of_operators; op++) {
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + shells->ord_sh[index_i]);
        p_num = symmetry->num_ij + op * (job->l_max + 2) + shells->ord_sh[index_i];
        for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i];i1++) vector_rot[i1] = k_zero;
        for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
          atm_rot = atom_p->K[(atm1 + i1) * symmetry->number_of_operators + op];
          offset = (atm_rot - atm1) * shells->shar[index_i];
          p_ind_i = symmetry->ind_i  + oppshift1;
          p_ind_j = symmetry->ind_j  + oppshift1;
          p_rot   = symmetry->rot    + oppshift1;
          for (i2 = 0; i2 < *p_num; i2++) {
            vector_rot[offset + *p_ind_i] += *p_rot * vector[i1 * shells->shar[index_i] + *p_ind_j];
            p_ind_i++;
            p_ind_j++;
            p_rot++;
           } // close loop on i2
          } // close loop on i1
          norm = k_zero;
          // Gram-Schmidt orthogonalise
          for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) vector_new[num_irp][i2] = k_zero;
          for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
            for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
              projection_operator[i1][i2] = k_zero;
             }
              projection_operator[i1][i1] = k_one;
             }
           for (i1 = 0; i1 < num_irp; i1++) {
             for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
               for (i3 = 0; i3 < atom_p->numb[atm] * shells->shar[index_i]; i3++) {
                 projection_operator[i2][i3] -= vector_new[i1][i2] * vector_new[i1][i3];
                }
               }
              }
           for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) {
             for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
               vector_new[num_irp][i1] += projection_operator[i1][i2] * vector_rot[i2];
             }
            }
          norm = k_zero;
          num_coef = 0;
          for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) 
          norm += vector_new[num_irp][i1] * vector_new[num_irp][i1];
          if (norm > 0.0000001) {
          for (i1 = 0; i1 < atom_p->numb[atm] * shells->shar[index_i]; i1++) vector_new[num_irp][i1] /= sqrt(norm); 
          //fprintf(file.out,"%3d %3d %3d\n",i,s,index_i);
          //for (i2 = 0; i2 < atom_p->numb[atm] * shells->shar[index_i]; i2++) {
          //fprintf(file.out,"%7.3lf",vector_new[num_irp][i2]); } fprintf(file.out,"\n"); 
          for (i2 = 0; i2 < shells->shar[index_i]; i2++) {
          for (i1 = 0; i1 < atom_p->numb[atm]; i1++) {
            salc->atm[total_coef] = i1; 
            salc->coeff->a[i1][total_atom_coef]  = vector_new[num_irp][i1 * shells->shar[index_i] + i2];  // FIX num_salc
            salc->bfn_posn->a[i1][total_atom_coef]  = bfposi + i2;
            total_coef++;
           }
            total_atom_coef++;
            num_coef++;
           }
            salc->num_coef[num_salc] = num_coef;
           (salc->num_irp[s])++;
            salc->irp[num_salc] = s;
            num_irp++;
            num_salc++;
          } // close if (norm
         } // close loop on op
        } // close loop on i
       //fprintf(file.out,"\n");
       salc->num_salc = num_salc;
       salc->total_coef = total_coef;

}
*/

	/*
  #include "SYMMETRY_ADAPTATION.h"
  SALC salc1;
  int r, t, i1, count_j, atm1, num_irrep_in_basis[symmetry->number_of_classes];
  double dot_product;
  DoubleMatrix *salc_rows, *salc_cols;
  IntMatrix *num_irrep_in_atom, *salc_count1;

  AllocateIntMatrix(&num_irrep_in_atom, &symmetry->number_of_classes, &atoms->number_of_unique_atoms, job);
  AllocateIntMatrix(&salc_count1, &symmetry->number_of_classes, &atoms->number_of_unique_atoms, job);
  ResetIntMatrix(salc_count1);
  count_atom_irrep(num_irrep_in_atom,atoms,atom_p,shells,symmetry,job,file);
  count_basis_irrep(num_irrep_in_basis,atoms,atom_p,shells,symmetry,job,file);

  if (job->verbosity >= 1)
  for (atm1 = 0; atm1 < atoms->number_of_unique_atoms; atm1++) {
    for (s = 0; s < symmetry->number_of_classes; s++) {
      count_atom_salc(s,atm1,&salc1,atoms,atom_p,pair_p,shells,symmetry,job,file);
      if (salc1.num_salc == 0) continue;
      allocate_SALC(&salc1,symmetry,job,file);
      fprintf(file.out,"\nsalc1 %3d %3d\n",atm1,s);
      generate_atom_salc(s,atm1,&salc1,atoms,atom_p,pair_p,shells,symmetry,job,file);
      bfposi = atoms->bfnposn_sh[atm1];
      for (k = 0; k < occupied[0]; k++) {
        count_j = 0;
        dot_product = k_zero;
        for (r = 0; r < salc1.num_irp[s]; r++){
          for (t = 0; t < salc1.num_coef[r]; t++){
            for (i1 = 0; i1 < atom_p->numb[atm1]; i1++) {
              fprintf(file.out,"wfn %3d irp %3d i1 %3d t %3d count %3d %3d %14.8f %10.4f\n",\
              k,r,i1,t,count_j+t,salc1.bfn_posn->a[i1][count_j+t],salc1.coeff->a[i1][count_j+t],eigenvectors->a[k][salc1.bfn_posn->a[k][count_j+t]]);
              if (i1 == 0) dot_product += salc1.coeff->a[i1][count_j+t] * eigenvectors->a[k][salc1.bfn_posn->a[k][count_j+t]];
             }
            }
           fprintf(file.out,"\n");
           count_j += salc1.num_coef[r];
          }
         fprintf(file.out,"dot product %3d %14.8f\n",k,dot_product);
        }
       free_SALC(&salc1,job);
      }
     }
  */

