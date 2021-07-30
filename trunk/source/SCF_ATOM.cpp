

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
#include "LINEAR_ALGEBRA_UTIL.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "INCOMPLETE_GAMMA.h"
#include "DENSITY_MATRIX_MOLECULE.h"
#include "SCF_ATOM.h"

using namespace std;

void atom_scf(ATOM *atoms, int atm, DoubleMatrix *eigvec, double *eigval, double *occupation, int occupied[2], SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int dim1, dim2, dim3, dim4;
  int iter;
  int i, j, k, s;
  int offset;
  int mix_hist = 12;
  int info = 0;
  double total_energy = k_zero;
  double alpha = k_one, beta = k_zero;
  double mix, tol = 1e-07, total_energy_change;
  double d_max;
  double time1, time2;
  double *eigenvalues;
  char trans = 'T', no_trans = 'N';
  char uplo = 'U';
  char jobz = 'V';
  DoubleMatrix *Over, *Overlap, *Kinetic, *Fock, *F_1e, *F_2e, *P_red, *P_rd1;
  DoubleMatrix *K_2e;
  DoubleMatrix *eigenvectors, *eigenvector1, *xtrn, *xtmp;
  DoubleMatrix *C_i, *C_o, *D;
  INTEGRAL_LIST integral_list;

    dim1 = atoms->bfnnumb_sh[atm];
    dim2 = dim1 * atoms->bfnnumb_sh[atm];
    dim3 = atoms->magnetic[atm] * dim1;
    dim4 = dim2 * dim2;
    allocate_integral_list(&integral_list,dim4,job,file);
    AllocateDoubleMatrix(&eigenvectors,&dim1,&dim1,job);
    AllocateDoubleMatrix(&eigenvector1,&dim1,&dim1,job);
    AllocateDoubleArray(&eigenvalues,&dim1,job);
    AllocateDoubleMatrix(&C_i,&mix_hist,&dim2,job);
    AllocateDoubleMatrix(&C_o,&mix_hist,&dim2,job);
    AllocateDoubleMatrix(&D,&mix_hist,&dim2,job);
    AllocateDoubleMatrix(&Overlap,&dim1,&dim1,job);
    AllocateDoubleMatrix(&Over,&dim1,&dim1,job);
    AllocateDoubleMatrix(&F_1e,&dim3,&dim1,job);
    AllocateDoubleMatrix(&F_2e,&dim3,&dim1,job);
    AllocateDoubleMatrix(&K_2e,&dim3,&dim1,job);
    AllocateDoubleMatrix(&Fock,&dim1,&dim1,job);
    AllocateDoubleMatrix(&P_red,&dim3,&dim1,job);
    AllocateDoubleMatrix(&P_rd1,&dim3,&dim1,job);
    AllocateDoubleMatrix(&xtrn,&dim1,&dim1,job);
    AllocateDoubleMatrix(&xtmp,&dim1,&dim1,job);
    ResetDoubleMatrix(eigenvectors);
    ResetDoubleMatrix(eigenvector1);
    ResetDoubleMatrix(C_i);
    ResetDoubleMatrix(C_o);
    ResetDoubleMatrix(D);
    ResetDoubleMatrix(Overlap);
    ResetDoubleMatrix(P_red);
    ResetDoubleMatrix(P_rd1);
    ResetDoubleMatrix(F_2e);
    ResetDoubleMatrix(K_2e);
    ResetDoubleMatrix(F_1e);
    ResetDoubleMatrix(Fock);
    ResetDoubleMatrix(xtrn);
    ResetDoubleMatrix(xtmp);

    initial_density_matrix_atom(P_red, atoms, &atm, occupation, occupied, shells, job, file);
 
    fock_element_1e_atm(F_1e, Overlap, &atm, atoms, shells, gaussians, job, file);
    time1 = MPI_Wtime();

    integral_list.num = 0;
    integrals_2e_atom(&integral_list, atm, atoms, shells, gaussians, job, file);
    time2 = MPI_Wtime();

    if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"Time for %d atom integrals %lf\n",integral_list.num,time2 - time1);
    for(i=0;i<integral_list.num;i++){
    fprintf(file.out,"dint list %d %d %d %d %e\n",integral_list.i[i],integral_list.j[i],integral_list.k[i],\
    integral_list.l[i],integral_list.value[i]);
    }
    fprintf(file.out,"Overlap\n");
    //PrintDoubleMatrix(Overlap, file);
    }

    for (i = 0; i < dim1; i++) {
      for (j = 0; j < dim1; j++) {
        Over->a[i][j] = Overlap->a[i][j];
      }
    }

    DiagonaliseSymmetrical(&Overlap, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);

    for (i = 0; i < dim1; i++) {
      for (j = 0; j < dim1; j++) {
        xtrn->a[i][j] = eigenvectors->a[i][j] / sqrt(*(eigenvalues + i));
      }
    }

    for (iter = 1; iter < 50; iter++) {

    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out,"beginning iteration %d\n",iter);
    fflush(file.out);

    ResetDoubleMatrix(F_2e);

    if (job->xc_num == 0) 
    fock_2e_matrix(F_2e, &integral_list, P_red, atm, atoms, job, file);

    if (job->xc_num > 0) 
    ////kohn_sham_2e_matrix(&atm, F_2e, &integral_list, P_red, atoms, shells, gaussians, job, file);

    //CHANGES2017 THIS ROUTINE IS CAUSING PROBELMS IN INITIAL DENSITY (CALC OF TOTAL ENEGY AND TOTAL ENERGY CHANGE
    ////total_energy_atom(F_1e, F_2e, P_red, atm, &total_energy, &total_energy_change, atoms, job, file);

    ResetDoubleMatrix(P_rd1);

    for (s = 0; s < atoms->magnetic[atm]; s++) {

    offset = s * atoms->bfnnumb_sh[atm];

    for (i = 0; i < dim1; i++) {
      for (j = 0; j < dim1; j++) {
        //Fock->a[i][j] = F_1e->a[i][j] + K_2e->a[offset + i][j] + F_2e->a[offset + i][j];
        Fock->a[i][j] = F_1e->a[i][j] + F_2e->a[offset + i][j];
       }
      }

    DoubleGEMM(&no_trans, &no_trans, &alpha, &xtrn, &Fock, &beta, &xtmp);
    ResetDoubleMatrix(Fock);
    DoubleGEMM(&no_trans, &trans,    &alpha, &xtmp, &xtrn, &beta, &Fock);

    DiagonaliseSymmetrical(&Fock, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);

    ResetDoubleMatrix(eigenvector1);
    DoubleGEMM(&no_trans, &no_trans, &alpha, &eigenvectors, &xtrn, &beta, &eigenvector1);

    for (i = 0; i < dim1; i++) {
      for (j = 0; j < dim1; j++) {
        for (k = 0; k < dim1; k++) {
          P_rd1->a[offset + i][j] += occupation[k] / atoms->magnetic[atm] * eigenvector1->a[k][i] * eigenvector1->a[k][j];
         }
        }
       }

    for (i = 0; i < dim1; i++) {
      for (j = 0; j < dim1; j++) {
        eigvec->a[offset + i][j] = eigenvector1->a[i][j];
       }
      }

    for (i = 0; i < dim1; i++) {
        eigval[offset + i] = eigenvalues[i];
      }

    } // close loop over up/down spins

    mix = 0.8;
    mix_hist = 5;

    //mix_density(C_i, C_o, D, P_red, &mix, &residue, &d_max, iter, &mix_hist, job, file);

    for (i = 0; i < atoms->magnetic[atm] * dim1; i++) {
      for (j = 0; j < dim1; j++) {
          P_red->a[i][j] = P_red->a[i][j] * mix + P_rd1->a[i][j] * (k_one - mix);
         }
        }
    
    if (fabs(total_energy_change) < tol) {
      ResetDoubleMatrix(F_2e);
      ////if (job->xc_num > 0) 
        ////kohn_sham_2e_matrix(&atm, F_2e, &integral_list, P_red, atoms, shells, gaussians, job, file);
      if (job->xc_num == 0) 
        fock_2e_matrix(F_2e, &integral_list, P_red, atm, atoms, job, file);
          if (job->taskid == 0)
          final_total_energy_atom(F_1e, F_2e, P_red, atm, &total_energy, &total_energy_change, atoms, job, file);
            break;
           }

  } // close loop on iter

    DestroyDoubleMatrix(&eigenvectors,job);
    DestroyDoubleMatrix(&eigenvector1,job);
    DestroyDoubleMatrix(&C_i,job);
    DestroyDoubleMatrix(&C_o,job);
    DestroyDoubleMatrix(&D,job);
    DestroyDoubleMatrix(&Overlap,job);
    DestroyDoubleMatrix(&Over,job);
    DestroyDoubleMatrix(&P_red,job);
    DestroyDoubleMatrix(&P_rd1,job);
    DestroyDoubleMatrix(&F_1e,job);
    DestroyDoubleMatrix(&F_2e,job);
    DestroyDoubleMatrix(&K_2e,job);
    DestroyDoubleMatrix(&Fock,job);
    DestroyDoubleMatrix(&xtrn,job);
    DestroyDoubleMatrix(&xtmp,job);
    free_integral_list(&integral_list,job);
    DestroyDoubleArray(&eigenvalues,&dim1,job);

}

void total_energy_atom(DoubleMatrix *F_1e, DoubleMatrix *F_2e, DoubleMatrix *P_red, int atm, double *current_total_energy, double *total_energy_change, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int i, j;
double oneelec_energy = k_zero;
double twoelec_energy = k_zero;
double kinetic_energy = k_zero;
double fock_energy = k_zero;
double coulomb_energy = k_zero;
double exchange_energy = k_zero;

  for (i = 0; i < atoms->magnetic[atm] * atoms->bfnnumb_sh[atm]; i++) {
   for (j = 0; j < atoms->bfnnumb_sh[atm]; j++) {
    oneelec_energy += F_1e->a[i][j] * P_red->a[i][j];
   }
  }

  for (i = 0; i < atoms->magnetic[atm] * atoms->bfnnumb_sh[atm]; i++) {
   for (j = 0; j < atoms->bfnnumb_sh[atm]; j++) {
    twoelec_energy += F_2e->a[i][j] * P_red->a[i][j];
   }
  }

  *total_energy_change  = (oneelec_energy + twoelec_energy / two) - *current_total_energy;
  *current_total_energy = (oneelec_energy + twoelec_energy / two);

  if (job->taskid == 0 && job->verbosity >= 1) {
  fprintf(file.out, "one electron energy     %16.8e\n", oneelec_energy);
  fprintf(file.out, "two electron energy     %16.8e\n", twoelec_energy / two);
  fprintf(file.out, "total energy            %16.8e\n", *current_total_energy);
  fprintf(file.out, "total energy change     %16.8e\n", *total_energy_change);
  fprintf(file.out,"\n");
  printf("one electron energy     %16.8e\n", oneelec_energy);
  printf("two electron energy     %16.8e\n", twoelec_energy / two);
  printf("total energy            %16.8e\n", *current_total_energy);
  printf("total energy change     %16.8e\n", *total_energy_change);
  printf("\n");
 }

}

void final_total_energy_atom(DoubleMatrix *F_1e, DoubleMatrix *F_2e, DoubleMatrix *P_red, int atm, double *current_total_energy, double *total_energy_change, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int i, j;
double oneelec_energy = k_zero;
double twoelec_energy = k_zero;

  for (i = 0; i < atoms->magnetic[atm] * atoms->bfnnumb_sh[atm]; i++) {
   for (j = 0; j < atoms->bfnnumb_sh[atm]; j++) {
    oneelec_energy += F_1e->a[i][j] * P_red->a[i][j];
   }
  }

  for (i = 0; i < atoms->magnetic[atm] * atoms->bfnnumb_sh[atm]; i++) {
   for (j = 0; j < atoms->bfnnumb_sh[atm]; j++) {
    twoelec_energy += F_2e->a[i][j] * P_red->a[i][j];
   }
  }

  *total_energy_change  = (oneelec_energy + twoelec_energy / two) - *current_total_energy;
  *current_total_energy = (oneelec_energy + twoelec_energy / two);

  if (job->taskid == 0) {
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"| ATOM CALCULATION          | ATOMIC NUMBER        %2d | TOTAL ENERGY %10.2e | ENERGY CHANGE %9.2e |\n", \
  atoms->atomic_number[atm],*current_total_energy,*total_energy_change);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 }

}

void initial_density_matrix_atom(DoubleMatrix *P_red, ATOM *atoms, int *atm, double *occupation, int occupied[2], SHELL *shells, JOB_PARAM *job, FILES file)

{

int i, j, k;
int count;
double *p_occupation, *p_occupation_up, *p_occupation_down, total_electron_count, electron_count_up, electron_count_down;

  switch (atoms->spin[*atm]) {

    case 0: //  atom->spin switch

     p_occupation = occupation;

      for (k = 0; k < atoms->bfnnumb_sh[*atm]; k++) {
       *p_occupation = k_zero;
        p_occupation++;
       }
        p_occupation = occupation;

    for (j = 0; j < atoms->nshel[*atm]; j++) {

      if (shells->nele[atoms->shelposn[*atm] + j] > k_zero) {

      switch (shells->type1[atoms->shelposn[*atm] + j]) {
           
      case 1:    // s shell

        if (shells->nele[atoms->shelposn[*atm] + j] <= two) {
           *p_occupation = shells->nele[atoms->shelposn[*atm] + j];
            p_occupation++;
           }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in s shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 3:     // p shell
        if (shells->nele[atoms->shelposn[*atm] + j] <= six) {
          for (k = 0; k < 3; k++) {
           *p_occupation = shells->nele[atoms->shelposn[*atm] + j] / three;
            p_occupation++;
           }
          }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in p shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 4:     // sp shell
        if (shells->nele[atoms->shelposn[*atm] + j] <= two) {
          *p_occupation = shells->nele[atoms->shelposn[*atm] + j];
           p_occupation++;
          for (k = 0; k < 3; k++) {
           *p_occupation = k_zero;
            p_occupation++;
           }
          }
        else if (shells->nele[atoms->shelposn[*atm] + j] > two && shells->nele[atoms->shelposn[*atm] + j] <= eight) {
          *p_occupation = two;
           p_occupation++;
          for (k = 0; k < 3; k++) {
           *p_occupation = (shells->nele[atoms->shelposn[*atm] + j] - two) / three;
            p_occupation++;
           }
          }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in sp shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 5:    // d shell
        if (shells->nele[atoms->shelposn[*atm] + j] <= ten) {
          for (k = 0; k < 5; k++) {
          *p_occupation = shells->nele[atoms->shelposn[*atm] + j] / five;
           p_occupation++;
           }
          }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in d shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 7:    // f shell
        if (shells->nele[atoms->shelposn[*atm] + j] <= fourteen) {
          for (k = 0; k < 7; k++) {
          *p_occupation = shells->nele[atoms->shelposn[*atm] + j] / seven;
           p_occupation++;
           }
          }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in f shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 9:    // g shell
        if (shells->nele[atoms->shelposn[*atm] + j] <= eighteen) {
          for (k = 0; k < 9; k++) {
          *p_occupation = shells->nele[atoms->shelposn[*atm] + j] / nine;
           p_occupation++;
           }
          }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in g shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

       } // close switch (shells->type1

      } // close if
   
      else { p_occupation +=  shells->type1[atoms->shelposn[*atm] + j]; } 

     } // close loop on j

    if (job->taskid == 0 && job->verbosity > 1) {
    count = 0;
    for (j = 0; j < atoms->nshel[*atm]; j++) {
     for (k = 0; k < shells->type1[atoms->shelposn[*atm] + j]; k++) {
       fprintf(file.out,"Occupation %d %d %lf\n",j,k,*(occupation + count));
         count++;
        }
       }
      }

 count = 0;
 total_electron_count = k_zero;
 for (i = 0; i < atoms->bfnnumb_sh[*atm]; i++) {
   P_red->a[i][i] = *(occupation + count);
   total_electron_count += *(occupation + count);
   count++;
  }
  
 occupied[0] = int((total_electron_count + 0.00001) / two);
 if (total_electron_count - (double)occupied[0] * two > 0.00001) 
 (occupied[0])++;              // corrects case of non-magnetic atoms with odd number of electrons
  occupied[1] = occupied[0];

         break;

    case -7 ... -1:   // atoms->spin switch
    case +1 ... +7: 

      p_occupation_up   = occupation;
      p_occupation_down = occupation + atoms->bfnnumb_sh[*atm];

      for (k = 0; k < atoms->bfnnumb_sh[*atm]; k++) {
       *p_occupation_up   = k_zero;
       *p_occupation_down = k_zero;
        p_occupation_up++;
        p_occupation_down++;
       }
      p_occupation_up   = occupation;
      p_occupation_down = occupation + atoms->bfnnumb_sh[*atm];

    for (j = 0; j < atoms->nshel[*atm]; j++) {

      if (shells->nele[atoms->shelposn[*atm] + j] > k_zero) {

  switch (shells->type1[atoms->shelposn[*atm] + j]) {
           
      case 1:   // s shell
        if (shells->nele[atoms->shelposn[*atm] + j] <= k_one) {
           *p_occupation_up   = shells->nele[atoms->shelposn[*atm] + j];
           *p_occupation_down = k_zero;
            p_occupation_up++;
            p_occupation_down++;
         } 
        else if (shells->nele[atoms->shelposn[*atm] + j] > k_one && shells->nele[atoms->shelposn[*atm] + j] <= two) {
           *p_occupation_up   = shells->nele[atoms->shelposn[*atm] + j] / two;
           *p_occupation_down = shells->nele[atoms->shelposn[*atm] + j] / two;
            p_occupation_up++;
            p_occupation_down++;
           }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in s shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 3:    // p shell
        if (shells->nele[atoms->shelposn[*atm] + j] <= three) {
          for (k = 0; k < 3; k++) {
           *p_occupation_up   = shells->nele[atoms->shelposn[*atm] + j] / three;
           *p_occupation_down = k_zero;
            p_occupation_up++;
            p_occupation_down++;
           }
          }
        else if (shells->nele[atoms->shelposn[*atm] + j] > three && shells->nele[atoms->shelposn[*atm] + j] <= six) {
          for (k = 0; k < 3; k++) {
           *p_occupation_up   = shells->nele[atoms->shelposn[*atm] + j] / six;
           *p_occupation_down = shells->nele[atoms->shelposn[*atm] + j] / six;
            p_occupation_up++;
            p_occupation_down++;
           }
          }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in p shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 4:      // sp shell
        //if (shells->nele[atoms->shelposn[*atm] + j] <= k_one) {
        if (shells->nele[atoms->shelposn[*atm] + j] <= two) {
          *p_occupation_up   = shells->nele[atoms->shelposn[*atm] + j];
          *p_occupation_down = k_zero;
           p_occupation_up++;
           p_occupation_down++;
          for (k = 0; k < 3; k++) {
           *p_occupation_up   = k_zero;
           *p_occupation_down = k_zero;
            p_occupation_up++;
            p_occupation_down++;
           }
          }
        else if (shells->nele[atoms->shelposn[*atm] + j] > k_one && shells->nele[atoms->shelposn[*atm] + j] <= two) {
          *p_occupation_up   = shells->nele[atoms->shelposn[*atm] + j] / two;
          *p_occupation_down = shells->nele[atoms->shelposn[*atm] + j] / two;
           p_occupation_up++;
           p_occupation_down++;
          for (k = 0; k < 3; k++) {
           *p_occupation_up   = k_zero;
           *p_occupation_down = k_zero;
            p_occupation_up++;
            p_occupation_down++;
           }
          }
        else if (shells->nele[atoms->shelposn[*atm] + j] > two && shells->nele[atoms->shelposn[*atm] + j] <= eight) {
          *p_occupation_up   = k_one;
          *p_occupation_down = k_one;
           p_occupation_up++;
           p_occupation_down++;
          for (k = 0; k < 3; k++) {
           *p_occupation_up   = (shells->nele[atoms->shelposn[*atm] + j] - two) / six;
           *p_occupation_down = (shells->nele[atoms->shelposn[*atm] + j] - two) / six;
            p_occupation_up++;
            p_occupation_down++;
           }
          }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in sp shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 5:  // d shell
        if (shells->nele[atoms->shelposn[*atm] + j] <= two) {
          *p_occupation_up = shells->nele[atoms->shelposn[*atm] + j] / two; // d_z^2 and d_x^2-y^2 occupation
           p_occupation_up++;
          *p_occupation_up = k_zero;
           p_occupation_up++;
          *p_occupation_up = k_zero;
           p_occupation_up++;
          *p_occupation_up = shells->nele[atoms->shelposn[*atm] + j] / two;
           p_occupation_up++;
          *p_occupation_up = k_zero;
           p_occupation_up++;
          for (k = 0; k < 5; k++) {
          *p_occupation_down = k_zero;
           p_occupation_down++;
           }
          }
        if (shells->nele[atoms->shelposn[*atm] + j] > two && shells->nele[atoms->shelposn[*atm] + j] <= three) {
          *p_occupation_up = k_zero;
           p_occupation_up++;
          *p_occupation_up = shells->nele[atoms->shelposn[*atm] + j] / three; // d_xz d_yz d_xy occupation
           p_occupation_up++;
          *p_occupation_up = shells->nele[atoms->shelposn[*atm] + j] / three;
           p_occupation_up++;
          *p_occupation_up = k_zero;
           p_occupation_up++;
          *p_occupation_up = shells->nele[atoms->shelposn[*atm] + j] / three;
           p_occupation_up++;
          for (k = 0; k < 5; k++) {
          *p_occupation_down = k_zero;
           p_occupation_down++;
           }
          }
        else if (shells->nele[atoms->shelposn[*atm] + j] > three && shells->nele[atoms->shelposn[*atm] + j] <= five) {
          *p_occupation_up = shells->nele[atoms->shelposn[*atm] + j] / five;
           p_occupation_up++;
          for (k = 0; k < 5; k++) {
          *p_occupation_down   = k_zero;
           p_occupation_down++;
           }
          }
        else if (shells->nele[atoms->shelposn[*atm] + j] > five && shells->nele[atoms->shelposn[*atm] + j] <= seven) {
          for (k = 0; k < 5; k++) {
          *p_occupation_up = k_one;
           p_occupation_up++;
          }
          *p_occupation_down = (shells->nele[atoms->shelposn[*atm] + j] - five) / two; // d_z^2 and d_x^2-y^2 occupation
           p_occupation_down++;
          *p_occupation_down = k_zero;
           p_occupation_down++;
          *p_occupation_down = k_zero;
           p_occupation_down++;
          *p_occupation_down = (shells->nele[atoms->shelposn[*atm] + j] - five) / two; 
           p_occupation_down++;
          *p_occupation_down = k_zero;
           p_occupation_down++;
         }
        else if (shells->nele[atoms->shelposn[*atm] + j] > seven && shells->nele[atoms->shelposn[*atm] + j] <= eight) {
          for (k = 0; k < 5; k++) {
          *p_occupation_up = k_one;
           p_occupation_up++;
          }
          *p_occupation_down = k_zero;
           p_occupation_down++;
          *p_occupation_down = (shells->nele[atoms->shelposn[*atm] + j] - five) / three; // d_xz d_yz d_xy occupation
           p_occupation_down++;
          *p_occupation_down = (shells->nele[atoms->shelposn[*atm] + j] - five) / three;
           p_occupation_down++;
          *p_occupation_down = k_zero;
           p_occupation_down++;
          *p_occupation_down = (shells->nele[atoms->shelposn[*atm] + j] - five)/ three;
           p_occupation_down++;
         }
        else if (shells->nele[atoms->shelposn[*atm] + j] > eight && shells->nele[atoms->shelposn[*atm] + j] <= ten) {
          for (k = 0; k < 5; k++) {
          *p_occupation_up = k_one;
           p_occupation_up++;
          *p_occupation_down = (shells->nele[atoms->shelposn[*atm] + j] - five) / five; 
           p_occupation_down++;
          }
         }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in d shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 7:
        if (shells->nele[atoms->shelposn[*atm] + j] <= fourteen) {
          for (k = 0; k < 7; k++) {
          *p_occupation_up   = shells->nele[atoms->shelposn[*atm] + j] / fourteen;
          *p_occupation_down = shells->nele[atoms->shelposn[*atm] + j] / fourteen;
           p_occupation_up++;
           p_occupation_down++;
           }
          }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in f shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

      case 9:
        if (shells->nele[atoms->shelposn[*atm] + j] <= eighteen) {
          for (k = 0; k < 9; k++) {
          *p_occupation_up   = shells->nele[atoms->shelposn[*atm] + j] / eighteen;
          *p_occupation_down = shells->nele[atoms->shelposn[*atm] + j] / eighteen;
           p_occupation_up++;
           p_occupation_down++;
           }
          }
        else {
           if (job->taskid == 0)
           fprintf(file.out,"Number of electrons %lf in g shell %d on atom %d is incorrect\n", \
           shells->nele[atoms->shelposn[*atm] + j],j,*atm);
           MPI_Finalize();
           exit(0);
          }
        break;

       } // close switch (shells->type1

      } // close if
   
      else {
            p_occupation_up   +=  shells->type1[atoms->shelposn[*atm] + j];  
            p_occupation_down +=  shells->type1[atoms->shelposn[*atm] + j]; 
           } 

     } // close loop on j

    if (job->taskid == 0 && job->verbosity > 1) {
    count = 0;
    for (j = 0; j < atoms->nshel[*atm]; j++) {
     for (k = 0; k < shells->type1[atoms->shelposn[*atm] + j]; k++) {
       fprintf(file.out,"Occupation %d %d %d %lf %lf\n",*atm,j,k,*(occupation + count),*(occupation + atoms->bfnnumb_sh[*atm] + count));
         count++;
        }
       }
      }

 count = 0;
 total_electron_count = k_zero;
 for (i = 0; i < atoms->magnetic[*atm]; i++) {
  electron_count_down = k_zero;
   for (j = 0; j < atoms->bfnnumb_sh[*atm]; j++) {
     P_red->a[i * atoms->bfnnumb_sh[*atm] + j][j] = *(occupation + count);
     total_electron_count += *(occupation + count);
     electron_count_down  += *(occupation + count);
     count++;
    }
   }

 electron_count_up = total_electron_count - electron_count_down;

 occupied[0] = int(electron_count_up   + 0.00001);
 occupied[1] = int(electron_count_down + 0.00001);

 if (occupied[0] - occupied[1] != abs(atoms->spin[*atm])) {
 fprintf(file.out," %3d %3d  %3d  %3d\n",\
 occupied[0], occupied[1],atoms->spin[*atm],*atm);
 fprintf(file.out,"Difference in number of up/down occupied states %3d incompatible with assigned spin %3d for atom %3d\n",\
 occupied[0] - occupied[1],atoms->spin[*atm],*atm);
 MPI_Finalize();
 exit(0);
}

         break;
 }

}

void fock_element_1e_atm(DoubleMatrix *F_1e, DoubleMatrix *Overlap, int *atm, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, mm;
  int ip, index_i, index_j, i4, j4;
  int bfposi, bfposj;
  int n1, n2, n3, n4, n5, n6;
  int t, u, v;
  int tmax, umax, vmax;
  int imax, jmax;
  int gausposi, gausposj, shelposi, shelposj;
  int bfposi1, bfposj1;
  int sheli, shelj;
  int dim, nd1;
  int nsheli, nshelj, sheli1, shelj1;
  int *p_i, *p_j, *p_i1, *p_j1;
  double SAB, ab, pab, p32, sabfac;
  double fn[55]; 
  double E1x[15][15][15], E1y[15][15][15], E1z[15][15][15];
  double *p_elecnuc, *p_overlap, *p_kinetic, *p_spin_orbit;
  double *kinetic, *elecnuc, *overlap, *spin_orbit;
  double *p_rot1, *p_rot2;

  ip = *atm;

   nd1 = atoms->bfnnumb[ip];
   dim = atoms->bfnnumb[ip] * atoms->bfnnumb[ip] ;

    kinetic = (double *) malloc(dim * sizeof(double));
    if (kinetic == NULL) { fprintf(stderr, "ERROR: not enough memory for overlap in\n"); exit(1); }
    p_kinetic = kinetic;
    for (i=0;i<dim;i++) { *p_kinetic = k_zero; p_kinetic++;}

    elecnuc = (double *) malloc(dim * sizeof(double));
    if (elecnuc == NULL) { fprintf(stderr, "ERROR: not enough memory for overlap in\n"); exit(1); }
    p_elecnuc = elecnuc;
    for (i=0;i<dim;i++) { *p_elecnuc = k_zero; p_elecnuc++;}

    spin_orbit = (double *) malloc(dim * sizeof(double));
    if (spin_orbit == NULL) { fprintf(stderr, "ERROR: not enough memory for overlap in\n"); exit(1); }
    p_spin_orbit = spin_orbit;
    for (i=0;i<dim;i++) { *p_spin_orbit = k_zero; p_spin_orbit++;}

    overlap = (double *) malloc(dim * sizeof(double));
    if (overlap == NULL) { fprintf(stderr, "ERROR: not enough memory for overlap in\n"); exit(1); }
    p_overlap = overlap;
    for (i=0;i<dim;i++) { *p_overlap = k_zero; p_overlap++;}

    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
    sheli = shells->type1_sh[index_i];
    imax = shells->imax_sh[index_i];

      shelposj = atoms->shelposn_sh[ip];
      gausposj = atoms->gausposn_sh[ip];
      bfposj = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[ip]; index_j++) {
      shelj = shells->type1_sh[index_j];
      jmax = shells->imax_sh[index_j];

      if (shelj != sheli) {
        bfposj += shelj;
        gausposj += shells->ng_sh[index_j];
        continue;
        }

        mm = imax + jmax;

        for (i4 = gausposi; i4 < gausposi + shells->ng_sh[index_i]; i4++) {
          for (j4 = gausposj; j4 < gausposj + shells->ng_sh[index_j]; j4++) {

            pab = gaussians->expo_sh[i4] + gaussians->expo_sh[j4];
            ab = gaussians->expo_sh[i4] * gaussians->expo_sh[j4];
            p32 = pab * sqrt(pab);
            SAB = pi32 / p32;
            sabfac = gaussians->c_sh[i4] * gaussians->c_sh[j4] * SAB;
            fn[0] = two * sqrt(pab / pi);
            for (i = 0; i < mm; i++) {
            fn[i+1] = fn[i] * -two * pab * (two * double(i) + k_one) / (two * double (i) + three);
            }

            for (i = 0; i <= imax + 2; i++) {
              for (j = 0; j <= jmax + 2; j++) {
                for (t = 0; t <= imax + jmax + 2; t++) {
                  E1x[t][i][j] = e_atom(i, j, t, pab);
                  E1y[t][i][j] = e_atom(i, j, t, pab);
                  E1z[t][i][j] = e_atom(i, j, t, pab);
                }
              }
            } // end loop to set up E factors

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

                 p_elecnuc = elecnuc + (bfposi + i) * nd1 + bfposj + j;
                  for (t = 0; t <= tmax; t++) {
                    for (u = 0; u <= umax; u++) {
                      for (v = 0; v <= vmax; v++) { 
                        *p_elecnuc -= E1x[t][n1][n4] * E1y[u][n2][n5] * E1z[v][n3][n6] * sabfac *
                        ftuvn_atom(t,u,v,0,fn) * atoms->atomic_number[ip];
                      }
                    }
                  } // end t u v loop

                 p_spin_orbit = spin_orbit + (bfposi + i) * nd1 + bfposj + j;
                  for (t = 0; t <= tmax - 1; t++) {
                    for (u = 0; u <= umax - 1; u++) {
                      for (v = 0; v <= vmax; v++) { 
                        if (n1 > 0 && n5 > 0)
                        *p_spin_orbit += double (n1 * n5) * E1x[t][n1 - 1][n4] * E1y[u][n2][n5 - 1] * E1z[v][n3][n6] * sabfac * \
                        ftuvn_atom(t,u,v,0,fn) * atoms->atomic_number[ip];
                        if (n2 > 0 && n4 > 0)
                        *p_spin_orbit -= E1x[t][n1][n4 - 1] * E1y[u][n2 - 1][n5] * E1z[v][n3][n6] * sabfac * \
                        ftuvn_atom(t,u,v,0,fn) * atoms->atomic_number[ip];
                      }
                    }
                  } // end t u v loop

                  for (t = 0; t <= tmax + 1; t++) {
                    for (u = 0; u <= umax - 1; u++) {
                      for (v = 0; v <= vmax; v++) { 
                        if (n5 > 0)
                       *p_spin_orbit -= two * double (n5) * gaussians->expo_sh[i4] * E1x[t][n1 + 1][n4] * E1y[u][n2][n5 - 1] * E1z[v][n3][n6] * \
                        sabfac * ftuvn_atom(t,u,v,0,fn) * atoms->atomic_number[ip];
                        if (n2 > 0)
                       *p_spin_orbit += two * double (n2) * gaussians->expo_sh[j4] * E1x[t][n1][n4 + 1] * E1y[u][n2 - 1][n5] * E1z[v][n3][n6] * \
                        sabfac * ftuvn_atom(t,u,v,0,fn) * atoms->atomic_number[ip];
                      }
                    }
                  } // end t u v loop

                  for (t = 0; t <= tmax - 1; t++) {
                    for (u = 0; u <= umax + 1; u++) {
                      for (v = 0; v <= vmax; v++) { 
                        if (n1 > 0)
                       *p_spin_orbit -= two * double (n1) * gaussians->expo_sh[j4] * E1x[t][n1 - 1][n4] * E1y[u][n2][n5 + 1] * E1z[v][n3][n6] * \
                        sabfac * ftuvn_atom(t,u,v,0,fn) * atoms->atomic_number[ip];
                        if (n4 > 0)
                       *p_spin_orbit += two * double (n4) * gaussians->expo_sh[i4] * E1x[t][n1][n4 - 1] * E1y[u][n2 + 1][n5] * E1z[v][n3][n6] * \
                        sabfac * ftuvn_atom(t,u,v,0,fn) * atoms->atomic_number[ip];
                      }
                    }
                  } // end t u v loop

                  for (t = 0; t <= tmax + 1; t++) {
                    for (u = 0; u <= umax + 1; u++) {
                      for (v = 0; v <= vmax; v++) { 
                      *p_spin_orbit += four * gaussians->expo_sh[i4] * gaussians->expo_sh[j4] * E1x[t][n1 + 1][n4] * E1y[u][n2][n5 + 1] * \
				       E1z[v][n3][n6] * \
                        sabfac * ftuvn_atom(t,u,v,0,fn) * atoms->atomic_number[ip];
                      *p_spin_orbit -= four * gaussians->expo_sh[i4] * gaussians->expo_sh[j4] * E1x[t][n1][n4 + 1] * E1y[u][n2 + 1][n5] * \
				       E1z[v][n3][n6] * \
                        sabfac * ftuvn_atom(t,u,v,0,fn) * atoms->atomic_number[ip];
                      }
                    }
                  } // end t u v loop

                  p_kinetic = kinetic + (bfposi + i) * nd1 + bfposj + j;
                 *p_kinetic += (- two * gaussians->expo_sh[j4] * double(2 * (n4 + n5 + n6) + 3)
                  * E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6] + four * gaussians->expo_sh[j4] * gaussians->expo_sh[j4]
                  *(E1x[0][n1][n4 + 2] * E1y[0][n2][n5] * E1z[0][n3][n6] + E1x[0][n1][n4] * E1y[0][n2][n5 + 2]
                  * E1z[0][n3][n6] + E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 + 2])) * sabfac;
                 if (n4 > 2)
                 *p_kinetic += double(n4 * (n4 - 1)) * E1x[0][n1][n4 - 2] * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac;
                 if (n5 > 2)
                 *p_kinetic += double(n5 * (n5 - 1)) * E1x[0][n1][n4] * E1y[0][n2][n5 - 2] * E1z[0][n3][n6] * sabfac;
                 if (n6 > 2)
                 *p_kinetic += double(n6 * (n6 - 1)) * E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6 - 2] * sabfac;

                 p_overlap = overlap + (bfposi + i) * nd1 + bfposj + j;
                *p_overlap += E1x[0][n1][n4] * E1y[0][n2][n5] * E1z[0][n3][n6] * sabfac;

               }
              }

          } // End j4 loop

        } // End i4 loop

        bfposj += shelj;
          gausposj += shells->ng_sh[index_j];

      } // close loop over index_j

      bfposi += sheli;
          gausposi += shells->ng_sh[index_i];

    } // close loop over index_i

    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli  = shells->type1[index_i];
    sheli1 = shells->type[index_i];
    nsheli = *(shells->num_ij + shells->ord[index_i]);
      shelposj = atoms->shelposn[ip];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[ip]; index_j++) {
      shelj  = shells->type1[index_j];      // 5 7
      shelj1 = shells->type[index_j];       // 6 10
      nshelj = *(shells->num_ij + shells->ord[index_j]);
        p_i1 = shells->ind_i + shells->opp[index_i];
        p_i  = shells->ind_j + shells->opp[index_i];
        p_rot1 = shells->rot + shells->opp[index_i];
        for (i = 0; i < nsheli; i++) {
          p_j1 = shells->ind_i + shells->opp[index_j];
          p_j  = shells->ind_j + shells->opp[index_j];
          p_rot2 = shells->rot + shells->opp[index_j];
          for (j = 0; j < nshelj; j++) {
            p_overlap = overlap + (bfposi1 + *p_i1) * nd1 + bfposj1 + *p_j1;
            p_kinetic = kinetic + (bfposi1 + *p_i1) * nd1 + bfposj1 + *p_j1;
            p_elecnuc = elecnuc + (bfposi1 + *p_i1) * nd1 + bfposj1 + *p_j1;
            p_spin_orbit = spin_orbit + (bfposi1 + *p_i1) * nd1 + bfposj1 + *p_j1;
            F_1e->a[bfposi + *p_i][bfposj + *p_j] += *p_rot1 * *p_rot2 * (-*p_kinetic/two + *p_elecnuc);
            Overlap->a[bfposi + *p_i][bfposj + *p_j] += *p_rot1 * *p_rot2 * (*p_overlap);
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj += shelj;
        bfposj1 += shelj1;
      } // close loop over index_j
      bfposi += sheli;
      bfposi1 += sheli1;
    } // close loop over index_i

    
  if (atoms->magnetic[*atm] == 2) {                        // double up F_1e for a spin polarised system
    for (i = 0; i < atoms->bfnnumb_sh[*atm]; i++) {
      for (j = 0; j < atoms->bfnnumb_sh[*atm]; j++) {
        F_1e->a[i + atoms->bfnnumb_sh[*atm]][j] = F_1e->a[i][j];
       }
      }
    if (job->spin_orb == 1) {                        // double up F_1e for a spin polarised system
    shelposi = atoms->shelposn[ip];
    bfposi = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[ip]; index_i++) {
    sheli  = shells->type1[index_i];
    sheli1 = shells->type[index_i];
    nsheli = *(shells->num_ij + shells->ord[index_i]);
      shelposj = atoms->shelposn[ip];
      bfposj = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel[ip]; index_j++) {
      shelj  = shells->type1[index_j];      // 5 7
      shelj1 = shells->type[index_j];       // 6 10
      nshelj = *(shells->num_ij + shells->ord[index_j]);
        p_i1 = shells->ind_i + shells->opp[index_i];
        p_i  = shells->ind_j + shells->opp[index_i];
        p_rot1 = shells->rot + shells->opp[index_i];
        for (i = 0; i < nsheli; i++) {
          p_j1 = shells->ind_i + shells->opp[index_j];
          p_j  = shells->ind_j + shells->opp[index_j];
          p_rot2 = shells->rot + shells->opp[index_j];
          for (j = 0; j < nshelj; j++) {
            p_spin_orbit = spin_orbit + (bfposi1 + *p_i1) * nd1 + bfposj1 + *p_j1;
            F_1e->a[bfposi + *p_i][bfposj + *p_j]                           += *p_rot1 * *p_rot2 * *p_spin_orbit/two;
            F_1e->a[bfposi + *p_i + atoms->bfnnumb_sh[*atm]][bfposj + *p_j] -= *p_rot1 * *p_rot2 * *p_spin_orbit/two;
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }
        bfposj += shelj;
        bfposj1 += shelj1;
      } // close loop over index_j
      bfposi += sheli;
      bfposi1 += sheli1;
    } // close loop over index_i

      } // close if spin_orb
     } // close if atoms->magnetic

  free(overlap);
  free(kinetic);
  free(elecnuc);
  free(spin_orbit);

}

void integrals_2e_atom(INTEGRAL_LIST *integral_list, int ip, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int index, index_i, index_j, index_k, index_l, i4, j4, k4, l4;
  int bfposi, bfposj, bfposk, bfposl;
  int bfposi1, bfposj1, bfposk1, bfposl1;
  int gausposi, gausposj, gausposk, gausposl;
  int sheli_lim, shelj_lim, shelk_lim, shell_lim;
  int i, j, k, l, m, n, p, q, mm, nd1;
  int n1, n2, n3, n4, n5, n6;
  int n7, n8, n9, n10, n11, n12;
  int t, u, v, tp, up, vp, tpupvp2;
  int tmax, umax, vmax, tpmax, upmax, vpmax;
  int imax, jmax, kmax, lmax;
  int shelposi, shelposj, shelposk, shelposl;
  int sheli, shelj, shelk, shell;
  int sheli1, shelj1, shelk1, shell1;
  int nsheli, nshelj, nshelk, nshell;
  int count;
  int dime1, dime2, dime3, dime4;
  int lim_ij, lim_kl;
  int slim_ij, slim_jk, slim_kl;
  int i1, j1, k1, l1;
  int *p_i, *p_j, *p_k, *p_l, *p_i1, *p_j1, *p_k1, *p_l1;
  double pab, pcd, p32, p_inv_ex;
  double SAB, SCD;
  double tpupvpsign;
  double en[55], em[55];
  double *F_temp, *p_F_temp, *F_sh, *p_F_sh;
  double fac, fac1, fac2;
  int tmp;
  double c1fac[9][9][9];
  double FTUV[17][17][17];
  double *p_rot1, *p_rot2, *p_rot3, *p_rot4;

  F_temp = (double *) malloc(50625 * sizeof(double));
  if (F_temp == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double F_temp! \n");
    exit(1);
  }

  F_sh = (double *) malloc(9 * 9 * 9 * 9 * sizeof(double));
  if (F_sh == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double F_sh! \n");
    exit(1);
  }

  count = 0;
  nd1 = atoms->bfnnumb_sh[ip];

    bfposi   = 0;
    bfposi1  = 0;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];

    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {

    sheli  = shells->type1_sh[index_i];
    sheli1 = shells->type_sh[index_i];
    imax   = shells->imax_sh[index_i];

    sheli_lim = sheli;

      bfposj   = 0;
      bfposj1  = 0;
      shelposj = atoms->shelposn_sh[ip];
      gausposj = atoms->gausposn_sh[ip];

      for (index = shelposi; index < index_i; index++) {
      bfposj   += shells->type1_sh[atoms->shelposn_sh[ip] + index_i - index - 1];
      bfposj1  += shells->type_sh[atoms->shelposn_sh[ip] + index_i - index - 1];
      gausposj += shells->ng_sh[atoms->shelposn_sh[ip] + index_i - index - 1];
      }

      for (index_j = index_i; index_j < shelposj + atoms->nshel_sh[ip]; index_j++) {

      shelj  = shells->type1_sh[index_j];
      shelj1 = shells->type_sh[index_j];
      jmax   = shells->imax_sh[index_j];

      shelj_lim = shelj;

       double sab[shells->ng_sh[index_i]][shells->ng_sh[index_j]];
       double C1x[shells->ng_sh[index_i]][shells->ng_sh[index_j]][imax + jmax + 1][imax + 1][jmax + 1];
       double pab_inv[shells->ng_sh[index_i]][shells->ng_sh[index_j]];

        for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
          lim_ij = 0;
            if (index_i == index_j) lim_ij = i4;
              for (j4 = lim_ij; j4 < shells->ng_sh[index_j]; j4++) {
                pab = gaussians->expo_sh[gausposi + i4] + gaussians->expo_sh[gausposj + j4];
                p32 = pab * sqrt(pab);
                SAB = pi32 / p32;
                pab_inv[i4][j4] = k_one / pab;
                sab[i4][j4] = gaussians->c_sh[gausposi + i4] * gaussians->c_sh[gausposj + j4] * SAB;

              for (t = 0; t <= imax + jmax; t++) {
                for (m = 0; m <= imax; m++) {
                  for (n = 0; n <= jmax; n++) {
                      C1x[i4][j4][t][m][n] = e_atom(m, n, t, pab);
                    }
                  }
                } // end loop to set up C1 factors

              }
            }

        bfposk   = 0;
        bfposk1  = 0;
        shelposk = atoms->shelposn_sh[ip];
        gausposk = atoms->gausposn_sh[ip];

        for (index = shelposi; index < index_i; index++) {
        bfposk   += shells->type1_sh[atoms->shelposn_sh[ip] + index_i - index - 1];
        bfposk1  += shells->type_sh[atoms->shelposn_sh[ip] + index_i - index - 1];
        gausposk += shells->ng_sh[atoms->shelposn_sh[ip] + index_i - index - 1];
        }

        for (index_k = index_i; index_k < shelposk + atoms->nshel_sh[ip]; index_k++) {

        shelk  = shells->type1_sh[index_k];
        shelk1 = shells->type_sh[index_k];
        kmax   = shells->imax_sh[index_k];

        shelk_lim = shelk;

          bfposl   = 0;
          bfposl1  = 0;
          shelposl = atoms->shelposn_sh[ip];
          gausposl = atoms->gausposn_sh[ip];

        for (index = shelposk; index < index_k; index++) {
        bfposl   += shells->type1_sh[atoms->shelposn_sh[ip] + index_k - index - 1];
        bfposl1  += shells->type_sh[atoms->shelposn_sh[ip] + index_k - index - 1];
        gausposl += shells->ng_sh[atoms->shelposn_sh[ip] + index_k - index - 1];
        }

          for (index_l = index_k; index_l < shelposl + atoms->nshel_sh[ip]; index_l++) {

          shell  = shells->type1_sh[index_l];
          shell1 = shells->type_sh[index_l];
          lmax   = shells->imax_sh[index_l];

        if ( (index_k == index_i && index_l < index_j) || 
             (((imax + jmax) / 2 ) *  2 == imax + jmax) && (((kmax + lmax + 1) / 2 ) * 2 == kmax + lmax + 1) ||
             (((imax + jmax + 1) / 2 ) *  2 == imax + jmax + 1) && (((kmax + lmax) / 2 ) * 2 == kmax + lmax)
           ) {
              bfposl   += shells->type1_sh[index_l];
              bfposl1  += shells->type_sh[index_l];
              gausposl += shells->ng_sh[index_l];
              continue;
             }

       shell_lim = shell;

       double scd[shells->ng_sh[index_k]][shells->ng_sh[index_l]];
       double C2x[shells->ng_sh[index_k]][shells->ng_sh[index_l]][kmax+lmax+1][kmax+1][lmax+1];
       double pcd_inv[shells->ng_sh[index_k]][shells->ng_sh[index_l]];

        for (k4 = 0; k4 < shells->ng_sh[index_k]; k4++) {
          lim_kl = 0;
            if (index_k == index_l) lim_kl = k4;
              for (l4 = lim_kl; l4 < shells->ng_sh[index_l]; l4++) {
                pcd = gaussians->expo_sh[gausposk + k4] + gaussians->expo_sh[gausposl + l4];
                p32 = pcd * sqrt(pcd);
                SCD = pi32 / p32;
                pcd_inv[k4][l4] = k_one / pcd;
                scd[k4][l4] = gaussians->c_sh[gausposk + k4] * gaussians->c_sh[gausposl + l4] * SCD;

                   for (t = 0; t <= kmax + lmax; t++) {
                    for (m = 0; m <= kmax; m++) {
                      for (n = 0; n <= lmax; n++) {
                      C2x[k4][l4][t][m][n] = e_atom(m, n, t, pcd);
                        }
                      }
                    } // end loop to set up C2 factors

                }
              }

            p_F_temp = F_temp;
            for (i = 0; i < sheli * shelj *shelk * shell; i++) {
             *p_F_temp = k_zero;
              p_F_temp++;
            }

                            dime4 = shell_lim ; 
                            dime3 = dime4*shelk_lim ; 
                            dime2 = dime3*shelj_lim ; 
                            dime1 = dime2*sheli_lim ;

            p_F_sh = F_sh;
            for (i = 0; i < sheli1 * shelj1 *shelk1 * shell1; i++) {
             *p_F_sh = k_zero;
              p_F_sh++;
            }

            mm = imax + jmax + kmax + lmax;
            en[0] = two / sqrt(pi);
            for (m = 0; m <= mm; m++) {
            en[m+1] = en[m] * -two * (two * double(m) + k_one) / (two * double (m) + three);
            }

                    for (i = 0; i < sheli_lim; i++) {
                      slim_ij = 0;
                       if (index_i == index_j) {
                         slim_ij = i;
                        }

                      for (j = slim_ij; j < shelj_lim; j++) {
                        tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
                        umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
                        vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
                        n1 = shells->tuv[imax][i][0];
                        n2 = shells->tuv[imax][i][1];
                        n3 = shells->tuv[imax][i][2];
                        n4 = shells->tuv[jmax][j][0];
                        n5 = shells->tuv[jmax][j][1];
                        n6 = shells->tuv[jmax][j][2];
                        slim_jk = 0;
                         for (k = 0; k < shelk_lim; k++) {
                          slim_kl = 0;
                            if (index_k == index_l) {
                              slim_kl = k;
                             }

                          for (l = slim_kl; l < shell_lim; l++) {
                            if (index_i == index_k && index_j == index_l && i * shelj_lim + j > k * shell_lim + l) 
                              continue;
                        tpmax = shells->tuv[kmax][k][0] + shells->tuv[lmax][l][0];
                        upmax = shells->tuv[kmax][k][1] + shells->tuv[lmax][l][1];
                        vpmax = shells->tuv[kmax][k][2] + shells->tuv[lmax][l][2];
                        n7 = shells->tuv[kmax][k][0];
                        n8 = shells->tuv[kmax][k][1];
                        n9 = shells->tuv[kmax][k][2];
                        n10 = shells->tuv[lmax][l][0];
                        n11 = shells->tuv[lmax][l][1];
                        n12 = shells->tuv[lmax][l][2];
        for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
          lim_ij = 0;
            if (index_i == index_j) lim_ij = i4;
              for (j4 = lim_ij; j4 < shells->ng_sh[index_j]; j4++) {

                for (t = 0; t <= tmax; t++) {
                for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                c1fac[t][u][v] = k_zero ; }}}

                  for (k4 = 0; k4 < shells->ng_sh[index_k]; k4++) {
                    lim_kl = 0;
                      if (index_k == index_l) lim_kl = k4;
                        for (l4 = lim_kl; l4 < shells->ng_sh[index_l]; l4++) {

                          if (index_i == index_k && index_j == index_l && \
                          i4 * shells->ng_sh[index_j] + j4 > k4 * shells->ng_sh[index_l] + l4) 
                            continue;

		            fac1 = sab[i4][j4] * scd[k4][l4];

		              if (index_i == index_j && i4 != j4) fac1 *= two;
		              if (index_k == index_l && k4 != l4) fac1 *= two; 
		              if (index_i == index_k && index_j == index_l && (i4 != k4 || j4 != l4)) fac1 *= two; 

                                p_inv_ex = pab_inv[i4][j4] + pcd_inv[k4][l4];

                                     f000m(em, k_zero, k_one / p_inv_ex, mm);
                                     //g000m(em, k_zero, k_one / p_inv_ex, mm)

                    for (m = 0; m <= mm; m++) {
                      for (n = 0; n <= mm; n++) {
                        for (p = 0; p <= mm; p++) {
                          if (m + n + p > mm) {
                            continue;
                          }
                          FTUV[m][n][p] = ftuvn_atom(m,n,p,0,em);
                        }
                      }
                    }


                            for (t = 0; t <= tmax; t++) {
                              if (fabs(C1x[i4][j4][t][n1][n4]) < 0.000000001)
                                continue;
                              for (u = 0; u <= umax; u++) {
                              if (fabs(C1x[i4][j4][u][n2][n5]) < 0.000000001)
                                  continue;
                                for (v = 0; v <= vmax; v++) {
                              if (fabs(C1x[i4][j4][v][n3][n6]) < 0.000000001)
                                    continue;
                                  for (tp = 0; tp <= tpmax; tp++) {
                              if (fabs(C2x[k4][l4][tp][n7][n10]) < 0.000000001)
                                      continue;
                                    for (up = 0; up <= upmax; up++) {
                              if (fabs(C2x[k4][l4][up][n8][n11]) < 0.000000001)
                                        continue;
                                      for (vp = 0; vp <= vpmax; vp++) {
                              if (fabs(C2x[k4][l4][vp][n9][n12]) < 0.000000001)
                                          continue;
                                          tpupvpsign = -k_one;
                                          tpupvp2 = tp + up + vp + 2;
                                          if ((tpupvp2 / 2) * 2 == tpupvp2)
                                            tpupvpsign = k_one;
                    c1fac[t][u][v] += fac1 * C2x[k4][l4][tp][n7][n10] * C2x[k4][l4][up][n8][n11] * C2x[k4][l4][vp][n9][n12] *
                       FTUV[t + tp][u + up][v + vp] * tpupvpsign;
                                }
                              }
                            } // end tp up vp loop
                                      }
                                    }
                                  } // end t u v loop

                  } // End l4 loop
                } // End k4 loop
   
                
                            for (t = 0; t <= tmax; t++) {
                              if (fabs(C1x[i4][j4][t][n1][n4]) < 0.000000001)
                                continue;
                              for (u = 0; u <= umax; u++) {
                              if (fabs(C1x[i4][j4][u][n2][n5]) < 0.000000001)
                                  continue;
                                for (v = 0; v <= vmax; v++) {
                              if (fabs(C1x[i4][j4][v][n3][n6]) < 0.000000001)
                                    continue;

                                p_F_temp = F_temp + dime2 * i + dime3 * j + dime4 * k + l ;
                               *p_F_temp += c1fac[t][u][v] * C1x[i4][j4][t][n1][n4] * C1x[i4][j4][u][n2][n5] * 
                                                       C1x[i4][j4][v][n3][n6];
                                     }}}
              } // End j4 loop
            } // End i4 loop

                          }
                        }
                      }
                    } // end ijkl loop

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
              i1 = *p_i1;
              j1 = *p_j1;
              k1 = *p_k1;
              l1 = *p_l1;
              fac2 = k_one;
              if (index_i == index_j && *p_i > *p_j) fac2 = k_zero;
              if (index_k == index_l && *p_k > *p_l) fac2 = k_zero;
              if (index_i == index_k && index_j == index_l && *p_i * shelj1 + *p_j > *p_k * shell1 + *p_l) fac2 = k_zero;
              if (index_i == index_j && *p_i1 > *p_j1) { j1 = *p_i1; i1 = *p_j1; }
              if (index_k == index_l && *p_k1 > *p_l1) { l1 = *p_k1; k1 = *p_l1; }
              if (index_i == index_k && index_j == index_l && i1 * shelj + j1 > k1 * shell + l1) { 
              tmp = i1; i1 = k1; k1 = tmp; tmp = j1; j1 = l1; l1 = tmp;}
                p_F_temp = F_temp + i1 * shelj * shelk * shell + j1 * shelk * shell + k1 * shell + l1;
                p_F_sh   = F_sh   + *p_i * shelj1 * shelk1 * shell1 + *p_j * shelk1 * shell1 + *p_k  * shell1 + *p_l;
               *p_F_sh   += fac2 * *p_F_temp * *p_rot1 * *p_rot2 * *p_rot3 * *p_rot4;
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

        p_F_sh = F_sh;
          for (i = 0; i < sheli1; i++) {
            for (j = 0; j < shelj1; j++) {
              for (k = 0; k < shelk1; k++) {
                for (l = 0; l < shell1; l++) {
                 fac = k_one;
                   if (index_i == index_j && i == j) fac /= two;
                   if (index_k == index_l && k == l) fac /= two;
                   if ((index_i == index_k && index_j == index_l) && (i == k && j == l)) fac /= two;
                   if (fabs(*p_F_sh) > 1e-09) { 
                     integral_list->value[count] = *p_F_sh * fac;
                     integral_list->i[count] = bfposi1 + i;
                     integral_list->j[count] = bfposj1 + j;
                     integral_list->k[count] = bfposk1 + k;
                     integral_list->l[count] = bfposl1 + l;
                     count++;
                    }
                 p_F_sh++;
                }
               }
              }
             }

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

   free(F_temp);
   free(F_sh);

}

void fock_2e_matrix(DoubleMatrix *F_2e, INTEGRAL_LIST *integrals, DoubleMatrix *P_reduced, int atm, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  int m, offset; 

 if (atoms->magnetic[atm] == 1) {
  for (m = 0; m < integrals->num; m++) {
    F_2e->a[integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->l[m]];
    F_2e->a[integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->k[m]];
    F_2e->a[integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->l[m]];
    F_2e->a[integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->k[m]];
    F_2e->a[integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->j[m]];
    F_2e->a[integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->i[m]];
    F_2e->a[integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->j[m]];
    F_2e->a[integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->i[m]];

    F_2e->a[integrals->i[m]][integrals->k[m]] -= integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->l[m]] / two;
    F_2e->a[integrals->i[m]][integrals->l[m]] -= integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->k[m]] / two;
    F_2e->a[integrals->j[m]][integrals->k[m]] -= integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->l[m]] / two;
    F_2e->a[integrals->j[m]][integrals->l[m]] -= integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->k[m]] / two;
    F_2e->a[integrals->k[m]][integrals->i[m]] -= integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->j[m]] / two;
    F_2e->a[integrals->k[m]][integrals->j[m]] -= integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->i[m]] / two;
    F_2e->a[integrals->l[m]][integrals->i[m]] -= integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->j[m]] / two;
    F_2e->a[integrals->l[m]][integrals->j[m]] -= integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->i[m]] / two;

   }
  }

 if (atoms->magnetic[atm] == 2) {
  offset = atoms->bfnnumb_sh[atm];
  for (m = 0; m < integrals->num; m++) {

    F_2e->a[integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->l[m]];
    F_2e->a[integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->k[m]];
    F_2e->a[integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->l[m]];
    F_2e->a[integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->k[m]];
    F_2e->a[integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->j[m]];
    F_2e->a[integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->i[m]];
    F_2e->a[integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->j[m]];
    F_2e->a[integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->i[m]];

    F_2e->a[integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[offset + integrals->k[m]][integrals->l[m]];
    F_2e->a[integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[offset + integrals->l[m]][integrals->k[m]];
    F_2e->a[integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[offset + integrals->k[m]][integrals->l[m]];
    F_2e->a[integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[offset + integrals->l[m]][integrals->k[m]];
    F_2e->a[integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[offset + integrals->i[m]][integrals->j[m]];
    F_2e->a[integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[offset + integrals->j[m]][integrals->i[m]];
    F_2e->a[integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[offset + integrals->i[m]][integrals->j[m]];
    F_2e->a[integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[offset + integrals->j[m]][integrals->i[m]];

    F_2e->a[offset + integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->l[m]];
    F_2e->a[offset + integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->k[m]];
    F_2e->a[offset + integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->l[m]];
    F_2e->a[offset + integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->k[m]];
    F_2e->a[offset + integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->j[m]];
    F_2e->a[offset + integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->i[m]];
    F_2e->a[offset + integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->j[m]];
    F_2e->a[offset + integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->i[m]];

    F_2e->a[offset + integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[offset + integrals->k[m]][integrals->l[m]];
    F_2e->a[offset + integrals->i[m]][integrals->j[m]] += integrals->value[m] * P_reduced->a[offset + integrals->l[m]][integrals->k[m]];
    F_2e->a[offset + integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[offset + integrals->k[m]][integrals->l[m]];
    F_2e->a[offset + integrals->j[m]][integrals->i[m]] += integrals->value[m] * P_reduced->a[offset + integrals->l[m]][integrals->k[m]];
    F_2e->a[offset + integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[offset + integrals->i[m]][integrals->j[m]];
    F_2e->a[offset + integrals->k[m]][integrals->l[m]] += integrals->value[m] * P_reduced->a[offset + integrals->j[m]][integrals->i[m]];
    F_2e->a[offset + integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[offset + integrals->i[m]][integrals->j[m]];
    F_2e->a[offset + integrals->l[m]][integrals->k[m]] += integrals->value[m] * P_reduced->a[offset + integrals->j[m]][integrals->i[m]];

    F_2e->a[integrals->i[m]][integrals->k[m]] -= integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->l[m]];
    F_2e->a[integrals->i[m]][integrals->l[m]] -= integrals->value[m] * P_reduced->a[integrals->j[m]][integrals->k[m]];
    F_2e->a[integrals->j[m]][integrals->k[m]] -= integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->l[m]];
    F_2e->a[integrals->j[m]][integrals->l[m]] -= integrals->value[m] * P_reduced->a[integrals->i[m]][integrals->k[m]];
    F_2e->a[integrals->k[m]][integrals->i[m]] -= integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->j[m]];
    F_2e->a[integrals->k[m]][integrals->j[m]] -= integrals->value[m] * P_reduced->a[integrals->l[m]][integrals->i[m]];
    F_2e->a[integrals->l[m]][integrals->i[m]] -= integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->j[m]];
    F_2e->a[integrals->l[m]][integrals->j[m]] -= integrals->value[m] * P_reduced->a[integrals->k[m]][integrals->i[m]];

    F_2e->a[offset + integrals->i[m]][integrals->k[m]] -= integrals->value[m] * P_reduced->a[offset + integrals->j[m]][integrals->l[m]];
    F_2e->a[offset + integrals->i[m]][integrals->l[m]] -= integrals->value[m] * P_reduced->a[offset + integrals->j[m]][integrals->k[m]];
    F_2e->a[offset + integrals->j[m]][integrals->k[m]] -= integrals->value[m] * P_reduced->a[offset + integrals->i[m]][integrals->l[m]];
    F_2e->a[offset + integrals->j[m]][integrals->l[m]] -= integrals->value[m] * P_reduced->a[offset + integrals->i[m]][integrals->k[m]];
    F_2e->a[offset + integrals->k[m]][integrals->i[m]] -= integrals->value[m] * P_reduced->a[offset + integrals->l[m]][integrals->j[m]];
    F_2e->a[offset + integrals->k[m]][integrals->j[m]] -= integrals->value[m] * P_reduced->a[offset + integrals->l[m]][integrals->i[m]];
    F_2e->a[offset + integrals->l[m]][integrals->i[m]] -= integrals->value[m] * P_reduced->a[offset + integrals->k[m]][integrals->j[m]];
    F_2e->a[offset + integrals->l[m]][integrals->j[m]] -= integrals->value[m] * P_reduced->a[offset + integrals->k[m]][integrals->i[m]];

   }
  }

}

double e_atom(int i, int ip, int t, double p)

{

  double recurrence;
  int flag;
  flag = (t == 0 && i == 0 && ip == 0) ? 1 : 0;
  if (flag == 1) {
    return k_one;
  }

  else if

  (i > ip)

  {
    recurrence = (t >= 0 && t <= i + ip) ? (e_atom(i - 1, ip, t - 1, p) / two / p + 
        (double)(t + 1) * e_atom(i - 1, ip, t + 1, p)) : k_zero;
  }

  else {
    recurrence = (t >= 0 && t <= i + ip) ? (e_atom(i, ip - 1, t - 1, p) / two / p + 
        + (double)(t + 1) * e_atom(i, ip - 1, t + 1, p)) : k_zero;
  }
  return recurrence;

}

double ftuvn_atom(int t, int u, int v, int n, double *f)

{

  double result;
  if (t < 0 || u < 0 || v < 0)
    return 0.0;
  if (t >= u && t >= v) {
    //fprintf(file_out,"t %d %d %d %d \n",t,u,v,n) ;
    if (t == 0 && u == 0 && v == 0)
      return f[n];
    result = (t > 0) ? (t - 1) * ftuvn_atom(t - 2, u, v, n + 1, f) : k_zero;
  } else if (u > t && u >= v) {
    result = (u > 0) ? (u - 1) * ftuvn_atom(t, u - 2, v, n + 1, f) : k_zero;
  } else {
    result = (v > 0) ? (v - 1) * ftuvn_atom(t, u, v - 2, n + 1, f) : k_zero;
  }
    //fprintf(file_out,"ftuvn_atom tuv n %d %d %d %d %lf \n",t,u,v,n,result) ;
  return result;

}

void g000m(int i, double *f, double x, int n)

{

  int i1, i2, j, type;
  double g[16], T;
  double i0;
  double fac = k_one;
  f[0] = -fac;
  for (i2 = 0; i2 <= n; i2++) {
    f[i2] = -two / rtpi / fac;
    fac *= (i2 + 0.5);
   }

  static double G[4338] = { 1.000000000000000, 0.333333333333333, 0.200000000000000, 0.142857142857143,
      0.111111111111111, 0.090909090909091, 0.076923076923077, 0.066666666666667, 0.058823529411765, 0.052631578947368,
      0.047619047619048, 0.043478260869565, 0.040000000000000, 0.037037037037037, 0.034482758620690, 0.032258064516129,
      0.030303030303030, 0.028571428571429, 0.983580385842959, 0.323509613422450, 0.192994157666346, 0.137413638310181,
      0.106660436705514, 0.087145058489144, 0.073662188798712, 0.063790298825366, 0.056250578797778, 0.050304150615128,
      0.045494371867164, 0.041523847097296, 0.038190587371005, 0.035352597744162, 0.032907145916707, 0.030778070837964,
      0.028907714761751, 0.027251626370837, 0.967643312635592, 0.314029472998161, 0.186255004792621, 0.132188029635739,
      0.102393947071065, 0.083540528018141, 0.070541950817984, 0.061039712989142, 0.053791384005819, 0.048080550314778,
      0.043465189724101, 0.039657830850816, 0.036463457664022, 0.033745117822999, 0.031403815925110, 0.029366218961129,
      0.027576848795228, 0.025992961032804, 0.952171930273449, 0.304879846161303, 0.179771873529508, 0.127171304074947,
      0.098303840331894, 0.080088621873286, 0.067556213936973, 0.058409349185289, 0.051440871180917, 0.045956112168412,
      0.041527182582549, 0.037876192694927, 0.034814851860861, 0.032211066988252, 0.029969440859114, 0.028019361630848,
      0.026307447104062, 0.024792593363273, 0.937150028797979, 0.296048189299992, 0.173534537054986, 0.122354830492368,
      0.094382650921487, 0.076782763038513, 0.064699100864166, 0.055893895390439, 0.049194194446507, 0.043926381281596,
      0.039676228180875, 0.036175096801002, 0.033241183362661, 0.030747077471355, 0.028600846621503, 0.026734497364043,
      0.025096663018395, 0.023647816322620, 0.922562012825585, 0.287522459508360, 0.167533190907351, 0.117730342930696,
      0.090623234886940, 0.073616661822110, 0.061964993943621, 0.053488276391324, 0.047046725596908, 0.041987104152047,
      0.037908391634972, 0.034550882525999, 0.031739030053163, 0.029349936515360, 0.027295005686608, 0.025508763680457,
      0.023941782045520, 0.022556048861491, 0.908392877032751, 0.279291093918388, 0.161758435122409, 0.113289924883880,
      0.087018755842400, 0.070584303166469, 0.059348523582399, 0.051187643149114, 0.044994044258325, 0.040134219516345,
      0.036219916881389, 0.033000056379079, 0.030305126728486, 0.028016579217385, 0.026049030312793, 0.024339430648799,
      0.022840215718441, 0.021514830044712, 0.894628182652567, 0.271342989905505, 0.156201257139716, 0.109025994256952,
      0.083562671542786, 0.067679934523372, 0.056844557197679, 0.048987362644441, 0.043031928496997, 0.038363849614615,
      0.034607218512822, 0.031519284357922, 0.028936357876419, 0.026744081702510, 0.024860166070068, 0.023223894733246,
      0.021789495731320, 0.020521813449768, 0.881254034939924, 0.263667486130357, 0.150853015444288, 0.104931288982249,
      0.080248721050131, 0.064898054269425, 0.054448188660050, 0.046883008181259, 0.041156345854049, 0.036672291853998,
      0.033066873987895, 0.030105384637691, 0.027629750789064, 0.025529654613710, 0.023725785668161, 0.022159672926270,
      0.020787268348402, 0.019574761827027, 0.868257061564465, 0.256254344380769, 0.145705423911703, 0.100998853263045,
      0.077070912466158, 0.062233400637391, 0.052154728210585, 0.044870350128707, 0.039363444787597, 0.035056010852631,
      0.031595616198023, 0.028755320596340, 0.026382468993385, 0.024370636903175, 0.022643383071061, 0.021144397154453,
      0.019831289073644, 0.018671542009427, 0.855624391892149, 0.249093732179515, 0.140750536825913, 0.097222024416930,
      0.074023511205878, 0.059680941140265, 0.049959692830285, 0.042945347081075, 0.037649546503491, 0.033511630846714,
      0.030190326374924, 0.027466194160764, 0.025191805984946, 0.023264489911013, 0.021610567884731, 0.020175808944565,
      0.018919417568867, 0.017810120059976, 0.843343637117978, 0.242176206124992, 0.135980734540445, 0.093594420292489,
      0.071101028788122, 0.057235862466009, 0.047858797041466, 0.041104137416879, 0.036011137156993, 0.032035928443995,
      0.028848027323112, 0.026235239458973, 0.024055179250812, 0.022208791718009, 0.020625060005242, 0.019251754337748,
      0.018049612808818, 0.016988556645913, 0.831402871214021, 0.235492695933329, 0.131388709754967, 0.090109927234008,
      0.068298212120022, 0.054893560821813, 0.045847944121597, 0.039343031238943, 0.034444860408437, 0.030625825707839,
      0.027565876962430, 0.025059816764170, 0.022970124568227, 0.021201231759716, 0.019684684515254, 0.018370179040275,
      0.017219928462084, 0.016205002628946, 0.819790612658427, 0.229034489151855, 0.126967454380422, 0.086762688570073,
      0.065610033253455, 0.052649632707753, 0.043923217710974, 0.037658502678186, 0.032947510316749, 0.029278383556699,
      0.026341162166356, 0.023937406717270, 0.021934290566308, 0.020239605689758, 0.018787366817268, 0.017529123799823,
      0.016428508487300, 0.015457694861437, 0.808495806912583, 0.222793216515124, 0.122710246967117, 0.083547093602981,
      0.063031679592470, 0.050499866100585, 0.042080873796450, 0.036047182544598, 0.031516024555401, 0.027990795464574,
      0.025171292882497, 0.022865604815018, 0.020945433538580, 0.019321810480782, 0.017931127992651, 0.016726719996769,
      0.015673582934602, 0.014744952178906, 0.797507809614983, 0.216760837915979, 0.118610640671282, 0.080457767076930,
      0.060558544531665, 0.048440232029311, 0.040317333054272, 0.034505851309679, 0.030147477936118, 0.026760381448656,
      0.024053796522294, 0.021842116151442, 0.020001412494763, 0.018445839752033, 0.017114080375922, 0.015961185440488,
      0.014953463942743, 0.014065171579673, 0.786816370461720, 0.210929628965312, 0.114662451736696, 0.077489559103912,
      0.058186218506350, 0.046466876524958, 0.038629173535821, 0.033031432405285, 0.028839076226285, 0.025584582331017,
      0.022986312607561, 0.020864750400973, 0.019100184440724, 0.017609779313041, 0.016334423334309, 0.015230820361083,
      0.014266541922722, 0.013416824582877, 0.776411617744798, 0.205292168115336, 0.110859748468989, 0.074637535527187,
      0.055910480436226, 0.044576112927827, 0.037013123680803, 0.031620985824540, 0.027588150246693, 0.024460954261802,
      0.021966587662061, 0.019931417032088, 0.018239799876060, 0.016811802913402, 0.015590439243014, 0.014534003587465,
      0.013611281919237, 0.012798453756530, 0.766284043520696, 0.199841324322276, 0.107196840681238, 0.071896968703107,
      0.053727289545085, 0.042764414536202, 0.035466055643125, 0.030271702011126, 0.026392150236829, 0.023387163491945,
      0.020992470336862, 0.019040120740836, 0.017418398499238, 0.016050168189089, 0.014880489647109, 0.013869188903083,
      0.012986220141648, 0.012208669407655, 0.756424489382800, 0.194570245225421, 0.103668269590400, 0.069263328682895,
      0.051632777539874, 0.041028407581246, 0.033984978915373, 0.028980896023868, 0.025248640475536, 0.022360981384005,
      0.020061906758730, 0.018188957094122, 0.016634205110690, 0.015323212796188, 0.014203011601361, 0.013234901571030,
      0.012389960656535, 0.011646146426928, 0.746824132812427, 0.189472345820492, 0.100268798145017, 0.066732274776822,
      0.049623241133157, 0.039364864513484, 0.032567034238442, 0.027746001964150, 0.024155294145404, 0.021380279650214,
      0.019172936091315, 0.017376108373082, 0.015885525704727, 0.014629350723368, 0.013556514179750, 0.012629735020648,
      0.011821172234322, 0.011109621280590, 0.737474474084268, 0.184541297606244, 0.096993401765513, 0.064299647484004,
      0.047695134893750, 0.037770697586949, 0.031209487783468, 0.026564567654249, 0.023109888429800, 0.020443025807358,
      0.018323686299358, 0.016599839607320, 0.015170743741529, 0.013967068774791, 0.012939575146764, 0.012052347688095,
      0.011278585342756, 0.010597889142759, 0.728367323703094, 0.179771018184097, 0.093837259479187, 0.061961460771753,
      0.045845064410997, 0.036242952727679, 0.029909725593813, 0.025434249555223, 0.022110299831936, 0.019547278838559,
      0.017512370106608, 0.015858494791219, 0.014488316590894, 0.013334923215579, 0.012350837782974, 0.011501460003709,
      0.010760989280409, 0.010109801161560, 0.719494790349553, 0.175155661291522, 0.090795745432831, 0.059713894689175,
      0.044069779758770, 0.034778803673860, 0.028665248275393, 0.024352807913504, 0.021154499705872, 0.018691185052509,
      0.016737281138529, 0.015150493273940, 0.013836772139813, 0.012731536572288, 0.011789007857708, 0.010975851519337,
      0.010267229443653, 0.009644261852824, 0.710849269313849, 0.170689607250686, 0.087864420766607, 0.057553288300347,
      0.042366169245927, 0.033375546375474, 0.027473665924173, 0.023318102125851, 0.020240549989818, 0.017872974131129,
      0.015996790241354, 0.014474326315099, 0.013214705556285, 0.012155594581222, 0.011252850742002, 0.010474358169107,
      0.009796204720884, 0.009200226615398, 0.702423431396967, 0.166367453814711, 0.085039025833577, 0.055476132923078,
      0.040731253440541, 0.032030593641872, 0.026332693280162, 0.022328086312765, 0.019366599132516, 0.017090955357030,
      0.015289341969355, 0.013828553798506, 0.012620776202180, 0.011605843277720, 0.010741188655302, 0.009995869657431,
      0.009346865008064, 0.008776699362369, 0.694210212260006, 0.162184007394613, 0.082315472749933, 0.053479065659867,
      0.039162179455790, 0.030741470026192, 0.025240145097729, 0.021380805090948, 0.018530878203924, 0.016343514012572,
      0.014613451232636, 0.013211801096670, 0.012053704688230, 0.011081086219899, 0.010252898039712, 0.009539326968316,
      0.008918208839917, 0.008372730262794, 0.686202802202985, 0.158134274650775, 0.079689838261642, 0.051558863208267,
      0.037656215485917, 0.029505806936059, 0.024193931722502, 0.020474389535789, 0.017731697181829, 0.015629107942668,
      0.013967700098076, 0.012622756079147, 0.011512270064622, 0.010580181840619, 0.009786907055856, 0.009103719990343,
      0.008509281131381, 0.007987413588766, 0.678394636355279, 0.154213454433455, 0.077158356913842, 0.049712435938430,
      0.036210745581215, 0.028321337960474, 0.023192054865572, 0.019607053325298, 0.016967441406379, 0.014946264273869,
      0.013350734736394, 0.012060166258092, 0.010995307140899, 0.010102040921734, 0.009342193194723, 0.008688085251913,
      0.008119171024173, 0.007619885662896, 0.670779385260629, 0.150416930057528, 0.074717414509926, 0.047936822226150,
      0.034823264651465, 0.027185894403238, 0.022232603566144, 0.018777089057267, 0.016236568194897, 0.014293576282572,
      0.012761262508642, 0.011522836064717, 0.010501703929208, 0.009645624184965, 0.008917781000089, 0.008291503761651,
      0.007747009833582, 0.007269322901516, 0.663350945840335, 0.146740261897302, 0.072363541847825, 0.046229183030232,
      0.033491373687731, 0.026097401013716, 0.021313750334149, 0.017982864731838, 0.015537603609712, 0.013669700405558,
      0.012198049185727, 0.011009624250612, 0.010030399205216, 0.009209939993986, 0.008512739896400, 0.007913098949057,
      0.007391969090783, 0.006934939949132, 0.656103432719009, 0.143179180287828, 0.070093408721529, 0.044586796703517,
      0.032212775192863, 0.025053871906138, 0.020433747464765, 0.017222820392002, 0.014869139372027, 0.013073353386359,
      0.011659916294865, 0.010519441408198, 0.009580380181232, 0.008794042162596, 0.008126182117209, 0.007552034700746,
      0.007053258676252, 0.006615987899866, 0.649031169897862, 0.139729578719752, 0.067903818176438, 0.043007054027354,
      0.030985268811507, 0.024053406659034, 0.019590923517100, 0.016495464914889, 0.014229829915213, 0.012503309551238,
      0.011145738587148, 0.010051247604830, 0.009150680286383, 0.008397027864040, 0.007757260729507, 0.007207513487829,
      0.006730125040014, 0.006311752601813, 0.642128682761166, 0.136387507315276, 0.065791701007599, 0.041487453459769,
      0.029806747150797, 0.023094186586794, 0.018783679949692, 0.015799372947044, 0.013618389571183, 0.011958398208895,
      0.010654441620680, 0.009604050125308, 0.008740377048890, 0.008018035636818, 0.007405167749497, 0.006878774580199,
      0.006421849504675, 0.006021553040458, 0.635390690402153, 0.133149166573358, 0.063754110490394, 0.040025596588011,
      0.028675191783335, 0.022174471175672, 0.018010487905780, 0.015133181977178, 0.013033589883805, 0.011437501168223,
      0.010184999453970, 0.009176901317836, 0.008348590075730, 0.007656243482505, 0.007069132345559, 0.006565092343673,
      0.006127746647395, 0.005744739797441, 0.628812098255143, 0.130010901372771, 0.061788217333676, 0.038619183776553,
      0.027588669424408, 0.021292594676922, 0.017269885141627, 0.014495589540202, 0.012474257043594, 0.010939550368757,
      0.009736432444556, 0.008768896538637, 0.007974479125203, 0.007310867051325, 0.006748419124377, 0.006265774616139,
      0.005847162757106, 0.005480693581160, 0.622387991021315, 0.126969195222147, 0.059891304845793, 0.037266010002049,
      0.026545328275766, 0.020446962850086, 0.016560473091488, 0.013885350546598, 0.011939269438161, 0.010463525618654,
      0.009307805148010, 0.008379172190729, 0.007617242268108, 0.006981157911424, 0.006442326496349, 0.005980161159035,
      0.005579474363475, 0.005228823825861, 0.616113625876012, 0.124020664746590, 0.058060764304363, 0.035963960867078,
      0.025543394528626, 0.019636049849733, 0.015880914063089, 0.013301274731496, 0.011427555313192, 0.010008452435306,
      0.008898224312755, 0.008006903852494, 0.007276114133444, 0.006666401897968, 0.006150185116623, 0.005707622180659,
      0.005324086834268, 0.004988567356006, 0.609984425946011, 0.121162054400888, 0.056294090521061, 0.034711008784913,
      0.024581169018883, 0.018858395249291, 0.015229928557782, 0.012742224218034, 0.010938090538914, 0.009573399983922,
      0.008506836966286, 0.007651304491941, 0.006950364234737, 0.006365917538368, 0.005871356398237, 0.005447556927960,
      0.005080433037927, 0.004759387112883, 0.603995974045666, 0.118390231399783, 0.054588877593034, 0.033505209327861,
      0.023657024027824, 0.018112601195873, 0.014606292709767, 0.012207111189863, 0.010469896477288, 0.009157479109586,
      0.008132828588623, 0.007311622762709, 0.006639295373279, 0.006079054550119, 0.005605231094025, 0.005199392343644,
      0.004847972068322, 0.004540770940540, 0.598144006661304, 0.115702180856173, 0.052942814832976, 0.032344697732067,
      0.022769400221965, 0.017397329690268, 0.014008835839083, 0.011694895667866, 0.010022037945344, 0.008759840458557,
      0.007775421368994, 0.006987141378066, 0.006342242114727, 0.005805192407889, 0.005351227944096, 0.004962581785540,
      0.004626188028780, 0.004332230428283, 0.592424408173695, 0.113095001118510, 0.051353682870177, 0.031227685551971,
      0.021916803725852, 0.016711299986552, 0.013436438113235, 0.011204583386404, 0.009593621270305, 0.008379672684726,
      0.007433872541951, 0.006677175559309, 0.006058569335684, 0.005543738976659, 0.005108792385852, 0.004736603805342,
      0.004414588872633, 0.004133299807096, 0.586833205308847, 0.110565899299016, 0.049819349915253, 0.030152457457925,
      0.021097803322022, 0.016053286106005, 0.012888028312635, 0.010735223764589, 0.009183792432348, 0.008016200737364,
      0.007107472799268, 0.006381071555154, 0.005787670837040, 0.005294129207862, 0.004877395323643, 0.004520960983971,
      0.004212705297645, 0.003943534897455, 0.581366561807785, 0.108112186984718, 0.048337768181548, 0.029117368170755,
      0.020311027772508, 0.015422114460250, 0.012362581695176, 0.010285907968323, 0.008791735291012, 0.007668684226444,
      0.006795544774171, 0.006098205228858, 0.005528968020983, 0.005055823895602, 0.004656531955290, 0.004315178820913,
      0.004020089691818, 0.003762512106161, 0.576020773306375, 0.105731276123646, 0.046906970456501, 0.028120839527311,
      0.019555163256556, 0.014816661578788, 0.011859117955530, 0.009855767058991, 0.008416669891485, 0.007336415862025,
      0.006497441594577, 0.005827980709952, 0.005281908628765, 0.004828308490180, 0.004445720652845, 0.004118804675042,
      0.003836315128175, 0.003589827469874, 0.570792262416601, 0.103420675078830, 0.045525066816584, 0.027161357671345,
      0.018828950919456, 0.014235851936276, 0.011376699274927, 0.009443970224930, 0.008057850847131, 0.007018719964304,
      0.006212545502202, 0.005569829107639, 0.005045965536406, 0.004611091966285, 0.004244501895075, 0.003931406754516,
      0.003660974406250, 0.003425095743197, 0.565677574001107, 0.101177984843109, 0.044190241479679, 0.026237470364259,
      0.018131184527611, 0.013678655875151, 0.010914428457361, 0.009049723091933, 0.007714565794824, 0.006714951041133,
      0.005940266534504, 0.005323207282997, 0.004820635605682, 0.004403705743313, 0.004052437249270, 0.003752573153485,
      0.003493679138097, 0.003267949529214, 0.560673370633112, 0.099000895407992, 0.042900749789027, 0.025347784410550,
      0.017460708225170, 0.013144087619358, 0.010471447148381, 0.008672266109234, 0.007386133919777, 0.006424492429929,
      0.005680041266617, 0.005087596677321, 0.004605438587838, 0.004205702655407, 0.003869108400092, 0.003581910933429,
      0.003334058876755, 0.003118038450500, 0.555776428234163, 0.096887182280156, 0.041654915323137, 0.024490963192973,
      0.016816414387792, 0.012631203375148, 0.010046934132754, 0.008310873007581, 0.007071904546731, 0.006146755001044,
      0.005431331610505, 0.004862502193998, 0.004399916077612, 0.004016655968932, 0.003694116223283, 0.003419045247043,
      0.003181760285193, 0.002975028358738, 0.550983631882507, 0.094834703139416, 0.040451127126301, 0.023665724312681,
      0.016197241569264, 0.012139099515104, 0.009640103707504, 0.007964849326159, 0.006771255794494, 0.005881175919801,
      0.005193623668744, 0.004647451131478, 0.004203630515231, 0.003836158445184, 0.003527079902161, 0.003263618502716,
      0.003036446343841, 0.002838600581102, 0.546291971785148, 0.092841394632250, 0.039287837054570, 0.022870837329790,
      0.015602172536927, 0.011666910841688, 0.009250204126935, 0.007633531005251, 0.006483593290973, 0.005627217464527,
      0.004966426640423, 0.004441992164996, 0.004016164234204, 0.003663821446239, 0.003367636084909, 0.003115289567695,
      0.002897795594927, 0.002708451201737, 0.541698539406959, 0.090905269295256, 0.038163557232277, 0.022105121600046,
      0.015030232391994, 0.011213808926822, 0.008876516116448, 0.007316283041700, 0.006208348945951, 0.005384365898041,
      0.004749271776791, 0.004245694374794, 0.003837118552768, 0.003499274081969, 0.003215438080784, 0.002973733008153,
      0.002765501421882, 0.002584290376660, 0.537200523750482, 0.089024412603105, 0.037076857614420, 0.021367444203416,
      0.014480486771073, 0.010779000524101, 0.008518351452072, 0.007012498204346, 0.005944979779010, 0.005152130390160,
      0.004541711384368, 0.004058146318729, 0.003666112907004, 0.003342162396303, 0.003070155093430, 0.002838638364447,
      0.002639271362218, 0.002465841680551, 0.532795207780330, 0.087196980135830, 0.036026363650389, 0.020656717960663,
      0.013952040125323, 0.010361726050468, 0.008175051602777, 0.006721595806730, 0.005692966800097, 0.004930041988909,
      0.004343317873366, 0.003878955147218, 0.003502784023695, 0.003192148590935, 0.002931471489586, 0.002709709459920,
      0.002518826452280, 0.002352841483925, 0.528479964986347, 0.085421194860481, 0.035010754044758, 0.019971899534082,
      0.013444034073856, 0.009961258134250, 0.007845986432779, 0.006443020534513, 0.005451813940361, 0.004717652638221,
      0.004153682849341, 0.003707745758594, 0.003346785131095, 0.003048910284745, 0.002799086101550, 0.002586663741704,
      0.002403900602422, 0.002245038359290, 0.524252256079918, 0.083695344522402, 0.034028758611000, 0.019311987608780,
      0.012955645828137, 0.009576900226641, 0.007530552961153, 0.006176241325143, 0.005221047030988, 0.004514534240015,
      0.003972416246105, 0.003544159993000, 0.003197785205871, 0.002912139807287, 0.002672711561825, 0.002469231652038,
      0.002294240001176, 0.002142192514925, 0.520109625818075, 0.082017779141582, 0.033079156214201, 0.018676021151034,
      0.012486086684289, 0.009207985273818, 0.007228174176211, 0.005920750297415, 0.005000212827859, 0.004320277758639,
      0.003799145498021, 0.003387855863075, 0.003055468254556, 0.002781543524766, 0.002552073668475, 0.002357156028670,
      0.002189602547066, 0.002044075254990, 0.516049699950252, 0.080386908609722, 0.032160772799005, 0.018063077740383,
      0.012034600580324, 0.008853874447031, 0.006938297902193, 0.005676061728713, 0.004788878079975, 0.004134492365743,
      0.003633514749873, 0.003238506819737, 0.002919532627914, 0.002656841197017, 0.002436910779760, 0.002250191533017,
      0.002089757306786, 0.001950468464753, 0.512070182282748, 0.078801200383852, 0.031272479499164, 0.017472271972312,
      0.011600462715479, 0.008513955928087, 0.006660395715956, 0.005441711077763, 0.004586628639662, 0.003956804623766,
      0.003475184102610, 0.003095801051450, 0.002789690366715, 0.002537765364045, 0.002326973236692, 0.002148104104767,
      0.001994483998514, 0.001861164119750, 0.508168851856189, 0.077259177272501, 0.030413190825266, 0.016902753928457,
      0.011182978228944, 0.008187643747807, 0.006393961911432, 0.005217254050879, 0.004393068612670, 0.003786857706263,
      0.003323828893316, 0.002959440815459, 0.002665666577477, 0.002424060760762, 0.002222022812229, 0.002050670441731,
      0.001903572499190, 0.001775963817770, 0.504343560231439, 0.075759415310596, 0.029581862927321, 0.016353707711456,
      0.010781480935389, 0.007874376675106, 0.006138512509716, 0.005002265709741, 0.004207819546375, 0.003624310653420,
      0.003179139007852, 0.002829141799504, 0.002547198836786, 0.002315483758633, 0.002121832185871, 0.001957677503734,
      0.001816822374650, 0.001694678332599, 0.500592228879628, 0.074300541719424, 0.028777491929038, 0.015824350041647,
      0.010395332114818, 0.007573617154463, 0.005893584312778, 0.004796339618849, 0.004030519654359, 0.003468837661141,
      0.003040818224679, 0.002704632512642, 0.002434036622891, 0.002211801832973, 0.002026184442480, 0.001868922039470,
      0.001734042431547, 0.001617127188510, 0.496912846672131, 0.072881232948157, 0.027999112330792, 0.015313928912968,
      0.010023919354391, 0.007284850289670, 0.005658733998841, 0.004599087030868, 0.003860823075720, 0.003320127402206,
      0.002908583588443, 0.002585653703831, 0.002325940773315, 0.002112793054727, 0.001934872594205, 0.001784210135222,
      0.001655050290050, 0.001543138254534, 0.493303467466497, 0.071500212793565, 0.027245795478358, 0.014821722305516,
      0.009666655439932, 0.007007582871801, 0.005433537257582, 0.004410136108179, 0.003698399167563, 0.003177882378020,
      0.002782164811959, 0.002471957806999, 0.002222682967292, 0.002018245605596, 0.001847699124452, 0.001703356784455,
      0.001579671976358, 0.001472547357584, 0.489762207784512, 0.070156250594710, 0.026516648094651, 0.014347036952326,
      0.009322977294987, 0.006741342449456, 0.005217587963382, 0.004229131179001, 0.003542931829163, 0.003041818299595,
      0.002661303705304, 0.002363308411410, 0.002124045231886, 0.001927957315434, 0.001764475552867, 0.001626185477309,
      0.001507741534096, 0.001405197913564, 0.486287244578707, 0.068848159499536, 0.025810810871829, 0.013889207158065,
      0.008992344965343, 0.006485676439441, 0.005010497384943, 0.004055732026545, 0.003394118856377, 0.002911663496414,
      0.002545753630793, 0.002259479756143, 0.002029819470701, 0.001841735220893, 0.001685022020367, 0.001552527809064,
      0.001439100653734, 0.001340940575615, 0.482876813083813, 0.067574794800390, 0.025127457121201, 0.013447593667389,
      0.008674240647043, 0.006240151276083, 0.004811893429647, 0.003889613209723, 0.003251671324939, 0.002787158351928,
      0.002435278982636, 0.002160256247594, 0.001939807014153, 0.001759395144332, 0.001609166893292, 0.001482223106703,
      0.001373598319173, 0.001279632898708, 0.479529204749764, 0.066335052335660, 0.024465791478527, 0.013021582580864,
      0.008368167756001, 0.006004351597487, 0.004621419921121, 0.003730463413990, 0.003115313001344, 0.002668054764477,
      0.002329654690182, 0.002065431998952, 0.001853818190307, 0.001680761292063, 0.001536746385798, 0.001415118072730,
      0.001311090470716, 0.001221139019817, 0.476242765253032, 0.065127866954810, 0.023825048662368, 0.012610584316399,
      0.008073650037422, 0.005777879467128, 0.004438735908541, 0.003577984830987, 0.002984779780070, 0.002554115632481,
      0.002228665743650, 0.001974810390635, 0.001771671915334, 0.001605665871031, 0.001467604199635, 0.001351066445456,
      0.001251439683649, 0.001165329352954, 0.473015892583175, 0.063952211044218, 0.023204492283274, 0.012214032614246,
      0.007790230713284, 0.005560353629201, 0.004263515006252, 0.003431892565682, 0.002859819145966, 0.002445114362805,
      0.002132106741337, 0.001888203650726, 0.001693195302698, 0.001533948723101, 0.001401591180529, 0.001289928674966,
      0.001194514861720, 0.001112080298361, 0.469847035201624, 0.062807093111329, 0.022603413701667, 0.011831383583717,
      0.007517471666243, 0.005351408796266, 0.004095444762373, 0.003291914069790, 0.002740189660648, 0.002340834401243,
      0.002039781457328, 0.001805432454510, 0.001618223290203, 0.001465456976108, 0.001338564990371, 0.001231571614065,
      0.001140190944813, 0.001061273965214, 0.466734690269854, 0.061691556424734, 0.022021130932389, 0.011462114789818,
      0.007254952658378, 0.005150694967770, 0.003934226055102, 0.003157788600295, 0.002625660471856, 0.002241068784130,
      0.001951502428764, 0.001726325542225, 0.001546598284076, 0.001400044710938, 0.001278389794520, 0.001175868223497,
      0.001088348630163, 0.001012797907205, 0.463677401944206, 0.060604677707905, 0.021456987593947, 0.011105724378117,
      0.007002270583268, 0.004957876778072, 0.003779572515486, 0.003029266701948, 0.002516010844712, 0.002145619710113,
      0.001867090561785, 0.001650719354193, 0.001478169819325, 0.001337572643864, 0.001220935963478, 0.001122697290775,
      0.001038874106489, 0.000966544870396, 0.460673759734724, 0.059545565884379, 0.020910351900600, 0.010761730236198,
      0.006759038749950, 0.004772632872711, 0.003631209975488, 0.002906109712684, 0.002411029713905, 0.002054298131183,
      0.001786374755303, 0.001578457682550, 0.001412794235601, 0.001277907823463, 0.001166079788322, 0.001071943161996,
      0.000991658800427, 0.000922412552785, 0.457722396925484, 0.058513360872317, 0.020380615695488, 0.010429669190149,
      0.006524886197392, 0.004594655311647, 0.003488875940241, 0.002788089290917, 0.002310515255868, 0.001966923362082,
      0.001709191541787, 0.001509391338809, 0.001350334367875, 0.001220923341425, 0.001113703209208, 0.001023495486040,
      0.000946599134715, 0.000880303375033, 0.454821989054011, 0.057507232426400, 0.019867193523092, 0.010109096234594,
      0.006299457038153, 0.004423648998316, 0.003352319083395, 0.002674986963750, 0.002214274480032, 0.001883322707271,
      0.001635384744285, 0.001443377836531, 0.001290659251229, 0.001166498056627, 0.001063693556390, 0.000977248970573,
      0.000903596297567, 0.000840124261827, 0.451971252447460, 0.056526379025170, 0.019369521739388, 0.009799583794838,
      0.006082409829961, 0.004259331133353, 0.003221298764568, 0.002566593695161, 0.002122122838321, 0.001803331104645,
      0.001564805148958, 0.001380281088415, 0.001233643839130, 0.001114516331852, 0.001015943303138, 0.000933103149319,
      0.000862556022728, 0.000801786433401, 0.449168942813318, 0.055570026801940, 0.018887057658110, 0.009500721019762,
      0.005873416974021, 0.004101430691923, 0.003095584567893, 0.002462709473278, 0.002033883852050, 0.001726790785257,
      0.001497310192403, 0.001319971117155, 0.001179168734569, 0.001064867782578, 0.000970349830031, 0.000890962160067,
      0.000823388379701, 0.000765205206722, 0.446413853882528, 0.054637428517529, 0.018419278731639, 0.009212113104152,
      0.005672164138880, 0.003949687923604, 0.002974955860749, 0.002363142915889, 0.001949388755453, 0.001653550948320,
      0.001432763663111, 0.001262323779427, 0.001127119933463, 0.001017447037278, 0.000926815200090, 0.000850734532923,
      0.000786007573694, 0.000730299805910, 0.443704816102947, 0.053727862573109, 0.017965681765056, 0.008933380639227,
      0.005478349708760, 0.003803853873854, 0.002859201371773, 0.002267710893376, 0.001868476155095, 0.001583467450788,
      0.001371035416403, 0.001207220502414, 0.001077388579765, 0.000972153508702, 0.000885245944250, 0.000812332988335,
      0.000750331754812, 0.000696993181443, 0.441040695381211, 0.052840632061560, 0.017525782161993, 0.008664158990154,
      0.005291684255293, 0.003663689926113, 0.002748118787313, 0.002176238168292, 0.001790991704456, 0.001516402510878,
      0.001312001102242, 0.001154548032294, 0.001029870731754, 0.000928891175641, 0.000845552856696, 0.000775674244433,
      0.000716282836086, 0.000665211837762, 0.438420391871084, 0.051975063855752, 0.017099113200958, 0.008404097699420,
      0.005111890031660, 0.003528967363635, 0.002641514365493, 0.002088557050853, 0.001716787793001, 0.001452224424879,
      0.001255541905334, 0.001104198194139, 0.000984467138976, 0.000887568374681, 0.000807650799615, 0.000740678833253,
      0.000683786319918, 0.000634885668864, 0.435842838806517, 0.051130507732287, 0.016685225340866, 0.008152859914947,
      0.004938700488155, 0.003399466950199, 0.002539202567126, 0.002004507069618, 0.001645723249087, 0.001390807296673,
      0.001201544296955, 0.001056067662718, 0.000941083029361, 0.000848097601497, 0.000771458516910, 0.000707270925442,
      0.000652771132553, 0.000605947801521, 0.433307001377653, 0.050306335529253, 0.016283685554567, 0.007910121841925,
      0.004771859808268, 0.003274978528862, 0.002441005702726, 0.001923934656697, 0.001577663056097, 0.001332030779372,
      0.001149899797979, 0.001010057743698, 0.000899627906049, 0.000810395321248, 0.000736898456487, 0.000675378163044,
      0.000623169466206, 0.000578334445775, 0.430811875648124, 0.049501940336625, 0.015894076689214, 0.007675572217332,
      0.004611122464387, 0.003155300637977, 0.002346753594912, 0.001846692846831, 0.001512478081190, 0.001275779828541,
      0.001100504752596, 0.000966074164765, 0.000860015353468, 0.000774381787644, 0.000703896600705, 0.000644931499994,
      0.000594916628494, 0.000551984752360, 0.428356487512047, 0.048716735718006, 0.015515996852355, 0.007448911806209,
      0.004456252792290, 0.003040240143719, 0.002256283255518, 0.001772640989734, 0.001450044816119, 0.001221944466474,
      0.001053260112235, 0.000924026876228, 0.000822162852263, 0.000739980870304, 0.000672382304611, 0.000615865049968,
      0.000567950898826, 0.000526840676735, 0.425939891689163, 0.047950154962437, 0.015149058822687, 0.007229852918748,
      0.004307024583609, 0.002929611888405, 0.002169438576773, 0.001701644475099, 0.001390245129568, 0.001170419557030,
      0.001008071229229, 0.000883829860652, 0.000785991602649, 0.000707119890004, 0.000642288141617, 0.000588115941243,
      0.000542213391432, 0.000502846849426, 0.423561170756673, 0.047201650365080, 0.014792889484455, 0.007018118947333,
      0.004163220695500, 0.002823238353910, 0.002086070035923, 0.001633574469715, 0.001332966030489, 0.001121104590555,
      0.000964847659778, 0.000845400951133, 0.000751426355813, 0.000675729461476, 0.000613549756249, 0.000561624178246,
      0.000517647924730, 0.000479950452371, 0.421219434216328, 0.046470692535598, 0.014447129284514, 0.006813443922671,
      0.004024632676776, 0.002720949339536, 0.002006034412707, 0.001568307666151, 0.001278099441953, 0.001073903478425,
      0.000923502975796, 0.000808661657801, 0.000718395252994, 0.000645743343386, 0.000586105723676, 0.000536332509491,
      0.000494200896723, 0.000458101100999, 0.418913817595444, 0.045756769733152, 0.014111431711130, 0.006615572088231,
      0.003891060409800, 0.002622581653686, 0.001929194519118, 0.001505726042488, 0.001225541985047, 0.001028724356787,
      0.000883954585232, 0.000773537002187, 0.000686829671900, 0.000617098295181, 0.000559897415673, 0.000512186301587,
      0.000471821166158, 0.000437250731776, 0.416643481580516, 0.045059387226919, 0.013785462793613, 0.006424257492202,
      0.003762311767464, 0.002527978818770, 0.001855418940915, 0.001445716632627, 0.001175194772352, 0.000985479399082,
      0.000846123560481, 0.000739955359094, 0.000656664080103, 0.000589733940482, 0.000534868872751, 0.000489133419059,
      0.000450459939178, 0.000417353494959, 0.414407611182190, 0.044378066681136, 0.013468900621929, 0.006239263596241,
      0.003638202284597, 0.002436990788739, 0.001784581790360, 0.001388171306684, 0.001126963210594, 0.000944084636950,
      0.000809934474527, 0.000707848305628, 0.000627835895119, 0.000563592636717, 0.000510966682144, 0.000467124109705,
      0.000430070661211, 0.000398365652319, 0.412205414930393, 0.043712345563669, 0.013161434885476, 0.006060362900298,
      0.003518554843202, 0.002349473678715, 0.001716562469699, 0.001332986561028, 0.001080756812042, 0.000904459789139,
      0.000775315244458, 0.000677150477064, 0.000600285350852, 0.000538619350724, 0.000488139861404, 0.000446110895227,
      0.000410608913848, 0.000380245479602, 0.410036124098465, 0.043061776577196, 0.012862766430220, 0.005887336582855,
      0.003403199370908, 0.002265289506197, 0.001651245444905, 0.001280063317551, 0.001036489014271, 0.000866528098072,
      0.000742196982024, 0.000647799429243, 0.000573955370112, 0.000514761540028, 0.000466339747335, 0.000426048466900,
      0.000392032316491, 0.000362953173512, 0.407898991955183, 0.042425927112073, 0.012572606833460, 0.005719974155915,
      0.003291972552097, 0.002184305943317, 0.001588520029233, 0.001229306731737, 0.000994077007911, 0.000830216173701,
      0.000710513850917, 0.000619735507188, 0.000548791442939, 0.000491969039551, 0.000445519890020, 0.000406893586051,
      0.000374300432541, 0.000346450762995, 0.405793293043624, 0.041804378720053, 0.012290677995478, 0.005558073134133,
      0.003184717551138, 0.002106396079697, 0.001528280176163, 0.001180626009158, 0.000953441572026, 0.000795453844350,
      0.000680202930476, 0.000592901719671, 0.000524741510453, 0.000470193953496, 0.000425635951712, 0.000388604989107,
      0.000357374679915, 0.000330702024638, 0.403718322485826, 0.041196726608001, 0.012016711747394, 0.005401438717495,
      0.003081283747234, 0.002031438195425, 0.001470424281319, 0.001133934230013, 0.000914506916789, 0.000762174014207,
      0.000651204085511, 0.000567243619450, 0.000501755853992, 0.000449390552166, 0.000406645610360, 0.000371143297023,
      0.000341218245698, 0.000315672401982, 0.401673395312296, 0.040602579150832, 0.011750449474564, 0.005249883486968,
      0.002981526480366, 0.001959315543719, 0.001414854992957, 0.001089148181362, 0.000877200533114, 0.000730312527186,
      0.000623459841991, 0.000542709188927, 0.000479786989288, 0.000429515173503, 0.000388508467563, 0.000354470928877,
      0.000325796004738, 0.000301328928570, 0.399657845815393, 0.040021557422905, 0.011491641754877, 0.005103227112598,
      0.002885306807884, 0.001889916142860, 0.001361479030667, 0.001046188196709, 0.000841453048950, 0.000699808036860,
      0.000596915268307, 0.000519248730971, 0.000458789565451, 0.000410526129116, 0.000371185960737, 0.000338552019433,
      0.000311074441988, 0.000287640154560, 0.397671026925699, 0.039453294747136, 0.011240048011349, 0.004961296073504,
      0.002792491271260, 0.001823132576999, 0.001310207011913, 0.001004978002608, 0.000807198091926, 0.000670601882190,
      0.000571517861875, 0.000496814764680, 0.000438720268544, 0.000392383614602, 0.000354641279312, 0.000323352340505,
      0.000297021578445, 0.000274576076732, 0.395712309610514, 0.038897436261143, 0.010995436178434, 0.004823923389309,
      0.002702951672607, 0.001758861805438, 0.001260953286073, 0.000965444571987, 0.000774372158072, 0.000642637968814,
      0.000547217440837, 0.000475361925850, 0.000419537729546, 0.000375049623957, 0.000338839284775, 0.000308839225938,
      0.000283606900499, 0.000262108071738, 0.393781082293616, 0.038353638499739, 0.010757582381481, 0.004690948362490,
      0.002616564860512, 0.001697004980035, 0.001213635775656, 0.000927517983889, 0.000742914486344, 0.000615862655603,
      0.000523966040614, 0.000454846871945, 0.000401202436482, 0.000358487867901, 0.000323746434364, 0.000294981500030,
      0.000270801292541, 0.000250208832435, 0.391876750295519, 0.037821568993138, 0.010526270628813, 0.004562216331231,
      0.002533212524814, 0.001637467270373, 0.001168175824371, 0.000891131289344, 0.000712766938691, 0.000590224646298,
      0.000501717815112, 0.000435228191356, 0.000383676650557, 0.000342663695922, 0.000309330708273, 0.000281749409255,
      0.000258576972687, 0.000238852307171, 0.389998735293426, 0.037300905880251, 0.010301292515914, 0.004437578432304,
      0.002452780999930, 0.001580157696362, 0.001124498051761, 0.000856220383112, 0.000683873885421, 0.000565674885954,
      0.000480428942357, 0.000416466316763, 0.000366924326087, 0.000327544021878, 0.000295561540180, 0.000269114557125,
      0.000246907431459, 0.000228013641873, 0.388146474800162, 0.036791337536481, 0.010082446941219, 0.004316891373590,
      0.002375161076382, 0.001524988967949, 0.001082530214104, 0.000822723881019, 0.000656182095627, 0.000542166462009,
      0.000460057534366, 0.000398523442397, 0.000350911034075, 0.000313097252992, 0.000282409750963, 0.000257049842035,
      0.000235767373300, 0.000217669124822, 0.386319421661361, 0.036292562215446, 0.009869539833062, 0.004200017215822,
      0.002300247820150, 0.001471877331635, 0.001042203071315, 0.000790583002658, 0.000629640632446, 0.000519654509753,
      0.000440563551060, 0.000381363445054, 0.000335603889244, 0.000299293222088, 0.000269847485448, 0.000245529397981,
      0.000225132660784, 0.000207796133970, 0.384517043570234, 0.035804287704087, 0.009662383887297, 0.004086823163168,
      0.002227940399553, 0.001420742423497, 0.001003450259581, 0.000759741459211, 0.000604200752949, 0.000498096122002,
      0.000421908718031, 0.000364951808655, 0.000320971480392, 0.000286103122915, 0.000257848152055, 0.000234528537990,
      0.000214980261394, 0.000198373086707, 0.382738822599254, 0.035326230990640, 0.009460798315188, 0.003977181362292,
      0.002158141919312, 0.001371507128448, 0.000966208169487, 0.000730145346161, 0.000579815812431, 0.000477450262795,
      0.000404056448009, 0.000349255552222, 0.000306983803897, 0.000273499448413, 0.000246386365209, 0.000224023700154,
      0.000205288196763, 0.000189379391942, 0.380984254748137, 0.034858117944956, 0.009264608601135, 0.003870968709543,
      0.002090759261499, 0.001324097445451, 0.000930415829385, 0.000701743040684, 0.000556441172930, 0.000457677684926,
      0.000386971765832, 0.000334243161098, 0.000293612200244, 0.000261455931804, 0.000235437890378, 0.000213992396145,
      0.000196035494247, 0.000180795404403, 0.379252849507499, 0.034399683010711, 0.009073646269851, 0.003768066665918,
      0.002025702933089, 0.001278442358445, 0.000896014793775, 0.000674485103515, 0.000534034115778, 0.000438740851144,
      0.000370621236787, 0.000319884521265, 0.000280829293425, 0.000249947490356, 0.000224979591638, 0.000204413162088,
      0.000187202140732, 0.000172602381045, 0.377544129437629, 0.033950668909015, 0.008887748662598, 0.003668361079502,
      0.001962886919823, 0.001234473712722, 0.000862949036498, 0.000648324185092, 0.000512553757992, 0.000420603858855,
      0.000354972898162, 0.000306150856630, 0.000268608933094, 0.000238950171717, 0.000214989381627, 0.000195265511701,
      0.000178769038569, 0.000164782439483, 0.375857629761790, 0.033510826353001, 0.008706758722120, 0.003571742015056,
      0.001902228546119, 0.001192126096540, 0.000831164848510, 0.000623214935780, 0.000491960972362, 0.000403232368169,
      0.000339996193850, 0.000293014669133, 0.000256926139332, 0.000228441102688, 0.000205446173792, 0.000186529891575,
      0.000170717963545, 0.000157318518333, 0.374192897973536, 0.033079913772951, 0.008530524785926, 0.003478103590460,
      0.001843648340780, 0.001151336727726, 0.000800610740045, 0.000599113920009, 0.000472218311041, 0.000386593533144,
      0.000325661911898, 0.000280449681551, 0.000245757049927, 0.000218398440330, 0.000196329836823, 0.000178187638517,
      0.000163031524780, 0.000150194339399, 0.372549493457515, 0.032657697051558, 0.008358900387589, 0.003387343819738,
      0.001787069908252, 0.001112045345076, 0.000771237346984, 0.000575979534149, 0.000453289932507, 0.000370655936080,
      0.000311942124832, 0.000268430782882, 0.000235078870032, 0.000208801325300, 0.000187621151167, 0.000170220938844,
      0.000155693126482, 0.000143394371592, 0.370926987123262, 0.032243949268929, 0.008191744065729, 0.003299364462384,
      0.001732419805194, 0.001074194104323, 0.000742997341235, 0.000553771927946, 0.000435141531729, 0.000355389524731,
      0.000298810132668, 0.000256933976189, 0.000224869824111, 0.000199629837307, 0.000179301767529, 0.000162612789549,
      0.000148686931453, 0.000136903796533, 0.369324961051508, 0.031838450456956, 0.008028919180379, 0.003214070878729,
      0.001679627422139, 0.001037727478499, 0.000715845344955, 0.000532452929383, 0.000417740273411, 0.000340765552303,
      0.000286240408456, 0.000245936328789, 0.000215109110055, 0.000190864952598, 0.000171354167272, 0.000155346961251,
      0.000141997826285, 0.000130708475729, 0.367743008152541, 0.031440987362687, 0.007870293736438, 0.003131371891105,
      0.001628624870031, 0.001002592162492, 0.000689737848451, 0.000511985972801, 0.000401054728158, 0.000326756520113,
      0.000274208546273, 0.000235415924686, 0.000205776855380, 0.000182488503374, 0.000163761624630, 0.000148407962836,
      0.000135611388151, 0.000124794919277, 0.366180731836188, 0.031051353220350, 0.007715740213928, 0.003051179650560,
      0.001579346871439, 0.000968736981612, 0.000664633131593, 0.000492336030139, 0.000385054811453, 0.000313336122787,
      0.000262691211527, 0.000225351819144, 0.000196854075403, 0.000174483139058, 0.000156508170640, 0.000141781007733,
      0.000129513853132, 0.000119150256005, 0.364637745692975, 0.030669347531712, 0.007565135404777, 0.002973409508908,
      0.001531730656236, 0.000936112804013, 0.000640491188592, 0.000473469545163, 0.000369711725311, 0.000300479195892,
      0.000251666093490, 0.000215723995298, 0.000188322633312, 0.000166832289325, 0.000149578558728, 0.000135451981725,
      0.000123692085992, 0.000113762204998, 0.363113673186079, 0.030294775854426, 0.007418360255871, 0.002897979895878,
      0.001485715861575, 0.000904672456787, 0.000617273655987, 0.000455354370540, 0.000354997902495, 0.000288161665883,
      0.000241111859946, 0.000206513322728, 0.000180165202045, 0.000159520128800, 0.000142958231866, 0.000129407412243,
      0.000118133551356, 0.000108619048433, 0.361608147353659, 0.029927449598083, 0.007275299718132, 0.002824812201166,
      0.001441244435958, 0.000874370645580, 0.000594943743726, 0.000437959707647, 0.000340886953170, 0.000276360502269,
      0.000231008113871, 0.000197701517885, 0.000172365227891, 0.000152531543383, 0.000136633291224, 0.000123634439068,
      0.000112826286204, 0.000103709605671, 0.360120810521208, 0.029567185827663, 0.007135842601363, 0.002753830661184,
      0.001398260547256, 0.000845163877588, 0.000573466169172, 0.000421256048987, 0.000327353613910, 0.000265053671899,
      0.000221335352033, 0.000189271106298, 0.000164906895732, 0.000145852098083, 0.000130590466261, 0.000118120786370,
      0.000107758873633, 0.000099023208547, 0.358651314023546, 0.029213807074103, 0.006999881434658, 0.002684962250312,
      0.001356710494483, 0.000817010387776, 0.000552807093938, 0.000405215123091, 0.000314373698925, 0.000254220095265,
      0.000212074925455, 0.000181205386481, 0.000157775095867, 0.000139468006327, 0.000124817086176, 0.000112854736035,
      0.000102920417835, 0.000094549677793, 0.357199317936115, 0.028867141151714, 0.006867312332122, 0.002618136576469,
      0.001316542623191, 0.000789870068185, 0.000532934063415, 0.000389809841817, 0.000301924053433, 0.000243839604741,
      0.000203209001630, 0.000173488395463, 0.000150955392326, 0.000133366100660, 0.000119301052664, 0.000107825102213,
      0.000098300520209, 0.000090279300558, 0.355764490815237, 0.028527020982173, 0.006738034863728, 0.002553285780834,
      0.001277707244318, 0.000763704400207, 0.000513815948869, 0.000375014249920, 0.000289982509073, 0.000233892904673,
      0.000194720528426, 0.000166104875859, 0.000144433992618, 0.000127533804771, 0.000114030813913, 0.000103021207032,
      0.000093889256588, 0.000086202808960, 0.354346509447012, 0.028193284424863, 0.006611951931069, 0.002490344441529,
      0.001240156556358, 0.000738476389680, 0.000495422892020, 0.000360803476802, 0.000278527841265, 0.000224361533222,
      0.000186593199599, 0.000159040244428, 0.000138197718850, 0.000121959106801, 0.000108995339792, 0.000098432857419,
      0.000089677155501, 0.000082311359623, 0.352945058604552, 0.027865774113296, 0.006488969647846, 0.002429249481127,
      0.001203844570703, 0.000714150504706, 0.000477726251964, 0.000347153690357, 0.000267539728438, 0.000215227825907,
      0.000178811421841, 0.000152280562035, 0.000132233980140, 0.000116630533848, 0.000104184098164, 0.000094050322993,
      0.000085655177442, 0.000078596514171, 0.351559830813235, 0.027544337297400, 0.006368997224877, 0.002369940077805,
      0.001168727040029, 0.000690692616063, 0.000460698554365, 0.000334042052794, 0.000256998713043, 0.000206474880746,
      0.000171360283300, 0.000145812504964, 0.000126530746292, 0.000111537127641, 0.000099587032287, 0.000089864314965,
      0.000081814695090, 0.000075050220607, 0.350190526123712, 0.027228825691448, 0.006251946859466, 0.002312357580012,
      0.001134761389603, 0.000668069940098, 0.000444313442805, 0.000321446678382, 0.000246886164262, 0.000198086524944,
      0.000164225523496, 0.000139623337520, 0.000121076522639, 0.000106668421313, 0.000095194539240, 0.000085865965998,
      0.000078147474435, 0.000071664795560, 0.348836851892378, 0.026919095327414, 0.006137733628954, 0.002256445424499,
      0.001101906651376, 0.000646250984004, 0.000428545632194, 0.000309346593007, 0.000237184242354, 0.000190047283045,
      0.000157393504583, 0.000133700885854, 0.000115860326039, 0.000102014417228, 0.000090997449334, 0.000082046811004,
      0.000074645656773, 0.000068432907347, 0.347498522569045, 0.026615006413544, 0.006026275388281, 0.002202149057571,
      0.001070123400771, 0.000625205493382, 0.000413370864171, 0.000297721695480, 0.000227875864556, 0.000182342346498,
      0.000150851183883, 0.000128033512967, 0.000110871661943, 0.000097565565815, 0.000086987006463, 0.000078398768803,
      0.000071301741533, 0.000065347559816, 0.346175259491556, 0.026316423197958, 0.005917492671417, 0.002149415859449,
      0.001039373696018, 0.000604904401973, 0.000398765864382, 0.000286552720508, 0.000218944672467, 0.000174957544551,
      0.000144586087646, 0.000122610094830, 0.000106100502502, 0.000093312745361, 0.000083154849347, 0.000074914124638,
      0.000068108569884, 0.000062402076935, 0.344866790687104, 0.026023213837088, 0.005811308596494, 0.002098195071585,
      0.001009621019966, 0.000585319783496, 0.000384708301584, 0.000275821203261, 0.000210375000861, 0.000167879316446,
      0.000138586285979, 0.000117419997576, 0.000101537265660, 0.000089247242723, 0.000079492993634, 0.000071585513473,
      0.000065059309101, 0.000059590088089, 0.343572850680004, 0.025735250268762, 0.005707648774498, 0.002048437726844,
      0.000980830224244, 0.000566424805468, 0.000371176748473, 0.000265509445464, 0.000202151847857, 0.000161094684820,
      0.000132840368881, 0.000112453055717, 0.000097172795189, 0.000085360734907, 0.000075993814814, 0.000068405904071,
      0.000062147437649, 0.000056905514056, 0.342293180305706, 0.025452408089779, 0.005606441221393, 0.002000096582410,
      0.000952967475686, 0.000548193684963, 0.000358150644173, 0.000255600482953, 0.000194260846390, 0.000154591230295,
      0.000127337423354, 0.000107699551325, 0.000092998341627, 0.000081645271492, 0.000072650031911, 0.000065368583787,
      0.000059366730949, 0.000054342553630, 0.341027526530814, 0.025174566437779, 0.005507616273522, 0.001953126055320,
      0.000926000204922, 0.000530601646186, 0.000345610258315, 0.000246078054615, 0.000186688236921, 0.000148357067175,
      0.000122067011522, 0.000103150194160, 0.000089005544066, 0.000078093257841, 0.000069454691921, 0.000062467144067,
      0.000056711247799, 0.000051895670857, 0.339775642278912, 0.024901607877284, 0.005411106506177, 0.001907482160518,
      0.000899897057043, 0.000513624879819, 0.000333536656630, 0.000236926572672, 0.000179420841341, 0.000142380820217,
      0.000117019149726, 0.000098796102677, 0.000085186412774, 0.000074697439082, 0.000066401154954, 0.000059695466604,
      0.000054175317418, 0.000049559582863, 0.338537286261987, 0.024633418289714, 0.005316846655200, 0.001863122451320,
      0.000874627844258, 0.000497240504045, 0.000321911667999, 0.000228131094239, 0.000172446038006, 0.000136651602423,
      0.000112184288550, 0.000094628785886, 0.000081533312589, 0.000071450884803, 0.000063483080057, 0.000057047710128,
      0.000051753527084, 0.000047329248239, 0.337312222817249, 0.024369886767261, 0.005224773541502, 0.001820005962208,
      0.000850163500464, 0.000481426527183, 0.000310717852896, 0.000219677294102, 0.000165751737860, 0.000131158993810,
      0.000107553293721, 0.000090640126021, 0.000078038947060, 0.000068346974452, 0.000060694411679, 0.000054518297800,
      0.000049440710343, 0.000045199855963, 0.336100221749171, 0.024110905510466, 0.005134825998396, 0.001778093153834,
      0.000826476037636, 0.000466161811873, 0.000299938473159, 0.000211551438660, 0.000159326361597, 0.000125893021112,
      0.000103117427873, 0.000086822361988, 0.000074696343312, 0.000065379383389, 0.000058029366749, 0.000052101905182,
      0.000047231935752, 0.000043166814829, 0.334901058176559, 0.023856369729357, 0.005046944801608, 0.001737345860178,
      0.000803538503978, 0.000451426040732, 0.000289557463035, 0.000203740360993, 0.000153158817810, 0.000120844138373,
      0.000098868333109, 0.000083168073553, 0.000071498837583, 0.000062542069573, 0.000055482422352, 0.000049793448760,
      0.000045122496144, 0.000041225743371, 0.333714512384474, 0.023606177548026, 0.004961072601899, 0.001697727237744,
      0.000781324943752, 0.000437199683437, 0.000279559401454, 0.000196231436992, 0.000147238482098, 0.000116003208394,
      0.000094798014339, 0.000079670166227, 0.000068440061413, 0.000059829260852, 0.000053048303951, 0.000047588074984,
      0.000043107898382, 0.000039372460235, 0.332540369680845, 0.023360229912506, 0.004877153860159, 0.001659201716727,
      0.000759810358720, 0.000423463965156, 0.000269929485454, 0.000189012562510, 0.000141555177069, 0.000111361484994,
      0.000090898823353, 0.000076321856841, 0.000065513928456, 0.000057235442827, 0.000050721974160, 0.000045481149814,
      0.000041183853590, 0.000037602975013, 0.331378420257612, 0.023118430501839, 0.004795134784910, 0.001621734954073,
      0.000738970671133, 0.000410200836286, 0.000260653504743, 0.000182072131493, 0.000136099153223, 0.000106910596047,
      0.000087163443608, 0.000073116659752, 0.000062714621887, 0.000054755347263, 0.000048498622012, 0.000043468248745,
      0.000039346267816, 0.000035913479484, 0.330228459056234, 0.022880685642212, 0.004714963272101, 0.001585293788342,
      0.000718782688196, 0.000397392943430, 0.000251717817316, 0.000175399015051, 0.000130861070652, 0.000102642527271,
      0.000083584875679, 0.000070048373672, 0.000060036582366, 0.000052383941026, 0.000046373652731, 0.000041545147279,
      0.000037591233145, 0.000034300339264, 0.329090285637414, 0.022646904224039, 0.004636588847122, 0.001549846196329,
      0.000699224067963, 0.000385023601574, 0.000243109326101, 0.000168982541412, 0.000125831981541, 0.000098549606714,
      0.000080156423360, 0.000067111069084, 0.000057474496558, 0.000050116415524, 0.000044342677952, 0.000039707811845,
      0.000035915019198, 0.000032760085840, 0.327963704054900, 0.022416997621898, 0.004559962608938, 0.001515361251335,
      0.000680273286586, 0.000373076767403, 0.000234815456591, 0.000162812476746, 0.000121003313425, 0.000094624489933,
      0.000076871680379, 0.000064299076217, 0.000055023286151, 0.000047948176615, 0.000042401506398, 0.000037952391120,
      0.000034314065035, 0.000031289408971, 0.326848522733229, 0.022190879617202, 0.004485037176276, 0.001481809083063,
      0.000661909606876, 0.000361537013710, 0.000226824135405, 0.000156879006796, 0.000116366853169, 0.000090860145818,
      0.000073724517700, 0.000061606973561, 0.000052678097383, 0.000045874834982, 0.000040546134976, 0.000036275207752,
      0.000032784971425, 0.000029885149434, 0.325744554349252, 0.021968466323522, 0.004411766635773, 0.001449160839036,
      0.000644113048116, 0.000350389504859, 0.000219123769752, 0.000151172719300, 0.000111914731647, 0.000087249843031,
      0.000070709071383, 0.000059029576886, 0.000050434291030, 0.000043892196934, 0.000038772740276, 0.000034672750465,
      0.000031324493465, 0.000028544292109, 0.324651615717350, 0.021749676114440, 0.004340106492022, 0.001417388647511,
      0.000626864357073, 0.000339619973252, 0.000211703227757, 0.000145684587158, 0.000107639409072, 0.000083787137049,
      0.000067819730989, 0.000056561928764, 0.000048287432843, 0.000041996255630, 0.000037077670458, 0.000033141666516,
      0.000029929533544, 0.000027263959379, 0.323569527678179, 0.021534429553869, 0.004270013619431, 0.001386465581800,
      0.000610144980164, 0.000329214696755, 0.000204551819610, 0.000140405952319, 0.000103533660976, 0.000080465857763,
      0.000065051128490, 0.000054199288543, 0.000046233284422, 0.000040183182694, 0.000035457437506, 0.000031678754502,
      0.000028597134627, 0.000026041404836, 0.322498114990849, 0.021322649328724, 0.004201446215832, 0.001356365625963,
      0.000593937036730, 0.000319160477050, 0.000197659279499, 0.000135328510356, 0.000099590564783, 0.000077280097628,
      0.000062398127670, 0.000051937122779, 0.000044267794490, 0.000038449320208, 0.000033908709834, 0.000030280957489,
      0.000027324473844, 0.000024874007267, 0.321437206228409, 0.021114260183880, 0.004134363757776, 0.001327063641805,
      0.000578223293368, 0.000309444618873, 0.000191015748300, 0.000130444295688, 0.000095803486966, 0.000074224200328,
      0.000059855813999, 0.000049771096088, 0.000042387090565, 0.000036791173072, 0.000032428305230, 0.000028945356463,
      0.000026108856376, 0.000023759264932, 0.320386633676531, 0.020909188859314, 0.004068726957446, 0.001298535337138,
      0.000562987139277, 0.000300054910096, 0.000184611756996, 0.000125745667446, 0.000092166070763, 0.000071292749931,
      0.000057419484945, 0.000047697062407, 0.000040587471003, 0.000035205401702, 0.000031013184114, 0.000027669164071,
      0.000024947709619, 0.000022694790093, 0.319346233235276, 0.020707364029376, 0.004004497721120, 0.001270757235241,
      0.000548212562584, 0.000290979602622, 0.000178438210776, 0.000121225295922, 0.000088672224412, 0.000068480560527,
      0.000055084640723, 0.000045711056642, 0.000038865397395, 0.000033688815064, 0.000029660443107, 0.000026449718650,
      0.000023838577615, 0.000021678303794, 0.318315844323860, 0.020508716244083, 0.003941639109140, 0.001243706645492,
      0.000533884127604, 0.000282207394056, 0.000172486373809, 0.000116876149609, 0.000085316109907, 0.000065782666313,
      0.000052846975443, 0.000043809286692, 0.000037217487312, 0.000032238364014, 0.000028367308891, 0.000025284478533,
      0.000022779115738, 0.000020707630890, 0.317295309788288, 0.020313177872392, 0.003880115297320, 0.001217361635104,
      0.000519986952997, 0.000273727410128, 0.000166747854643, 0.000112691482780, 0.000082092132226, 0.000063194312107,
      0.000050702368657, 0.000041988125821, 0.000035640507366, 0.000030851134947, 0.000027131132345, 0.000024171016604,
      0.000021767085624, 0.000019780695295, 0.316284475811793, 0.020120683047359, 0.003819891539734, 0.001191701001937,
      0.000506506690784, 0.000265529187822, 0.000161214592216, 0.000108664823602, 0.000078994929030, 0.000060710944284,
      0.000048646877271, 0.000040244105378, 0.000034131366588, 0.000029524343724, 0.000025949382948, 0.000023107015111,
      0.000020800350330, 0.000018895515454, 0.315283191827960, 0.019931167613133, 0.003760934132844, 0.001166704248339,
      0.000493429506197, 0.000257602659197, 0.000155878842450, 0.000104789962746, 0.000076019360806, 0.000058328202104,
      0.000046676727818, 0.000038573907830, 0.000032687110103, 0.000028255329880, 0.000024819643437, 0.000022090260711,
      0.000019876869719, 0.000018050200018, 0.314291310436461, 0.019744569073700, 0.003703210380910, 0.001142351555970,
      0.000480742058312, 0.000249938135860, 0.000150733165398, 0.000101060942488, 0.000073160501427, 0.000056041909410,
      0.000044788309067, 0.000036974360112, 0.000031304913079, 0.000027041551091, 0.000023739604706, 0.000021118639739,
      0.000018994696055, 0.000017242943722, 0.313308687321307, 0.019560826543338, 0.003646688562632, 0.001118623761579,
      0.000468431481447, 0.000242526294070, 0.000145770412929, 0.000097472046261, 0.000070413629126, 0.000053848066702,
      0.000042978164965, 0.000035442427273, 0.000029982074961, 0.000025880577882, 0.000022707060932, 0.000020190133695,
      0.000018151969790, 0.000016472023449, 0.312335181171541, 0.019379880698699, 0.003591337898982, 0.001095502333691,
      0.000456485367289, 0.000235358160442, 0.000140983716918, 0.000094017788662, 0.000067774217854, 0.000051742843543,
      0.000041242987877, 0.000033975206394, 0.000028716013949, 0.000024770088587, 0.000021719904930, 0.000019302814937,
      0.000017346915563, 0.000015735794478, 0.311370653604281, 0.019201673732482, 0.003537128522183, 0.001072969350175,
      0.000444891747725, 0.000228425098237, 0.000136366477934, 0.000090692905868, 0.000065237929009, 0.000049722571301,
      0.000039579612138, 0.000032569920791, 0.000027504261736, 0.000023707864526, 0.000020776123711, 0.000018454842564,
      0.000016577838359, 0.000015032686899, 0.310414969090047, 0.019026149308615, 0.003484031445780, 0.001051007476642,
      0.000433639078336, 0.000221718794196, 0.000131912354388, 0.000087492346467, 0.000062800603520, 0.000047783736210,
      0.000037985007875, 0.000031223914462, 0.000026344458467, 0.000022691785401, 0.000019873794246, 0.000017644458497,
      0.000015843119863, 0.000014361202199, 0.309467994880286, 0.018853252518922, 0.003432018535780, 0.001029599945666,
      0.000422716222542, 0.000215231245909, 0.000127615252136, 0.000084411262670, 0.000060458254273, 0.000045922972724,
      0.000036456275109, 0.000029934646787, 0.000025234347935, 0.000021719824903, 0.000019011079415, 0.000016869983724,
      0.000015141214967, 0.000013719909996, 0.308529600937036, 0.018682929841196, 0.003381062482805, 0.001008730536771,
      0.000412112436369, 0.000208954749698, 0.000123469314508, 0.000081445001897, 0.000058207058858, 0.000044137057169,
      0.000034990638115, 0.000028699687462, 0.000024171772989, 0.000020790046510, 0.000018186224150, 0.000016129814731,
      0.000014470648445, 0.000013107444933, 0.307599659864639, 0.018515129098646, 0.003331136775231, 0.000988383557171,
      0.000401817353801, 0.000202881888982, 0.000119468912753, 0.000078589098715, 0.000056043352628, 0.000042422901661,
      0.000033585440027, 0.000027516711654, 0.000023154671140, 0.000019900599482, 0.000017397551741, 0.000015422420084,
      0.000013830011776, 0.000012522503705, 0.306678046843466, 0.018349799420661, 0.003282215673267, 0.000968543823228,
      0.000391820972705, 0.000197005523109, 0.000115608636872, 0.000075839267120, 0.000053963622058, 0.000040777548297,
      0.000032238137678, 0.000026383495378, 0.000022181070375, 0.000019049715026, 0.000016643460316, 0.000014746337170,
      0.000013217960110, 0.000011963842228, 0.305764639565554, 0.018186891204844, 0.003234274183950, 0.000949196642593,
      0.000382113641308, 0.000191318776642, 0.000111883286839, 0.000073191393138, 0.000051964498387, 0.000039198163595,
      0.000030946296674, 0.000025297911070, 0.000021249085145, 0.000018235702645, 0.000015922419478, 0.000014100169089,
      0.000012633209380, 0.000011430272938, 0.304859318172125, 0.018026356080282, 0.003187288037010, 0.000930327797023,
      0.000372686045195, 0.000185815029076, 0.000108287864170, 0.000070641527742, 0.000050042751524, 0.000037682033168,
      0.000029707586676, 0.000024257923363, 0.000020356912544, 0.000017456946640, 0.000015232967096, 0.000013482581681,
      0.000012074533536, 0.000010920662208, 0.303961965192900, 0.017868146871994, 0.003141233661586, 0.000911923525819,
      0.000363529194807, 0.000180487904956, 0.000104817563853, 0.000068185880064, 0.000048195284232, 0.000036226556643,
      0.000028519776895, 0.000023261585046, 0.000019502828650, 0.000016711902778, 0.000014573706235, 0.000012892300695,
      0.000011540761913, 0.000010433927889, 0.303072465487172, 0.017712217566527, 0.003096088163744, 0.000893970509893,
      0.000354634413432, 0.000175331264405, 0.000101467766602, 0.000065820810893, 0.000046419126540, 0.000034829242785,
      0.000027380731775, 0.000022307033200, 0.000018685185032, 0.000015999095102, 0.000013943302234, 0.000012328109078,
      0.000011030776716, 0.000009969036964, 0.302190706186563, 0.017558523278657, 0.003051829304785, 0.000876455856410,
      0.000345993325651, 0.000170339194023, 0.000098234031429, 0.000063542826449, 0.000044711430410, 0.000033487704844,
      0.000026288406881, 0.000021392485514, 0.000017902405416, 0.000015317112892, 0.000013340479907, 0.000011788844398,
      0.000010543510620, 0.000009525003303, 0.301316576639425, 0.017407020219148, 0.003008435480306, 0.000859367084010,
      0.000337597846238, 0.000165505998149, 0.000095112088525, 0.000061348572409, 0.000043069464623, 0.000032199656094,
      0.000025240844950, 0.000020516236751, 0.000017152982492, 0.000014664607753, 0.000012764020878, 0.000011273396375,
      0.000010077944478, 0.000009100885528, 0.300449968356824, 0.017257665663554, 0.002965885699981, 0.000842692108559,
      0.000329440169478, 0.000160826190478, 0.000092097832429, 0.000059234828195, 0.000041490609881, 0.000030962905575,
      0.000024236172122, 0.000019676655378, 0.000016435474869, 0.000014040290846, 0.000012212761027, 0.000010780704528,
      0.000009633105137, 0.000008695784972, 0.299590774960061, 0.017110417921997, 0.002924159568050, 0.000826419229436,
      0.000321512758911, 0.000156294486005, 0.000089187315482, 0.000057198501494, 0.000039972354131, 0.000029775354007,
      0.000023272594344, 0.000018872180348, 0.000015748504166, 0.000013442930230, 0.000011685588063, 0.000010309755926,
      0.000009208063347, 0.000008308843734, 0.298738892129676, 0.016965236309919, 0.002883237264469, 0.000810537116321,
      0.000313808337452, 0.000151905793294, 0.000086376741542, 0.000055236623008, 0.000038512288070, 0.000028634989892,
      0.000022348393912, 0.000018101318014, 0.000015090752229, 0.000012871348331, 0.000011181439196, 0.000009859583045,
      0.000008801931778, 0.000007939242818, 0.297894217555898, 0.016822081119753, 0.002843099526722, 0.000795034796460,
      0.000306319877906, 0.000147655207055, 0.000083662459962, 0.000053346341426, 0.000037108100854, 0.000027539885776,
      0.000021461926185, 0.000017362639186, 0.000014460958475, 0.000012324419526, 0.000010699298922, 0.000009429261721,
      0.000008413863111, 0.000007586200361, 0.297056650890482, 0.016680913593483, 0.002803727632245, 0.000779901642401,
      0.000299040593839, 0.000143538001015, 0.000081040959814, 0.000051524918601, 0.000035757575985, 0.000026488194678,
      0.000020611616426, 0.000016654776319, 0.000013857917351, 0.000011801067826, 0.000010238196899, 0.000009017909191,
      0.000008043048227, 0.000007248969945, 0.296226093699900, 0.016541695896081, 0.002765103381464, 0.000765127360183,
      0.000291963930795, 0.000139549621068, 0.000078508864341, 0.000049769724938, 0.000034458587375, 0.000025478146667,
      0.000019795956795, 0.000015976420818, 0.000013280475909, 0.000011300264674, 0.000009797205931, 0.000008624682232,
      0.000007688714477, 0.000006926838977, 0.295402449419841, 0.016404391089764, 0.002727209081400, 0.000750701977940,
      0.000285083557861, 0.000135685678703, 0.000076062925647, 0.000048078234963, 0.000033209095575, 0.000024508045594,
      0.000019013503455, 0.000015326320471, 0.000012727531486, 0.000010821026837, 0.000009375440028, 0.000008248775374,
      0.000007350124027, 0.000006619127157, 0.294585623310968, 0.016268963109063, 0.002690027529849, 0.000736615834935,
      0.000278393359537, 0.000131941944692, 0.000073700019596, 0.000046448023084, 0.000032007144163, 0.000023576265959,
      0.000018262873828, 0.000014703276988, 0.000012198029492, 0.000010362414385, 0.000008972052568, 0.000007889419202,
      0.000007026572283, 0.000006325184996, 0.293775522415920, 0.016135376736672, 0.002653542000084, 0.000722859570985,
      0.000271887427934, 0.000128314343018, 0.000071417140918, 0.000044876759530, 0.000030850856289, 0.000022681249913,
      0.000017542743948, 0.000014106143650, 0.000011690961289, 0.000009923528776, 0.000008586234531, 0.000007545878731,
      0.000006717386384, 0.000006044392427, 0.292972055517484, 0.016003597580040, 0.002617736226086, 0.000709424116267,
      0.000265560055253, 0.000124798945056, 0.000069211398518, 0.000043362206447, 0.000029738431364, 0.000021821504396,
      0.000016851845944, 0.000013533823064, 0.000011205362171, 0.000009503511012, 0.000008217212819, 0.000007217451858,
      0.000006421923762, 0.000005776157454, 0.292175133097928, 0.015873592048701, 0.002582594388275, 0.000696300681497,
      0.000259405726557, 0.000121391963977, 0.000067080010973, 0.000041902214175, 0.000028668141890, 0.000020995598386,
      0.000016188965625, 0.000012985265015, 0.000010740309430, 0.000009101539887, 0.000007864248647, 0.000006903467876,
      0.000006139570779, 0.000005519914884, 0.291384667299453, 0.015745327332291, 0.002548101099720, 0.000683480748461,
      0.000253419112820, 0.000118089749369, 0.000065020302216, 0.000040494717662, 0.000027638330425, 0.000020202160279,
      0.000015552940175, 0.000012459464412, 0.000010294920509, 0.000008716830306, 0.000007526636011, 0.000006603286068,
      0.000005869741408, 0.000005275125102, 0.290600571885719, 0.015618771379260, 0.002514241392812, 0.000670956060891,
      0.000247595064235, 0.000114888782077, 0.000063029697385, 0.000039137733039, 0.000026647406673, 0.000019439875374,
      0.000014942655948, 0.000011955459324, 0.000009868351235, 0.000008348631688, 0.000007203700221, 0.000006316294350,
      0.000005611875988, 0.000005041272911, 0.289822762204426, 0.015493892876221, 0.002481000706394, 0.000658718615662,
      0.000241928603772, 0.000111785669244, 0.000061105718855, 0.000037829354333, 0.000025693844700, 0.000018707483466,
      0.000014357046350, 0.000011472329105, 0.000009459794135, 0.000007996226432, 0.000006894796501, 0.000006041907991,
      0.000005365440034, 0.000004817866419, 0.289051155150906, 0.015370661227944, 0.002448364873313, 0.000646760654309,
      0.000236414920990, 0.000108777139540, 0.000059245982414, 0.000036567750310, 0.000024776180271, 0.000018003776545,
      0.000013795089830, 0.000011009192599, 0.000009068476821, 0.000007658928455, 0.000006599308656, 0.000005779568376,
      0.000005129923091, 0.000004604435981, 0.288285669132705, 0.015249046537959, 0.002416320108392, 0.000635074654851,
      0.000231049366068, 0.000105860038589, 0.000057448193609, 0.000035351161463, 0.000023893008289, 0.000017327596595,
      0.000013255807938, 0.000010565206425, 0.000008693660450, 0.000007336081798, 0.000006316647792, 0.000005528741836,
      0.000004904837653, 0.000004400533189, 0.287526224035115, 0.015129019589749, 0.002384852996808, 0.000623653323903,
      0.000225827444075, 0.000103031324568, 0.000055710144229, 0.000034177897110, 0.000023042980356, 0.000016677833482,
      0.000012738263488, 0.000010139563335, 0.000008334638254, 0.000007027059287, 0.000006046251098, 0.000005288918524,
      0.000004689718124, 0.000004205729906, 0.286772741187638, 0.015010551828512, 0.002353950482852, 0.000612489589068,
      0.000220744809445, 0.000100288063976, 0.000054029708936, 0.000033046332623, 0.000022224802423, 0.000016053422938,
      0.000012241558783, 0.000009731490648, 0.000007990734128, 0.000006731261263, 0.000005787580683, 0.000005059611347,
      0.000004484119827, 0.000004019617342, 0.286025143331342, 0.014893615343482, 0.002323599859060, 0.000601576591605,
      0.000215797260664, 0.000097627427577, 0.000052404842024, 0.000031954906761, 0.000021437232548, 0.000015453344627,
      0.000011764833928, 0.000009340248750, 0.000007661301289, 0.000006448114362, 0.000005540122462, 0.000004840354943,
      0.000004287618060, 0.000003841805181, 0.285283354587090, 0.014778182850780, 0.002293788755715, 0.000590907679355,
      0.000210980735153, 0.000095046686494, 0.000050833574321, 0.000030902119127, 0.000020679078749, 0.000014876620296,
      0.000011307265212, 0.000008965129659, 0.000007345720990, 0.000006177070348, 0.000005303385097, 0.000004630704703,
      0.000004099807190, 0.000003671920738, 0.284547300424615, 0.014664227676782, 0.002264505130685, 0.000580476399916,
      0.000206291304348, 0.000092543208458, 0.000049314010205, 0.000029886527712, 0.000019949196940, 0.000014322312006,
      0.000010868063560, 0.000008605455655, 0.000007043401289, 0.000005917605009, 0.000005076898979, 0.000004430235843,
      0.000003920299797, 0.000003509608157, 0.283816907632413, 0.014551723741984, 0.002235737259601, 0.000570276494065,
      0.000201725168960, 0.000090114454201, 0.000047844324742, 0.000028906746560, 0.000019246488957, 0.000013789520438,
      0.000010446473056, 0.000008260577972, 0.000006753775878, 0.000005669217088, 0.000004860215264, 0.000004238542509,
      0.000003748725848, 0.000003354527654, 0.283092104288427, 0.014440645545357, 0.002207473726359, 0.000560301889406,
      0.000197278654409, 0.000087757973992, 0.000046422760943, 0.000027961443512, 0.000018569900676, 0.000013277383268,
      0.000010041769523, 0.000007929875539, 0.000006476302957, 0.000005431427270, 0.000004652904942, 0.000004055236933,
      0.000003584731913, 0.000003206354780, 0.282372819731508, 0.014330968149162, 0.002179703413933, 0.000550546694246,
      0.000192948206432, 0.000085471404308, 0.000045047627126, 0.000027049338052, 0.000017918420191, 0.000012785073616,
      0.000009653259172, 0.000007612753786, 0.000006210464166, 0.000005203777215, 0.000004454557959, 0.000003879948619,
      0.000003427980414, 0.000003064779729, 0.281658984533620, 0.014222667164220, 0.002152415495485, 0.000541005191685,
      0.000188730386856, 0.000083252464629, 0.000043717294387, 0.000026169199241, 0.000017291076089, 0.000012311798558,
      0.000009280277301, 0.000007308643498, 0.000005955763555, 0.000004985828627, 0.000004264782368, 0.000003712323566,
      0.000003278148914, 0.000002929506672, 0.280950530472770, 0.014115718735620, 0.002125599425775, 0.000531671833916,
      0.000184621869522, 0.000081098954360, 0.000042430194172, 0.000025319843734, 0.000016686935782, 0.000011856797707,
      0.000008922187060, 0.000007016999714, 0.000005711726608, 0.000004777162375, 0.000004083203527, 0.000003552023535,
      0.000003134929429, 0.000002800253126, 0.280247390506643, 0.014010099528844, 0.002099244932838, 0.000522541236721,
      0.000180619436364, 0.000079008749876, 0.000041184815944, 0.000024500133875, 0.000016105103918, 0.000011419341842,
      0.000008578378262, 0.000006737300687, 0.000005477899302, 0.000004577377640, 0.000003909463325, 0.000003398725333,
      0.000002998027779, 0.000002676749347, 0.279549498746922, 0.013905786716312, 0.002073342009944, 0.000513608174162,
      0.000176719973636, 0.000076979801675, 0.000039979704943, 0.000023708975881, 0.000015544720853, 0.000010998731614,
      0.000008248266249, 0.000006469046874, 0.000005253847216, 0.000004386091112, 0.000003743219451, 0.000003252120150,
      0.000002867162964, 0.000002558737749, 0.278856790434270, 0.013802757964309, 0.002047880907807, 0.000504867573456,
      0.000172920468266, 0.000075010131642, 0.000038813460041, 0.000022945318085, 0.000015004961190, 0.000010594296288,
      0.000007931290804, 0.000006211759983, 0.000005039154673, 0.000004202936217, 0.000003584144685, 0.000003111912903,
      0.000002742066573, 0.000002445972361, 0.278169201913952, 0.013700991420309, 0.002022852127057, 0.000496314510030,
      0.000169218004361, 0.000073097830422, 0.000037684731673, 0.000022208149265, 0.000014485032379, 0.000010205392550,
      0.000007626915112, 0.000005964982053, 0.000004833423919, 0.000004027562379, 0.000003431926234, 0.000002977821632,
      0.000002622482212, 0.000002338218292, 0.277486670612087, 0.013600465700657, 0.001998246410945, 0.000487944202747,
      0.000165609759831, 0.000071241054892, 0.000036592219859, 0.000021496497034, 0.000013984173373, 0.000009831403366,
      0.000007334624768, 0.000005728274581, 0.000004636274347, 0.000003859634313, 0.000003286265086, 0.000002849576903,
      0.000002508164966, 0.000002235251235, 0.276809135012507, 0.013501159878617, 0.001974054738293, 0.000479752009298,
      0.000162093003134, 0.000069438025724, 0.000035534672298, 0.000020809426294, 0.000013501653344, 0.000009471736878,
      0.000007053926821, 0.000005501217678, 0.000004447341744, 0.000003698831353, 0.000003146875399, 0.000002726921255,
      0.000002398880880, 0.000002136856986, 0.276136534634197, 0.013403053472768, 0.001950268316656, 0.000471733421752,
      0.000158665090149, 0.000067687025056, 0.000034510882545, 0.000020146037762, 0.000013036770449, 0.000009125825361,
      0.000006784348868, 0.000005283409274, 0.000004266277578, 0.000003544846809, 0.000003013483916, 0.000002609608659,
      0.000002294406468, 0.000002042830984, 0.275468810009318, 0.013306126435727, 0.001926878575720, 0.000463884062262,
      0.000155323461162, 0.000065986394236, 0.000033519688252, 0.000019505466546, 0.000012588850648, 0.000008793124214,
      0.000006525438178, 0.000005074464346, 0.000004092748318, 0.000003397387352, 0.000002885829407, 0.000002497404010,
      0.000002194528239, 0.000001952977879, 0.274805902661778, 0.013210359143210, 0.001903877160895, 0.000456199678917,
      0.000152065637952, 0.000064334531662, 0.000032559969484, 0.000018886880783, 0.000012157246574, 0.000008473110997,
      0.000006276760863, 0.000004874014188, 0.000003926434775, 0.000003256172421, 0.000002763662133, 0.000002390082636,
      0.000002099042246, 0.000001867111111, 0.274147755086342, 0.013115732383399, 0.001881255927124, 0.000448676141738,
      0.000148889220998, 0.000062729890703, 0.000031630647102, 0.000018289480336, 0.000011741336447, 0.000008165284512,
      0.000006037901073, 0.000004681705705, 0.000003767031480, 0.000003120933670, 0.000002646743339, 0.000002287429834,
      0.000002007753662, 0.000001785052515, 0.273494310728272, 0.013022227346615, 0.001859006932881, 0.000441309438812,
      0.000145791886778, 0.000061170977698, 0.000030730681206, 0.000017712495539, 0.000011340523035, 0.000007869163917,
      0.000005808460241, 0.000004497200748, 0.000003614246089, 0.000002991414422, 0.000002534844765, 0.000002189240420,
      0.000001920476365, 0.000001706631938, 0.272845513963469, 0.012929825615289, 0.001837122434360, 0.000434095672556,
      0.000142771385167, 0.000059656350030, 0.000029859069643, 0.000017155185995, 0.000010954232657, 0.000007584287881,
      0.000005588056345, 0.000004320175468, 0.000003467798810, 0.000002867369165, 0.000002427748180, 0.000002095318310,
      0.000001837032547, 0.000001631686881, 0.272201310079113, 0.012838509154216, 0.001815594879855, 0.000427031056110,
      0.000139825536936, 0.000058184614278, 0.000029014846572, 0.000016616839420, 0.000010581914225, 0.000007310213778,
      0.000005376323210, 0.000004150319702, 0.000003327421859, 0.000002748563054, 0.000002325244940, 0.000002005476109,
      0.000001757252342, 0.000001560062149, 0.271561645254783, 0.012748260301089, 0.001794416904304, 0.000420111909855,
      0.000136952231333, 0.000056754424432, 0.000028197081087, 0.000016096770541, 0.000010223038334, 0.000007046516906,
      0.000005172909841, 0.000003987336388, 0.000003192858935, 0.000002634771445, 0.000002227135560, 0.000001919534725,
      0.000001680973468, 0.000001491609523, 0.270926466544045, 0.012659061757296, 0.001773581324017, 0.000413334658046,
      0.000134149423757, 0.000055364480182, 0.000027404875892, 0.000015594320032, 0.000009877096382, 0.000006792789746,
      0.000004977479779, 0.000003830941002, 0.000003063864726, 0.000002525779450, 0.000002133229310, 0.000001837322995,
      0.000001608040884, 0.000001426187443, 0.270295721856488, 0.012570896578985, 0.001753081131565, 0.000406695825558,
      0.000131415133516, 0.000054013525271, 0.000026637366031, 0.000015108853492, 0.000009543599724, 0.000006548641256,
      0.000004789710490, 0.000003680861022, 0.000002940204428, 0.000002421381506, 0.000002043343830, 0.000001758677337,
      0.000001538306471, 0.000001363660708, 0.269669359940206, 0.012483748168370, 0.001732909490819, 0.000400192034754,
      0.000128747441660, 0.000052700345904, 0.000025893717664, 0.000014639760473, 0.000009222078874, 0.000006313696178,
      0.000004609292777, 0.000003536835412, 0.000002821653290, 0.000002321380965, 0.000001957304752, 0.000001683441407,
      0.000001471628713, 0.000001303900191, 0.269047330364708, 0.012397600265291, 0.001713059732165, 0.000393820002439,
      0.000126144488903, 0.000051423769227, 0.000025173126895, 0.000014186453536, 0.000008912082726, 0.000006087594399,
      0.000004435930216, 0.000003398614131, 0.000002707996181, 0.000002225589707, 0.000001874945354, 0.000001611465778,
      0.000001407872407, 0.000001246782559, 0.268429583504242, 0.012312436939002, 0.001693525347845, 0.000387576536937,
      0.000123604473604, 0.000050182661855, 0.000024474818643, 0.000013748367356, 0.000008613177814, 0.000005869990314,
      0.000004269338622, 0.000003265957664, 0.000002599027168, 0.000002133827762, 0.000001796106219, 0.000001542607632,
      0.000001346908375, 0.000001192190018, 0.267816070521524, 0.012228242580193, 0.001674299987456, 0.000381458535250,
      0.000121125649838, 0.000048975928458, 0.000023798045559, 0.000013324957852, 0.000008324947602, 0.000005660552234,
      0.000004109245532, 0.000003138636569, 0.000002494549123, 0.000002045922958, 0.000001720634912, 0.000001476730463,
      0.000001288613198, 0.000001140010057, 0.267206743351862, 0.012145001893231, 0.001655377453586, 0.000375462980325,
      0.000118706325522, 0.000047802510405, 0.000023142086985, 0.000012915701364, 0.000008046991803, 0.000005458961812,
      0.000003955389711, 0.000003016431052, 0.000002394373336, 0.000001961710574, 0.000001648385668, 0.000001413703799,
      0.000001232868954, 0.000001090135213, 0.266601554687662, 0.012062699888617, 0.001636751697585, 0.000369586938402,
      0.000116344860615, 0.000046661384448, 0.000022506247950, 0.000012520093852, 0.000007778925726, 0.000005264913489,
      0.000003807520682, 0.000002899130547, 0.000002298319156, 0.000001881033016, 0.000001579219105, 0.000001353402932,
      0.000001179562972, 0.000001042462843, 0.266000457963302, 0.011981321875652, 0.001618416815457, 0.000363827556462,
      0.000114039665379, 0.000045551561468, 0.000021889858213, 0.000012137650133, 0.000007520379647, 0.000005078113972,
      0.000003665398273, 0.000002786533329, 0.000002206213637, 0.000001803739510, 0.000001513001934, 0.000001295708661,
      0.000001128587598, 0.000000996894906, 0.265403407340377, 0.011900853455304, 0.001600367043893, 0.000358182059748,
      0.000111789198706, 0.000044472085258, 0.000021292271334, 0.000011767903150, 0.000007270998214, 0.000004898281728,
      0.000003528792186, 0.000002678446135, 0.000002117891212, 0.000001729685794, 0.000001449606697, 0.000001240507046,
      0.000001079839970, 0.000000953337755, 0.264810357693287, 0.011821280513269, 0.001582596756417, 0.000352647749376,
      0.000109591966498, 0.000043422031358, 0.000020712863789, 0.000011410403267, 0.000007030439864, 0.000004725146501,
      0.000003397481583, 0.000002574683802, 0.000002033193365, 0.000001658733842, 0.000001388911503, 0.000001187689176,
      0.000001033221803, 0.000000911701938, 0.264221264595168, 0.011742589213223, 0.001565100459645, 0.000347222000025,
      0.000107446520112, 0.000042400505927, 0.000020151034113, 0.000011064717592, 0.000006798376277, 0.000004558448847,
      0.000003271254686, 0.000002475068922, 0.000001951968336, 0.000001590751588, 0.000001330799788, 0.000001137150942,
      0.000000988639180, 0.000000871902013, 0.263636084304161, 0.011664765990261, 0.001547872789669, 0.000341902257703,
      0.000105351454852, 0.000041406644664, 0.000019606202082, 0.000010730429327, 0.000006574491843, 0.000004397939693,
      0.000003149908403, 0.000002379431513, 0.000001874070821, 0.000001525612666, 0.000001275160076, 0.000001088792830,
      0.000000946002365, 0.000000833856363, 0.263054773749999, 0.011587797544510, 0.001530908508542, 0.000336686037581,
      0.000103305408522, 0.000040439611765, 0.000019077807920, 0.000010407137152, 0.000006358483156, 0.000004243379913,
      0.000003033247962, 0.000002287608703, 0.000001799361698, 0.000001463196160, 0.000001221885758, 0.000001042519707,
      0.000000905225608, 0.000000797487023, 0.262477290520914, 0.011511670834915, 0.001514202500873, 0.000331570921909,
      0.000101307060023, 0.000039498598919, 0.000018565311547, 0.000010094454617, 0.000006150058527, 0.000004094539917,
      0.000002921086561, 0.000002199444427, 0.000001727707762, 0.000001403386369, 0.000001170874874, 0.000000998240636,
      0.000000866226968, 0.000000762719521, 0.261903592850847, 0.011436373073198, 0.001497749770524, 0.000326554557992,
      0.000099355128006, 0.000038582824343, 0.000018068191841, 0.000009792009578, 0.000005948937514, 0.000003951199263,
      0.000002813245039, 0.000002114789136, 0.000001658981464, 0.000001346072576, 0.000001122029911, 0.000000955868682,
      0.000000828928146, 0.000000729482716, 0.261333639606953, 0.011361891717972, 0.001481545437405, 0.000321634656236,
      0.000097448369568, 0.000037691531848, 0.000017585945945, 0.000009499443638, 0.000005754850477, 0.000003813146283,
      0.000002709551556, 0.000002033499525, 0.000001593060673, 0.000001291148836, 0.000001075257607, 0.000000915320740,
      0.000000793254319, 0.000000697708649, 0.260767390277409, 0.011288214469018, 0.001465584734370, 0.000316808988257,
      0.000095585578991, 0.000036823989949, 0.000017118088586, 0.000009216411624, 0.000005567538150, 0.000003680177725,
      0.000002609841284, 0.000001955438262, 0.000001529828439, 0.000001238513760, 0.000001030468763, 0.000000876517360,
      0.000000759133985, 0.000000667332399, 0.260204804959486, 0.011215329261713, 0.001449863004198, 0.000312075385054,
      0.000093765586536, 0.000035979490995, 0.000016664151432, 0.000008942581072, 0.000005386751222, 0.000003552098407,
      0.000002513956120, 0.000001880473741, 0.000001469172772, 0.000001188070326, 0.000000987578063, 0.000000839382590,
      0.000000726498817, 0.000000638291949, 0.259645844347909, 0.011143224261611, 0.001434375696673, 0.000307431735236,
      0.000091987257266, 0.000035157350345, 0.000016223682465, 0.000008677631746, 0.000005212249950, 0.000003428720892,
      0.000002421744403, 0.000001808479832, 0.000001410986428, 0.000001139725680, 0.000000946503906, 0.000000803843816,
      0.000000695283516, 0.000000610528054, 0.259090469723486, 0.011071887859159, 0.001419118365740, 0.000302875983314,
      0.000090249489916, 0.000034356905559, 0.000015796245384, 0.000008421255164, 0.000005043803773, 0.000003309865170,
      0.000002333060644, 0.000001739335657, 0.000001355166707, 0.000001093390962, 0.000000907168240, 0.000000769831616,
      0.000000665425680, 0.000000583984115, 0.258538642941982, 0.011001308664560, 0.001404086666760, 0.000298406128041,
      0.000088551215807, 0.000033577515636, 0.000015381419028, 0.000008173154151, 0.000004881190955, 0.000003195358357,
      0.000002247765275, 0.000001672925360, 0.000001301615253, 0.000001048981127, 0.000000869496409, 0.000000737279619,
      0.000000636865674, 0.000000558606061, 0.257990326423266, 0.010931475502769, 0.001389276353830, 0.000294020220815,
      0.000086891397788, 0.000032818560262, 0.000014978796820, 0.000007933042405, 0.000004724198227, 0.000003085034404,
      0.000002165724397, 0.000001609137899, 0.000001250237871, 0.000001006414783, 0.000000833417001, 0.000000706124369,
      0.000000609546505, 0.000000534342235, 0.257445483140691, 0.010862377408624, 0.001374683277197, 0.000289716364121,
      0.000085269029221, 0.000032079439101, 0.000014587986236, 0.000007700644083, 0.000004572620461, 0.000002978733820,
      0.000002086809548, 0.000001547866844, 0.000001200944350, 0.000000965614027, 0.000000798861706, 0.000000676305197,
      0.000000583413704, 0.000000511143288, 0.256904076610722, 0.010794003622097, 0.001360303380735, 0.000285492710036,
      0.000083683133002, 0.000031359571101, 0.000014208608290, 0.000007475693402, 0.000004426260345, 0.000002876303406,
      0.000002010897475, 0.000001489010175, 0.000001153648287, 0.000000926504301, 0.000000765765181, 0.000000647764095,
      0.000000558415215, 0.000000488962071, 0.256366070882789, 0.010726343583679, 0.001346132699505, 0.000281347458773,
      0.000082132760609, 0.000030658393834, 0.000013840297042, 0.000007257934255, 0.000004284928074, 0.000002777596000,
      0.000001937869922, 0.000001432470105, 0.000001108266927, 0.000000889014239, 0.000000734064921, 0.000000620445600,
      0.000000534501287, 0.000000467753540, 0.255831430529383, 0.010659386929876, 0.001332167357386, 0.000277278857274,
      0.000080616991190, 0.000029975362848, 0.000013482699124, 0.000007047119844, 0.000004148441055, 0.000002682470232,
      0.000001867613419, 0.000001378152894, 0.000001064721008, 0.000000853075536, 0.000000703701130, 0.000000594296684,
      0.000000511624368, 0.000000447474658 };

  if (x <= 0.000000000000001)
    type = 0;
  else if (x >= 30.000000000000000)
    type = 3;
  else if (x > 0.000000000000001 && x < 12.000000000000000)
    type = 1;
  else
    type = 2;

  switch (type) {
    case 0:
      for (i2 = 0; i2 <= n; i2++) {
        f[i2] *= G[i2];
        f[i2] += k_one;
      }
      break;

    case 1:

      T = -0.05 * modf(x * 20.0, &i0); // 0.05 is stepsize  20.0 its reciprocal
      i1 = 18 * int(i0); // 18 is number of tabulated Fn

      for (i2 = 0; i2 <= n; i2++) {
        f[i2] *= ((((((G[i1 + i2 + 6] * T * sixth + G[i1 + i2 + 5]) * T * fifth + G[i1 + i2 + 4]) * T * fourth + G[i1
            + i2 + 3]) * T * third + G[i1 + i2 + 2]) * T * half + G[i1 + i2 + 1]) * T + G[i1 + i2]);
        f[i2] += k_one;
      }
      break;

    case 2:
      double g1, x1, expx;
      x1 = k_one / x;
      expx = exp(-x);
      if (x >= 24.000000000000000 && x < 30.000000000000000) {
        g1 = 0.490;
      } else if (x >= 18.000000000000000 && x < 24.000000000000000) {
        g1 = 0.499093162 - 0.2152832 * x1;
      } else if (x >= 15.000000000000000 && x < 18.000000000000000) {
        g1 = 0.4998436875 - 0.24249438 * x1 + 0.24642845 * x1 * x1;
      } else {
        g1 = (((-1.186732847845 * x1 + k_one) * -1.298418478457 * x1 + k_one) * -0.2473631686 * x1 + 0.4999489092);
      }
      g[0] = sqrt(fourth * pi * x1) - expx * x1 * g1;
      f[0] *= g[0];
      f[0] += k_one;
      for (j = 0; j < n; j++) {
        g[j + 1] = ((2 * j + 1) * g[j] - expx) * half * x1;
        f[j + 1] *= g[j + 1];
        f[j + 1] += k_one;
      }
      break;

    case 3:
      x1 = k_one / x;
      expx = exp(-x);
      g[0] = sqrt(fourth * pi * x1);
      f[0] *= g[0];
      f[0] += k_one;
      if (x < 36.000000000000000) {
        for (j = 0; j < n; j++) {
          g[j + 1] = ((2 * j + 1) * g[j] - expx) * half * x1;
          f[j + 1] *= g[j + 1];
          f[j + 1] += k_one;
        }
      } else {
        for (j = 0; j < n; j++) {
          g[j + 1] = (2 * j + 1) * g[j] * half * x1;
          f[j + 1] *= g[j + 1];
          f[j + 1] += k_one;
        }
      }
      break;

  }

}

