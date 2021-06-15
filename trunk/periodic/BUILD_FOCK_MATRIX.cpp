/*
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include "INTEGRALS_TWO_CENTRE.h"
*/

#include <cstring>
#include "mycomplex.h"
#include "myconstants.h"
#include "conversion_factors.h"
#include <xc.h>
#include "USER_DATA.h"
#include "DFT.h"
#include "LIMITS.h"
#include "PAIRS_QUADS.h"
#include "PRINT_UTIL.h"
#include "MATRIX_UTIL.h"
#include "PARALLEL.h"
#include "ROTATIONS_MOLECULE.h"
#include "ALLOCATE_MEMORY.h"
#include "DENSITY_MATRIX.h"
#include "INTEGRALS_4C_CRYSTAL.h"
#include "INTEGRALS_4C_MOLECULE.h"
#include "BUILD_FOCK_MATRIX.h"

using namespace std;

  // ******************************************************************************************
  // * This routine builds the Fock Matrix                                                    *
  // ******************************************************************************************
/*
void build_fock_matrix_molecule(double *Fock, INT_1E *one_ints, FERMI *fermi, double *total_electrons, double *S1, double *P, double *F, double *delta_P, double *delta_F, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry,REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k, l, m, p, q, s, count;
int dimp, dimf, dimp_spin, dimf_spin, dimm_spin;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
double time1, time2;
double *Fock_2c, *Fock_2e, *Kohn_2e, *F_up_down;

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;

  AllocateDoubleArray(&Fock_2c,&dimp_spin,job);
  AllocateDoubleArray(&Fock_2e,&dimp_spin,job);
  AllocateDoubleArray(&Kohn_2e,&dimp_spin,job);
  ResetDoubleArray(Fock_2c, &dimp_spin);
  ResetDoubleArray(Fock_2e, &dimp_spin);
  ResetDoubleArray(Kohn_2e, &dimp_spin);

  PAIR_TRAN pair_q;
  //count_density_pairs4(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  //allocate_PAIR_TRAN(&pair_q,atoms,symmetry,R_tables,job,file);
  //generate_density_pairs4(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  //pair_q.cutoff = (job->itol1 > job->itol2) ? job->itol1 : job->itol2;
  pair_q.cutoff = R->cutoff;
  count_range_selected_pairs(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_q,atoms,symmetry,R_tables,job,file);
  generate_range_selected_pairs(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->iter == 1) print_pairs(&pair_q,atoms,R,job,file);
  */

  /*
  switch (crystal->type[0]) {

  case 'C':
  case 'S':
  case 'P':

  if (job->scf_coulomb == 1) {

  if (job->scf_direct == 0) {

    AllocateDoubleArray(&F_up_down,&job->dimf,job);
    if (job->spin_dim == 1) {
      for (m = 0; m < job->dimf; m++) F_up_down[m] = F[m];
     }
    else if (job->spin_dim == 2) {
      ResetDoubleArray(F_up_down,&job->dimf);
      for (s = 0; s < job->spin_dim; s++) {
        for (m = 0; m < job->dimf; m++) F_up_down[m] += F[s * job->dimf + m];
        }
       }

  if (job->int_exist == 0)

    coulomb_matrix_crystal_compute_integrals(Fock_2c,S1,F_up_down,pair_p,&pair_q,atom_p,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,G,job,file);

  else

    coulomb_matrix_crystal_read_integrals(Fock_2c,F_up_down,pair_p,atoms,shells,symmetry,R_tables,job,file);

  }

  else if (job->scf_direct == 1) {

    AllocateDoubleArray(&F_up_down,&job->dimf,job);
    if (job->spin_dim == 1) {
      for (m = 0; m < job->dimf; m++) F_up_down[m] = delta_F[m];
     }
    else if (job->spin_dim == 2) {
      ResetDoubleArray(F_up_down,&job->dimf);
      for (s = 0; s < job->spin_dim; s++) {
        for (m = 0; m < job->dimf; m++) F_up_down[m] += delta_F[s * job->dimf + m];
        }
       }

  if (job->int_exist == 0)

    coulomb_matrix_crystal_compute_integrals(Fock_2c,S1,F_up_down,pair_p,&pair_q,atom_p,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,G,job,file);

  else

    coulomb_matrix_crystal_read_integrals(Fock_2c,F_up_down,pair_p,atoms,shells,symmetry,R_tables,job,file);

  }

    DestroyDoubleArray(&F_up_down,&job->dimf,job);

 } // close if (job->scf_coulomb == 1)

 
  if (job->scf_exchange == 1) {

  if (job->scf_direct == 0) {

    exchange_matrix_crystal_compute_integrals(Fock_2e,S2,F,pair_p,&pair_q,atom_p,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,G,job,file);
 
  } 

  else if (job->scf_direct == 1) {

    exchange_matrix_crystal_compute_integrals(Fock_2e,S2,delta_F,pair_p,&pair_q,atom_p,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,G,job,file);

  } 

 } // close if (job->scf_exchange == 1)
 
  break;

  case 'M':
  */
/*
  if (job->scf_coulomb == 1 || job->scf_exchange == 1) {

  if (job->scf_direct == 0) {

  if (job->int_exist == 0)

    fock_matrix_molecule_compute_integrals(Fock_2c,Fock_2e,S1,F,pair_p,&pair_q,atom_p,atoms,shells,gaussians,\
    crystal,symmetry,R,R_tables,job,file);

  else

     fock_matrix_molecule_read_integrals(Fock_2c,Fock_2e,F,pair_p,atoms,shells,symmetry,job,file);

  }

  else if (job->scf_direct == 1) {

     fock_matrix_molecule_compute_integrals(Fock_2c,Fock_2e,S1,delta_F,pair_p,&pair_q,atom_p,atoms,shells,gaussians,\
     crystal,symmetry,R,R_tables,job,file);
 
  }

 } // close if (job->scf_coulomb == 1 || job->scf_exchange == 1)
 */

  /*
  break;

 }
 */

  ////double *Fock_EN;
  ////AllocateDoubleArray(&Fock_EN,&dimp_spin,job);
  ////ResetDoubleArray(Fock_EN, &dimp_spin);
  //fock_element_elecnuc(Fock_2c,&pair_q,R,G,atoms,shells,gaussians,crystal,job,file);
  ////fock_element_elecnuc(Fock_EN,pair_p,R,G,atoms,shells,gaussians,crystal,job,file);
  //DestroyDoubleArray(&Fock_EN,&dimp_spin,job);

  //for (i = 0; i < dimp; i++) { //fprintf(file.out,"EN %3d %10.4f %10.4f\n",i,Fock_EN[i],one_ints->ElecNuc[i]);
  //if (fabs(Fock_EN[i]-one_ints->ElecNuc[i]) > k_zero) fprintf(file.out,"diff "); }
  /*  
  if (job->scf_direct == 0) {
  ResetDoubleArray(Fock, &dimp_spin);
  count = 0;
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < dimp; i++) {
      Fock[count] += one_ints->Fock[i] + Fock_2c[count] + (double)job->spin_dim * Fock_2e[count] + two * Kohn_2e[count];
      count++;
     }
    }
   }

  if (job->scf_direct == 1) {
  if (job->iter == 1) {
  ResetDoubleArray(Fock, &dimp_spin);
  count = 0;
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < dimp; i++) {
      Fock[count] += one_ints->Fock[i];
      count++;
     }
    }
   }
  count = 0;
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < dimp; i++) {
      Fock[count] += Fock_2c[count] + (double) job->spin_dim * Fock_2e[count] + two * Kohn_2e[count];
      count++;
     }
    }
   }

  if (job->taskid == 0 && job->scf_direct == 0)
  total_energy_final(one_ints, Fock_2c, Fock_2e, Kohn_2e, P, pair_p, atoms,job,file);
  else if (job->taskid == 0 && job->scf_direct == 1)
  total_energy_direct(one_ints, Fock, P, pair_p, atoms,job,file);

  DestroyDoubleArray(&Fock_2c,&dimp_spin,job);
  DestroyDoubleArray(&Fock_2e,&dimp_spin,job);
  DestroyDoubleArray(&Kohn_2e,&dimp_spin,job);

}
*/

void build_fock_matrix(double *Fock, INT_1E *one_ints, FERMI *fermi, double *total_electrons, double *S1, double *S2, double *P, double *F, double *delta_P, double *delta_F, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry,REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k, l, m, p, q, s, count;
int dimp, dimf, dimp_spin, dimf_spin, dimm_spin;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
double time1, time2;
double *Fock_2c, *Fock_2e, *Kohn_2e, *F_up_down;

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;

  AllocateDoubleArray(&Fock_2c,&dimp_spin,job);
  AllocateDoubleArray(&Fock_2e,&dimp_spin,job);
  AllocateDoubleArray(&Kohn_2e,&dimp_spin,job);
  ResetDoubleArray(Fock_2c, &dimp_spin);
  ResetDoubleArray(Fock_2e, &dimp_spin);
  ResetDoubleArray(Kohn_2e, &dimp_spin);

  PAIR_TRAN pair_q;
  //count_density_pairs4(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  //allocate_PAIR_TRAN(&pair_q,atoms,symmetry,R_tables,job,file);
  //generate_density_pairs4(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  //pair_q.cutoff = (job->itol1 > job->itol2) ? job->itol1 : job->itol2;
  pair_q.cutoff = R->cutoff;
  count_range_selected_pairs(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_q,atoms,symmetry,R_tables,job,file);
  generate_range_selected_pairs(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->iter == 1) print_pairs(&pair_q,atoms,R,job,file);

  switch (crystal->type[0]) {

  case 'C':
  case 'S':
  case 'P':

  if (job->scf_coulomb == 1) {

  if (job->scf_direct == 0) {

    AllocateDoubleArray(&F_up_down,&job->dimf,job);
    if (job->spin_dim == 1) {
      for (m = 0; m < job->dimf; m++) F_up_down[m] = F[m];
     }
    else if (job->spin_dim == 2) {
      ResetDoubleArray(F_up_down,&job->dimf);
      for (s = 0; s < job->spin_dim; s++) {
        for (m = 0; m < job->dimf; m++) F_up_down[m] += F[s * job->dimf + m];
        }
       }

  if (job->int_exist == 0)

    fock_matrix_crystal_compute_coulomb_integrals(Fock_2c,S1,F_up_down,pair_p,&pair_q,atom_p,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,G,job,file);

  else

    fock_matrix_crystal_read_coulomb_integrals(Fock_2c,F_up_down,pair_p,atoms,shells,symmetry,R_tables,job,file);

  }

  else if (job->scf_direct == 1) {

    AllocateDoubleArray(&F_up_down,&job->dimf,job);
    if (job->spin_dim == 1) {
      for (m = 0; m < job->dimf; m++) F_up_down[m] = delta_F[m];
     }
    else if (job->spin_dim == 2) {
      ResetDoubleArray(F_up_down,&job->dimf);
      for (s = 0; s < job->spin_dim; s++) {
        for (m = 0; m < job->dimf; m++) F_up_down[m] += delta_F[s * job->dimf + m];
        }
       }
  if (job->int_exist == 0)

    fock_matrix_crystal_compute_coulomb_integrals(Fock_2c,S1,F_up_down,pair_p,&pair_q,atom_p,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,G,job,file);

  else

    fock_matrix_crystal_read_coulomb_integrals(Fock_2c,F_up_down,pair_p,atoms,shells,symmetry,R_tables,job,file);

  }

    DestroyDoubleArray(&F_up_down,&job->dimf,job);

 } // close if (job->scf_coulomb == 1)

 
  if (job->scf_exchange == 1) {

  if (job->scf_direct == 0) {

    fock_matrix_crystal_compute_exchange_integrals(Fock_2e,S2,F,pair_p,&pair_q,atom_p,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,G,job,file);
 
  } 

  else if (job->scf_direct == 1) {

    fock_matrix_crystal_compute_exchange_integrals(Fock_2e,S2,delta_F,pair_p,&pair_q,atom_p,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,G,job,file);

  } 

 } // close if (job->scf_exchange == 1)
 
  break;

  case 'M':

  if (job->scf_coulomb == 1 || job->scf_exchange == 1) {

  if (job->scf_direct == 0) {

  if (job->int_exist == 0)

    fock_matrix_molecule_compute_integrals(Fock_2c,Fock_2e,S1,F,pair_p,&pair_q,atom_p,atoms,shells,gaussians,\
    crystal,symmetry,R,R_tables,job,file);

  else

     fock_matrix_molecule_read_integrals(Fock_2c,Fock_2e,F,pair_p,atoms,shells,symmetry,job,file);

  }

  else if (job->scf_direct == 1) {

     fock_matrix_molecule_compute_integrals(Fock_2c,Fock_2e,S1,delta_F,pair_p,&pair_q,atom_p,atoms,shells,gaussians,\
     crystal,symmetry,R,R_tables,job,file);
 
  }

 } // close if (job->scf_coulomb == 1 || job->scf_exchange == 1)

  break;

 }

  ////double *Fock_EN;
  ////AllocateDoubleArray(&Fock_EN,&dimp_spin,job);
  ////ResetDoubleArray(Fock_EN, &dimp_spin);
  //fock_element_elecnuc(Fock_2c,&pair_q,R,G,atoms,shells,gaussians,crystal,job,file);
  ////fock_element_elecnuc(Fock_EN,pair_p,R,G,atoms,shells,gaussians,crystal,job,file);
  //DestroyDoubleArray(&Fock_EN,&dimp_spin,job);

  //for (i = 0; i < dimp; i++) { //fprintf(file.out,"EN %3d %10.4f %10.4f\n",i,Fock_EN[i],one_ints->ElecNuc[i]);
  //if (fabs(Fock_EN[i]-one_ints->ElecNuc[i]) > k_zero) fprintf(file.out,"diff "); }
    
  if (job->scf_direct == 0) {
  ResetDoubleArray(Fock, &dimp_spin);
  count = 0;
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < dimp; i++) {
      Fock[count] += one_ints->Fock[i] + Fock_2c[count] + (double)job->spin_dim * Fock_2e[count] + two * Kohn_2e[count];
      count++;
     }
    }
   }

  if (job->scf_direct == 1) {
  if (job->iter == 1) {
  ResetDoubleArray(Fock, &dimp_spin);
  count = 0;
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < dimp; i++) {
      Fock[count] += one_ints->Fock[i];
      count++;
     }
    }
   }
  count = 0;
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < dimp; i++) {
      Fock[count] += Fock_2c[count] + (double) job->spin_dim * Fock_2e[count] + two * Kohn_2e[count];
      count++;
     }
    }
   }

  if (job->taskid == 0 && job->scf_direct == 0)
  total_energy_final(one_ints, Fock_2c, Fock_2e, Kohn_2e, P, pair_p, atoms,job,file);
  else if (job->taskid == 0 && job->scf_direct == 1)
  total_energy_direct(one_ints, Fock, P, pair_p, atoms,job,file);

  DestroyDoubleArray(&Fock_2c,&dimp_spin,job);
  DestroyDoubleArray(&Fock_2e,&dimp_spin,job);
  DestroyDoubleArray(&Kohn_2e,&dimp_spin,job);

}

void build_fock_matrix_no_sym(double *Fock, DFT_GRID *dft_grid, INT_1E *one_ints, FERMI *fermi, double *total_electrons, double *S1, double *S2, double *P, double *F, double *delta_P, double *delta_F, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry,REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k, l, p, q, s, count;
int dimp, dimf, dimp_spin, dimf_spin, dimm_spin;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
double time1, time2;
double *Fock_2c, *Fock_2e, *Kohn_2e;

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;

  AllocateDoubleArray(&Fock_2c,&dimp_spin,job);
  AllocateDoubleArray(&Fock_2e,&dimp_spin,job);
  AllocateDoubleArray(&Kohn_2e,&dimp_spin,job);
  ResetDoubleArray(Fock_2c, &dimp_spin);
  ResetDoubleArray(Fock_2e, &dimp_spin);
  ResetDoubleArray(Kohn_2e, &dimp_spin);

  PAIR_TRAN pair_q;
  pair_q.cutoff = (job->itol1 > job->itol2) ? job->itol1 : job->itol2;
  pair_q.cutoff = R->cutoff;
  count_range_selected_pairs(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_q,atoms,symmetry,R_tables,job,file);
  generate_range_selected_pairs(&pair_q,atoms,atom_p,symmetry,R,R_tables,job,file);
  //if (job->iter == 1) print_pairs(&pair_q,atoms,R,job,file);

  switch (crystal->type[0]) {

  case 'C':
  case 'S':
  case 'P':

//job->int_exist_no_sym = 1;

  if (job->scf_coulomb == 1) {

  if (job->int_exist_no_sym == 0 || job->scf_direct == 1)

    fock_matrix_crystal_compute_coulomb_integrals_no_sym(Fock_2c,S1,F,pair_p,&pair_q,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,G,job,file);

  else

    fock_matrix_crystal_read_coulomb_integrals_no_sym(Fock_2c,F,pair_p,atoms,job,file);

 } // close if (job->scf_coulomb == 1)

 
  if (job->scf_exchange == 1) {

    fock_matrix_crystal_compute_exchange_integrals_no_sym(Fock_2e,S2,F,pair_p,&pair_q,atoms,shells,gaussians,crystal,symmetry,R,\
    R_tables,job,file);

 } // close if (job->scf_exchange == 1)
 
  break;

  case 'M':

  if (job->scf_coulomb == 1 || job->scf_exchange == 1) {

  if (job->int_exist_no_sym == 0 || job->scf_direct == 1) 

   fock_matrix_molecule_compute_integrals_no_sym(Fock_2c,Fock_2e,S2,F,pair_p,&pair_q,atoms,shells,gaussians,\
   crystal,symmetry,R,job,file);

  else

   fock_matrix_molecule_read_integrals_no_sym(Fock_2c,Fock_2e,F,pair_p,atoms,shells,symmetry,job,file);

 } // close if (job->scf_coulomb == 1 || job->scf_exchange == 1)

  break;

 }

  count = 0;
  ResetDoubleArray(Fock, &dimp_spin);
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < dimp; i++) {
      //16/09/20 Fock[count] += one_ints->Fock[count] + Fock_2c[count] + Fock_2e[count] + two * Kohn_2e[count];
      Fock[count] += one_ints->Fock[i] + Fock_2c[count] + (double)job->spin_dim * Fock_2e[count] + two * Kohn_2e[count];
      count++;
     }
    }

  //double total_energy = job->total_energy;
  if (job->taskid == 0)
  total_energy_final(one_ints, Fock_2c, Fock_2e, Kohn_2e, P, pair_p, atoms,job,file);
  //job->total_energy = total_energy;

  //print_Fock_matrix(Fock, pair_p, atoms, job, file);

  DestroyDoubleArray(&Fock_2c,&dimp_spin,job);
  DestroyDoubleArray(&Fock_2e,&dimp_spin,job);
  DestroyDoubleArray(&Kohn_2e,&dimp_spin,job);

}

/*
void total_energy(INT_1E *one_ints, double *Fock, double *P, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int j, p, q, s, count, count1;
double current_total_energy;
double total_energy = k_zero;

    current_total_energy = job->total_energy;

    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      count1 = 0;
      for (p = 0; p < pair_p->nump; p++) {
	q = pair_p->posn[p];
	for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
	  //total_energy += (one_ints->Kinetic[count1] + one_ints->ElecNuc[count1]) * P[count] * pair_p->numb[p];
	  total_energy += (Fock[count1] + one_ints->Kinetic[count1] + one_ints->ElecNuc[count1]) * P[count] * pair_p->numb[p];
	  count++;
	  count1++;
	 }
	}
       }

//May2013
    ////job->energy_change = total_energy / two + job->nuc_nuc - current_total_energy;
    ////job->total_energy  = total_energy / two + job->nuc_nuc;
//May2013
    //job->total_energy  = total_energy ;

    //if (job->taskid == 0) 
    //fprintf(file.out, "Total Energy %17.9e Total Energy Change %17.9e\n\n", job->total_energy,job->energy_change);
    //if (job->taskid == 0) 

}

void total_energy_final(INT_1E *one_ints, double *F_2c, double *F_2e, double *K_2e, double *P, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int i, j, p, q, s, count, count1;
double w;
double current_total_energy;
double oneelec_energy  = k_zero;
double twoelec_energy  = k_zero;
double kinetic_energy  = k_zero;
double elecnuc_energy  = k_zero;
double fock_energy     = k_zero;
double coulomb_energy  = k_zero;
double exchange_energy = k_zero;
double exc_corr_energy = k_zero;

  w = k_one;
 // if (job->xc_typ[0] > 0) {
 //   xc_func_type func[job->xc_num];
 //   w = k_zero;
 //   for (i = 0; i < job->xc_num; i++) {
 //     xc_func_init(&func[i], job->xc_typ[i], job->spin_dim);
 //       if (func[i].info->family == 32)
 //         if (job->taskid == 0) fprintf(file.out,"Fix call to xc_func at line 1190 in SCF.cpp\n");
 //         //w = xc_hyb_gga_exx_coef(func[i].gga);
 //        //fprintf(file.out,"weight %lf\n",w);
 //       }
 //      xc_func_end(&func[0]);
 //     }

    current_total_energy = job->total_energy;

    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      count1 = 0;
      for (p = 0; p < pair_p->nump; p++) {
	q = pair_p->posn[p];
	for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
	  kinetic_energy += one_ints->Kinetic[count1] * P[count] * pair_p->numb[p];
	  count++;
	  count1++;
	 }
	}
       }

    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      count1 = 0;
      for (p = 0; p < pair_p->nump; p++) {
	q = pair_p->posn[p];
	for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
	  elecnuc_energy += one_ints->ElecNuc[count1] * P[count] * pair_p->numb[p];
	  count++;
	  count1++;
	 }
	}
       }

    oneelec_energy = elecnuc_energy + kinetic_energy;

    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      for (p = 0; p < pair_p->nump; p++) {
	q = pair_p->posn[p];
	for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
	  coulomb_energy += F_2c[count] * P[count] * pair_p->numb[p];
	  //fprintf(file.out,"coulomb %3d %3d %3d %3d %12.5e %12.5e %12.5e %3d\n",\
          s,p,q,count,coulomb_energy,F_2c[count] , P[count] , pair_p->numb[p]);
	  count++;
	 }
	}
       }

    if (job->xc_hfx == 1) {
    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      for (p = 0; p < pair_p->nump; p++) {
	q = pair_p->posn[p];
	for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
	  exchange_energy += (double) job->spin_dim * F_2e[count] * P[count] * pair_p->numb[p] * w;
	  //fprintf(file.out,"exchange %4d %4d %4d %4d %12.5e %12.5e %12.5e %4d\n",\
          s, p, q, count, exchange_energy, F_2e[count], P[count], pair_p->numb[p]);
	  count++;
	 }
	}
       }
      }

    if (job->xc_num > 0) {
    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      for (p = 0; p < pair_p->nump; p++) {
	q = pair_p->posn[p];
	for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
	  exc_corr_energy += K_2e[count] * P[count] * pair_p->numb[p];
	  count++;
	 }
	}
       }
      }

    twoelec_energy = coulomb_energy + exchange_energy + two * exc_corr_energy;

    job->energy_change = (oneelec_energy + twoelec_energy / two) + job->nuc_nuc - current_total_energy;
    job->total_energy  = (oneelec_energy + twoelec_energy / two) + job->nuc_nuc;

    if (job->taskid == 0) {
    printf(           "kinetic energy               %17.9e\n", kinetic_energy);
    printf(           "elecnuc energy               %17.9e\n", elecnuc_energy);
    printf(           "coulomb energy               %17.9e\n", coulomb_energy / two);
    if (job->xc_num  > 0) 
    printf(           "exchange correlation energy  %17.9e\n", exc_corr_energy);
    if (job->xc_hfx == 1) 
    printf(           "exchange energy              %17.9e\n", exchange_energy / two);
    printf(           "nuclear repulsion energy     %17.9e\n", job->nuc_nuc);
    printf(           "total energy                 %17.9e\n", job->total_energy);
    printf(           "total energy change          %17.9e\n", job->energy_change);
   }

    if (job->taskid == 0 && fabs(job->energy_change) < job->scf_tol) {

    fprintf(file.out,  "kinetic energy               %17.9e\n", kinetic_energy);
    fprintf(file.out,  "elecnuc energy               %17.9e\n", elecnuc_energy);
    fprintf(file.out,  "coulomb energy               %17.9e\n", coulomb_energy / two);
    if (job->xc_num  > 0) 
    fprintf(file.out,  "exchange correlation energy  %17.9e\n", exc_corr_energy);
    if (job->xc_hfx == 1) 
    fprintf(file.out,  "exchange energy              %17.9e\n", exchange_energy / two);
    fprintf(file.out,  "nuclear repulsion energy     %17.9e\n", job->nuc_nuc);
    fprintf(file.out,  "total energy                 %17.9e\n", job->total_energy);
    fprintf(file.out,  "total energy change          %17.9e\n", job->energy_change);
   }

}

void total_energy_direct(INT_1E *one_ints, double *Fock, double *P, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int i, j, p, q, s, count, count1;
double w;
double current_total_energy;
double oneelec_energy  = k_zero;
double twoelec_energy  = k_zero;
double kinetic_energy  = k_zero;
double elecnuc_energy  = k_zero;

  job->twoe_energy = k_zero;

//  w = k_one;
//  if (job->xc_typ[0] > 0) {
//    xc_func_type func[job->xc_num];
//    w = k_zero;
//    for (i = 0; i < job->xc_num; i++) {
//      xc_func_init(&func[i], job->xc_typ[i], job->spin_dim);
//	if (func[i].info->family == 32)
//	  if (job->taskid == 0) fprintf(file.out,"Fix call to xc_func at line 1190 in SCF.cpp\n");
//	  //w = xc_hyb_gga_exx_coef(func[i].gga);
//	 //fprintf(file.out,"weight %lf\n",w);
//	}
//       xc_func_end(&func[0]);
//      }

    current_total_energy = job->total_energy;

    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      count1 = 0;
      for (p = 0; p < pair_p->nump; p++) {
	q = pair_p->posn[p];
	for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
	  kinetic_energy += one_ints->Kinetic[count1] * P[count] * pair_p->numb[p];
	  count++;
	  count1++;
	 }
	}
       }

    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      count1 = 0;
      for (p = 0; p < pair_p->nump; p++) {
	q = pair_p->posn[p];
	for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
	  elecnuc_energy += one_ints->ElecNuc[count1] * P[count] * pair_p->numb[p];
	  count++;
	  count1++;
	 }
	}
       }

    oneelec_energy = elecnuc_energy + kinetic_energy;

    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
      count1 = 0;
      for (p = 0; p < pair_p->nump; p++) {
	q = pair_p->posn[p];
	for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
	  job->twoe_energy += (Fock[count] - one_ints->Fock[count1]) * P[count] * pair_p->numb[p];
	  count++;
	  count1++;
	 }
	}
       }

    //twoelec_energy = coulomb_energy + exchange_energy + two * exc_corr_energy;
    twoelec_energy = job->twoe_energy;

    job->energy_change = (oneelec_energy + twoelec_energy / two) + job->nuc_nuc - current_total_energy;
    job->total_energy  = (oneelec_energy + twoelec_energy / two) + job->nuc_nuc;

    if (job->taskid == 0) {
    printf(           "\n");
    printf(           "kinetic energy               %17.9e\n", kinetic_energy);
    printf(           "elecnuc energy               %17.9e\n", elecnuc_energy);
    printf(           "twoelec energy               %17.9e\n", job->twoe_energy / two);
    printf(           "nuclear repulsion energy     %17.9e\n", job->nuc_nuc);
    printf(           "total energy                 %17.9e\n", job->total_energy);
    printf(           "total energy change          %17.9e\n", job->energy_change);
   }

    if (job->taskid == 0 && fabs(job->energy_change) < job->scf_tol) {

    fprintf(file.out,  "kinetic energy               %17.9e\n", kinetic_energy);
    fprintf(file.out,  "elecnuc energy               %17.9e\n", elecnuc_energy);
    fprintf(file.out,  "twoelec energy               %17.9e\n", job->twoe_energy / two);
    fprintf(file.out,  "nuclear repulsion energy     %17.9e\n", job->nuc_nuc);
    fprintf(file.out,  "total energy                 %17.9e\n", job->total_energy);
    fprintf(file.out,  "total energy change          %17.9e\n", job->energy_change);
   }

}
*/

void fock_matrix_crystal_compute_coulomb_integrals(double *Fock_2c, double *S1, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i1;
int dimp, dimf, dimp_spin, dimf_spin;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
int dimi, nshells;
int atm_n1, atm_n2, atm_n3, atm_n4;
int lat_n1, lat_n2, lat_n3, lat_n4;
int nd1, nd2, nd3, nd4;
int all_quads, unique_quads, total_quads, largest;
int winRank;
long i, j, q;
long total_integrals;
long begin_task, end_task, total_tasks;
long *counter, myCounter, increment;
const long izero = 0;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, time11, time12;
char xx[4], yy[24] = "scf_integrals_2C.";
double *Fock_2c_buffer;
FILE *integrals_2c; 
MPI_Win win;

  sprintf(xx, "%d", job->taskid);

  //if (job->scf_direct == 0 && (job->taskid > 0 || job->numtasks == 1)) integrals_2c = fopen(strcat(yy,xx), "wb");
  if (job->taskid > 0 || job->numtasks == 1) integrals_2c = fopen(strcat(yy,xx), "wb");

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;
  AllocateDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2c_buffer,&dimp_spin);

  total_tasks = (long) pair_q->tot * (long) pair_q->tot;

  time1 = k_zero;
  time3 = k_zero;
  time5 = k_zero;
  time7 = k_zero;
  time9 = k_zero;
  time11 = k_zero;
  all_quads = 0;
  unique_quads = 0;
  total_quads = 0;
  total_integrals = 0;

  winRank = 0;
  CreateLongCounter(job->taskid, winRank, &counter, izero, &win, MPI_COMM_WORLD);
  myCounter = 0;

  time2 = MPI_Wtime();
  if (job->taskid > 0 || job->numtasks == 1) while (myCounter < total_tasks) {

    increment = (total_tasks - myCounter) / job->numtasks / 2;
    if (increment <  1) increment =  1;
    if (increment > 64) increment = 64;
    increment = 1; // set small for small systems with large basis sets
    myCounter = GetLongCounter(winRank, increment , &win);
    begin_task = (myCounter             < total_tasks) ? myCounter : total_tasks;
    end_task   = (myCounter + increment < total_tasks) ? myCounter + increment : total_tasks;
 
    for (q = begin_task; q < end_task; q++) {

      i = q / pair_q->tot;
      j = q - i * pair_q->tot;
      atm_n1 = pair_q->cell1[i];
      atm_n2 = pair_q->cell2[i];
      lat_n1 = pair_q->latt1[i];
      lat_n2 = pair_q->latt2[i];
      atm_n3 = pair_q->cell1[j];
      atm_n4 = pair_q->cell2[j];
      lat_n3 = pair_q->latt1[j];
      lat_n4 = pair_q->latt2[j];

      if (lat_n4 >= R_tables->last_vector || lat_n2 >= R_tables->last_vector) continue;

      if (pair_q->dist[i] > job->itol1 || pair_q->dist[j] > job->itol1) continue;

      nd1 = atoms->bfnnumb_sh[atm_n1];
      nd2 = atoms->bfnnumb_sh[atm_n2];
      nd3 = atoms->bfnnumb_sh[atm_n3];
      nd4 = atoms->bfnnumb_sh[atm_n4];
      dimi = nd1 * nd2 * nd3 * nd4;
      QUAD_TRAN Quad;
      INTEGRAL_LIST integral_list;
      Quad.tot = 8 * symmetry->number_of_operators;
      allocate_QUAD_TRAN(&Quad,job,file);
      allocate_integral_list(&integral_list, dimi, job, file);

      time4 = MPI_Wtime();
      generate_c_quads(pair_p,&Quad,atm_n1,atm_n2,atm_n3,atm_n4,lat_n1,lat_n2,lat_n3,lat_n4,atoms,atom_p,symmetry,R_tables,job,file);
      time3 += MPI_Wtime() - time4;

      if (Quad.tot > 0)   { 
      time6 = MPI_Wtime();
      integral_list.num = 0;
      unique_quads++;
      all_quads += Quad.tot; 
      total_quads += Quad.tot; 
      if (pair_p->uniq[Quad.latt2[0] * dim2 + Quad.cell1[0] * dim1 + Quad.cell2[0]] == 0) total_quads--;
      nshells = atoms->nshel_sh[Quad.cell1[0]] * atoms->nshel_sh[Quad.cell2[0]] * \
                atoms->nshel_sh[Quad.cell3[0]] * atoms->nshel_sh[Quad.cell4[0]];
      int start_index[nshells];

      for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
      //shell_overlap(start_index,&Quad,atoms,shells,R,job,file);
      //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 1;
      //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
      if (job->scf_direct == 0)
      shell_screen1(start_index,S1,pair_p,&Quad,atoms,shells,job,file);
      //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 1;
      //if (job->taskid == 1) {printf("%2d %2d  ",lat_n2,lat_n4);for(i1 = 0; i1 < nshells; i1++) printf("%2d",start_index[i1]);printf("\n");}
      //else if (job->scf_direct == 1)
      //shell_screen_crystal_coulomb_direct(start_index,S1,F,pair_p,&Quad,atoms,shells,R_tables,symmetry,job,file);
      //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 1;
      time5 += MPI_Wtime() - time6;

      time8 = MPI_Wtime();
      //count = 0; 
      integrals_crystal_coulomb_ijkl(&integral_list,&Quad,start_index,R,G,atoms,shells,gaussians,symmetry,crystal,job,file);
      //integrals_crystal_coulomb_ijkl(&integral_list,&Quad,start_index,&count,R,G,atoms,shells,gaussians,symmetry,crystal,job,file);
      total_integrals += integral_list.num;
      time7 += MPI_Wtime() - time8;
      time10 = MPI_Wtime();
      contract_integrals_crystal_coulomb_ijkl(Fock_2c_buffer,&integral_list,F,pair_p,&Quad,atoms,shells,symmetry,R_tables,job,file);
      time9 += MPI_Wtime() - time10;

      time12 = MPI_Wtime();
      //if (job->scf_direct == 0 && integral_list.num > 0 && (job->taskid > 0 || job->numtasks == 1))
      if (integral_list.num > 0 && (job->taskid > 0 || job->numtasks == 1))
      pack_write_integrals_crystal_coulomb_ijkl(&integral_list, &Quad, integrals_2c, job, file);
      time11 += MPI_Wtime() - time12;
     }

      free_QUAD_TRAN(&Quad,job);
      free_integral_list(&integral_list,job);
    } // close loop on q
   } // while (myCounter
      time1 += MPI_Wtime() - time2;

      DestroyLongCounter(job->taskid, winRank, &win, counter);

      MPI_Allreduce(Fock_2c_buffer,Fock_2c,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      DestroyDoubleArray(&Fock_2c_buffer,&dimp_spin,job);

      if (job->taskid > 0 || job->numtasks == 1) fclose(integrals_2c); 
      job->int_exist = 1;
      //if (job->scf_direct == 0 && (job->taskid > 0 || job->numtasks == 1)) { 
      //fclose(integrals_2c); 
     //}
      //job->int_exist = 1;

      if (job->verbosity >= 1 && job->iter >= 1) \
      printf("%2d c_quads uniq %4d tot %5d int %8li quad %8.2e scrn %8.2e comp %8.2e contr %8.2e write %8.2e tot %8.2e\n",\
      job->taskid,unique_quads,total_quads,total_integrals,time3,time5,time7,time9,time11,time1);

      int all_quads_sum = 0, total_quads_sum = 0, unique_quads_sum = 0;
      MPI_Reduce(&all_quads,&all_quads_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&total_quads,&total_quads_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&unique_quads,&unique_quads_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
      if (job->taskid == 0) printf("total c_quads = %6d %6d %6d\n",unique_quads_sum,total_quads_sum,all_quads_sum);

}

void fock_matrix_crystal_compute_exchange_integrals(double *Fock_2e, double *S2, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int h, i1;
int Quad_max, h_max;
int dimp, dimf, dimp_spin, dimf_spin;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
int ip, jp, kp, lp, gj, gk, gl;
int count, dimi, nshells;
int atm_n1, atm_n2, atm_n3, atm_n4;
int lat_n1, lat_n2, lat_n3, lat_n4;
int nd1, nd2, nd3, nd4;
int unique_quads, total_quads, all_quads, largest;
int winRank;
long i, j, q;
long total_integrals;
long begin_task, end_task, total_tasks;
long *counter, myCounter, increment;
const long izero = 0;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, time11, time12;
double *Fock_2e_buffer;
double R_AB, R_CD, R_AC, R_BD;
VECTOR_DOUBLE R_AB_1e, R_CD_1e, R_AC_1e, R_BD_1e;
MPI_Win win;

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;
  AllocateDoubleArray(&Fock_2e_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2e_buffer,&dimp_spin);

  total_tasks = (long) pair_q->tot * (long) pair_q->tot;

  time1 = k_zero;
  time3 = k_zero;
  time5 = k_zero;
  time7 = k_zero;
  time9 = k_zero;
  time11 = k_zero;
  all_quads = 0;
  total_quads = 0;
  unique_quads = 0;
  total_integrals = 0;
  largest = 0;

  h_max = 43;
  h_max = 79;
  if (h_max > R_tables->last_vector) h_max = R_tables->last_vector;
  Quad_max = 8 * symmetry->number_of_operators;

  winRank = 0;
  CreateLongCounter(job->taskid, winRank, &counter, izero, &win, MPI_COMM_WORLD);
  myCounter = 0;

  time2 = MPI_Wtime();
  if (job->taskid > 0 || job->numtasks == 1) while (myCounter < total_tasks) {

    increment = (total_tasks - myCounter) / job->numtasks / 2;
    if (increment <  1) increment =  1;
    if (increment > 64) increment = 64;
    myCounter = GetLongCounter(winRank, increment , &win);
    begin_task = (myCounter             < total_tasks) ? myCounter : total_tasks;
    end_task   = (myCounter + increment < total_tasks) ? myCounter + increment : total_tasks;
 
    for (q = begin_task; q < end_task; q++) {
      for (h = 0; h < h_max; h++) {
	i = q / pair_q->tot;
	j = q - i * pair_q->tot;

	atm_n1 = pair_q->cell1[i];
	atm_n3 = pair_q->cell2[i];
	lat_n1 = 0;
	lat_n3 = pair_q->latt2[i];

	atm_n2 = pair_q->cell1[j];
	atm_n4 = pair_q->cell2[j];
	lat_n2 = h;
	lat_n4 = R_tables->sumvec[h * R_tables->margin_vector + pair_q->latt2[j]];

        if (lat_n2 >= R_tables->last_vector || lat_n3 >= R_tables->last_vector || lat_n4 >= R_tables->last_vector) continue;
        //if (lat_n2 > 400 || pair_q->latt2[j] > 400) printf("%3d %3d %3d\n",lat_n2,pair_q->latt2[j],lat_n4);
        //if (lat_n2 >= 400 || lat_n3 >= 400 || lat_n4 >= 400) continue;

  // ******************************************************************************************
  // * Initial filter of quads by distance limits and generation of e_quads                   *
  // ******************************************************************************************

        //if (pair_q->dist[i] > job->itol4) continue; // R_AC condition
        //printf("AC %li %10.4lf %10.4lf\n",i,pair_q->dist[i], job->itol4);

	//(0|h)
	R_AB_1e.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
	R_AB_1e.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
	R_AB_1e.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;
	R_AB = R_AB_1e.comp1 * R_AB_1e.comp1 + R_AB_1e.comp2 * R_AB_1e.comp2 + R_AB_1e.comp3 * R_AB_1e.comp3;

	if (R_AB > job->itol3 * job->itol3) continue; // R_AB condition
        //printf("AB %li %10.4lf %10.4lf\n",i,sqrt(R_AB), job->itol3);

	//(2g|4h+l)
	R_CD_1e.comp1 = atoms->cell_vector[atm_n3].comp1 + R->vec_ai[lat_n3].comp1 - atoms->cell_vector[atm_n4].comp1 - \
	R->vec_ai[lat_n4].comp1;
	R_CD_1e.comp2 = atoms->cell_vector[atm_n3].comp2 + R->vec_ai[lat_n3].comp2 - atoms->cell_vector[atm_n4].comp2 - \
	R->vec_ai[lat_n4].comp2;
	R_CD_1e.comp3 = atoms->cell_vector[atm_n3].comp3 + R->vec_ai[lat_n3].comp3 - atoms->cell_vector[atm_n4].comp3 - \
	R->vec_ai[lat_n4].comp3;
	R_CD = R_CD_1e.comp1 * R_CD_1e.comp1 + R_CD_1e.comp2 * R_CD_1e.comp2 + R_CD_1e.comp3 * R_CD_1e.comp3;

	if (R_CD > job->itol3 * job->itol3) continue; // R_CD condition
        //printf("AC AB CD %li %10.4lf %10.4lf %10.4lf %10.4lf\n",i,pair_q->dist[i],sqrt(R_AB),sqrt(R_CD), job->itol3);

        time4 = MPI_Wtime();
	QUAD_TRAN Quad;
	Quad.tot = Quad_max;
	allocate_QUAD_TRAN(&Quad,job,file);
	generate_e_quads(pair_p,&Quad,atm_n1,atm_n2,atm_n3,atm_n4,lat_n1,lat_n2,lat_n3,lat_n4,atoms,atom_p,symmetry,\
	R_tables,job,file);
        if (Quad.tot > largest) largest = Quad.tot;
	if (Quad.tot < 1) {
	free_QUAD_TRAN(&Quad,job);
	continue; 
       }
        time3 += MPI_Wtime() - time4;

  // ******************************************************************************************
  // * Select quads in Quad from generate_e_quads which meet distance limits                  *
  // ******************************************************************************************

  time6 = MPI_Wtime();
  int needed[Quad.tot];
  for (i1 = 0; i1 < Quad.tot; i1++) needed[i1] = 0;

  count = 0;
  for (i1 = 0; i1 < Quad.tot; i1++) {

    ip = Quad.cell1[i1];
    jp = Quad.cell2[i1];
    kp = Quad.cell3[i1];
    lp = Quad.cell4[i1];

    gj = Quad.latt2[i1];
    gk = Quad.latt3[i1];
    gl = Quad.latt4[i1];

    //(10|3h) 
    R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
    R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
    R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;

    //(10|2g)
    R_AC_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[kp].comp1 - R->vec_ai[gk].comp1;
    R_AC_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[kp].comp2 - R->vec_ai[gk].comp2;
    R_AC_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[kp].comp3 - R->vec_ai[gk].comp3;

    //(2g|4h+l)
    R_CD_1e.comp1 = atoms->cell_vector[kp].comp1 + R->vec_ai[gk].comp1 - atoms->cell_vector[lp].comp1 - R->vec_ai[gl].comp1;
    R_CD_1e.comp2 = atoms->cell_vector[kp].comp2 + R->vec_ai[gk].comp2 - atoms->cell_vector[lp].comp2 - R->vec_ai[gl].comp2;
    R_CD_1e.comp3 = atoms->cell_vector[kp].comp3 + R->vec_ai[gk].comp3 - atoms->cell_vector[lp].comp3 - R->vec_ai[gl].comp3;

    //(3h|4h+l) ITOL5
    R_BD_1e.comp1 = atoms->cell_vector[jp].comp1 + R->vec_ai[gj].comp1 - atoms->cell_vector[lp].comp1 - R->vec_ai[gl].comp1;
    R_BD_1e.comp2 = atoms->cell_vector[jp].comp2 + R->vec_ai[gj].comp2 - atoms->cell_vector[lp].comp2 - R->vec_ai[gl].comp2;
    R_BD_1e.comp3 = atoms->cell_vector[jp].comp3 + R->vec_ai[gj].comp3 - atoms->cell_vector[lp].comp3 - R->vec_ai[gl].comp3;

    R_AB = sqrt(double_vec_dot(&R_AB_1e, &R_AB_1e));
    R_CD = sqrt(double_vec_dot(&R_CD_1e, &R_CD_1e));
    R_AC = sqrt(double_vec_dot(&R_AC_1e, &R_AC_1e));
    R_BD = sqrt(double_vec_dot(&R_BD_1e, &R_BD_1e));

    if (R_AB < job->itol3 && R_CD < job->itol3 && R_AC < job->itol4 && R_BD < job->itol5 && \
    pair_p->uniq[gk * dim2 + ip * dim1 + kp] == -1) {
    needed[i1] = 1;
    count++;
   }
  }

    total_quads += count; // count all quads which are needed for unique part of Fock matrix (pair_p->uniq == -1)

  // ******************************************************************************************
  // * Transfer selected quads from Quad to quad                                              *
  // ******************************************************************************************

  QUAD_TRAN quad;
  quad.tot = Quad.tot;
  allocate_QUAD_TRAN(&quad,job,file);
  quad.tot = 0;

  for (i1 = 0; i1 < Quad.tot; i1++) {
  if ((i1 == 0 && count > 0) || needed[i1] == 1) {
  quad.p[quad.tot] = Quad.p[i1];
  quad.k[quad.tot] = Quad.k[i1];
  quad.cell1[quad.tot] = Quad.cell1[i1];
  quad.cell2[quad.tot] = Quad.cell2[i1];
  quad.cell3[quad.tot] = Quad.cell3[i1];
  quad.cell4[quad.tot] = Quad.cell4[i1];
  quad.latt1[quad.tot] = Quad.latt1[i1];
  quad.latt2[quad.tot] = Quad.latt2[i1];
  quad.latt3[quad.tot] = Quad.latt3[i1];
  quad.latt4[quad.tot] = Quad.latt4[i1];
  quad.tot++;
 }
}

  time5 += MPI_Wtime() - time6;

  // ******************************************************************************************
  // * Compute integrals and contract                                                         *
  // ******************************************************************************************

  if (quad.tot > 0) {
  time8 = MPI_Wtime();
  nd1 = atoms->bfnnumb_sh[atm_n1];
  nd2 = atoms->bfnnumb_sh[atm_n2];
  nd3 = atoms->bfnnumb_sh[atm_n3];
  nd4 = atoms->bfnnumb_sh[atm_n4];
  dimi = nd1 * nd2 * nd3 * nd4;
  INTEGRAL_LIST integral_list;
  allocate_integral_list(&integral_list, dimi, job, file);
  integral_list.num = 0;
  all_quads += quad.tot; // include all quads whether (pair_p->uniq == 0/-1)
  unique_quads++;
  nshells = atoms->nshel_sh[quad.cell1[0]] * atoms->nshel_sh[quad.cell2[0]] * \
            atoms->nshel_sh[quad.cell3[0]] * atoms->nshel_sh[quad.cell4[0]];
  int start_index[nshells];

  //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
  //shell_overlap(start_index,&quad,atoms,shells,R,job,file);
  //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 1;
  for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
  shell_screen2(start_index,S2,pair_p,&quad,R_tables,atoms,shells,job,file);
  time7 += MPI_Wtime() - time8;

  time10 = MPI_Wtime();
  count = 0; 
  integrals_crystal_exchange_ijkl(&integral_list,&quad,start_index,&count,R,atoms,shells,gaussians,symmetry,crystal,job,file);
  //integrals_exchange_ijkl1(&integral_list,&quad,start_index,&count,R,atoms,shells,gaussians,symmetry,crystal,job,file);
  total_integrals += integral_list.num;
  time9 += MPI_Wtime() - time10;
  time12 = MPI_Wtime();
  contract_integrals_crystal_exchange_ijkl(Fock_2e_buffer,&integral_list,F,pair_p,&quad,atoms,shells,symmetry,R_tables,job,file);
  time11 += MPI_Wtime() - time12;
  free_integral_list(&integral_list,job);
  } // close if (quad.tot > 0)

  free_QUAD_TRAN(&quad,job);
  free_QUAD_TRAN(&Quad,job);

  } // close loop on h
 } // close loop on q
} // while (myCounter
   time1 += MPI_Wtime() - time2;

   DestroyLongCounter(job->taskid, winRank, &win, counter);

   MPI_Allreduce(Fock_2e_buffer,Fock_2e,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   DestroyDoubleArray(&Fock_2e_buffer,&dimp_spin,job);

   if (job->verbosity > 1 && job->iter == 1) \
   printf("%2d e_quads uniq %4d tot %5d int %8li select %8.2e quad %8.2e scrn %8.2e comp %8.2e contr %8.2e tot %8.2e\n",\
   job->taskid,unique_quads,total_quads,total_integrals,time3,time5,time7,time9,time11,time1);

   int unique_quads_sum = 0, total_quads_sum = 0, all_quads_sum = 0;
   MPI_Reduce(&unique_quads,&unique_quads_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(&total_quads,&total_quads_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(&all_quads,&all_quads_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
   //if (job->taskid == 0 && job->iter == 1) 
   if (job->taskid == 0)
   printf("total e_quads = %6d %6d %6d\n",unique_quads_sum,total_quads_sum,all_quads_sum);

}

/*
void fock_matrix_molecule_compute_integrals(double *Fock_2c, double *Fock_2e, double *S1, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

int i1;
int dimp, dimf, dimp_spin, dimf_spin;
int count, dimi, nshells;
int atm_n1, atm_n2, atm_n3, atm_n4;
int lat_n1, lat_n2, lat_n3, lat_n4;
int nd1, nd2, nd3, nd4;
int all_quads, unique_quads, total_quads, largest;
int winRank;
long i, j, q;
long total_integrals;
long begin_task,end_task, total_tasks;
long *counter, myCounter, increment;
const long izero = 0;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, time11, time12;
char xx[4], yy[24] = "scf_integrals_2M.";
double *Fock_2c_buffer, *Fock_2e_buffer;
FILE *integrals_2m; 
MPI_Win win;

  sprintf(xx, "%d", job->taskid);

  if (job->scf_direct == 0 && (job->taskid > 0 || job->numtasks == 1)) integrals_2m = fopen(strcat(yy,xx), "wb");

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;
  AllocateDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  AllocateDoubleArray(&Fock_2e_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2c_buffer,&dimp_spin);
  ResetDoubleArray(Fock_2e_buffer,&dimp_spin);

  total_tasks = (long) pair_q->tot * (long) pair_q->tot;

  time1 = k_zero;
  time3 = k_zero;
  time5 = k_zero;
  time7 = k_zero;
  time9 = k_zero;
  time11 = k_zero;
  all_quads = 0;
  unique_quads = 0;
  total_quads = 0;
  total_integrals = 0;

  winRank = 0;
  CreateLongCounter(job->taskid, winRank, &counter, izero, &win, MPI_COMM_WORLD);
  myCounter = 0;

  time2 = MPI_Wtime();
  if (job->taskid > 0 || job->numtasks == 1) while (myCounter < total_tasks) {

    increment = (total_tasks - myCounter) / job->numtasks / 2;
    if (increment <  1) increment =  1;
    if (increment > 64) increment = 64;
    increment = 1; // set small for small systems with large basis sets
    myCounter = GetLongCounter(winRank, increment , &win);
    begin_task = (myCounter             < total_tasks) ? myCounter : total_tasks;
    end_task   = (myCounter + increment < total_tasks) ? myCounter + increment : total_tasks;
 
    for (q = begin_task; q < end_task; q++) {

      i = q / pair_q->tot;
      j = q - i * pair_q->tot;
      atm_n1 = pair_q->cell1[i];
      atm_n2 = pair_q->cell2[i];
      atm_n3 = pair_q->cell1[j];
      atm_n4 = pair_q->cell2[j];
      lat_n1 = 0;
      lat_n2 = 0;
      lat_n3 = 0;
      lat_n4 = 0;

      if (job->pms == 1 && atm_n1 * atoms->number_of_atoms_in_unit_cell + atm_n2 > \
      atm_n3 * atoms->number_of_atoms_in_unit_cell + atm_n4) continue;

      nd1 = atoms->bfnnumb_sh[atm_n1];
      nd2 = atoms->bfnnumb_sh[atm_n2];
      nd3 = atoms->bfnnumb_sh[atm_n3];
      nd4 = atoms->bfnnumb_sh[atm_n4];
      dimi = nd1 * nd2 * nd3 * nd4;
      QUAD_TRAN Quad;
      INTEGRAL_LIST integral_list;
      Quad.tot = 8 * symmetry->number_of_operators;
      allocate_QUAD_TRAN(&Quad,job,file);
      allocate_integral_list(&integral_list, dimi, job, file);

      time4 = MPI_Wtime();
      generate_molecule_quads(pair_p,&Quad,atm_n1,atm_n2,atm_n3,atm_n4,lat_n1,lat_n2,lat_n3,lat_n4,atoms,atom_p,symmetry,\
      R_tables,job,file);
      time3 += MPI_Wtime() - time4;

      if (Quad.tot > 0)   { 
      time6 = MPI_Wtime();
      integral_list.num = 0;
      total_quads += Quad.tot; 
      unique_quads++;
      nshells = atoms->nshel_sh[Quad.cell1[0]] * atoms->nshel_sh[Quad.cell2[0]] * \
                atoms->nshel_sh[Quad.cell3[0]] * atoms->nshel_sh[Quad.cell4[0]];
      int start_index[nshells];
      for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 1;
      ////for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
      ////if (job->scf_direct == 0 || job->iter == 1) 
      ////shell_screen1(start_index,S1,pair_p,&Quad,atoms,shells,job,file);
      ////else if (job->scf_direct == 1 && job->iter > 1) 
      ////shell_screen_direct(start_index,S1,F,pair_p,&Quad,atoms,shells,symmetry,job,file);
      time5 += MPI_Wtime() - time6;

      time8 = MPI_Wtime();
      count = 0; 
      //integrals_molecule_ijkl(&integral_list,atm_n1,atm_n2,atm_n3,atm_n4,start_index,&count,R,atoms,shells,gaussians,\
      symmetry,crystal,job,file);
      integrals_molecule_ijkl(&integral_list,atm_n1,atm_n2,atm_n3,atm_n4,start_index,R,atoms,shells,gaussians,symmetry,job,file);
      total_integrals += integral_list.num;
      time7 += MPI_Wtime() - time8;
      time10 = MPI_Wtime();
      contract_integrals_molecule_ijkl(Fock_2c_buffer,Fock_2e_buffer,&integral_list,F,pair_p,&Quad,atoms,shells,symmetry,job,file);
      time9 += MPI_Wtime() - time10;

      time12 = MPI_Wtime();
      if (job->scf_direct == 0 && integral_list.num > 0 && (job->taskid > 0 || job->numtasks == 1))
      pack_write_molecule_2e_integrals(&integral_list, &Quad, integrals_2m, job, file);
      time11 += MPI_Wtime() - time12;
     }

      free_QUAD_TRAN(&Quad,job);
      free_integral_list(&integral_list,job);
    } // close loop on q
   } // while (myCounter
      time1 += MPI_Wtime() - time2;

      DestroyLongCounter(job->taskid, winRank, &win, counter);

      MPI_Allreduce(Fock_2c_buffer,Fock_2c,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      DestroyDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
      MPI_Allreduce(Fock_2e_buffer,Fock_2e,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      DestroyDoubleArray(&Fock_2e_buffer,&dimp_spin,job);

      if (job->scf_direct == 0 && (job->taskid > 0 || job->numtasks == 1)) { fclose(integrals_2m); job->int_exist == 1; }
      //if (job->taskid > 0 || job->numtasks == 1) fclose(integrals_2m); 
      //job->int_exist = 1;

      if (job->taskid == 0 && job->iter == 1) \
      printf("%2d m_quads uniq %4d tot %5d int %8li quad %8.2e scrn %8.2e comp %8.2e contr %8.2e write %8.2e tot %8.2e\n",\
      job->taskid,unique_quads,total_quads,total_integrals,time3,time5,time7,time9,time11,time1);

}
*/

void fock_matrix_crystal_read_coulomb_integrals(double *Fock_2c, double *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * This routine reads crystal Coulomb integrals from disk and contracts them              *
  // ******************************************************************************************

int dimf, dimp, dimp_spin;
double *Fock_2c_buffer;
double time1, time2;
size_t result;
FILE *integrals_2c; 
char xx[4], yy[24] = "scf_integrals_2C.";

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file);
  dimp_spin = dimp * job->spin_dim;
  AllocateDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2c_buffer,&dimp_spin);

  sprintf(xx, "%d", job->taskid);

  if (job->taskid > 0 || job->numtasks == 1) { 

  if ((integrals_2c = fopen(strcat(yy,xx), "rb")) == NULL) {
  printf("Failed to open file %s in coulomb_matrix_crystal_read_integrals\n",yy);
  fprintf(file.out,"Failed to open file %s in coulomb_matrix_crystal_read_integrals\n",yy);
  MPI_Finalize();
  exit(1);
 }

  while (!feof(integrals_2c)) {
  QUAD_TRAN Quad;
  INTEGRAL_LIST integral_list;
  result = fread(&Quad.tot,sizeof(int),1,integrals_2c);
  result = fread(&integral_list.num,sizeof(int),1,integrals_2c);
  if (feof(integrals_2c)) break;
  allocate_QUAD_TRAN(&Quad,job,file);
  allocate_integral_list(&integral_list, integral_list.num, job, file);
  read_unpack_integrals_crystal_coulomb_ijkl(&integral_list, &Quad, integrals_2c, job, file);
  contract_integrals_crystal_coulomb_ijkl(Fock_2c_buffer,&integral_list,F,pair_p,&Quad,atoms,shells,symmetry,R_tables,job,file);
  free_QUAD_TRAN(&Quad,job);
  free_integral_list(&integral_list,job);
 } //while (!feof

  fclose(integrals_2c);
}

  MPI_Allreduce(Fock_2c_buffer,Fock_2c,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&Fock_2c_buffer,&dimp_spin,job);

}

void fock_matrix_crystal_read_exchange_integrals(double *Fock_2e, double *F, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * This routine reads crystal exchange integrals from disk and contracts them             *
  // * routine is not used and not tested                                                     *
  // ******************************************************************************************

int dimf, dimp, dimp_spin;
double *Fock_2e_buffer;
double time1, time2;
size_t result;
FILE *integrals_2e; 
char xx[4], yy[24] = "scf_integrals_2E.";

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file);
  dimp_spin = dimp * job->spin_dim;
  AllocateDoubleArray(&Fock_2e_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2e_buffer,&dimp_spin);
	  
  sprintf(xx, "%d", job->taskid);

  if (job->taskid > 0 || job->numtasks == 1) {

  if ((integrals_2e = fopen(strcat(yy,xx), "rb")) == NULL) {
  printf("Failed to open file %s in exchange_matrix_crystal_read_integrals\n",yy);
  fprintf(file.out,"Failed to open file %s in exchange_matrix_crystal_read_integrals\n",yy);
  MPI_Finalize();
  exit(1);
 }

  while (!feof(integrals_2e)) {

  QUAD_TRAN Quad;
  INTEGRAL_LIST integral_list;
  result = fread(&Quad.tot,sizeof(int),1,integrals_2e);
  result = fread(&integral_list.num,sizeof(int),1,integrals_2e);
  if (feof(integrals_2e)) break;
  allocate_QUAD_TRAN(&Quad,job,file);
  allocate_integral_list(&integral_list, integral_list.num, job, file);
  read_unpack_integrals_crystal_exchange_ijkl(&integral_list, &Quad, integrals_2e, job, file);
  contract_integrals_crystal_exchange_ijkl(Fock_2e_buffer,&integral_list,F,pair_p,&Quad,atoms,shells,symmetry,R_tables,job,file);
  free_QUAD_TRAN(&Quad,job);
  free_integral_list(&integral_list,job);

 } //while (!feof

  fclose(integrals_2e);

 }

  MPI_Allreduce(Fock_2e_buffer,Fock_2e,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&Fock_2e_buffer,&dimp_spin,job);

}

/*
void fock_matrix_molecule_read_integrals(double *Fock_2c, double *Fock_2e, double *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * This routine reads molecule integrals from disk and contracts them                     *
  // ******************************************************************************************

int dimf, dimp, dimp_spin;
double *Fock_2e_buffer, *Fock_2c_buffer;
double time1, time2;
FILE *integrals_2m; 
char xx[4], yy[24] = "scf_integrals_2M.";

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file);
  dimp_spin = dimp * job->spin_dim;
  AllocateDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2c_buffer,&dimp_spin);
  AllocateDoubleArray(&Fock_2e_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2e_buffer,&dimp_spin);
	  
  sprintf(xx, "%d", job->taskid);

  if (job->taskid > 0 || job->numtasks == 1) { 

  if ((integrals_2m = fopen(strcat(yy,xx), "rb")) == NULL) {
  printf("Failed to open file %s in coulomb_exchange_matrix_molecule_read_integrals\n",yy);
  fprintf(file.out,"Failed to open file %s in coulomb_exchange_matrix_molecule_read_integrals\n",yy);
  MPI_Finalize();
  exit(1);
 }

  while (!feof(integrals_2m)) {

  QUAD_TRAN Quad;
  size_t result;
  INTEGRAL_LIST integral_list;
  result = fread(&Quad.tot,sizeof(int),1,integrals_2m);
  result = fread(&integral_list.num,sizeof(int),1,integrals_2m);
  if (feof(integrals_2m)) break;
  allocate_QUAD_TRAN(&Quad,job,file);
  allocate_integral_list(&integral_list, integral_list.num, job, file);
  read_unpack_molecule_2e_integrals(&integral_list, &Quad, integrals_2m, job, file);
  contract_integrals_molecule_ijkl(Fock_2c_buffer,Fock_2e_buffer,&integral_list,F,pair_p,&Quad,atoms,shells,symmetry,job,file);
  free_QUAD_TRAN(&Quad,job);
  free_integral_list(&integral_list,job);

 } //while (!feof

  fclose(integrals_2m);
}

  MPI_Allreduce(Fock_2c_buffer,Fock_2c,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  MPI_Allreduce(Fock_2e_buffer,Fock_2e,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&Fock_2e_buffer,&dimp_spin,job);

}
*/

void fock_matrix_crystal_compute_coulomb_integrals_no_sym(double *Fock_2c, double *S1, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i1, m;
int dimp, dimf, dimp_spin, dimf_spin;
int dimi, nshells;
int atm_n1, atm_n2, atm_n3, atm_n4;
int lat_n1, lat_n2, lat_n3, lat_n4;
int nd1, nd2, nd3, nd4;
int total_quads;
int Fock_temp_offset, Density_matrix_offset;
int F_ptr, D_ptr;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
int winRank;
long i, j, q;
long total_integrals;
long begin_task,end_task, total_tasks;
long *counter, myCounter, increment;
const long izero = 0;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, time11, time12;
char xx[4], yy[30] = "scf_integrals_2C_no_sym.";
double *Fock_2c_buffer;
FILE *integrals_2c; 
MPI_Win win;

  sprintf(xx, "%d", job->taskid);

  if (job->scf_direct == 0 && (job->taskid > 0 || job->numtasks == 1)) integrals_2c = fopen(strcat(yy,xx), "wb");

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;
  AllocateDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2c_buffer,&dimp_spin);

  total_tasks = (long) pair_q->tot * (long) pair_q->tot;

  time1 = k_zero;
  time3 = k_zero;
  time5 = k_zero;
  time7 = k_zero;
  time9 = k_zero;
  time11 = k_zero;
  total_quads = 0;
  total_integrals = 0;

  winRank = 0;
  CreateLongCounter(job->taskid, winRank, &counter, izero, &win, MPI_COMM_WORLD);
  myCounter = 0;

  time2 = MPI_Wtime();
  if (job->taskid > 0 || job->numtasks == 1) while (myCounter < total_tasks) {

    increment = (total_tasks - myCounter) / job->numtasks / 2;
    if (increment <  1) increment =  1;
    if (increment > 64) increment = 64;
    increment = 1; // set small for small systems with large basis sets
    myCounter = GetLongCounter(winRank, increment , &win);
    begin_task = (myCounter             < total_tasks) ? myCounter : total_tasks;
    end_task   = (myCounter + increment < total_tasks) ? myCounter + increment : total_tasks;
 
    for (q = begin_task; q < end_task; q++) {

      i = q / pair_q->tot;
      j = q - i * pair_q->tot;
      atm_n1 = pair_q->cell1[i];
      atm_n2 = pair_q->cell2[i];
      lat_n1 = pair_q->latt1[i];
      lat_n2 = pair_q->latt2[i];
      atm_n3 = pair_q->cell1[j];
      atm_n4 = pair_q->cell2[j];
      lat_n3 = pair_q->latt1[j];
      lat_n4 = pair_q->latt2[j];

      if (lat_n4 >= R_tables->last_vector || lat_n2 >= R_tables->last_vector) continue;

      if (pair_q->dist[i] > job->itol1 || pair_q->dist[j] > job->itol1) continue;

      if (pair_p->uniq[lat_n2 * dim2 + atm_n1 * dim1 + atm_n2] == -1) {

      QUAD_TRAN Quad;
      Quad.tot = 1;
      allocate_QUAD_TRAN(&Quad,job,file);

      Quad.cell1[0] = atm_n1;
      Quad.cell2[0] = atm_n2;
      Quad.cell3[0] = atm_n3;
      Quad.cell4[0] = atm_n4;
      Quad.latt1[0] = lat_n1;
      Quad.latt2[0] = lat_n2;
      Quad.latt3[0] = lat_n3;
      Quad.latt4[0] = lat_n4;

      time6 = MPI_Wtime();
      nd1 = atoms->bfnnumb_sh[atm_n1];
      nd2 = atoms->bfnnumb_sh[atm_n2];
      nd3 = atoms->bfnnumb_sh[atm_n3];
      nd4 = atoms->bfnnumb_sh[atm_n4];
      dimi = nd1 * nd2 * nd3 * nd4;
      INTEGRAL_LIST integral_list;
      allocate_integral_list(&integral_list, dimi, job, file);
      integral_list.num = 0;
      total_quads++; 
      nshells = atoms->nshel_sh[Quad.cell1[0]] * atoms->nshel_sh[Quad.cell2[0]] * \
                atoms->nshel_sh[Quad.cell3[0]] * atoms->nshel_sh[Quad.cell4[0]];
      int start_index[nshells];

      //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
      //shell_overlap(start_index,&Quad,atoms,shells,R,job,file);
      for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 1;
      //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
      //shell_screen1(start_index,S1,pair_p,&Quad,atoms,shells,job,file);
      time5 += MPI_Wtime() - time6;

      time8 = MPI_Wtime();
      //count = 0; 
      integrals_crystal_coulomb_ijkl(&integral_list,&Quad,start_index,R,G,atoms,shells,gaussians,symmetry,crystal,job,file);
      //integrals_crystal_coulomb_ijkl(&integral_list,&Quad,start_index,&count,R,G,atoms,shells,gaussians,symmetry,crystal,job,file);
      total_integrals += integral_list.num;
      time7 += MPI_Wtime() - time8;

      time10 = MPI_Wtime();
      F_ptr = pair_p->Off[pair_p->Ptr[lat_n2 * dim2 + atm_n1 * dim1 + atm_n2]];
      D_ptr = pair_p->off[pair_p->ptr[lat_n4 * dim2 + atm_n3 * dim1 + atm_n4]];
      for (m = 0; m < integral_list.num; m++) {
      Fock_temp_offset      =         nd2 * integral_list.i[m] + integral_list.j[m];
      Density_matrix_offset = D_ptr + nd4 * integral_list.k[m] + integral_list.l[m];
      Fock_2c_buffer[F_ptr + Fock_temp_offset] += integral_list.value[m] * F[Density_matrix_offset];
     }

      time12 = MPI_Wtime();
      if (job->scf_direct == 0 && integral_list.num > 0 && (job->taskid > 0 || job->numtasks == 1))
      pack_write_integrals_crystal_coulomb_ijkl(&integral_list, &Quad, integrals_2c, job, file);
      time11 += MPI_Wtime() - time12;

      free_integral_list(&integral_list,job);
      free_QUAD_TRAN(&Quad,job);
     } // close if (pair_q->uniq == -1)
    } // close loop on q
   } // while (myCounter
      time1 += MPI_Wtime() - time2;

      DestroyLongCounter(job->taskid, winRank, &win, counter);

      MPI_Allreduce(Fock_2c_buffer,Fock_2c,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      DestroyDoubleArray(&Fock_2c_buffer,&dimp_spin,job);

      if (job->taskid > 0 || job->numtasks == 1) fclose(integrals_2c); 
      job->int_exist_no_sym = 1;

      if (job->iter == 1) \
      printf("%2d c_quads no sym tot %5d int %8li quad %8.2e scrn %8.2e comp %8.2e contr %8.2e write %8.2e tot %8.2e\n",\
      job->taskid,total_quads,total_integrals,time3,time5,time7,time9,time11,time1);

      int total_quads_sum = 0;
      MPI_Reduce(&total_quads,&total_quads_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
      if (job->taskid == 0) printf("total c_quads no sym = %6d\n",total_quads_sum);

}

void fock_matrix_crystal_compute_exchange_integrals_no_sym(double *Fock_2e, double *S2, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

int gl0, h, m, i1;
int h_max;
int dimp, dimf, dimp_spin, dimf_spin;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
int count, dimi, nshells;
int atm_n1, atm_n2, atm_n3, atm_n4;
int lat_n1, lat_n2, lat_n3, lat_n4;
int nd1, nd2, nd3, nd4;
int total_quads;
int Fock_temp_offset, Density_matrix_offset;
int F_ptr, D_ptr;
int winRank;
long i, j, q;
long total_integrals;
long begin_task,end_task, total_tasks;
long *counter, myCounter, increment;
const long izero = 0;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10;
double *Fock_2e_buffer;
double R_AB, R_CD, R_AC, R_BD;
VECTOR_DOUBLE R_AB_1e, R_CD_1e, R_AC_1e, R_BD_1e;
MPI_Win win;

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;
  AllocateDoubleArray(&Fock_2e_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2e_buffer,&dimp_spin);

  total_tasks = (long) pair_q->tot * (long) pair_q->tot;

  time1 = k_zero;
  time3 = k_zero;
  time5 = k_zero;
  time7 = k_zero;
  time9 = k_zero;
  total_quads = 0;
  total_integrals = 0;

  h_max = R_tables->last_vector;

  winRank = 0;
  CreateLongCounter(job->taskid, winRank, &counter, izero, &win, MPI_COMM_WORLD);
  myCounter = 0;

  time2 = MPI_Wtime();
  if (job->taskid > 0 || job->numtasks == 1) while (myCounter < total_tasks) {

   increment = (total_tasks - myCounter) / job->numtasks / 2;
   if (increment <  1) increment =  1;
   if (increment > 64) increment = 64;
   myCounter = GetLongCounter(winRank, increment , &win);
   begin_task = (myCounter             < total_tasks) ? myCounter : total_tasks;
   end_task   = (myCounter + increment < total_tasks) ? myCounter + increment : total_tasks;

   for (q = begin_task; q < end_task; q++) {
     for (h = 0; h < h_max; h++) {
       i = q / pair_q->tot;
       j = q - i * pair_q->tot;

       atm_n1 = pair_q->cell1[i];
       atm_n2 = pair_q->cell2[i];
       lat_n1 = 0;
       lat_n2 = pair_q->latt2[i];

       atm_n3 = pair_q->cell1[j];
       atm_n4 = pair_q->cell2[j];
       lat_n3 = h;
       lat_n4 = R_tables->sumvec[h * R_tables->margin_vector + pair_q->latt2[j]];

       if (lat_n2 >= R_tables->last_vector || lat_n3 >= R_tables->last_vector || lat_n4 >= R_tables->last_vector) continue;

  // ******************************************************************************************
  // * Filter quads by distance limits                                                       *
  // ******************************************************************************************

       time4 = MPI_Wtime();

       //(10|3h) 
       R_AB_1e.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n2].comp1 - R->vec_ai[lat_n2].comp1;
       R_AB_1e.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n2].comp2 - R->vec_ai[lat_n2].comp2;
       R_AB_1e.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n2].comp3 - R->vec_ai[lat_n2].comp3;

       //(10|2g)
       R_AC_1e.comp1 = atoms->cell_vector[atm_n1].comp1 - atoms->cell_vector[atm_n3].comp1 - R->vec_ai[lat_n3].comp1;
       R_AC_1e.comp2 = atoms->cell_vector[atm_n1].comp2 - atoms->cell_vector[atm_n3].comp2 - R->vec_ai[lat_n3].comp2;
       R_AC_1e.comp3 = atoms->cell_vector[atm_n1].comp3 - atoms->cell_vector[atm_n3].comp3 - R->vec_ai[lat_n3].comp3;

       //(2g|4h+l)
       R_CD_1e.comp1 = atoms->cell_vector[atm_n3].comp1 + R->vec_ai[lat_n3].comp1 - atoms->cell_vector[atm_n4].comp1 - \
       R->vec_ai[lat_n4].comp1;
       R_CD_1e.comp2 = atoms->cell_vector[atm_n3].comp2 + R->vec_ai[lat_n3].comp2 - atoms->cell_vector[atm_n4].comp2 - \
       R->vec_ai[lat_n4].comp2;
       R_CD_1e.comp3 = atoms->cell_vector[atm_n3].comp3 + R->vec_ai[lat_n3].comp3 - atoms->cell_vector[atm_n4].comp3 - \
       R->vec_ai[lat_n4].comp3;

       //(3h|4h+l)
       R_BD_1e.comp1 = atoms->cell_vector[atm_n2].comp1 + R->vec_ai[lat_n2].comp1 - atoms->cell_vector[atm_n4].comp1 - \
       R->vec_ai[lat_n4].comp1;
       R_BD_1e.comp2 = atoms->cell_vector[atm_n2].comp2 + R->vec_ai[lat_n2].comp2 - atoms->cell_vector[atm_n4].comp2 - \
       R->vec_ai[lat_n4].comp2;
       R_BD_1e.comp3 = atoms->cell_vector[atm_n2].comp3 + R->vec_ai[lat_n2].comp3 - atoms->cell_vector[atm_n4].comp3 - \
       R->vec_ai[lat_n4].comp3;

       R_AB = sqrt(double_vec_dot(&R_AB_1e, &R_AB_1e));
       R_CD = sqrt(double_vec_dot(&R_CD_1e, &R_CD_1e));
       R_AC = sqrt(double_vec_dot(&R_AC_1e, &R_AC_1e));
       R_BD = sqrt(double_vec_dot(&R_BD_1e, &R_BD_1e));

       time3 += MPI_Wtime() - time4;

       if ((R_AB > job->itol3 || R_CD > job->itol3 || R_AC > job->itol4 || R_BD > job->itol5) || \
       pair_p->uniq[lat_n3 * dim2 + atm_n1 * dim1 + atm_n3] != -1) continue;

       time6 = MPI_Wtime();
       nd1 = atoms->bfnnumb_sh[atm_n1];
       nd2 = atoms->bfnnumb_sh[atm_n2];
       nd3 = atoms->bfnnumb_sh[atm_n3];
       nd4 = atoms->bfnnumb_sh[atm_n4];
       dimi = nd1 * nd2 * nd3 * nd4;
       QUAD_TRAN Quad;
       INTEGRAL_LIST integral_list;

       Quad.tot = 1;
       allocate_QUAD_TRAN(&Quad,job,file);
       allocate_integral_list(&integral_list, dimi, job, file);
       integral_list.num = 0;

       Quad.cell1[0] = atm_n1;
       Quad.cell2[0] = atm_n2;
       Quad.cell3[0] = atm_n3;
       Quad.cell4[0] = atm_n4;
       Quad.latt1[0] = lat_n1;
       Quad.latt2[0] = lat_n2;
       Quad.latt3[0] = lat_n3;
       Quad.latt4[0] = lat_n4;

       total_quads++; 
       nshells = atoms->nshel_sh[Quad.cell1[0]] * atoms->nshel_sh[Quad.cell2[0]] * \
                 atoms->nshel_sh[Quad.cell3[0]] * atoms->nshel_sh[Quad.cell4[0]];
       int start_index[nshells];

       //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
       //shell_overlap(start_index,&quad,atoms,shells,R,job,file);
       for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 1;
       //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
       //shell_screen2(start_index,S2,pair_p,&quad,atoms,shells,job,file);
       time5 += MPI_Wtime() - time6;
       count = 0; 
       time8 = MPI_Wtime();
       integrals_crystal_exchange_ijkl(&integral_list,&Quad,start_index,&count,R,atoms,shells,gaussians,symmetry,crystal,job,file);
       total_integrals += integral_list.num;
       time7 += MPI_Wtime() - time8;

       time10 = MPI_Wtime();
       F_ptr = pair_p->Off[pair_p->Ptr[lat_n3 * dim2 + atm_n1 * dim1 + atm_n3]];
       gl0 = R_tables->diffvec[lat_n4 * R_tables->margin_vector + lat_n2];
       D_ptr = pair_p->off[pair_p->ptr[gl0 * dim2 + atm_n2 * dim1 + atm_n4]];
       for (m = 0; m < integral_list.num; m++) {
       Fock_temp_offset      =         nd3 * integral_list.i[m] + integral_list.k[m];
       Density_matrix_offset = D_ptr + nd4 * integral_list.j[m] + integral_list.l[m];
       Fock_2e_buffer[F_ptr + Fock_temp_offset] -= integral_list.value[m] * F[Density_matrix_offset] / two;
      }
       time9 += MPI_Wtime() - time10;
       free_QUAD_TRAN(&Quad,job);
       free_integral_list(&integral_list,job);

       } // close loop on h
      } // close loop on q
     } // while (myCounter
       time1 += MPI_Wtime() - time2;

   DestroyLongCounter(job->taskid, winRank, &win, counter);

   MPI_Allreduce(Fock_2e_buffer,Fock_2e,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   DestroyDoubleArray(&Fock_2e_buffer,&dimp_spin,job);

   if (job->verbosity >= 1 && job->iter == 1) \
   printf("%2d e_quads nosym tot %5d int %8li filter %8.2e scrn %8.2e comp %8.2e contr %8.2e tot %8.2e\n",\
   job->taskid,total_quads,total_integrals,time3,time5,time7,time9,time1);

   int total_quads_sum = 0, all_quads_sum = 0;
   MPI_Reduce(&total_quads,&total_quads_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
   if (job->taskid == 0 && job->iter == 1) printf("total e_quads no sym = %6d\n",total_quads_sum);
}

/*
void fock_matrix_molecule_compute_integrals_no_sym(double *Fock_2c, double *Fock_2e, double *S1, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, JOB_PARAM *job, FILES file)

{

int i1, m;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dimp, dimf, dimp_spin, dimf_spin;
int count, dimi, nshells;
int atm_n1, atm_n2, atm_n3, atm_n4;
int nd1, nd2, nd3, nd4;
int total_quads;
int Fock_temp_offset, Density_matrix_offset;
int F_ptr, D_ptr;
int winRank;
long i, j, q;
long total_integrals;
long begin_task,end_task, total_tasks;
long *counter, myCounter, increment;
const long izero = 0;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10;
char xx[4], yy[30] = "scf_integrals_2M_no_sym.";
double *Fock_2c_buffer, *Fock_2e_buffer;
FILE *integrals_2m; 
MPI_Win win;

  sprintf(xx, "%d", job->taskid);

  if (job->scf_direct == 0 && (job->taskid > 0 || job->numtasks == 1)) integrals_2m = fopen(strcat(yy,xx), "wb");

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;
  AllocateDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  AllocateDoubleArray(&Fock_2e_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2c_buffer,&dimp_spin);
  ResetDoubleArray(Fock_2e_buffer,&dimp_spin);

  total_tasks = (long) pair_q->tot * (long) pair_q->tot;

  time1 = k_zero;
  time3 = k_zero;
  time5 = k_zero;
  time7 = k_zero;
  time9 = k_zero;
  total_quads = 0;
  total_integrals = 0;

  winRank = 0;
  CreateLongCounter(job->taskid, winRank, &counter, izero, &win, MPI_COMM_WORLD);
  myCounter = 0;

  time2 = MPI_Wtime();
  if (job->taskid > 0 || job->numtasks == 1) while (myCounter < total_tasks) {

    increment = (total_tasks - myCounter) / job->numtasks / 2;
    if (increment <  1) increment =  1;
    if (increment > 64) increment = 64;
    increment = 1; // set small for small systems with large basis sets
    myCounter = GetLongCounter(winRank, increment , &win);
    begin_task = (myCounter             < total_tasks) ? myCounter : total_tasks;
    end_task   = (myCounter + increment < total_tasks) ? myCounter + increment : total_tasks;
 
    for (q = begin_task; q < end_task; q++) {

      i = q / pair_q->tot;
      j = q - i * pair_q->tot;
      atm_n1 = pair_q->cell1[i];
      atm_n2 = pair_q->cell2[i];
      atm_n3 = pair_q->cell1[j];
      atm_n4 = pair_q->cell2[j];

      if (pair_p->uniq[atm_n1 * dim1 + atm_n2] == -1 || pair_p->uniq[atm_n1 * dim1 + atm_n3] == -1) {

      time4 = MPI_Wtime();
      QUAD_TRAN Quad;
      Quad.tot = 1;
      allocate_QUAD_TRAN(&Quad,job,file);

      Quad.cell1[0] = atm_n1;
      Quad.cell2[0] = atm_n2;
      Quad.cell3[0] = atm_n3;
      Quad.cell4[0] = atm_n4;

      nd1 = atoms->bfnnumb_sh[atm_n1];
      nd2 = atoms->bfnnumb_sh[atm_n2];
      nd3 = atoms->bfnnumb_sh[atm_n3];
      nd4 = atoms->bfnnumb_sh[atm_n4];

      dimi = nd1 * nd2 * nd3 * nd4;
      INTEGRAL_LIST integral_list;
      allocate_integral_list(&integral_list, dimi, job, file);
      integral_list.num = 0;
      total_quads++; 
      nshells = atoms->nshel_sh[atm_n1] * atoms->nshel_sh[atm_n2] * atoms->nshel_sh[atm_n3] * atoms->nshel_sh[atm_n4];
      int start_index[nshells];
      for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 1;
      //for (i1 = 0; i1 < nshells; i1++) start_index[i1] = 0;
      //shell_screen1(start_index,S1,pair_p,&Quad,atoms,shells,job,file);
      time3 += MPI_Wtime() - time4;

      time6 = MPI_Wtime();
      count = 0; 
      //integrals_molecule_ijkl(&integral_list,atm_n1,atm_n2,atm_n3,atm_n4,start_index,&count,R,atoms,shells,gaussians,\
      symmetry,crystal,job,file);
      integrals_molecule_ijkl(&integral_list,atm_n1,atm_n2,atm_n3,atm_n4,start_index,R,atoms,shells,gaussians,symmetry,job,file);
      total_integrals += integral_list.num;
      time5 += MPI_Wtime() - time6;

      time8 = MPI_Wtime();
      if (pair_p->uniq[atm_n1 * dim1 + atm_n2] == -1) {
      F_ptr = pair_p->Off[pair_p->Ptr[atm_n1 * dim1 + atm_n2]];
      D_ptr = pair_p->off[pair_p->ptr[atm_n3 * dim1 + atm_n4]];
      for (m = 0; m < integral_list.num; m++) {
        Fock_temp_offset      =         nd2 * integral_list.i[m] + integral_list.j[m];
        Density_matrix_offset = D_ptr + nd4 * integral_list.k[m] + integral_list.l[m];
        Fock_2c_buffer[F_ptr + Fock_temp_offset] += integral_list.value[m] * F[Density_matrix_offset];
       }
      }
      if (pair_p->uniq[atm_n1 * dim1 + atm_n3] == -1) {
      F_ptr = pair_p->Off[pair_p->Ptr[atm_n1 * dim1 + atm_n3]];
      D_ptr = pair_p->off[pair_p->ptr[atm_n2 * dim1 + atm_n4]];
      for (m = 0; m < integral_list.num; m++) {
        Fock_temp_offset      =         nd3 * integral_list.i[m] + integral_list.k[m];
        Density_matrix_offset = D_ptr + nd4 * integral_list.j[m] + integral_list.l[m];
        Fock_2e_buffer[F_ptr + Fock_temp_offset] -= integral_list.value[m] * F[Density_matrix_offset] / two;
       }
      }
      time7 += MPI_Wtime() - time8;

      time10 = MPI_Wtime();
      if (job->scf_direct == 0 && integral_list.num > 0 && (job->taskid > 0 || job->numtasks == 1))
      pack_write_molecule_2e_integrals(&integral_list, &Quad, integrals_2m, job, file);
      time9 += MPI_Wtime() - time10;

      free_integral_list(&integral_list,job);
      free_QUAD_TRAN(&Quad,job);
    } // close if (pair->uniq == -1)
   } // close loop on q
  } // while (myCounter
      time1 += MPI_Wtime() - time2;

      DestroyLongCounter(job->taskid, winRank, &win, counter);

      MPI_Allreduce(Fock_2c_buffer,Fock_2c,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      DestroyDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
      MPI_Allreduce(Fock_2e_buffer,Fock_2e,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      DestroyDoubleArray(&Fock_2e_buffer,&dimp_spin,job);

      //if (job->scf_direct == 0 && (job->taskid > 0 || job->numtasks == 1)) { fclose(integrals_2m); job->int_exist == 1; }
      if (job->taskid > 0 || job->numtasks == 1) fclose(integrals_2m); 
      job->int_exist_no_sym = 1;

      if (job->iter == 1) \
      printf("%2d m_quads no sym tot %5d int %8li scrn %8.2e comp %8.2e contr %8.2e write %8.2e tot %8.2e\n",\
      job->taskid,total_quads,total_integrals,time3,time5,time7,time9,time1);

}
*/

void fock_matrix_crystal_read_coulomb_integrals_no_sym(double *Fock_2c, double *F, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * This routine reads crystal Coulomb integrals from disk and contracts them - no sym     *
  // ******************************************************************************************

int m;
int dimf, dimp, dimp_spin;
int Fock_temp_offset, Density_matrix_offset;
int F_ptr, D_ptr;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
int atm_n1, atm_n2, atm_n3, atm_n4;
int lat_n2, lat_n4;
int nd2, nd4;
double *Fock_2c_buffer;
double time1, time2;
size_t result;
FILE *integrals_2c; 
char xx[4], yy[30] = "scf_integrals_2C_no_sym.";

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file);
  dimp_spin = dimp * job->spin_dim;
  AllocateDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2c_buffer,&dimp_spin);
	  
  sprintf(xx, "%d", job->taskid);

  if (job->taskid > 0 || job->numtasks == 1) {

  if ((integrals_2c = fopen(strcat(yy,xx), "rb")) == NULL) {
  printf("Failed to open file %s in coulomb_matrix_crystal_read_integrals_no_sym\n",yy);
  fprintf(file.out,"Failed to open file %s in coulomb_matrix_crystal_read_integrals_no_sym\n",yy);
  MPI_Finalize();
  exit(1);
 }

  while (!feof(integrals_2c)) {

  QUAD_TRAN Quad;
  INTEGRAL_LIST integral_list;
  result = fread(&Quad.tot,sizeof(int),1,integrals_2c);
  result = fread(&integral_list.num,sizeof(int),1,integrals_2c);
  if (feof(integrals_2c)) break;
  allocate_QUAD_TRAN(&Quad,job,file);
  allocate_integral_list(&integral_list, integral_list.num, job, file);
  read_unpack_integrals_crystal_coulomb_ijkl(&integral_list, &Quad, integrals_2c, job, file);
  atm_n1 = Quad.cell1[0];
  atm_n2 = Quad.cell2[0];
  atm_n3 = Quad.cell3[0];
  atm_n4 = Quad.cell4[0];
  lat_n2 = Quad.latt2[0];
  lat_n4 = Quad.latt4[0];
  nd2 = atoms->bfnnumb_sh[atm_n2];
  nd4 = atoms->bfnnumb_sh[atm_n4];
  F_ptr = pair_p->Off[pair_p->Ptr[lat_n2 * dim2 + atm_n1 * dim1 + atm_n2]];
  D_ptr = pair_p->off[pair_p->ptr[lat_n4 * dim2 + atm_n3 * dim1 + atm_n4]];
  for (m = 0; m < integral_list.num; m++) {
  Fock_temp_offset      =         nd2 * integral_list.i[m] + integral_list.j[m];
  Density_matrix_offset = D_ptr + nd4 * integral_list.k[m] + integral_list.l[m];
  Fock_2c_buffer[F_ptr + Fock_temp_offset] += integral_list.value[m] * F[Density_matrix_offset];
 }
  free_QUAD_TRAN(&Quad,job);
  free_integral_list(&integral_list,job);

 } //while (!feof

  fclose(integrals_2c);

 }

  MPI_Allreduce(Fock_2c_buffer,Fock_2c,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&Fock_2c_buffer,&dimp_spin,job);

}

void fock_matrix_crystal_read_exchange_integrals_no_sym(double *Fock_2e, double *F, PAIR_TRAN *pair_p, ATOM *atoms, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * This routine reads crystal exchange integrals from disk and contracts them - no sym    *
  // * not used and not tested                                                                *
  // ******************************************************************************************

int m, gl0;
int dimf, dimp, dimp_spin;
int Fock_temp_offset, Density_matrix_offset;
int F_ptr, D_ptr;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
int atm_n1, atm_n2, atm_n3, atm_n4;
int lat_n2, lat_n3, lat_n4;
int nd3, nd4;
double *Fock_2e_buffer;
double time1, time2;
size_t result;
FILE *integrals_2e; 
char xx[4], yy[24] = "scf_integrals_2E.";

  sprintf(xx, "%d", job->taskid);

  if (job->taskid > 0 || job->numtasks == 1) {

  if ((integrals_2e = fopen(strcat(yy,xx), "rb")) == NULL) {
  printf("Failed to open file %s in exchange_matrix_crystal_read_integrals\n",yy);
  fprintf(file.out,"Failed to open file %s in exchange_matrix_crystal_read_integrals\n",yy);
  MPI_Finalize();
  exit(1);
 }

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file);
  dimp_spin = dimp * job->spin_dim;
  AllocateDoubleArray(&Fock_2e_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2e_buffer,&dimp_spin);
	  
  while (!feof(integrals_2e)) {

  QUAD_TRAN Quad;
  INTEGRAL_LIST integral_list;
  result = fread(&Quad.tot,sizeof(int),1,integrals_2e);
  result = fread(&integral_list.num,sizeof(int),1,integrals_2e);
  if (feof(integrals_2e)) break;
  allocate_QUAD_TRAN(&Quad,job,file);
  allocate_integral_list(&integral_list, integral_list.num, job, file);
  read_unpack_integrals_crystal_exchange_ijkl(&integral_list, &Quad, integrals_2e, job, file);
  atm_n1 = Quad.cell1[0];
  atm_n2 = Quad.cell2[0];
  atm_n3 = Quad.cell3[0];
  atm_n4 = Quad.cell4[0];
  lat_n2 = Quad.latt2[0];
  lat_n3 = Quad.latt3[0];
  lat_n4 = Quad.latt4[0];
  nd3 = atoms->bfnnumb_sh[atm_n3];
  nd4 = atoms->bfnnumb_sh[atm_n4];
  F_ptr = pair_p->Off[pair_p->Ptr[lat_n3 * dim2 + atm_n1 * dim1 + atm_n3]];
  gl0 = R_tables->diffvec[lat_n4 * R_tables->margin_vector + lat_n2];
  D_ptr = pair_p->off[pair_p->ptr[gl0 * dim2 + atm_n2 * dim1 + atm_n4]];
  for (m = 0; m < integral_list.num; m++) {
  Fock_temp_offset      =         nd3 * integral_list.i[m] + integral_list.k[m];
  Density_matrix_offset = D_ptr + nd4 * integral_list.j[m] + integral_list.l[m];
  Fock_2e_buffer[F_ptr + Fock_temp_offset] -= integral_list.value[m] * F[Density_matrix_offset] / two;
 }
  free_QUAD_TRAN(&Quad,job);
  free_integral_list(&integral_list,job);

 } //while (!feof

  fclose(integrals_2e);

 }

  MPI_Allreduce(Fock_2e_buffer,Fock_2e,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&Fock_2e_buffer,&dimp_spin,job);

}

/*
void fock_matrix_molecule_read_integrals_no_sym(double *Fock_2c, double *Fock_2e, double *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * This routine reads molecule integrals from disk and contracts them                     *
  // ******************************************************************************************

int m;
int Fock_temp_offset, Density_matrix_offset;
int F_ptr, D_ptr;
int atm_n1, atm_n2, atm_n3, atm_n4;
int nd2, nd3, nd4;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dimf, dimp, dimp_spin;
double *Fock_2e_buffer, *Fock_2c_buffer;
double time1, time2;
FILE *integrals_2m; 
char xx[4], yy[30] = "scf_integrals_2M_no_sym.";

  sh_array_dimensions(&dimp, &dimf, pair_p, atoms, job, file);
  dimp_spin = dimp * job->spin_dim;
  AllocateDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2c_buffer,&dimp_spin);
  AllocateDoubleArray(&Fock_2e_buffer,&dimp_spin,job);
  ResetDoubleArray(Fock_2e_buffer,&dimp_spin);
	  
  sprintf(xx, "%d", job->taskid);

  if (job->taskid > 0 || job->numtasks == 1) { 

  if ((integrals_2m = fopen(strcat(yy,xx), "rb")) == NULL) {
  printf("Failed to open file %s in coulomb_exchange_matrix_molecule_read_integrals\n",yy);
  fprintf(file.out,"Failed to open file %s in coulomb_exchange_matrix_molecule_read_integrals\n",yy);
  MPI_Finalize();
  exit(1);
 }

  while (!feof(integrals_2m)) {

  size_t result;
  QUAD_TRAN Quad;
  INTEGRAL_LIST integral_list;
  result = fread(&Quad.tot,sizeof(int),1,integrals_2m);
  result = fread(&integral_list.num,sizeof(int),1,integrals_2m);
  if (feof(integrals_2m)) break;
  allocate_QUAD_TRAN(&Quad,job,file);
  allocate_integral_list(&integral_list, integral_list.num, job, file);
  read_unpack_molecule_2e_integrals(&integral_list, &Quad, integrals_2m, job, file);
  atm_n1 = Quad.cell1[0];
  atm_n2 = Quad.cell2[0];
  atm_n3 = Quad.cell3[0];
  atm_n4 = Quad.cell4[0];
  nd2 = atoms->bfnnumb_sh[atm_n2];
  nd3 = atoms->bfnnumb_sh[atm_n3];
  nd4 = atoms->bfnnumb_sh[atm_n4];
  if (pair_p->uniq[atm_n1 * dim1 + atm_n2] == -1) {
  F_ptr = pair_p->Off[pair_p->Ptr[atm_n1 * dim1 + atm_n2]];
  D_ptr = pair_p->off[pair_p->ptr[atm_n3 * dim1 + atm_n4]];
  for (m = 0; m < integral_list.num; m++) {
    Fock_temp_offset      =         nd2 * integral_list.i[m] + integral_list.j[m];
    Density_matrix_offset = D_ptr + nd4 * integral_list.k[m] + integral_list.l[m];
    Fock_2c_buffer[F_ptr + Fock_temp_offset] += integral_list.value[m] * F[Density_matrix_offset];
   }
  }
  if (pair_p->uniq[atm_n1 * dim1 + atm_n3] == -1) {
  F_ptr = pair_p->Off[pair_p->Ptr[atm_n1 * dim1 + atm_n3]];
  D_ptr = pair_p->off[pair_p->ptr[atm_n2 * dim1 + atm_n4]];
  for (m = 0; m < integral_list.num; m++) {
    Fock_temp_offset      =         nd3 * integral_list.i[m] + integral_list.k[m];
    Density_matrix_offset = D_ptr + nd4 * integral_list.j[m] + integral_list.l[m];
    Fock_2e_buffer[F_ptr + Fock_temp_offset] -= integral_list.value[m] * F[Density_matrix_offset] / two;
   }
  }
  free_QUAD_TRAN(&Quad,job);
  free_integral_list(&integral_list,job);

 } //while (!feof

  fclose(integrals_2m);
}

  MPI_Allreduce(Fock_2c_buffer,Fock_2c,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&Fock_2c_buffer,&dimp_spin,job);
  MPI_Allreduce(Fock_2e_buffer,Fock_2e,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&Fock_2e_buffer,&dimp_spin,job);

}
*/

void contract_integrals_crystal_coulomb_ijkl(double *Fock_2c_buffer, INTEGRAL_LIST *integral_list_coulomb, double *F, PAIR_TRAN *pair_p, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{
int i, m, s, s1, op, pm, dimc, F_ptr, D_ptr;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
int ip0, jp0, kp0, lp0, gj0, gk0, gl0, ip, jp, kp, lp, gj, gl, gl0_gk0;
int nd1, nd2, nd3, nd4;
int Fock_temp_offset, Density_matrix_offset, coulomb_ints_num; 
int *p_coulomb_ints_i, *p_coulomb_ints_j, *p_coulomb_ints_k, *p_coulomb_ints_l;
double *p_coulomb_ints_value, *Fock_2c_temp;

      for (i = 0; i < quad->tot; i++) {

        pm = quad->p[i];
        op = quad->k[i];

        for (s = 0; s < job->spin_dim; s++) {

          switch (pm) {

          case 0:

          ip0 = quad->cell1[0];
          jp0 = quad->cell2[0];
          kp0 = quad->cell3[0];
          lp0 = quad->cell4[0];

          gl0 = quad->latt4[0];

          p_coulomb_ints_i = integral_list_coulomb->i;
          p_coulomb_ints_j = integral_list_coulomb->j;
          p_coulomb_ints_k = integral_list_coulomb->k;
          p_coulomb_ints_l = integral_list_coulomb->l;

          break;
         
          case 1:

          ip0 = quad->cell1[0];
          jp0 = quad->cell2[0];
          kp0 = quad->cell4[0];
          lp0 = quad->cell3[0];

          gl0 = R_tables->diffvec[quad->latt4[0]];

          p_coulomb_ints_i = integral_list_coulomb->i;
          p_coulomb_ints_j = integral_list_coulomb->j;
          p_coulomb_ints_k = integral_list_coulomb->l;
          p_coulomb_ints_l = integral_list_coulomb->k;

          break;
         
          case 2:

          ip0 = quad->cell2[0];
          jp0 = quad->cell1[0];
          kp0 = quad->cell3[0];
          lp0 = quad->cell4[0];

          gl0 = quad->latt4[0];

          p_coulomb_ints_i = integral_list_coulomb->j;
          p_coulomb_ints_j = integral_list_coulomb->i;
          p_coulomb_ints_k = integral_list_coulomb->k;
          p_coulomb_ints_l = integral_list_coulomb->l;

          break;
         
          case 3:

          ip0 = quad->cell2[0];
          jp0 = quad->cell1[0];
          kp0 = quad->cell4[0];
          lp0 = quad->cell3[0];

          gl0 = R_tables->diffvec[quad->latt4[0]];

          p_coulomb_ints_i = integral_list_coulomb->j;
          p_coulomb_ints_j = integral_list_coulomb->i;
          p_coulomb_ints_k = integral_list_coulomb->l;
          p_coulomb_ints_l = integral_list_coulomb->k;

          break;
         
          case 4:

          ip0 = quad->cell3[0];
          jp0 = quad->cell4[0];
          kp0 = quad->cell1[0];
          lp0 = quad->cell2[0];

          gl0 = quad->latt2[0];

          p_coulomb_ints_i = integral_list_coulomb->k;
          p_coulomb_ints_j = integral_list_coulomb->l;
          p_coulomb_ints_k = integral_list_coulomb->i;
          p_coulomb_ints_l = integral_list_coulomb->j;

          break;
         
          case 5:

          ip0 = quad->cell3[0];
          jp0 = quad->cell4[0];
          kp0 = quad->cell2[0];
          lp0 = quad->cell1[0];

          gl0 = R_tables->diffvec[quad->latt2[0]];

          p_coulomb_ints_i = integral_list_coulomb->k;
          p_coulomb_ints_j = integral_list_coulomb->l;
          p_coulomb_ints_k = integral_list_coulomb->j;
          p_coulomb_ints_l = integral_list_coulomb->i;

          break;
         
          case 6:

          ip0 = quad->cell4[0];
          jp0 = quad->cell3[0];
          kp0 = quad->cell1[0];
          lp0 = quad->cell2[0];

          gl0 = quad->latt2[0];

          p_coulomb_ints_i = integral_list_coulomb->l;
          p_coulomb_ints_j = integral_list_coulomb->k;
          p_coulomb_ints_k = integral_list_coulomb->i;
          p_coulomb_ints_l = integral_list_coulomb->j;

          break;
         
          case 7:

          ip0 = quad->cell4[0];
          jp0 = quad->cell3[0];
          kp0 = quad->cell2[0];
          lp0 = quad->cell1[0];

          gl0 = R_tables->diffvec[quad->latt2[0]];

          p_coulomb_ints_i = integral_list_coulomb->l;
          p_coulomb_ints_j = integral_list_coulomb->k;
          p_coulomb_ints_k = integral_list_coulomb->j;
          p_coulomb_ints_l = integral_list_coulomb->i;

          break;
         
         } // close switch (pm)

          ip = quad->cell1[i];
          jp = quad->cell2[i];
          kp = quad->cell3[i];
          lp = quad->cell4[i];

          gj = quad->latt2[i];
          gl = quad->latt4[i];

          nd1 = atoms->bfnnumb_sh[ip];
          nd2 = atoms->bfnnumb_sh[jp];
          nd4 = atoms->bfnnumb_sh[lp0];

          dimc = nd1 * nd2;

          AllocateDoubleArray(&Fock_2c_temp,&dimc,job);
          ResetDoubleArray(Fock_2c_temp,&dimc);

          coulomb_ints_num     = integral_list_coulomb->num;
          p_coulomb_ints_value = integral_list_coulomb->value;

          if (pair_p->uniq[gj * dim2 + ip * dim1 + jp] == -1) {

          F_ptr = s * job->dimp + pair_p->Off[pair_p->Ptr[gj  * dim2 + ip  * dim1 + jp]];
          D_ptr = pair_p->off[pair_p->ptr[gl0 * dim2 + kp0 * dim1 + lp0]];
          for (m = 0; m < coulomb_ints_num; m++) {
          Fock_temp_offset      =         nd2 * *p_coulomb_ints_i + *p_coulomb_ints_j;
          Density_matrix_offset = D_ptr + nd4 * *p_coulomb_ints_k + *p_coulomb_ints_l;
          Fock_2c_temp[Fock_temp_offset] += *p_coulomb_ints_value * F[Density_matrix_offset];
          p_coulomb_ints_value++;
          p_coulomb_ints_i++;
          p_coulomb_ints_j++;
          p_coulomb_ints_k++;
          p_coulomb_ints_l++;
         }
          //for (s1 = 0; s1 < job->spin_dim; s1++) {
          //D_ptr = pair_p->off[pair_p->ptr[kp0 * dim1 + lp0]];
          //for (m = 0; m < coulomb_ints_num; m++) {
          //Fock_temp_offset      =         nd2 * *p_coulomb_ints_i + *p_coulomb_ints_j;
          //Density_matrix_offset = D_ptr + nd4 * *p_coulomb_ints_k + *p_coulomb_ints_l;
          //Fock_2c_temp[Fock_temp_offset] += *p_coulomb_ints_value * F[Density_matrix_offset];
          //p_coulomb_ints_value++;
          //p_coulomb_ints_i++;
          //p_coulomb_ints_j++;
          //p_coulomb_ints_k++;
          //p_coulomb_ints_l++;
         //}
          //p_coulomb_ints_value -= coulomb_ints_num;
          //p_coulomb_ints_i     -= coulomb_ints_num;
          //p_coulomb_ints_j     -= coulomb_ints_num;
          //p_coulomb_ints_k     -= coulomb_ints_num;
          //p_coulomb_ints_l     -= coulomb_ints_num;
         //}
          rotate_sum_block(Fock_2c_temp,&Fock_2c_buffer[F_ptr],ip,jp,op,atoms,shells,symmetry,job,file);
         }
          DestroyDoubleArray(&Fock_2c_temp,&dimc,job);
 
         } // close loop on s

        } // close loop on i

}

void contract_integrals_crystal_exchange_ijkl(double *Fock_2e_buffer, INTEGRAL_LIST *integral_list_exchange, double *F, PAIR_TRAN *pair_p, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

int i, m, s, op, pm, dime, F_ptr, D_ptr;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
int jp0,kp0,lp0,gj0,gk0,gl0,ip,jp,kp,lp,gj,gk,gl,gl0_gj0;
int nd1, nd2, nd3, nd4;
int Fock_temp_offset, Density_matrix_offset, exchange_ints_num; 
int *p_exchange_ints_i, *p_exchange_ints_j, *p_exchange_ints_k, *p_exchange_ints_l;
double *p_exchange_ints_value, *Fock_2e_temp;

      for (i = 0; i < quad->tot; i++) {

        pm = quad->p[i];
        op = quad->k[i];

        for (s = 0; s < job->spin_dim; s++) {

          switch (pm) {

          case 0:

          jp0 = quad->cell2[0];
          kp0 = quad->cell3[0];
          lp0 = quad->cell4[0];

          gj0 = quad->latt2[0];
          gl0 = quad->latt4[0];

          p_exchange_ints_i = integral_list_exchange->i;
          p_exchange_ints_j = integral_list_exchange->j;
          p_exchange_ints_k = integral_list_exchange->k;
          p_exchange_ints_l = integral_list_exchange->l;

          break;
         
          case 1:

          jp0 = quad->cell2[0];
          kp0 = quad->cell4[0];
          lp0 = quad->cell3[0];

          gj0 = quad->latt2[0];
          gl0 = quad->latt3[0];

          p_exchange_ints_i = integral_list_exchange->i;
          p_exchange_ints_j = integral_list_exchange->j;
          p_exchange_ints_k = integral_list_exchange->l;
          p_exchange_ints_l = integral_list_exchange->k;

          break;
         
          case 2:

          jp0 = quad->cell1[0];
          kp0 = quad->cell3[0];
          lp0 = quad->cell4[0];

          gj0 = quad->latt1[0];
          gl0 = quad->latt4[0];

          p_exchange_ints_i = integral_list_exchange->j;
          p_exchange_ints_j = integral_list_exchange->i;
          p_exchange_ints_k = integral_list_exchange->k;
          p_exchange_ints_l = integral_list_exchange->l;

          break;
         
          case 3:

          jp0 = quad->cell1[0];
          kp0 = quad->cell4[0];
          lp0 = quad->cell3[0];

          gj0 = quad->latt1[0];
          gl0 = quad->latt3[0];

          p_exchange_ints_i = integral_list_exchange->j;
          p_exchange_ints_j = integral_list_exchange->i;
          p_exchange_ints_k = integral_list_exchange->l;
          p_exchange_ints_l = integral_list_exchange->k;

          break;
         
          case 4:

          jp0 = quad->cell4[0];
          kp0 = quad->cell1[0];
          lp0 = quad->cell2[0];

          gj0 = quad->latt4[0];
          gl0 = quad->latt2[0];

          p_exchange_ints_i = integral_list_exchange->k;
          p_exchange_ints_j = integral_list_exchange->l;
          p_exchange_ints_k = integral_list_exchange->i;
          p_exchange_ints_l = integral_list_exchange->j;

          break;
         
          case 5:

          jp0 = quad->cell4[0];
          kp0 = quad->cell2[0];
          lp0 = quad->cell1[0];

          gj0 = quad->latt4[0];
          gl0 = quad->latt1[0];

          p_exchange_ints_i = integral_list_exchange->k;
          p_exchange_ints_j = integral_list_exchange->l;
          p_exchange_ints_k = integral_list_exchange->j;
          p_exchange_ints_l = integral_list_exchange->i;

          break;
         
          case 6:

          jp0 = quad->cell3[0];
          kp0 = quad->cell1[0];
          lp0 = quad->cell2[0];

          gj0 = quad->latt3[0];
          gl0 = quad->latt2[0];

          p_exchange_ints_i = integral_list_exchange->l;
          p_exchange_ints_j = integral_list_exchange->k;
          p_exchange_ints_k = integral_list_exchange->i;
          p_exchange_ints_l = integral_list_exchange->j;

          break;
         
          case 7:

          jp0 = quad->cell3[0];
          kp0 = quad->cell2[0];
          lp0 = quad->cell1[0];

          gj0 = quad->latt3[0];
          gl0 = quad->latt1[0];

          p_exchange_ints_i = integral_list_exchange->l;
          p_exchange_ints_j = integral_list_exchange->k;
          p_exchange_ints_k = integral_list_exchange->j;
          p_exchange_ints_l = integral_list_exchange->i;

          break;
         
         } // close switch (pm)

          ip = quad->cell1[i];
          jp = quad->cell2[i];
          kp = quad->cell3[i];
          lp = quad->cell4[i];

          gj = quad->latt2[i];
          gk = quad->latt3[i];
          gl = quad->latt4[i];

          gl0_gj0 = R_tables->diffvec[gl0 * R_tables->margin_vector + gj0];
 
          nd1 = atoms->bfnnumb_sh[ip];
          nd3 = atoms->bfnnumb_sh[kp];
          nd4 = atoms->bfnnumb_sh[lp0];
          dime = nd1 * nd3;
          AllocateDoubleArray(&Fock_2e_temp,&dime,job);
          ResetDoubleArray(Fock_2e_temp,&dime);
          exchange_ints_num     = integral_list_exchange->num;
          p_exchange_ints_value = integral_list_exchange->value;

          if (job->scf_exchange == 1 && job->xc_hfx == 1 && pair_p->uniq[gk * dim2 + ip * dim1 + kp] == -1) {
          F_ptr = s * job->dimp + pair_p->Off[pair_p->Ptr[gk * dim2 + ip * dim1 + kp]];
          D_ptr = s * job->dimf + pair_p->off[pair_p->ptr[gl0_gj0  * dim2 + jp0  * dim1 + lp0]];
          for (m = 0; m < exchange_ints_num; m++) {
          Fock_temp_offset      =         nd3 * *p_exchange_ints_i + *p_exchange_ints_k;
          Density_matrix_offset = D_ptr + nd4 * *p_exchange_ints_j + *p_exchange_ints_l;
          Fock_2e_temp[Fock_temp_offset] -= *p_exchange_ints_value * F[Density_matrix_offset] / two;
          p_exchange_ints_value++;
          p_exchange_ints_i++;
          p_exchange_ints_j++;
          p_exchange_ints_k++;
          p_exchange_ints_l++;
         }
          rotate_sum_block(Fock_2e_temp,&Fock_2e_buffer[F_ptr],ip,kp,op,atoms,shells,symmetry,job,file);
         }

          DestroyDoubleArray(&Fock_2e_temp,&dime,job);

        } // close loop on s

       } // close loop on i

}

/*
void contract_integrals_molecule_ijkl(double *Fock_2c_buffer, double *Fock_2e_buffer, INTEGRAL_LIST *integral_list_molecule, double *F, PAIR_TRAN *pair_p, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{
int i, m, pm, op, s, s1;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1;
int ip, jp, kp, lp;
int ip0,jp0,kp0,lp0;
int nd1, nd2, nd3, nd4;
int Fock_temp_offset, Density_matrix_offset; 
int F_ptr, D_ptr;
int dimc, dime;
int molecule_ints_num;
int *p_coulomb_ints_i, *p_coulomb_ints_j, *p_coulomb_ints_k, *p_coulomb_ints_l;
int *p_exchange_ints_i, *p_exchange_ints_j, *p_exchange_ints_k, *p_exchange_ints_l;
double *p_coulomb_ints_value, *p_exchange_ints_value;
double *p_molecule_ints_value;
double *Fock_2c_temp, *Fock_2e_temp;

        for (i = 0; i < quad->tot; i++) {

          pm = quad->p[i];
          op = quad->k[i];

        for (s = 0; s < job->spin_dim; s++) {

          switch (pm) {

          case 0:

          jp0 = quad->cell2[0];
          kp0 = quad->cell3[0];
          lp0 = quad->cell4[0];

          p_coulomb_ints_i = integral_list_molecule->i;
          p_coulomb_ints_j = integral_list_molecule->j;
          p_coulomb_ints_k = integral_list_molecule->k;
          p_coulomb_ints_l = integral_list_molecule->l;

          p_exchange_ints_i = integral_list_molecule->i;
          p_exchange_ints_j = integral_list_molecule->j;
          p_exchange_ints_k = integral_list_molecule->k;
          p_exchange_ints_l = integral_list_molecule->l;

          break;
         
          case 1:

          jp0 = quad->cell2[0];
          kp0 = quad->cell4[0];
          lp0 = quad->cell3[0];

          p_coulomb_ints_i = integral_list_molecule->i;
          p_coulomb_ints_j = integral_list_molecule->j;
          p_coulomb_ints_k = integral_list_molecule->l;
          p_coulomb_ints_l = integral_list_molecule->k;

          p_exchange_ints_i = integral_list_molecule->i;
          p_exchange_ints_j = integral_list_molecule->j;
          p_exchange_ints_k = integral_list_molecule->l;
          p_exchange_ints_l = integral_list_molecule->k;

          break;
         
          case 2:

          jp0 = quad->cell1[0];
          kp0 = quad->cell3[0];
          lp0 = quad->cell4[0];

          p_coulomb_ints_i = integral_list_molecule->j;
          p_coulomb_ints_j = integral_list_molecule->i;
          p_coulomb_ints_k = integral_list_molecule->k;
          p_coulomb_ints_l = integral_list_molecule->l;

          p_exchange_ints_i = integral_list_molecule->j;
          p_exchange_ints_j = integral_list_molecule->i;
          p_exchange_ints_k = integral_list_molecule->k;
          p_exchange_ints_l = integral_list_molecule->l;

          break;
         
          case 3:

          jp0 = quad->cell1[0];
          kp0 = quad->cell4[0];
          lp0 = quad->cell3[0];

          p_coulomb_ints_i = integral_list_molecule->j;
          p_coulomb_ints_j = integral_list_molecule->i;
          p_coulomb_ints_k = integral_list_molecule->l;
          p_coulomb_ints_l = integral_list_molecule->k;

          p_exchange_ints_i = integral_list_molecule->j;
          p_exchange_ints_j = integral_list_molecule->i;
          p_exchange_ints_k = integral_list_molecule->l;
          p_exchange_ints_l = integral_list_molecule->k;

          break;
         
          case 4:

          jp0 = quad->cell4[0];
          kp0 = quad->cell1[0];
          lp0 = quad->cell2[0];

          p_coulomb_ints_i = integral_list_molecule->k;
          p_coulomb_ints_j = integral_list_molecule->l;
          p_coulomb_ints_k = integral_list_molecule->i;
          p_coulomb_ints_l = integral_list_molecule->j;

          p_exchange_ints_i = integral_list_molecule->k;
          p_exchange_ints_j = integral_list_molecule->l;
          p_exchange_ints_k = integral_list_molecule->i;
          p_exchange_ints_l = integral_list_molecule->j;

          break;
         
          case 5:

          jp0 = quad->cell4[0];
          kp0 = quad->cell2[0];
          lp0 = quad->cell1[0];

          p_coulomb_ints_i = integral_list_molecule->k;
          p_coulomb_ints_j = integral_list_molecule->l;
          p_coulomb_ints_k = integral_list_molecule->j;
          p_coulomb_ints_l = integral_list_molecule->i;

          p_exchange_ints_i = integral_list_molecule->k;
          p_exchange_ints_j = integral_list_molecule->l;
          p_exchange_ints_k = integral_list_molecule->j;
          p_exchange_ints_l = integral_list_molecule->i;

          break;
         
          case 6:

          jp0 = quad->cell3[0];
          kp0 = quad->cell1[0];
          lp0 = quad->cell2[0];

          p_coulomb_ints_i = integral_list_molecule->l;
          p_coulomb_ints_j = integral_list_molecule->k;
          p_coulomb_ints_k = integral_list_molecule->i;
          p_coulomb_ints_l = integral_list_molecule->j;

          p_exchange_ints_i = integral_list_molecule->l;
          p_exchange_ints_j = integral_list_molecule->k;
          p_exchange_ints_k = integral_list_molecule->i;
          p_exchange_ints_l = integral_list_molecule->j;

          break;
         
          case 7:

          jp0 = quad->cell3[0];
          kp0 = quad->cell2[0];
          lp0 = quad->cell1[0];

          p_coulomb_ints_i = integral_list_molecule->l;
          p_coulomb_ints_j = integral_list_molecule->k;
          p_coulomb_ints_k = integral_list_molecule->j;
          p_coulomb_ints_l = integral_list_molecule->i;

          p_exchange_ints_i = integral_list_molecule->l;
          p_exchange_ints_j = integral_list_molecule->k;
          p_exchange_ints_k = integral_list_molecule->j;
          p_exchange_ints_l = integral_list_molecule->i;

          break;
         
         } // close switch (pm)

          ip = quad->cell1[i];
          jp = quad->cell2[i];
          kp = quad->cell3[i];
          lp = quad->cell4[i];

          nd1 = atoms->bfnnumb_sh[ip];
          nd2 = atoms->bfnnumb_sh[jp];
          nd3 = atoms->bfnnumb_sh[kp];
          nd4 = atoms->bfnnumb_sh[lp0];

//fprintf(file.out,"%3d %3d   %3d %3d %3d %3d   %3d %3d %3d %3d\n",pm,op,ip0,jp0,kp0,lp0,ip,jp,kp0,lp0);

          dimc = nd1 * nd2;
          dime = nd1 * nd3;

          AllocateDoubleArray(&Fock_2c_temp,&dimc,job);
          ResetDoubleArray(Fock_2c_temp,&dimc);
          AllocateDoubleArray(&Fock_2e_temp,&dime,job);
          ResetDoubleArray(Fock_2e_temp,&dime);

          molecule_ints_num     = integral_list_molecule->num;
          p_coulomb_ints_value  = integral_list_molecule->value;
          p_exchange_ints_value = integral_list_molecule->value;

          if (job->scf_coulomb == 1 && pair_p->uniq[ip * dim1 + jp] == -1) {
          //HERE 22/05/20 F_ptr = pair_p->Off[pair_p->Ptr[ip  * dim1 + jp]];
          //D_ptr = pair_p->off[pair_p->ptr[kp0 * dim1 + lp0]];
          F_ptr = s * job->dimp + pair_p->Off[pair_p->Ptr[ip  * dim1 + jp]];
          for (s1 = 0; s1 < job->spin_dim; s1++) {
          D_ptr = s1 * job->dimf + pair_p->off[pair_p->ptr[kp0 * dim1 + lp0]];
          for (m = 0; m < molecule_ints_num; m++) {
          Fock_temp_offset      =         nd2 * *p_coulomb_ints_i + *p_coulomb_ints_j;
          Density_matrix_offset = D_ptr + nd4 * *p_coulomb_ints_k + *p_coulomb_ints_l;
          Fock_2c_temp[Fock_temp_offset] += *p_coulomb_ints_value * F[Density_matrix_offset];
//fprintf(file.out,"INT %3d %3d %3d %3d  %9.2e %9.2e %9.2e\n",*p_coulomb_ints_i,*p_coulomb_ints_j,*p_coulomb_ints_k,*p_coulomb_ints_l,\
*p_coulomb_ints_value,F[Density_matrix_offset],*p_coulomb_ints_value * F[Density_matrix_offset]);
          p_coulomb_ints_value++;
          p_coulomb_ints_i++;
          p_coulomb_ints_j++;
          p_coulomb_ints_k++;
          p_coulomb_ints_l++;
         }
          p_coulomb_ints_value -= molecule_ints_num;
          p_coulomb_ints_i     -= molecule_ints_num;
          p_coulomb_ints_j     -= molecule_ints_num;
          p_coulomb_ints_k     -= molecule_ints_num;
          p_coulomb_ints_l     -= molecule_ints_num;
         }
          rotate_sum_block(Fock_2c_temp,&Fock_2c_buffer[F_ptr],ip,jp,op,atoms,shells,symmetry,job,file);
         }
          DestroyDoubleArray(&Fock_2c_temp,&dimc,job);

          if (job->scf_exchange == 1 && job->xc_hfx == 1 && pair_p->uniq[ip * dim1 + kp] == -1) {
          F_ptr = s * job->dimp + pair_p->Off[pair_p->Ptr[ip  * dim1 + kp]];
          D_ptr = s * job->dimf + pair_p->off[pair_p->ptr[jp0 * dim1 + lp0]];
          //HERE 22/05/20 F_ptr = s * job->dimp + pair_p->Off[pair_p->Ptr[ip  * dim1 + kp]];
          //D_ptr = s * job->dimf + pair_p->off[pair_p->ptr[jp0 * dim1 + lp0]];
          for (m = 0; m < molecule_ints_num; m++) {
          Fock_temp_offset      =         nd3 * *p_exchange_ints_i + *p_exchange_ints_k;
          Density_matrix_offset = D_ptr + nd4 * *p_exchange_ints_j + *p_exchange_ints_l;
          Fock_2e_temp[Fock_temp_offset] -= *p_exchange_ints_value * F[Density_matrix_offset] / two;
//fprintf(file.out,"INT %3d %3d %3d %3d  %9.2e %9.2e %9.2e\n",*p_exchange_ints_i,*p_exchange_ints_j,*p_exchange_ints_k,*p_exchange_ints_l,\
*p_exchange_ints_value,F[Density_matrix_offset],*p_exchange_ints_value * F[Density_matrix_offset]);
          p_exchange_ints_value++;
          p_exchange_ints_i++;
          p_exchange_ints_j++;
          p_exchange_ints_k++;
          p_exchange_ints_l++;
         }
          rotate_sum_block(Fock_2e_temp,&Fock_2e_buffer[F_ptr],ip,kp,op,atoms,shells,symmetry,job,file);
         }

          DestroyDoubleArray(&Fock_2e_temp,&dime,job);

        } // close loop on s

        } // close loop on i

}

void shell_screen_molecule_compute_integrals(double *S1, PAIR_TRAN *pair_p, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int i, j, k, p, q;
int ip, jp, gi, gj;
int nd1, nd2, nd3, nd4, nd12, nd34;
int dim;
int count;

double *S0;

  dim = 0;
  AllocateDoubleArray(&S0,&job->dimp,job);
  for (i = 0; i < job->dimp; i++) S0[i] = k_zero;
  for (i = 0; i < job->dimf; i++) S1[i] = k_zero;
  for (p = 0; p < pair_p->nump; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    gj = 0;
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    integrals_molecule_screen(&S0[dim],ip,jp,R,atoms,shells,gaussians,symmetry,crystal,job,file);
    dim += nd34;
   }
    expand_screening_integral_matrix(S0,S1,pair_p,atoms,shells,symmetry,job,file);
    DestroyDoubleArray(&S0,&job->dimp,job);

    if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"screened integrals matrix %3d %3d %3d %3d\n",pair_p->nump,dim,job->dimp,job->dimf);
    count = 0;
    for (p = 0; p < pair_p->nump; p++) {
      q = pair_p->posn[p];
      fprintf(file.out,"pair %d [%3d %3d] [%3d]\n",p,pair_p->cell1[q],pair_p->cell2[q],pair_p->latt2[q]);
      for (k = 0; k < pair_p->numb[p]; k++) {
        for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
          for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
            fprintf(file.out,"%10.2e ",S1[count]);
            count++;
           }
          fprintf(file.out,"\n");
         }
        fprintf(file.out,"\n");
       }
      fprintf(file.out,"\n\n");
     }
    }

}
*/

void shell_screen_direct(int *start_index, double *S1, double *F, PAIR_TRAN *pair_p, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int index_i, index_j, index_k, index_l;
int shelposi, shelposj, shelposk, shelposl;
int bfposi1, bfposj1, bfposk1, bfposl1;
int sheli1, shelj1, shelk1, shell1;
int shell_count, flag;
int i, j, k, l, ip, jp, kp, lp, gi, gj, gk, gl;
int i1, m, pm, op, p, s, s1;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1, dimi;
int dimc, dime;
int ip0, jp0, kp0, lp0, gl0;
int ip1, jp1, kp1, lp1;
int nd1, nd2, nd3, nd4;
int nd11, nd12, nd13, nd14;
int Fock_temp_offset, Density_matrix_offset; 
int F_ptr, D_ptr;
int molecule_ints_num;
int *p_coulomb_ints_i, *p_coulomb_ints_j, *p_coulomb_ints_k, *p_coulomb_ints_l;
int *p_exchange_ints_i, *p_exchange_ints_j, *p_exchange_ints_k, *p_exchange_ints_l;
double *Fock_2c, *Fock_2e, *Fock_2c_temp, *Fock_2e_temp;
double *p_coulomb_ints_value, *p_exchange_ints_value;
double *p_molecule_ints_value;
double largest_integral, testint, integral_rejection_threshold_sqrd;
INTEGRAL_LIST integral_list_molecule;

  integral_rejection_threshold_sqrd = integral_rejection_threshold * integral_rejection_threshold;
  //integral_rejection_threshold_sqrd = 1e-16;

  ip1 = quad->cell1[0];
  jp1 = quad->cell2[0];
  kp1 = quad->cell3[0];
  lp1 = quad->cell4[0];

  nd11 = atoms->bfnnumb_sh[ip1];
  nd12 = atoms->bfnnumb_sh[jp1];
  nd13 = atoms->bfnnumb_sh[kp1];
  nd14 = atoms->bfnnumb_sh[lp1];

  shell_count = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip1];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip1]; index_i++) {
    sheli1 = shells->type_sh[index_i];
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp1];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp1]; index_j++) {
      shelj1 = shells->type_sh[index_j];
      bfposk1  = 0;
      shelposk = atoms->shelposn_sh[kp1];
      for (index_k = shelposk; index_k < shelposk + atoms->nshel_sh[kp1]; index_k++) {
        shelk1 = shells->type_sh[index_k];
        bfposl1  = 0;
        shelposl = atoms->shelposn_sh[lp1];
        for (index_l = shelposl; index_l < shelposl + atoms->nshel_sh[lp1]; index_l++) {
          shell1 = shells->type_sh[index_l];
          dimi = nd11 * nd12 * nd13 * nd14;
          F_ptr = pair_p->off[pair_p->ptr[ip1 * dim1 + jp1]];
          D_ptr = pair_p->off[pair_p->ptr[kp1 * dim1 + lp1]];
          allocate_integral_list(&integral_list_molecule, dimi, job, file);
          integral_list_molecule.num = 0;
          largest_integral = k_zero;
          for (i = 0; i < sheli1; i++) {
            for (j = 0; j < shelj1; j++) {
              for (k = 0; k < shelk1; k++) {
                for (l = 0; l < shell1; l++) {
                  testint = fabs(S1[F_ptr + (bfposi1 + i) * nd12 + bfposj1 + j] * S1[D_ptr + (bfposk1 + k) * nd14 + bfposl1 + l]);
                  if (testint > integral_rejection_threshold_sqrd) {
                    if (largest_integral < testint) largest_integral = testint;
                    integral_list_molecule.value[integral_list_molecule.num] = sqrt(testint);
                    integral_list_molecule.i[integral_list_molecule.num] = bfposi1 + i;
                    integral_list_molecule.j[integral_list_molecule.num] = bfposj1 + j;
                    integral_list_molecule.k[integral_list_molecule.num] = bfposk1 + k;
                    integral_list_molecule.l[integral_list_molecule.num] = bfposl1 + l;
                    integral_list_molecule.num++;
                   }
                  }
                 }
                }
               }

          if (integral_list_molecule.num == 0) { 
            shell_count++; 
            bfposl1  += shell1;
            free_integral_list(&integral_list_molecule,job);
            continue; 
           }

          else if (largest_integral > 1.0e-04) { 
            start_index[shell_count] = 1;
            shell_count++; 
            bfposl1  += shell1;
            free_integral_list(&integral_list_molecule,job);
            continue; 
           }

            //start_index[shell_count] = 1;

          for (p = 0; p < quad->tot; p++) {

            if (start_index[shell_count] == 1) break;

            pm = quad->p[p];
            //op = quad->k[p];

            for (s = 0; s < job->spin_dim; s++) {

              switch (pm) {

              case 0:

              jp0 = quad->cell2[0];
              kp0 = quad->cell3[0];
              lp0 = quad->cell4[0];
   
              p_coulomb_ints_i = integral_list_molecule.i;
              p_coulomb_ints_j = integral_list_molecule.j;
              p_coulomb_ints_k = integral_list_molecule.k;
              p_coulomb_ints_l = integral_list_molecule.l;
   
              p_exchange_ints_i = integral_list_molecule.i;
              p_exchange_ints_j = integral_list_molecule.j;
              p_exchange_ints_k = integral_list_molecule.k;
              p_exchange_ints_l = integral_list_molecule.l;

              break;
   
              case 1:
 
              jp0 = quad->cell2[0];
              kp0 = quad->cell4[0];
              lp0 = quad->cell3[0];
           
              p_coulomb_ints_i = integral_list_molecule.i;
              p_coulomb_ints_j = integral_list_molecule.j;
              p_coulomb_ints_k = integral_list_molecule.l;
              p_coulomb_ints_l = integral_list_molecule.k;
           
              p_exchange_ints_i = integral_list_molecule.i;
              p_exchange_ints_j = integral_list_molecule.j;
              p_exchange_ints_k = integral_list_molecule.l;
              p_exchange_ints_l = integral_list_molecule.k;
             
              break;
               
              case 2:

              jp0 = quad->cell1[0];
              kp0 = quad->cell3[0];
              lp0 = quad->cell4[0];
           
              p_coulomb_ints_i = integral_list_molecule.j;
              p_coulomb_ints_j = integral_list_molecule.i;
              p_coulomb_ints_k = integral_list_molecule.k;
              p_coulomb_ints_l = integral_list_molecule.l;
           
              p_exchange_ints_i = integral_list_molecule.j;
              p_exchange_ints_j = integral_list_molecule.i;
              p_exchange_ints_k = integral_list_molecule.k;
              p_exchange_ints_l = integral_list_molecule.l;
           
              break;
             
              case 3:
 
              jp0 = quad->cell1[0];
              kp0 = quad->cell4[0];
              lp0 = quad->cell3[0];
           
              p_coulomb_ints_i = integral_list_molecule.j;
              p_coulomb_ints_j = integral_list_molecule.i;
              p_coulomb_ints_k = integral_list_molecule.l;
              p_coulomb_ints_l = integral_list_molecule.k;
           
              p_exchange_ints_i = integral_list_molecule.j;
              p_exchange_ints_j = integral_list_molecule.i;
              p_exchange_ints_k = integral_list_molecule.l;
              p_exchange_ints_l = integral_list_molecule.k;
           
              break;
             
              case 4:
           
              jp0 = quad->cell4[0];
              kp0 = quad->cell1[0];
              lp0 = quad->cell2[0];
 
              p_coulomb_ints_i = integral_list_molecule.k;
              p_coulomb_ints_j = integral_list_molecule.l;
              p_coulomb_ints_k = integral_list_molecule.i;
              p_coulomb_ints_l = integral_list_molecule.j;
 
              p_exchange_ints_i = integral_list_molecule.k;
              p_exchange_ints_j = integral_list_molecule.l;
              p_exchange_ints_k = integral_list_molecule.i;
              p_exchange_ints_l = integral_list_molecule.j;
    
              break;
             
              case 5:
           
              jp0 = quad->cell4[0];
              kp0 = quad->cell2[0];
              lp0 = quad->cell1[0];
  
              p_coulomb_ints_i = integral_list_molecule.k;
              p_coulomb_ints_j = integral_list_molecule.l;
              p_coulomb_ints_k = integral_list_molecule.j;
              p_coulomb_ints_l = integral_list_molecule.i;
           
              p_exchange_ints_i = integral_list_molecule.k;
              p_exchange_ints_j = integral_list_molecule.l;
              p_exchange_ints_k = integral_list_molecule.j;
              p_exchange_ints_l = integral_list_molecule.i;
           
              break;
             
              case 6:
           
              jp0 = quad->cell3[0];
              kp0 = quad->cell1[0];
              lp0 = quad->cell2[0];
           
              p_coulomb_ints_i = integral_list_molecule.l;
              p_coulomb_ints_j = integral_list_molecule.k;
              p_coulomb_ints_k = integral_list_molecule.i;
              p_coulomb_ints_l = integral_list_molecule.j;
           
              p_exchange_ints_i = integral_list_molecule.l;
              p_exchange_ints_j = integral_list_molecule.k;
              p_exchange_ints_k = integral_list_molecule.i;
              p_exchange_ints_l = integral_list_molecule.j;
           
              break;
             
              case 7:
           
              jp0 = quad->cell3[0];
              kp0 = quad->cell2[0];
              lp0 = quad->cell1[0];
           
              p_coulomb_ints_i = integral_list_molecule.l;
              p_coulomb_ints_j = integral_list_molecule.k;
              p_coulomb_ints_k = integral_list_molecule.j;
              p_coulomb_ints_l = integral_list_molecule.i;
           
              p_exchange_ints_i = integral_list_molecule.l;
              p_exchange_ints_j = integral_list_molecule.k;
              p_exchange_ints_k = integral_list_molecule.j;
              p_exchange_ints_l = integral_list_molecule.i;
           
              break;
             
             } // close switch (pm)
           
              ip = quad->cell1[p];
              jp = quad->cell2[p];
              kp = quad->cell3[p];
              lp = quad->cell4[p];
 
              nd1 = atoms->bfnnumb_sh[ip];
              nd2 = atoms->bfnnumb_sh[jp];
              nd3 = atoms->bfnnumb_sh[kp];
              nd4 = atoms->bfnnumb_sh[lp0];

              //dimc = nd1 * nd2;
              //dime = nd1 * nd3;

              molecule_ints_num     = integral_list_molecule.num;
              p_coulomb_ints_value  = integral_list_molecule.value;
              p_exchange_ints_value = integral_list_molecule.value;
             
              if (job->scf_coulomb == 1 && pair_p->uniq[ip * dim1 + jp] == -1) {
              //AllocateDoubleArray(&Fock_2c,&dimc,job);
              //AllocateDoubleArray(&Fock_2c_temp,&dimc,job);
              for (s1 = 0; s1 < job->spin_dim; s1++) {
              flag = 0;
              //ResetDoubleArray(Fock_2c,&dimc);
              //ResetDoubleArray(Fock_2c_temp,&dimc);
              D_ptr = s1 * job->dimf + pair_p->off[pair_p->ptr[kp0 * dim1 + lp0]];
              for (m = 0; m < molecule_ints_num; m++) {
              //Fock_temp_offset      =         nd2 * *p_coulomb_ints_i + *p_coulomb_ints_j;
              Density_matrix_offset = D_ptr + nd4 * *p_coulomb_ints_k + *p_coulomb_ints_l;
              //Fock_2c_temp[Fock_temp_offset] += *p_coulomb_ints_value * F[Density_matrix_offset];
              if (fabs(*p_coulomb_ints_value * F[Density_matrix_offset]) > 1e-10) {
                start_index[shell_count] = 1; 
                flag = 1;
                //break;
               }
              if (flag == 1) break;
              p_coulomb_ints_value++;
              p_coulomb_ints_i++;
              p_coulomb_ints_j++;
              p_coulomb_ints_k++;
              p_coulomb_ints_l++;
              }
              //rotate_sum_block(Fock_2c_temp,Fock_2c,ip,jp,op,atoms,shells,symmetry,job,file);
              //for (i1 = 0; i1 < dimc; i1++) {
                //if (fabs(Fock_2c[i1]) > 1e-10) {
                //if (fabs(Fock_2c_temp[i1]) > 1e-10) {
                //start_index[shell_count] = 1; 
                //flag = 1;
                //break;
               //}
                //if (flag == 1) break;
              //}
              p_coulomb_ints_value -= molecule_ints_num;
              p_coulomb_ints_i     -= molecule_ints_num;
              p_coulomb_ints_j     -= molecule_ints_num;
              p_coulomb_ints_k     -= molecule_ints_num;
              p_coulomb_ints_l     -= molecule_ints_num;
             }
              //DestroyDoubleArray(&Fock_2c,&dimc,job);
              //DestroyDoubleArray(&Fock_2c_temp,&dimc,job);
             }

              if (job->scf_exchange == 1 && job->xc_hfx == 1 && pair_p->uniq[ip * dim1 + kp] == -1 && start_index[shell_count] == 0) {
              flag = 0;
              //AllocateDoubleArray(&Fock_2e,&dime,job);
              //AllocateDoubleArray(&Fock_2e_temp,&dime,job);
              //ResetDoubleArray(Fock_2e,&dime);
              //ResetDoubleArray(Fock_2e_temp,&dime);
              D_ptr = s * job->dimf + pair_p->off[pair_p->ptr[jp0 * dim1 + lp0]];
              for (m = 0; m < molecule_ints_num; m++) {
              //Fock_temp_offset      =         nd3 * *p_exchange_ints_i + *p_exchange_ints_k;
              Density_matrix_offset = D_ptr + nd4 * *p_exchange_ints_j + *p_exchange_ints_l;
              //Fock_2e_temp[Fock_temp_offset] -= *p_exchange_ints_value * F[Density_matrix_offset] / two;
              if (fabs(*p_exchange_ints_value * F[Density_matrix_offset] / two) > 1e-10) {
              start_index[shell_count] = 1; 
              flag = 1;
             }
              if (flag == 1) break;
              p_exchange_ints_value++;
              p_exchange_ints_i++;
              p_exchange_ints_j++;
              p_exchange_ints_k++;
              p_exchange_ints_l++;
              }
              //rotate_sum_block(Fock_2e_temp,Fock_2e,ip,kp,op,atoms,shells,symmetry,job,file);
              //for (i1 = 0; i1 < dime; i1++) {
                //if (fabs(Fock_2e[i1]) > 1e-10) {
                //if (fabs(Fock_2e_temp[i1]) > 1e-10) {
                //start_index[shell_count] = 1; 
                //flag = 1;
                //break;
               //}
                //if (flag == 1) break;
              //}
              //DestroyDoubleArray(&Fock_2e,&dime,job);
              //DestroyDoubleArray(&Fock_2e_temp,&dime,job);
             }

            } // close loop on s
           } // close loop on p
 
          shell_count++;
          free_integral_list(&integral_list_molecule,job);
          bfposl1  += shell1;
         }
        bfposk1  += shelk1;
       }
      bfposj1  += shelj1;
     }
    bfposi1  += sheli1;
   }

}

void shell_screen_crystal_coulomb_direct(int *start_index, double *S1, double *F, PAIR_TRAN *pair_p, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, REAL_LATTICE_TABLES *R_tables, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int index_i, index_j, index_k, index_l;
int shelposi, shelposj, shelposk, shelposl;
int bfposi1, bfposj1, bfposk1, bfposl1;
int sheli1, shelj1, shelk1, shell1;
int shell_count, flag;
int i, j, k, l, ip, jp, kp, lp, gi, gj, gk, gl;
int i1, m, pm, op, p, s, s1;
int dim1 = atoms->number_of_atoms_in_unit_cell, dim2 = dim1 * dim1, dimi;
int dimc, dime;
int ip0, jp0, kp0, lp0, gl0;
int ip1, jp1, kp1, lp1;
int nd1, nd2, nd3, nd4;
int nd11, nd12, nd13, nd14;
int Fock_temp_offset, Density_matrix_offset; 
int F_ptr, D_ptr;
int coulomb_ints_num;
int *p_coulomb_ints_i, *p_coulomb_ints_j, *p_coulomb_ints_k, *p_coulomb_ints_l;
//int *p_exchange_ints_i, *p_exchange_ints_j, *p_exchange_ints_k, *p_exchange_ints_l;
//double *Fock_2c, *Fock_2e, *Fock_2c_temp, *Fock_2e_temp;
double *p_coulomb_ints_value;
//double *p_molecule_ints_value;
double largest_integral, testint, integral_rejection_threshold_sqrd;
INTEGRAL_LIST integral_list_coulomb;

  integral_rejection_threshold_sqrd = integral_rejection_threshold * integral_rejection_threshold;
  //integral_rejection_threshold_sqrd = 1e-16;

  ip1 = quad->cell1[0];
  jp1 = quad->cell2[0];
  kp1 = quad->cell3[0];
  lp1 = quad->cell4[0];

  nd11 = atoms->bfnnumb_sh[ip1];
  nd12 = atoms->bfnnumb_sh[jp1];
  nd13 = atoms->bfnnumb_sh[kp1];
  nd14 = atoms->bfnnumb_sh[lp1];

  shell_count = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip1];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip1]; index_i++) {
    sheli1 = shells->type_sh[index_i];
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp1];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp1]; index_j++) {
      shelj1 = shells->type_sh[index_j];
      bfposk1  = 0;
      shelposk = atoms->shelposn_sh[kp1];
      for (index_k = shelposk; index_k < shelposk + atoms->nshel_sh[kp1]; index_k++) {
        shelk1 = shells->type_sh[index_k];
        bfposl1  = 0;
        shelposl = atoms->shelposn_sh[lp1];
        for (index_l = shelposl; index_l < shelposl + atoms->nshel_sh[lp1]; index_l++) {
          shell1 = shells->type_sh[index_l];
          dimi = nd11 * nd12 * nd13 * nd14;
          F_ptr = pair_p->off[pair_p->ptr[ip1 * dim1 + jp1]];
          D_ptr = pair_p->off[pair_p->ptr[kp1 * dim1 + lp1]];
          allocate_integral_list(&integral_list_coulomb, dimi, job, file);
          integral_list_coulomb.num = 0;
          largest_integral = k_zero;
          for (i = 0; i < sheli1; i++) {
            for (j = 0; j < shelj1; j++) {
              for (k = 0; k < shelk1; k++) {
                for (l = 0; l < shell1; l++) {
                  testint = fabs(S1[F_ptr + (bfposi1 + i) * nd12 + bfposj1 + j] * S1[D_ptr + (bfposk1 + k) * nd14 + bfposl1 + l]);
                  if (testint > integral_rejection_threshold_sqrd) {
                    if (largest_integral < testint) largest_integral = testint;
                    integral_list_coulomb.value[integral_list_coulomb.num] = sqrt(testint);
                    integral_list_coulomb.i[integral_list_coulomb.num] = bfposi1 + i;
                    integral_list_coulomb.j[integral_list_coulomb.num] = bfposj1 + j;
                    integral_list_coulomb.k[integral_list_coulomb.num] = bfposk1 + k;
                    integral_list_coulomb.l[integral_list_coulomb.num] = bfposl1 + l;
                    integral_list_coulomb.num++;
                   }
                  }
                 }
                }
               }

          if (integral_list_coulomb.num == 0) { 
            shell_count++; 
            bfposl1  += shell1;
            free_integral_list(&integral_list_coulomb,job);
            continue; 
           }

          else if (largest_integral > 1.0e-04) { 
            start_index[shell_count] = 1;
            shell_count++; 
            bfposl1  += shell1;
            free_integral_list(&integral_list_coulomb,job);
            continue; 
           }

            //start_index[shell_count] = 1;

          for (p = 0; p < quad->tot; p++) {

            if (start_index[shell_count] == 1) break;

            pm = quad->p[p];
            //op = quad->k[p];

            for (s = 0; s < job->spin_dim; s++) {

              switch (pm) {

              case 0:

              jp0 = quad->cell2[0];
              kp0 = quad->cell3[0];
              lp0 = quad->cell4[0];
              gl0 = quad->latt4[0];

              p_coulomb_ints_i = integral_list_coulomb.i;
              p_coulomb_ints_j = integral_list_coulomb.j;
              p_coulomb_ints_k = integral_list_coulomb.k;
              p_coulomb_ints_l = integral_list_coulomb.l;
   
              break;
   
              case 1:
 
              jp0 = quad->cell2[0];
              kp0 = quad->cell4[0];
              lp0 = quad->cell3[0];
              gl0 = R_tables->diffvec[quad->latt4[0]];

              p_coulomb_ints_i = integral_list_coulomb.i;
              p_coulomb_ints_j = integral_list_coulomb.j;
              p_coulomb_ints_k = integral_list_coulomb.l;
              p_coulomb_ints_l = integral_list_coulomb.k;
             
              break;
               
              case 2:

              jp0 = quad->cell1[0];
              kp0 = quad->cell3[0];
              lp0 = quad->cell4[0];
              gl0 = quad->latt4[0];

              p_coulomb_ints_i = integral_list_coulomb.j;
              p_coulomb_ints_j = integral_list_coulomb.i;
              p_coulomb_ints_k = integral_list_coulomb.k;
              p_coulomb_ints_l = integral_list_coulomb.l;
           
              break;
             
              case 3:
 
              jp0 = quad->cell1[0];
              kp0 = quad->cell4[0];
              lp0 = quad->cell3[0];
              gl0 = R_tables->diffvec[quad->latt4[0]];

              p_coulomb_ints_i = integral_list_coulomb.j;
              p_coulomb_ints_j = integral_list_coulomb.i;
              p_coulomb_ints_k = integral_list_coulomb.l;
              p_coulomb_ints_l = integral_list_coulomb.k;
           
              break;
             
              case 4:
           
              jp0 = quad->cell4[0];
              kp0 = quad->cell1[0];
              lp0 = quad->cell2[0];
              gl0 = quad->latt2[0];

              p_coulomb_ints_i = integral_list_coulomb.k;
              p_coulomb_ints_j = integral_list_coulomb.l;
              p_coulomb_ints_k = integral_list_coulomb.i;
              p_coulomb_ints_l = integral_list_coulomb.j;
    
              break;
             
              case 5:
           
              jp0 = quad->cell4[0];
              kp0 = quad->cell2[0];
              lp0 = quad->cell1[0];
              gl0 = R_tables->diffvec[quad->latt2[0]];

              p_coulomb_ints_i = integral_list_coulomb.k;
              p_coulomb_ints_j = integral_list_coulomb.l;
              p_coulomb_ints_k = integral_list_coulomb.j;
              p_coulomb_ints_l = integral_list_coulomb.i;
           
              break;
             
              case 6:
           
              jp0 = quad->cell3[0];
              kp0 = quad->cell1[0];
              lp0 = quad->cell2[0];
              gl0 = quad->latt2[0];

              p_coulomb_ints_i = integral_list_coulomb.l;
              p_coulomb_ints_j = integral_list_coulomb.k;
              p_coulomb_ints_k = integral_list_coulomb.i;
              p_coulomb_ints_l = integral_list_coulomb.j;
           
              break;
             
              case 7:
           
              jp0 = quad->cell3[0];
              kp0 = quad->cell2[0];
              lp0 = quad->cell1[0];
              gl0 = R_tables->diffvec[quad->latt2[0]];

              p_coulomb_ints_i = integral_list_coulomb.l;
              p_coulomb_ints_j = integral_list_coulomb.k;
              p_coulomb_ints_k = integral_list_coulomb.j;
              p_coulomb_ints_l = integral_list_coulomb.i;
           
              break;
             
             } // close switch (pm)
           
              ip = quad->cell1[p];
              jp = quad->cell2[p];
              kp = quad->cell3[p];
              lp = quad->cell4[p];
              gj = quad->latt2[p];

              nd1 = atoms->bfnnumb_sh[ip];
              nd2 = atoms->bfnnumb_sh[jp];
              nd3 = atoms->bfnnumb_sh[kp];
              nd4 = atoms->bfnnumb_sh[lp0];

              //dimc = nd1 * nd2;
              //dime = nd1 * nd3;

              coulomb_ints_num     = integral_list_coulomb.num;
              p_coulomb_ints_value = integral_list_coulomb.value;
             
              //if (job->scf_coulomb == 1 && pair_p->uniq[ip * dim1 + jp] == -1) {
              if (job->scf_coulomb == 1 && pair_p->uniq[gj * dim2 + ip * dim1 + jp] == -1) {
              //AllocateDoubleArray(&Fock_2c,&dimc,job);
              //AllocateDoubleArray(&Fock_2c_temp,&dimc,job);
              for (s1 = 0; s1 < job->spin_dim; s1++) {
              flag = 0;
              //ResetDoubleArray(Fock_2c,&dimc);
              //ResetDoubleArray(Fock_2c_temp,&dimc);
              D_ptr = s1 * job->dimf + pair_p->off[pair_p->ptr[gl0 * dim2 + kp0 * dim1 + lp0]];
              //D_ptr = s1 * job->dimf + pair_p->off[pair_p->ptr[kp0 * dim1 + lp0]];
              for (m = 0; m < coulomb_ints_num; m++) {
              //Fock_temp_offset      =         nd2 * *p_coulomb_ints_i + *p_coulomb_ints_j;
              Density_matrix_offset = D_ptr + nd4 * *p_coulomb_ints_k + *p_coulomb_ints_l;
              //Fock_2c_temp[Fock_temp_offset] += *p_coulomb_ints_value * F[Density_matrix_offset];
              if (fabs(*p_coulomb_ints_value * F[Density_matrix_offset]) > 1e-10) {
                start_index[shell_count] = 1; 
                flag = 1;
                //break;
               }
              if (flag == 1) break;
              p_coulomb_ints_value++;
              p_coulomb_ints_i++;
              p_coulomb_ints_j++;
              p_coulomb_ints_k++;
              p_coulomb_ints_l++;
              }
              //rotate_sum_block(Fock_2c_temp,Fock_2c,ip,jp,op,atoms,shells,symmetry,job,file);
              //for (i1 = 0; i1 < dimc; i1++) {
                //if (fabs(Fock_2c[i1]) > 1e-10) {
                //if (fabs(Fock_2c_temp[i1]) > 1e-10) {
                //start_index[shell_count] = 1; 
                //flag = 1;
                //break;
               //}
                //if (flag == 1) break;
              //}
              p_coulomb_ints_value -= coulomb_ints_num;
              p_coulomb_ints_i     -= coulomb_ints_num;
              p_coulomb_ints_j     -= coulomb_ints_num;
              p_coulomb_ints_k     -= coulomb_ints_num;
              p_coulomb_ints_l     -= coulomb_ints_num;
             }
              //DestroyDoubleArray(&Fock_2c,&dimc,job);
              //DestroyDoubleArray(&Fock_2c_temp,&dimc,job);
             }
/*
              if (job->scf_exchange == 1 && job->xc_hfx == 1 && pair_p->uniq[ip * dim1 + kp] == -1 && start_index[shell_count] == 0) {
              flag = 0;
              //AllocateDoubleArray(&Fock_2e,&dime,job);
              //AllocateDoubleArray(&Fock_2e_temp,&dime,job);
              //ResetDoubleArray(Fock_2e,&dime);
              //ResetDoubleArray(Fock_2e_temp,&dime);
              D_ptr = s * job->dimf + pair_p->off[pair_p->ptr[jp0 * dim1 + lp0]];
              for (m = 0; m < molecule_ints_num; m++) {
              //Fock_temp_offset      =         nd3 * *p_exchange_ints_i + *p_exchange_ints_k;
              Density_matrix_offset = D_ptr + nd4 * *p_exchange_ints_j + *p_exchange_ints_l;
              //Fock_2e_temp[Fock_temp_offset] -= *p_exchange_ints_value * F[Density_matrix_offset] / two;
              if (fabs(*p_exchange_ints_value * F[Density_matrix_offset] / two) > 1e-10) {
              start_index[shell_count] = 1; 
              flag = 1;
             }
              if (flag == 1) break;
              p_exchange_ints_value++;
              p_exchange_ints_i++;
              p_exchange_ints_j++;
              p_exchange_ints_k++;
              p_exchange_ints_l++;
              }
              //rotate_sum_block(Fock_2e_temp,Fock_2e,ip,kp,op,atoms,shells,symmetry,job,file);
              //for (i1 = 0; i1 < dime; i1++) {
                //if (fabs(Fock_2e[i1]) > 1e-10) {
                //if (fabs(Fock_2e_temp[i1]) > 1e-10) {
                //start_index[shell_count] = 1; 
                //flag = 1;
                //break;
               //}
                //if (flag == 1) break;
              //}
              //DestroyDoubleArray(&Fock_2e,&dime,job);
              //DestroyDoubleArray(&Fock_2e_temp,&dime,job);
             }
*/

            } // close loop on s
           } // close loop on p
 
          shell_count++;
          free_integral_list(&integral_list_coulomb,job);
          bfposl1  += shell1;
         }
        bfposk1  += shelk1;
       }
      bfposj1  += shelj1;
     }
    bfposi1  += sheli1;
   }

}

void shell_screen_compute_integrals(double *S1, double *S2, PAIR_TRAN *pair_p, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int i, j, k, p, q;
int ip, jp, gi, gj;
int nd1, nd2, nd3, nd4, nd12, nd34;
int dim;
int count;

double *S0;

  dim = 0;
  AllocateDoubleArray(&S0,&job->dimp,job);
  for (i = 0; i < job->dimp; i++) S0[i] = k_zero;
  for (i = 0; i < job->dimf; i++) S1[i] = k_zero;

  switch (crystal->type[0]) {

  case 'C':

  for (i = 0; i < job->dimp; i++) S0[i] = k_zero;
  for (i = 0; i < job->dimf; i++) S1[i] = k_zero;

  dim = 0;
  for (p = 0; p < pair_p->nump; p++) {
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
    integrals_crystal_coulomb_screen_ijij(&S0[dim],ip,jp,gj,R,G,atoms,shells,gaussians,symmetry,crystal,job,file);
    dim += nd34;
   }
  
    expand_screening_integral_matrix(S0,S1,pair_p,atoms,shells,symmetry,job,file);

  for (i = 0; i < job->dimp; i++) S0[i] = k_zero;
  for (i = 0; i < job->dimf; i++) S2[i] = k_zero;

  dim = 0;
  for (p = 0; p < pair_p->nump; p++) {
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
    integrals_crystal_exchange_screen_ijij(&S0[dim],ip,jp,gj,R,atoms,shells,gaussians,symmetry,crystal,job,file);
    dim += nd34;
   }
  
    expand_screening_integral_matrix(S0,S2,pair_p,atoms,shells,symmetry,job,file);

  break;

  case 'S':
  case 'P':

  if (job->taskid == 0)
  fprintf(file.out,"shell_screen routine is for 3-D systems only\n");
  MPI_Finalize();
  exit(1);

  break;

  case 'M':

  for (i = 0; i < job->dimp; i++) S0[i] = k_zero;
  for (i = 0; i < job->dimf; i++) S1[i] = k_zero;

  dim = 0;
  for (p = 0; p < pair_p->nump; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    gj = 0;
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    integrals_molecule_screen(&S0[dim],ip,jp,R,atoms,shells,gaussians,symmetry,crystal,job,file);
    dim += nd34;
   }
  
    expand_screening_integral_matrix(S0,S1,pair_p,atoms,shells,symmetry,job,file);

  break;

   } // close switch

    DestroyDoubleArray(&S0,&job->dimp,job);

    if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"screened integrals matrix %3d %3d %3d %3d\n",pair_p->nump,dim,job->dimp,job->dimf);
    count = 0;
    for (p = 0; p < pair_p->nump; p++) {
      q = pair_p->posn[p];
      fprintf(file.out,"pair %d [%3d %3d] [%3d]\n",p,pair_p->cell1[q],pair_p->cell2[q],pair_p->latt2[q]);
      for (k = 0; k < pair_p->numb[p]; k++) {
        for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
          for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
            fprintf(file.out,"%10.2e ",S2[count]);
            count++;
           }
          fprintf(file.out,"\n");
         }
        fprintf(file.out,"\n");
       }
      fprintf(file.out,"\n\n");
     }
    }

}

void shell_screen_complex(Complex *S1, PAIR_TRAN *pair_p, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, Q_LATTICE *q_G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int i, j, k, p, q;
int ip, jp, gi, gj;
int nd1, nd2, nd3, nd4, nd12, nd34;
int dim;
int count;
Complex *S0;

  AllocateComplexArray(&S0,&job->dimp,job);

  switch (crystal->type[0]) {

  case 'C':

  for (i = 0; i < job->dimp; i++) S0[i] = k_zero;
  for (i = 0; i < job->dimf; i++) S1[i] = k_zero;

  dim = 0;
  for (p = 0; p < pair_p->nump; p++) {
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
    //fprintf(file.out,"pair %3d   %3d %3d  %3d\n",p,ip,jp,gj);
    integrals_crystal_screen_complex(&S0[dim],ip,jp,gj,R,G,q_G,atoms,shells,gaussians,symmetry,crystal,job,file);
    dim += nd34;
   }
  
    expand_screening_integral_matrix_complex(S0,S1,pair_p,atoms,shells,symmetry,job,file);

  break;

  case 'S':
  case 'P':
  case 'M':

  if (job->taskid == 0)
  fprintf(file.out,"shell_screen_coulomb routine is for 3-D systems only\n");
  MPI_Finalize();
  exit(1);

  break;

   } // close switch

    DestroyComplexArray(&S0,&job->dimp,job);

    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out,"screened complex integrals matrix %3d %3d %3d\n",pair_p->nump,dim,job->dimp);
    count = 0;
    if (job->taskid == 0 && job->verbosity > 1)
    for (p = 0; p < pair_p->nump; p++) {
      q = pair_p->posn[p];
      fprintf(file.out,"pair %d [%3d %3d] [%3d]\n",p,pair_p->cell1[q],pair_p->cell2[q],pair_p->latt2[q]);
      for (k = 0; k < pair_p->numb[p]; k++) {
        for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
          for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
            fprintf(file.out,"%9.2e  ", (S1[count]).real());
            //fprintf(file.out,"%9.2e %9.2e  ", (S1[count]).real(), (S1[count]).imag());
            count++;
           }
          fprintf(file.out,"\n");
         }
        fprintf(file.out,"\n");
       }
      fprintf(file.out,"\n\n");
     }

}

void shell_screen1(int *start_index, double *S1, PAIR_TRAN *pair_p, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, JOB_PARAM *job, FILES file)

{

int index_i, index_j, index_k, index_l;
int shelposi, shelposj, shelposk, shelposl;
int bfposi1, bfposj1, bfposk1, bfposl1;
int sheli1, shelj1, shelk1, shell1;
int shell_count, flag;
int i, j, k, l, ip, jp, kp, lp, gi, gj, gk, gl;
int ptr1, ptr2, i1, i2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int nd2, nd4;
double test;
double integral_rejection_threshold_sqrd;

  integral_rejection_threshold_sqrd = integral_rejection_threshold * integral_rejection_threshold * 1000000.0;

  ip = quad->cell1[0];
  jp = quad->cell2[0];
  kp = quad->cell3[0];
  lp = quad->cell4[0];
  gi = quad->latt1[0];
  gj = quad->latt2[0];
  gk = quad->latt3[0];
  gl = quad->latt4[0];
  ptr1 = pair_p->off[pair_p->ptr[gj * dim2 + ip * dim1 + jp]];
  ptr2 = pair_p->off[pair_p->ptr[gl * dim2 + kp * dim1 + lp]];
  nd2 = atoms->bfnnumb_sh[jp];
  nd4 = atoms->bfnnumb_sh[lp];

  shell_count = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
    sheli1 = shells->type_sh[index_i];
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
      shelj1 = shells->type_sh[index_j];
      bfposk1  = 0;
      shelposk = atoms->shelposn_sh[kp];
      for (index_k = shelposk; index_k < shelposk + atoms->nshel_sh[kp]; index_k++) {
        shelk1 = shells->type_sh[index_k];
        bfposl1  = 0;
        shelposl = atoms->shelposn_sh[lp];
        for (index_l = shelposl; index_l < shelposl + atoms->nshel_sh[lp]; index_l++) {
          shell1 = shells->type_sh[index_l];
          flag = 0;
          for (i = 0; i < sheli1; i++) {
            for (j = 0; j < shelj1; j++) {
              for (k = 0; k < shelk1; k++) {
                for (l = 0; l < shell1; l++) {
                  if (flag == 1) continue;
                  test = fabs(S1[ptr1 + (bfposi1 + i) * nd2 + bfposj1 + j] * S1[ptr2 + (bfposk1 + k) * nd4 + (bfposl1 + l)]);
                  if (test > integral_rejection_threshold_sqrd) {
                  start_index[shell_count] = 1; 
                  flag = 1;
                 }
                  //fprintf(file.out,"%3d %3d %3d %3d  %3d %3d %3d %3d %3d %9.2e %9.2e  %9.2e\n",\
                  ip,jp,kp,lp,bfposi1+i,bfposj1+j,bfposk1+k,bfposl1+l,start_index[shell_count],\
                  S1[ptr1 + (bfposi1 + i) * nd2 + bfposj1 + j],S1[ptr2 + (bfposk1 + k) * nd4 + (bfposl1 + l)],test);
                }
               }
              }
             }
            shell_count++;
            bfposl1  += shell1;
           }
          bfposk1  += shelk1;
         }
        bfposj1  += shelj1;
       }
      bfposi1  += sheli1;
     }

}

void shell_screen2(int *start_index, double *S2, PAIR_TRAN *pair_p, QUAD_TRAN *quad, REAL_LATTICE_TABLES *R_tables, ATOM *atoms, SHELL *shells, JOB_PARAM *job, FILES file)

{

int index_i, index_j, index_k, index_l;
int shelposi, shelposj, shelposk, shelposl;
int bfposi1, bfposj1, bfposk1, bfposl1;
int sheli1, shelj1, shelk1, shell1;
int shell_count, flag;
int i, j, k, l, ip, jp, kp, lp, gi, gj, gk, gl;
int ptr1, ptr2, i1, i2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int nd2, nd4;
double test;
double integral_rejection_threshold_sqrd;

  integral_rejection_threshold_sqrd = integral_rejection_threshold * integral_rejection_threshold * 1000000000000.0;
  //integral_rejection_threshold_sqrd = integral_rejection_threshold * integral_rejection_threshold * 100000000.0;

  ip = quad->cell1[0];
  jp = quad->cell2[0];
  kp = quad->cell3[0];
  lp = quad->cell4[0];
  gi = quad->latt1[0];
  gj = quad->latt2[0];
  gk = quad->latt3[0];
  gl = quad->latt4[0];

  gl = R_tables->diffvec[quad->latt4[0] * R_tables->margin_vector + quad->latt3[0]];
  if (gl > R_tables->last_vector) return;

  ptr1 = pair_p->off[pair_p->ptr[gj * dim2 + ip * dim1 + jp]];
  ptr2 = pair_p->off[pair_p->ptr[gl * dim2 + kp * dim1 + lp]];
  nd2 = atoms->bfnnumb_sh[jp];
  nd4 = atoms->bfnnumb_sh[lp];

  shell_count = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
    sheli1 = shells->type_sh[index_i];
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
      shelj1 = shells->type_sh[index_j];
      bfposk1  = 0;
      shelposk = atoms->shelposn_sh[kp];
      for (index_k = shelposk; index_k < shelposk + atoms->nshel_sh[kp]; index_k++) {
        shelk1 = shells->type_sh[index_k];
        bfposl1  = 0;
        shelposl = atoms->shelposn_sh[lp];
        for (index_l = shelposl; index_l < shelposl + atoms->nshel_sh[lp]; index_l++) {
          shell1 = shells->type_sh[index_l];
          flag = 0;
          for (i = 0; i < sheli1; i++) {
            for (j = 0; j < shelj1; j++) {
              for (k = 0; k < shelk1; k++) {
                for (l = 0; l < shell1; l++) {
                  if (flag == 1) continue;
                  test = fabs(S2[ptr1 + (bfposi1 + i) * nd2 + bfposj1 + j] * S2[ptr2 + (bfposk1 + k) * nd4 + (bfposl1 + l)]);
                  if (test > integral_rejection_threshold_sqrd) {
                  start_index[shell_count] = 1; 
                  flag = 1;
                 }
                }
               }
              }
             }
            shell_count++;
            bfposl1  += shell1;
           }
          bfposk1  += shelk1;
         }
        bfposj1  += shelj1;
       }
      bfposi1  += sheli1;
     }
}

void integral_screen_crystal_ija(int *start_index, ComplexMatrix *V_screen, Complex *S1, PAIR_TRAN *pair_p, TRIPLE_TRAN *triple, ATOM *atoms, SHELL *shells, ATOM *atoms_ax, SHELL *shells_ax, JOB_PARAM *job, FILES file)
//void shell_screen3(int *start_index, ComplexMatrix *V_screen, Complex *S1, PAIR_TRAN *pair_p, TRIPLE_TRAN *triple, ATOM *atoms, SHELL *shells, ATOM *atoms_ax, SHELL *shells_ax, JOB_PARAM *job, FILES file)

{

int index_i, index_j, index_k;
int shelposi, shelposj, shelposk;
int bfposi1, bfposj1, bfposk1;
int sheli1, shelj1, shelk1;
int shell_count, flag;
int i, j, k, ip, jp, kp, gi, gj, gk;
int ptr1, ptr2;
int dim1 = atoms->number_of_atoms_in_unit_cell;
int dim2 = dim1 * dim1;
int nd2;
double mod_test;
Complex test;

  ip = triple->cell1[0];
  jp = triple->cell2[0];
  kp = triple->cell3[0];
  gi = triple->latt1[0];
  gj = triple->latt2[0];
  gk = triple->latt3[0];
  ptr1 = pair_p->off[pair_p->ptr[gj * dim2 + ip * dim1 + jp]];
  ptr2 = atoms_ax->bfnposn_sh[kp];
  nd2 = atoms->bfnnumb_sh[jp];

  shell_count = 0;
  bfposi1  = 0;
  shelposi = atoms->shelposn_sh[ip];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
    sheli1 = shells->type_sh[index_i];
    bfposj1  = 0;
    shelposj = atoms->shelposn_sh[jp];
    for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
      shelj1 = shells->type_sh[index_j];
      bfposk1  = 0;
      shelposk = atoms_ax->shelposn_sh[kp];
      for (index_k = shelposk; index_k < shelposk + atoms_ax->nshel_sh[kp]; index_k++) {
        shelk1 = shells_ax->type_sh[index_k];
        flag = 0;
        for (i = 0; i < sheli1; i++) {
          for (j = 0; j < shelj1; j++) {
            for (k = 0; k < shelk1; k++) {
              if (flag == 1) continue;
              test = S1[ptr1 + (bfposi1 + i) * nd2 + (bfposj1 + j)] * V_screen->a[ptr2 + bfposk1 + k][ptr2 + bfposk1 + k];
              mod_test = sqrt(test.real() * test.real() + test.imag() * test.imag());
              //if (mod_test > 1.0e-12) {
              ////if (mod_test > 1.0e-10) {
              if (mod_test > 1.0e-07) {
              start_index[shell_count] = 1; 
              flag = 1;
              }
             }
            }
           }
          shell_count++;
          bfposk1  += shelk1;
         }
        bfposj1  += shelj1;
       }
      bfposi1  += sheli1;
     }

}

/*
void shell_overlap(int *start_index, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, REAL_LATTICE *R, JOB_PARAM *job, FILES file)

{

int index_i, index_j, index_k, index_l;
int shelposi, shelposj, shelposk, shelposl;
int shell_count;
int i, ip, jp, kp, lp, gi, gj, gk, gl, qp;
double ab, pab, p32, cd, pcd, overlap_ab, overlap_cd;
double AB, CD, SAB, SCD, KAB, KCD, R_AB_1esqrd, R_CD_1esqrd;
VECTOR_DOUBLE R_AB_1e, R_CD_1e;

  ip = quad->cell1[0];
  jp = quad->cell2[0];
  kp = quad->cell3[0];
  lp = quad->cell4[0];
  gi = quad->latt1[0];
  gj = quad->latt2[0];
  gk = quad->latt3[0];
  gl = quad->latt4[0];

  shell_count = 0;

  shelposi = atoms->shelposn_sh[ip];
  for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      shelposj = atoms->shelposn_sh[jp];
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        pab = shells->min_expo_sh[index_i] + shells->min_expo_sh[index_j];
        ab =  shells->min_expo_sh[index_i] * shells->min_expo_sh[index_j];
        AB = pow(four * ab, three_quarters);
        KAB = ab * R_AB_1esqrd / pab;
        p32 = pab * sqrt(pab);
        //SAB = pi32 * exp(-KAB) / p32;
        SAB = AB * exp(-KAB) / p32;
        //overlap_ab = shells->min_coef_sh[index_i] * shells->min_coef_sh[index_j] * SAB;
        shelposk = atoms->shelposn_sh[kp];
        for (index_k = shelposk; index_k < shelposk + atoms->nshel_sh[kp]; index_k++) {
          shelposl = atoms->shelposn_sh[lp];
          for (index_l = shelposl; index_l < shelposl + atoms->nshel_sh[lp]; index_l++) {
            R_CD_1e.comp1 = atoms->cell_vector[kp].comp1 + R->vec_ai[gk].comp1 - atoms->cell_vector[lp].comp1 - R->vec_ai[gl].comp1;
            R_CD_1e.comp2 = atoms->cell_vector[kp].comp2 + R->vec_ai[gk].comp1 - atoms->cell_vector[lp].comp2 - R->vec_ai[gl].comp2;
            R_CD_1e.comp3 = atoms->cell_vector[kp].comp3 + R->vec_ai[gk].comp1 - atoms->cell_vector[lp].comp3 - R->vec_ai[gl].comp3;
            R_CD_1esqrd = double_vec_dot(&R_CD_1e, &R_CD_1e);
            pcd = shells->min_expo_sh[index_k] + shells->min_expo_sh[index_l];
            cd = shells->min_expo_sh[index_k] * shells->min_expo_sh[index_l];
            CD = pow(four * cd, three_quarters);
            KCD = cd * R_CD_1esqrd / pcd;
            p32 = pcd * sqrt(pcd);
            //SCD = pi32 * exp(-KCD) / p32;
            SCD = CD * exp(-KCD) / p32;
            //overlap_cd = shells->min_coef_sh[index_k] * shells->min_coef_sh[index_l] * SCD;
            if (SAB > 1.0e-6 && SCD > 1.0e-6) start_index[shell_count] = 1;
            //if (SAB > 1.0e-11 || SCD > 1.0e-11) start_index[shell_count] = 1;
            //if (overlap_ab > 1.0e-13 || overlap_cd > 1.0e-13 || overlap_ab * overlap_cd > 1.0e-10) start_index[shell_count] = 1;
            //if (overlap_ab * overlap_cd > 1.0e-4) start_index[shell_count] = 1;
            //if (overlap_ab * overlap_cd > 1.0e-06) \
            fprintf(file.out,"shell_overlap %3d %3d %3d %3d   %3d %3d %3d %3d     %3d %3d %3d %3d   %16.8e\n",\
            index_i,index_j,index_k,index_l,gi,gj,gk,gl,ip,jp,kp,lp,overlap_ab * overlap_cd);
            //fprintf(file.out,"%3d %3d ex %9.2e %9.2e sab %9.2e ab %9.2e ex %9.2e %9.2e scd %9.2e cd %9.2e abcd %e %3d %3d %3d %3d %3d\n",\
            qp,shell_count,\
            shells->min_expo_sh[index_i], \
            shells->min_coef_sh[index_j],SAB,overlap_ab,shells->min_expo_sh[index_k],shells->min_coef_sh[index_l],SCD,overlap_cd, \
            overlap_ab * overlap_cd,start_index[shell_count],ip,jp,kp,lp);
            shell_count++;
           }
          }
         }
        }

}

void print_Fock_matrix(double *Fock, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file)

{

int i, j, p, q, s;
int count;

  fprintf(file.out,"Fock matrix unique pairs [%d]\n",pair_p->nump);
  count = 0;
  for (s = 0; s < job->spin_dim; s++) {
   for (p = 0; p < pair_p->nump; p++) {
    q = pair_p->posn[p];
    fprintf(file.out,"Pair %3d   %3d %3d   %3d\n",p,pair_p->cell1[q],pair_p->cell2[q],pair_p->latt2[q]);
    for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
      for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
        //fprintf(file.out,"%6.2f ",Fock[count]);
        fprintf(file.out,"%10.6f ",Fock[count]);
        count++;
       }
      fprintf(file.out,"\n");
     }
    fprintf(file.out,"\n\n");
   }
  fprintf(file.out,"\n\n");
 }
 fflush(file.out);

}
*/

