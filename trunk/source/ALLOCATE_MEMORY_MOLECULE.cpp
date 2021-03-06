
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson and Alin-Marin Elena        *
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
#include "ALLOCATE_MEMORY_MOLECULE.h"

using namespace std;

void allocate_JOB_PARAM(JOB_PARAM *job, FILES file)

{

  int i;

  job->verbosity   = 1;
  job->print_pairs = 0;

  job->C09 = 0;

  job->vectors = 0;
  job->values  = 0;

  job->memory = (int *) malloc(job->numtasks * sizeof(int));
  if (job->memory == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for memory array \n");
  MPI_Finalize(); 
  exit(1);;
  }

  job->max_memory = (int *) malloc(job->numtasks * sizeof(int));
  if (job->max_memory == NULL) {
  fprintf(stderr, "ERROR: There is not enough max_memory for max_memory array \n");
  MPI_Finalize(); 
  exit(1);;
  }

  for (i = 0; i < job->numtasks; i++)
  job->memory[job->taskid] = 0;

  for (i = 0; i < job->numtasks; i++)
  job->max_memory[i] = 0;

}

void free_JOB_PARAM(JOB_PARAM *job, FILES file)

{

  free(job->memory);
  free(job->max_memory);

}

void allocate_SYMMETRY(SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, lmax;

  lmax = 26;
  for (i = 2; i <= job->l_max; i++) 
  lmax += (2 * i + 1) * (2 * i + 1);

  symmetry->memory = 0;

  symmetry->num_ij = (int *) malloc((job->l_max + 2) * symmetry->number_of_operators * sizeof(int)); // 6 is s + p + sp + d + f + g
  symmetry->memory += (job->l_max + 2) * symmetry->number_of_operators * sizeof(int);
  if (symmetry->num_ij == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for num_ij array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->op_shift = (int *) malloc((job->l_max + 2) * symmetry->number_of_operators * sizeof(int)); // 6 = s + sp + p + d + f + g
  symmetry->memory += (job->l_max + 2) * symmetry->number_of_operators * sizeof(int);
  if (symmetry->op_shift == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for op_shift array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->inverse = (int *) malloc(symmetry->number_of_operators * sizeof(int));
  symmetry->memory += symmetry->number_of_operators * sizeof(int);
  if (symmetry->inverse == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for inverse array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->ind_i = (int *) malloc(lmax * symmetry->number_of_operators * sizeof(int)); // 100 = 1 + 9 + 16 + 25 + 49
  symmetry->memory += lmax * symmetry->number_of_operators * sizeof(int);
  if (symmetry->ind_i == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for ind_i array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->ind_j = (int *) malloc(lmax * symmetry->number_of_operators * sizeof(int));
  symmetry->memory += lmax * symmetry->number_of_operators * sizeof(int);
  if (symmetry->ind_j == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for ind_j array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->rot = (double *) malloc(lmax * symmetry->number_of_operators * sizeof(double));
  symmetry->memory += lmax * symmetry->number_of_operators * sizeof(double);
  if (symmetry->rot == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for rot array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->irr = (double *) malloc(9 * symmetry->number_of_operators * sizeof(double));
  symmetry->memory += 9 * symmetry->number_of_operators * sizeof(double);
  if (symmetry->irr == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for irr array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->character_table = (double *) malloc(symmetry->number_of_operators * symmetry->number_of_operators * sizeof(double));
  symmetry->memory += symmetry->number_of_operators * symmetry->number_of_operators * sizeof(double);
  if (symmetry->character_table == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for character_table array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->inr = (int *) malloc(9 * symmetry->number_of_operators * sizeof(int));
  symmetry->memory += 9 * symmetry->number_of_operators * sizeof(int);
  if (symmetry->inr == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for inr array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->grp_pm = (int *) malloc(64 * sizeof(int));
  symmetry->memory += 64 * sizeof(int);
  if (symmetry->grp_pm == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for grp_pm array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->irp_dim_k = (int *) malloc(symmetry->number_of_operators * symmetry->number_of_operators * sizeof(int));
  symmetry->memory += symmetry->number_of_operators * symmetry->number_of_operators * sizeof(int);
  if (symmetry->irp_dim_k == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for irp_dim_k array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->grp_k = (int *) malloc(symmetry->number_of_operators * symmetry->number_of_operators * sizeof(int));
  symmetry->memory += symmetry->number_of_operators * symmetry->number_of_operators * sizeof(int);
  if (symmetry->grp_k == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for grp_k array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->cls_num_pm = (int *) malloc(64 * sizeof(int));
  symmetry->memory += 64 * sizeof(int);
  if (symmetry->cls_num_pm == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for cls_num_pm array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->cls_num_k = (int *) malloc(symmetry->number_of_operators * sizeof(int));
  symmetry->memory += symmetry->number_of_operators * sizeof(int);
  if (symmetry->cls_num_k == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for cls_num_k array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->cls_pos_pm = (int *) malloc(64 * sizeof(int));
  symmetry->memory += 64 * sizeof(int);
  if (symmetry->cls_pos_pm == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for cls_pos_pm array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->cls_pos_k = (int *) malloc(symmetry->number_of_operators * sizeof(int));
  symmetry->memory += symmetry->number_of_operators * sizeof(int);
  if (symmetry->cls_pos_k == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for cls_pos_k array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->cls_pm = (int *) malloc(64 * sizeof(int));
  symmetry->memory += 64 * sizeof(int);
  if (symmetry->cls_pm == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for cls_pos_pm array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->cls_k = (int *) malloc(symmetry->number_of_operators * sizeof(int));
  symmetry->memory += symmetry->number_of_operators * sizeof(int);
  if (symmetry->cls_k == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for cls_pos_k array \n");
  MPI_Finalize(); 
  exit(1);
  }

  symmetry->taur = (VECTOR_DOUBLE *) malloc(symmetry->number_of_operators * sizeof(VECTOR_DOUBLE));
  symmetry->memory += symmetry->number_of_operators * sizeof(VECTOR_DOUBLE);
  if (symmetry->taur == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for taur array \n");
  MPI_Finalize(); 
  exit(1);
  }

  job->memory[job->taskid] += symmetry->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"%d SYMMETRY OPERATORS ALLOCATED %d MEMORY ALLOCATED %d\n",symmetry->number_of_operators,symmetry->memory,job->memory[job->taskid]);

}

void free_SYMMETRY(SYMMETRY *symmetry, JOB_PARAM *job)

{

  job->memory[job->taskid] -= symmetry->memory;

  free(symmetry->num_ij);
  free(symmetry->op_shift);
  free(symmetry->inverse);
  free(symmetry->ind_i);
  free(symmetry->ind_j);
  free(symmetry->rot);
  free(symmetry->irr);
  free(symmetry->character_table);
  free(symmetry->inr);
  free(symmetry->grp_pm);
  free(symmetry->irp_dim_k);
  free(symmetry->grp_k);
  free(symmetry->cls_num_pm);
  free(symmetry->cls_num_k);
  free(symmetry->cls_pos_pm);
  free(symmetry->cls_pos_k);
  free(symmetry->cls_pm);
  free(symmetry->cls_k);
  free(symmetry->taur);

}

void allocate_ATOM(ATOM *atoms, JOB_PARAM *job, FILES file)

{

  int i;

  atoms->cell_vector = (VECTOR_DOUBLE *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(VECTOR_DOUBLE));
  if (atoms->cell_vector == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->atomic_number = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->atomic_number == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->nshel = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->nshel == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->shelposn = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->shelposn == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->gausposn = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->gausposn == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->nshel_sh = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->nshel_sh == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->shelposn_sh = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->shelposn_sh == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->gausposn_sh = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->gausposn_sh == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->shlposn_sh = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->shlposn_sh == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->bfnposn_sh = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->bfnposn_sh == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->bfnnumb_sh = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->bfnnumb_sh == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->shlposn = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->shlposn == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->bfnposn = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->bfnposn == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->bfnnumb = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->bfnnumb == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->basis_set = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->basis_set == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->magnetic = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->magnetic == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }
  for (i=0;i<atoms->number_of_atoms_in_unit_cell;i++) {
    atoms->magnetic[i] = 1; // initialise to non-magnetic state as default
   }

  atoms->spin = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->spin == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->uniq = (int *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (atoms->uniq == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);
  }

  atoms->pop = (double *) malloc(2 * atoms->number_of_atoms_in_unit_cell * sizeof(double));
  if (atoms->pop == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atoms in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atoms->memory = 16 * atoms->number_of_atoms_in_unit_cell * sizeof(int) + atoms->number_of_atoms_in_unit_cell * (2 * sizeof(double) + sizeof(VECTOR_DOUBLE));
  job->memory[job->taskid] += atoms->memory;
  
  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"%d ATOMS ALLOCATED %d MEMORY ALLOCATED\n",atoms->number_of_atoms_in_unit_cell,atoms->memory);

}

void free_ATOM(ATOM *atoms, JOB_PARAM *job)

{

  job->memory[job->taskid] -= atoms->memory;

  free(atoms->cell_vector);
  free(atoms->nshel);
  free(atoms->shelposn);
  free(atoms->gausposn);
  free(atoms->bfnposn);
  free(atoms->bfnnumb);
  free(atoms->atomic_number);
  free(atoms->nshel_sh);
  free(atoms->shelposn_sh);
  free(atoms->gausposn_sh);
  free(atoms->bfnposn_sh);
  free(atoms->bfnnumb_sh);
  free(atoms->shlposn);
  free(atoms->shlposn_sh);
  free(atoms->basis_set);
  free(atoms->magnetic);
  free(atoms->spin);
  free(atoms->uniq);
  free(atoms->pop);

}

void allocate_ATOM_TRAN(ATOM_TRAN *atom_p, ATOM *atoms, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  atom_p->K = (int *) malloc(atoms->number_of_atoms_in_unit_cell * symmetry->number_of_operators * sizeof(int));
  if (atom_p->K == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atom_p in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atom_p->O = (int *) malloc(atoms->number_of_atoms_in_unit_cell * symmetry->number_of_operators * sizeof(int));
  if (atom_p->O == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atom_p in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atom_p->P = (int *) malloc(atoms->number_of_atoms_in_unit_cell * symmetry->number_of_operators * sizeof(int));
  if (atom_p->P == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atom_p in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atom_p->expo = (double *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(double));
  if (atom_p->expo == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atom_p in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atom_p->coef = (double *) malloc(atoms->number_of_atoms_in_unit_cell * sizeof(double));
  if (atom_p->coef == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atom_p in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atom_p->numb = (int *) malloc(atoms->number_of_unique_atoms * sizeof(int));
  if (atom_p->P == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atom_p in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atom_p->posn = (int *) malloc(atoms->number_of_unique_atoms * sizeof(int));
  if (atom_p->P == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for atom_p in main() \n");
  MPI_Finalize(); exit(1);;
  }

  atom_p->memory = 4 * symmetry->number_of_operators * atoms->number_of_atoms_in_unit_cell * sizeof(int) + \
  2 * atoms->number_of_atoms_in_unit_cell * sizeof(double);
  job->memory[job->taskid] += atom_p->memory;
  
  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"%d ATOM_TRAN ATOMS ALLOCATED %d MEMORY ALLOCATED\n",atoms->number_of_atoms_in_unit_cell,atom_p->memory);

}

void free_ATOM_TRAN(ATOM_TRAN *atom_p, JOB_PARAM *job)

{

  job->memory[job->taskid] -= atom_p->memory;

  free(atom_p->K);
  free(atom_p->O);
  free(atom_p->P);
  free(atom_p->expo);
  free(atom_p->coef);
  free(atom_p->numb);
  free(atom_p->posn);

}

void allocate_SHELL_GAUSSIAN(SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  shells->ng = (int *) malloc(atoms->number_of_shells_in_unit_cell * sizeof(int));
  if (shells->ng == NULL) {
  fprintf(file.out, "ERROR: There is not enough memory for shells in main()\n");
  MPI_Finalize(); 
  exit(1);
  }

  shells->imax = (int *) malloc(atoms->number_of_shells_in_unit_cell * sizeof(int));
  shells->type = (int *) malloc(atoms->number_of_shells_in_unit_cell * sizeof(int));
  shells->type1 = (int *) malloc(atoms->number_of_shells_in_unit_cell * sizeof(int));
  shells->ord = (int *) malloc(atoms->number_of_shells_in_unit_cell * sizeof(int));
  shells->opp = (int *) malloc(atoms->number_of_shells_in_unit_cell * sizeof(int));
  shells->ind_i = (int *) malloc(56 * sizeof(int));
  shells->ind_j = (int *) malloc(56 * sizeof(int));
  shells->num_ij = (int *) malloc(6 * sizeof(int));
  shells->rot = (double *) malloc(56 * sizeof(double));
  shells->nele = (double *) malloc(atoms->number_of_shells_in_unit_cell * sizeof(double));
  shells->pop = (double *) malloc(2 * atoms->number_of_shells_in_unit_cell * sizeof(double));
  shells->pop_sh = (double *) malloc(2 * atoms->number_of_sh_shells_in_unit_cell * sizeof(double));
  shells->min_coef_sh = (double *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(double));
  shells->min_expo_sh = (double *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(double));
  shells->ng_sh = (int *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(int));
  shells->imax_sh = (int *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(int));
  shells->type_sh = (int *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(int));
  shells->type1_sh = (int *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(int));
  shells->cart = (int *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(int));
  shells->shar = (int *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(int));
  shells->ord_sh = (int *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(int));
  shells->opp_sh = (int *) malloc(atoms->number_of_sh_shells_in_unit_cell * sizeof(int));
  gaussians->expo = (double *) malloc(atoms->number_of_gaussians_in_unit_cell * sizeof(double));
  gaussians->sc = (double *) malloc(atoms->number_of_gaussians_in_unit_cell * sizeof(double));
  gaussians->pc = (double *) malloc(atoms->number_of_gaussians_in_unit_cell * sizeof(double));
  gaussians->dc = (double *) malloc(atoms->number_of_gaussians_in_unit_cell * sizeof(double));
  gaussians->fc = (double *) malloc(atoms->number_of_gaussians_in_unit_cell * sizeof(double));
  gaussians->gc = (double *) malloc(atoms->number_of_gaussians_in_unit_cell * sizeof(double));
  gaussians->expo_sh = (double *) malloc(atoms->number_of_sh_gaussians_in_unit_cell * sizeof(double));
  gaussians->c_sh = (double *) malloc(atoms->number_of_sh_gaussians_in_unit_cell * sizeof(double));

  shells->memory = (6 * atoms->number_of_shells_in_unit_cell + 8 * atoms->number_of_sh_shells_in_unit_cell + 118) * sizeof(int) + \
  (5 * atoms->number_of_shells_in_unit_cell    + 2 * atoms->number_of_sh_shells_in_unit_cell + 56) * sizeof(double);

  gaussians->memory = (6 * atoms->number_of_gaussians_in_unit_cell + 2 * atoms->number_of_sh_gaussians_in_unit_cell) * sizeof(double);

  job->memory[job->taskid] += shells->memory;
  job->memory[job->taskid] += gaussians->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"%d SHELLS %d GAUSSIANS ALLOCATED MEMORY ALLOCATED %5d + %5d\n", \
  atoms->number_of_shells_in_unit_cell,atoms->number_of_gaussians_in_unit_cell,shells->memory,gaussians->memory);

}

void free_SHELL_GAUSSIAN(SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job)

{

  job->memory[job->taskid] -= shells->memory;
  job->memory[job->taskid] -= gaussians->memory;

  free(shells->ng);
  free(shells->imax);
  free(shells->type);
  free(shells->type1);
  free(shells->ord);
  free(shells->opp);
  free(shells->ind_i);
  free(shells->ind_j);
  free(shells->num_ij);
  free(shells->rot);
  free(shells->nele);
  free(shells->pop);
  free(shells->pop_sh);
  free(shells->min_coef_sh);
  free(shells->min_expo_sh);
  free(shells->ng_sh);
  free(shells->imax_sh);
  free(shells->type_sh);
  free(shells->type1_sh);
  free(shells->shar);
  free(shells->cart);
  free(shells->ord_sh);
  free(shells->opp_sh);
  free(gaussians->expo);
  free(gaussians->sc);
  free(gaussians->pc);
  free(gaussians->dc);
  free(gaussians->fc);
  free(gaussians->gc);
  free(gaussians->expo_sh);
  free(gaussians->c_sh);

}

void allocate_SALC(SALC *salc, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  salc->atm = (int *) malloc(salc->total_coef * sizeof(int));
  salc->memory +=  salc->total_coef * sizeof(int);
  if (salc->atm == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for salc in main() \n");
  MPI_Finalize();
  exit(1);
  }

  salc->irp = (int *) malloc(salc->num_salc * sizeof(int));
  salc->memory +=  salc->num_salc * sizeof(int);
  if (salc->irp == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for salc in main() \n");
  MPI_Finalize();
  exit(1);
  }

  salc->num_irp = (int *) malloc(symmetry->number_of_classes * sizeof(int));
  salc->memory +=  symmetry->number_of_classes * sizeof(int);
  if (salc->num_irp == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for salc in main() \n");
  MPI_Finalize();
  exit(1);
  }

  salc->num_coef = (int *) malloc(salc->num_salc * sizeof(int));
  salc->memory +=  salc->num_salc * sizeof(int);
  if (salc->num_coef == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for salc in main() \n");
  MPI_Finalize();
  exit(1);
  }

  AllocateDoubleMatrix(&salc->coeff,&salc->num_atom,&salc->total_coef,job);
  AllocateIntMatrix(&salc->bfn_posn,&salc->num_atom,&salc->total_coef,job);
  salc->memory +=  salc->num_atom * salc->total_coef * sizeof(double);
  salc->memory +=  salc->num_atom * salc->total_coef * sizeof(int);

}

void free_SALC(SALC *salc, JOB_PARAM *job)

{

  free(salc->atm);
  free(salc->irp);
  free(salc->num_irp);
  free(salc->num_coef);
  DestroyDoubleMatrix(&salc->coeff,job);
  DestroyIntMatrix(&salc->bfn_posn,job);

  job->memory[job->taskid] -= salc->memory;

}

void allocate_PAIR_TRAN(PAIR_TRAN *pair_p, ATOM *atoms, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file)

{

  int i, last_vector;
  int dim  = atoms->number_of_atoms_in_unit_cell * atoms->number_of_atoms_in_unit_cell * R_tables->last_vector;
  int dim1 = atoms->number_of_unique_atoms * atoms->number_of_unique_atoms;

  pair_p->cell1 = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->cell2 = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->latt1 = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->latt2 = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->lattd = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->k     = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->p     = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->off   = (int *) malloc((pair_p->tot + 1) * sizeof(int));
  pair_p->Off   = (int *) malloc((pair_p->nump + 1) * sizeof(int));
  pair_p->O     = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->P     = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->ptr   = (int *) malloc(dim * sizeof(int));
  pair_p->Ptr   = (int *) malloc(dim * sizeof(int));
  pair_p->uniq  = (int *) malloc(dim * sizeof(int));
  pair_p->posn  = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->numb  = (int *) malloc(pair_p->tot * sizeof(int));
  pair_p->K     = (int *) malloc(pair_p->nump * symmetry->number_of_operators * 8 * sizeof(int));
  pair_p->rot   = (int *) malloc(pair_p->nump * symmetry->number_of_operators * 2 * sizeof(int));
  pair_p->dist  = (double *) malloc(pair_p->tot * sizeof(double));

  if (pair_p->cell1 == NULL) { fprintf(file.err,"ERROR: There is not enough memory for pair_p"); exit(1); }

  pair_p->memory = ((12 * pair_p->tot + (10 * symmetry->number_of_operators + 1) * pair_p->nump + 3 * dim) * sizeof(int));

  job->memory[job->taskid] += pair_p->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"UNIQUE %5d TOTAL %5d PAIRS ALLOCATED %6d MEMORY ALLOCATED\n",pair_p->nump,pair_p->tot,pair_p->memory);

}

void free_PAIR_TRAN(PAIR_TRAN *pair_p, JOB_PARAM *job)

{

  job->memory[job->taskid] -= pair_p->memory;

  free(pair_p->cell1);
  free(pair_p->cell2);
  free(pair_p->latt1);
  free(pair_p->latt2);
  free(pair_p->lattd);
  free(pair_p->k);
  free(pair_p->p);
  free(pair_p->off);
  free(pair_p->Off);
  free(pair_p->O);
  free(pair_p->P);
  free(pair_p->ptr);
  free(pair_p->Ptr);
  free(pair_p->uniq);
  free(pair_p->posn);
  free(pair_p->numb);
  free(pair_p->K);
  free(pair_p->rot);
  free(pair_p->dist);

}

void allocate_TRIPLE_TRAN(TRIPLE_TRAN *triple, JOB_PARAM *job, FILES file)

{

  triple->numb  = (int *) malloc(triple->tot * sizeof(int));
  if (triple->numb == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->cell1 = (int *) malloc(triple->tot * sizeof(int));
  if (triple->cell1 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->cell2 = (int *) malloc(triple->tot * sizeof(int));
  if (triple->cell2 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->cell3 = (int *) malloc(triple->tot * sizeof(int));
  if (triple->cell3 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->latt1 = (int *) malloc(triple->tot * sizeof(int));
  if (triple->latt1 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->latt2 = (int *) malloc(triple->tot * sizeof(int));
  if (triple->latt2 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->latt3 = (int *) malloc(triple->tot * sizeof(int));
  if (triple->latt3 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->posn = (int *) malloc(triple->tot * sizeof(int));
  if (triple->posn == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->k     = (int *) malloc(triple->tot * sizeof(int));
  if (triple->k == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->p     = (int *) malloc(triple->tot * sizeof(int));
  if (triple->p == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for triple"); 
  MPI_Finalize();
  exit(1); 
 }

  triple->memory = 10 * triple->tot * sizeof(int);

  job->memory[job->taskid] += triple->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"TRIPLES ALLOCATED %6d MEMORY ALLOCATED %6d\n",triple->tot,triple->memory);

}

void free_TRIPLE_TRAN(TRIPLE_TRAN *triple, JOB_PARAM *job)

{

  job->memory[job->taskid] -= triple->memory;

  free(triple->numb);
  free(triple->cell1);
  free(triple->cell2);
  free(triple->cell3);
  free(triple->latt1);
  free(triple->latt2);
  free(triple->latt3);
  free(triple->posn);
  free(triple->k);
  free(triple->p);

}

void allocate_QUAD_TRAN(QUAD_TRAN *quad, JOB_PARAM *job, FILES file)

{

  quad->cell1 = (int *) malloc(quad->tot * sizeof(int));
  if (quad->cell1 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->cell2 = (int *) malloc(quad->tot * sizeof(int));
  if (quad->cell2 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->cell3 = (int *) malloc(quad->tot * sizeof(int));
  if (quad->cell3 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->cell4 = (int *) malloc(quad->tot * sizeof(int));
  if (quad->cell4 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->latt1 = (int *) malloc(quad->tot * sizeof(int));
  if (quad->latt1 == NULL) { 

  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->latt2 = (int *) malloc(quad->tot * sizeof(int));
  if (quad->latt2 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->latt3 = (int *) malloc(quad->tot * sizeof(int));
  if (quad->latt3 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->latt4 = (int *) malloc(quad->tot * sizeof(int));
  if (quad->latt4 == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->k     = (int *) malloc(quad->tot * sizeof(int));
  if (quad->k == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->p     = (int *) malloc(quad->tot * sizeof(int));
  if (quad->p == NULL) { 
  fprintf(file.err,"ERROR: There is not enough memory for quad"); 
  MPI_Finalize();
  exit(1); 
 }

  quad->memory = 10 * quad->tot * sizeof(int);

  job->memory[job->taskid] += quad->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"QUADS ALLOCATED %6d MEMORY ALLOCATED %6d\n",quad->tot,quad->memory);

}

void free_QUAD_TRAN(QUAD_TRAN *quad, JOB_PARAM *job)

{

  job->memory[job->taskid] -= quad->memory;

  free(quad->cell1);
  free(quad->cell2);
  free(quad->cell3);
  free(quad->cell4);
  free(quad->latt1);
  free(quad->latt2);
  free(quad->latt3);
  free(quad->latt4);
  free(quad->k);
  free(quad->p);

}

void allocate_fermi(FERMI *fermi, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  int i;

  fermi->occupied = (int *) malloc(job->spin_dim * fermi->nkunique * sizeof(int));
  if (fermi->occupied == NULL) { fprintf(file.out,"ERROR: There is not enough memory for fermi"); MPI_Finalize(); exit(1); }

  fermi->occupation = (double *) malloc(job->spin_dim * fermi->nkunique * atoms->number_of_sh_bfns_in_unit_cell * sizeof(double));
  if (fermi->occupation == NULL) { fprintf(file.out,"ERROR: There is not enough memory for fermi"); MPI_Finalize(); exit(1); }

  for (i = 0; i < job->spin_dim * fermi->nkunique * atoms->number_of_sh_bfns_in_unit_cell; i++)
  fermi->occupation[i] = k_zero;

  for (i = 0; i < job->spin_dim * fermi->nkunique; i++)
  fermi->occupied[i] = 0;

  fermi->memory =  job->spin_dim * fermi->nkunique * sizeof(int) + job->spin_dim * fermi->nkunique * atoms->number_of_sh_bfns_in_unit_cell * sizeof(double);

  job->memory[job->taskid] += fermi->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"%d KPOINTS IN FERMI %d MEMORY ALLOCATED\n", fermi->nkunique, fermi->memory);

}

void free_fermi(FERMI *fermi, JOB_PARAM *job)

{

  job->memory[job->taskid] -= fermi->memory;

  free(fermi->occupied);
  free(fermi->occupation);

}

void allocate_INT_1E(INT_1E *one_ints, int dim, int Function[], JOB_PARAM *job, FILES file)

{
 
  int i;
 
  one_ints->memory = 0;

  if (Function[0] == 1) {
   one_ints->Fock = (double *) malloc(dim * sizeof(double));
   for (i = 0; i < dim; i++) {
    one_ints->Fock[i] = k_zero;
   }
  one_ints->memory += dim * sizeof(double);
 }
  else one_ints->Fock = (double *) malloc(1);

  if (Function[1] == 1) {
   one_ints->Kinetic = (double *) malloc(dim * sizeof(double));
   for (i = 0; i < dim; i++) {
    one_ints->Kinetic[i] = k_zero;
   }
  one_ints->memory += dim * sizeof(double);
 }
  else one_ints->Kinetic = (double *) malloc(1);

  if (Function[2] == 1) {
   one_ints->ElecNuc = (double *) malloc(dim * sizeof(double));
   for (i = 0; i < dim; i++) {
    one_ints->ElecNuc[i] = k_zero;
   }
  one_ints->memory += dim * sizeof(double);
 }
  else one_ints->ElecNuc = (double *) malloc(1);

  if (Function[3] == 1) {
   one_ints->Momentum = (double *) malloc(3 * dim * sizeof(double));
   for (i = 0; i < 3 * dim; i++) {
    one_ints->Momentum[i] = k_zero;
   }
  one_ints->memory += 3 * dim * sizeof(double);
 }
  else one_ints->Momentum = (double *) malloc(1);

  if (Function[4] == 1) {
   one_ints->Overlap = (double *) malloc(dim * sizeof(double));
   for (i = 0; i < dim; i++) {
    one_ints->Overlap[i] = k_zero;
   }
  one_ints->memory += dim * sizeof(double);
 }
  else one_ints->Overlap = (double *) malloc(1);

  if (Function[5] == 1) {
   one_ints->Grad_Grad = (double *) malloc(9 * dim * sizeof(double));
   for (i = 0; i < 9 * dim; i++) {
    one_ints->Grad_Grad[i] = k_zero;
   }
  one_ints->memory += 9 * dim * sizeof(double);
 }
  else one_ints->Grad_Grad = (double *) malloc(1);

  if (Function[6] == 1) {
   one_ints->Dipole = (double *) malloc(3 * dim * sizeof(double));
   for (i = 0; i < 3 * dim; i++) {
    one_ints->Dipole[i] = k_zero;
   }
  one_ints->memory += 3 * dim * sizeof(double);
 }
  else one_ints->Dipole = (double *) malloc(1);

  if (Function[7] == 1) {
   one_ints->Coulomb = (double *) malloc(dim * sizeof(double));
   for (i = 0; i < dim; i++) {
    one_ints->Coulomb[i] = k_zero;
   }
  one_ints->memory += dim * sizeof(double);
 }
  else one_ints->Coulomb = (double *) malloc(1);

  job->memory[job->taskid] += one_ints->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"INT_1E %d MEMORY ALLOCATED\n", one_ints->memory);

}

void free_INT_1E(INT_1E *one_ints, int Function[], JOB_PARAM *job, FILES file)

{

  job->memory[job->taskid] -= one_ints->memory;

  if (Function[0] == 1)
  free(one_ints->Fock);

  if (Function[1] == 1)
  free(one_ints->Kinetic);

  if (Function[2] == 1)
  free(one_ints->ElecNuc);

  if (Function[3] == 1)
  free(one_ints->Momentum);

  if (Function[4] == 1)
  free(one_ints->Overlap);

  if (Function[5] == 1)
  free(one_ints->Grad_Grad);

  if (Function[6] == 1)
  free(one_ints->Dipole);

  if (Function[7] == 1)
  free(one_ints->Coulomb);

}

void allocate_integral_list(INTEGRAL_LIST *integral_list, int num, JOB_PARAM *job, FILES file)

{

  integral_list->i = (int *) malloc(num * sizeof(int));
  if (integral_list->i == NULL) { fprintf(file.out,"ERROR: There is not enough memory for integral_list"); MPI_Finalize(); exit(1); }

  integral_list->j = (int *) malloc(num * sizeof(int));
  if (integral_list->j == NULL) { fprintf(file.out,"ERROR: There is not enough memory for integral_list"); MPI_Finalize(); exit(1); }

  integral_list->k = (int *) malloc(num * sizeof(int));
  if (integral_list->k == NULL) { fprintf(file.out,"ERROR: There is not enough memory for integral_list"); MPI_Finalize(); exit(1); }

  integral_list->l = (int *) malloc(num * sizeof(int));
  if (integral_list->l == NULL) { fprintf(file.out,"ERROR: There is not enough memory for integral_list"); MPI_Finalize(); exit(1); }

  integral_list->gj = (int *) malloc(num * sizeof(int));
  if (integral_list->gj == NULL) { fprintf(file.out,"ERROR: There is not enough memory for integral_list"); MPI_Finalize(); exit(1); }

  integral_list->gk = (int *) malloc(num * sizeof(int));
  if (integral_list->gk == NULL) { fprintf(file.out,"ERROR: There is not enough memory for integral_list"); MPI_Finalize(); exit(1); }

  integral_list->gl = (int *) malloc(num * sizeof(int));
  if (integral_list->gl == NULL) { fprintf(file.out,"ERROR: There is not enough memory for integral_list"); MPI_Finalize(); exit(1); }

  integral_list->offset1 = (int *) malloc(num * sizeof(int));
  if (integral_list->offset1 == NULL) { fprintf(file.out,"ERROR: There is not enough memory for integral_list"); MPI_Finalize(); exit(1); }

  integral_list->offset2 = (int *) malloc(num * sizeof(int));
  if (integral_list->offset2 == NULL) { fprintf(file.out,"ERROR: There is not enough memory for integral_list"); MPI_Finalize(); exit(1); }

  integral_list->value = (double *) malloc(num * sizeof(double));
  if (integral_list->value == NULL) { fprintf(file.err,"ERROR: There is not enough memory for integral_list"); exit(1); }

  integral_list->memory = (9 * num * sizeof(int) + num * sizeof(double));
  job->memory[job->taskid] += integral_list->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"INT_2E %d MEMORY ALLOCATED\n", integral_list->memory);

}

void free_integral_list(INTEGRAL_LIST *integral_list, JOB_PARAM *job)

{

  job->memory[job->taskid] -= integral_list->memory;

  free(integral_list->i);
  free(integral_list->j);
  free(integral_list->k);
  free(integral_list->l);
  free(integral_list->gj);
  free(integral_list->gk);
  free(integral_list->gl);
  free(integral_list->offset1);
  free(integral_list->offset2);
  free(integral_list->value);

}

void allocate_dft_grid(DFT_GRID *dft_grid, int ngrid_points, ATOM *atoms, JOB_PARAM *job, FILES file)

{
 
  dft_grid->posn   = (int *)    malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (dft_grid->posn == NULL) { fprintf(file.out,"ERROR: There is not enough memory for dft_grid"); MPI_Finalize(); exit(1); }
  dft_grid->numb   = (int *)    malloc(atoms->number_of_atoms_in_unit_cell * sizeof(int));
  if (dft_grid->numb == NULL) { fprintf(file.out,"ERROR: There is not enough memory for dft_grid"); MPI_Finalize(); exit(1); }
  dft_grid->weight = (double *) malloc(ngrid_points * sizeof(double));
  if (dft_grid->weight == NULL) { fprintf(file.out,"ERROR: There is not enough memory for dft_grid"); MPI_Finalize(); exit(1); }
  dft_grid->r      = (double *) malloc(ngrid_points * sizeof(double));
  if (dft_grid->r == NULL) { fprintf(file.out,"ERROR: There is not enough memory for dft_grid"); MPI_Finalize(); exit(1); }
  dft_grid->y      = (double *) malloc(ngrid_points * sizeof(double));
  if (dft_grid->y == NULL) { fprintf(file.out,"ERROR: There is not enough memory for dft_grid"); MPI_Finalize(); exit(1); }
  dft_grid->x      = (VECTOR_DOUBLE *) malloc(ngrid_points * sizeof(VECTOR_DOUBLE));
  if (dft_grid->x == NULL) { fprintf(file.out,"ERROR: There is not enough memory for dft_grid"); MPI_Finalize(); exit(1); }

  dft_grid->memory = ngrid_points * (sizeof(VECTOR_DOUBLE) + 3 * sizeof(double)) + 2 * atoms->number_of_atoms_in_unit_cell * sizeof(int);
  job->max_memory[job->taskid] = dft_grid->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"dft_grid %d MEMORY ALLOCATED\n", dft_grid->memory);

}

void free_dft_grid(DFT_GRID *dft_grid, JOB_PARAM *job)

{

  free(dft_grid->r);
  free(dft_grid->x);
  free(dft_grid->y);
  free(dft_grid->weight);
  free(dft_grid->posn);
  free(dft_grid->numb);

  job->memory[job->taskid] -= dft_grid->memory;

}
 
void allocate_REAL_LATTICE(REAL_LATTICE *R, JOB_PARAM *job, FILES file)

{

  R->vec_ai = (VECTOR_DOUBLE *) malloc(R->max_vector * sizeof(VECTOR_DOUBLE));
  if (R->vec_ai == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for Rvec_ai in main() \n");
  MPI_Finalize(); 
  exit(1);;
  }

  R->vec_ai_ord = (int *) malloc(R->max_vector * sizeof(int));
  if (R->vec_ai_ord == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for R->vec_ai_ord in main() \n");
  MPI_Finalize(); 
  exit(1);;
  }

  R->vec_ai_inv = (int *) malloc(R->max_vector * sizeof(int));
  if (R->vec_ai_inv == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for R->vec_ai_inv in main() \n");
  MPI_Finalize(); 
  exit(1);;
  }

  R->num = (int *) malloc(R->max_vector * sizeof(int));
  if (R->num == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for R->num in main() \n");
  MPI_Finalize(); 
  exit(1);;
  }

  R->mag = (double *) malloc(R->max_vector * sizeof(double));
  if (R->mag == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for R->mag in main() \n");
  MPI_Finalize(); 
  exit(1);;
  }

  R->memory = 3 * R->max_vector * sizeof(int) + R->max_vector * sizeof(double) + R->max_vector * sizeof(VECTOR_DOUBLE);

  job->memory[job->taskid] += R->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"%d REAL LATTICE VECTORS ALLOCATED %d MEMORY ALLOCATED \n",R->max_vector,R->memory);

}

void free_REAL_LATTICE(REAL_LATTICE *R, JOB_PARAM *job)

{

  job->memory[job->taskid] -= R->memory;

  free(R->vec_ai);
  free(R->vec_ai_ord);
  free(R->vec_ai_inv);
  free(R->num);
  free(R->mag);

}

void allocate_REAL_LATTICE_TABLES(REAL_LATTICE_TABLES *R_tables, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  R_tables->sumvec = (int *) malloc(R_tables->margin_vector * R_tables->margin_vector * sizeof(int));
  if (R_tables->sumvec == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for sumvec in main() \n");
  MPI_Finalize(); exit(1);;
  }

  R_tables->diffvec = (int *) malloc(R_tables->margin_vector * R_tables->margin_vector * sizeof(int));
  if (R_tables->diffvec == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for diffvec in main() \n");
  MPI_Finalize(); exit(1);;
  }

  R_tables->lattvec = (int *) malloc(R_tables->margin_vector * symmetry->number_of_operators * sizeof(int));
  if (R_tables->lattvec == NULL) {
  fprintf(stderr, "ERROR: There is not enough memory for lattvec in main() \n");
  MPI_Finalize(); exit(1);;
  }

  R_tables->memory = 3 * R_tables->margin_vector * symmetry->number_of_operators * sizeof(int);
  job->memory[job->taskid] += R_tables->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"%d REAL_LATTICE_TABLES ALLOCATED %d MEMORY ALLOCATED\n",R_tables->margin_vector * R_tables->margin_vector,R_tables->memory);

}

void free_REAL_LATTICE_TABLES(REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job)

{

  job->memory[job->taskid] -= R_tables->memory;

  free(R_tables->sumvec);
  free(R_tables->diffvec);
  free(R_tables->lattvec);

}

void allocate_k_points(KPOINT_TRAN *knet, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int i;

  knet->bz  = (int *) malloc(knet->nktot * sizeof(int));
  knet->fbz = (int *) malloc(knet->nktot * sizeof(int));
  knet->opr = (int *) malloc(knet->nktot * sizeof(int));
  knet->trs = (int *) malloc(knet->nktot * sizeof(int));
  knet->ibz = (int *) malloc(knet->nktot * sizeof(int));
  knet->num = (int *) malloc(knet->unique * sizeof(int));
  knet->weight = (double *) malloc(knet->nktot * sizeof(double));
  knet->oblique = (VECTOR_INT *) malloc(knet->nktot * sizeof(VECTOR_INT));
  knet->cart = (VECTOR_DOUBLE *) malloc(knet->nktot * sizeof(VECTOR_DOUBLE));

   for (i = 0; i < knet->unique; i++)
    knet->num[i] = 0;
   for (i = 0; i < knet->nktot; i++)
    knet->trs[i] = -1;

  knet->memory = (4 * knet->nktot + 2 * knet->unique) * sizeof(int) + knet->unique * sizeof(double) + \
  knet->nktot * (sizeof(VECTOR_INT) + sizeof(VECTOR_DOUBLE)) ;
  job->memory[job->taskid] += knet->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"%d KPOINTS %d MEMORY ALLOCATED\n", knet->nktot,knet->memory);

}

void free_k_points(KPOINT_TRAN *knet, JOB_PARAM *job)

{

  job->memory[job->taskid] -= knet->memory;

  free(knet->bz);
  free(knet->fbz);
  free(knet->opr);
  free(knet->trs);
  free(knet->ibz);
  free(knet->num);
  free(knet->weight);
  free(knet->oblique);
  free(knet->cart);

}

void AllocateIntArray(int **a, int *n, JOB_PARAM *job)

{

  *a = (int *) malloc(*n * sizeof(int));
  if (*a == NULL) {
    fprintf(stderr, "ERROR: not enough memory for int vector! \n");
    exit(1);
   }

   job->memory[job->taskid] += *n * sizeof(int);
 
  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

}

void ResetIntArray(int *a,int *n)

{

  for (int i = 0; i < *n; i++) {
      a[i] = 0;
     }
}

void DestroyIntArray(int **a, int *n, JOB_PARAM *job)

{

   job->memory[job->taskid] -= *n * sizeof(int);

   free(*a);
  
}

void AllocateDoubleArray(double **a, int *n, JOB_PARAM *job)

{

  *a = (double *) malloc(*n * sizeof(double));
  if (*a == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double vector! \n");
    exit(1);
   }

   job->memory[job->taskid] += *n * sizeof(double);
 
  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

}

void ResetDoubleArray(double *a,int *n)

{

  for (int i = 0; i < *n; i++) {
      a[i] = k_zero;
     }
}

void DestroyDoubleArray(double **a, int *n, JOB_PARAM *job)

{

   job->memory[job->taskid] -= *n * sizeof(double);

   free(*a);
  
}

void AllocateIntMatrix(IntMatrix** g, int *iRows, int *iCols, JOB_PARAM *job)

{

  int i;

  *g=(IntMatrix *)malloc(sizeof(IntMatrix));
  if (*g==NULL){
    fprintf(stderr,"out of memory int matrix!\n");
    exit(1);
  }

  (*g)->a=(int **)malloc(*iRows*sizeof(int *));
  if ((*g)->a==NULL){
    fprintf(stderr,"out of memory for rows\n");
    exit(1);
  }

  (*g)->a[0]=(int *) malloc(*iRows * *iCols*sizeof(int));
  if ((*g)->a[0]==NULL){
    fprintf(stderr,"out of memory for columns\n");
    exit(1);
  }

  for(i=1;i<*iRows;i++){
    (*g)->a[i]=(*g)->a[i-1]+*iCols;
  }
  (*g)->iRows=*iRows;
  (*g)->iCols=*iCols;

  (*g)->memory = *iRows * *iCols * sizeof(int);

  job->memory[job->taskid] += (*g)->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

}

void ResetIntMatrix(IntMatrix* g)

{

  int i,j;

  for (i=0;i<g->iRows;i++){
    for(j=0;j<g->iCols;j++){
      g->a[i][j] = 0;
    }
  }
}

void DestroyIntMatrix(IntMatrix **g, JOB_PARAM *job)

{

  job->memory[job->taskid] -= (*g)->memory;

  free((*g)->a[0]);
  free((*g)->a);
  free((*g));

}

void AllocateDoubleMatrix(DoubleMatrix** g, int *iRows, int *iCols, JOB_PARAM *job)

{

  int i;
  long long memsize;

  *g=(DoubleMatrix *)malloc(sizeof(DoubleMatrix));
  if (*g==NULL){
    fprintf(stderr,"out of memory  double matrix!\n");
    exit(1);
  }

  (*g)->a=(double **)malloc(*iRows*sizeof(double *));
  if ((*g)->a==NULL){
    fprintf(stderr,"out of memory  for rows %3d\n",*iRows);
    exit(1);
  }

  memsize = *iRows * (*iCols * sizeof(double));

  (*g)->a[0]=(double *) malloc(memsize);
  if ((*g)->a[0]==NULL){
    fprintf(stderr,"out of memory for columns %d %d %lld\n",*iRows,*iCols,memsize);
    exit(1);
  }

  for(i=1;i<*iRows;i++){
    (*g)->a[i]=(*g)->a[i-1]+*iCols;
  }
  (*g)->iRows=*iRows;
  (*g)->iCols=*iCols;

  (*g)->memory = memsize;

  job->memory[job->taskid] += (*g)->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

}

void ResetDoubleMatrix(DoubleMatrix* g)

{

  int i,j;

  for (i=0;i<g->iRows;i++){
    for(j=0;j<g->iCols;j++){
      g->a[i][j]=k_zero;
    }
  }
}

void DestroyDoubleMatrix(DoubleMatrix **g, JOB_PARAM *job)

{

  job->memory[job->taskid] -= (*g)->memory;

  free((*g)->a[0]);
  free((*g)->a);
  free((*g));

}

void AllocateComplexMatrix(ComplexMatrix** g, int *iRows, int *iCols, JOB_PARAM *job)

{

  int i,j;

  *g = (ComplexMatrix *) malloc(sizeof(ComplexMatrix));
  if (*g == NULL){
    fprintf(stderr,"out of memory  complex matrix!\n");
    exit(1);
  }

  (*g)->a = (Complex **)malloc(*iRows*sizeof(Complex *));
  if ((*g)->a == NULL){
    fprintf(stderr,"out of memory  for rows\n");
    exit(1);
  }

  (*g)->a[0] = (Complex *) malloc(*iRows * *iCols*sizeof(Complex));
  if ((*g)->a[0] == NULL){
    fprintf(stderr,"out of memory for columns\n");
    exit(1);
  }

  for(i = 1;i < *iRows; i++){
    (*g)->a[i] = (*g)->a[i - 1] + *iCols;
  }

  (*g)->iRows = *iRows;
  (*g)->iCols = *iCols;

  (*g)->memory = *iRows * *iCols * sizeof(Complex);

  job->memory[job->taskid] += (*g)->memory;

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

}

void ResetComplexMatrix(ComplexMatrix* g)

{

int i,j;

  for (i = 0; i < g->iRows; i++){
    for(j = 0; j < g->iCols; j++){
      g->a[i][j] = Complex(k_zero, k_zero);
     }
    }

}

void DestroyComplexMatrix(ComplexMatrix **g, JOB_PARAM *job)

{

  job->memory[job->taskid] -= (*g)->memory;

  free((*g)->a[0]);
  free((*g)->a);
  free((*g));

}

void AllocateComplexArray(Complex **a, int *n, JOB_PARAM *job)

{

  *a = (Complex *) malloc(*n * sizeof(Complex));
  if (*a == NULL) {
    fprintf(stderr, "ERROR: not enough memory for Complex vector! \n");
    exit(1);
  }

  job->memory[job->taskid] += *n * sizeof(Complex);

  if (job->max_memory[job->taskid] < job->memory[job->taskid]) 
  job->max_memory[job->taskid] = job->memory[job->taskid];

}

void ResetComplexArray(Complex *a,int *n) {

  for (int i = 0; i < *n; i++) {
      a[i] = Complex(k_zero, k_zero);
    }

}

void DestroyComplexArray(Complex **a, int *n, JOB_PARAM *job)

{
 
  job->memory[job->taskid] -= *n * sizeof(Complex);

  free((*a));

}

void array_dimensions(int *dim, int *dim1, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file)

{
  int i, q, tmp;

  *dim = 0;
  *dim1 = 0;

   for (i = 0; i < pair_p->nump; i++) {
    q = pair_p->posn[i];
    tmp = atoms->bfnnumb[pair_p->cell1[q]] * atoms->bfnnumb[pair_p->cell2[q]] ;
    *dim  += tmp ;
    *dim1 += tmp * pair_p->numb[i] ;
   }

   if (job->taskid == 0 && job->verbosity > 1)
   fprintf(file.out,"Matrix dimension reduced %d full %d\n", *dim, *dim1);

}

void sh_array_dimensions(int *dim, int *dim1, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file)

{
  int i, q, tmp;

  *dim = 0;
  *dim1 = 0;

   for (i = 0; i < pair_p->nump; i++) {
    q = pair_p->posn[i];
    tmp = atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]] ;
    *dim  += tmp ;
    *dim1 += tmp * pair_p->numb[i] ;
   }

    job->dimp = *dim;
    job->dimf = *dim1;

   if (job->taskid == 0 && job->verbosity > 1)
   fprintf(file.out,"Reduced (Full) sh matrix dimension %d (%d)\n",*dim,*dim1);

}

void triple_array_dimensions(int *dim, int *dim1, TRIPLE_TRAN *triple, ATOM *atoms, JOB_PARAM *job, FILES file)

{
  int i, q, tmp;

  *dim = 0;
  *dim1 = 0;

   for (i = 0; i < triple->nump; i++) {
    q = triple->posn[i];
    tmp = atoms->bfnnumb[triple->cell1[q]] * atoms->bfnnumb[triple->cell2[q]] * atoms->bfnnumb[triple->cell3[q]];
    *dim  += tmp ;
    *dim1 += tmp * triple->numb[i] ;
   }

   if (job->taskid == 0 && job->verbosity > 1)
   fprintf(file.out,"Triples matrix dimension reduced %d full %d\n", *dim, *dim1);

}

void sh_triple_array_dimensions(int *dim, int *dim1, TRIPLE_TRAN *triple, ATOM *atoms, JOB_PARAM *job, FILES file)

{
  int i, q, tmp;

  *dim = 0;
  *dim1 = 0;

   for (i = 0; i < triple->nump; i++) {
    q = triple->posn[i];
    tmp = atoms->bfnnumb_sh[triple->cell1[q]] * atoms->bfnnumb_sh[triple->cell2[q]] * atoms->bfnnumb_sh[triple->cell3[q]] ;
    *dim  += tmp ;
    *dim1 += tmp * triple->numb[i] ;
   }

   if (job->taskid == 0 && job->verbosity > 1)
   fprintf(file.out,"Reduced (Full) sh triple matrix dimension %d (%d)\n",*dim,*dim1);

}
