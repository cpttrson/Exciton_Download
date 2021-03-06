
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
#include "conversion_factors.h"
#include "myconstants.h"
#include "mylogical.h"
#include "USER_DATA.h"
#include "TOOLS.h"
#include "HEADER.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "SETUP_SYMMETRY.h"
#include "SETUP_ATOMS.h"
#include "IVANIC_RUEDENBERG.h"
#include "INPUT_MOLECULE.h"
#include "PAIRS_QUADS.h"

using namespace std;

int main(int argc, char *argv[]) {

  if (argc < 3) {
  printf("Usage: exciton input_file output_file\n");
  exit(1);
  }

int rank;

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

int i;
int print;
char jobname[17];
char title[100];

  // ******************************************************************************************
  // * BEGIN Initialisation Section                                                           *
  // ******************************************************************************************
 
FILES   file ;
JOB_PARAM job;
CRYSTAL crystal ;
REAL_LATTICE R;
REAL_LATTICE_TABLES R_tables;
RECIPROCAL_LATTICE G;
SYMMETRY symmtemp ;
SYMMETRY symmetry ;
ATOM atoms, atoms_ax;
ATOM_TRAN atom_p;
ATOM_TRAN atom_i;
SHELL shells, shells_ax;
GAUSSIAN gaussians, gaussians_ax;

  // ******************************************************************************************
  // * JOB_PARAM structure contains switches to control job flow and monitors memory use      *
  // ******************************************************************************************
 
  file.job = fopen(argv[1], "r");

  MPI_Comm_size(MPI_COMM_WORLD, &job.numtasks);
  job.taskid = rank;

  allocate_JOB_PARAM(&job,file);

  char xx[4], yy[4] = { '.' };
  sprintf(xx, "%d", rank);
  strcat(yy, xx);

  file.out = fopen(strcat(argv[2], yy), "w");
  if (file.out == NULL) {
  fprintf(file.out, "Cannot open output file. Job terminating.");
  MPI_Finalize();
  exit(1);
 }

  file.job = fopen(argv[1], "r");
  if (file.job == NULL) {
  fprintf(file.out, "Cannot open input file. Job terminating.");
  MPI_Finalize();
  exit(1);
 }

  file.in = fopen("INPUT", "r");
  if (file.in == NULL) {
  fprintf(file.out, "Cannot open file 'INPUT'. Job terminating.");
  MPI_Finalize();
  exit(1);
 }

  if (job.taskid == 0)
  print_header(file);

  job.mxr = 1;
  job.mxg = 1;
  job.sym_adapt = 1;
  job.kss = 1;
  job.rss = 1;
  job.sgs = 1;
  job.pms = 1;
  job.trs = 1;
  job.l_max = 2;
  job.verbosity = 1;
  symmtemp.number_of_operators = 48;
  allocate_SYMMETRY(&symmtemp,&job,file);

  print = 0;
  for (i = 0; i < 3; i++) {
  crystal.conventional_cell[i].comp1 = k_zero;
  crystal.conventional_cell[i].comp2 = k_zero;
  crystal.conventional_cell[i].comp3 = k_zero;
 }
  read_symmetry_group(&crystal,&symmtemp,print,&job,file);
  count_unique_atoms(&atoms,&job,file);
  ATOMTYPE atms[atoms.number_of_unique_atoms];
  read_unique_atoms(&atoms,atms,&crystal,&job,file);
  count_basis_sets(&atoms,0,&job,file);
  rewind(file.in);

  print = 0;
  read_symmetry_group(&crystal,&symmtemp,print,&job,file);
  count_unique_atoms(&atoms_ax,&job,file);
  read_unique_atoms(&atoms_ax,atms,&crystal,&job,file);
  count_basis_sets(&atoms_ax,1,&job,file);
  rewind(file.in);

  do {

     read_line(file.job, title, 99);
     sscanf(title, "%s", jobname);

    // *****JOB: SET MXR ********************************************************** 

    if (!strcmp(jobname, "SET_MXR")) {

      read_line(file.job, title, 99);
      sscanf(title, "%d", &job.mxr);

    if (job.taskid == 0 && job.verbosity > 1)
      fprintf(file.out,"MXR SET TO %d\n\n",job.mxr);

    }

     // *****JOB: SET MXG ********************************************************** 

    if (!strcmp(jobname, "SET_MXG")) {

      read_line(file.job, title, 99);
      sscanf(title, "%d", &job.mxg);

    if (job.taskid == 0 && job.verbosity > 1)
      fprintf(file.out,"MXG SET TO %d\n\n",job.mxg);

    }

    // *****JOB: PRINTING LEVEL *************************************************** 

    if (!strcmp(jobname, "PRINTING_LEVEL")) {

      read_line(file.job, title, 99);
      sscanf(title, "%d", &job.verbosity);

    if (job.taskid == 0 && job.verbosity > 1)
      fprintf(file.out,"PRINTING LEVEL SET TO %d\n\n",job.verbosity);

    }

    // *****JOB: PRINT PAIRS ****************************************************** 

    if (!strcmp(jobname, "PRINT_PAIRS")) {

      job.print_pairs = 1;

    if (job.taskid == 0 && job.verbosity > 1)
      fprintf(file.out,"PRINTING PAIRS and QUADS\n\n");

    }

    // *****JOB: SCF COULOMB INTEGRALS **********************************************

    if (!strcmp(jobname, "COU_OFF")) {

      job.scf_coulomb = 0;

    if (job.taskid == 0)
      fprintf(file.out,"2E COULOMB INTEGRALS SWITCHED OFF\n\n");

    }

    // *****JOB: SCF EXCHANGE INTEGRALS **********************************************

    if (!strcmp(jobname, "EXC_OFF")) {

      job.scf_exchange = 0;

    if (job.taskid == 0)
      fprintf(file.out,"2E EXCHANGE INTEGRALS SWITCHED OFF\n\n");

    }

    // *****JOB: SYMMETRY ADAPTATION **********************************************

    if (!strcmp(jobname, "SYM_ADAPT_OFF")) {

      job.sym_adapt = 0;

    if (job.taskid == 0)
      fprintf(file.out,"SYMMETRY ADAPTED LINEAR COMBINATIONS OF BASIS FUNCTIONS SWITCHED OFF\n\n");

    }

    // *****JOB: PERMUTATION SYMMETRY ********************************************** 

    if (!strcmp(jobname, "PMS_OFF")) {

      job.pms = 0;

    if (job.taskid == 0 && job.verbosity > 1)
      fprintf(file.out,"PERMUTATION SYMMETRY SWITCHED OFF\n\n");

    }

    // *****JOB: TIME REVERSAL SYMMETRY ********************************************* 

    if (!strcmp(jobname, "TRS_OFF")) {

      job.trs = 0;

  if (job.taskid == 0)
      fprintf(file.out,"TIME REVERSAL SYMMETRY SWITCHED OFF\n\n");

    }

    // *****JOB: SPACE GROUP SYMMETRY ********************************************** 

    if (!strcmp(jobname, "SGS_OFF")) {

      job.sgs = 0;
      job.sym_adapt = 0;

  if (job.taskid == 0)
      fprintf(file.out,"SPACE GROUP SYMMETRY AND SYMMETRY ADAPTATION SWITCHED OFF\n\n");

    }

    // *****JOB: REAL SPACE SYMMETRY **********************************************

    if (!strcmp(jobname, "RSS_OFF")) {

      job.rss = 0;

  if (job.taskid == 0)
      fprintf(file.out,"REAL SPACE SYMMETRY SWITCHED OFF\n\n");

    }

    // *****JOB: K SPACE SYMMETRY ************************************************** 

    if (!strcmp(jobname, "KSS_OFF")) {

      job.kss = 0;

  if (job.taskid == 0)
      fprintf(file.out,"RECIPROCAL SPACE SYMMETRY SWITCHED OFF\n\n");

    }

  } while (strcmp(jobname, "BEGIN")); fflush(file.out);

  // ******************************************************************************************
  // * END Initialisation Section                                                             *
  // ******************************************************************************************
 
  symmetry.number_of_operators = symmtemp.number_of_operators ;
  symmetry.number_of_permutations = 8;

  allocate_SYMMETRY(&symmetry,&job,file);

  for(i=0;i < 9 * symmetry.number_of_operators;i++) {
  symmetry.irr[i] = symmtemp.irr[i];}
  for(i=0;i < symmetry.number_of_operators;i++) {
  symmetry.taur[i] = symmtemp.taur[i];}
  for(i=0;i < 9 * symmetry.number_of_operators;i++) {
  symmetry.inr[i] = symmtemp.inr[i];}

  //conventional_unit_cell(&crystal,file);
  //primitive_unit_cell(&crystal,&symmetry,&job,file);
  generate_cartesian_symmetry_operators(&crystal,&symmetry,&job,file);
  print = 0;
  int flag;
  flag = 0;
  ATOMTYPE basis[atoms.number_of_basis_sets + 1];
  ATOMTYPE basis_ax[atoms_ax.number_of_basis_sets + 1];
  rewind(file.in);

  read_symmetry_group(&crystal,&symmtemp,print,&job,file);
  count_unique_atoms(&atoms,&job,file);
  read_unique_atoms(&atoms,atms,&crystal,&job,file);
  read_basis_sets(&atoms,flag,basis,&job,file);
  rewind(file.in);
  read_symmetry_group(&crystal,&symmtemp,print,&job,file);
  count_unique_atoms(&atoms_ax,&job,file);
  read_unique_atoms(&atoms_ax,atms,&crystal,&job,file);
  flag = 1;
  read_basis_sets(&atoms_ax,flag,basis_ax,&job,file);

  free_SYMMETRY(&symmtemp,&job);

  switch (crystal.type[0]) {

  case 'C':
  case 'S':
  case 'P':

  break;

  case 'M':
  
  G.max_vector = 1;
  G.last_vector = 1;
  G.cutoff = k_zero;
  //allocate_RECIPROCAL_LATTICE(&G,&job,file);
  R.max_vector =  1;
  R.cutoff = 9999.0;
  allocate_REAL_LATTICE(&R,&job,file);
  R.vec_ai[0].comp1 = k_zero;
  R.vec_ai[0].comp2 = k_zero;
  R.vec_ai[0].comp3 = k_zero;
  R.last_ewald_vector = 1;

  break;

  } // end switch

  count_all_atoms(&crystal,&atoms,atms,&symmetry,&job,file);
  atoms_ax.number_of_atoms_in_unit_cell = atoms.number_of_atoms_in_unit_cell;
  allocate_ATOM(&atoms,&job,file);
  allocate_ATOM(&atoms_ax,&job,file);
  allocate_ATOM_TRAN(&atom_p,&atoms,&symmetry,&job,file);
  allocate_ATOM_TRAN(&atom_i,&atoms,&symmetry,&job,file);
  generate_all_atoms(&crystal,&R,&symmetry,atms,basis,&atoms,&atom_p,&atom_i,&job,file);
  generate_all_atoms(&crystal,&R,&symmetry,atms,basis_ax,&atoms_ax,&atom_p,&atom_i,&job,file);

  if (job.sgs == 0) {
    for (i = 0; i < atoms.number_of_atoms_in_unit_cell; i++) {
      atom_p.O[i] = 0; 
      atom_p.K[i] = i; 
     }
    symmetry.number_of_operators = 1;
   }

  if (job.sym_adapt == 0) symmetry.number_of_classes = 1;

  PAIR_TRAN pair_p;
  count_pairs4(&pair_p,&atoms,&atom_p,&symmetry,&R,&job,file);

  //generate_real_lattice_shells(&crystal,&R,&job,file);
  R_tables.last_vector = R.last_vector;
  R_tables.margin_vector = R.margin_vector;
  R_tables.max_vector = R.max_vector;
  allocate_REAL_LATTICE_TABLES(&R_tables,&symmetry,&job,file);
  R_tables.diffvec[0] = 0;
  R_tables.sumvec[0] = 0;
  for (i = 0; i < symmetry.number_of_operators; i++)
  R_tables.lattvec[i] = 0;
  allocate_SHELL_GAUSSIAN(&shells,&gaussians,&atoms,&job,file);
  allocate_SHELL_GAUSSIAN(&shells_ax,&gaussians_ax,&atoms_ax,&job,file);
  generate_rotation_operators_ivanic_ruedenberg(&symmetry,&job,file);
  generate_runtime_arrays(&atoms,&atom_p,basis,&shells,&gaussians,&job,file);
  generate_runtime_arrays(&atoms_ax,&atom_p,basis_ax,&shells_ax,&gaussians_ax,&job,file);
  test_rotation_operators(&atoms,&shells,&symmetry,&job,file);
  count_density_pairs4(&pair_p, &atoms, &atom_p, &symmetry, &R, &R_tables, &job, file);
  allocate_PAIR_TRAN(&pair_p,&atoms,&symmetry,&R_tables,&job,file);
  generate_density_pairs4(&pair_p,&atoms,&atom_p,&symmetry,&R,&R_tables,&job,file); 
  print_pairs(&pair_p,&atoms,&R,&job,file);

  startjob(&atoms,&atoms_ax,&atom_p,&shells,&gaussians,&shells_ax,&gaussians_ax,&crystal,&symmetry,&R,&R_tables,&G,&job,file);

  free_PAIR_TRAN(&pair_p,&job);
  free_SHELL_GAUSSIAN(&shells,&gaussians,&job);
  free_SHELL_GAUSSIAN(&shells_ax,&gaussians_ax,&job);
  free_REAL_LATTICE_TABLES(&R_tables,&job);
  free_ATOM_TRAN(&atom_i,&job);
  free_ATOM_TRAN(&atom_p,&job);
  free_ATOM(&atoms,&job);
  free_ATOM(&atoms_ax,&job);
  free_REAL_LATTICE(&R,&job);
  free_SYMMETRY(&symmetry,&job);
  //printf("residual memory %5d\n",job.memory[job.taskid]);
  //printf("max memory %5d\n",job.max_memory[job.taskid]);
  free_JOB_PARAM(&job,file);

  MPI_Finalize();

  fclose(file.in);
  //if (job.taskid == 0)
  fclose(file.out);
  fclose(file.job);

  return 0;

}
