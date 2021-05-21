  /*! \mainpage Exciton08 Manual
 this is what appears on the main page of the Exciton08 manual
 you can create sections
 \section first text 1
 just another thing
 \section second text 2
 is quite simple to use inline mathematics with latex \f$ \alpha=0\f$
 or to write proper formulae
 \f[
 E=mc^2
 \f]
 more info here http://www.stack.nl/~dimitri/doxygen/manual.html
 */

/*! \file MAIN.cpp
 \brief main file
 \details a proper description
 should be added later
 \remarks  here just remarks

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include "mycomplex.h"
*/
#include <mpi.h>
#include <cstring>
#include "conversion_factors.h"
#include "myconstants.h"
#include "mylogical.h"
#include "USER_DATA.h"
#include "TOOLS.h"
#include "HEADER.h"
#include "ALLOCATE_MEMORY.h"
#include "SETUP_SYMMETRY.h"
#include "SETUP_ATOMS.h"
#include "SETUP_CRYSTAL.h"
#include "SETUP_REAL_LATTICE.h"
#include "SETUP_RECIPROCAL_LATTICE.h"
#include "IVANIC_RUEDENBERG.h"
#include "INPUT_ALL.h"
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
  int read_crystal;
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

  read_crystal = 0;
  job.mxr = 10; // real space lattice is 10 x 10 x 10 for C, 10 x 10 for S 10 for P
  job.C09 = 0;
  //job.mxr = 1;
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

  if (job.taskid == 0)
      fprintf(file.out,"SPACE GROUP SYMMETRY SWITCHED OFF\n\n");

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

        // *****JOB: LATTICE SUM LIMITS *************************************************** /

    if (!strcmp(jobname, "LATTICE_SUM_LIMITS")) {

  if (job.taskid == 0 && job.verbosity > 1)
      fprintf(file.out,"LATTICE SUM LIMIT DEFAULT VALUES\nREAL SPACE LIMIT (Angs)    %4.1lf\nRECIPROCAL SPACE LIMIT (1/Bohr) %4.1lf\n\n",R.cutoff, G.cutoff);

      read_line(file.job, title, 99);
      sscanf(title, "%lf %lf", &R.cutoff, &G.cutoff);

  if (job.taskid == 0 && job.verbosity > 1)
      fprintf(file.out,"LATTICE SUM LIMITS RESET TO\nREAL SPACE LIMIT (Angs)    %4.1lf\nRECIPROCAL SPACE LIMIT (1/Bohr) %4.1lf\n\n",R.cutoff, G.cutoff);

    }

        // *****JOB: READ FROM CRYSTAL09 OUTPUT ****************************************

    if (!strcmp(jobname, "READ_CRYSTAL_09")) {

  if (job.taskid == 0 && job.trs == 0)
      fprintf(file.out,"MUST USE TIME REVERSAL SYMMETRY WHEN READING FROM CRYSTAL09\n\n");
  if (job.trs == 0) {
      MPI_Finalize();
      exit(0);
     }

  if (job.taskid == 0 && job.verbosity > 1)
      fprintf(file.out,"READING ION POSITIONS FROM CRYSTAL09 OUTPUT FILE\n\n");
      crystal.type[1] = '1';
      read_crystal = 1;
      job.C09 = 1;     // Read geometry from INPUT by default

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

  conventional_unit_cell(&crystal,file);
  primitive_unit_cell(&crystal,&symmetry,&job,file);
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

  R.max_vector = (2 * job.mxr + 1) * (2 * job.mxr + 1) * (2 * job.mxr + 1);
  G.cutoff = 10.0;

  break;

  case 'S':

  R.max_vector = (2 * job.mxr + 1) * (2 * job.mxr + 1);
  G.cutoff = 10.0;

  break;

  case 'P':

  R.max_vector = (2 * job.mxr + 1);
  G.cutoff = 10.0;

  break;

  case 'M':
 
  job.mxr = 1;
  G.max_vector = 1;
  G.last_vector = 1;
  G.cutoff = k_zero;
  allocate_RECIPROCAL_LATTICE(&G,&job,file);
  R.max_vector =  1;
  R.cutoff = 9999.0;
  allocate_REAL_LATTICE(&R,&job,file);
  R.vec_ai[0].comp1 = k_zero;
  R.vec_ai[0].comp2 = k_zero;
  R.vec_ai[0].comp3 = k_zero;
  R.last_ewald_vector = 1;

  break;

  } // end switch

  count_reciprocal_lattice_vectors(crystal,&G,&job,file);
  allocate_RECIPROCAL_LATTICE(&G,&job,file);
  generate_reciprocal_lattice(crystal,&G,&job,file);

  R.cutoff /= bohr_to_AA;
  allocate_REAL_LATTICE(&R,&job,file);
  generate_real_lattice(&crystal,&R,&G,&job,file);

    switch (read_crystal) {

    case 0:

  count_all_atoms(&crystal,&atoms,atms,&symmetry,&job,file);
  atoms_ax.number_of_atoms_in_unit_cell = atoms.number_of_atoms_in_unit_cell;
  allocate_ATOM(&atoms,&job,file);
  allocate_ATOM(&atoms_ax,&job,file);
  allocate_ATOM_TRAN(&atom_p,&atoms,&symmetry,&job,file);
  allocate_ATOM_TRAN(&atom_i,&atoms,&symmetry,&job,file);
  generate_all_atoms(&crystal,&R,&symmetry,atms,basis,&atoms,&atom_p,&atom_i,&job,file);
  generate_all_atoms(&crystal,&R,&symmetry,atms,basis_ax,&atoms_ax,&atom_p,&atom_i,&job,file);

    break;

    case 1:

  count_all_crystal_atoms(&atoms,&job,file);
  atoms_ax.number_of_atoms_in_unit_cell = atoms.number_of_atoms_in_unit_cell;
  allocate_ATOM(&atoms,&job,file);
  allocate_ATOM(&atoms_ax,&job,file);
  allocate_ATOM_TRAN(&atom_p,&atoms,&symmetry,&job,file);
  allocate_ATOM_TRAN(&atom_i,&atoms,&symmetry,&job,file);
  generate_all_crystal_atoms(&crystal,&R,&symmetry,basis,&atoms,&atom_p,&atom_i,&job,file);
  generate_all_crystal_atoms(&crystal,&R,&symmetry,basis_ax,&atoms_ax,&atom_p,&atom_i,&job,file);

  if (job.taskid == 0 && job.verbosity > 1)
  fprintf(file.out,"READ ION POSITIONS FROM CRYSTAL OUTPUT FILE\n\n");

  break;

  }

  if (job.sgs == 0) {
    for (i = 0; i < atoms.number_of_atoms_in_unit_cell; i++) {
      atom_p.O[i] = 0; 
      atom_p.K[i] = i; 
     }
    symmetry.number_of_operators = 1;
   }

  PAIR_TRAN pair_p;
  count_pairs4(&pair_p,&atoms,&atom_p,&symmetry,&R,&job,file);
  generate_real_lattice_shells(&crystal,&R,&job,file);
  R_tables.last_vector = R.last_vector;
  R_tables.margin_vector = R.margin_vector;
  R_tables.max_vector = R.max_vector;
  allocate_REAL_LATTICE_TABLES(&R_tables,&symmetry,&job,file);
  R_tables.diffvec[0] = 0;
  R_tables.sumvec[0] = 0;
  for (i = 0; i < symmetry.number_of_operators; i++)
  R_tables.lattvec[i] = 0;
  generate_real_lattice_tables(&crystal,&symmetry,&R,&R_tables,&job,file);
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
  free_RECIPROCAL_LATTICE(&G,&job);
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
