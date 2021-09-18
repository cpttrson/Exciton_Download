

  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

#include <cstring>
#include <cstdlib>
#include <mkl_scalapack.h>
#include "mycomplex.h"
#include "conversion_factors.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "myconstants.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "PARALLEL.h"
#include "PAIRS_QUADS.h"
#include "ROTATIONS_MOLECULE.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "INTEGRALS_2C_MOLECULE.h"
#include "INTEGRALS_3C_MOLECULE.h"
#include "PRINT_MOLECULE.h"
#include "SCALAPACK.h"
#include "DENSITY_FITTING_MOLECULE.h"
#include "GW_BSE_MOLECULE.h"
extern "C" void  Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern "C" int   Cblacs_pnum(int, int, int);

using namespace std;

void gw_molecule(FERMI* fermi, ATOM* atoms, ATOM_TRAN* atom_p, SHELL* shells, GAUSSIAN* gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int i, j1, I1, I2, il, il1, jl, jl1, mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int buffer_size1, buffer_size2;
int ictxt, nbsize; 
int block_cyclic;
int *dim_send, *offset_j;
int begin_j[job->numtasks], end_j[job->numtasks];
int num_proc = job->numtasks < dim1 ? job->numtasks : dim1;
int offset;
double time1, time2;
double time3, time4;
double *dSigma_dE;
double *Ham_buffer1, *Ham_buffer2;
double *integral_buffer1, *integral_buffer2;
DoubleMatrix *V_inv;
DoubleMatrix *Sigma;

  time1 = MPI_Wtime();

  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in gw. HOMO level not in range of bands. ntransitions = %d %d %d\n",ntransitions,nocc,nvir);
  MPI_Finalize();
  exit(0);
 }

  initialise_spk_grid(&ntransitions, &ictxt, &nbsize, job, file);
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, &nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, &nbsize, &mycol, &izero, &npcol);
  if (job->taskid == 0) printf("mpA %7d nqA %7d blocksize %7d ntransitions %7d\n",mpA,nqA,nbsize,ntransitions);

  // ******************************************************************************************
  // * Generate three centre integrals needed for matrix elements                             *
  // * Generate inverse of Coulomb matrix (if needed)                                         *
  // ******************************************************************************************
  
  AllocateIntArray(&dim_send,&job->numtasks,job);
  count_integral_buffer_sizes(dim_send,fermi,atoms_ax,job,file);
  AllocateDoubleArray(&integral_buffer1,&dim_send[job->taskid],job);
  AllocateDoubleArray(&integral_buffer2,&dim_send[job->taskid],job);
  generate_integral_buffers_molecule_ija(integral_buffer1,integral_buffer2,&ictxt,&nbsize,fermi,atoms,atom_p,shells, \
  gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  MPI_Barrier(MPI_COMM_WORLD);

  // ******************************************************************************************
  // * Generate inverse of Coulomb matrix                                                     *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  AllocateDoubleMatrix(&V_inv,&dim1ax,&dim1ax,job);
  generate_coulomb_matrix_inverse(V_inv,atom_p,atoms,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("generate_coulomb_matrix_inverse            %10.2f\n",time4);

  // ******************************************************************************************
  // * Allocate memory for BSE eigenvectors and eigenvalues                                   *
  // ******************************************************************************************

  buffer_size1 = mpA * nqA;
  buffer_size2 = 1;
  AllocateDoubleArray(&Ham_buffer1,&buffer_size1,job);
  AllocateDoubleArray(&Ham_buffer2,&buffer_size2,job);
  ResetDoubleArray(Ham_buffer1,&buffer_size1);
  ResetDoubleArray(Ham_buffer2,&buffer_size2);

  // ******************************************************************************************
  // * Generate and diagonalise Casida matrix                                                 *
  // ******************************************************************************************
  
  bse_hamiltonian(&ictxt,&nbsize,Ham_buffer1,Ham_buffer2,integral_buffer1,integral_buffer2,fermi,atom_p,atoms,\
  atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  DestroyIntArray(&dim_send,&job->numtasks,job);
  time3 = MPI_Wtime();
  diagonalise_rpa_hamiltonian(Ham_buffer1,&ictxt,&nbsize,fermi,atoms,atoms_ax,job,file);
  DestroyDoubleArray(&Ham_buffer1,&buffer_size1,job);
  DestroyDoubleArray(&Ham_buffer2,&buffer_size2,job);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("diagonalise_rpa_hamiltonian                %10.2f\n",time4);

  // ******************************************************************************************
  // * Calculate self energy and quasi-particle energy shifts                                 *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  AllocateIntArray(&offset_j,&dim1,job);
  mpi_begin_end(begin_j, end_j, dim1, num_proc, job, file);
  for (i = num_proc; i < job->numtasks; i++) begin_j[i] = 0;
  for (i = num_proc; i < job->numtasks; i++) end_j[i] = 0;
  for (i = 0; i < dim1; i++) offset_j[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1] * nbands * nbands; }}
  self_energy_diagonal(&ictxt,&nbsize,begin_j,end_j,offset_j,integral_buffer1,integral_buffer2,\
  atoms_ax,fermi,job,file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("self_energy                                %10.2f\n",time4);
  DestroyDoubleArray(&integral_buffer1,&dim_send[job->taskid],job);
  DestroyDoubleArray(&integral_buffer2,&dim_send[job->taskid],job);
  DestroyIntArray(&offset_j,&dim1,job);
  DestroyDoubleMatrix(&V_inv,job);
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("gw                                         %10.2f\n",time2);

}

void bse_molecule(FERMI* fermi, ATOM* atoms, ATOM_TRAN* atom_p, int *numfrag, int *natom, int nat[][2], SHELL* shells, GAUSSIAN* gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates BSE and BSE-TDA excitations                                         *
  // ******************************************************************************************
  
  // ******************************************************************************************
  // * Allocate fermi structure                                                               *
  // ******************************************************************************************

int i, j;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int tda, count;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int ictxt, nbsize, cblacs_taskid, itemp, nprow, npcol, myrow, mycol, mpA, nqA, izero = 0, ione = 1;
int descA[9];
int info = 0;
int buffer_size1, buffer_size2;
int *dim_send;
double *integral_buffer1, *integral_buffer2;
double *bse_eigenvalues, *Ham_buffer1, *Ham_buffer2;
double time1, time2, time3, time4, time5, time6;
double *GW_eigenvalues, *scf_eigenvalues;
char yy[20] = "./bse_eigenvalues";
FILE *bse_evalues;
DoubleMatrix *V_inv;

  time1 = MPI_Wtime();
  initialise_spk_grid(&ntransitions, &ictxt, &nbsize, job, file);
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, &nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, &nbsize, &mycol, &izero, &npcol);
  itemp = max(1, mpA);
  descinit_(descA, &ntransitions, &ntransitions, &nbsize, &nbsize, &izero, &izero, &ictxt, &itemp, &info);
  if (job->taskid == 0) printf("mpA %7d nqA %7d blocksize %7d ntransitions %7d\n",mpA,nqA,nbsize,ntransitions);

  // ******************************************************************************************
  // * Allocate memory for BSE eigenvectors and eigenvalues                                   *
  // ******************************************************************************************

  buffer_size1 = mpA * nqA;
  if      (job->bse_tda == 0) buffer_size2 = buffer_size1;
  else if (job->bse_tda == 1) buffer_size2 = 1;
  AllocateDoubleArray(&bse_eigenvalues,&ntransitions,job);
  AllocateDoubleArray(&Ham_buffer1,&buffer_size1,job);
  AllocateDoubleArray(&Ham_buffer2,&buffer_size2,job);
  ResetDoubleArray(bse_eigenvalues,&ntransitions);
  ResetDoubleArray(Ham_buffer1,&buffer_size1);
  ResetDoubleArray(Ham_buffer2,&buffer_size2);

  // ******************************************************************************************
  // * Generate three centre integrals needed for matrix elements                             *
  // * Generate RPA, BSE or BSE-TDA Hamiltonian                                               *
  // ******************************************************************************************
  
  AllocateIntArray(&dim_send,&job->numtasks,job);
  count_integral_buffer_sizes(dim_send,fermi,atoms_ax,job,file);
  AllocateDoubleArray(&integral_buffer1,&dim_send[job->taskid],job);
  AllocateDoubleArray(&integral_buffer2,&dim_send[job->taskid],job);
  generate_integral_buffers_molecule_ija(integral_buffer1,integral_buffer2,&ictxt,&nbsize,fermi,atoms,atom_p,shells, \
  gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);

  // ******************************************************************************************
  // * Generate and diagonalise Casida matrix for BSE or BSE-TDA calculation                  *
  // ******************************************************************************************

  if (job->bse_ham == 1) { // RPA calculation needed for GW self-energy in GW/BSE calculation
    tda = job->bse_tda;
    job->bse_ham = 2;
    job->bse_tda = 1;
    bse_hamiltonian(&ictxt,&nbsize,Ham_buffer1,Ham_buffer2,integral_buffer1,integral_buffer2,fermi,atom_p,atoms,\
    atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
    time3 = MPI_Wtime();
    diagonalise_rpa_hamiltonian(Ham_buffer1,&ictxt,&nbsize,fermi,atoms,atoms_ax,job,file);
    job->bse_ham = 1;
    job->bse_tda = tda;
    time4 = MPI_Wtime() - time3;
    if (job->taskid == 0) printf("diagonalise_rpa_hamiltonian                %10.2f\n",time4);
   }

  // ******************************************************************************************
  // * Generate TDHF, TDHF-TDA, BSE or BSE-TDA Hamiltonian                                    *
  // ******************************************************************************************

  bse_hamiltonian(&ictxt,&nbsize,Ham_buffer1,Ham_buffer2,integral_buffer1,integral_buffer2,fermi,atom_p,atoms, \
  atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  DestroyDoubleArray(&integral_buffer1,&dim_send[job->taskid],job);
  DestroyDoubleArray(&integral_buffer2,&dim_send[job->taskid],job);
  DestroyIntArray(&dim_send,&job->numtasks,job);
  time2 = MPI_Wtime() - time1;

  time3 = MPI_Wtime();
  switch (job->bse_tda) {

  case 0: // Full-matrix BSE or TDHF

  diagonalise_bse_hamiltonian(&ictxt,&nbsize,Ham_buffer1,Ham_buffer2,fermi,job,file);

  break;

  case 1: // Tamm-Dancoff Approximation: TDHF or BSE

   if (job->taskid == 0 && job->verbosity > 1) {
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"|                                     GW BSE HAMILTONIAN MATRIX (eV)                                      |\n");
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   for (i = 0; i < mpA; i++) {
     for (j = 0; j < nqA; j++) {
       fprintf(file.out,"%10.4f ",Ham_buffer1[i + mpA * j]);
      }
       fprintf(file.out,"\n");
     }
    }

  diagonalise_bse_tda_hamiltonian(&ictxt,&nbsize,Ham_buffer1,fermi,job,file);

  break;

 } // close switch (job->bse_tda)
  time4 = MPI_Wtime() - time3;

  char zz6[24] = "/bse_eigenvalues";
  read_SCF_GW_eigenvalues(bse_eigenvalues, 0, ntransitions, zz6, job, file);

  if (job->taskid == 0) {
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|                                GW BSE HAMILTONIAN EIGENVALUES (eV)                                      |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  count = 0;
  int numb = ntransitions / 4; int rem = ntransitions - 4 * numb;
  for (i = 0; i < numb; i++) {
  fprintf(file.out,"| %4d           %10.4lf | %4d         %10.4lf | %4d         %10.4lf | %4d         %10.4lf |\n",\
  count+1,bse_eigenvalues[count+0]*au_to_eV,count+2,bse_eigenvalues[count+1]*au_to_eV,\
  count+3,bse_eigenvalues[count+2]*au_to_eV,count+4,bse_eigenvalues[count+3]*au_to_eV);
  count += 4;
 }
  if (rem == 1) 
  fprintf(file.out,"| %4d           %10.4lf |                         |                         |                         |\n",\
  count+1,bse_eigenvalues[count+0]*au_to_eV);
  if (rem == 2) 
  fprintf(file.out,"| %4d           %10.4lf | %4d         %10.4lf |                         |                         |\n",\
  count+1,bse_eigenvalues[count+0]*au_to_eV,count+2,bse_eigenvalues[count+1]*au_to_eV);
  if (rem == 3) 
  fprintf(file.out,"| %4d           %10.4lf | %4d         %10.4lf | %4d         %10.4lf |                         |\n",\
  count+1,bse_eigenvalues[count+0]*au_to_eV,count+2,bse_eigenvalues[count+1]*au_to_eV,count+3,bse_eigenvalues[count+2]*au_to_eV);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"| HAM GEN TIME    %9.2e | HAM DIAG TIME %9.2e |                         |                         |\n",\
  time2,time4);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 }

  DestroyDoubleArray(&bse_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&Ham_buffer1,&buffer_size1,job);
  DestroyDoubleArray(&Ham_buffer2,&buffer_size2,job);
  if (job->taskid == 0) printf("bse                                        %10.2f\n",time2+time4);

}

void rpa_molecule(FERMI* fermi, ATOM* atoms, ATOM_TRAN* atom_p, SHELL* shells, GAUSSIAN* gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates RPA excitations                                                     *
  // ******************************************************************************************
  
  // ******************************************************************************************
  // * Allocate fermi structure                                                               *
  // ******************************************************************************************

int i, j;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int count;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int ictxt, nbsize, cblacs_taskid, itemp, nprow, npcol, myrow, mycol, mpA, nqA, izero = 0, ione = 1;
int descA[9];
int info = 0;
int buffer_size1, buffer_size2;
int *dim_send;
double *integral_buffer1, *integral_buffer2;
double *Ham_buffer1, *Ham_buffer2;
double time1, time2, time3, time4, time5, time6;
double *scf_eigenvalues;
DoubleMatrix *V_inv;

  time1 = MPI_Wtime();
  initialise_spk_grid(&ntransitions, &ictxt, &nbsize, job, file);
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, &nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, &nbsize, &mycol, &izero, &npcol);
  itemp = max(1, mpA);
  descinit_(descA, &ntransitions, &ntransitions, &nbsize, &nbsize, &izero, &izero, &ictxt, &itemp, &info);
  //if (job->taskid == 0) printf("mpA %7d nqA %7d blocksize %7d ntransitions %7d\n",mpA,nqA,nbsize,ntransitions);

  // ******************************************************************************************
  // * Allocate memory for RPA eigenvectors and eigenvalues                                   *
  // ******************************************************************************************

  buffer_size1 = mpA * nqA;
  buffer_size2 = 1;
  AllocateDoubleArray(&Ham_buffer1,&buffer_size1,job);
  AllocateDoubleArray(&Ham_buffer2,&buffer_size2,job);
  ResetDoubleArray(Ham_buffer1,&buffer_size1);

  // ******************************************************************************************
  // * Generate three centre integrals needed for matrix elements                             *
  // * Generate RPA Hamiltonian                                                               *
  // ******************************************************************************************
  
  AllocateIntArray(&dim_send,&job->numtasks,job);
  count_integral_buffer_sizes(dim_send,fermi,atoms_ax,job,file);
  AllocateDoubleArray(&integral_buffer1,&dim_send[job->taskid],job);
  AllocateDoubleArray(&integral_buffer2,&dim_send[job->taskid],job);
  generate_integral_buffers_molecule_ija(integral_buffer1,integral_buffer2,&ictxt,&nbsize,fermi,atoms,atom_p,shells, \
  gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  bse_hamiltonian(&ictxt,&nbsize,Ham_buffer1,Ham_buffer2,integral_buffer1,integral_buffer2,fermi,atom_p,atoms, \
  atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  DestroyDoubleArray(&integral_buffer1,&dim_send[job->taskid],job);
  DestroyDoubleArray(&integral_buffer2,&dim_send[job->taskid],job);
  DestroyIntArray(&dim_send,&job->numtasks,job);
  time2 = MPI_Wtime() - time1;

  time3 = MPI_Wtime();
  diagonalise_rpa_hamiltonian(Ham_buffer1,&ictxt,&nbsize,fermi,atoms,atoms_ax,job,file);
  time4 = MPI_Wtime() - time3;

  DestroyDoubleArray(&Ham_buffer1,&buffer_size1,job);
  DestroyDoubleArray(&Ham_buffer2,&buffer_size2,job);
  if (job->taskid == 0) printf("rpa                                        %10.2f\n",time2+time4);

}

void bse_hamiltonian(int *ictxt, int *nbsize, double *Ham_buffer1, double *Ham_buffer2, double *integral_buffer1, double *integral_buffer2, FERMI* fermi, ATOM_TRAN *atom_p, ATOM *atoms, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k, j1;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int mpA, nqA, row, col, nprow, npcol, myrow, mycol, izero = 0;
int mpa, nqa, mpa_nqa, MPA[job->numtasks],NQA[job->numtasks];;
int cblacs_taskid, nprow_nbsize, nprow_myrow, npcol_nbsize, npcol_mycol;
int begin_j[job->numtasks], end_j[job->numtasks];
int *dim_ham, *offset_j;
int offset;
int num_proc = job->numtasks < dim1 ? job->numtasks : dim1;
double factor;
double time1, time2, time3, time4;
double time5, time6, time7, time8;
double time9, time10;
double *Hamiltonian_buffer;
double *Ham_temp1, *Ham_temp2;
double *GW_eigenvalues, *scf_eigenvalues;

  time1  = MPI_Wtime();
  time4  = k_zero;
  time6  = k_zero;
  time8  = k_zero;
  time10 = k_zero;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  AllocateDoubleArray(&scf_eigenvalues,&nbands,job);
  AllocateDoubleArray(&GW_eigenvalues,&nbands,job);
  ResetDoubleArray(scf_eigenvalues,&nbands);

  char zz6[24] = "/evalfile";
  read_SCF_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz6, job, file);

  AllocateIntArray(&dim_ham,&job->numtasks,job);
  for (row = 0; row < nprow; row++) {
    mpa = numroc_(&ntransitions, nbsize, &row, &izero, &nprow);
    for (col = 0; col < npcol; col++) {
      nqa = numroc_(&ntransitions, nbsize, &col, &izero, &npcol);
      cblacs_taskid = Cblacs_pnum(*ictxt,row,col);
      dim_ham[cblacs_taskid] = mpa * nqa;
      MPA[cblacs_taskid] = mpa;
      NQA[cblacs_taskid] = nqa;
     }
    }
  AllocateIntArray(&offset_j,&dim1,job);
  mpi_begin_end(begin_j, end_j, dim1, num_proc, job, file);
  for (i = num_proc; i < job->numtasks; i++) begin_j[i] = 0;
  for (i = num_proc; i < job->numtasks; i++) end_j[i] = 0;
  for (i = 0; i < dim1; i++) offset_j[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1] * nbands * nbands; }}

  for (j = 0; j < job->numtasks; j++) {
    row = j / npcol; // fix for variable j1 case
    col = j - row * npcol;
    cblacs_taskid = Cblacs_pnum(*ictxt,row,col);
    mpa = MPA[cblacs_taskid];
    nqa = NQA[cblacs_taskid];
    nprow_nbsize = nprow * *nbsize;
    nprow_myrow = ((nprow + row) % nprow) * *nbsize;
    npcol_nbsize = npcol * *nbsize;
    npcol_mycol = ((npcol + col) % npcol) * *nbsize;
    AllocateDoubleArray(&Hamiltonian_buffer,&dim_ham[j],job);

  if (job->bse_ham == 0 || job->bse_ham == 1) {

  switch (job->bse_tda) {

  case 0: // Full-matrix BSE or TDHF

  ResetDoubleArray(Hamiltonian_buffer,&dim_ham[j]);

  if (job->bse_spin == 0) { // singlet excited state: unscreened ladders
  time3 = MPI_Wtime();
  factor = two;
  hamiltonian_in_core_ia_jb(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  &factor,Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time4 += MPI_Wtime() - time3;
  time5 = MPI_Wtime();
  hamiltonian_in_core_ij_ab(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time6 += MPI_Wtime() - time5;
 }
  if (job->bse_spin == 1) { // triplet excited state: unscreened ladders
  time5 = MPI_Wtime();
  hamiltonian_in_core_ij_ab(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time6 += MPI_Wtime() - time5;
 }

  MPI_Reduce(Hamiltonian_buffer,Ham_buffer1,dim_ham[j],MPI_DOUBLE,MPI_SUM,j,MPI_COMM_WORLD);

  ResetDoubleArray(Hamiltonian_buffer,&dim_ham[j]);

  if (job->bse_spin == 0) { // singlet excited state: unscreened ladders and rings
  time9 = MPI_Wtime();
  factor = two;
  hamiltonian_in_core_ia_jb(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  &factor,Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  hamiltonian_in_core_ib_ja(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time10 += MPI_Wtime() - time9;
 }
  if (job->bse_spin == 1) { // triplet excited state: unscreened ladders
  hamiltonian_in_core_ib_ja(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
 }

  MPI_Reduce(Hamiltonian_buffer,Ham_buffer2,dim_ham[j],MPI_DOUBLE,MPI_SUM,j,MPI_COMM_WORLD);
  DestroyDoubleArray(&Hamiltonian_buffer,&dim_ham[j],job);

  break;

  case 1: // Tamm-Dancoff Approximation: TDHF or BSE

  ResetDoubleArray(Hamiltonian_buffer,&dim_ham[j]);

  if (job->bse_spin == 0) { // singlet excited state: unscreened ladders and rings
  time3 = MPI_Wtime();
  factor = two;
  hamiltonian_in_core_ia_jb(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  &factor,Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time4 += MPI_Wtime() - time3;
  time5 = MPI_Wtime();
  hamiltonian_in_core_ij_ab(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time6 += MPI_Wtime() - time5;
 }
  if (job->bse_spin == 1) { // triplet excited state: unscreened ladders
  time3 = MPI_Wtime();
  hamiltonian_in_core_ij_ab(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time4 += MPI_Wtime() - time3;
 }

  MPI_Reduce(Hamiltonian_buffer,Ham_buffer1,dim_ham[j],MPI_DOUBLE,MPI_SUM,j,MPI_COMM_WORLD);
  DestroyDoubleArray(&Hamiltonian_buffer,&dim_ham[j],job);

  break;

 } // close switch (job->bse_tda)

 } // close if (job->bse_ham

  else if (job->bse_ham == 2) {

  ResetDoubleArray(Hamiltonian_buffer,&dim_ham[j]);

  factor = four;
  time3 = MPI_Wtime();
  hamiltonian_in_core_ia_jb(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  &factor,Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time4 += MPI_Wtime() - time3;
  time9 = MPI_Wtime();
  MPI_Reduce(Hamiltonian_buffer,Ham_buffer1,dim_ham[j],MPI_DOUBLE,MPI_SUM,j,MPI_COMM_WORLD);
  DestroyDoubleArray(&Hamiltonian_buffer,&dim_ham[j],job);
  time10 += MPI_Wtime() - time9;
 } // close else if (job->bse_ham

 } // close loop over j

  if (job->bse_cou == 1)
  hamiltonian_in_core_coulomb_exchange_energy(ictxt,nbsize,begin_j,end_j,offset_j,integral_buffer1,integral_buffer2, \
  atoms_ax,fermi,job,file);

  if (job->bse_ham == 1 && job->bse_tda == 0) { // BSE screened ladders
  time7 = MPI_Wtime();
  hamiltonian_in_core_screened_ij_ab_and_ib_ja(ictxt,nbsize,begin_j,end_j,offset_j,GW_eigenvalues,Ham_buffer1,Ham_buffer2,\
  integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time8 = MPI_Wtime() - time7;
 }

  if (job->bse_tda == 0) { // A + B and A - B matrices
  AllocateDoubleArray(&Ham_temp1,&nqA,job);
  AllocateDoubleArray(&Ham_temp2,&nqA,job);
  //fprintf(file.out,"Performing TDA in TDHF\n");  Uncomment next two lines to do TDA by setting B = 0
  //int buffer_size2 = mpA * nqA;
  //ResetDoubleArray(Ham_buffer2,&buffer_size2);
  for (i = 0; i < mpA; i++) {
    for (k = 0; k < nqA; k++) Ham_temp1[k] = Ham_buffer1[i * nqA + k] + Ham_buffer2[i * nqA + k];
    for (k = 0; k < nqA; k++) Ham_temp2[k] = Ham_buffer1[i * nqA + k] - Ham_buffer2[i * nqA + k];
    for (k = 0; k < nqA; k++) Ham_buffer1[i * nqA + k] = Ham_temp1[k]; // A + B matrix
    for (k = 0; k < nqA; k++) Ham_buffer2[i * nqA + k] = Ham_temp2[k]; // A - B matrix
   } 
  DestroyDoubleArray(&Ham_temp1,&nqA,job);
  DestroyDoubleArray(&Ham_temp2,&nqA,job);
 }

  if (job->bse_ham == 1 && job->bse_tda == 1) { // TDA-BSE screened ladders
  time7 = MPI_Wtime();
  hamiltonian_in_core_screened_ij_ab(ictxt,nbsize,begin_j,end_j,offset_j,GW_eigenvalues,Ham_buffer1,\
  integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time8 = MPI_Wtime() - time7;
 }

  if (job->bse_ham == 0 || job->bse_ham == 2) {
  casida_diagonal_molecule_mpi_spk(&ntransitions,ictxt,nbsize,&nbands,fermi,Ham_buffer1,scf_eigenvalues,job,file);
  if (job->bse_tda == 0) {
  casida_diagonal_molecule_mpi_spk(&ntransitions,ictxt,nbsize,&nbands,fermi,Ham_buffer2,scf_eigenvalues,job,file);
 }
 }
  else if (job->bse_ham == 1) {
  casida_diagonal_molecule_mpi_spk(&ntransitions,ictxt,nbsize,&nbands,fermi,Ham_buffer1,GW_eigenvalues,job,file);
  if (job->bse_tda == 0) {
  casida_diagonal_molecule_mpi_spk(&ntransitions,ictxt,nbsize,&nbands,fermi,Ham_buffer2,GW_eigenvalues,job,file);
 }
 }

  DestroyDoubleArray(&scf_eigenvalues,&nbands,job);
  DestroyDoubleArray(&GW_eigenvalues,&nbands,job);
  DestroyIntArray(&dim_ham,&job->numtasks,job);
  DestroyIntArray(&offset_j,&dim1,job);

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("bse_hamiltonian                            %10.2f\n",time2);
  if (job->taskid == 0) printf("hamiltonian_in_core_ia_jb                  %10.2f\n",time4);
  if (job->taskid == 0) printf("hamiltonian_in_core_ij_ab                  %10.2f\n",time6);
  if (job->taskid == 0) printf("hamiltonian_in_core_screened_ij_ab         %10.2f\n",time8);
  if (job->taskid == 0) printf("MPI_Reduce                                 %10.2f\n",time10);

}

void diagonalise_bse_hamiltonian(int *ictxt, int *nbsize, double *Ham_buffer1, double *Ham_buffer2, FERMI* fermi, JOB_PARAM *job, FILES file)

{

int i, j;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int nt = ntransitions;
int buffer_size;
int nprow, npcol, myrow, mycol, mpA, nqA, izero = 0, ione = 1;
int itemp, descA[9];
int info = 0;
int lwork = -1;
double time1, time2, time3, time4;
double alpha = k_one, beta = k_zero;
double *eigenvalues, *bse_eigenvalues;
double *L2TL1, *U, *VT, *X1, *X2;
double *work;
char joba = 'A';
char jobv = 'V';
char jobz = 'V';
char side = 'L';
char diag = 'N';
char uplo0 = 'L';
char uplo1 = 'U';
char uplo = 'U';
char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
char xx[20] = "/bse_eigenvectors";
char yy[20] = "bse_eigenvalues";
FILE *bse_evalues;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
  itemp = max(1, mpA);
  descinit_(descA, &ntransitions, &ntransitions, nbsize, nbsize, &izero, &izero, ictxt, &itemp, &info);

  buffer_size = mpA * nqA;

  AllocateDoubleArray(&eigenvalues,&nt,job);
  AllocateDoubleArray(&bse_eigenvalues,&nt,job);
  AllocateDoubleArray(&U ,&buffer_size,job);
  AllocateDoubleArray(&VT,&buffer_size,job);
  AllocateDoubleArray(&X1,&buffer_size,job);
  AllocateDoubleArray(&X2,&buffer_size,job);
  AllocateDoubleArray(&L2TL1,&buffer_size,job);

  ResetDoubleArray(eigenvalues,&nt);
  ResetDoubleArray(bse_eigenvalues,&nt);
  ResetDoubleArray(U ,&buffer_size);
  ResetDoubleArray(VT,&buffer_size);
  ResetDoubleArray(X1,&buffer_size);
  ResetDoubleArray(X2,&buffer_size);
  ResetDoubleArray(L2TL1,&buffer_size);

  time1 = MPI_Wtime();
  pdpotrf_(&uplo0,&nt,Ham_buffer1,&ione,&ione,descA,&info);
  pdpotrf_(&uplo0,&nt,Ham_buffer2,&ione,&ione,descA,&info);
  block_cyclic_zero_triangle(&uplo,&nt,ictxt,nbsize,Ham_buffer1); 
  block_cyclic_zero_triangle(&uplo,&nt,ictxt,nbsize,Ham_buffer2);
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("pdpotrf spk                                %10.2lf\n",time2);
  for (i = 0; i < buffer_size; i++) L2TL1[i] = Ham_buffer1[i]; 
  pdgemm_(&Trans,&NoTrans,&nt,&nt,&nt,&alpha,Ham_buffer2,&ione,&ione,descA,Ham_buffer1,&ione,&ione,descA,&beta,L2TL1,&ione, \
  &ione,descA);
  //pdtrmm_(&side,&uplo,&Trans,&diag,&nt,&nt,&alpha,Ham_buffer2,&ione,&ione,descA,L2TL1,&ione,&ione,descA);

  time1 = MPI_Wtime();
  work = (double *) calloc(2, sizeof(double)) ;
  pdgesvd_(&jobv,&jobv,&nt,&nt,L2TL1,&ione,&ione,descA,eigenvalues,U,&ione,&ione,descA,VT,&ione,&ione,descA,work,&lwork,&info);
  lwork= (int) work[0];
  free(work);
  work = (double *) calloc(lwork, sizeof(double)) ;
  pdgesvd_(&jobv,&jobv,&nt,&nt,L2TL1,&ione,&ione,descA,eigenvalues,U,&ione,&ione,descA,VT,&ione,&ione,descA,work,&lwork,&info);
  free(work);
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("pdgesvd spk                                %10.2lf\n",time2);

/* uncomment for T amplitudes
 
DoubleMatrix *X, *Y, *T, *I, *XCHECK, *YCHECK;
int *ipiv;
AllocateIntArray(&ipiv,&ntransitions,job);
AllocateDoubleMatrix(&X,&ntransitions,&ntransitions,job);
AllocateDoubleMatrix(&Y,&ntransitions,&ntransitions,job);
AllocateDoubleMatrix(&T,&ntransitions,&ntransitions,job);
AllocateDoubleMatrix(&I,&ntransitions,&ntransitions,job);
AllocateDoubleMatrix(&XCHECK,&ntransitions,&ntransitions,job);
AllocateDoubleMatrix(&YCHECK,&ntransitions,&ntransitions,job);
ResetDoubleMatrix(X);
ResetDoubleMatrix(Y);
ResetDoubleMatrix(T);
ResetDoubleMatrix(I);
ResetDoubleMatrix(XCHECK);
ResetDoubleMatrix(YCHECK);

char equed='N', fact = 'N';
int lda = X->iRows, n = X->iRows, ldaf = X->iRows, ldb = X->iRows, ldx = X->iRows, nrhs = X->iRows, iwork[lda];
double *R, *C, rcond, FERR[lda], BERR[lda];
DoubleMatrix *AF;
AllocateDoubleMatrix(&AF,&ntransitions,&ntransitions,job);

for (i = 0; i < ntransitions; i++) { bse_eigenvalues[i] = eigenvalues[ntransitions - 1 - i]; I->a[i][i] = 1.0; }
int count = ntransitions - 1;
count = 0;
for (i = 0; i < ntransitions; i++) {
for (j = 0; j < ntransitions; j++) {
X->a[i][j] = X1[count] / two / sqrt(eigenvalues[i]);
Y->a[i][j] = X2[count] / two / sqrt(eigenvalues[i]);
XCHECK->a[i][j] = X1[count] / two / sqrt(eigenvalues[i]);
YCHECK->a[i][j] = X2[count] / two / sqrt(eigenvalues[i]);
 ////fprintf(file.out,"%3d X %10.4lf Y %10.4lf\n",i,X->a[i][j],Y->a[i][j]); // print X, Y
 //fprintf(file.out,"%3d X %10.4lf Y %10.4lf %10.4lf\n",\
 i,X1[count]/two/sqrt(eigenvalues[i]),X2[count]/two/sqrt(eigenvalues[i]),eigenvalues[i]); // print X, Y
//count--;
count++;
}fprintf(file.out,"\n");}
dgesv_(&ntransitions,&ntransitions,X->a[0],&X->iRows,ipiv,Y->a[0],&Y->iCols,&info);
//dgesvx_(&fact,&Trans,&n,&nrhs,X->a[0],&lda,AF->a[0],&ldaf,ipiv,&equed,R,C,Y->a[0],&ldb,Y->a[0],&ldx,&rcond,FERR,BERR,work,iwork,&info);

for (i = 0; i < ntransitions; i++) {
for (j = 0; j < ntransitions; j++) {
 T->a[i][j] = Y->a[i][j];
 fprintf(file.out,"%10.4lf",T->a[i][j]); // print X, Y
}fprintf(file.out,"\n");}
//dgemm_(&NoTrans,&NoTrans,&dim1ax,&nqA,&dim1ax,&alph,&V_inv_k1->a[0][0],&dim1ax,temp2b_buffer,&dim1ax,&bet,temp3_buffer,&dim1ax,&ione,&ione);
ResetDoubleMatrix(X);
ResetDoubleMatrix(Y);
for (i = 0; i < ntransitions; i++) {
for (j = 0; j < ntransitions; j++) {
for (int k = 0; k < ntransitions; k++) {
Y->a[i][j] += T->a[i][k] * XCHECK->a[k][j];
X->a[i][j] += (I->a[i][k] + T->a[i][k]) * XCHECK->a[k][j];
}}}

for (i = 0; i < ntransitions; i++) fprintf(file.out,"%3d\n",ipiv[i]);
//for (i = 0; i < ntransitions; i++) {
//for (j = 0; j < ntransitions; j++) {
//fprintf(file.out,"%10.4lf %10.4lf",X->a[i][j], YCHECK->a[i][j]); 
//}fprintf(file.out,"\n");}
for (i = 0; i < ntransitions; i++) {
for (j = 0; j < ntransitions; j++) {
 fprintf(file.out,"%3d X %10.4lf Y %10.4lf X + Y %10.4lf   Y = T . X  %10.4f X + Y  = (1 + T) . X %10.4lf\n",\
 i,XCHECK->a[i][j],YCHECK->a[i][j],XCHECK->a[i][j]+YCHECK->a[i][j],Y->a[i][j],X->a[i][j]); // print X, Y
}fprintf(file.out,"\n");}
*/

  time1 = MPI_Wtime();
  ResetDoubleArray(L2TL1,&buffer_size);
  pdgemm_(&NoTrans,&NoTrans,&nt,&nt,&nt,&alpha,Ham_buffer2,&ione,&ione,descA,U,&ione,&ione,descA,&beta,X1,&ione,&ione,descA);
  for (i = 0; i < buffer_size; i++)  X2[i] = X1[i]; // X1 = L2.U, X2 = L2.U
  pdgemm_(&NoTrans,  &Trans,&nt,&nt,&nt,&alpha,Ham_buffer1,&ione,&ione,descA,VT,&ione,&ione,descA,&beta,L2TL1,&ione,&ione,descA);
  for (i = 0; i < buffer_size; i++)  X1[i] += L2TL1[i]; // X1 = L2.U + L1.VT
  for (i = 0; i < buffer_size; i++)  X2[i] -= L2TL1[i]; // X2 = L2.U - L1.VT
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("pdgemm spk                                 %10.2lf\n",time2);

 // fprintf(file.out,"|                                     X1   HAMILTONIAN MATRIX (eV)                                      |\n");
 // for (i = 0; i < nt; i++) { for (j = 0; j < nt; j++) { fprintf(file.out,"%10.4f ",X1[i + nt * j]); } fprintf(file.out,"\n"); }
 // fprintf(file.out,"|                                     X2   HAMILTONIAN MATRIX (eV)                                      |\n");
 // for (i = 0; i < nt; i++) { for (j = 0; j < nt; j++) { fprintf(file.out,"%10.4f ",X2[i + nt * j]); } fprintf(file.out,"\n"); }
 // int ntnt = nt * nt;
 // double *In;
 // AllocateDoubleArray(&In,&ntnt,job);
 // ResetDoubleArray(In,&ntnt);
 // beta = k_zero;
 // pdgemm_(&Trans,&NoTrans,&nt,&nt,&nt,&alpha,X2,&ione,&ione,descA,X2,&ione,&ione,descA,&beta,In,&ione,&ione,descA);
 // beta = -k_one;
 // pdgemm_(&Trans,&NoTrans,&nt,&nt,&nt,&alpha,X1,&ione,&ione,descA,X1,&ione,&ione,descA,&beta,In,&ione,&ione,descA);
 // fprintf(file.out,"|                                    NORM  HAMILTONIAN MATRIX (eV)                                      |\n");
 // for (i = 0; i < nt; i++) {
 // for (j = 0; j < nt; j++) {
 // fprintf(file.out,"%10.4f ",In[i + nt * j]);
 // }
 // fprintf(file.out,"\n");
 // }
  //print_real_matrix2(In,0,6,k_one,file);
  //DoubleMatrix *tmp, *LAMBDA;
  //AllocateDoubleMatrix(&LAMBDA,&ntransitions,&ntransitions,job);
  //AllocateDoubleMatrix(&tmp,&ntransitions,&ntransitions,job);
  //ResetDoubleMatrix(LAMBDA);
  //ResetDoubleMatrix(tmp);
  //for (i = 0; i < ntransitions; i++) LAMBDA->a[i][i] = k_one / sqrt(bse_eigenvalues[i]);
  //ResetDoubleMatrix(L2TL1);
  //beta = k_zero;
  //DoubleGEMM3(&Trans,&NoTrans,&nt,&nt,&nt,&alpha,&U->a[0],&nt,&LAMBDA->a[0],&nt,&beta,&tmp->a[0],&nt); // U . L+
  //beta = k_zero;
  //DoubleGEMM3(&NoTrans,&Trans,&nt,&nt,&nt,&alpha,&tmp->a[0],&nt,&VT->a[0],&nt,&beta,&L2TL1->a[0],&nt);  // L2T . L1 = U. L+ .VT

  time3 = MPI_Wtime();

  for (i = 0; i < buffer_size; i++) X2[i] += X1[i]; 
  //block_cyclic_to_linear(&nt,ictxt,nbsize,X2,xx,job,file);
  block_cyclic_to_linear_limit(&nt,ictxt,nbsize,job->bse_lim,X2,xx,job,file);

  for (i = 0; i < ntransitions; i++) bse_eigenvalues[i] = eigenvalues[ntransitions - 1 - i];
  //for (i = 0; i < ntransitions; i++) bse_eigenvalues[i] = eigenvalues[i + nt - job->bse_lim];

  char zz8[24] = "/bse_eigenvalues";
  write_SCF_GW_eigenvalues(bse_eigenvalues, 0, ntransitions, zz8, job, file); 

  DestroyDoubleArray(&U,&buffer_size,job);
  DestroyDoubleArray(&VT,&buffer_size,job);
  DestroyDoubleArray(&X1,&buffer_size,job);
  DestroyDoubleArray(&X2,&buffer_size,job);
  DestroyDoubleArray(&L2TL1,&buffer_size,job);
  DestroyDoubleArray(&eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&bse_eigenvalues,&ntransitions,job);

}

void diagonalise_bse_tda_hamiltonian(int *ictxt, int *nbsize, double *Ham_buffer1, FERMI* fermi, JOB_PARAM *job, FILES file)

{

int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int nt = ntransitions;
int buffer_size;
int nprow, npcol, myrow, mycol, mpA, nqA, izero = 0, ione = 1;
int itemp, descA[9];
int info = 0;
int lwork = -1;
double time1, time2, time3, time4;
double *bse_eigenvalues, *bse_eigvec;
double *work;
char jobz = 'V';
char uplo = 'U';
char xx[20] = "/bse_eigenvectors";
char yy[20] = "bse_eigenvalues";
FILE *bse_evalues;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
  itemp = max(1, mpA);
  descinit_(descA, &ntransitions, &ntransitions, nbsize, nbsize, &izero, &izero, ictxt, &itemp, &info);

  buffer_size = mpA * nqA;
  AllocateDoubleArray(&bse_eigenvalues,&nt,job);
  AllocateDoubleArray(&bse_eigvec,&buffer_size,job);
  ResetDoubleArray(bse_eigenvalues,&nt);
  ResetDoubleArray(bse_eigvec,&buffer_size);

  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0)) {

  time1 = MPI_Wtime();
  int lwork;
  double *work;
  work = (double *) calloc(2, sizeof(double)) ;
  if (work==NULL) { 
  fprintf(file.out,"ERROR: allocation of work array on core %3d\n",job->taskid); 
  MPI_Finalize();
  exit(1); 
 }
  lwork=-1;
  pdsyev_(&jobz,&uplo,&ntransitions,Ham_buffer1,&ione,&ione,descA,bse_eigenvalues,bse_eigvec,&ione,&ione,descA,work,&lwork,&info);
  lwork= (int) work[0];
  free(work);
  work = (double *) calloc(lwork, sizeof(double)) ;
  if (work==NULL) { 
  fprintf(file.out,"ERROR: allocation of work array on core %3d\n",job->taskid); 
  MPI_Finalize();
  exit(1); 
 }
  pdsyev_(&jobz,&uplo,&ntransitions,Ham_buffer1,&ione,&ione,descA,bse_eigenvalues,bse_eigvec,&ione,&ione,descA,work,&lwork,&info);
  free(work);
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("pdsyev spk                                 %10.2f\n",time2);
  } // close if (

  block_cyclic_to_linear_limit(&nt,ictxt,nbsize,job->bse_lim,bse_eigvec,xx,job,file);

  if (job->taskid == 0) {
  bse_evalues  = fopen(yy, "wb");
  fwrite(bse_eigenvalues, sizeof(double), ntransitions, bse_evalues);
  fflush(bse_evalues);
  fclose(bse_evalues);
 }

  char zz8[24] = "/bse_eigenvalues";
  write_SCF_GW_eigenvalues(bse_eigenvalues, 0, ntransitions, zz8, job, file); 

  DestroyDoubleArray(&bse_eigenvalues,&nt,job);
  DestroyDoubleArray(&bse_eigvec,&buffer_size,job);

}

void diagonalise_rpa_hamiltonian(double *Hamiltonian, int *ictxt, int *nbsize, FERMI *fermi, ATOM *atoms, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int i, j, k, m, n, s;
int count;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int I1, I2, il, jl, mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int info = 0, ione = 1;
char uplo = 'U';
char jobz = 'V';
double time1, time2, time3, time4;
double *scf_eigenvalues;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  AllocateDoubleArray(&scf_eigenvalues,&nbands,job);
  ResetDoubleArray(scf_eigenvalues,&nbands);

  char zz6[24] = "/evalfile";
  read_SCF_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz6, job, file);

  char xc[22] = "/cas_eigenvectors_mpi";
  char bufcas[120];
  strcpy(bufcas,file.scf_eigvec);
  strcat(bufcas,xc);

  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0)) {
  double root_scf_transition[ntransitions];
  count = 0;
  for (m = 0; m < fermi->occupied[0]; m++) {
    for (n = fermi->occupied[0]; n < nbands; n++) {
      root_scf_transition[count] = sqrt(scf_eigenvalues[n] - scf_eigenvalues[m]);
      count++;
     }
    }
  for (jl = 0; jl < nqA; jl++) {
    I2 = npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + mycol) % npcol) * *nbsize;
    for (il = 0; il < mpA; il++) {
      Hamiltonian[il + mpA * jl] *= root_scf_transition[I2];
     }
    }
  for (jl = 0; jl < nqA; jl++) {
    for (il = 0; il < mpA; il++) {
    I1 = nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize;
    Hamiltonian[il + mpA * jl] *= root_scf_transition[I1];
   }
  }
 }

  if (job->verbosity > 1 && job->taskid == 0) {
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|                                            CASIDA MATRIX (eV)                                           |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  for (il = 0; il < mpA; il++) { for (jl = 0; jl < nqA; jl++) { fprintf(file.out,"%7.3f ",Hamiltonian[il + mpA * jl] * au_to_eV);
  } fprintf(file.out,"\n"); }
 }

  // ******************************************************************************************
  // * Diagonalise Hamiltonian                                                                *
  // ******************************************************************************************

  int buffer_size;
  double *cas_eigenvalues, *eigvec_buffer;
  buffer_size = mpA * nqA;
  if (buffer_size < 1 || buffer_size > 10000 * 10000) buffer_size = 1;
  AllocateDoubleArray(&cas_eigenvalues,&ntransitions,job);
  AllocateDoubleArray(&eigvec_buffer,&buffer_size,job);
  ResetDoubleArray(eigvec_buffer,&buffer_size);

  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0)) {

  int descA[9];
  int itemp = max(1, mpA);
  descinit_(descA, &ntransitions, &ntransitions, nbsize, nbsize, &izero, &izero, ictxt, &itemp, &info);

  time1 = MPI_Wtime();
  int lwork;
  double *work;
  work = (double *) calloc(2, sizeof(double)) ;
  if (work==NULL) { 
  fprintf(file.out,"ERROR: allocation of work array on core %3d\n",job->taskid); 
  MPI_Finalize();
  exit(1); 
 }
  lwork=-1;
  pdsyev_(&jobz,&uplo,&ntransitions,Hamiltonian,&ione,&ione,descA,cas_eigenvalues,eigvec_buffer,&ione,&ione,descA,work,&lwork,&info);
  lwork= (int) work[0];
  free(work);
  work = (double *) calloc(lwork, sizeof(double)) ;
  if (work==NULL) { 
  fprintf(file.out,"ERROR: allocation of work array on core %3d\n",job->taskid); 
  MPI_Finalize();
  exit(1); 
 }
  pdsyev_(&jobz,&uplo,&ntransitions,Hamiltonian,&ione,&ione,descA,cas_eigenvalues,eigvec_buffer,&ione,&ione,descA,work,&lwork,&info);
  free(work);
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("Diagonalise rpa matrix                     %10.2f\n",time2);
  for (i = 0 ; i < ntransitions; i++) cas_eigenvalues[i] = sqrt(cas_eigenvalues[i]); // eigenvalues of Ham are omega^2
 }

  time1 = MPI_Wtime();
  block_cyclic_to_linear_limit(&ntransitions,ictxt,nbsize,job->rpa_lim,eigvec_buffer,xc,job,file);
  DestroyDoubleArray(&eigvec_buffer,&buffer_size,job);
  MPI_Bcast(cas_eigenvalues,ntransitions,MPI_DOUBLE,0,MPI_COMM_WORLD);
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("block_cyclic_to_linear                     %10.2f\n",time2);

  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD,bufcas,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;

  if (job->taskid == 0) {
  double norm, root_scf_transition[ntransitions];
  count = 0;
  for (m = 0; m < fermi->occupied[0]; m++) {
    for (n = fermi->occupied[0]; n < nbands; n++) {
      root_scf_transition[count] = sqrt(scf_eigenvalues[n] - scf_eigenvalues[m]);
      count++;
     }
    }
  DoubleMatrix *cas_eigenvectors;
  if (job->rpa_lim == 0 || job->rpa_lim > ntransitions) job->rpa_lim = ntransitions;

  int block_size = job->rpa_lim, div, rem;

  if (ntransitions > 20000) {

  //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|                                   CASIDA MATRIX EIGENVALUES (eV)                                        |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  count = 0;
  int numb = ntransitions / 4; int rem = ntransitions - 4 * numb;
  for (i = 0; i < numb; i++) {
  fprintf(file.out,"| %4d           %10.4lf | %4d         %10.4lf | %4d         %10.4lf | %4d         %10.4lf |\n",\
  count+1,cas_eigenvalues[count+0]*au_to_eV,count+2,cas_eigenvalues[count+1]*au_to_eV,\
  count+3,cas_eigenvalues[count+2]*au_to_eV,count+4,cas_eigenvalues[count+3]*au_to_eV);
  count += 4;
 }
  if (rem == 1) 
  fprintf(file.out,"| %4d           %10.4lf |                         |                         |                         |\n",\
  count+1,cas_eigenvalues[count+0]*au_to_eV);
  if (rem == 2) 
  fprintf(file.out,"| %4d           %10.4lf | %4d         %10.4lf |                         |                         |\n",\
  count+1,cas_eigenvalues[count+0]*au_to_eV,count+2,cas_eigenvalues[count+1]*au_to_eV);
  if (rem == 3) 
  fprintf(file.out,"| %4d           %10.4lf | %4d         %10.4lf | %4d         %10.4lf |                         |\n",\
  count+1,cas_eigenvalues[count+0]*au_to_eV,count+2,cas_eigenvalues[count+1]*au_to_eV,count+3,cas_eigenvalues[count+2]*au_to_eV);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"| HAM GEN TIME    %9.2e | HAM DIAG TIME %9.2e |                         |                         |\n",\
  time1,time3);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

  time1 = MPI_Wtime();
  block_size = 1000;
  div = job->rpa_lim / block_size;
  rem = job->rpa_lim - (job->rpa_lim / block_size) * block_size;
  long long global_pointer;
  AllocateDoubleMatrix(&cas_eigenvectors,&block_size,&ntransitions,job);
  for (i = 0 ; i < div; i++) {
  ResetDoubleMatrix(cas_eigenvectors);
  global_pointer = (i * block_size) * (ntransitions * sizeof(double));
  //fprintf(file.out,"%d %lld %d\n",i,global_pointer,(i * block_size) * (ntransitions * sizeof(double)));
  MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
  MPI_File_read(fh, &cas_eigenvectors->a[0][0], block_size * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  for (j = 0 ; j < block_size; j++) {
    norm = k_zero;
    for (k = 0 ; k < ntransitions; k++) {
    norm += (root_scf_transition[k] * root_scf_transition[k] / cas_eigenvalues[i * block_size + j] + \
    cas_eigenvalues[i * block_size + j] / root_scf_transition[k] / root_scf_transition[k]) * \
    cas_eigenvectors->a[j][k] * cas_eigenvectors->a[j][k] / two;
   }
    for (k = 0 ; k < ntransitions; k++) {
    cas_eigenvectors->a[j][k] *= root_scf_transition[k] / sqrt(cas_eigenvalues[i * block_size + j] * norm);
   }
  }
  MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
  MPI_File_write(fh, &cas_eigenvectors->a[0][0], block_size * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
 }
  ResetDoubleMatrix(cas_eigenvectors);
  global_pointer = (div * block_size) * (ntransitions * sizeof(double));
  MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
  MPI_File_read(fh, &cas_eigenvectors->a[0][0], rem * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  for (j = 0 ; j < rem; j++) {
    norm = k_zero;
    for (k = 0 ; k < ntransitions; k++) {
    norm += (root_scf_transition[k] * root_scf_transition[k] / cas_eigenvalues[div * block_size + j] + \
    cas_eigenvalues[div * block_size + j] / root_scf_transition[k] / root_scf_transition[k]) * \
    cas_eigenvectors->a[j][k] * cas_eigenvectors->a[j][k] / two;
   }
    for (k = 0 ; k < ntransitions; k++) {
    cas_eigenvectors->a[j][k] *= root_scf_transition[k] / sqrt(cas_eigenvalues[div * block_size + j] * norm);
   }
  }
  MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
  MPI_File_write(fh, &cas_eigenvectors->a[0][0], rem * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  DestroyDoubleMatrix(&cas_eigenvectors,job);
  time2 = MPI_Wtime() - time1;
  printf("transform cas eigenvectors                 %10.2f\n",time2);
 } // close if (job->rpa_lim

  else {

  AllocateDoubleMatrix(&cas_eigenvectors,&job->rpa_lim,&ntransitions,job);
  ResetDoubleMatrix(cas_eigenvectors);
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  MPI_File_read(fh, &cas_eigenvectors->a[0][0], job->rpa_lim * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  for (i = 0 ; i < job->rpa_lim; i++) {
    norm = k_zero;
    for (j = 0 ; j < ntransitions; j++) {
    norm += (root_scf_transition[j] * root_scf_transition[j] / cas_eigenvalues[i] + cas_eigenvalues[i] / root_scf_transition[j] / \
    root_scf_transition[j]) * cas_eigenvectors->a[i][j] * cas_eigenvectors->a[i][j] / two;
   }
    for (j = 0 ; j < ntransitions; j++) {
    cas_eigenvectors->a[i][j] *= root_scf_transition[j] / sqrt(cas_eigenvalues[i] * norm);
    //fprintf(file.out,"CAS %3d %3d %10.4lf %10.4lf  %10.4lf %10.4lf %10.4lf\n",i,j,root_scf_transition[j]*root_scf_transition[j],\
    cas_eigenvalues[i],cas_eigenvectors->a[i][j],norm,sqrt(norm));
   }
  }
 
  if (job->verbosity > 1) {
  //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|                        CASIDA MATRIX EIGENVECTORS AND sqrt[EIGENVALUES] (eV)                            |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  //print_real_eigenvector_matrix2(cas_eigenvectors, cas_eigenvalues, 8, 35, au_to_eV, file);
  fprintf(file.out,"| CAS GEN TIME    %9.2e | CAS DIAG TIME %9.2e |                         |                         |\n", \
  time2,time4);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 }
  else if (job->verbosity == 1) {
  //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|                                   CASIDA MATRIX EIGENVALUES (eV)                                        |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  count = 0;
  int numb = ntransitions / 4; int rem = ntransitions - 4 * numb;
  for (i = 0; i < numb; i++) {
  fprintf(file.out,"| %4d           %10.4lf | %4d         %10.4lf | %4d         %10.4lf | %4d         %10.4lf |\n",\
  count+1,cas_eigenvalues[count+0]*au_to_eV,count+2,cas_eigenvalues[count+1]*au_to_eV,\
  count+3,cas_eigenvalues[count+2]*au_to_eV,count+4,cas_eigenvalues[count+3]*au_to_eV);
  count += 4;
 }
  if (rem == 1) 
  fprintf(file.out,"| %4d           %10.4lf |                         |                         |                         |\n",\
  count+1,cas_eigenvalues[count+0]*au_to_eV);
  if (rem == 2) 
  fprintf(file.out,"| %4d           %10.4lf | %4d         %10.4lf |                         |                         |\n",\
  count+1,cas_eigenvalues[count+0]*au_to_eV,count+2,cas_eigenvalues[count+1]*au_to_eV);
  if (rem == 3) 
  fprintf(file.out,"| %4d           %10.4lf | %4d         %10.4lf | %4d         %10.4lf |                         |\n",\
  count+1,cas_eigenvalues[count+0]*au_to_eV,count+2,cas_eigenvalues[count+1]*au_to_eV,count+3,cas_eigenvalues[count+2]*au_to_eV);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"| HAM GEN TIME    %9.2e | HAM DIAG TIME %9.2e |                         |                         |\n",\
  time1,time3);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 }
  int vector_size = job->rpa_lim * ntransitions;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  MPI_File_write(fh, &cas_eigenvectors->a[0][0], vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
  DestroyDoubleMatrix(&cas_eigenvectors,job);

 } // close else

} //if (job->taskid == 0) {

  MPI_File_close(&fh);

  char buf1[110];
  char xx[14] = "/cas_evalues";
  strcpy(buf1,file.scf_eigvec);
  strcat(buf1,xx);
  MPI_File eh;
  MPI_File_open(MPI_COMM_WORLD,buf1,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&eh) ;
  MPI_File_seek(eh, 0, MPI_SEEK_SET) ;
  for (s = 0; s < job->spin_dim; s++) MPI_File_write(eh, &cas_eigenvalues[s * dim1], ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_close(&eh);


  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&scf_eigenvalues,&nbands,job);

}

void casida_diagonal_molecule_mpi_spk(int *nt, int *ictxt, int *nbsize, int *nbands, FERMI *fermi, double *Ham_buffer, double *eigenvalues, JOB_PARAM *job, FILES file)

{

int m, n, I1, I2, il, jl;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0) != 1) return; 
  mpA = numroc_(nt, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(nt, nbsize, &mycol, &izero, &npcol);

  for (m = 0; m < fermi->occupied[0]; m++) {
    for (n = fermi->occupied[0]; n < *nbands; n++) {
      I1 = m * (*nbands - fermi->occupied[0]) + n - fermi->occupied[0];
      if ((I1 / *nbsize) % nprow == myrow && (I1 / *nbsize) % npcol == mycol) {
      il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
      jl = *nbsize * (I1 / (*nbsize * npcol)) + I1 % *nbsize;
      Ham_buffer[il + mpA * jl] += eigenvalues[n] - eigenvalues[m];
     }
    }
   }

}

void hamiltonian_in_core_ia_jb(int *begin_j, int *end_j, int *offset_j, int *nbsize, int *mpA, int *nqA, int *nprow_nbsize, int *npcol_nbsize, int* nprow_myrow, int *npcol_mycol, double *factor, double *Hamiltonian_buffer, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int a1, il, jl, j1, I1, I2, m, n, r, s;
int nd6, offset1;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];

  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = offset_j[j1];
    for (jl = 0; jl < *nqA; jl++) {
      I2 = *npcol_nbsize * (jl / *nbsize) + jl % *nbsize + *npcol_mycol;
      r = I2 / nvir;
      s = I2  - r  * nvir + nocc;
      for (il = 0; il < *mpA; il++) {
        I1 = *nprow_nbsize * (il / *nbsize) + il % *nbsize + *nprow_myrow;
        m = I1 / nvir;
        n = I1  - m  * nvir + nocc;
        for (a1 = 0; a1 < nd6; a1++) {
          Hamiltonian_buffer[il + *mpA * jl] += *factor * integral_buffer1[offset1 + m * nbands * nd6 + n * nd6 + a1] * \
                                                          integral_buffer2[offset1 + r * nbands * nd6 + s * nd6 + a1];
         }
        }
       }
      }

}

void hamiltonian_in_core_ij_ab(int *begin_j, int *end_j, int *offset_j, int *nbsize, int *mpA, int *nqA, int *nprow_nbsize, int *npcol_nbsize, int* nprow_myrow, int *npcol_mycol, double *Hamiltonian_buffer, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int a1, il, jl, j1, I1, I2, m, n, r, s;
int nd6, offset1;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];

  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = offset_j[j1];
    for (jl = 0; jl < *nqA; jl++) {
      I2 = *npcol_nbsize * (jl / *nbsize) + jl % *nbsize + *npcol_mycol;
      r = I2 / nvir;
      s = I2  - r  * nvir + nocc;
      for (il = 0; il < *mpA; il++) {
        I1 = *nprow_nbsize * (il / *nbsize) + il % *nbsize + *nprow_myrow;
        m = I1 / nvir;
        n = I1  - m  * nvir + nocc;
        for (a1 = 0; a1 < nd6; a1++) {
          Hamiltonian_buffer[il + *mpA * jl] -= integral_buffer1[offset1 + m * nbands * nd6 + r * nd6 + a1] * \
                                                integral_buffer2[offset1 + n * nbands * nd6 + s * nd6 + a1];
         }
        }
       }
      }

}

void hamiltonian_in_core_ib_ja(int *begin_j, int *end_j, int *offset_j, int *nbsize, int *mpA, int *nqA, int *nprow_nbsize, int *npcol_nbsize, int* nprow_myrow, int *npcol_mycol, double *Hamiltonian_buffer, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int a1, il, jl, j1, I1, I2, m, n, r, s;
int nd6, offset1;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];

  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = offset_j[j1];
    for (jl = 0; jl < *nqA; jl++) {
      I2 = *npcol_nbsize * (jl / *nbsize) + jl % *nbsize + *npcol_mycol;
      r = I2 / nvir;
      s = I2  - r  * nvir + nocc;
      for (il = 0; il < *mpA; il++) {
        I1 = *nprow_nbsize * (il / *nbsize) + il % *nbsize + *nprow_myrow;
        m = I1 / nvir;
        n = I1  - m  * nvir + nocc;
        for (a1 = 0; a1 < nd6; a1++) {
          Hamiltonian_buffer[il + *mpA * jl] -= integral_buffer1[offset1 + m * nbands * nd6 + s * nd6 + a1] * \
                                                integral_buffer2[offset1 + n * nbands * nd6 + r * nd6 + a1];
         }
        }
       }
      }

}

void hamiltonian_in_core_screened_ij_ab(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *GW_eigenvalues, double *Ham_buffer1, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, k, l, a1, il, jl, j1, I1, I2, m, n, r, s;
int nd6, offset, offset1, offset2;
int *offset_j1;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int dim4, mdim;
int mdim1, mdim2, l1, l2;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int itr, rem;
double time1, time2, time3, time4, time5, time6, time7, time8;
double *temp4;
double *M, *M_buffer;
double *cas_eigenvalues, *scf_eigenvalues;
double Sigma[nbands], dSigma_dE[nbands], sigma_factor, denom;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->rpa_lim == 0 || job->rpa_lim > ntransitions) job->rpa_lim = ntransitions;
  if (job->taskid == 0) printf("RPA vector limit %3d\n",job->rpa_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->rpa_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->rpa_lim / nvir;
  rem = job->rpa_lim - itr * nvir;
  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4(temp4,integral_buffer2,fermi,atoms_ax,job,file);
  time4 = MPI_Wtime() - time3;

  AllocateIntArray(&offset_j1,&dim1,job);
  for (i = 0; i < dim1; i++) offset_j1[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j1[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1]; }}

  // ******************************************************************************************
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  AllocateDoubleArray(&scf_eigenvalues,&nbands,job);
  ResetDoubleArray(scf_eigenvalues,&nbands);
  AllocateDoubleArray(&cas_eigenvalues,&job->rpa_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->rpa_lim);
  time3 = MPI_Wtime();

  char zz6[24] = "/evalfile";
  read_SCF_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz6, job, file);

  char zz7[24] = "/cas_evalues";
  read_SCF_GW_eigenvalues(cas_eigenvalues, 0, job->rpa_lim, zz7, job, file);
  time4 = MPI_Wtime() - time3;

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with <alpha|occ-occ>                             *
  // ******************************************************************************************
  
  time2 = k_zero;
  time4 = k_zero;
  time6 = k_zero;
  time8 = k_zero;

  for (i = 0; i < nbands; i++) Sigma[i] = k_zero;
  for (i = 0; i < nbands; i++) dSigma_dE[i] = k_zero;

  for (l1 = 0; l1 < itr; l1++) {

  time1 = MPI_Wtime();
  mdim = nbands * nbands * nvir;
  AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->rpa_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l2 = 0; l2 < nvir; l2++) {
      for (m = 0; m < nbands; m++) {
        for (r = 0; r < nbands; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer[m * nbands * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
           }
          }
         }
        }
       }
  time2 += MPI_Wtime() - time1;

  AllocateDoubleArray(&M,&mdim,job);
  ResetDoubleArray(M,&mdim);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  DestroyDoubleArray(&M_buffer,&mdim,job);

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

  for (i = 0; i < nbands; i++) {
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }
  time6 += MPI_Wtime() - time5;

  time7 = MPI_Wtime();
  for (il = 0; il < mpA; il++) {
    I1 = nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize;
    m = I1 / nvir;
    n = I1  - m  * nvir + nocc;
    for (jl = 0; jl < nqA; jl++) {
      I2 = npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + mycol) % npcol) * *nbsize;
      r = I2 / nvir;
      s = I2  - r  * nvir + nocc;
      for (l2 = 0; l2 < nvir; l2++) {
        Ham_buffer1[il + mpA * jl] += two * M[m * nbands * nvir + r * nvir + l2] * M[n * nbands * nvir + s * nvir + l2] / \
        cas_eigenvalues[l1 * nvir + l2];
       }
      }
     }
  time8 += MPI_Wtime() - time7;
  DestroyDoubleArray(&M,&mdim,job);

 } // close loop on l1

  if (rem > 0) {

  time1 = MPI_Wtime();
  mdim = nbands * nbands * nvir;
  AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->rpa_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l2 = 0; l2 < rem; l2++) {
      for (m = 0; m < nbands; m++) {
        for (r = 0; r < nbands; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer[m * nbands * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
           }
          }
         }
        }
       }
  time2 += MPI_Wtime() - time1;

  AllocateDoubleArray(&M,&mdim,job);
  ResetDoubleArray(M,&mdim);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  DestroyDoubleArray(&M_buffer,&mdim,job);

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < rem; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

  for (i = 0; i < nbands; i++) {
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < rem; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }
  time6 += MPI_Wtime() - time5;

  time7 = MPI_Wtime();
  for (il = 0; il < mpA; il++) {
    I1 = nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize;
    m = I1 / nvir;
    n = I1  - m  * nvir + nocc;
    for (jl = 0; jl < nqA; jl++) {
      I2 = npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + mycol) % npcol) * *nbsize;
      r = I2 / nvir;
      s = I2  - r  * nvir + nocc;
      for (l2 = 0; l2 < rem; l2++) {
        Ham_buffer1[il + mpA * jl] += two * M[m * nbands * nvir + r * nvir + l2] * M[n * nbands * nvir + s * nvir + l2] / \
        cas_eigenvalues[l1 * nvir + l2];
       }
      }
     }
  time8 += MPI_Wtime() - time7;
  DestroyDoubleArray(&M,&mdim,job);

 } // if (rem > 0)

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);
  if (job->taskid == 0 ) printf("Contract                                   %10.2f\n",time8);

  if (job->taskid == 0) {

  //fprintf(file.out,"GW Eigenvalues: %4d RPA states used\n",job->rpa_lim);
  fprintf(file.out,"\n\n===========================================================================================================\n");
  fprintf(file.out,"|                                GW EIGENVALUES AND SELF-ENERGY (eV)                                      |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|    i    |      E[i]     |   E[i] + S[i] | E[i] + Z S[i] |      S[i]     |    Z S[i]     |       Z       |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  for (m = 0; m < nbands; m++) {
    fprintf(file.out,"|   %3d   | %10.4lf    | %10.4lf    | %10.4lf    | %10.4lf    | %10.4lf    | %10.4lf    |\n", m + 1, \
    scf_eigenvalues[m] * au_to_eV, (scf_eigenvalues[m] + Sigma[m]) * au_to_eV, (scf_eigenvalues[m] + Sigma[m] / \
    (k_one - dSigma_dE[m])) * au_to_eV, Sigma[m] * au_to_eV, Sigma[m] / (k_one - dSigma_dE[m]) * au_to_eV, \
     k_one / (k_one - dSigma_dE[m]));
   }
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fflush(file.out);

 }

  DestroyIntArray(&offset_j1,&dim1,job);
  DestroyDoubleArray(&scf_eigenvalues,&nbands,job);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&temp4,&dim4,job);

}

void hamiltonian_in_core_screened_ij_ab_and_ib_ja(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *GW_eigenvalues, double *Ham_buffer1, double *Ham_buffer2, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, k, l, a1, il, jl, j1, I1, I2, m, n, r, s;
int nd6, offset, offset1, offset2;
int *offset_j1;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int dim4, mdim;
int mdim1, mdim2, l1, l2;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int itr, rem;
double time1, time2, time3, time4, time5, time6, time7, time8;
double *temp4;
double *M, *M_buffer;
double *cas_eigenvalues, *scf_eigenvalues;
double Sigma[nbands], dSigma_dE[nbands], sigma_factor, denom;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->rpa_lim == 0 || job->rpa_lim > ntransitions) job->rpa_lim = ntransitions;
  if (job->taskid == 0) printf("RPA vector limit %3d\n",job->rpa_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->rpa_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->rpa_lim / nvir;
  rem = job->rpa_lim - itr * nvir;

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4(temp4,integral_buffer2,fermi,atoms_ax,job,file);
  time4 = MPI_Wtime() - time3;

  AllocateIntArray(&offset_j1,&dim1,job);
  for (i = 0; i < dim1; i++) offset_j1[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j1[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1]; }}

  // ******************************************************************************************
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  AllocateDoubleArray(&scf_eigenvalues,&nbands,job);
  ResetDoubleArray(scf_eigenvalues,&nbands);
  AllocateDoubleArray(&cas_eigenvalues,&job->rpa_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->rpa_lim);
  time3 = MPI_Wtime();

  char zz6[24] = "/evalfile";
  read_SCF_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz6, job, file);

  char zz7[24] = "/cas_evalues";
  read_SCF_GW_eigenvalues(cas_eigenvalues, 0, job->rpa_lim, zz7, job, file);
  time4 = MPI_Wtime() - time3;

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with <alpha|occ-occ>                             *
  // ******************************************************************************************
  
  time2 = k_zero;
  time4 = k_zero;
  time6 = k_zero;
  time8 = k_zero;

  for (i = 0; i < nbands; i++) Sigma[i] = k_zero;
  for (i = 0; i < nbands; i++) dSigma_dE[i] = k_zero;
  mdim = nbands * nbands * nvir;
  AllocateDoubleArray(&M,&mdim,job);
  AllocateDoubleArray(&M_buffer,&mdim,job);

  for (l1 = 0; l1 < itr; l1++) {

  time1 = MPI_Wtime();
  ResetDoubleArray(M_buffer,&mdim);
  ResetDoubleArray(M,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->rpa_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l2 = 0; l2 < nvir; l2++) {
      for (m = 0; m < nbands; m++) {
        for (r = 0; r < nbands; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer[m * nbands * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
           }
          }
         }
        }
       }
  time2 += MPI_Wtime() - time1;

  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

  for (i = 0; i < nbands; i++) {
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }
  time6 += MPI_Wtime() - time5;

  time7 = MPI_Wtime();
  for (il = 0; il < mpA; il++) {
    I1 = nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize;
    m = I1 / nvir;
    n = I1  - m  * nvir + nocc;
    for (jl = 0; jl < nqA; jl++) {
      I2 = npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + mycol) % npcol) * *nbsize;
      r = I2 / nvir;
      s = I2  - r  * nvir + nocc;
      for (l2 = 0; l2 < nvir; l2++) {
        Ham_buffer1[il + mpA * jl] += two * M[m * nbands * nvir + r * nvir + l2] * M[n * nbands * nvir + s * nvir + l2] / \
        cas_eigenvalues[l1 * nvir + l2];
       }
      }
     }

  for (il = 0; il < mpA; il++) {
    I1 = nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize;
    m = I1 / nvir;
    n = I1  - m  * nvir + nocc;
    for (jl = 0; jl < nqA; jl++) {
      I2 = npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + mycol) % npcol) * *nbsize;
      r = I2 / nvir;
      s = I2  - r  * nvir + nocc;
      for (l2 = 0; l2 < nvir; l2++) {
        Ham_buffer2[il + mpA * jl] += two * M[m * nbands * nvir + s * nvir + l2] * M[n * nbands * nvir + r * nvir + l2] / \
        cas_eigenvalues[l1 * nvir + l2];
       }
      }
     }
  time8 += MPI_Wtime() - time7;

 } // close loop on l1

  if (rem > 0) {

  time1 = MPI_Wtime();
  mdim = nbands * nbands * nvir;
  AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->rpa_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l2 = 0; l2 < rem; l2++) {
      for (m = 0; m < nbands; m++) {
        for (r = 0; r < nbands; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer[m * nbands * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
           }
          }
         }
        }
       }
  time2 += MPI_Wtime() - time1;

  AllocateDoubleArray(&M,&mdim,job);
  ResetDoubleArray(M,&mdim);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  DestroyDoubleArray(&M_buffer,&mdim,job);

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < rem; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

  for (i = 0; i < nbands; i++) {
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < rem; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }
  time6 += MPI_Wtime() - time5;

  time7 = MPI_Wtime();
  for (il = 0; il < mpA; il++) {
    I1 = nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize;
    m = I1 / nvir;
    n = I1  - m  * nvir + nocc;
    for (jl = 0; jl < nqA; jl++) {
      I2 = npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + mycol) % npcol) * *nbsize;
      r = I2 / nvir;
      s = I2  - r  * nvir + nocc;
      for (l2 = 0; l2 < rem; l2++) {
        Ham_buffer1[il + mpA * jl] += two * M[m * nbands * nvir + r * nvir + l2] * M[n * nbands * nvir + s * nvir + l2] / \
        cas_eigenvalues[l1 * nvir + l2];
       }
      }
     }
  time8 += MPI_Wtime() - time7;
  DestroyDoubleArray(&M,&mdim,job);

 } // if (rem > 0)

  DestroyDoubleArray(&M_buffer,&mdim,job);
  DestroyDoubleArray(&M,&mdim,job);

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);
  if (job->taskid == 0 ) printf("Contract                                   %10.2f\n",time8);

  if (job->taskid == 0) {

  fprintf(file.out,"\n\n===========================================================================================================\n");
  fprintf(file.out,"|                                GW EIGENVALUES AND SELF-ENERGY (eV)                                      |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|    i    |      E[i]     |   E[i] + S[i] | E[i] + Z S[i] |      S[i]     |    Z S[i]     |       Z       |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  for (m = 0; m < nbands; m++) {
    fprintf(file.out,"|   %3d   | %10.4lf    | %10.4lf    | %10.4lf    | %10.4lf    | %10.4lf    | %10.4lf    |\n", m + 1, \
    scf_eigenvalues[m] * au_to_eV, (scf_eigenvalues[m] + Sigma[m]) * au_to_eV, (scf_eigenvalues[m] + Sigma[m] / \
    (k_one - dSigma_dE[m])) * au_to_eV, Sigma[m] * au_to_eV, Sigma[m] / (k_one - dSigma_dE[m]) * au_to_eV, \
     k_one / (k_one - dSigma_dE[m]));
   }
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fflush(file.out);

 }

  DestroyDoubleArray(&scf_eigenvalues,&nbands,job);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&temp4,&dim4,job);

}

void hamiltonian_in_core_coulomb_exchange_energy(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int a1, il, jl, j1, I1, I2, m, n, r, s;
int nd6, offset1;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
double coulomb_buffer,  coulomb_energy  = k_zero;
double exchange_buffer, exchange_energy = k_zero;

  coulomb_buffer  = k_zero;
  exchange_buffer = k_zero;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = offset_j[j1];
    for (m = 0; m < nocc; m++) {
      for (r = 0; r < nocc; r++) {
        for (a1 = 0; a1 < nd6; a1++) {
          coulomb_buffer += two * integral_buffer1[offset1 + m * nbands * nd6 + m * nd6 + a1] * \
                                  integral_buffer2[offset1 + r * nbands * nd6 + r * nd6 + a1];
          exchange_buffer      -= integral_buffer1[offset1 + m * nbands * nd6 + r * nd6 + a1] * \
                                  integral_buffer2[offset1 + m * nbands * nd6 + r * nd6 + a1];
         }
        }
       }
      }

  MPI_Reduce(&coulomb_buffer,&coulomb_energy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&exchange_buffer,&exchange_energy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if (job->taskid == 0) printf("coulomb energy                       %16.8f\n",coulomb_energy);
  if (job->taskid == 0) printf("exchange energy                      %16.8f\n",exchange_energy);

}

void self_energy_diagonal(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, k, l, a1, il, jl, j1, I1, I2, m, n, r, s;
int nd6, offset, offset1, offset2;
int *offset_j1;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int dim4, mdim;
int mdim1, mdim2, l1, l2;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int itr, rem;
int ii, begin_j2, end_j2;
double time1, time2, time3, time4, time5, time6, time7, time8;
double *temp4;
double *M, *M_buffer;
double *cas_eigenvalues, *scf_eigenvalues, *GW_eigenvalues;
double Sigma[nbands], dSigma_dE[nbands], sigma_factor, denom;
double EE, start_energy, energy_increment;
FILE *file_prt[job->nspectra];
struct filename { char name[20]; };
struct filename *filearray;
DoubleMatrix *Sigma_plot_buffer, *Sigma_plot, *dSigma_dE_plot_buffer, *dSigma_dE_plot;

  start_energy = job->energy_range[0];
  energy_increment = (job->energy_range[1] - job->energy_range[0]) / (double) job->npoints;

  // ******************************************************************************************
  // * Open files and allocate memory for self-energy diagonal matrix element plots           *
  // ******************************************************************************************

  if (job->self_plot == 1) {

  AllocateDoubleMatrix(&Sigma_plot,&job->nspectra,&job->npoints,job);
  AllocateDoubleMatrix(&Sigma_plot_buffer,&job->nspectra,&job->npoints,job);
  AllocateDoubleMatrix(&dSigma_dE_plot,&job->nspectra,&job->npoints,job);
  AllocateDoubleMatrix(&dSigma_dE_plot_buffer,&job->nspectra,&job->npoints,job);
  ResetDoubleMatrix(Sigma_plot);
  ResetDoubleMatrix(Sigma_plot_buffer);
  ResetDoubleMatrix(dSigma_dE_plot);
  ResetDoubleMatrix(dSigma_dE_plot_buffer);

  filearray = (struct filename *) malloc(10 * sizeof(struct filename));
  if (filearray == NULL) {
    fprintf(file.out, "CANNOT OPEN MEMORY FOR FILEARRAY IN GW_self_energy\n");
    exit(1);
  }

  strcpy(filearray[0].name, "self_energy_01.dat");
  strcpy(filearray[1].name, "self_energy_02.dat");
  strcpy(filearray[2].name, "self_energy_03.dat");
  strcpy(filearray[3].name, "self_energy_04.dat");
  strcpy(filearray[4].name, "self_energy_05.dat");
  strcpy(filearray[5].name, "self_energy_06.dat");
  strcpy(filearray[6].name, "self_energy_07.dat");
  strcpy(filearray[7].name, "self_energy_08.dat");
  strcpy(filearray[8].name, "self_energy_09.dat");
  strcpy(filearray[9].name, "self_energy_10.dat");

  for (k = 0; k < job->nspectra; k++) {
    file_prt[k] = fopen(filearray[k].name, "w");
    if (file_prt[k] == NULL) {
      fprintf(file.out, "CANNOT OPEN FILES FOR PLOTTING IN GW_self_energy\n");
      exit(1);
    }
   }

  } // close if (job->self_plot

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->rpa_lim == 0 || job->rpa_lim > ntransitions) job->rpa_lim = ntransitions;
  if (job->taskid == 0) printf("RPA vector limit %3d\n",job->rpa_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->rpa_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->rpa_lim / nvir;
  rem = job->rpa_lim - itr * nvir;

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4(temp4,integral_buffer2,fermi,atoms_ax,job,file);
  time4 = MPI_Wtime() - time3;

  AllocateIntArray(&offset_j1,&dim1,job);
  for (i = 0; i < dim1; i++) offset_j1[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j1[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1]; }}

  // ******************************************************************************************
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  AllocateDoubleArray(&GW_eigenvalues,&nbands,job);
  AllocateDoubleArray(&scf_eigenvalues,&nbands,job);
  ResetDoubleArray(scf_eigenvalues,&nbands);
  AllocateDoubleArray(&cas_eigenvalues,&job->rpa_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->rpa_lim);
  time3 = MPI_Wtime();
  char zz6[24] = "/evalfile";
  read_SCF_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz6, job, file);
  char zz7[24] = "/cas_evalues";
  read_SCF_GW_eigenvalues(cas_eigenvalues, 0, job->rpa_lim, zz7, job, file);
  time4 = MPI_Wtime() - time3;

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with <alpha|occ-occ>                             *
  // ******************************************************************************************
  
  time2 = k_zero;
  time4 = k_zero;
  time6 = k_zero;
  time8 = k_zero;

  for (i = 0; i < nbands; i++) Sigma[i] = k_zero;
  for (i = 0; i < nbands; i++) dSigma_dE[i] = k_zero;

  mdim = nbands * nbands * nvir;
  AllocateDoubleArray(&M,&mdim,job);
  AllocateDoubleArray(&M_buffer,&mdim,job);

  for (l1 = 0; l1 < itr; l1++) {

  time1 = MPI_Wtime();
  ResetDoubleArray(M_buffer,&mdim);
  ResetDoubleArray(M,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->rpa_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l2 = 0; l2 < nvir; l2++) {
      for (m = 0; m < nbands; m++) {
        for (r = 0; r < nbands; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer[m * nbands * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
           }
          }
         }
        }
       }
  time2 += MPI_Wtime() - time1;

  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;

  if (job->self_plot == 1) {

  for (l = 0; l < job->nspectra; l++) {
    i = fermi->plot_bands[l];
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        for (ii = 0; ii < job->npoints; ii++) {
          EE = start_energy + ii * energy_increment;
          denom = EE - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
          sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
          Sigma_plot_buffer->a[l][ii] += sigma_factor * denom * denom / (denom * denom + 1.0e-04);
          dSigma_dE_plot_buffer->a[l][ii] -= sigma_factor / denom;
         }
        }
       }
      }

  for (l = 0; l < job->nspectra; l++) {
    i = fermi->plot_bands[l];
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        for (ii = 0; ii < job->npoints; ii++) {
          EE = start_energy + ii * energy_increment;
          denom = EE - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
          sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
          Sigma_plot_buffer->a[l][ii] += sigma_factor * denom * denom / (denom * denom + 1.0e-04);
          dSigma_dE_plot_buffer->a[l][ii] -= sigma_factor / denom;
         }
        }
       }
      }

  } // close if job->self_plot

 //fprintf(file.out,"%10.4f %10.4f %3d\n",(start_energy+ 508 * energy_increment)*au_to_eV,scf_eigenvalues[5]*au_to_eV,fermi->plot_bands[5]);
 //fprintf(file.out,"plot %3d %10.4f %10.4f\n", l1, (start_energy + 508 * energy_increment) * au_to_eV, Sigma_plot_buffer->a[5][508] * au_to_eV);
 //fprintf(file.out,"%10.4f %10.4f %3d\n",(start_energy+ 540 * energy_increment)*au_to_eV,scf_eigenvalues[17]*au_to_eV,fermi->plot_bands[0]);
 //fprintf(file.out,"plot %3d %10.4f %10.4f\n", l1, (start_energy + 540 * energy_increment) * au_to_eV, Sigma_plot_buffer->a[0][540] * au_to_eV);

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

  for (i = 0; i < nbands; i++) {
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

 //fprintf(file.out,"not  %3d %10.4f %10.4f\n", l1, scf_eigenvalues[5] * au_to_eV, Sigma[5] * au_to_eV);
 //fprintf(file.out,"not  %3d %10.4f %10.4f\n", l1, scf_eigenvalues[17] * au_to_eV, Sigma[17] * au_to_eV);

  time6 += MPI_Wtime() - time5;

 } // close loop on l1

  if (rem > 0) {

  time1 = MPI_Wtime();
  ResetDoubleArray(M_buffer,&mdim);
  ResetDoubleArray(M,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->rpa_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l2 = 0; l2 < rem; l2++) {
      for (m = 0; m < nbands; m++) {
        for (r = 0; r < nbands; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer[m * nbands * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
           }
          }
         }
        }
       }
  time2 += MPI_Wtime() - time1;

  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;

  if (job->self_plot == 1) {

  for (l = 0; l < job->nspectra; l++) {
    i = fermi->plot_bands[l];
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < rem; l2++) {
        for (ii = 0; ii < job->npoints; ii++) {
          EE = start_energy + ii * energy_increment;
          denom = EE - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
          sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
          Sigma_plot_buffer->a[l][ii] += sigma_factor * denom * denom / (denom * denom + 1.0e-04);
          dSigma_dE_plot_buffer->a[l][ii] -= sigma_factor / denom;
         }
        }
       }
      }

  for (l = 0; l < job->nspectra; l++) {
    i = fermi->plot_bands[l];
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < rem; l2++) {
        for (ii = 0; ii < job->npoints; ii++) {
          EE = start_energy + ii * energy_increment;
          denom = EE - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
          sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
          Sigma_plot_buffer->a[l][ii] += sigma_factor * denom * denom / (denom * denom + 1.0e-04);
          dSigma_dE_plot_buffer->a[l][ii] -= sigma_factor / denom;
         }
        }
       }
      }

  } // close if job->self_plot

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < rem; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

  for (i = 0; i < nbands; i++) {
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < rem; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

  time6 += MPI_Wtime() - time5;

 } // if rem > 0

  if (job->self_plot == 1) {
  MPI_Reduce(&Sigma_plot_buffer->a[0][0],&Sigma_plot->a[0][0],job->nspectra*job->npoints,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&dSigma_dE_plot_buffer->a[0][0],&dSigma_dE_plot->a[0][0],job->nspectra*job->npoints,MPI_DOUBLE,MPI_SUM,0,\
  MPI_COMM_WORLD);
  if (job->taskid == 0)
  for (i = 0; i < job->nspectra; i++) for (ii = 0; ii < job->npoints; ii++) fprintf(file_prt[i],"%10.4lf %10.4lf %10.4lf\n", \
  (start_energy + ii * energy_increment) * au_to_eV, Sigma_plot->a[i][ii] * au_to_eV, \
  Sigma_plot->a[i][ii] / (k_one - dSigma_dE_plot->a[i][ii]) * au_to_eV);
  DestroyDoubleMatrix(&dSigma_dE_plot,job);
  DestroyDoubleMatrix(&dSigma_dE_plot_buffer,job);
  DestroyDoubleMatrix(&Sigma_plot,job);
  DestroyDoubleMatrix(&Sigma_plot_buffer,job);
 } // close if (job->self_plot
  DestroyDoubleArray(&M_buffer,&mdim,job);
  DestroyDoubleArray(&M,&mdim,job);

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);

  if (job->taskid == 0) {

  fprintf(file.out,"\n\n===========================================================================================================\n");
  fprintf(file.out,"|                                GW EIGENVALUES AND SELF-ENERGY (eV)                                      |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|    i    |      E[i]     |   E[i] + S[i] | E[i] + Z S[i] |      S[i]     |    Z S[i]     |       Z       |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  for (m = 0; m < nbands; m++) {
    fprintf(file.out,"|   %3d   | %10.4lf    | %10.4lf    | %10.4lf    | %10.4lf    | %10.4lf    | %10.4lf    |\n", m + 1, \
    scf_eigenvalues[m] * au_to_eV, (scf_eigenvalues[m] + Sigma[m]) * au_to_eV, (scf_eigenvalues[m] + Sigma[m] / \
    (k_one - dSigma_dE[m])) * au_to_eV, Sigma[m] * au_to_eV, Sigma[m] / (k_one - dSigma_dE[m]) * au_to_eV, \
     k_one / (k_one - dSigma_dE[m]));
   }
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fflush(file.out);

 }

  char zz8[14] = "/GW_evalues";
  write_SCF_GW_eigenvalues(GW_eigenvalues, 0, job->spin_dim * nbands, zz8, job, file); 

  DestroyDoubleArray(&GW_eigenvalues,&nbands,job);
  DestroyDoubleArray(&scf_eigenvalues,&nbands,job);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&temp4,&dim4,job);

}

void count_integral_buffer_sizes(int *dim_send, FERMI *fermi, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int i, j1;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int num_proc = job->numtasks < dim1 ? job->numtasks : dim1;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int begin_j[job->numtasks], end_j[job->numtasks];

  mpi_begin_end(begin_j, end_j, dim1, num_proc, job, file);
  for (i = num_proc; i < job->numtasks; i++) begin_j[i] = 0;
  for (i = num_proc; i < job->numtasks; i++) end_j[i] = 0;
  for (i = 0; i < job->numtasks; i++) { dim_send[i] = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) {
  dim_send[i] += atoms_ax->bfnnumb_sh[j1] * nbands * nbands; }}

}

void read_write_scf_eigenvectors(FERMI *fermi, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Read in SCF eigenvectors from disk and write to MPI file scf_evec_spk                  *
  // ******************************************************************************************

int i, j;
int vector_size;
int dimk, dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
DoubleMatrix *eigvec;
ComplexMatrix *scf_eigenvectors;
double time1, time2;
char zz2[24] = "scf_evectors";
FILE *scf_evectors;

  time1 = MPI_Wtime();
  MPI_File fh;
  char buf2[110], xy[14] = "/scf_evec_sp";
  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xy);

  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  if (job->taskid == 0) {
  dimk = job->spin_dim * nbands;
  vector_size = dimk * dim1;
  AllocateComplexMatrix(&scf_eigenvectors,&dimk,&dim1,job);
  AllocateDoubleMatrix(&eigvec,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  scf_evectors = fopen(zz2, "rb");
  fseek(scf_evectors, dim1 * (fermi->bands[0] - 1) * sizeof(Complex),SEEK_SET);
  size_t result = fread(&scf_eigenvectors->a[0][0],sizeof(Complex),vector_size,scf_evectors);
  fclose(scf_evectors);
  for (i = 0; i < nbands; i++) {
    for (j = 0; j < dim1; j++) {
      eigvec->a[i][j] = (scf_eigenvectors->a[i][j]).real();
      //fprintf(file.out,"read write eigvec %3d %3d %10.4lf\n",i,j,eigvec->a[i][j]);
     }
    }
  DestroyComplexMatrix(&scf_eigenvectors,job);
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  MPI_File_write(fh, &eigvec->a[0][0], vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
  DestroyDoubleMatrix(&eigvec,job);
 }
  MPI_File_close(&fh);
  time2 = MPI_Wtime() - time1;
  //if (job->taskid == 0) printf("end read/write evecs %3d %f\n",job->taskid,time2);

}

void read_scf_GW_eigenvalues(double *eigenvalues, int seekpoint, int datasize, char *zz, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Read GW eigenvalues from disk                                                          *
  // ******************************************************************************************

int i, s;
double *evalues_temp;
size_t result;
FILE *evals;

  AllocateDoubleArray(&evalues_temp,&datasize,job);
  ResetDoubleArray(evalues_temp,&datasize);
  if (job->taskid == 0) {
  evals = fopen(zz, "rb");
  fseek(evals,seekpoint*sizeof(double),SEEK_SET);
  result = fread(evalues_temp,sizeof(double),job->spin_dim * datasize,evals);
  fclose(evals);
 } 
  MPI_Allreduce(evalues_temp,eigenvalues,job->spin_dim * datasize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&evalues_temp,&datasize,job);

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Eigenvalues\n");
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < datasize; i++) {
      fprintf(file.out,"%5d %5d %10.4lf\n",s,i,eigenvalues[s * datasize + i] * au_to_eV);
     }
    }
   }

}

void write_SCF_GW_eigenvalues(double *eigenvalues, int seekpoint, int size, char *zz, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Write SCF or GW eignvalues to MPI file SCF_eval or GW_eval                             *
  // ******************************************************************************************

int i;
char buf1[110];
MPI_File fh;

  strcpy(buf1,file.scf_eigvec);
  strcat(buf1,zz);

  MPI_File_open(MPI_COMM_WORLD,buf1, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh) ;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  MPI_File_write(fh, eigenvalues, size, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"Writing %s at %d\n", buf1, seekpoint);
  for (i = 0; i < size; i++)
  fprintf(file.out,"%5d %10.4lf\n",i,eigenvalues[i] * au_to_eV);
 }

}

void read_SCF_GW_eigenvalues(double *eigenvalues, int seekpoint, int size, char *zz, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Read SCF, GW or CAS eigenvalues from disk                                              *
  // ******************************************************************************************

int i;
char buf1[110];
MPI_File eh;

  strcpy(buf1, file.scf_eigvec);
  strcat(buf1, zz);

  MPI_File_open(MPI_COMM_WORLD, buf1, MPI_MODE_RDONLY, MPI_INFO_NULL, &eh) ;
  MPI_File_seek(eh, seekpoint * sizeof(double), MPI_SEEK_SET) ;
  MPI_File_read(eh, eigenvalues, size, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_close(&eh);

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"reading %s at %d\n",buf1,seekpoint);
  for (i = 0; i < size; i++)
  fprintf(file.out,"%5d %10.4lf\n",i,eigenvalues[i] * au_to_eV);
 }

}

void generate_temp4(double *temp4, double *integral_buffer2, FERMI *fermi, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int i, m, n;
int a1, nd6;
int j1, l1, l2;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int itr, rem;
int begin_j[job->numtasks], end_j[job->numtasks];
int num_proc = job->numtasks < dim1 ? job->numtasks : dim1;
int offset, offset1, offset2;
int *offset_j1;
char xc[22] = "/cas_eigenvectors_mpi";
char bufcas[120];
double time1, time2;
double time3, time4;
DoubleMatrix *cas_eigenvectors;
MPI_File fh;

  strcpy(bufcas,file.scf_eigvec);
  strcat(bufcas,xc);

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with job->rpa_lim CASIDA eigenvectors            *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  time3 = MPI_Wtime();

  itr = job->rpa_lim / nvir;
  rem = job->rpa_lim - itr * nvir;

  AllocateIntArray(&offset_j1,&dim1,job);
  mpi_begin_end(begin_j, end_j, dim1, num_proc, job, file);
  for (i = num_proc; i < job->numtasks; i++) begin_j[i] = 0;
  for (i = num_proc; i < job->numtasks; i++) end_j[i] = 0;
  for (i = 0; i < dim1; i++) offset_j1[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j1[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1]; }}
  if (end_j[job->taskid] == begin_j[job->taskid]) { itr = 0; rem = 0; }

  MPI_File_open(MPI_COMM_WORLD,bufcas,MPI_MODE_RDWR,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;

  for (l1 = 0; l1 < itr; l1++) {
    AllocateDoubleMatrix(&cas_eigenvectors,&nvir,&ntransitions,job);
    ResetDoubleMatrix(cas_eigenvectors);
    //MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
    MPI_File_read(fh, &cas_eigenvectors->a[0][0], nvir * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
    for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
      nd6 = atoms_ax->bfnnumb_sh[j1];
      offset1 = job->rpa_lim * offset_j1[j1];
      offset2 = offset_j1[j1] * nbands * nbands;
      for (l2 = 0; l2 < nvir; l2++) {
        for (m = 0; m < nocc; m++) {
          for (n = 0; n < nvir; n++) {
            for (a1 = 0; a1 < nd6; a1++) {
              temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1] += \
              integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1] * cas_eigenvectors->a[l2][m * nvir + n];
             }
            }
           }
          }
         }
        DestroyDoubleMatrix(&cas_eigenvectors,job);
       }

  if (rem > 0) {
  AllocateDoubleMatrix(&cas_eigenvectors,&nvir,&ntransitions,job);
  ResetDoubleMatrix(cas_eigenvectors);
  MPI_File_read(fh, &cas_eigenvectors->a[0][0], rem * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = offset_j1[j1] * job->rpa_lim;
    offset2 = offset_j1[j1] * nbands * nbands;
    for (l2 = 0; l2 < rem; l2++) {
      for (m = 0; m < nocc; m++) {
        for (n = 0; n < nvir; n++) {
          for (a1 = 0; a1 < nd6; a1++) {
            temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1] += \
            integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1] * cas_eigenvectors->a[l2][m * nvir + n];
           }
          }
         }
        }
       }
      DestroyDoubleMatrix(&cas_eigenvectors,job);
     }

  MPI_File_close(&fh);
  DestroyIntArray(&offset_j1,&dim1,job);
  time4 = MPI_Wtime() - time3; 
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("generate_temp4                             %10.2f\n",time2);

}
