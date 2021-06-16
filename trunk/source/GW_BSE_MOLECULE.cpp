#include <cstring>
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

fermi->nkunique = 1;
allocate_fermi(fermi,atoms,job,file);
fermi->occupied[0] = fermi->homo[0] - fermi->bands[0] + 1;

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

  ////read_write_scf_eigenvectors(fermi,atoms,job,file);

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
  //MPI_Bcast(&V_inv->a[0][0],dim1ax*dim1ax,MPI_DOUBLE,0,MPI_COMM_WORLD);
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
  
  bse_hamiltonian_in_core1(&ictxt,&nbsize,Ham_buffer1,Ham_buffer2,integral_buffer1,integral_buffer2,fermi,atom_p,atoms,\
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
  self_energy_diagonal_in_core1(&ictxt,&nbsize,begin_j,end_j,offset_j,integral_buffer1,integral_buffer2,\
  atoms_ax,fermi,job,file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("self_energy                                %10.2f\n",time4);
  if (job->gw_int == 2) {
  DestroyDoubleArray(&integral_buffer1,&dim_send[job->taskid],job);
  DestroyDoubleArray(&integral_buffer2,&dim_send[job->taskid],job);
 }
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

fermi->nkunique = 1;
allocate_fermi(fermi,atoms,job,file);
fermi->occupied[0] = fermi->homo[0] - fermi->bands[0] + 1;

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

  ////read_write_scf_eigenvectors(fermi,atoms,job,file);

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
  //gw_bse_molecule_integrals_in_core(integral_buffer1,integral_buffer2,&ictxt,&nbsize,fermi,atoms,atom_p,shells, \
  gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  generate_integral_buffers_molecule_ija(integral_buffer1,integral_buffer2,&ictxt,&nbsize,fermi,atoms,atom_p,shells, \
  gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  bse_hamiltonian_in_core1(&ictxt,&nbsize,Ham_buffer1,Ham_buffer2,integral_buffer1,integral_buffer2,fermi,atom_p,atoms, \
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
  //bse_evalues  = fopen(yy, "rb");
  //size_t result = fread(bse_eigenvalues, sizeof(double), ntransitions, bse_evalues);
  //fclose(bse_evalues);
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

void rpa_molecule1(FERMI* fermi, ATOM* atoms, ATOM_TRAN* atom_p, SHELL* shells, GAUSSIAN* gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates RPA excitations                                                     *
  // ******************************************************************************************
  
  // ******************************************************************************************
  // * Allocate fermi structure                                                               *
  // ******************************************************************************************

fermi->nkunique = 1;
allocate_fermi(fermi,atoms,job,file);
fermi->occupied[0] = fermi->homo[0] - fermi->bands[0] + 1;

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
  if (job->taskid == 0) printf("mpA %7d nqA %7d blocksize %7d ntransitions %7d\n",mpA,nqA,nbsize,ntransitions);

  ////read_write_scf_eigenvectors(fermi,atoms,job,file);

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
  bse_hamiltonian_in_core1(&ictxt,&nbsize,Ham_buffer1,Ham_buffer2,integral_buffer1,integral_buffer2,fermi,atom_p,atoms, \
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

void bse_hamiltonian_in_core1(int *ictxt, int *nbsize, double *Ham_buffer1, double *Ham_buffer2, double *integral_buffer1, double *integral_buffer2, FERMI* fermi, ATOM_TRAN *atom_p, ATOM *atoms, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

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
//char zz4[24] = "scf_evalues";
//char zz5[24] = "GW_evalues";

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

  ////read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);
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
  hamiltonian_in_core_screened_ij_ab3(ictxt,nbsize,begin_j,end_j,offset_j,GW_eigenvalues,Ham_buffer1,\
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
  if (job->taskid == 0) printf("bse_hamiltonian_in_core                    %10.2f\n",time2);
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
  block_cyclic_zero_triangle(&uplo,&nt,ictxt,nbsize,Ham_buffer1); // FIX uplo
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
  block_cyclic_to_linear(&nt,ictxt,nbsize,X2,xx,job,file);

  for (i = 0; i < ntransitions; i++) bse_eigenvalues[i] = eigenvalues[ntransitions - 1 - i];
  //for (i = 0; i < ntransitions; i++) bse_eigenvalues[i] = eigenvalues[i + nt - job->bse_lim];
  //if (job->taskid == 0) {
  //bse_evalues  = fopen(yy, "wb");
  //fwrite(bse_eigenvalues, sizeof(double), ntransitions, bse_evalues);
  //fflush(bse_evalues);
  //fclose(bse_evalues);
 //}

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
//char zz4[24] = "scf_evalues";
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

  ////read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);

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
////CHP
  if (buffer_size < 1 || buffer_size > 10000 * 10000) buffer_size = 1;
////CHP
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
////CHP
  DestroyDoubleArray(&eigvec_buffer,&buffer_size,job);
////CHP
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

// *
  int block_size = job->rpa_lim, div, rem;
////CHP
  //if (job->rpa_lim * ntransitions > 2000) {

  if (ntransitions > 20000) {
////CHP

  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
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
////CHP
long long global_pointer;
////CHP
  AllocateDoubleMatrix(&cas_eigenvectors,&block_size,&ntransitions,job);
  for (i = 0 ; i < div; i++) {
  ResetDoubleMatrix(cas_eigenvectors);
////CHP
  global_pointer = (i * block_size) * (ntransitions * sizeof(double));
  //fprintf(file.out,"%d %lld %d\n",i,global_pointer,(i * block_size) * (ntransitions * sizeof(double)));
  MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
  //MPI_File_seek(fh, i * block_size * ntransitions * sizeof(double), MPI_SEEK_SET) ;
////CHP
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
////CHP
  //MPI_File_seek(fh, i * block_size * ntransitions * sizeof(double), MPI_SEEK_SET) ;
  MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
  MPI_File_write(fh, &cas_eigenvectors->a[0][0], block_size * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
////CHP
 }
  ResetDoubleMatrix(cas_eigenvectors);
////CHP
  global_pointer = (div * block_size) * (ntransitions * sizeof(double));
  MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
  //MPI_File_seek(fh, div * block_size * ntransitions * sizeof(double), MPI_SEEK_SET) ;
////CHP
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
////CHP
  //global_pointer = (div * block_size) * (ntransitions * sizeof(double));
  //MPI_File_seek(fh, div * block_size * ntransitions * sizeof(double), MPI_SEEK_SET) ;
  //MPI_File_write(fh, &cas_eigenvectors->a[0][0], block_size * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
  MPI_File_write(fh, &cas_eigenvectors->a[0][0], rem * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
////CHP
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
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|                        CASIDA MATRIX EIGENVECTORS AND sqrt[EIGENVALUES] (eV)                            |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  ////print_real_eigenvector_matrix2(cas_eigenvectors, cas_eigenvalues, 8, 35, au_to_eV, file);
  fprintf(file.out,"| CAS GEN TIME    %9.2e | CAS DIAG TIME %9.2e |                         |                         |\n", \
  time2,time4);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 }
  else if (job->verbosity == 1) {
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
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

  //FILE *cas_evals;
  //char zz6[24] = "cas_evalues";
  //cas_evals = fopen(zz6, "wb");
  //for (s = 0; s < job->spin_dim; s++) {
    //fwrite(&cas_eigenvalues[s * dim1], sizeof(double), ntransitions, cas_evals);
    //fflush(cas_evals);
    //fclose(cas_evals);
   //}
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
  ////CHPDestroyDoubleArray(&eigvec_buffer,&buffer_size,job);
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
      //printf("m n I1 %3d   %3d %3d %3d   mpA %3d nqA %3d\n",job->taskid,m,n,I1,mpA,nqA);
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

void hamiltonian_in_core_screened_ij_ab3(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *GW_eigenvalues, double *Ham_buffer1, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

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
//char zz4[24] = "scf_evalues";
//char zz5[24] = "cas_evalues";

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;
  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
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
  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();

  char zz6[24] = "/evalfile";
  read_SCF_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz6, job, file);

  ////read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);

  char zz7[24] = "/cas_evalues";
  read_SCF_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz7, job, file);
  //for (i = 0; i < job->bse_lim; i++) fprintf(file.out,"1419 %5d %10.4lf\n",i,cas_eigenvalues[i] * au_to_eV);
  ////read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz5, job, file);
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

////CHP
//int kk  = 0;
////CHP
  time1 = MPI_Wtime();
  mdim = nbands * nbands * nvir;
  //fprintf(file.out,"hamiltonian_in_core_screened_ij_ab3  allocating nbands^2 %7d * nvir %7d = %7d\n", \
  nbands*nbands,nvir, nbands*nbands*nvir); fflush(file.out);
  AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l2 = 0; l2 < nvir; l2++) {
      for (m = 0; m < nbands; m++) {
        for (r = 0; r < nbands; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer[m * nbands * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
////CHP
      //      if (kk < 1000 && (M_buffer[m * nbands * nvir + r * nvir + l2] != M_buffer[m * nbands * nvir + r * nvir + l2] || \
      //          integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] !=  integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] || \
      //          temp4[(l1 * nvir + l2) * nd6 + offset1 + a1] != temp4[(l1 * nvir + l2) * nd6 + offset1 + a1])){
      //          fprintf(file.out,"%3d %3d %3d %3d %3d %f %f %f\n",j1,l2,m,r,a1,M_buffer[m * nbands * nvir + r * nvir + l2],\
      //          integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1],  temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]);
      //          kk++;
      //     }
      //
           }
          }
         }
        }
       }
  time2 += MPI_Wtime() - time1;

  AllocateDoubleArray(&M,&mdim,job);
  //fprintf(file.out,"hamiltonian_in_core_screened_ij_ab3  allocating mdim %7d\n", mdim); fflush(file.out);
  ResetDoubleArray(M,&mdim);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  DestroyDoubleArray(&M_buffer,&mdim,job);

////CHP
//int ii = 0;

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
  //      if (M[i * nbands * nvir + k * nvir + l2] != M[i * nbands * nvir + k * nvir + l2] || \
  //          M[k * nbands * nvir + i * nvir + l2] != M[k * nbands * nvir + i * nvir + l2] || 
  //          cas_eigenvalues[l1 * nvir + l2] != cas_eigenvalues[l1 * nvir + l2]) {
  //          sigma_factor = 0.0;
  //          if (ii < 1000) {
  //          fprintf(file.out,"i %3d k %3d l2 %3d l1 %3d Mik %f Mki %f cas_eval %f\n",\
  //          i,k,l2,l1,M[i * nbands * nvir + k * nvir + l2],M[k * nbands * nvir + i * nvir + l2],cas_eigenvalues[l1 * nvir + l2]);
  //          ii++;
  //       }
  //      }
////CHP
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

////CHP
//int jj = 0;

  for (i = 0; i < nbands; i++) {
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
 //       if (M[i * nbands * nvir + k * nvir + l2] != M[i * nbands * nvir + k * nvir + l2] || \
 //           M[k * nbands * nvir + i * nvir + l2] != M[k * nbands * nvir + i * nvir + l2] || 
 //           cas_eigenvalues[l1 * nvir + l2] != cas_eigenvalues[l1 * nvir + l2]) {
 //           sigma_factor = 0.0;
 //           if (jj < 1000) {
 //           fprintf(file.out,"i %3d k %3d l2 %3d l1 %3d Mik %f Mki %f cas_eval %f\n",\
 //           i,k,l2,l1,M[i * nbands * nvir + k * nvir + l2],M[k * nbands * nvir + i * nvir + l2],cas_eigenvalues[l1 * nvir + l2]);
 //           jj++;
 //        }
 //       }
////CHP
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }
  //for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);
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

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);
  if (job->taskid == 0 ) printf("Contract                                   %10.2f\n",time8);
  if (job->taskid == 0) {
  fprintf(file.out,"GW Eigenvalues: %4d RPA states used\n",job->bse_lim);
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
//char zz4[24] = "scf_evalues";
//char zz5[24] = "cas_evalues";

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
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
  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();

  char zz6[24] = "/evalfile";
  read_SCF_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz6, job, file);

  ////read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);

  char zz7[24] = "/cas_evalues";
  read_SCF_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz7, job, file);
  //for (i = 0; i < job->bse_lim; i++) fprintf(file.out,"1649 %5d %10.4lf\n",i,cas_eigenvalues[i] * au_to_eV);
  ////read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz5, job, file);
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
  //mdim = nbands * nbands * nvir;
  //AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  ResetDoubleArray(M,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
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

  //AllocateDoubleArray(&M,&mdim,job);
  //ResetDoubleArray(M,&mdim);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  //DestroyDoubleArray(&M_buffer,&mdim,job);

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;

        //fprintf(file.out,"%3d %3d %3d %10.3e %10.3e %10.3e %10.3e\n",i,k,l2,\
        sigma_factor,M[i * nbands * nvir + k * nvir + l2], M[k * nbands * nvir + i * nvir + l2], denom);

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
  //for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);
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

  //DestroyDoubleArray(&M,&mdim,job);

 } // close loop on l1
  DestroyDoubleArray(&M_buffer,&mdim,job);
  DestroyDoubleArray(&M,&mdim,job);

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);
  if (job->taskid == 0 ) printf("Contract                                   %10.2f\n",time8);
  if (job->taskid == 0) {
  fprintf(file.out,"GW Eigenvalues: %4d RPA states used\n",job->bse_lim);
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

void self_energy_diagonal_in_core1(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

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
double *cas_eigenvalues, *scf_eigenvalues, *GW_eigenvalues;
double Sigma[nbands], dSigma_dE[nbands], sigma_factor, denom;
//char zz4[24] = "scf_evalues";
//char zz5[24] = "cas_evalues";

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
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
  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();
  char zz6[24] = "/evalfile";
  read_SCF_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz6, job, file);
  char zz7[24] = "/cas_evalues";
  read_SCF_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz7, job, file);
  ////read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);
  ////read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz5, job, file);
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
  //mdim = nbands * nbands * nvir;
  //AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  ResetDoubleArray(M,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
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

  //AllocateDoubleArray(&M,&mdim,job);
  //ResetDoubleArray(M,&mdim);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  //DestroyDoubleArray(&M_buffer,&mdim,job);

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;

        //fprintf(file.out,"%3d %3d %3d %10.3e %10.3e %10.3e %10.3e\n",i,k,l2,\
        sigma_factor,M[i * nbands * nvir + k * nvir + l2], M[k * nbands * nvir + i * nvir + l2], denom);

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
  //for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);
  time6 += MPI_Wtime() - time5;

 } // close loop on l1
  DestroyDoubleArray(&M_buffer,&mdim,job);
  DestroyDoubleArray(&M,&mdim,job);

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);
  if (job->taskid == 0) {
  fprintf(file.out,"GW Eigenvalues: %4d RPA states used\n",job->bse_lim);
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
  //char zz4[24] = "GW_evalues";
  //FILE *evals = fopen(zz4, "wb");
  //fseek(evals,0,SEEK_SET);
  //fwrite(GW_eigenvalues,sizeof(double),job->spin_dim * nbands,evals);
  //fclose(evals);
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
//fprintf(file.out,"Writing %3d eigenvectors to %s\n",fermi->bands[0]-1,buf2);


  time1 = MPI_Wtime();
  MPI_File fh;
  char buf2[110], xy[14] = "/scf_evec_spk";
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

void generate_temp4_in_core(double *temp4, double *integral_buffer2, FERMI *fermi, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

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
int num_eigenvectors;
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
  // * Contract V <alpha|beta><beta|occ-vir> with CASIDA eigenvectors                         *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  time3 = MPI_Wtime();

  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  num_eigenvectors = job->bse_lim;
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;

  AllocateIntArray(&offset_j1,&dim1,job);
  mpi_begin_end(begin_j, end_j, dim1, num_proc, job, file);
  for (i = num_proc; i < job->numtasks; i++) begin_j[i] = 0;
  for (i = num_proc; i < job->numtasks; i++) end_j[i] = 0;
  for (i = 0; i < dim1; i++) offset_j1[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j1[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1]; }}
  if (end_j[job->taskid] == begin_j[job->taskid]) { itr = 0; rem = 0; }

  //printf("task %3d %3d %3d %3d %3d\n",job->taskid,nvir, ntransitions,itr,rem);
  MPI_File_open(MPI_COMM_WORLD,bufcas,MPI_MODE_RDWR,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;

////CHP
//int kk = 0;
////CHP
  for (l1 = 0; l1 < itr; l1++) {
    //fprintf(file.out,"generate_temp4_in_core %3d allocating nvir %7d * ntransitions %7d = %7d\n",l1, \
    nvir,ntransitions, nvir * ntransitions); fflush(file.out);
    AllocateDoubleMatrix(&cas_eigenvectors,&nvir,&ntransitions,job);
    ResetDoubleMatrix(cas_eigenvectors);
    //MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
    MPI_File_read(fh, &cas_eigenvectors->a[0][0], nvir * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
    for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
      nd6 = atoms_ax->bfnnumb_sh[j1];
      offset1 = job->bse_lim * offset_j1[j1];
      offset2 = offset_j1[j1] * nbands * nbands;
      for (l2 = 0; l2 < nvir; l2++) {
        for (m = 0; m < nocc; m++) {
          for (n = 0; n < nvir; n++) {
            for (a1 = 0; a1 < nd6; a1++) {
              temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1] += \
              integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1] * cas_eigenvectors->a[l2][m * nvir + n];
////CHP
 //           if (kk < 1000 && (temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1] != temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1] || \
 //             integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1] != integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1] || \
 //              cas_eigenvectors->a[l2][m * nvir + n] != cas_eigenvectors->a[l2][m * nvir + n])) {
 //              fprintf(file.out,"main %3d %3d %3d %3d %3d %3d %f %f %f\n",l1,j1,l2,m,n,a1,temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1],\
 //              integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1],cas_eigenvectors->a[l2][m * nvir + n]);
 //              kk++;
 //             }
//////CHP
             }
            }
           }
          }
         }
        DestroyDoubleMatrix(&cas_eigenvectors,job);
       }

  if (rem > 0) {
  //fprintf(file.out,"generate_temp4_in_core  allocating nvir %7d * ntransitions %7d = %7d\n", \
  nvir,ntransitions, nvir * ntransitions); fflush(file.out);
  AllocateDoubleMatrix(&cas_eigenvectors,&nvir,&ntransitions,job);
  ResetDoubleMatrix(cas_eigenvectors);
  MPI_File_read(fh, &cas_eigenvectors->a[0][0], rem * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = offset_j1[j1] * job->bse_lim;
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

      //for (a1 = 0; a1 < nocc * nvir * 6; a1++) fprintf(file.out,"TEMP4a %3d %10.4lf\n", a1, temp4[a1]);


  MPI_File_close(&fh);
  DestroyIntArray(&offset_j1,&dim1,job);
  time4 = MPI_Wtime() - time3; 
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("generate_temp4_in_core                     %10.2f\n",time2);

}

/*
void optical_spectrum_molecule(FERMI* fermi, ATOM* atoms, ATOM_TRAN* atom_p, int *numfrag, int *natom, int nat[][2], SHELL* shells, GAUSSIAN* gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

fermi->nkunique = 1;
allocate_fermi(fermi,atoms,job,file);
fermi->occupied[0] = fermi->homo[0] - fermi->bands[0] + 1;

int i, j, k, i1, j1, i3, j3, l, m, n, s, t;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int nt = ntransitions;
int dim4 = nbands * 3;
int seek, count;
double E, E_tot, energy, range = job->energy_range[1] - job->energy_range[0];
double increment = range / (double) (job->npoints + 1);
double *scf_eigenvalues, *GW_eigenvalues, *bse_eigenvalues;
//double del = 0.002;
//double del = 0.00367; // 0.1 eV
char xx[20] = "/bse_eigenvectors";
char yy[20] = "bse_eigenvalues";
char zz4[24] = "scf_evalues";
char zz5[24] = "GW_evalues";
char buf2[110];
FILE *bse_evalues, *scf_evalues, *spect;
DoubleMatrix *bse_eigenvectors, *M_k, *M_x, *tmp;
MPI_File fh;

  // ******************************************************************************************
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************
 
  AllocateDoubleArray(&scf_eigenvalues,&nbands,job);
  ResetDoubleArray(scf_eigenvalues,&nbands);
  read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);
  AllocateDoubleArray(&bse_eigenvalues,&ntransitions,job);
  ResetDoubleArray(bse_eigenvalues,&ntransitions);
  read_scf_GW_eigenvalues(bse_eigenvalues, 0, ntransitions, yy, job, file);
  //AllocateDoubleArray(&GW_eigenvalues,&nbands,job);
  //ResetDoubleArray(GW_eigenvalues,&nbands);
  //read_scf_GW_eigenvalues(GW_eigenvalues, fermi->bands[0] - 1, nbands, zz5, job, file);

  // ******************************************************************************************
  // * Allocate memory for BSE eigenvectors                                                   *
  // * Read eigenvectors from disk                                                            *
  // ******************************************************************************************

  if (job->bse_lim == 0) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  if (job->taskid == 0)
  AllocateDoubleMatrix(&bse_eigenvectors,&job->bse_lim,&ntransitions,job);
  //AllocateDoubleMatrix(&bse_eigenvectors,&ntransitions,&ntransitions,job);
  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xx);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  if (job->taskid == 0 && job->bse_tda == 0) {
  //printf("bse_eigenvectors directory: %s %s\n",buf2,file.scf_eigvec);
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  //for (i = 0; i < nt; i++) 
  for (i = 0; i < job->bse_lim; i++) 
  //MPI_File_read(fh, &bse_eigenvectors->a[nt - 1 - i][0], nt, MPI_DOUBLE, MPI_STATUS_IGNORE);
  //MPI_File_read(fh, &bse_eigenvectors->a[i][0], nt, MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_read(fh, &bse_eigenvectors->a[job->bse_lim - 1 - i][0], nt, MPI_DOUBLE, MPI_STATUS_IGNORE);
  for (j = 0; j < bse_eigenvectors->iRows; j++) { 
    for (k = 0; k < bse_eigenvectors->iCols; k++) { 
      //fprintf(file.out,"%3d %3d %10.4f\n",j,k,bse_eigenvectors->a[j][k]);
      bse_eigenvectors->a[j][k] /= two * sqrt(bse_eigenvalues[j]); 
     }
    } 
   }
  else if (job->taskid == 0 && job->bse_tda == 1) {
//CHP2018
long long lim, memsize;
lim = job->bse_lim;
if (lim > 25000) lim = 25000;
//long long lim, nt1;
//lim = 25000;
memsize = lim * nt;
//printf("%d %d %lld %lld %lld\n",job->bse_lim,nt,lim,nt1,memsize);
//CHP2018
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  MPI_File_read(fh, &bse_eigenvectors->a[0][0], memsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
  //CHP2018MPI_File_read(fh, &bse_eigenvectors->a[0][0], job->bse_lim * nt, MPI_DOUBLE, MPI_STATUS_IGNORE);
  //MPI_File_read(fh, &bse_eigenvectors->a[0][0], nt * nt, MPI_DOUBLE, MPI_STATUS_IGNORE);
//CHP2018
for (i = lim; i < job->bse_lim; i++) {
  for (j = 0; j < nt; j++) { 
    bse_eigenvectors->a[i][j] = k_zero;
   }
  }
//CHP2018
  }
  MPI_File_close(&fh);
  //fprintf(file.out,"|                                    BSE EIGENVECTOR MATRIX (eV)                                      |\n");
  //if (job->taskid == 0)  
  //print_real_eigenvector_matrix2(bse_eigenvectors, bse_eigenvalues, 8, 50, au_to_eV, file);

  // ******************************************************************************************
  // * Calculate dipole matrix elements                                                       *
  // ******************************************************************************************

  //if (job->taskid == 0)
  //for (i = 0; i < bse_eigenvectors->iRows; i++) { 
  //fprintf(file.out,"%3d %10.4f ",i,bse_eigenvalues[i]);
  //for (j = 0; j < bse_eigenvectors->iCols; j++) { 
  //fprintf(file.out,"%8.3f",bse_eigenvectors->a[i][j]); 
  //}
  //fprintf(file.out,"\n");
  //}

  if (job->taskid == 0) {
  AllocateDoubleMatrix(&M_x,&dim4,&nbands,job);
  ResetDoubleMatrix(M_x);
 }
  dipole_matrix_elements_molecule(M_x,fermi,atoms,atom_p,shells,gaussians,crystal,symmetry,R,R_tables,G,job,file);

  // ******************************************************************************************
  // * Allocate memory for optical matrix elements                                            *
  // ******************************************************************************************

  if (job->taskid == 0) {
  //double bse_matrix_element[job->field_dirs][ntransitions];
  double free_particle_matrix_element[job->field_dirs][ntransitions];
  double bse_matrix_element[job->field_dirs][job->bse_lim];
  double bse_spectrum[job->spin_dim][job->field_dirs][job->npoints + 1];
  double free_particle_spectrum[job->spin_dim][job->field_dirs][job->npoints + 1];
  double E_trans[ntransitions];
  double TRK_sum_rule[job->field_dirs];
  double TRK_sum_rule1[job->field_dirs];
  char field[80];

  for (i = 0; i < job->spin_dim; i++) {
    for (j3 = 0; j3 < job->field_dirs; j3++) { 
      for (j = 0; j < job->npoints + 1; j++) {
        bse_spectrum[i][j3][j] = k_zero;
        free_particle_spectrum[i][j3][j] = k_zero;
       }
      }
     }

    for (j3 = 0; j3 < job->field_dirs; j3++) TRK_sum_rule[j3] = k_zero;
    for (j3 = 0; j3 < job->field_dirs; j3++) TRK_sum_rule1[j3] = k_zero;

  // ******************************************************************************************
  // * Generate Bethe-Salpeter and single-particle dielectric functions                       *
  // ******************************************************************************************
  
  //fprintf(file.out,"BSE Eigenvalues\n");
  //for (s = 0; s < job->spin_dim; s++) {
  //for (i = 0; i < nbands; i++) {
  //fprintf(file.out,"%5d %5d %10.4lf\n",s,i,bse_eigenvalues[s * nbands + i] * au_to_eV);
  //}
  //}
  //print_real_matrix2(M_x,0,8,1.0,file);

   sprintf(field," %d     %4.1lf %4.1lf %4.1lf    %4.1lf %4.1lf %4.1lf   %4.1lf %4.1lf %4.1lf",job->field_dirs, \
   job->e_field[0].comp1, job->e_field[0].comp2, job->e_field[0].comp3, \
   job->e_field[1].comp1, job->e_field[1].comp2, job->e_field[1].comp3, \
   job->e_field[2].comp1, job->e_field[2].comp2, job->e_field[2].comp3);
   fprintf(file.out,"|                                      LARGEST OSCILLATOR STRENGTHS                                       |\n");
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"| MODE NUMBER   ENERGY (eV) | FIELDS   %66s |\n", field);
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

  for (s = 0; s < job->spin_dim; s++) {

  for (j3 = 0; j3 < job->field_dirs; j3++) {
    for (t = 0; t < ntransitions; t++) {
      free_particle_matrix_element[j3][t] = k_zero;
     }
    }

  for (j3 = 0; j3 < job->field_dirs; j3++) {
    for (t = 0; t < job->bse_lim; t++) {
      bse_matrix_element[j3][t] = k_zero;
     }
    }

    count = 0;
    for (m = 0; m < fermi->occupied[0]; m++) {
      for (n = fermi->occupied[0]; n < nbands; n++) {
        for (j3 = 0; j3 < job->field_dirs; j3++) {
          free_particle_matrix_element[j3][count] = \
         (job->e_field[j3].comp1 * M_x->a[m][n] + \
          job->e_field[j3].comp2 * M_x->a[nbands + m][n] + \
          job->e_field[j3].comp3 * M_x->a[2 * nbands + m][n]);
          E_trans[count] = scf_eigenvalues[n] - scf_eigenvalues[m];
          //if (fabs(free_particle_matrix_element[j3][count]) > 1.5)
          //fprintf(file.out,"%3d %3d %3d %3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",m,n,count,j3,\
          free_particle_matrix_element[j3][count], bse_matrix_element[0][54],bse_matrix_element[1][54],bse_matrix_element[2][54],\
          scf_eigenvalues[m] * au_to_eV, scf_eigenvalues[n] * au_to_eV, E_trans[count] * au_to_eV,bse_eigenvectors->a[5][count],\
          bse_eigenvectors->a[6][count], bse_eigenvectors->a[7][count]);
          //for (t = 0; t < ntransitions; t++) {
          for (t = 0; t < job->bse_lim; t++) {
            bse_matrix_element[j3][t] += free_particle_matrix_element[j3][count] * bse_eigenvectors->a[t][count];
         }
        }
        count++; // counter for all transitions
       } // close loop on n
      } // close loop on m
      //for (t = 0; t < job->bse_lim; t++) {
      //fprintf(file.out,"bsemat %3d %10.4f %10.4f %10.4f\n",\
      t,bse_matrix_element[0][t],bse_matrix_element[1][t],bse_matrix_element[2][t]);
      //}

      for (j1 = 0; j1 < job->npoints; j1++) {
        E = increment * (double) j1 + job->energy_range[0];
        for (j3 = 0; j3 < job->field_dirs; j3++) {
          //for (t = 0; t < ntransitions; t++) {
          for (t = 0; t < job->bse_lim; t++) {
            //bse_spectrum[s][j3][j1] += bse_matrix_element[j3][t] * bse_matrix_element[j3][t] * \
            pi / del / E / E * exp(-(*(bse_eigenvalues + t) - E) * (*(bse_eigenvalues + t) - E) / del / del);
            //free_particle_spectrum[s][j3][j1] += free_particle_matrix_element[j3][t] * free_particle_matrix_element[j3][t] * \
            pi / del / E / E * exp(-(E_trans[t] - E) * (E_trans[t] - E) / del / del);
            //bse_spectrum[s][j3][j1] += job->spin_fac * two * bse_matrix_element[j3][t] * bse_matrix_element[j3][t] * \
            bse_eigenvalues[t] * \
            job->linewidth / pi / ((E - bse_eigenvalues[t]) * (E - bse_eigenvalues[t]) + job->linewidth * job->linewidth);
            bse_spectrum[s][j3][j1] += job->spin_fac * two * bse_matrix_element[j3][t] * bse_matrix_element[j3][t] * \
            bse_eigenvalues[t] * \
            job->linewidth / ((E - bse_eigenvalues[t]) * (E - bse_eigenvalues[t]) + job->linewidth * job->linewidth);
            free_particle_spectrum[s][j3][j1] += job->spin_fac * two * free_particle_matrix_element[j3][t] * \
            free_particle_matrix_element[j3][t] * E_trans[t] * \
            job->linewidth / ((E - E_trans[t]) * (E - E_trans[t]) + job->linewidth * job->linewidth);
            //free_particle_spectrum[s][j3][j1] += job->spin_fac * two * free_particle_matrix_element[j3][t] * \
            free_particle_matrix_element[j3][t] * E_trans[t] * \
            job->linewidth / pi / ((E - E_trans[t]) * (E - E_trans[t]) + job->linewidth * job->linewidth);
           }
          }
         }

        //for (t = 0; t < ntransitions; t++) {
        for (t = 0; t < job->bse_lim; t++) {
          for (j3 = 0; j3 < job->field_dirs; j3++) {
            TRK_sum_rule[j3] +=  job->spin_fac * two / three * bse_matrix_element[j3][t] * bse_matrix_element[j3][t] * \
            bse_eigenvalues[t];
           }
          }

        //for (t = 0; t < ntransitions; t++) {
        for (t = 0; t < job->bse_lim; t++) {
          if ((bse_matrix_element[0][t] * bse_matrix_element[0][t] * bse_eigenvalues[t] > 0.05) || \
              (bse_matrix_element[1][t] * bse_matrix_element[1][t] * bse_eigenvalues[t] > 0.05) || \
              (bse_matrix_element[2][t] * bse_matrix_element[2][t] * bse_eigenvalues[t] > 0.05)) {
          //if ((bse_matrix_element[0][t] * bse_matrix_element[0][t] * bse_eigenvalues[t] > 0.00) || \
              (bse_matrix_element[1][t] * bse_matrix_element[1][t] * bse_eigenvalues[t] > 0.00) || \
              (bse_matrix_element[2][t] * bse_matrix_element[2][t] * bse_eigenvalues[t] > 0.00)) {
          fprintf(file.out,"| %4d           %10.4lf ", t + 1, bse_eigenvalues[t] * au_to_eV);
          for (j3 = 0; j3 < job->field_dirs; j3++) {
            fprintf(file.out,"|              %10.4lf ", \
            job->spin_fac * two / three * bse_matrix_element[j3][t] * bse_matrix_element[j3][t] * bse_eigenvalues[t]);
           }
          fprintf(file.out,"|\n");
         }
        }

   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"| TRK SUM RULE   %10.4lf |              %10.4lf |              %10.4lf |              %10.4lf |\n", \
   TRK_sum_rule[0] + TRK_sum_rule[1] + TRK_sum_rule[2], TRK_sum_rule[0], TRK_sum_rule[1], TRK_sum_rule[2]);
   fprintf(file.out,"===========================================================================================================\n");
   if (job->taskid == 0 && job->field_dirs < 3)
   fprintf(file.out,"WARNING: Thomas-Reiche-Kun Sum Rule: %3d out of 3 components calculated.\n", job->field_dirs);
   fflush(file.out);

     } // close loop on s

    spect = fopen("Free_particle_optical_spectrum.dat", "w");
    if (spect == NULL) { fprintf(file.out, "cannot open file Free_particle_optical_spectrum.dat in bethe_salpeter\n"); exit(1); }
        for (j = 0; j < job->npoints; j++) {
         energy = job->energy_range[0] + range * (double) j / (double) job->npoints;
          fprintf(spect, "%10.4e   ", energy * au_to_eV);
           for (l = 0; l < job->spin_dim; l++) {
            for (j3 = 0; j3 < job->field_dirs; j3++) {
            fprintf(spect, "%12.4e", free_particle_spectrum[l][j3][j]);
           }
          }
         fprintf(spect,"\n");
        }
       fflush(spect);
       fclose(spect);

    if      (job->bse_tda == 0 && job->bse_ham == 0) spect = fopen("TDHF_optical_spectrum.dat", "w");
    else if (job->bse_tda == 0 && job->bse_ham == 1) spect = fopen("BSE_optical_spectrum.dat", "w");
    else if (job->bse_tda == 1 && job->bse_ham == 0) spect = fopen("TDHF_TDA_optical_spectrum.dat", "w");
    else if (job->bse_tda == 1 && job->bse_ham == 1) spect = fopen("BSE_TDA_optical_spectrum.dat", "w");
    if (spect == NULL) { fprintf(file.out, "cannot open spectrum file in bethe_salpeter routine\n"); exit(1); }
        for (j = 0; j < job->npoints; j++) {
         energy = job->energy_range[0] + range * (double) j / (double) job->npoints;
          fprintf(spect, "%10.4e   ", energy * au_to_eV);
           for (l = 0; l < job->spin_dim; l++) {
            for (j3 = 0; j3 < job->field_dirs; j3++) {
            fprintf(spect, "%12.4e", bse_spectrum[l][j3][j]);
           }
          }
         fprintf(spect,"\n");
        }
       fflush(spect);
       fclose(spect);
       DestroyDoubleArray(&bse_eigenvalues,&ntransitions,job);
       DestroyDoubleArray(&scf_eigenvalues,&nbands,job);
       DestroyDoubleMatrix(&M_x,job);
       DestroyDoubleMatrix(&bse_eigenvectors,job);
      } // close if (job->taskid == 0)

}

void gw_bse_molecule_integrals(int *ictxt, int *nbsize, FERMI* fermi, ATOM* atoms, ATOM_TRAN* atom_p, SHELL* shells, GAUSSIAN* gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
double time1, time2;
double time3, time4;

  time1 = MPI_Wtime();

  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in gw_bse_molecule_integrals. HOMO level not in range of bands. ntransitions = %d %d %d\n",\
  ntransitions,nocc,nvir);
  MPI_Finalize();
  exit(0);
 }

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate three centre integrals needed for matrix elements                             *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  generate_density_fitting_intermediates_spk(ictxt,nbsize,fermi,atom_p,atoms,shells,gaussians,atoms_ax,shells_ax,\
  gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  time4 = MPI_Wtime();
  //if (job->taskid == 0) printf("generate_density_fitting_intermediates_spk %10.2f\n",time4 - time3);

  time3 = MPI_Wtime();
  integrals_occ_vir_occ_vir(ictxt,nbsize,0,fermi,atom_p,atoms,shells,gaussians,atoms_ax,shells_ax,\
  gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  time4 = MPI_Wtime();
  //if (job->taskid == 0) printf("block_cyclic %3d integrals_occ_vir_occ_vir %10.2f\n",0,time4 - time3);

  //time3 = MPI_Wtime();
  //integrals_occ_vir_occ_vir(ictxt,nbsize,1,fermi,atom_p,atoms,shells,gaussians,atoms_ax,shells_ax,\
  gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  //time4 = MPI_Wtime();
  //if (job->taskid == 0) printf("block_cyclic %3d integrals_occ_vir_occ_vir %10.2f\n",1,time4 - time3);

  time3 = MPI_Wtime();
  integrals_occ_occ_vir_vir(ictxt,nbsize,0,fermi,atom_p,atoms,shells,gaussians,atoms_ax,shells_ax,\
  gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  time4 = MPI_Wtime();
  //if (job->taskid == 0) printf("block_cyclic %3d integrals_occ_occ_vir_vir %10.2f\n",0,time4 - time3);

  time3 = MPI_Wtime();
  integrals_occ_occ_vir_vir(ictxt,nbsize,1,fermi,atom_p,atoms,shells,gaussians,atoms_ax,shells_ax,\
  gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  time4 = MPI_Wtime();
  //if (job->taskid == 0) printf("block_cyclic %3d integrals_occ_occ_vir_vir %10.2f\n",1,time4 - time3);

  time3 = MPI_Wtime();
  integrals_self_energy(ictxt,nbsize,fermi,atom_p,atoms,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,\
  symmetry,R,R_tables,G,job,file);
  time4 = MPI_Wtime();
  //if (job->taskid == 0) printf("integrals_self_energy                      %10.2f\n",time4 - time3);

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("integrals                                  %10.2f\n",time2);

}

void gw_bse_molecule_integrals_in_core(double *integral_buffer1, double *integral_buffer2, int *ictxt, int *nbsize, FERMI *fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int i, j, k, p;
int j1, j2, jj, kk;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int begin_j[job->numtasks], end_j[job->numtasks];
int *dim_send, *dim_recv, *offset_j;
int offset, offset1, dim3;
int num_proc = job->numtasks < dim1 ? job->numtasks : dim1;
double time1, time2, time3, time4;
double *integral_buffer;
char buf2[110], xy[14] = "/scf_evec_spk";
DoubleMatrix *V_inv;
MPI_File fh;

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xy);
  
  // ******************************************************************************************
  // * Generate inverse of Coulomb matrix                                                     *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  AllocateDoubleMatrix(&V_inv,&dim1ax,&dim1ax,job);
  //allocate_INT_1E(&one_ints, dim, Function, job, file);
  //fock_element_1e1(&one_ints, dim, &pair_p, pair_p.nump, Function, R, G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  //if (job->taskid == 0) 
  generate_coulomb_matrix_inverse(V_inv,atom_p,atoms,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  //free_INT_1E(&one_ints, Function, job, file);
  //MPI_Bcast(&V_inv->a[0][0],dim1ax*dim1ax,MPI_DOUBLE,0,MPI_COMM_WORLD);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("generate_coulomb_matrix_inverse            %10.2f\n",time4);

  // ******************************************************************************************
  // * Generate three centre integrals <ij|alpha> and <beta|V|kl> in core                     *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  AllocateIntArray(&dim_send,&job->numtasks,job);
  AllocateIntArray(&dim_recv,&job->numtasks,job);
  AllocateIntArray(&offset_j,&dim1,job);
  mpi_begin_end(begin_j, end_j, dim1, num_proc, job, file);
  for (i = num_proc; i < job->numtasks; i++) begin_j[i] = 0;
  for (i = num_proc; i < job->numtasks; i++) end_j[i] = 0;
  for (i = 0; i < dim1; i++) offset_j[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1] * nbands * nbands; }}
  for (i = 0; i < job->numtasks; i++) { dim_send[i] = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) {
  dim_send[i] += atoms_ax->bfnnumb_sh[j1] * nbands * nbands; }}
  AllocateDoubleArray(&integral_buffer,&dim_send[job->taskid],job);
  ResetDoubleArray(integral_buffer,&dim_send[job->taskid]);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
  contract_integrals_molecule_ij_alpha(&j1,fh,&integral_buffer[offset_j[j1]],fermi,atom_p,atoms,shells,gaussians,atoms_ax, \
  shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  //integrals_occ_occ_vir_vir1(&j1,fh,&integral_buffer[offset_j[j1]],fermi,atom_p,atoms,shells,gaussians,atoms_ax, \
  shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  } // close loop on j1
  MPI_File_close(&fh);
  MPI_Barrier(MPI_COMM_WORLD);
  initialise_ring_rotate(num_proc,&dim3,dim_send,dim_recv,job,file);
  for (i = 0; i < dim_send[job->taskid]; i++) integral_buffer1[i] = integral_buffer[i]; 
  ResetDoubleArray(integral_buffer2,&dim_send[job->taskid]);
  for (i = 0; i < num_proc; i++) {
    jj = (job->taskid - i) % num_proc;
    kk = jj < 0 ? jj + num_proc : jj;
    for (j1 = begin_j[kk]; j1 < end_j[kk]; j1++) {
      for (j2 = begin_j[job->taskid]; j2 < end_j[job->taskid]; j2++) {
        contract_coulomb_integrals(&j1,&j2,V_inv,&integral_buffer[offset_j[j1]],&integral_buffer2[offset_j[j2]],fermi,atom_p,\
        atoms,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
       }
      }
     ring_rotate_once1(&integral_buffer,i,num_proc,&dim3,dim_send,dim_recv,job,file);
    } // close loop over i
     DestroyDoubleArray(&integral_buffer,&dim_send[job->taskid],job);
     DestroyDoubleMatrix(&V_inv,job);
     DestroyIntArray(&dim_send,&job->numtasks,job);
     DestroyIntArray(&dim_recv,&job->numtasks,job);
     DestroyIntArray(&offset_j,&dim1,job);
     time2 = MPI_Wtime() - time1;
     if (job->taskid == 0) printf("integrals in core                          %10.2f\n",time2);

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
      //printf("m n I1 %3d   %3d %3d %3d   mpA %3d nqA %3d\n",job->taskid,m,n,I1,mpA,nqA);
      il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
      jl = *nbsize * (I1 / (*nbsize * npcol)) + I1 % *nbsize;
      Ham_buffer[il + mpA * jl] += eigenvalues[n] - eigenvalues[m];
     }
    }
   }

}

*/
/*
void hamiltonian_ia_jb(double *Hamiltonian, double factor, DoubleMatrix *V_inv, FERMI *fermi, int *ictxt, int *nbsize, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int il, jl, kl;
int j1, nd6, bfposa1, dim2, count;
int dim2a, dim2b;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
double time1, time2;
double *temp2a_buffer, *temp2b_buffer, *temp3_buffer;
char xx[4], yz[24] = "integrals_bc_ov1.", zz[24] = "integrals_bc_ov2.";
FILE *integrals_ov1, *integrals_ov2;

  sprintf(xx, "%d", job->taskid);
  strcat(yz,xx);
  strcat(zz,xx);

  time1 = MPI_Wtime();

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  dim2b = nqA * dim1ax;
  AllocateDoubleArray(&temp3_buffer,&dim2b,job);
  ResetDoubleArray(temp3_buffer,&dim2b);
  integrals_ov2 = fopen(zz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2b = nqA * nd6;
    AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
    ResetDoubleArray(temp2b_buffer,&dim2b);
    fread(temp2b_buffer,sizeof(double),dim2b,integrals_ov2);
    for (il = 0; il < dim1ax; il++) {
      for (jl = 0; jl < nd6; jl++) {
        for (kl = 0; kl < nqA; kl++) {
          temp3_buffer[il * nqA + kl] += factor * V_inv->a[il][bfposa1 + jl] * temp2b_buffer[jl * nqA + kl];
         }
        }
       }
      DestroyDoubleArray(&temp2b_buffer,&dim2b,job);
     }
      fclose(integrals_ov2);

  integrals_ov1 = fopen(yz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2a = nd6 * mpA;
    AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
    ResetDoubleArray(temp2a_buffer,&dim2a);
    fread(temp2a_buffer,sizeof(double),dim2a,integrals_ov1);
    for (jl = 0; jl < nd6; jl++) {
      for (il = 0; il < mpA; il++) {
        for (kl = 0; kl < nqA; kl++) {
          Hamiltonian[il + mpA * kl] += temp2a_buffer[jl * mpA + il] * temp3_buffer[(bfposa1 + jl) * nqA + kl];
         }
        }
       }
      DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
     }
      fclose(integrals_ov1);
      dim2b = nqA * dim1ax;
      DestroyDoubleArray(&temp3_buffer,&dim2b,job);
      time2 = MPI_Wtime() - time1;
      if (job->taskid == 0) printf("hamiltonian_ia_jb                          %10.2f\n",time2);

}

void hamiltonian_ij_ab(double *Hamiltonian, DoubleMatrix *V_inv, FERMI *fermi, int *ictxt, int *nbsize, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int a1, a2, il, jl, kl;
int i, j, j1, k, l, m, n, r, s, nd6, bfposa1, bfposa2;
int dim2, dim2a, dim2b, dim3, dim3a, dim6;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int begin_occ_row[job->numtasks], end_occ_row[job->numtasks];
int begin_occ_array[job->numtasks], end_occ_array[job->numtasks];
int sendtotal, sendcounts[job->numtasks], senddispls[job->numtasks];
int recvtotal, recvcounts[job->numtasks], recvdispls[job->numtasks];
int *sendindex1, *recvindex1;
int *sendindex2, *recvindex2;
char xx[4], xy[24] = "integrals_bc_oo.";
double time1, time2;
double time3, time4;
double *temp2a_buffer, *temp3_buffer;
double *send_buffer[job->numtasks];
double *recv_buffer[job->numtasks];
double *senddata, *recvdata;
FILE *integrals_temp2a;
IntMatrix *array_sizes;

  // ******************************************************************************************
  // * Set up Cblacs parameters                                                               *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  sprintf(xx, "%d", job->taskid);
  strcat(xy,xx);

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  mpi_begin_end_lower_triangle(begin_occ_row, end_occ_row, begin_occ_array, end_occ_array, nocc, nprow, job, file);

  // ******************************************************************************************
  // * Read <alpha|occ-occ> integrals and multiply in V <alpha|beta>                          *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();

  time3 = MPI_Wtime();
  integrals_temp2a = fopen(xy, "rb");
  dim3a = nocc * nocc * dim1ax;
  AllocateDoubleArray(&temp3_buffer,&dim3a,job);
  ResetDoubleArray(temp3_buffer,&dim3a);
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2a = nd6 * (end_occ_array[myrow] - begin_occ_array[myrow]);
    AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
    ResetDoubleArray(temp2a_buffer,&dim2a);
    fread(temp2a_buffer,sizeof(double),dim2a,integrals_temp2a);
    for (i = 0; i < (end_occ_array[myrow] - begin_occ_array[myrow]); i++) {
      for (a1 = 0; a1 < nd6; a1++) {
        for (a2 = 0; a2 < dim1ax; a2++) {
          temp3_buffer[i * dim1ax + a2] += temp2a_buffer[i * nd6 + a1] * V_inv->a[bfposa1 + a1][a2];
         }
        }
       }
      DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
     } // close loop on j1
      fclose(integrals_temp2a);

  time3 = MPI_Wtime();
  AllocateIntMatrix(&array_sizes,&job->numtasks,&job->numtasks,job);
  count_ij_ab_array_size(array_sizes,&sendtotal,sendcounts,senddispls,&recvtotal,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);
  AllocateDoubleArray(&senddata,&sendtotal,job);
  AllocateIntArray(&sendindex1,&sendtotal,job);
  AllocateIntArray(&sendindex2,&sendtotal,job);

  time3 = MPI_Wtime();
  senddata_ij_ab(temp3_buffer,senddata,&sendtotal,sendindex1,sendindex2,senddispls,ictxt,nbsize,fermi,atoms_ax,job,file);

  time3 = MPI_Wtime();
  recvdata_buffer(Hamiltonian,array_sizes,&mpA,&sendtotal,senddata,sendindex1,sendindex2,sendcounts,senddispls,&recvtotal,\
  recvdata,recvindex1,recvindex2,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);
  DestroyIntMatrix(&array_sizes,job);
  DestroyDoubleArray(&temp3_buffer,&dim3a,job);

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("hamiltonian_ij_ab                          %10.2f\n",time2);

}

void hamiltonian_ia_jb1(double *Hamiltonian, double factor, DoubleMatrix *V_inv, FERMI *fermi, int *ictxt, int *nbsize, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int il, jl, kl;
int j1, nd6, bfposa1, dim2, count;
int dim2a, dim2b;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int sendtotal, sendcounts[job->numtasks], senddispls[job->numtasks];
int recvtotal, recvcounts[job->numtasks], recvdispls[job->numtasks];
int *sendindex1, *recvindex1;
int *sendindex2, *recvindex2;
double time1, time2;
double *temp2a_buffer, *temp2b_buffer, *temp3_buffer;
double *send_buffer[job->numtasks];
double *recv_buffer[job->numtasks];
double *senddata, *recvdata;
char xx[4], yz[24] = "integrals_ov1.", zz[24] = "integrals_ov2.";
FILE *integrals_ov1, *integrals_ov2;
IntMatrix *array_sizes;

  sprintf(xx, "%d", job->taskid);
  strcat(yz,xx);
  strcat(zz,xx);

  time1 = MPI_Wtime();

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  dim2b = nqA * dim1ax;
  AllocateDoubleArray(&temp3_buffer,&dim2b,job);
  ResetDoubleArray(temp3_buffer,&dim2b);
  integrals_ov2 = fopen(zz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2b = nqA * nd6;
    AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
    ResetDoubleArray(temp2b_buffer,&dim2b);
    fread(temp2b_buffer,sizeof(double),dim2b,integrals_ov2);
    for (il = 0; il < dim1ax; il++) {
      for (jl = 0; jl < nd6; jl++) {
        for (kl = 0; kl < nqA; kl++) {
          temp3_buffer[il * nqA + kl] += factor * V_inv->a[il][bfposa1 + jl] * temp2b_buffer[jl * nqA + kl];
         }
        }
       }
      DestroyDoubleArray(&temp2b_buffer,&dim2b,job);
     }
      fclose(integrals_ov2);

  AllocateIntMatrix(&array_sizes,&job->numtasks,&job->numtasks,job);
  count_ia_jb_array_size(array_sizes,&sendtotal,sendcounts,senddispls,&recvtotal,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);
  AllocateDoubleArray(&senddata,&sendtotal,job);
  AllocateIntArray(&sendindex1,&sendtotal,job);
  AllocateIntArray(&sendindex2,&sendtotal,job);

  senddata_ia_jb(temp3_buffer,senddata,&sendtotal,sendindex1,sendindex2,senddispls,ictxt,nbsize,fermi,atoms_ax,job,file);
  //time4 = MPI_Wtime() - time3; printf("senddata %10.4f\n",time4);

  //time3 = MPI_Wtime();
  recvdata_buffer(Hamiltonian,array_sizes,&mpA,&sendtotal,senddata,sendindex1,sendindex2,sendcounts,senddispls,&recvtotal,\
  recvdata,recvindex1,recvindex2,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);

  DestroyIntMatrix(&array_sizes,job);
  DestroyDoubleArray(&temp3_buffer,&dim2b,job);

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("hamiltonian_ia_jb1                         %10.2f\n",time2);

}

void hamiltonian_ib_ja(double *Hamiltonian, DoubleMatrix *V_inv, FERMI *fermi, int *ictxt, int *nbsize, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int il, jl, kl;
int j1, nd6, bfposa1, dim2, count;
int dim2a, dim2b;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int sendtotal, sendcounts[job->numtasks], senddispls[job->numtasks];
int recvtotal, recvcounts[job->numtasks], recvdispls[job->numtasks];
int *sendindex1, *recvindex1;
int *sendindex2, *recvindex2;
double time1, time2;
double *temp2a_buffer, *temp2b_buffer, *temp3_buffer;
double *send_buffer[job->numtasks];
double *recv_buffer[job->numtasks];
double *senddata, *recvdata;
char xx[4], yz[24] = "integrals_ov1.", zz[24] = "integrals_ov2.";
FILE *integrals_ov1, *integrals_ov2;
IntMatrix *array_sizes;

  sprintf(xx, "%d", job->taskid);
  strcat(yz,xx);
  strcat(zz,xx);

  time1 = MPI_Wtime();

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  dim2b = nqA * dim1ax;
  AllocateDoubleArray(&temp3_buffer,&dim2b,job);
  ResetDoubleArray(temp3_buffer,&dim2b);
  integrals_ov2 = fopen(zz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2b = nqA * nd6;
    AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
    ResetDoubleArray(temp2b_buffer,&dim2b);
    fread(temp2b_buffer,sizeof(double),dim2b,integrals_ov2);
    for (il = 0; il < dim1ax; il++) {
      for (jl = 0; jl < nd6; jl++) {
        for (kl = 0; kl < nqA; kl++) {
          temp3_buffer[il * nqA + kl] += V_inv->a[il][bfposa1 + jl] * temp2b_buffer[jl * nqA + kl];
         }
        }
       }
      DestroyDoubleArray(&temp2b_buffer,&dim2b,job);
     }
      fclose(integrals_ov2);

  AllocateIntMatrix(&array_sizes,&job->numtasks,&job->numtasks,job);
  count_ib_ja_array_size(array_sizes,&sendtotal,sendcounts,senddispls,&recvtotal,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);
  AllocateDoubleArray(&senddata,&sendtotal,job);
  AllocateIntArray(&sendindex1,&sendtotal,job);
  AllocateIntArray(&sendindex2,&sendtotal,job);

  senddata_ib_ja(temp3_buffer,senddata,&sendtotal,sendindex1,sendindex2,senddispls,ictxt,nbsize,fermi,atoms_ax,job,file);
  //time4 = MPI_Wtime() - time3; printf("senddata %10.4f\n",time4);

  //time3 = MPI_Wtime();
  recvdata_buffer(Hamiltonian,array_sizes,&mpA,&sendtotal,senddata,sendindex1,sendindex2,sendcounts,senddispls,&recvtotal,\
  recvdata,recvindex1,recvindex2,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);

  DestroyIntMatrix(&array_sizes,job);
  DestroyDoubleArray(&temp3_buffer,&dim2b,job);

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("hamiltonian_ib_ja                          %10.2f\n",time2);

}

void hamiltonian_screened_ij_ab(double *Hamiltonian, DoubleMatrix *V_inv, FERMI *fermi, int *ictxt, int *nbsize, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int a1, a2, il, jl, kl;
int i, j, j1, k, l, m, n, r, s, nd6, bfposa1, bfposa2;
int dim2, dim2a, dim2b, dim3, dim4, dim6;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int begin_occ[job->numtasks], end_occ[job->numtasks];
double time1, time2;
double time3, time4;
double *temp4, *temp2a_buffer;
char xx[4], yz[24] = "integrals_self_energy.";
char xy[24] = "integrals_bc_oo.", xz[24] = "integrals_bc_vv.";
FILE *integrals_occvir, *integrals_occ, *integrals_vir;

  // ******************************************************************************************
  // * Set up Cblacs parameters                                                               *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  sprintf(xx, "%d", job->taskid);
  strcat(xy,xx);
  strcat(xz,xx);
  strcat(yz,xx);

  mpi_begin_end(begin_occ, end_occ, nocc, job->numtasks, job, file);

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Read in SCF or SCF_GW eigenvalues from disk                                            *
  // ******************************************************************************************

  FILE *evals;
  double *eigenvalues, *evalues_temp; // read these in
  char zz4[24] = "GW_evalues";
  AllocateDoubleArray(&eigenvalues,&nbands,job);
  ResetDoubleArray(eigenvalues,&nbands);
  if (job->taskid == 0) {
  evals = fopen(zz4, "rb");
  fseek(evals,(fermi->bands[0] - 1) * sizeof(double),SEEK_SET);
  fread(eigenvalues,sizeof(double),job->spin_dim * nbands,evals);
  fclose(evals);
 }
  MPI_Bcast(eigenvalues,job->spin_dim * nbands,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // ******************************************************************************************
  // * Read CASIDA matrix eigenvalues                                                         *
  // ******************************************************************************************
  
  double *cas_eigenvalues, *sqrt_cas_eigenvalues;
  FILE *cas_evals;
  char zz6[24] = "cas_evalues";

  AllocateDoubleArray(&sqrt_cas_eigenvalues,&ntransitions,job);
  if (job->taskid == 0) {
  AllocateDoubleArray(&cas_eigenvalues,&ntransitions,job);
  cas_evals = fopen(zz6, "rb");
  fread(&cas_eigenvalues[0], sizeof(double), ntransitions, cas_evals);
  if (job->verbosity > 1) for (i = 0; i < ntransitions; i++) fprintf(file.out,"%3d %10.4f\n",i, cas_eigenvalues[i]);
  for (i = 0; i < ntransitions; i++) sqrt_cas_eigenvalues[i] = sqrt(cas_eigenvalues[i]);
  fclose(cas_evals);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
 }
  MPI_Bcast(sqrt_cas_eigenvalues,ntransitions,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->rpa_lim == 0 || job->rpa_lim > ntransitions) job->rpa_lim = ntransitions;
  if (job->taskid == 0) printf("RPA vector limit %3d\n",job->rpa_lim);
  dim4 = job->rpa_lim * dim1ax;
  //dim4 = ntransitions * dim1ax;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  generate_temp4(temp4,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);

 //for(i = 0; i < dim4; i++) fprintf(file.out,"TEMP4 %3d %10.4f\n",i,temp4[i]);

  // ******************************************************************************************
  // * Set up M_buffer arrays                                                                 *
  // ******************************************************************************************

int I1, I2, cblacs_taskid, row_proc, col_proc;
int *integral_indices[job->numtasks], *integral_indices_buffer[job->numtasks];
int count1[job->numtasks];
int begin_vir_row[job->numtasks], end_vir_row[job->numtasks];
int begin_vir_array[job->numtasks], end_vir_array[job->numtasks];
int begin_occ_row[job->numtasks], end_occ_row[job->numtasks];
int begin_occ_array[job->numtasks], end_occ_array[job->numtasks];
int occ_pair, vir_pair;
int sendtotal, sendcounts[job->numtasks], senddispls[job->numtasks];
int recvtotal, recvcounts[job->numtasks], recvdispls[job->numtasks];
int *sendindex1, *recvindex1;
int *sendindex2, *recvindex2;
double *send_buffer[job->numtasks];
double *recv_buffer[job->numtasks];
double *senddata, *recvdata;
IntMatrix *array_sizes_buffer, *array_sizes;
DoubleMatrix *M_buffer1, *M_buffer2;

  mpi_begin_end_lower_triangle(begin_occ_row, end_occ_row, begin_occ_array, end_occ_array, nocc, nprow, job, file);
  mpi_begin_end_lower_triangle(begin_vir_row, end_vir_row, begin_vir_array, end_vir_array, nvir, npcol, job, file);

  // ******************************************************************************************
  // * Read from disk and contract V <alpha|beta><beta|occ-vir> with <alpha|occ-occ>          *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  dim6 = end_occ_array[myrow] - begin_occ_array[myrow]; 
  AllocateDoubleMatrix(&M_buffer1,&dim6,&job->rpa_lim,job);
  //AllocateDoubleMatrix(&M_buffer1,&dim6,&ntransitions,job);
  ResetDoubleMatrix(M_buffer1);
  integrals_occ = fopen(xy, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2a = nd6 * (end_occ_array[myrow] - begin_occ_array[myrow]);
    AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
    ResetDoubleArray(temp2a_buffer,&dim2a);
    fread(temp2a_buffer,sizeof(double),dim2a,integrals_occ);
    for (i = begin_occ_row[myrow]; i < end_occ_row[myrow]; i++) {
      occ_pair = (i * (i + 1)) / 2 - begin_occ_array[myrow];
      for (j = 0; j <= i; j++) {
        for (l = 0; l < job->rpa_lim; l++) {
        //for (l = 0; l < ntransitions; l++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer1->a[occ_pair + j][l] += temp2a_buffer[(occ_pair + j) * nd6 + a1] * temp4[l * dim1ax + bfposa1 + a1] / \
            sqrt_cas_eigenvalues[l];
           }
          }
         }
        }
      DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
     } // close loop on j1
      fclose(integrals_occ);
      time4 = MPI_Wtime() - time3; 
      //printf("M_buffer1 %10.4f\n",time4);

  // ******************************************************************************************
  // * Read from disk and contract V <alpha|beta><beta|occ-vir> with <alpha|vir-vir>          *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  dim6 = end_vir_array[mycol] - begin_vir_array[mycol]; 
  AllocateDoubleMatrix(&M_buffer2,&dim6,&job->rpa_lim,job);
  //AllocateDoubleMatrix(&M_buffer2,&dim6,&ntransitions,job);
  ResetDoubleMatrix(M_buffer2);
  integrals_vir = fopen(xz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2a = nd6 * (end_vir_array[mycol] - begin_vir_array[mycol]);
    AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
    ResetDoubleArray(temp2a_buffer,&dim2a);
    fread(temp2a_buffer,sizeof(double),dim2a,integrals_vir);
    for (i = begin_vir_row[mycol]; i < end_vir_row[mycol]; i++) {
      vir_pair = (i * (i + 1)) / 2 - begin_vir_array[mycol];
      for (j = 0; j <= i; j++) {
        for (l = 0; l < job->rpa_lim; l++) {
        //for (l = 0; l < ntransitions; l++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer2->a[vir_pair + j][l] += temp2a_buffer[(vir_pair + j) * nd6 + a1] * temp4[l * dim1ax + bfposa1 + a1] / \
            sqrt_cas_eigenvalues[l];
           }
          }
         }
        }
      DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
     } // close loop on j1
      fclose(integrals_vir);
      time4 = MPI_Wtime() - time3; 
      //printf("M_buffer2 %10.4f\n",time4);

  DestroyDoubleArray(&temp4,&dim4,job);

  time3 = MPI_Wtime();
  AllocateIntMatrix(&array_sizes,&job->numtasks,&job->numtasks,job);
  count_ij_ab_array_size(array_sizes,&sendtotal,sendcounts,senddispls,&recvtotal,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);
  AllocateDoubleArray(&senddata,&sendtotal,job);
  AllocateIntArray(&sendindex1,&sendtotal,job);
  AllocateIntArray(&sendindex2,&sendtotal,job);
  //time4 = MPI_Wtime() - time3; printf("countdata %10.4f\n",time4);

  time3 = MPI_Wtime();
  senddata_screened_ij_ab(M_buffer1,M_buffer2,senddata,&sendtotal,sendindex1,sendindex2,senddispls,ictxt,nbsize,fermi,job,file);
  DestroyDoubleMatrix(&M_buffer1,job);
  DestroyDoubleMatrix(&M_buffer2,job);
  //time4 = MPI_Wtime() - time3; printf("senddata %10.4f\n",time4);

  time3 = MPI_Wtime();
  recvdata_buffer(Hamiltonian,array_sizes,&mpA,&sendtotal,senddata,sendindex1,sendindex2,sendcounts,senddispls,&recvtotal,\
  recvdata,recvindex1,recvindex2,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);
  DestroyIntMatrix(&array_sizes,job);
  DestroyDoubleArray(&sqrt_cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&eigenvalues,&nbands,job);
  //time4 = MPI_Wtime() - time3; printf("recvdata %10.4f\n",time4);

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("hamiltonian_screened_ij_ab                 %10.2f\n",time2);

}

void hamiltonian_screened_ib_ja(double *Hamiltonian, DoubleMatrix *V_inv, FERMI *fermi, int *ictxt, int *nbsize, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int a1, a2, il, jl, kl;
int i, j, j1, k, l, m, n, r, s, nd6, bfposa1, bfposa2;
int dim2, dim2a, dim2b, dim3, dim4, dim6;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
double time1, time2;
double time3, time4;
double *temp4, *temp2a_buffer, *temp2b_buffer;
char xx[4], yz[24] = "integrals_ov1.", zz[24] = "integrals_ov2.";
FILE *integrals_ov1, *integrals_ov2;
DoubleMatrix *M_buffer1, *M_buffer2;

  // ******************************************************************************************
  // * Set up Cblacs parameters                                                               *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  sprintf(xx, "%d", job->taskid);
  strcat(yz,xx);
  strcat(zz,xx);

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  int global_row_index, global_col_index;
  global_row_index = 0;
  global_col_index = 0;
  for (i = 0; i < myrow; i++) global_row_index += numroc_(&ntransitions, nbsize, &i, &izero, &nprow);
  for (i = 0; i < mycol; i++) global_col_index += numroc_(&ntransitions, nbsize, &i, &izero, &npcol);
  //printf("%3d   %3d %3d  %3d %3d  %3d %3d\n",job->taskid,myrow,mycol,global_row_index,global_col_index,\
  numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow),numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol));

  // ******************************************************************************************
  // * Read in SCF or SCF_GW eigenvalues from disk                                            *
  // ******************************************************************************************

  FILE *evals;
  double *eigenvalues, *evalues_temp; // read these in
  char zz4[24] = "GW_evalues";
  AllocateDoubleArray(&eigenvalues,&nbands,job);
  ResetDoubleArray(eigenvalues,&nbands);
  if (job->taskid == 0) {
  evals = fopen(zz4, "rb");
  fseek(evals,(fermi->bands[0] - 1) * sizeof(double),SEEK_SET);
  fread(eigenvalues,sizeof(double),job->spin_dim * nbands,evals);
  fclose(evals);
 }
  MPI_Bcast(eigenvalues,job->spin_dim * nbands,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // ******************************************************************************************
  // * Read CASIDA matrix eigenvalues                                                         *
  // ******************************************************************************************
  
  double *cas_eigenvalues, *sqrt_cas_eigenvalues;
  FILE *cas_evals;
  char zz6[24] = "cas_evalues";

  AllocateDoubleArray(&sqrt_cas_eigenvalues,&ntransitions,job);
  if (job->taskid == 0) {
  AllocateDoubleArray(&cas_eigenvalues,&ntransitions,job);
  cas_evals = fopen(zz6, "rb");
  fread(&cas_eigenvalues[0], sizeof(double), ntransitions, cas_evals);
  if (job->verbosity > 1) for (i = 0; i < ntransitions; i++) fprintf(file.out,"%3d %10.4f\n",i, cas_eigenvalues[i]);
  for (i = 0; i < ntransitions; i++) sqrt_cas_eigenvalues[i] = sqrt(cas_eigenvalues[i]);
  fclose(cas_evals);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
 }
  MPI_Bcast(sqrt_cas_eigenvalues,ntransitions,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->rpa_lim == 0 || job->rpa_lim > ntransitions) job->rpa_lim = ntransitions;
  if (job->taskid == 0) printf("RPA vector limit %3d\n",job->rpa_lim);
  dim4 = job->rpa_lim * dim1ax;
  //dim4 = ntransitions * dim1ax;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  generate_temp4(temp4,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);

  // ******************************************************************************************
  // * Set up M_buffer arrays                                                                 *
  // ******************************************************************************************

  // ******************************************************************************************
  // * Read from disk and contract V <alpha|beta><beta|occ-vir> with <alpha|occ-vir>          *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  AllocateDoubleMatrix(&M_buffer1,&mpA,&job->rpa_lim,job);
  //AllocateDoubleMatrix(&M_buffer1,&mpA,&ntransitions,job);
  ResetDoubleMatrix(M_buffer1);
  integrals_ov1 = fopen(yz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2a = nd6 * mpA;
    AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
    ResetDoubleArray(temp2a_buffer,&dim2a);
    fread(temp2a_buffer,sizeof(double),dim2a,integrals_ov1);
    for (il = 0; il < mpA; il++) {
      for (l = 0; l < job->rpa_lim; l++) {
      //for (l = 0; l < ntransitions; l++) {
        for (a1 = 0; a1 < nd6; a1++) {
          M_buffer1->a[il][l] += temp2a_buffer[il + a1 * mpA] * temp4[l * dim1ax + bfposa1 + a1] / \
          sqrt_cas_eigenvalues[l];
         }
        }
       }
       DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
      } // close loop on j1
       fclose(integrals_ov1);
       time4 = MPI_Wtime() - time3; 
       //printf("M_buffer1 %10.4f\n",time4);

  // ******************************************************************************************
  // * Read from disk and contract V <alpha|beta><beta|occ-vir> with <alpha|occ-vir>          *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  AllocateDoubleMatrix(&M_buffer2,&nqA,&job->rpa_lim,job);
  //AllocateDoubleMatrix(&M_buffer2,&nqA,&ntransitions,job);
  ResetDoubleMatrix(M_buffer2);
  integrals_ov2 = fopen(zz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2b = nd6 * nqA;
    AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
    ResetDoubleArray(temp2b_buffer,&dim2b);
    fread(temp2b_buffer,sizeof(double),dim2b,integrals_ov2);
    for (il = 0; il < nqA; il++) {
      for (a1 = 0; a1 < nd6; a1++) {
        for (l = 0; l < job->rpa_lim; l++) {
        //for (l = 0; l < ntransitions; l++) {
          M_buffer2->a[il][l] += temp2b_buffer[il + a1 * nqA] * temp4[l * dim1ax + bfposa1 + a1] / \
          sqrt_cas_eigenvalues[l];
         }
        }
       }
       DestroyDoubleArray(&temp2b_buffer,&dim2b,job);
      } // close loop on j1
       fclose(integrals_ov2);
       time4 = MPI_Wtime() - time3; 
       //printf("M_buffer2 %10.4f\n",time4);

  DestroyDoubleArray(&temp4,&dim4,job);

int sendtotal, sendcounts[job->numtasks], senddispls[job->numtasks];
int recvtotal, recvcounts[job->numtasks], recvdispls[job->numtasks];
int *sendindex1, *recvindex1;
int *sendindex2, *recvindex2;
double *send_buffer[job->numtasks];
double *recv_buffer[job->numtasks];
double *senddata, *recvdata;
IntMatrix *array_sizes_buffer, *array_sizes;

  time3 = MPI_Wtime();
  AllocateIntMatrix(&array_sizes,&job->numtasks,&job->numtasks,job);
  count_ib_ja_array_size(array_sizes,&sendtotal,sendcounts,senddispls,&recvtotal,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);
  AllocateDoubleArray(&senddata,&sendtotal,job);
  AllocateIntArray(&sendindex1,&sendtotal,job);
  AllocateIntArray(&sendindex2,&sendtotal,job);
  //time4 = MPI_Wtime() - time3; printf("countdata %10.4f\n",time4);

  time3 = MPI_Wtime();
  senddata_screened_ib_ja(M_buffer1,M_buffer2,senddata,&sendtotal,sendindex1,sendindex2,senddispls,ictxt,nbsize,fermi,job,file);
  DestroyDoubleMatrix(&M_buffer1,job);
  DestroyDoubleMatrix(&M_buffer2,job);
  //time4 = MPI_Wtime() - time3; printf("senddata %10.4f\n",time4);

  time3 = MPI_Wtime();
  recvdata_buffer(Hamiltonian,array_sizes,&mpA,&sendtotal,senddata,sendindex1,sendindex2,sendcounts,senddispls,&recvtotal,\
  recvdata,recvindex1,recvindex2,recvcounts,recvdispls,ictxt,nbsize,fermi,job,file);
  DestroyIntMatrix(&array_sizes,job);
  DestroyDoubleArray(&sqrt_cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&eigenvalues,&nbands,job);
  //time4 = MPI_Wtime() - time3; printf("recvdata %10.4f\n",time4);

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("hamiltonian_screened_ib_ja                 %10.2f\n",time2);

}

void generate_temp4(double *temp4, DoubleMatrix *V_inv, FERMI *fermi, int *ictxt, int *nbsize, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int i, j1, k;
int a1, a2, nd6;
int bfposa1, bfposa2;
int dim2, dim2b, dim3;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int begin_occ[job->numtasks], end_occ[job->numtasks];
double time1, time2;
double time3, time4;
double *temp2b_buffer, *temp3_buffer;
char xx[4], yz[24] = "integrals_self_energy.";
FILE *integrals_occvir;

  // ******************************************************************************************
  // * Set up Cblacs parameters                                                               *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  sprintf(xx, "%d", job->taskid);
  strcat(yz,xx);

  mpi_begin_end(begin_occ, end_occ, nocc, job->numtasks, job, file);

//  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
//  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
//  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Read <alpha|occ-vir> integrals and multiply in V <alpha|beta>                          *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  dim2 = (end_occ[job->taskid] - begin_occ[job->taskid]) * nvir;
  dim3 = dim2 * dim1ax;
  AllocateDoubleArray(&temp3_buffer,&dim3,job);
  ResetDoubleArray(temp3_buffer,&dim3);
  integrals_occvir = fopen(yz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa2 = atoms_ax->bfnposn_sh[j1];
    dim2b = dim2 * nd6;
    AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
    ResetDoubleArray(temp2b_buffer,&dim2b);
    fread(temp2b_buffer,sizeof(double),dim2b,integrals_occvir);
    for (a1 = 0; a1 < dim1ax; a1++) {
      for (k = 0; k < dim2; k++) {
        for (a2 = 0; a2 < nd6; a2++) {
          temp3_buffer[a1 * dim2 + k] += V_inv->a[a1][bfposa2 + a2] * temp2b_buffer[k * nd6 + a2];
         }
        }
       }
      DestroyDoubleArray(&temp2b_buffer,&dim2b,job);
     }
  fclose(integrals_occvir);
  time4 = MPI_Wtime() - time3; 
  //printf("temp3 %10.4f\n",time4);

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with CASIDA eigenvectors                         *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  if (job->rpa_lim == 0 || job->rpa_lim > ntransitions) job->rpa_lim = ntransitions;
  if (job->taskid == 0) printf("RPA vector limit %3d\n",job->rpa_lim);
  int dim4 = job->rpa_lim * dim1ax;
  //int dim4 = ntransitions * dim1ax;
  double *temp4_buffer;
  AllocateDoubleArray(&temp4_buffer,&dim4,job);
  ResetDoubleArray(temp4_buffer,&dim4);

  DoubleMatrix *cas_eigenvectors;
  char xc[22] = "/cas_eigenvectors_mpi";
  char bufcas[120];
  MPI_File fh;
  strcpy(bufcas,file.scf_eigvec);
  strcat(bufcas,xc);
  int ione = 1;

  int itr = job->rpa_lim / nvir;
  int rem = job->rpa_lim - itr * nvir;
  AllocateDoubleMatrix(&cas_eigenvectors,&nvir,&ntransitions,job);
  MPI_File_open(MPI_COMM_WORLD,bufcas,MPI_MODE_RDWR,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  int l1, l2, k1 = begin_occ[job->taskid] * nvir;
  for (l1 = 0; l1 < itr; l1++) {
  //for (l1 = 0; l1 < nocc; l1++) {
    MPI_File_read(fh, &cas_eigenvectors->a[0][0], nvir * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
//for(i = 0; i < nvir; i++) {\
//for(int j = 0; j < nocc; j++){ \
//for(k = 0; k < nvir; k++) {\
//fprintf(file.out,"CAS %3d %3d %3d %10.4f \n", i,j,k, cas_eigenvectors->a[i][j * nvir + k]); }}}

    for (l2 = 0; l2 < nvir; l2++) {
      for (a1 = 0; a1 < dim1ax; a1++) {
        for (k = 0; k < (end_occ[job->taskid] - begin_occ[job->taskid]) * nvir; k++) {
          temp4_buffer[((l1 * nvir) + l2) * dim1ax + a1] += temp3_buffer[a1 * dim2 + k] * cas_eigenvectors->a[l2][k1 + k];
          //temp4_buffer[l1  * dim1ax + a1] += temp3_buffer[a1 * dim2 + k] * cas_eigenvectors->a[l1][k1 + k];
//fprintf(file.out,"%3d %3d %3d %3d %10.4f %10.4f %10.4f\n", l1,l2,k,a1,\
          temp4_buffer[((l1 * nvir) + l2) * dim1ax + a1], temp3_buffer[a1 * dim2 + k], cas_eigenvectors->a[l2][k1 + k]);
         }
        }
       }
      }

    MPI_File_read(fh, &cas_eigenvectors->a[0][0], rem * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
    for (l2 = 0; l2 < rem; l2++) {
      for (a1 = 0; a1 < dim1ax; a1++) {
        for (k = 0; k < (end_occ[job->taskid] - begin_occ[job->taskid]) * nvir; k++) {
          temp4_buffer[((itr * nvir) + l2) * dim1ax + a1] += temp3_buffer[a1 * dim2 + k] * cas_eigenvectors->a[l2][k1 + k];
         }
        }
       }

  MPI_File_close(&fh);
  DestroyDoubleMatrix(&cas_eigenvectors,job);
  DestroyDoubleArray(&temp3_buffer,&dim2b,job);
  ResetDoubleArray(temp4,&dim4);
  MPI_Allreduce(temp4_buffer,temp4,dim4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  DestroyDoubleArray(&temp4_buffer,&dim4,job);
  time4 = MPI_Wtime() - time3; 

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("generate_temp4                             %10.2f\n",time2);

}
void generate_temp4_in_core(double *temp4, double *integral_buffer2, FERMI *fermi, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

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
int num_eigenvectors;
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
  // * Contract V <alpha|beta><beta|occ-vir> with CASIDA eigenvectors                         *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  time3 = MPI_Wtime();

  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  num_eigenvectors = job->bse_lim;
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;

  AllocateIntArray(&offset_j1,&dim1,job);
  mpi_begin_end(begin_j, end_j, dim1, num_proc, job, file);
  for (i = num_proc; i < job->numtasks; i++) begin_j[i] = 0;
  for (i = num_proc; i < job->numtasks; i++) end_j[i] = 0;
  for (i = 0; i < dim1; i++) offset_j1[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j1[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1]; }}
  if (end_j[job->taskid] == begin_j[job->taskid]) { itr = 0; rem = 0; }

  //printf("task %3d %3d %3d %3d %3d\n",job->taskid,nvir, ntransitions,itr,rem);
  MPI_File_open(MPI_COMM_WORLD,bufcas,MPI_MODE_RDWR,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;

////CHP
//int kk = 0;
////CHP
  for (l1 = 0; l1 < itr; l1++) {
    //fprintf(file.out,"generate_temp4_in_core %3d allocating nvir %7d * ntransitions %7d = %7d\n",l1, \
    nvir,ntransitions, nvir * ntransitions); fflush(file.out);
    AllocateDoubleMatrix(&cas_eigenvectors,&nvir,&ntransitions,job);
    ResetDoubleMatrix(cas_eigenvectors);
    //MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
    MPI_File_read(fh, &cas_eigenvectors->a[0][0], nvir * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
    for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
      nd6 = atoms_ax->bfnnumb_sh[j1];
      offset1 = job->bse_lim * offset_j1[j1];
      offset2 = offset_j1[j1] * nbands * nbands;
      for (l2 = 0; l2 < nvir; l2++) {
        for (m = 0; m < nocc; m++) {
          for (n = 0; n < nvir; n++) {
            for (a1 = 0; a1 < nd6; a1++) {
              temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1] += \
              integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1] * cas_eigenvectors->a[l2][m * nvir + n];
////CHP
 //           if (kk < 1000 && (temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1] != temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1] || \
 //             integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1] != integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1] || \
 //              cas_eigenvectors->a[l2][m * nvir + n] != cas_eigenvectors->a[l2][m * nvir + n])) {
 //              fprintf(file.out,"main %3d %3d %3d %3d %3d %3d %f %f %f\n",l1,j1,l2,m,n,a1,temp4[((l1 * nvir) + l2) * nd6 + offset1 + a1],\
 //              integral_buffer2[offset2 + m * nbands * nd6 + (nocc + n) * nd6 + a1],cas_eigenvectors->a[l2][m * nvir + n]);
 //              kk++;
 //             }
//////CHP
             }
            }
           }
          }
         }
        DestroyDoubleMatrix(&cas_eigenvectors,job);
       }

  if (rem > 0) {
  //fprintf(file.out,"generate_temp4_in_core  allocating nvir %7d * ntransitions %7d = %7d\n", \
  nvir,ntransitions, nvir * ntransitions); fflush(file.out);
  AllocateDoubleMatrix(&cas_eigenvectors,&nvir,&ntransitions,job);
  ResetDoubleMatrix(cas_eigenvectors);
  MPI_File_read(fh, &cas_eigenvectors->a[0][0], rem * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = offset_j1[j1] * job->bse_lim;
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

      //for (a1 = 0; a1 < nocc * nvir * 6; a1++) fprintf(file.out,"TEMP4a %3d %10.4lf\n", a1, temp4[a1]);


  MPI_File_close(&fh);
  DestroyIntArray(&offset_j1,&dim1,job);
  time4 = MPI_Wtime() - time3; 
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("generate_temp4_in_core                     %10.2f\n",time2);

}

*/
/*
void self_energy(DoubleMatrix *Sigma, double *dSigma_dE, DoubleMatrix *V_inv, FERMI *fermi, int *ictxt, int *nbsize, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int a1, a2, il, jl, kl;
int i, j, j1, k, l, m, n, r, s, nd6, bfposa1, bfposa2;
int dim2, dim2a, dim2b, dim3;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int global_row_index, global_col_index;
int begin_occ[job->numtasks], end_occ[job->numtasks];
double sigma_factor, denom;
double time1, time2;
double *temp2a_buffer, *temp2b_buffer, *temp3_buffer;
char xx[4], yz[24] = "integrals_self_energy.";
//char xx[4], yz[24] = "integrals_occvir.", zz[24] = "integrals_rpa2.";
//char xy[24] = "integrals_2a.", xz[24] = "integrals_2b.";
char xy[24] = "integrals_oo.", xz[24] = "integrals_vv.";
FILE *integrals_rpa1, *integrals_rpa2, *integrals_occvir, *integrals_occ, *integrals_vir;

  // ******************************************************************************************
  // * Set up Cblacs parameters                                                               *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  sprintf(xx, "%d", job->taskid);
  strcat(xy,xx);
  strcat(xz,xx);
  strcat(yz,xx);
  //strcat(zz,xx);

  mpi_begin_end(begin_occ, end_occ, nocc, job->numtasks, job, file);

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
  global_row_index = 0;
  global_col_index = 0;
  for (i = 0; i < myrow; i++) global_row_index += numroc_(&ntransitions, nbsize, &i, &izero, &nprow);
  for (i = 0; i < mycol; i++) global_col_index += numroc_(&ntransitions, nbsize, &i, &izero, &npcol);

  // ******************************************************************************************
  // * Read in SCF or SCF_GW eigenvalues from disk                                            *
  // ******************************************************************************************

  FILE *evals, *GW_evals;
  double *eigenvalues, *evalues_temp; // read these in
  char zz4[24] = "scf_evalues";
  //char zz4[24] = "GW_evalues";
  AllocateDoubleArray(&eigenvalues,&nbands,job);
  ResetDoubleArray(eigenvalues,&nbands);
  if (job->taskid == 0) {
  evals = fopen(zz4, "rb");
  fseek(evals,(fermi->bands[0] - 1) * sizeof(double),SEEK_SET);
  fread(eigenvalues,sizeof(double),job->spin_dim * nbands,evals);
  fclose(evals);
 }
  MPI_Bcast(eigenvalues,job->spin_dim * nbands,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //for (m = 0; m < nbands; m++) fprintf(file.out,"GW1 evals %3d %10.4lf\n",m,eigenvalues[m]*au_to_eV);

  // ******************************************************************************************
  // * Read CASIDA eigenvectors and eigenvalues                                               *
  // ******************************************************************************************
  
  double *cas_eigenvalues;
  FILE *cas_evals;
  char zz6[24] = "cas_evalues";

  AllocateDoubleArray(&cas_eigenvalues,&ntransitions,job);

  if (job->taskid == 0) {
  cas_evals = fopen(zz6, "rb");
  fread(&cas_eigenvalues[0], sizeof(double), ntransitions, cas_evals);
  if (job->verbosity > 1) for (i = 0; i < ntransitions; i++) fprintf(file.out,"%3d %10.4f\n",i, cas_eigenvalues[i]);
  fclose(cas_evals);
 }
  MPI_Bcast(cas_eigenvalues,ntransitions,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // ******************************************************************************************
  // * Read <alpha|occ-vir> integrals and multiply in V <alpha|beta>                          *
  // ******************************************************************************************
  
  dim2 = (end_occ[job->taskid] - begin_occ[job->taskid]) * nvir;
  dim3 = dim2 * dim1ax;
  AllocateDoubleArray(&temp3_buffer,&dim3,job);
  ResetDoubleArray(temp3_buffer,&dim3);
  integrals_occvir = fopen(yz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa2 = atoms_ax->bfnposn_sh[j1];
    dim2b = dim2 * nd6;
    AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
    ResetDoubleArray(temp2b_buffer,&dim2b);
    fread(temp2b_buffer,sizeof(double),dim2b,integrals_occvir);
    for (a1 = 0; a1 < dim1ax; a1++) {
      for (k = 0; k < dim2; k++) {
        for (a2 = 0; a2 < nd6; a2++) {
          temp3_buffer[a1 * dim2 + k] += V_inv->a[a1][bfposa2 + a2] * temp2b_buffer[k * nd6 + a2];
          //fprintf(file.out,"%3d %3d %3d %3d %10.4f %10.4f %10.4f\n",\
          a1*dim2+k,a1,bfposa2+a2,k,temp3_buffer[a1 * dim2 + k],V_inv->a[a1][bfposa2 + a2],temp2b_buffer[k * nd6 + a2]);
         }
        }
       }
      DestroyDoubleArray(&temp2b_buffer,&dim2b,job);
     }
  fclose(integrals_occvir);

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with CASIDA eigenvectors                         *
  // ******************************************************************************************
  
  if (job->rpa_lim == 0 || job->rpa_lim > ntransitions) job->rpa_lim = ntransitions;
  if (job->taskid == 0) printf("RPA vector limit %3d\n",job->rpa_lim);
  int dim4 = job->rpa_lim * dim1ax;
  //int dim4 = ntransitions * dim1ax;
  double *temp4_buffer, *temp4;
  AllocateDoubleArray(&temp4_buffer,&dim4,job);
  ResetDoubleArray(temp4_buffer,&dim4);

  DoubleMatrix *cas_eigenvectors;
  char xc[22] = "/cas_eigenvectors_mpi";
  char bufcas[120];
  MPI_File fh;
  strcpy(bufcas,file.scf_eigvec);
  strcat(bufcas,xc);
  int ione = 1;

  int itr = job->rpa_lim / nvir;
  int rem = job->rpa_lim - itr * nvir;
  AllocateDoubleMatrix(&cas_eigenvectors,&nvir,&ntransitions,job);
  MPI_File_open(MPI_COMM_WORLD,bufcas,MPI_MODE_RDWR,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  int l1, l2, k1 = begin_occ[job->taskid] * nvir;
  for (l1 = 0; l1 < itr; l1++) {
  //for (l1 = 0; l1 < nocc; l1++) {
    MPI_File_read(fh, &cas_eigenvectors->a[0][0], nvir * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
    for (l2 = 0; l2 < nvir; l2++) {
      for (a1 = 0; a1 < dim1ax; a1++) {
        for (k = 0; k < (end_occ[job->taskid] - begin_occ[job->taskid]) * nvir; k++) {
          temp4_buffer[((l1 * nvir) + l2) * dim1ax + a1] += temp3_buffer[a1 * dim2 + k] * cas_eigenvectors->a[l2][k1 + k];
          //temp4_buffer[l1  * dim1ax + a1] += temp3_buffer[a1 * dim2 + k] * cas_eigenvectors->a[l1][k1 + k];
         }
        }
       }
      }

    MPI_File_read(fh, &cas_eigenvectors->a[0][0], rem * ntransitions, MPI_DOUBLE, MPI_STATUS_IGNORE);
    for (l2 = 0; l2 < rem; l2++) {
      for (a1 = 0; a1 < dim1ax; a1++) {
        for (k = 0; k < (end_occ[job->taskid] - begin_occ[job->taskid]) * nvir; k++) {
          temp4_buffer[((itr * nvir) + l2) * dim1ax + a1] += temp3_buffer[a1 * dim2 + k] * cas_eigenvectors->a[l2][k1 + k];
         }
        }
       }

  MPI_File_close(&fh);
  DestroyDoubleMatrix(&cas_eigenvectors,job);
  DestroyDoubleArray(&temp3_buffer,&dim2b,job);
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  MPI_Allreduce(temp4_buffer,temp4,dim4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //for(l=0;l<job->rpa_lim;l++){for(i=0;i<dim1ax;i++)fprintf(file.out,"T4 %3d %3d %10.4f\n",l,i,temp4_buffer[l*dim1ax+i]);}
  DestroyDoubleArray(&temp4_buffer,&dim4,job);
  //for(l=0;l<job->rpa_lim;l++){for(i=0;i<dim1ax;i++)fprintf(file.out,"t4 %d %d %d %f\n",job->taskid,l,i,temp4[l * dim1ax + i]);}

  // ******************************************************************************************
  // * Set up M_buffer arrays                                                                 *
  // ******************************************************************************************
  
  double *dSigma_dE_buffer;
  DoubleMatrix *Sigma_buffer, *M_buffer;
  AllocateDoubleArray(&dSigma_dE_buffer,&nbands,job);
  AllocateDoubleMatrix(&Sigma_buffer,&nbands,&nbands,job);
  ResetDoubleMatrix(Sigma_buffer);
  ResetDoubleArray(dSigma_dE_buffer,&nbands);

  // ******************************************************************************************
  // * Read from disk and contract V <alpha|beta><beta|occ-vir> with <alpha|occ-vir>          *
  // ******************************************************************************************
 
  AllocateDoubleMatrix(&M_buffer,&job->rpa_lim,&dim2,job);
  //AllocateDoubleMatrix(&M_buffer,&ntransitions,&dim2,job);
  ResetDoubleMatrix(M_buffer);
  integrals_occvir = fopen(yz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2a = dim2 * nd6;
    AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
    ResetDoubleArray(temp2a_buffer,&dim2a);
    fread(temp2a_buffer,sizeof(double),dim2a,integrals_occvir);
    //for (l = 0; l < ntransitions; l++) {
    for (l = 0; l < job->rpa_lim; l++) {
      for (k = 0; k < dim2; k++) {
        for (a1 = 0; a1 < nd6; a1++) {
          M_buffer->a[l][k] += temp2a_buffer[k * nd6 + a1] * temp4[l * dim1ax + bfposa1 + a1];
          //fprintf(file.out,"%3d %3d %3d %3d %10.4f %10.4f %10.4f\n", \
          j1,l,k,a1,M_buffer->a[l][k],temp2a_buffer[k * nd6 + a1], temp4[l * dim1ax + bfposa1 + a1]);
         }
        }
       }
      DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
     }
  fclose(integrals_occvir);

  //for (l = 0; l < ntransitions; l++) {
  for (l = 0; l < job->rpa_lim; l++) {
    for (i = 0; i < dim2; i++) {
      m = (begin_occ[job->taskid] * nvir + i) / nvir;
      n = (begin_occ[job->taskid] * nvir + i) - m  * nvir + nocc; 
      denom = eigenvalues[m] - eigenvalues[n] - cas_eigenvalues[l];
      sigma_factor = two * M_buffer->a[l][i] * M_buffer->a[l][i] / denom;
      Sigma_buffer->a[m][m] += sigma_factor;
      dSigma_dE_buffer[m] -= sigma_factor / denom;
      //fprintf(file.out,"occvir %3d %3d %3d %3d %10.4f %10.4f    %10.4f %10.4f %10.4f %10.4f \n",\
      l,m,n,i,M_buffer->a[l][i], M_buffer->a[l][i], eigenvalues[n],- eigenvalues[m], cas_eigenvalues[l],denom);
      denom = eigenvalues[n] - eigenvalues[m] + cas_eigenvalues[l];
      sigma_factor = two * M_buffer->a[l][i] * M_buffer->a[l][i] / denom;
      Sigma_buffer->a[n][n] += sigma_factor;
      dSigma_dE_buffer[n] -= sigma_factor / denom;
      //fprintf(file.out,"occvi1 %3d %3d %3d %3d %10.4f %10.4f    %10.4f %10.4f %10.4f %10.4f \n",\
      l,m,n,i,M_buffer->a[l][i], M_buffer->a[l][i], eigenvalues[m],- eigenvalues[n], -cas_eigenvalues[l],denom);
     }
    }
   DestroyDoubleMatrix(&M_buffer,job);

  // ******************************************************************************************
  // * Read from disk and contract V <alpha|beta><beta|occ-vir> with <alpha|occ-occ>          *
  // ******************************************************************************************
  
int dim6;
int occ_pair, vir_pair;
int begin_vir_row[job->numtasks], end_vir_row[job->numtasks];
int begin_vir_array[job->numtasks], end_vir_array[job->numtasks];
int begin_occ_row[job->numtasks], end_occ_row[job->numtasks];
int begin_occ_array[job->numtasks], end_occ_array[job->numtasks];

  mpi_begin_end_lower_triangle(begin_occ_row, end_occ_row, begin_occ_array, end_occ_array, nocc, job->numtasks, job, file);
  mpi_begin_end_lower_triangle(begin_vir_row, end_vir_row, begin_vir_array, end_vir_array, nvir, job->numtasks, job, file);

  dim6 = end_occ_array[job->taskid] - begin_occ_array[job->taskid]; 
  AllocateDoubleMatrix(&M_buffer,&job->rpa_lim,&dim6,job);
  //AllocateDoubleMatrix(&M_buffer,&ntransitions,&dim6,job);
  ResetDoubleMatrix(M_buffer);
  integrals_occ = fopen(xy, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2a = nd6 * (end_occ_array[job->taskid] - begin_occ_array[job->taskid]);
    AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
    ResetDoubleArray(temp2a_buffer,&dim2a);
    fread(temp2a_buffer,sizeof(double),dim2a,integrals_occ);
    for (a1 = 0; a1 < nd6; a1++) {
      //for (l = 0; l < ntransitions; l++) {
      for (l = 0; l < job->rpa_lim; l++) {
        for (i = begin_occ_row[job->taskid]; i < end_occ_row[job->taskid]; i++) {
          occ_pair = (i * (i + 1)) / 2 - begin_occ_array[job->taskid];
          for (j = 0; j <= i; j++) {
            M_buffer->a[l][occ_pair + j] += temp2a_buffer[(occ_pair + j) * nd6 + a1] * temp4[l * dim1ax + bfposa1 + a1];
            //M_buffer->a[l][occ_pair + j] += temp2a_buffer[(occ_pair + j) * nd6 + a1] * temp4[(bfposa1 + a1) * ntransitions + l];
            //fprintf(file.out,"%3d %3d %3d %3d %3d %3d %10.4f %10.4f %10.4f\n", \
            j1,a1,l,i,j,occ_pair+j,temp2a_buffer[(occ_pair + j) * nd6 + a1],temp4[l * dim1ax + bfposa1 + a1], M_buffer->a[l][i]);
           }
          }
         }
        }
      DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
     } // close loop on j1
      fclose(integrals_occ);
  //for(l = 0; l < ntransitions; l++) {for(il = 0; il < dim6; il++) { \
  fprintf(file.out,"M_buffer %3d %3d %10.4f\n",l,il,M_buffer->a[l][il]); }}
  //for (l = 0; l < ntransitions; l++) {
  for (l = 0; l < job->rpa_lim; l++) {
    for (m = begin_occ_row[job->taskid]; m < end_occ_row[job->taskid]; m++) {
      occ_pair = (m * (m + 1)) / 2 - begin_occ_array[job->taskid];
      for (r = 0; r <= m; r++) {
        //if (job->taskid > 0) continue;
        denom = eigenvalues[m] - eigenvalues[r] + cas_eigenvalues[l];
        sigma_factor = two * M_buffer->a[l][occ_pair + r] * M_buffer->a[l][occ_pair + r] / denom;
        Sigma_buffer->a[m][m] += sigma_factor;
        dSigma_dE_buffer[m] -= sigma_factor / denom;
        //fprintf(file.out,"occocc %3d %3d %3d %3d %3d %10.4f %10.4f    %10.4f %10.4f %10.4f %10.4f\n",\
        l,m,r,occ_pair+r,m,M_buffer->a[l][occ_pair+r], M_buffer->a[l][occ_pair+r], eigenvalues[m], - eigenvalues[r], \
        cas_eigenvalues[l],denom);
        if (r == m) continue;
        denom = eigenvalues[r] - eigenvalues[m] + cas_eigenvalues[l];
        sigma_factor = two * M_buffer->a[l][occ_pair + r] * M_buffer->a[l][occ_pair + r] / denom;
        Sigma_buffer->a[r][r] += sigma_factor;
        dSigma_dE_buffer[r] -= sigma_factor / denom;
        //fprintf(file.out,"occoc1 %3d %3d %3d %3d %3d %10.4f %10.4f    %10.4f %10.4f %10.4f %10.4f\n",\
        l,m,r,occ_pair+r,r,M_buffer->a[l][occ_pair+r], M_buffer->a[l][occ_pair+r], eigenvalues[m], - eigenvalues[r], \
        cas_eigenvalues[l],denom);
       }
      }
     }
    DestroyDoubleMatrix(&M_buffer,job);

  // ******************************************************************************************
  // * Read from disk and contract V <alpha|beta><beta|occ-vir> with <alpha|vir-vir>          *
  // ******************************************************************************************
  
  dim6 = end_vir_array[job->taskid] - begin_vir_array[job->taskid]; 
  //AllocateDoubleMatrix(&M_buffer,&ntransitions,&dim6,job);
  AllocateDoubleMatrix(&M_buffer,&job->rpa_lim,&dim6,job);
  ResetDoubleMatrix(M_buffer);
  integrals_vir = fopen(xz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    bfposa1 = atoms_ax->bfnposn_sh[j1];
    dim2a = nd6 * (end_vir_array[job->taskid] - begin_vir_array[job->taskid]);
    AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
    ResetDoubleArray(temp2a_buffer,&dim2a);
    fread(temp2a_buffer,sizeof(double),dim2a,integrals_vir);
    for (a1 = 0; a1 < nd6; a1++) {
      //for (l = 0; l < ntransitions; l++) {
      for (l = 0; l < job->rpa_lim; l++) {
        for (i = begin_vir_row[job->taskid]; i < end_vir_row[job->taskid]; i++) {
          vir_pair = (i * (i + 1)) / 2 - begin_vir_array[job->taskid];
          for (j = 0; j <= i; j++) {
            M_buffer->a[l][vir_pair + j] += temp2a_buffer[(vir_pair + j) * nd6 + a1] * temp4[l * dim1ax + bfposa1 + a1];
            //fprintf(file.out,"%3d %3d %3d %3d %3d %3d %10.4f %10.4f %10.4f\n", \
            j1,a1,l,i,j,vir_pair+j,temp2a_buffer[(vir_pair+j) * nd6 + a1],temp4[(bfposa1+a1)*ntransitions+l], M_buffer->a[l][i]);
           }
          }
         }
        }
      DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
     } // close loop on j1
      fclose(integrals_vir);

  //for (l = 0; l < ntransitions; l++) {
  for (l = 0; l < job->rpa_lim; l++) {
    for (n = begin_vir_row[job->taskid]; n < end_vir_row[job->taskid]; n++) {
      vir_pair = (n * (n + 1)) / 2 - begin_vir_array[job->taskid];
      for (s = 0; s <= n; s++) {
        denom = eigenvalues[nocc + n] - eigenvalues[nocc + s] - cas_eigenvalues[l];
        sigma_factor = two * M_buffer->a[l][vir_pair + s] * M_buffer->a[l][vir_pair + s] / denom;
        Sigma_buffer->a[nocc + n][nocc + n] += sigma_factor;
        dSigma_dE_buffer[nocc + n] -= sigma_factor / denom;
        //fprintf(file.out,"virvir %3d %3d %3d %3d %3d %10.4f %10.4f    %10.4f %10.4f %10.4f %10.4f\n",\
        l,n,s,vir_pair+s,nocc+n,M_buffer->a[l][vir_pair+s],M_buffer->a[l][vir_pair+s],eigenvalues[nocc+n],-eigenvalues[nocc+s],\
        -cas_eigenvalues[l],denom);
        if (n == s) continue;
        denom = eigenvalues[nocc + s] - eigenvalues[nocc + n] - cas_eigenvalues[l];
        sigma_factor = two * M_buffer->a[l][vir_pair + s] * M_buffer->a[l][vir_pair + s] / denom;
        Sigma_buffer->a[nocc + s][nocc + s] += sigma_factor;
        dSigma_dE_buffer[nocc + s] -= sigma_factor / denom;
        //fprintf(file.out,"virvi1 %3d %3d %3d %3d %3d %10.4f %10.4f    %10.4f %10.4f %10.4f %10.4f\n",\
        l,n,s,vir_pair+s,nocc+s,M_buffer->a[l][vir_pair+s],M_buffer->a[l][vir_pair+s],eigenvalues[nocc+s],-eigenvalues[nocc+n],\
        -cas_eigenvalues[l],denom);
       }
      }
     }
    DestroyDoubleMatrix(&M_buffer,job);

  MPI_Reduce(&Sigma_buffer->a[0][0],&Sigma->a[0][0],nbands*nbands,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(dSigma_dE_buffer,dSigma_dE,nbands,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  DestroyDoubleArray(&dSigma_dE_buffer,&nbands,job);
  DestroyDoubleMatrix(&Sigma_buffer,job);
  DestroyDoubleArray(&temp4,&dim4,job);

  // ******************************************************************************************
  // * Generate Self-energy for each state required                                           *
  // ******************************************************************************************
  
  fprintf(file.out,"GW Eigenvalues: %4d RPA states used\n",job->rpa_lim);
  for (m = 0; m < nbands; m++) {
    fprintf(file.out,"%3d %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n", m + 1, eigenvalues[m] * au_to_eV, \
    (eigenvalues[m] + Sigma->a[m][m]) * au_to_eV, (eigenvalues[m] + Sigma->a[m][m] / (k_one - dSigma_dE[m])) * \
    au_to_eV, Sigma->a[m][m] * au_to_eV, Sigma->a[m][m] / (k_one - dSigma_dE[m]) * au_to_eV, k_one / (k_one - dSigma_dE[m]));
   }
  fflush(file.out);

  if (job->taskid == 0) {
  AllocateDoubleArray(&evalues_temp,&nbands,job);
  for (i = 0; i < nbands; i++) {
  evalues_temp[i] = eigenvalues[i] + Sigma->a[i][i] / (k_one - dSigma_dE[i]);
  //fprintf(file.out,"GW0 evals %3d %10.4lf\n",i,evalues_temp[i]*au_to_eV);
 }
  char zz5[24] = "GW_evalues";
  GW_evals = fopen(zz5, "wb");
  fseek(GW_evals,(fermi->bands[0] - 1) * sizeof(double),SEEK_SET);
  fwrite(evalues_temp,sizeof(double),job->spin_dim * nbands,GW_evals);
  fclose(GW_evals);
  DestroyDoubleArray(&evalues_temp,&nbands,job);
 }
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&eigenvalues,&nbands,job);

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("end Sigma                                  %10.2f\n",time2);

}

void count_ij_ab_array_size(IntMatrix *array_sizes, int *sendtotal, int *sendcounts, int *senddispls, int *recvtotal, int *recvcounts, int *recvdispls, int *ictxt, int *nbsize, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, m, n, r, s;
int I1, I2, row_proc, col_proc, cblacs_taskid;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int count1[job->numtasks];
int begin_vir_row[job->numtasks],   end_vir_row[job->numtasks];
int begin_vir_array[job->numtasks], end_vir_array[job->numtasks];
int begin_occ_row[job->numtasks],   end_occ_row[job->numtasks];
int begin_occ_array[job->numtasks], end_occ_array[job->numtasks];
double time1, time2;
IntMatrix *array_sizes_buffer;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
  mpi_begin_end_lower_triangle(begin_occ_row, end_occ_row, begin_occ_array, end_occ_array, nocc, nprow, job, file);
  mpi_begin_end_lower_triangle(begin_vir_row, end_vir_row, begin_vir_array, end_vir_array, nvir, npcol, job, file);

  time1 = MPI_Wtime();
  for (i = 0; i < job->numtasks; i++) count1[i] = 0;
  AllocateIntMatrix(&array_sizes_buffer,&job->numtasks,&job->numtasks,job);
  ResetIntMatrix(array_sizes_buffer);
  ResetIntMatrix(array_sizes);
  for (m = begin_occ_row[myrow]; m < end_occ_row[myrow]; m++) {
    for (r = 0; r <= m; r++) {
      for (n = begin_vir_row[mycol]; n < end_vir_row[mycol]; n++) {
        for (s = 0; s <= n; s++) {
          I1 = m * nvir + n;
          I2 = r * nvir + s;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          (count1[cblacs_taskid])++;
          (array_sizes_buffer->a[job->taskid][cblacs_taskid])++;
          if (r < m) {
          I1 = r * nvir + n;
          I2 = m * nvir + s;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          (count1[cblacs_taskid])++;
          (array_sizes_buffer->a[job->taskid][cblacs_taskid])++;
         }
          if (s < n) {
          I1 = m * nvir + s;
          I2 = r * nvir + n;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          (count1[cblacs_taskid])++;
          (array_sizes_buffer->a[job->taskid][cblacs_taskid])++;
         }
          if (r< m && s < n) {
          I1 = r * nvir + s;
          I2 = m * nvir + n;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          (count1[cblacs_taskid])++;
          (array_sizes_buffer->a[job->taskid][cblacs_taskid])++;
         }
        }
       }
      }
     }
    MPI_Allreduce(&array_sizes_buffer->a[0][0],&array_sizes->a[0][0],job->numtasks * job->numtasks,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    DestroyIntMatrix(&array_sizes_buffer,job);

    *sendtotal = 0;
    *recvtotal = 0;
    senddispls[0] = 0;
    recvdispls[0] = 0;
    for (i = 0; i < job->numtasks; i++) sendcounts[i] = count1[i];
    for (i = 0; i < job->numtasks - 1; i++) { senddispls[i+1] = sendcounts[i] + *sendtotal; *sendtotal += sendcounts[i]; }
    *sendtotal += sendcounts[job->numtasks - 1];
    for (i = 0; i < job->numtasks; i++) recvcounts[i] = array_sizes->a[i][job->taskid];
    for (i = 0; i < job->numtasks - 1; i++) { recvdispls[i+1] = recvcounts[i] + *recvtotal; *recvtotal += recvcounts[i]; }
    *recvtotal += recvcounts[job->numtasks - 1];

    time2 = MPI_Wtime() - time1;
    if (job->taskid == 0) printf("count_ij_ab_array_size                     %10.2f\n",time2);

}

void count_ia_jb_array_size(IntMatrix *array_sizes, int *sendtotal, int *sendcounts, int *senddispls, int *recvtotal, int *recvcounts, int *recvdispls, int *ictxt, int *nbsize, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, m, n, r, s;
int il, jl;
int I1, I2, I3, I4, row_proc, col_proc, cblacs_taskid;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int count1[job->numtasks];
int global_row_index, global_col_index;
double time1, time2;
IntMatrix *array_sizes_buffer;

  time1 = MPI_Wtime();
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  global_row_index = 0;
  global_col_index = 0;
  for (i = 0; i < myrow; i++) global_row_index += numroc_(&ntransitions, nbsize, &i, &izero, &nprow);
  for (i = 0; i < mycol; i++) global_col_index += numroc_(&ntransitions, nbsize, &i, &izero, &npcol);

  for (i = 0; i < job->numtasks; i++) count1[i] = 0;
  AllocateIntMatrix(&array_sizes_buffer,&job->numtasks,&job->numtasks,job);
  ResetIntMatrix(array_sizes_buffer);
  ResetIntMatrix(array_sizes);
  for (il = 0; il < mpA; il++) {
    I1 = global_row_index + il;
    m = I1 / (nbands - fermi->occupied[0]);
    //s = I1  - m  * (nbands - fermi->occupied[0]);
    n = I1  - m  * (nbands - fermi->occupied[0]);
    for (jl = 0; jl < nqA; jl++) {
      I2 = global_col_index + jl;
      r = I2 / (nbands - fermi->occupied[0]);
      s = I2  - r  * (nbands - fermi->occupied[0]);
      //n = I2  - r  * (nbands - fermi->occupied[0]);
      I3 = m * nvir + n;
      I4 = r * nvir + s;
      row_proc = (I3 / *nbsize) % nprow;
      col_proc = (I4 / *nbsize) % npcol;
      //fprintf(file.out,"il jl %3d %3d glob rc %3d %3d msrn %3d %3d %3d %3d mnrs %3d %3d %3d %3d rc %3d %3d\n",\
      il,jl,global_row_index,global_col_index,m,s,r,n,m,n,r,s,row_proc,col_proc);
      cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
      (count1[cblacs_taskid])++;
      (array_sizes_buffer->a[job->taskid][cblacs_taskid])++;
     }
    }
    MPI_Allreduce(&array_sizes_buffer->a[0][0],&array_sizes->a[0][0],job->numtasks * job->numtasks,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    DestroyIntMatrix(&array_sizes_buffer,job);

    *sendtotal = 0;
    *recvtotal = 0;
    senddispls[0] = 0;
    recvdispls[0] = 0;
    for (i = 0; i < job->numtasks; i++) sendcounts[i] = count1[i];
    for (i = 0; i < job->numtasks - 1; i++) { senddispls[i+1] = sendcounts[i] + *sendtotal; *sendtotal += sendcounts[i]; }
    *sendtotal += sendcounts[job->numtasks - 1];
    for (i = 0; i < job->numtasks; i++) recvcounts[i] = array_sizes->a[i][job->taskid];
    for (i = 0; i < job->numtasks - 1; i++) { recvdispls[i+1] = recvcounts[i] + *recvtotal; *recvtotal += recvcounts[i]; }
    *recvtotal += recvcounts[job->numtasks - 1];

    time2 = MPI_Wtime() - time1;
    if (job->taskid == 0) printf("count_ia_jb_array_size                     %10.2f\n",time2);

}

void count_ib_ja_array_size(IntMatrix *array_sizes, int *sendtotal, int *sendcounts, int *senddispls, int *recvtotal, int *recvcounts, int *recvdispls, int *ictxt, int *nbsize, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, m, n, r, s;
int il, jl;
int I1, I2, I3, I4, row_proc, col_proc, cblacs_taskid;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int count1[job->numtasks];
int global_row_index, global_col_index;
double time1, time2;
IntMatrix *array_sizes_buffer;

  time1 = MPI_Wtime();
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  global_row_index = 0;
  global_col_index = 0;
  for (i = 0; i < myrow; i++) global_row_index += numroc_(&ntransitions, nbsize, &i, &izero, &nprow);
  for (i = 0; i < mycol; i++) global_col_index += numroc_(&ntransitions, nbsize, &i, &izero, &npcol);

  for (i = 0; i < job->numtasks; i++) count1[i] = 0;
  AllocateIntMatrix(&array_sizes_buffer,&job->numtasks,&job->numtasks,job);
  ResetIntMatrix(array_sizes_buffer);
  ResetIntMatrix(array_sizes);
  for (il = 0; il < mpA; il++) {
    I1 = global_row_index + il;
    m = I1 / (nbands - fermi->occupied[0]);
    s = I1  - m  * (nbands - fermi->occupied[0]);
    for (jl = 0; jl < nqA; jl++) {
      I2 = global_col_index + jl;
      r = I2 / (nbands - fermi->occupied[0]);
      n = I2  - r  * (nbands - fermi->occupied[0]);
      I3 = m * nvir + n;
      I4 = r * nvir + s;
      row_proc = (I3 / *nbsize) % nprow;
      col_proc = (I4 / *nbsize) % npcol;
      //fprintf(file.out,"il jl %3d %3d glob rc %3d %3d msrn %3d %3d %3d %3d mnrs %3d %3d %3d %3d rc %3d %3d\n",\
      il,jl,global_row_index,global_col_index,m,s,r,n,m,n,r,s,row_proc,col_proc);
      cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
      (count1[cblacs_taskid])++;
      (array_sizes_buffer->a[job->taskid][cblacs_taskid])++;
     }
    }
    MPI_Allreduce(&array_sizes_buffer->a[0][0],&array_sizes->a[0][0],job->numtasks * job->numtasks,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    DestroyIntMatrix(&array_sizes_buffer,job);

    *sendtotal = 0;
    *recvtotal = 0;
    senddispls[0] = 0;
    recvdispls[0] = 0;
    for (i = 0; i < job->numtasks; i++) sendcounts[i] = count1[i];
    for (i = 0; i < job->numtasks - 1; i++) { senddispls[i+1] = sendcounts[i] + *sendtotal; *sendtotal += sendcounts[i]; }
    *sendtotal += sendcounts[job->numtasks - 1];
    for (i = 0; i < job->numtasks; i++) recvcounts[i] = array_sizes->a[i][job->taskid];
    for (i = 0; i < job->numtasks - 1; i++) { recvdispls[i+1] = recvcounts[i] + *recvtotal; *recvtotal += recvcounts[i]; }
    *recvtotal += recvcounts[job->numtasks - 1];

    time2 = MPI_Wtime() - time1;
    if (job->taskid == 0) printf("count_ib_ja_array_size                     %10.2f\n",time2);

}

void senddata_ij_ab(double *temp3_buffer, double *senddata, int *sendtotal, int *sendindex1, int *sendindex2, int *senddispls, int *ictxt, int *nbsize, FERMI *fermi, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int i, il, jl, m, n, r, s;
int bfposa1, a1, j1, nd6;
int dim2b;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int I1, I2, row_proc, col_proc, cblacs_taskid;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int occ_pair, vir_pair;
int count1[job->numtasks];
int begin_vir_row[job->numtasks], end_vir_row[job->numtasks];
int begin_vir_array[job->numtasks], end_vir_array[job->numtasks];
int begin_occ_row[job->numtasks], end_occ_row[job->numtasks];
int begin_occ_array[job->numtasks], end_occ_array[job->numtasks];
double time1, time2;
double *temp2b_buffer;
FILE *integrals_temp2b;
char xx[4], zz[24] = "integrals_bc_vv.";

  sprintf(xx, "%d", job->taskid);
  strcat(zz,xx);

  time1 = MPI_Wtime();
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
  mpi_begin_end_lower_triangle(begin_occ_row, end_occ_row, begin_occ_array, end_occ_array, nocc, nprow, job, file);
  mpi_begin_end_lower_triangle(begin_vir_row, end_vir_row, begin_vir_array, end_vir_array, nvir, npcol, job, file);

  ResetDoubleArray(senddata,sendtotal);
  ResetIntArray(sendindex1,sendtotal);
  ResetIntArray(sendindex2,sendtotal);

  time1 = MPI_Wtime();
  integrals_temp2b = fopen(zz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
  nd6 = atoms_ax->bfnnumb_sh[j1];
  bfposa1 = atoms_ax->bfnposn_sh[j1];
  dim2b = nd6 * (end_vir_array[mycol] - begin_vir_array[mycol]);
  AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
  ResetDoubleArray(temp2b_buffer,&dim2b);
  fread(temp2b_buffer,sizeof(double),dim2b,integrals_temp2b);
  for (i = 0; i < job->numtasks; i++) {
    count1[i] = 0;
   }
  for (m = begin_occ_row[myrow]; m < end_occ_row[myrow]; m++) {
    for (r = 0; r <= m; r++) {
      occ_pair = (m * (m + 1)) / 2 + r - begin_occ_array[myrow];
      for (n = begin_vir_row[mycol]; n < end_vir_row[mycol]; n++) {
        for (s = 0; s <= n; s++) {
          vir_pair = (n * (n + 1)) / 2 + s - begin_vir_array[mycol];
          I1 = m * nvir + n;
          I2 = r * nvir + s;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          for (a1 = 0; a1 < nd6; a1++) {
            senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] -= temp3_buffer[occ_pair * dim1ax + bfposa1 + a1] * \
            temp2b_buffer[vir_pair * nd6 + a1];
           }
          (count1[cblacs_taskid])++;
          if (r < m) {
          I1 = r * nvir + n;
          I2 = m * nvir + s;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          for (a1 = 0; a1 < nd6; a1++) {
            senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] -= temp3_buffer[occ_pair * dim1ax + bfposa1 + a1] * \
            temp2b_buffer[vir_pair * nd6 + a1];
           }
          (count1[cblacs_taskid])++;
          }
          if (s < n) {
          I1 = m * nvir + s;
          I2 = r * nvir + n;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          for (a1 = 0; a1 < nd6; a1++) {
            senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] -= temp3_buffer[occ_pair * dim1ax + bfposa1 + a1] * \
            temp2b_buffer[vir_pair * nd6 + a1];
           }
          (count1[cblacs_taskid])++;
          }
          if (r < m && s < n) {
          I1 = r * nvir + s;
          I2 = m * nvir + n;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          for (a1 = 0; a1 < nd6; a1++) {
            senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] -= temp3_buffer[occ_pair * dim1ax + bfposa1 + a1] * \
            temp2b_buffer[vir_pair * nd6 + a1];
           }
          (count1[cblacs_taskid])++;
          }
         }
        }
       }
      }
    DestroyDoubleArray(&temp2b_buffer,&dim2b,job);
     }
    fclose(integrals_temp2b);

    for (i = 0; i < job->numtasks; i++) {
      count1[i] = 0;
     }
  for (m = begin_occ_row[myrow]; m < end_occ_row[myrow]; m++) {
    for (r = 0; r <= m; r++) {
      occ_pair = (m * (m + 1)) / 2 + r - begin_occ_array[myrow];
      for (n = begin_vir_row[mycol]; n < end_vir_row[mycol]; n++) {
        for (s = 0; s <= n; s++) {
          vir_pair = (n * (n + 1)) / 2 + s - begin_vir_array[mycol];
          I1 = m * nvir + n;
          I2 = r * nvir + s;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
          jl = *nbsize * (I2 / (*nbsize * npcol)) + I2 % *nbsize;
          sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il;
          sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl;
          (count1[cblacs_taskid])++;
          if (r < m) {
          I1 = r * nvir + n;
          I2 = m * nvir + s;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
          jl = *nbsize * (I2 / (*nbsize * npcol)) + I2 % *nbsize;
          sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il;
          sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl;
          (count1[cblacs_taskid])++;
          }
          if (s < n) {
          I1 = m * nvir + s;
          I2 = r * nvir + n;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
          jl = *nbsize * (I2 / (*nbsize * npcol)) + I2 % *nbsize;
          sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il;
          sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl;
          (count1[cblacs_taskid])++;
          }
          if (r < m && s < n) {
          I1 = r * nvir + s;
          I2 = m * nvir + n;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
          jl = *nbsize * (I2 / (*nbsize * npcol)) + I2 % *nbsize;
          sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il;
          sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl;
          (count1[cblacs_taskid])++;
          }
         }
        }
       }
      }
    time2 = MPI_Wtime() - time1;
    if (job->taskid == 0) printf("senddata_ij_ab                             %10.2f\n",time2);

}

void senddata_ia_jb(double *temp3_buffer, double *senddata, int *sendtotal, int *sendindex1, int *sendindex2, int *senddispls, int *ictxt, int *nbsize, FERMI *fermi, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int i, il, jl, il1, jl1, m, n, r, s;
int bfposa1, a1, j1, nd6;
int dim2a;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int I1, I2, I3, I4, row_proc, col_proc, cblacs_taskid;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int count1[job->numtasks];
int global_row_index, global_col_index;
double time1, time2;
double *temp2a_buffer;
FILE *integrals_ov1;
char xx[4], zz[24] = "integrals_ov1.";

  sprintf(xx, "%d", job->taskid);
  strcat(zz,xx);

  time1 = MPI_Wtime();
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  global_row_index = 0;
  global_col_index = 0;
  for (i = 0; i < myrow; i++) global_row_index += numroc_(&ntransitions, nbsize, &i, &izero, &nprow);
  for (i = 0; i < mycol; i++) global_col_index += numroc_(&ntransitions, nbsize, &i, &izero, &npcol);

  ResetDoubleArray(senddata,sendtotal);
  ResetIntArray(sendindex1,sendtotal);
  ResetIntArray(sendindex2,sendtotal);

  time1 = MPI_Wtime();
  integrals_ov1 = fopen(zz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
  nd6 = atoms_ax->bfnnumb_sh[j1];
  bfposa1 = atoms_ax->bfnposn_sh[j1];
  dim2a = nd6 * mpA;
  AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
  ResetDoubleArray(temp2a_buffer,&dim2a);
  fread(temp2a_buffer,sizeof(double),dim2a,integrals_ov1);
  for (i = 0; i < job->numtasks; i++) {
    count1[i] = 0;
   }
  for (il = 0; il < mpA; il++) {
    I1 = global_row_index + il;
    m = I1 / (nbands - fermi->occupied[0]);
    n = I1  - m  * (nbands - fermi->occupied[0]);
    //s = I1  - m  * (nbands - fermi->occupied[0]);
    for (jl = 0; jl < nqA; jl++) {
      I2 = global_col_index + jl;
      r = I2 / (nbands - fermi->occupied[0]);
      s = I2  - r  * (nbands - fermi->occupied[0]);
      //n = I2  - r  * (nbands - fermi->occupied[0]);
      I3 = m * nvir + n;
      I4 = r * nvir + s;
      row_proc = (I3 / *nbsize) % nprow;
      col_proc = (I4 / *nbsize) % npcol;
      cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
      for (a1 = 0; a1 < nd6; a1++) {
        senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] += temp2a_buffer[a1 * mpA + il] * \
        temp3_buffer[(bfposa1 + a1) * nqA + jl];
       }
      (count1[cblacs_taskid])++;
      }
     }
    DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
   }
    fclose(integrals_ov1);
    for (i = 0; i < job->numtasks; i++) {
      count1[i] = 0;
     }
  for (il = 0; il < mpA; il++) {
    I1 = global_row_index + il;
    m = I1 / (nbands - fermi->occupied[0]);
    n = I1  - m  * (nbands - fermi->occupied[0]);
    //s = I1  - m  * (nbands - fermi->occupied[0]);
    for (jl = 0; jl < nqA; jl++) {
      I2 = global_col_index + jl;
      r = I2 / (nbands - fermi->occupied[0]);
      s = I2  - r  * (nbands - fermi->occupied[0]);
      //n = I2  - r  * (nbands - fermi->occupied[0]);
      I3 = m * nvir + n;
      I4 = r * nvir + s;
      row_proc = (I3 / *nbsize) % nprow;
      col_proc = (I4 / *nbsize) % npcol;
      cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
      il1 = *nbsize * (I3 / (*nbsize * nprow)) + I3 % *nbsize;
      jl1 = *nbsize * (I4 / (*nbsize * npcol)) + I4 % *nbsize;
      sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il1;
      sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl1;
      (count1[cblacs_taskid])++;
      }
     }
    time2 = MPI_Wtime() - time1;
    if (job->taskid == 0) printf("senddata_ia_jb                             %10.2f\n",time2);

}

void senddata_ib_ja(double *temp3_buffer, double *senddata, int *sendtotal, int *sendindex1, int *sendindex2, int *senddispls, int *ictxt, int *nbsize, FERMI *fermi, ATOM *atoms_ax, JOB_PARAM *job, FILES file)

{

int i, il, jl, il1, jl1, m, n, r, s;
int bfposa1, a1, j1, nd6;
int dim2a;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int I1, I2, I3, I4, row_proc, col_proc, cblacs_taskid;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int count1[job->numtasks];
int global_row_index, global_col_index;
double time1, time2;
double *temp2a_buffer;
FILE *integrals_ov1;
char xx[4], zz[24] = "integrals_ov1.";

  sprintf(xx, "%d", job->taskid);
  strcat(zz,xx);

  time1 = MPI_Wtime();
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  global_row_index = 0;
  global_col_index = 0;
  for (i = 0; i < myrow; i++) global_row_index += numroc_(&ntransitions, nbsize, &i, &izero, &nprow);
  for (i = 0; i < mycol; i++) global_col_index += numroc_(&ntransitions, nbsize, &i, &izero, &npcol);

  ResetDoubleArray(senddata,sendtotal);
  ResetIntArray(sendindex1,sendtotal);
  ResetIntArray(sendindex2,sendtotal);

  time1 = MPI_Wtime();
  integrals_ov1 = fopen(zz, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
  nd6 = atoms_ax->bfnnumb_sh[j1];
  bfposa1 = atoms_ax->bfnposn_sh[j1];
  dim2a = nd6 * mpA;
  AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
  ResetDoubleArray(temp2a_buffer,&dim2a);
  fread(temp2a_buffer,sizeof(double),dim2a,integrals_ov1);
  for (i = 0; i < job->numtasks; i++) {
    count1[i] = 0;
   }
  for (il = 0; il < mpA; il++) {
    I1 = global_row_index + il;
    m = I1 / (nbands - fermi->occupied[0]);
    s = I1  - m  * (nbands - fermi->occupied[0]);
    for (jl = 0; jl < nqA; jl++) {
      I2 = global_col_index + jl;
      r = I2 / (nbands - fermi->occupied[0]);
      n = I2  - r  * (nbands - fermi->occupied[0]);
      I3 = m * nvir + n;
      I4 = r * nvir + s;
      row_proc = (I3 / *nbsize) % nprow;
      col_proc = (I4 / *nbsize) % npcol;
      cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
      for (a1 = 0; a1 < nd6; a1++) {
        senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] -= temp2a_buffer[a1 * mpA + il] * \
        temp3_buffer[(bfposa1 + a1) * nqA + jl];
       }
      (count1[cblacs_taskid])++;
      }
     }
    DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
   }
    fclose(integrals_ov1);
    for (i = 0; i < job->numtasks; i++) {
      count1[i] = 0;
     }
  for (il = 0; il < mpA; il++) {
    I1 = global_row_index + il;
    m = I1 / (nbands - fermi->occupied[0]);
    s = I1  - m  * (nbands - fermi->occupied[0]);
    for (jl = 0; jl < nqA; jl++) {
      I2 = global_col_index + jl;
      r = I2 / (nbands - fermi->occupied[0]);
      n = I2  - r  * (nbands - fermi->occupied[0]);
      I3 = m * nvir + n;
      I4 = r * nvir + s;
      row_proc = (I3 / *nbsize) % nprow;
      col_proc = (I4 / *nbsize) % npcol;
      cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
      il1 = *nbsize * (I3 / (*nbsize * nprow)) + I3 % *nbsize;
      jl1 = *nbsize * (I4 / (*nbsize * npcol)) + I4 % *nbsize;
      sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il1;
      sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl1;
      (count1[cblacs_taskid])++;
      }
     }
    time2 = MPI_Wtime() - time1;
    if (job->taskid == 0) printf("senddata_ib_ja                             %10.2f\n",time2);

}

void senddata_screened_ij_ab(DoubleMatrix *M_buffer1, DoubleMatrix *M_buffer2, double *senddata, int *sendtotal, int *sendindex1, int *sendindex2, int *senddispls, int *ictxt, int *nbsize, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, il, jl, l, m, n, r, s;
int I1, I2, row_proc, col_proc, cblacs_taskid;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int occ_pair, vir_pair;
int count1[job->numtasks];
int begin_vir_row[job->numtasks], end_vir_row[job->numtasks];
int begin_vir_array[job->numtasks], end_vir_array[job->numtasks];
int begin_occ_row[job->numtasks], end_occ_row[job->numtasks];
int begin_occ_array[job->numtasks], end_occ_array[job->numtasks];
double time1, time2;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
  mpi_begin_end_lower_triangle(begin_occ_row, end_occ_row, begin_occ_array, end_occ_array, nocc, nprow, job, file);
  mpi_begin_end_lower_triangle(begin_vir_row, end_vir_row, begin_vir_array, end_vir_array, nvir, npcol, job, file);

  ResetDoubleArray(senddata,sendtotal);
  ResetIntArray(sendindex1,sendtotal);
  ResetIntArray(sendindex2,sendtotal);

  time1 = MPI_Wtime();
  for (i = 0; i < job->numtasks; i++) {
    count1[i] = 0;
   }
  for (m = begin_occ_row[myrow]; m < end_occ_row[myrow]; m++) {
    for (r = 0; r <= m; r++) {
      occ_pair = (m * (m + 1)) / 2 + r - begin_occ_array[myrow];
      for (n = begin_vir_row[mycol]; n < end_vir_row[mycol]; n++) {
        for (s = 0; s <= n; s++) {
          vir_pair = (n * (n + 1)) / 2 + s - begin_vir_array[mycol];
          I1 = m * nvir + n;
          I2 = r * nvir + s;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
          jl = *nbsize * (I2 / (*nbsize * npcol)) + I2 % *nbsize;
          sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il;
          sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl;
          //for (l = 0; l < ntransitions; l++) {
          for (l = 0; l < job->rpa_lim; l++) {
          senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] += two * M_buffer1->a[occ_pair][l] * M_buffer2->a[vir_pair][l];
         }
          (count1[cblacs_taskid])++;
          if (r < m) {
          I1 = r * nvir + n;
          I2 = m * nvir + s;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
          jl = *nbsize * (I2 / (*nbsize * npcol)) + I2 % *nbsize;
          sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il;
          sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl;
          //for (l = 0; l < ntransitions; l++) {
          for (l = 0; l < job->rpa_lim; l++) {
          senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] += two * M_buffer1->a[occ_pair][l] * M_buffer2->a[vir_pair][l];
         }
          (count1[cblacs_taskid])++;
          }
          if (s < n) {
          I1 = m * nvir + s;
          I2 = r * nvir + n;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
          jl = *nbsize * (I2 / (*nbsize * npcol)) + I2 % *nbsize;
          sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il;
          sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl;
          //for (l = 0; l < ntransitions; l++) {
          for (l = 0; l < job->rpa_lim; l++) {
          senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] += two * M_buffer1->a[occ_pair][l] * M_buffer2->a[vir_pair][l];
         }
          (count1[cblacs_taskid])++;
          }
          if (r < m && s < n) {
          I1 = r * nvir + s;
          I2 = m * nvir + n;
          row_proc = (I1 / *nbsize) % nprow;
          col_proc = (I2 / *nbsize) % npcol;
          cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
          il = *nbsize * (I1 / (*nbsize * nprow)) + I1 % *nbsize;
          jl = *nbsize * (I2 / (*nbsize * npcol)) + I2 % *nbsize;
          sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il;
          sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl;
          //for (l = 0; l < ntransitions; l++) {
          for (l = 0; l < job->rpa_lim; l++) {
          senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] += two * M_buffer1->a[occ_pair][l] * M_buffer2->a[vir_pair][l];
         }
          (count1[cblacs_taskid])++;
          }
         }
        }
       }
      }
    time2 = MPI_Wtime() - time1;
    if (job->taskid == 0) printf("senddata_screened_ij_ab                    %10.2f\n",time2);

}

void senddata_screened_ib_ja(DoubleMatrix *M_buffer1, DoubleMatrix *M_buffer2, double *senddata, int *sendtotal, int *sendindex1, int *sendindex2, int *senddispls, int *ictxt, int *nbsize, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, il, jl, il1, jl1, l, m, n, r, s;
int I1, I2, I3, I4, row_proc, col_proc, cblacs_taskid;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int count1[job->numtasks];
int global_row_index, global_col_index;
double time1, time2;

  time1 = MPI_Wtime();
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  global_row_index = 0;
  global_col_index = 0;
  for (i = 0; i < myrow; i++) global_row_index += numroc_(&ntransitions, nbsize, &i, &izero, &nprow);
  for (i = 0; i < mycol; i++) global_col_index += numroc_(&ntransitions, nbsize, &i, &izero, &npcol);

  ResetDoubleArray(senddata,sendtotal);
  ResetIntArray(sendindex1,sendtotal);
  ResetIntArray(sendindex2,sendtotal);

  time1 = MPI_Wtime();
  for (i = 0; i < job->numtasks; i++) {
    count1[i] = 0;
   }
  for (il = 0; il < mpA; il++) {
    I1 = global_row_index + il;
    m = I1 / (nbands - fermi->occupied[0]);
    s = I1  - m  * (nbands - fermi->occupied[0]);
    for (jl = 0; jl < nqA; jl++) {
      I2 = global_col_index + jl;
      r = I2 / (nbands - fermi->occupied[0]);
      n = I2  - r  * (nbands - fermi->occupied[0]);
      I3 = m * nvir + n;
      I4 = r * nvir + s;
      row_proc = (I3 / *nbsize) % nprow;
      col_proc = (I4 / *nbsize) % npcol;
      cblacs_taskid = Cblacs_pnum(*ictxt,row_proc,col_proc);
      il1 = *nbsize * (I3 / (*nbsize * nprow)) + I3 % *nbsize;
      jl1 = *nbsize * (I4 / (*nbsize * npcol)) + I4 % *nbsize;
      sendindex1[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = il1;
      sendindex2[senddispls[cblacs_taskid] + count1[cblacs_taskid]] = jl1;
      //for (l = 0; l < ntransitions; l++) {
      for (l = 0; l < job->rpa_lim; l++) {
        senddata[senddispls[cblacs_taskid] + count1[cblacs_taskid]] += two * M_buffer1->a[il][l] * M_buffer2->a[jl][l];
      }
     (count1[cblacs_taskid])++;
    }
   }
    time2 = MPI_Wtime() - time1;
    if (job->taskid == 0) printf("senddata_screened_ib_ja                    %10.2f\n",time2);

}

void recvdata_buffer(double *Hamiltonian, IntMatrix *array_sizes, int *mpA, int *sendtotal, double *senddata, int *sendindex1, int *sendindex2, int *sendcounts, int *senddispls, int *recvtotal, double* recvdata, int *recvindex1, int *recvindex2, int *recvcounts, int *recvdispls, int *ictxt, int *nbsize, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, il, j, jl;
int count;
double time1, time2;

  time1 = MPI_Wtime();

  AllocateDoubleArray(&recvdata,recvtotal,job);
  ResetDoubleArray(recvdata,recvtotal);
  MPI_Alltoallv(senddata,sendcounts,senddispls,MPI_DOUBLE,recvdata,recvcounts,recvdispls,MPI_DOUBLE,MPI_COMM_WORLD);
  DestroyDoubleArray(&senddata,sendtotal,job);
  //if (job->taskid == 0) printf("senddata\n");

  AllocateIntArray(&recvindex1,recvtotal,job);
  ResetIntArray(recvindex1,recvtotal);
  MPI_Alltoallv(sendindex1,sendcounts,senddispls,MPI_INT,recvindex1,recvcounts,recvdispls,MPI_INT,MPI_COMM_WORLD);
  DestroyIntArray(&sendindex1,sendtotal,job);
  //if (job->taskid == 0) printf("sendindex1\n");

  AllocateIntArray(&recvindex2,recvtotal,job);
  ResetIntArray(recvindex2,recvtotal);
  MPI_Alltoallv(sendindex2,sendcounts,senddispls,MPI_INT,recvindex2,recvcounts,recvdispls,MPI_INT,MPI_COMM_WORLD);
  DestroyIntArray(&sendindex2,sendtotal,job);
  //if (job->taskid == 0) printf("sendindex2\n");

  count = 0;
  for (i = 0; i < job->numtasks; i++) {
    for (j = 0; j <array_sizes->a[i][job->taskid]; j++) {
      il = recvindex1[count];
      jl = recvindex2[count];
      Hamiltonian[il + *mpA * jl] += recvdata[count];
      count++;
     }
    }

  DestroyDoubleArray(&recvdata,recvtotal,job);
  DestroyIntArray(&recvindex1,recvtotal,job);
  DestroyIntArray(&recvindex2,recvtotal,job);

  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("recvdata_buffer                            %10.2f\n",time2);

}

void integrals_occ_vir_occ_vir(int *ictxt, int *nbsize, int block_cyclic, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int a1, i1, j1, k1, l1, j2, kp, lp;
int j, l, m, n, q, r, s;
int nd4, nd5, nd6;
int op, pm;
int dim1a, dim456, dim2a, dim2b, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int bfposa1, bfposk1, bfposl1;
int count;
int ione = 1;
int begin_evec[job->numtasks], end_evec[job->numtasks];
int I1, I2, il, jl, mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int active_proc[atoms_ax->number_of_atoms_in_unit_cell];
int temp1[atoms_ax->number_of_atoms_in_unit_cell];
int begin_ax[job->numtasks], end_ax[job->numtasks];
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
double *three_centre_integrals, *reduced_three_centre_integrals;
double time1, time2, time3, time4;
double *temp1a_buffer, *temp1b_buffer, *temp2a_buffer, *temp2b_buffer;
double *eigvec, *eigvec5;
FILE *integrals, *integrals_ov1, *integrals_ov2;
TRIPLE_TRAN triple;
INTEGRAL_LIST integral_list;
MPI_File fh;

char buf2[110], xy[14] = "/scf_evec_spk";
char xx[4], yy[24] = "integrals_3c.", bc[4] = "bc_", yz[24] = "integrals_", zz[24] = "integrals_";
char xo[5] = "ov1.", xv[5] = "ov2.";

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xy);
  sprintf(xx, "%d", job->taskid);
  if (block_cyclic) {
  strcat(yz,bc);
  strcat(zz,bc);
 }
  strcat(yy,xx);
  strcat(yz,xo);
  strcat(zz,xv);
  strcat(yz,xx);
  strcat(zz,xx);
  integrals_ov1 = fopen(yz, "wb");
  fclose(integrals_ov1);
  integrals_ov2 = fopen(zz, "wb");
  fclose(integrals_ov2);

  time1 = MPI_Wtime();
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
  mpi_begin_end(begin_ax, end_ax, atoms_ax->number_of_atoms_in_unit_cell, job->numtasks, job, file);
  for (i1 = 0; i1 < atoms_ax->number_of_atoms_in_unit_cell; i1++) { temp1[i1] = 0; active_proc[i1] = 0;
  if (i1 >= begin_ax[job->taskid] && i1 < end_ax[job->taskid]) temp1[i1] = job->taskid; }
  MPI_Allreduce(temp1,active_proc,atoms_ax->number_of_atoms_in_unit_cell,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  mpi_begin_end(begin_evec, end_evec, nbands, job->numtasks, job, file);
  //printf("begin_ax %3d end_ax %3d task %3d\n",begin_ax[job->taskid],end_ax[job->taskid],job->taskid);
  //printf("begin_evec %3d end_evec %3d task %3d\n",begin_evec[job->taskid],end_evec[job->taskid],job->taskid);
  //printf("mpA %3d nqA %3d %3d %3d %3d %3d  %3d %8d\n",mpA,nqA,myrow,mycol,nprow,npcol,nbands,ntransitions);

  int global_row_index, global_col_index, i;
  global_row_index = 0;
  global_col_index = 0;
  if (!block_cyclic) {
  for (i = 0; i < myrow; i++) global_row_index += numroc_(&ntransitions, nbsize, &i, &izero, &nprow);
  for (i = 0; i < mycol; i++) global_col_index += numroc_(&ntransitions, nbsize, &i, &izero, &npcol);
  //printf("%3d   %3d %3d  %3d %3d  %3d %3d\n",job->taskid,myrow,mycol,global_row_index,global_col_index,\
  numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow),numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol));
 }

  time1 = MPI_Wtime();
  //dim2a = dim1ax * mpA;
  //dim2b = nqA * dim1ax;
  //AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
  //AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
  //ResetDoubleArray(temp2a_buffer,&dim2a);
  //ResetDoubleArray(temp2b_buffer,&dim2b);
  //int dim5 = atoms->number_of_sh_bfns_in_unit_cell * (end_evec[job->taskid] - begin_evec[job->taskid]);
  //AllocateDoubleArray(&eigvec5,&dim5,job);
  AllocateDoubleArray(&eigvec,&dim1,job);
  //AllocateDoubleArray(&eigvec,&atoms->number_of_sh_bfns_in_unit_cell,job);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  integrals = fopen(yy, "rb");
  time1 = MPI_Wtime();
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    time3 = MPI_Wtime();
    nd6 = atoms_ax->bfnnumb_sh[j1];
    dim1a = dim1 * nd6 * nbands; 
    AllocateDoubleArray(&temp1a_buffer,&dim1a,job);
    AllocateDoubleArray(&temp1b_buffer,&dim1a,job);
    ResetDoubleArray(temp1a_buffer,&dim1a);
    ResetDoubleArray(temp1b_buffer,&dim1a);
    //if (job->taskid == 0) printf("%3d ",j1);
    count_triples1_reversed(j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
    allocate_TRIPLE_TRAN(&triple, job, file);
    generate_triples1_reversed(j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
    for (j = 0; j < triple.nump; j++) {
      q = triple.posn[j];
      nd4 = atoms->bfnnumb_sh[triple.cell1[q]];
      nd5 = atoms->bfnnumb_sh[triple.cell2[q]];
      nd6 = atoms_ax->bfnnumb_sh[triple.cell3[q]];
      dim456 = nd4 * nd5 * nd6;
      AllocateDoubleArray(&reduced_three_centre_integrals,&dim456,job);
      ResetDoubleArray(reduced_three_centre_integrals,&dim456);
      AllocateDoubleArray(&three_centre_integrals,&dim456,job);
      if (active_proc[j1] == job->taskid) {
      fread(&integral_list.num,sizeof(int),1,integrals);
      allocate_integral_list(&integral_list, integral_list.num, job, file);
      read_unpack_molecule_3c_integrals(&integral_list, integrals, job, file);
      for (jl = 0; jl < integral_list.num; jl++) {
        k1 = integral_list.i[jl];
        l1 = integral_list.j[jl];
        a1 = integral_list.k[jl];
        reduced_three_centre_integrals[k1 * nd5 * nd6 + l1 * nd6 + a1] = integral_list.value[jl];
       }
      MPI_Bcast(reduced_three_centre_integrals,dim456,MPI_DOUBLE,active_proc[j1],MPI_COMM_WORLD);
      free_integral_list(&integral_list,job);
     }
      else if (active_proc[j1] != job->taskid)
      MPI_Bcast(reduced_three_centre_integrals,dim456,MPI_DOUBLE,active_proc[j1],MPI_COMM_WORLD);
        //MPI_File_seek(fh, begin_evec[job->taskid] * dim1 * sizeof(double), MPI_SEEK_SET) ;
        //MPI_File_seek(fh, begin_evec[job->taskid] * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET) ;
        //MPI_File_read(fh, eigvec, dim5, MPI_DOUBLE, MPI_STATUS_IGNORE);
      for (il = begin_evec[job->taskid]; il < end_evec[job->taskid]; il++) {
        //MPI_File_seek(fh, il * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET) ;
        MPI_File_seek(fh, il * dim1 * sizeof(double), MPI_SEEK_SET) ;
        MPI_File_read(fh, eigvec, dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
        for (j2 = 0; j2 < triple.numb[j]; j2++) {
          kp = triple.cell1[q + j2];
          lp = triple.cell2[q + j2];
          op = triple.k[q + j2];
          pm = triple.p[q + j2];
          nd4 = atoms->bfnnumb_sh[kp];
          nd5 = atoms->bfnnumb_sh[lp];
          nd6 = atoms_ax->bfnnumb_sh[triple.cell3[q + j2]];
          bfposk1 = atoms->bfnposn_sh[kp]; 
          bfposl1 = atoms->bfnposn_sh[lp];
          bfposa1 = atoms_ax->bfnposn_sh[j1];
          rotate_permute_triple_ax_reversed2(&kp,&lp,&j1,&op,&pm,reduced_three_centre_integrals, \
          three_centre_integrals,atom_p,atoms,shells,atoms_ax,shells_ax,symmetry,job,file);
          count = 0;
          for (k1 = 0; k1 < nd4; k1++) { 
            for (l1 = 0; l1 < nd5; l1++) {
              for (a1 = 0; a1 < nd6; a1++) {
                //temp1b_buffer[il * dim1 * nd6 + (bfposl1 + l1) * nd6 + a1] += three_centre_integrals[count] * \
                eigvec[(il - begin_evec[job->taskid]) * dim1 + bfposk1 + k1];
                temp1b_buffer[il * dim1 * nd6 + (bfposl1 + l1) * nd6 + a1] += three_centre_integrals[count] * eigvec[bfposk1 + k1];
                count++;
               }
              }
             }
            } // close loop on j2
           } // close loop on il
            DestroyDoubleArray(&reduced_three_centre_integrals,&dim456,job);
            DestroyDoubleArray(&three_centre_integrals,&dim456,job);
           } // close loop on j
      MPI_Allreduce(temp1b_buffer,temp1a_buffer,dim1a,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      DestroyDoubleArray(&temp1b_buffer,&dim1a,job);

      dim2a = nd6 * mpA;
      AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
      ResetDoubleArray(temp2a_buffer,&dim2a);

      for (il = 0; il < mpA; il++) {
        if (block_cyclic) I1 = nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        else I1 = global_row_index + il;
        m = I1 / (nbands - fermi->occupied[0]);
        n = I1  - m  * (nbands - fermi->occupied[0]) + fermi->occupied[0];
	//MPI_File_seek(fh, n * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET) ;
        MPI_File_seek(fh, n * dim1 * sizeof(double), MPI_SEEK_SET) ;
        MPI_File_read(fh, eigvec, dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
          for (k1 = 0; k1 < dim1; k1++) {
            for (a1 = 0; a1 < nd6; a1++) {
              //temp2a_buffer[il * dim1ax + bfposa1 + a1] += temp1a_buffer[m * dim1 * nd6 + k1 * nd6 + a1] * eigvec[k1];
              temp2a_buffer[a1 * mpA + il] += temp1a_buffer[m * dim1 * nd6 + k1 * nd6 + a1] * eigvec[k1];
              //temp2a_buffer[(bfposa1 + a1) * mpA + il] += temp1a_buffer[m * dim1 * nd6 + k1 * nd6 + a1] * eigvec[k1];
             }
            }
           }

        integrals_ov1 = fopen(yz, "ab");
        fwrite(temp2a_buffer, sizeof(double), dim2a, integrals_ov1);
        fclose(integrals_ov1);
        DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
        dim2b = nqA * nd6;
        AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
        ResetDoubleArray(temp2b_buffer,&dim2b);

      for (il = 0; il < nqA; il++) {
        if (block_cyclic) I2 = npcol * *nbsize * (il / *nbsize) + il % *nbsize + ((npcol + mycol) % npcol) * *nbsize;
        else I2 = global_col_index + il;
        r = I2 / (nbands - fermi->occupied[0]);
        s = I2  - r  * (nbands - fermi->occupied[0]) + fermi->occupied[0];
        //MPI_File_seek(fh, s * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET) ;
        MPI_File_seek(fh, s * dim1 * sizeof(double), MPI_SEEK_SET) ;
        MPI_File_read(fh, eigvec, dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
        for (k1 = 0; k1 < dim1; k1++) {
          for (a1 = 0; a1 < nd6; a1++) {
            temp2b_buffer[a1 * nqA + il] += temp1a_buffer[r * dim1 * nd6 + k1 * nd6 + a1] * eigvec[k1];
            //temp2b_buffer[(bfposa1 + a1) * nqA + il] += temp1a_buffer[r * dim1 * nd6 + k1 * nd6 + a1] * eigvec[k1];
           }
          }
         }

        integrals_ov2 = fopen(zz, "ab");
        fwrite(temp2b_buffer, sizeof(double), dim2b, integrals_ov2);
        fclose(integrals_ov2);
        DestroyDoubleArray(&temp2b_buffer,&dim2b,job);

        free_TRIPLE_TRAN(&triple,job);
        DestroyDoubleArray(&temp1a_buffer,&dim1a,job);
        time4 = MPI_Wtime() - time3;
        //if (job->taskid == 0) printf("%f\n",time4);
       } // close loop on j1
        //DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
        //DestroyDoubleArray(&temp2b_buffer,&dim2b,job);
        DestroyDoubleArray(&eigvec,&dim1,job);
        //DestroyDoubleArray(&eigvec,&atoms->number_of_sh_bfns_in_unit_cell,job);
        fclose(integrals);
        MPI_File_close(&fh);
        //if (job->taskid == 0) printf("end intermediates %3d %f\n",job->taskid,time2);
        time2 = MPI_Wtime() - time1;
        if (job->taskid == 0) printf("integrals_occ_vir_occ_vir                  %10.2f\n",time2);

}

void integrals_occ_occ_vir_vir(int *ictxt, int *nbsize, int block_cyclic, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int a1, i1, j1, k1, l1, j2, kp, lp;
int i, j, l, m, n, q, r, s;
int nd4, nd5, nd6;
int op, pm;
int dim1a, dim1b, dim1c, dim2a, dim2b, dim456, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int bfposa1, bfposk1, bfposl1;
int count, task;
int ione = 1;
int begin_vir_row[job->numtasks], end_vir_row[job->numtasks];
int begin_vir_array[job->numtasks], end_vir_array[job->numtasks];
int begin_occ_row[job->numtasks], end_occ_row[job->numtasks];
int begin_occ_array[job->numtasks], end_occ_array[job->numtasks];
int I1, I2, il, jl, mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
int active_proc[atoms_ax->number_of_atoms_in_unit_cell];
int temp1[atoms_ax->number_of_atoms_in_unit_cell];
int begin_ax[job->numtasks], end_ax[job->numtasks];
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int occ_pair, vir_pair;
int ntransitions = nocc * nvir;
double *three_centre_integrals, *reduced_three_centre_integrals;
double time1, time2;
double time3, time4;
double *temp1a_buffer, *temp1b_buffer, *temp1c_buffer;
double *temp2a_buffer, *temp2b_buffer;
double *eigvec;
char buf2[110], xy[14] = "/scf_evec_spk";
char xx[4], yy[24] = "integrals_3c.", bc[4] = "bc_", yz[24] = "integrals_", zz[24] = "integrals_";
char xo[4] = "oo.", xv[4] = "vv.";
FILE *integrals, *integrals_temp2a, *integrals_temp2b;
TRIPLE_TRAN triple;
INTEGRAL_LIST integral_list;
MPI_File fh;

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xy);
  sprintf(xx, "%d", job->taskid);
  if (block_cyclic) {
  strcat(yz,bc);
  strcat(zz,bc);
 }
  strcat(yy,xx);
  strcat(yz,xo);
  strcat(zz,xv);
  strcat(yz,xx);
  strcat(zz,xx);
  integrals_temp2a = fopen(yz, "wb");
  fclose(integrals_temp2a);
  integrals_temp2b = fopen(zz, "wb");
  fclose(integrals_temp2b);

  AllocateDoubleArray(&eigvec,&dim1,job);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;

  time1 = MPI_Wtime();
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
  mpi_begin_end(begin_ax, end_ax, atoms_ax->number_of_atoms_in_unit_cell, job->numtasks, job, file);

  if (!block_cyclic) {
    nprow = job->numtasks;
    npcol = job->numtasks;
    myrow = job->taskid;
    mycol = job->taskid;
   }

  for (i1 = 0; i1 < atoms_ax->number_of_atoms_in_unit_cell; i1++) { temp1[i1] = 0; active_proc[i1] = 0;
  if (i1 >= begin_ax[job->taskid] && i1 < end_ax[job->taskid]) temp1[i1] = job->taskid; }
  MPI_Allreduce(temp1,active_proc,atoms_ax->number_of_atoms_in_unit_cell,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  //printf("mpA %3d nqA %3d %3d %3d %3d %3d  %3d %8d\n",mpA,nqA,myrow,mycol,nprow,npcol,nbands,ntransitions);
  mpi_begin_end_lower_triangle(begin_occ_row, end_occ_row, begin_occ_array, end_occ_array, nocc, nprow, job, file);
  mpi_begin_end_lower_triangle(begin_vir_row, end_vir_row, begin_vir_array, end_vir_array, nvir, npcol, job, file);
  //mpi_begin_end(begin_evec, end_evec, nbands, job->numtasks, job, file);
  //printf("begin_ax %3d end_ax %3d task %3d\n",begin_ax[job->taskid],end_ax[job->taskid],job->taskid);
  //printf("begin_evec %3d end_evec %3d task %3d\n",begin_evec[job->taskid],end_evec[job->taskid],job->taskid);
  //printf("mpA %3d nqA %3d %3d %3d %3d %3d  %3d %8d\n",mpA,nqA,myrow,mycol,nprow,npcol,nbands,ntransitions);

//HERE
  int dim5 = dim1 * (end_occ_row[myrow] - begin_occ_row[myrow]);
  int dim6 = dim1 * (end_vir_row[mycol] - begin_vir_row[mycol]);
  double *eigvec1a, *eigvec1b;
  //printf("%3d\n",dim5);
  AllocateDoubleArray(&eigvec1a,&dim5,job);
  AllocateDoubleArray(&eigvec1b,&dim6,job);
//HERE

  integrals = fopen(yy, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    time3 = MPI_Wtime();
    nd6 = atoms_ax->bfnnumb_sh[j1];
    dim1a = (end_occ_row[myrow] - begin_occ_row[myrow]) * dim1 * nd6; 
    dim1b = (end_vir_row[mycol] - begin_vir_row[mycol]) * dim1 * nd6; 
    //if (job->taskid == 0) printf("%7d %7d %7d %7d  %7d %6d\n",dim1a,dim1b,dim2a,dim2b,dim1a+dim1b+dim2a+dim2b,job->taskid);
    AllocateDoubleArray(&temp1a_buffer,&dim1a,job);
    ResetDoubleArray(temp1a_buffer,&dim1a);
    AllocateDoubleArray(&temp1b_buffer,&dim1b,job);
    ResetDoubleArray(temp1b_buffer,&dim1b);
    //if (job->taskid == 0) printf("%3d ",j1);
    count_triples1_reversed(j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
    allocate_TRIPLE_TRAN(&triple, job, file);
    generate_triples1_reversed(j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
    for (j = 0; j < triple.nump; j++) {
      q = triple.posn[j];
      nd4 = atoms->bfnnumb_sh[triple.cell1[q]];
      nd5 = atoms->bfnnumb_sh[triple.cell2[q]];
      nd6 = atoms_ax->bfnnumb_sh[triple.cell3[q]];
      dim456 = nd4 * nd5 * nd6;
      AllocateDoubleArray(&reduced_three_centre_integrals,&dim456,job);
      ResetDoubleArray(reduced_three_centre_integrals,&dim456);
      AllocateDoubleArray(&three_centre_integrals,&dim456,job);
      if (active_proc[j1] == job->taskid) {
      fread(&integral_list.num,sizeof(int),1,integrals);
      allocate_integral_list(&integral_list, integral_list.num, job, file);
      read_unpack_molecule_3c_integrals(&integral_list, integrals, job, file);
      for (jl = 0; jl < integral_list.num; jl++) {
        k1 = integral_list.i[jl];
        l1 = integral_list.j[jl];
        a1 = integral_list.k[jl];
        reduced_three_centre_integrals[k1 * nd5 * nd6 + l1 * nd6 + a1] = integral_list.value[jl];
       }
      MPI_Bcast(reduced_three_centre_integrals,dim456,MPI_DOUBLE,active_proc[j1],MPI_COMM_WORLD);
      free_integral_list(&integral_list,job);
     }
      else if (active_proc[j1] != job->taskid)
      MPI_Bcast(reduced_three_centre_integrals,dim456,MPI_DOUBLE,active_proc[j1],MPI_COMM_WORLD);
      //for (il = 0; il < end_vir[job->taskid] - begin_vir[job->taskid]; il++) {
      for (j2 = 0; j2 < triple.numb[j]; j2++) {
        kp = triple.cell1[q + j2];
        lp = triple.cell2[q + j2];
        op = triple.k[q + j2];
        pm = triple.p[q + j2];
        nd4 = atoms->bfnnumb_sh[kp];
        nd5 = atoms->bfnnumb_sh[lp];
        nd6 = atoms_ax->bfnnumb_sh[triple.cell3[q + j2]];
        bfposk1 = atoms->bfnposn_sh[kp]; 
        bfposl1 = atoms->bfnposn_sh[lp];
        bfposa1 = atoms_ax->bfnposn_sh[j1];
        rotate_permute_triple_ax_reversed2(&kp,&lp,&j1,&op,&pm,reduced_three_centre_integrals, \
        three_centre_integrals,atom_p,atoms,shells,atoms_ax,shells_ax,symmetry,job,file);
    //MPI_File_seek(fh, begin_occ_row[myrow] * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET);
    MPI_File_seek(fh, begin_occ_row[myrow] * dim1 * sizeof(double), MPI_SEEK_SET);
    MPI_File_read(fh, eigvec1a, (end_occ_row[myrow] - begin_occ_row[myrow]) * dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
        for (i = 0; i < end_occ_row[myrow] - begin_occ_row[myrow]; i++) {
        //for (il = 0; il < end_occ_row[job->taskid] - begin_occ_row[job->taskid]; il++) {
        //for (il = 0; il < end_occ[job->taskid] - begin_occ[job->taskid]; il++) {
        //MPI_File_seek(fh, (begin_occ_row[myrow] + i) * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET);
        //MPI_File_seek(fh, (begin_occ_row[myrow] + i) * dim1 * sizeof(double), MPI_SEEK_SET);
        //MPI_File_read(fh, eigvec, dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
        count = 0;
        for (k1 = 0; k1 < nd4; k1++) { 
          for (l1 = 0; l1 < nd5; l1++) {
            for (a1 = 0; a1 < nd6; a1++) {
              temp1a_buffer[i * dim1 * nd6 + (bfposl1 + l1) * nd6 + a1] += three_centre_integrals[count] * \
              eigvec1a[i * dim1 + bfposk1 + k1];
              //temp1a_buffer[i * dim1 * nd6 + (bfposl1 + l1) * nd6 + a1] += three_centre_integrals[count] * eigvec[bfposk1 + k1];
              count++;
             }
            }
           }
          } // close loop on i
    //MPI_File_seek(fh, (nocc + begin_vir_row[mycol]) * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET);
    MPI_File_seek(fh, (nocc + begin_vir_row[mycol]) * dim1 * sizeof(double), MPI_SEEK_SET);
    MPI_File_read(fh, eigvec1b, dim6, MPI_DOUBLE, MPI_STATUS_IGNORE);
        for (i = 0; i < end_vir_row[mycol] - begin_vir_row[mycol]; i++) {
          ////MPI_File_seek(fh, (nocc + begin_vir_row[mycol] + i) * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET);
          ////MPI_File_seek(fh, (nocc + begin_vir_row[mycol] + i) * dim1 * sizeof(double), MPI_SEEK_SET);
          ////MPI_File_read(fh, eigvec, dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
          count = 0;
          for (k1 = 0; k1 < nd4; k1++) { 
            for (l1 = 0; l1 < nd5; l1++) {
              for (a1 = 0; a1 < nd6; a1++) {
                temp1b_buffer[i * dim1 * nd6 + (bfposl1 + l1) * nd6 + a1] += three_centre_integrals[count] * \
                eigvec1b[i * dim1 + bfposk1 + k1];
                //temp1b_buffer[i * dim1 * nd6 + (bfposl1 + l1) * nd6 + a1] += three_centre_integrals[count] * eigvec[bfposk1 + k1];
                count++;
               }
              }
             }
            } // close loop on i
           } // close loop on j2
          DestroyDoubleArray(&reduced_three_centre_integrals,&dim456,job);
          DestroyDoubleArray(&three_centre_integrals,&dim456,job);
         } // close loop on j
          dim2a = nd6 * (end_occ_array[myrow] - begin_occ_array[myrow]);
          AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
          ResetDoubleArray(temp2a_buffer,&dim2a);
          for (i = begin_occ_row[myrow]; i < end_occ_row[myrow]; i++) {
            for (j = 0; j <= i; j++) {
             occ_pair = (i * (i + 1)) / 2 + j - begin_occ_array[myrow];
             //MPI_File_seek(fh, j * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET) ;
             MPI_File_seek(fh, j * dim1 * sizeof(double), MPI_SEEK_SET) ;
             MPI_File_read(fh, eigvec, dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
             for (k1 = 0; k1 < dim1; k1++) {
               for (a1 = 0; a1 < nd6; a1++) {
                 temp2a_buffer[occ_pair * nd6 + a1] += temp1a_buffer[(i - begin_occ_row[myrow]) * dim1 * nd6 + k1 * nd6 + a1] * \
                 eigvec[k1];
                }
               }
              }
             }
          integrals_temp2a = fopen(yz, "ab");
          fwrite(temp2a_buffer, sizeof(double), dim2a, integrals_temp2a);
          fclose(integrals_temp2a);
          DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
          dim2b = nd6 * (end_vir_array[mycol] - begin_vir_array[mycol]);
          AllocateDoubleArray(&temp2b_buffer,&dim2b,job);
          ResetDoubleArray(temp2b_buffer,&dim2b);
          for (i = begin_vir_row[mycol]; i < end_vir_row[mycol]; i++) {
            for (j = 0; j <= i; j++) {
              vir_pair = (i * (i + 1)) / 2 + j - begin_vir_array[mycol];
              //MPI_File_seek(fh, (nocc + j) * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET) ;
              MPI_File_seek(fh, (nocc + j) * dim1 * sizeof(double), MPI_SEEK_SET) ;
              MPI_File_read(fh, eigvec, dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
              for (k1 = 0; k1 < dim1; k1++) {
                for (a1 = 0; a1 < nd6; a1++) {
                  temp2b_buffer[vir_pair * nd6 + a1] += temp1b_buffer[(i - begin_vir_row[mycol]) * dim1 * nd6 + k1 * nd6 + a1] * \
                  eigvec[k1];
                 }
                }
               }
              }
          //for (i= 0; i < dim2b; i++) fprintf(file.out,"%3d %10.4f\n",i,temp2b_buffer[i]);
          integrals_temp2b = fopen(zz, "ab");
          fwrite(temp2b_buffer, sizeof(double), dim2b, integrals_temp2b);
          fclose(integrals_temp2b);
          DestroyDoubleArray(&temp2b_buffer,&dim2b,job);

          free_TRIPLE_TRAN(&triple,job);
          DestroyDoubleArray(&temp1a_buffer,&dim1a,job);
          DestroyDoubleArray(&temp1b_buffer,&dim1b,job);
          time4 = MPI_Wtime() - time3;
          //if (job->taskid == 0) printf("%f\n",time4);
         } // close loop on j1
          fclose(integrals);
          DestroyDoubleArray(&eigvec1a,&dim5,job);
          DestroyDoubleArray(&eigvec1b,&dim6,job);
          DestroyDoubleArray(&eigvec,&dim1,job);
          MPI_File_close(&fh);
          time2 = MPI_Wtime() - time1;
          if (job->taskid == 0) printf("integrals_occ_occ_vir_vir                  %10.2f\n",time2);
          //if (job->taskid == 0) printf("block_cyclic %3d integrals_occ_occ_vir_vir %10.2f\n",block_cyclic,time2);

}

void integrals_self_energy(int *ictxt, int *nbsize, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int a1, i1, j1, jl, k1, l1, j2, kp, lp;
int j, l, m, n, q, r, s;
int nd4, nd5, nd6;
int op, pm;
int dim1a, dim456, dim2a, dim5, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int bfposa1, bfposk1, bfposl1;
int count;
int ione = 1;
int active_proc[atoms_ax->number_of_atoms_in_unit_cell];
int temp1[atoms_ax->number_of_atoms_in_unit_cell];
int begin_ax[job->numtasks], end_ax[job->numtasks];
int begin_occ[job->numtasks], end_occ[job->numtasks];
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
double *three_centre_integrals, *reduced_three_centre_integrals;
double time1, time2, time3, time4;
double *temp1a_buffer, *temp2a_buffer;
double *eigvec, *eigve1;
TRIPLE_TRAN triple;
INTEGRAL_LIST integral_list;

  mpi_begin_end(begin_ax, end_ax, atoms_ax->number_of_atoms_in_unit_cell, job->numtasks, job, file);
  for (i1 = 0; i1 < atoms_ax->number_of_atoms_in_unit_cell; i1++) { temp1[i1] = 0; active_proc[i1] = 0;
  if (i1 >= begin_ax[job->taskid] && i1 < end_ax[job->taskid]) temp1[i1] = job->taskid; }
  MPI_Allreduce(temp1,active_proc,atoms_ax->number_of_atoms_in_unit_cell,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  mpi_begin_end(begin_occ, end_occ, nocc, job->numtasks, job, file);
  //printf("begin_ax %3d end_ax %3d task %3d\n",begin_ax[job->taskid],end_ax[job->taskid],job->taskid);
  //printf("mpA %3d nqA %3d %3d %3d %3d %3d  %3d %8d\n",mpA,nqA,myrow,mycol,nprow,npcol,nbands,ntransitions);

  MPI_File fh;
  char buf2[110], xy[14] = "/scf_evec_spk";
  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xy);
  FILE *integrals, *integrals_occvir;
  char xx[4], yy[24] = "integrals_3c.", xz[28] = "integrals_self_energy.";
  sprintf(xx, "%d", job->taskid);
  integrals = fopen(strcat(yy,xx), "rb");
  fclose(integrals);
  strcat(xz,xx);
  integrals_occvir = fopen(xz, "wb");
  fclose(integrals_occvir);

  time1 = MPI_Wtime();
  dim5 = dim1 * (end_occ[job->taskid] - begin_occ[job->taskid]);
  AllocateDoubleArray(&eigvec,&dim5,job);
  AllocateDoubleArray(&eigve1,&dim1,job);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  integrals = fopen(yy, "rb");
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    time3 = MPI_Wtime();
    nd6 = atoms_ax->bfnnumb_sh[j1];
    dim1a = dim5 * nd6; 
    //dim1a = dim1 * nd6 * nbands; 
    AllocateDoubleArray(&temp1a_buffer,&dim1a,job);
    //AllocateDoubleArray(&temp1b_buffer,&dim1a,job);
    ResetDoubleArray(temp1a_buffer,&dim1a);
    //ResetDoubleArray(temp1b_buffer,&dim1a);
    //if (job->taskid == 0) printf("%3d ",j1);
    count_triples1_reversed(j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
    allocate_TRIPLE_TRAN(&triple, job, file);
    generate_triples1_reversed(j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
    for (j = 0; j < triple.nump; j++) {
      q = triple.posn[j];
      nd4 = atoms->bfnnumb_sh[triple.cell1[q]];
      nd5 = atoms->bfnnumb_sh[triple.cell2[q]];
      nd6 = atoms_ax->bfnnumb_sh[triple.cell3[q]];
      dim456 = nd4 * nd5 * nd6;
      AllocateDoubleArray(&reduced_three_centre_integrals,&dim456,job);
      ResetDoubleArray(reduced_three_centre_integrals,&dim456);
      AllocateDoubleArray(&three_centre_integrals,&dim456,job);
      if (active_proc[j1] == job->taskid) {
      fread(&integral_list.num,sizeof(int),1,integrals);
      allocate_integral_list(&integral_list, integral_list.num, job, file);
      read_unpack_molecule_3c_integrals(&integral_list, integrals, job, file);
      for (jl = 0; jl < integral_list.num; jl++) {
        k1 = integral_list.i[jl];
        l1 = integral_list.j[jl];
        a1 = integral_list.k[jl];
        reduced_three_centre_integrals[k1 * nd5 * nd6 + l1 * nd6 + a1] = integral_list.value[jl];
       }
      MPI_Bcast(reduced_three_centre_integrals,dim456,MPI_DOUBLE,active_proc[j1],MPI_COMM_WORLD);
      free_integral_list(&integral_list,job);
     }
      else if (active_proc[j1] != job->taskid)
      MPI_Bcast(reduced_three_centre_integrals,dim456,MPI_DOUBLE,active_proc[j1],MPI_COMM_WORLD);
      //MPI_File_seek(fh, begin_occ[job->taskid] * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET) ;
      MPI_File_seek(fh, begin_occ[job->taskid] * dim1 * sizeof(double), MPI_SEEK_SET) ;
      MPI_File_read(fh, eigvec, dim5, MPI_DOUBLE, MPI_STATUS_IGNORE);
      for (m = 0; m < end_occ[job->taskid] - begin_occ[job->taskid]; m++) {
      //for (il = begin_occ[job->taskid]; il < end_occ[job->taskid]; il++) {
        //MPI_File_seek(fh, il * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET) ;
        //MPI_File_seek(fh, il * dim1 * sizeof(double), MPI_SEEK_SET) ;
        //MPI_File_read(fh, eigvec, dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
        for (j2 = 0; j2 < triple.numb[j]; j2++) {
          kp = triple.cell1[q + j2];
          lp = triple.cell2[q + j2];
          op = triple.k[q + j2];
          pm = triple.p[q + j2];
          nd4 = atoms->bfnnumb_sh[kp];
          nd5 = atoms->bfnnumb_sh[lp];
          nd6 = atoms_ax->bfnnumb_sh[triple.cell3[q + j2]];
          bfposk1 = atoms->bfnposn_sh[kp]; 
          bfposl1 = atoms->bfnposn_sh[lp];
          bfposa1 = atoms_ax->bfnposn_sh[j1];
          rotate_permute_triple_ax_reversed2(&kp,&lp,&j1,&op,&pm,reduced_three_centre_integrals, \
          three_centre_integrals,atom_p,atoms,shells,atoms_ax,shells_ax,symmetry,job,file);
          count = 0;
          for (k1 = 0; k1 < nd4; k1++) { 
            for (l1 = 0; l1 < nd5; l1++) {
              for (a1 = 0; a1 < nd6; a1++) {
                temp1a_buffer[m * dim1 * nd6 + (bfposl1 + l1) * nd6 + a1] += three_centre_integrals[count] * \
                eigvec[m * dim1 + bfposk1 + k1];
                count++;
               }
              }
             }
            } // close loop on j2
           } // close loop on m
            DestroyDoubleArray(&reduced_three_centre_integrals,&dim456,job);
            DestroyDoubleArray(&three_centre_integrals,&dim456,job);
           } // close loop on j
      dim2a = (end_occ[job->taskid] - begin_occ[job->taskid]) * nvir * nd6;
      AllocateDoubleArray(&temp2a_buffer,&dim2a,job);
      ResetDoubleArray(temp2a_buffer,&dim2a);
      //MPI_File_seek(fh, nocc * dim1 * sizeof(MPI_DOUBLE), MPI_SEEK_SET) ;
      MPI_File_seek(fh, nocc * dim1 * sizeof(double), MPI_SEEK_SET) ;
      for (n = 0; n < nvir; n++) {
        MPI_File_read(fh, eigve1, dim1, MPI_DOUBLE, MPI_STATUS_IGNORE);
        for (m = 0; m < end_occ[job->taskid] - begin_occ[job->taskid]; m++) {
          for (k1 = 0; k1 < dim1; k1++) {
            for (a1 = 0; a1 < nd6; a1++) {
              //temp2a_buffer[a1 * ntransitions + m * nvir + n] += temp1a_buffer[m * dim1 * nd6 + k1 * nd6 + a1] * eigvec[k1];
              temp2a_buffer[m * nvir * nd6 + n * nd6 + a1] += temp1a_buffer[m * dim1 * nd6 + k1 * nd6 + a1] * eigve1[k1];
              //temp2a_buffer[a1 * mpA + il] += temp1a_buffer[m * dim1 * nd6 + k1 * nd6 + a1] * eigvec[k1];
             }
            }
           }
          }
        integrals_occvir = fopen(xz, "ab");
        fwrite(temp2a_buffer, sizeof(double), dim2a, integrals_occvir);
        fclose(integrals_occvir);
        DestroyDoubleArray(&temp2a_buffer,&dim2a,job);
        free_TRIPLE_TRAN(&triple,job);
        DestroyDoubleArray(&temp1a_buffer,&dim1a,job);
        time4 = MPI_Wtime() - time3;
        //if (job->taskid == 0) printf("%f\n",time4);
       } // close loop on j1
        DestroyDoubleArray(&eigvec,&dim5,job);
        DestroyDoubleArray(&eigve1,&dim1,job);
        fclose(integrals);
        MPI_File_close(&fh);
        time2 = MPI_Wtime() - time1;
        if (job->taskid == 0) printf("integrals_self_energy                      %10.2f\n",time2);

}


void integrals_occ_occ_vir_vir1(int *j1, MPI_File fh, double *integral_buffer, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int i, j, k, l, q;
int a1, k1, l1, kp, lp, j2;
int nd4, nd5, nd6;
int op, pm;
int dim1a, dim2, dim456, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int bfposa1, bfposk1, bfposl1;
int count;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
double *three_centre_integrals, *reduced_three_centre_integrals;
double time1, time2, time3, time4;
double *temp1_buffer;
double *eigvec;
char buf2[110], xy[14] = "/scf_evec_spk";
TRIPLE_TRAN triple;

  time1 = MPI_Wtime();
  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xy);
  dim2 = nbands * dim1;
  AllocateDoubleArray(&eigvec,&dim2,job);
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ; //read_write_eigenvectors writes from (fermi->bands[0] - 1)
  MPI_File_read(fh, eigvec, dim2, MPI_DOUBLE, MPI_STATUS_IGNORE);
  nd6 = atoms_ax->bfnnumb_sh[*j1];
  dim1a = nbands * dim1 * nd6; 
  AllocateDoubleArray(&temp1_buffer,&dim1a,job);
  ResetDoubleArray(temp1_buffer,&dim1a);
  if (job->taskid == 0) printf("J1 %3d ",*j1);
  count_triples1_reversed(*j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
  allocate_TRIPLE_TRAN(&triple, job, file);
  generate_triples1_reversed(*j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
  for (j = 0; j < triple.nump; j++) {
    q = triple.posn[j];
    nd4 = atoms->bfnnumb_sh[triple.cell1[q]];
    nd5 = atoms->bfnnumb_sh[triple.cell2[q]];
    dim456 = nd4 * nd5 * nd6;
    AllocateDoubleArray(&reduced_three_centre_integrals,&dim456,job);
    AllocateDoubleArray(&three_centre_integrals,&dim456,job);
    ResetDoubleArray(reduced_three_centre_integrals,&dim456);
    //three_centre_coulomb1_reversed2(j,&triple,reduced_three_centre_integrals,R,G,atoms,shells, \
    gaussians,atoms_ax, shells_ax,gaussians_ax,crystal,job,file);
    integrals_ij_alpha(j,&triple,reduced_three_centre_integrals,R,G,atoms,shells, \
    gaussians,atoms_ax, shells_ax,gaussians_ax,crystal,job,file);
    for (j2 = 0; j2 < triple.numb[j]; j2++) {
      kp = triple.cell1[q + j2];
      lp = triple.cell2[q + j2];
      op = triple.k[q + j2];
      pm = triple.p[q + j2];
      nd4 = atoms->bfnnumb_sh[kp];
      nd5 = atoms->bfnnumb_sh[lp];
      bfposk1 = atoms->bfnposn_sh[kp];
      bfposl1 = atoms->bfnposn_sh[lp];
      ResetDoubleArray(three_centre_integrals,&dim456);
      rotate_permute_triple_ax_reversed2(&kp,&lp,j1,&op,&pm,reduced_three_centre_integrals, \
      three_centre_integrals,atom_p,atoms,shells,atoms_ax,shells_ax,symmetry,job,file);
      count = 0;
      for (k1 = 0; k1 < nd4; k1++) {
        for (l1 = 0; l1 < nd5; l1++) {
          for (a1 = 0; a1 < nd6; a1++) {
            for (l = 0; l < nbands; l++) {
              temp1_buffer[(bfposk1 + k1) * nbands * nd6 + l * nd6 + a1] += three_centre_integrals[count] * \
              eigvec[l * dim1 + bfposl1 + l1];
              }
             count++;
            }
           }
          }
         } // close loop on j2
       DestroyDoubleArray(&reduced_three_centre_integrals,&dim456,job);
       DestroyDoubleArray(&three_centre_integrals,&dim456,job);
      } // close loop on j
      free_TRIPLE_TRAN(&triple,job);
      for (k = 0; k < nbands; k++) {
        for (k1 = 0; k1 < dim1; k1++) {
          for (l = 0; l < nbands; l++) {
            for (a1 = 0; a1 < nd6; a1++) {
              integral_buffer[k * nbands * nd6 + l * nd6 + a1] += temp1_buffer[k1 * nbands * nd6 + l * nd6 + a1] * \
              eigvec[k * dim1 + k1];
             }
            }
           }
          }
          DestroyDoubleArray(&temp1_buffer,&dim1a,job);
          DestroyDoubleArray(&eigvec,&dim2,job);
          time2 = MPI_Wtime() - time1;
          printf("integrals_occ_occ_vir_vir1                 %10.2f\n",time2);

}

void contract_coulomb_integrals(int *j1, int *j2, DoubleMatrix *V_inv, double *integral_buffer, double *integral_buffer2, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Multiply inverse square root of coulomb matrix into three centre integrals             *
  // ******************************************************************************************

int i, a1, a2;
int bfposa1, bfposa2, nd6a, nd6b;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;

  nd6a = atoms_ax->bfnnumb_sh[*j1];
  nd6b = atoms_ax->bfnnumb_sh[*j2];
  bfposa1 = atoms_ax->bfnposn_sh[*j1];
  bfposa2 = atoms_ax->bfnposn_sh[*j2];
  for (i = 0; i < nbands * nbands; i++) {
    for (a1 = 0; a1 < nd6a; a1++) {
      for (a2 = 0; a2 < nd6b; a2++) {
        //integral_buffer2[i * nd6a + a1] += V_inv->a[bfposa1 + a1][bfposa2 + a2] * integral_buffer[i * nd6b + a2];
        integral_buffer2[i * nd6b + a2] += V_inv->a[bfposa2 + a2][bfposa1 + a1] * integral_buffer[i * nd6a + a1]; // 1proc
       }
      }
     }

}
*/

/*
void rpa_bse_hamiltonian_in_core(double *Ham_buffer1, int *ictxt, int *nbsize, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int i, j, k, p;
int m, n, r, s, nd6, a1;
int I1, I2, il, jl;
int j1, j2, jj, kk;
int dima, dimb, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
int mpA, nqA, row, col, nprow, npcol, myrow, mycol, izero = 0;
int mpa, nqa, mpa_nqa, MPA[job->numtasks],NQA[job->numtasks];;
int cblacs_taskid, nprow_nbsize, nprow_myrow, npcol_nbsize, npcol_mycol;
int begin_j[job->numtasks], end_j[job->numtasks];
int *dim_send, *dim_recv, *dim_ham, *offset_j;
int *dim_send1;
int offset, offset1, dim3;
int num_proc = job->numtasks < dim1 ? job->numtasks : dim1;
int begin_h[job->numtasks], end_h[job->numtasks];
int sendtotal, sendcounts[job->numtasks], senddispls[job->numtasks];
int recvtotal, recvcounts[job->numtasks], recvdispls[job->numtasks];
double factor1;
double time1, time2, time3, time4;
double time5, time6, time7, time8, time9, time10;
double *integral_buffer, *integral_buffer1, *integral_buffer2, *Hamiltonian_buffer;
double *integrals1, *integrals2;
char buf2[110], xy[14] = "/scf_evec_spk";
DoubleMatrix *V_inv;
MPI_File fh;

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xy);
  
  // ******************************************************************************************
  // * Generate inverse of Coulomb matrix                                                     *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  AllocateDoubleMatrix(&V_inv,&dim1ax,&dim1ax,job);
  if (job->taskid == 0) 
  generate_coulomb_matrix_inverse(V_inv,fermi,atom_p,atoms,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  MPI_Bcast(&V_inv->a[0][0],dim1ax*dim1ax,MPI_DOUBLE,0,MPI_COMM_WORLD);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("generate_coulomb_matrix_inverse            %10.2f\n",time4);

  // ******************************************************************************************
  // * Generate three centre integrals <ij|alpha> and <beta|V|kl> in core                     *
  // ******************************************************************************************
  
  time1 = MPI_Wtime();
  AllocateIntArray(&dim_send,&job->numtasks,job);
  AllocateIntArray(&dim_send1,&job->numtasks,job);
  AllocateIntArray(&dim_ham,&job->numtasks,job);
  AllocateIntArray(&dim_recv,&job->numtasks,job);
  AllocateIntArray(&offset_j,&dim1,job);
  if (job->taskid == 0) printf("NTRANSITIONS %3d\n",ntransitions);
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
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
  mpi_begin_end(begin_j, end_j, dim1, num_proc, job, file);
  for (i = num_proc; i < job->numtasks; i++) begin_j[i] = 0;
  for (i = num_proc; i < job->numtasks; i++) end_j[i] = 0;

  for (i = 0; i < dim1; i++) offset_j[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1] * nbands * nbands; }}

  for (i = 0; i < job->numtasks; i++) { dim_send[i] = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) {
  dim_send[i] += atoms_ax->bfnnumb_sh[j1] * nbands * nbands; }}
  AllocateDoubleArray(&integral_buffer,&dim_send[job->taskid],job);
  AllocateDoubleArray(&integral_buffer1,&dim_send[job->taskid],job);
  AllocateDoubleArray(&integral_buffer2,&dim_send[job->taskid],job);
  ResetDoubleArray(integral_buffer,&dim_send[job->taskid]);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
  integrals_occ_occ_vir_vir1(&j1,fh,&integral_buffer[offset_j[j1]],fermi,atom_p,atoms,shells,gaussians,atoms_ax, \
  shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  } // close loop on j1
  MPI_File_close(&fh);
  //count_arrays(&sendtotal,senddispls,sendcounts,&recvtotal,recvdispls,recvcounts,atoms_ax,fermi,job,file);
  MPI_Barrier(MPI_COMM_WORLD);

  initialise_ring_rotate(num_proc,&dim3,dim_send,dim_recv,job,file);
  for (i = 0; i < dim_send[job->taskid]; i++) integral_buffer1[i] = integral_buffer[i]; 
  time3 = MPI_Wtime();
  ResetDoubleArray(integral_buffer2,&dim_send[job->taskid]);
  for (i = 0; i < num_proc; i++) {
    jj = (job->taskid - i) % num_proc;
    kk = jj < 0 ? jj + num_proc : jj;
    //for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
      //for (j2 = begin_j[kk]; j2 < end_j[kk]; j2++) {
      for (j1 = begin_j[kk]; j1 < end_j[kk]; j1++) {
    for (j2 = begin_j[job->taskid]; j2 < end_j[job->taskid]; j2++) {
        contract_coulomb_integrals(&j1,&j2,V_inv,&integral_buffer[offset_j[j1]],&integral_buffer2[offset_j[j2]],fermi,atom_p,\
        atoms,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
       }
      }
     ring_rotate_once1(&integral_buffer,i,num_proc,&dim3,dim_send,dim_recv,job,file);
    } // close loop over i
     //AllocateDoubleArray(&integrals1,&recvtotal,job);
     //AllocateDoubleArray(&integrals2,&recvtotal,job);
     //MPI_Alltoallv(integral_buffer1,sendcounts,senddispls,MPI_DOUBLE,integrals1,recvcounts,recvdispls,MPI_DOUBLE,MPI_COMM_WORLD);
     //MPI_Alltoallv(integral_buffer2,sendcounts,senddispls,MPI_DOUBLE,integrals2,recvcounts,recvdispls,MPI_DOUBLE,MPI_COMM_WORLD);
     DestroyDoubleArray(&integral_buffer,&dim_send[job->taskid],job);
     time4 = MPI_Wtime() - time3;
     if (job->taskid == 0) printf("i loop                                     %10.2f\n",time4);

  // ******************************************************************************************
  // * Contract integrals <ij|alpha> and <beta|V|kl> to generate Hamiltonian                  *
  // ******************************************************************************************
  
     time3 = MPI_Wtime();
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
       ResetDoubleArray(Hamiltonian_buffer,&dim_ham[j]);

  switch (job->bse_tda) {

  case 0: // Full-matrix BSE or TDHF

 // time3 = MPI_Wtime();
 // if (job->bse_spin == 0) { // singlet excited state: unscreened ladders
 // hamiltonian_ij_ab(Ham_buffer1,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);
 //}
 // if (job->bse_spin == 1) { // triplet excited state: unscreened ladders
 // hamiltonian_ij_ab(Ham_buffer1,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);
 //}
 // time4 = MPI_Wtime() - time3;
 // if (job->taskid == 0) printf("hamiltonian_ij_ab                          %10.2f\n",time4);
 //
 // time3 = MPI_Wtime();
 // if (job->bse_ham == 1) { // Full-matrix BSE screened ladders
 // hamiltonian_screened_ij_ab(Ham_buffer1,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);
 //}
 //
 // time4 = MPI_Wtime() - time3;
 // if (job->taskid == 0) printf("hamiltonian_screened_ij_ab                 %10.2f\n",time4);
 //
 // time3 = MPI_Wtime();
 // if (job->bse_spin == 0) { // singlet excited state: unscreened ladders and rings
 // hamiltonian_ib_ja(Ham_buffer2,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);
 //}
 // if (job->bse_spin == 1) { // triplet excited state: unscreened ladders
 // hamiltonian_ib_ja(Ham_buffer2,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);
 //}
 // time4 = MPI_Wtime() - time3;
 // if (job->taskid == 0) printf("hamiltonian_ib_ja                          %10.2f\n",time4);
 //
 // time3 = MPI_Wtime();
 // if (job->bse_ham == 1)  { // Full matrix-BSE screened ladders
 // hamiltonian_screened_ib_ja(Ham_buffer2,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);
 //}
 // time4 = MPI_Wtime() - time3;
 // if (job->taskid == 0) printf("hamiltonian_screened_ib_ja                 %10.2f\n",time4);
 //
 // double *Ham_temp1, *Ham_temp2;
 // AllocateDoubleArray(&Ham_temp1,&nqA,job);
 // AllocateDoubleArray(&Ham_temp2,&nqA,job);
 // for (i = 0; i < mpA; i++) {
 //   for (j = 0; j < nqA; j++) Ham_temp1[j] = (Ham_buffer1[i * nqA + j] + Ham_buffer2[i * nqA + j]);
 //   for (j = 0; j < nqA; j++) Ham_temp2[j] = (Ham_buffer1[i * nqA + j] - Ham_buffer2[i * nqA + j]);
 //   for (j = 0; j < nqA; j++) Ham_buffer1[i * nqA + j] = Ham_temp1[j]; // A + B matrix
 //   for (j = 0; j < nqA; j++) Ham_buffer2[i * nqA + j] = Ham_temp2[j]; // A - B matrix
 //  } 
 // DestroyDoubleArray(&Ham_temp1,&nqA,job);
 // DestroyDoubleArray(&Ham_temp2,&nqA,job);
 //
 // if (job->bse_spin == 0) { // singlet excited state: rings
 // time3 = MPI_Wtime();
 // //hamiltonian_ia_jb(Ham_buffer1,four,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);
 // hamiltonian_ia_jb1(Ham_buffer1,four,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);
 // time4 = MPI_Wtime() - time3;
 // if (job->taskid == 0) printf("hamiltonian_ia_jb1                         %10.2f\n",time4);
 //}
 //
 // if (job->bse_ham == 0) {
 // casida_diagonal_molecule_mpi_spk(&nt,ictxt,nbsize,&nbands,fermi,Ham_buffer1,scf_eigenvalues,job,file);
 // casida_diagonal_molecule_mpi_spk(&nt,ictxt,nbsize,&nbands,fermi,Ham_buffer2,scf_eigenvalues,job,file);
 //}
 //
 // else if (job->bse_ham == 1) {
 // casida_diagonal_molecule_mpi_spk(&nt,ictxt,nbsize,&nbands,fermi,Ham_buffer1,GW_eigenvalues,job,file);
 // casida_diagonal_molecule_mpi_spk(&nt,ictxt,nbsize,&nbands,fermi,Ham_buffer2,GW_eigenvalues,job,file);
 //}

  time2 = MPI_Wtime() - time1;
  //if (job->taskid == 0) printf("Hamiltonian setup spk %10.2lf\n",time2);

  break;

  case 1: // Tamm-Dancoff Approximation: TDHF or BSE

  // ******************************************************************************************
  // * Generate inverse of Coulomb matrix                                                     *
  // ******************************************************************************************
  
  time3 = MPI_Wtime();
  if (job->bse_spin == 0) { // singlet excited state: unscreened ladders and rings
  time3 = MPI_Wtime();
  factor1 = two;
  hamiltonian_in_core_ia_jb(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  &factor1,Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("hamiltonian_in_core_ia_jb                  %10.2f\n",time4);
  time3 = MPI_Wtime();
  hamiltonian_in_core_ij_ab(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("hamiltonian_in_core_ij_ab                  %10.2f\n",time4);
 }
  if (job->bse_spin == 1) { // triplet excited state: unscreened ladders
  time3 = MPI_Wtime();
  hamiltonian_in_core_ij_ab(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
  Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("hamiltonian_in_core_ij_ab                  %10.2f\n",time4);
 }
  if (job->bse_ham == 1) { // TDA-BSE screened ladders
  time3 = MPI_Wtime();
  ////hamiltonian_screened_ij_ab(Ham_buffer1,V_inv,fermi,ictxt,nbsize,atoms_ax,job,file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("hamiltonian_in_core_screened_ij_ab         %10.2f\n",time4);
 }

//  if (job->bse_ham == 0) {
//  casida_diagonal_molecule_mpi_spk(&nt,ictxt,nbsize,&nbands,fermi,Ham_buffer1,scf_eigenvalues,job,file);
// }
//  else if (job->bse_ham == 1) {
//  casida_diagonal_molecule_mpi_spk(&nt,ictxt,nbsize,&nbands,fermi,Ham_buffer1,GW_eigenvalues,job,file);
// }

  break;

  case 2: // Random Phase Approximation: RPA

       factor1 = four;
       hamiltonian_in_core_ia_jb(begin_j,end_j,offset_j,nbsize,&mpa,&nqa,&nprow_nbsize,&npcol_nbsize,&nprow_myrow,&npcol_mycol, \
       &factor1,Hamiltonian_buffer,integral_buffer1,integral_buffer2,atoms_ax,fermi,job,file);

 } // close switch (job->bse_tda)
  //     for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
  //       nd6 = atoms_ax->bfnnumb_sh[j1];
  //       offset1 = offset_j[j1];
  //       for (il = 0; il < mpa; il++) {
  //         I1 = nprow_nbsize * (il / *nbsize) + il % *nbsize + nprow_myrow;
  //         //I1 = nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + row) % nprow) * *nbsize;
  //         m = I1 / nvir;
  //         n = I1  - m  * nvir + nocc;
  //         for (a1 = 0; a1 < nd6; a1++) {
  //           for (jl = 0; jl < nqa; jl++) {
  //             I2 = npcol_nbsize * (jl / *nbsize) + jl % *nbsize + npcol_mycol;
  //             //I2 = npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + col) % npcol) * *nbsize;
  //             r = I2 / nvir;
  //             s = I2  - r  * nvir + nocc;
  //             Hamiltonian_buffer[il + mpa * jl] += four * integral_buffer1[offset1 + m * nbands * nd6 + n * nd6 + a1] * \
  //                                                         integral_buffer2[offset1 + r * nbands * nd6 + s * nd6 + a1];
  //            }
  //           }
  //          }
  //         } // close loop on j1
            MPI_Reduce(Hamiltonian_buffer,Ham_buffer1,dim_ham[j],MPI_DOUBLE,MPI_SUM,j,MPI_COMM_WORLD);
            //MPI_Reduce(Hamiltonian_buffer,Hamiltonian,dim_ham[j],MPI_DOUBLE,MPI_SUM,j,MPI_COMM_WORLD);
            DestroyDoubleArray(&Hamiltonian_buffer,&dim_ham[j],job);
           } // close loop over j
            DestroyDoubleArray(&integral_buffer1,&dim_send[job->taskid],job);
            DestroyDoubleArray(&integral_buffer2,&dim_send[job->taskid],job);
  DestroyDoubleMatrix(&V_inv,job);
            //DestroyDoubleArray(&integrals1,&recvtotal,job);
            //DestroyDoubleArray(&integrals2,&recvtotal,job);
            DestroyIntArray(&dim_send,&job->numtasks,job);
            DestroyIntArray(&dim_ham,&job->numtasks,job);
            DestroyIntArray(&dim_recv,&job->numtasks,job);
            DestroyIntArray(&offset_j,&dim1,job);
            time4 = MPI_Wtime() - time3;
            if (job->taskid == 0) printf("j loop                                     %10.2f\n",time4);
            time2 = MPI_Wtime() - time1;
            if (job->taskid == 0) printf("hamiltonian_ia_jb4                         %10.2f\n",time2);

}

void count_arrays(int *sendtotal, int *senddispls, int *sendcounts, int *recvtotal, int *recvdispls, int *recvcounts, ATOM *atoms_ax, FERMI *fermi,  JOB_PARAM *job, FILES file)

{

int i, j, j1;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int begin_ax1[job->numtasks], end_ax1[job->numtasks];
int begin_ax2[job->numtasks], end_ax2[job->numtasks];
int begin_j[job->numtasks], end_j[job->numtasks];
int pointer;
int num_proc = job->numtasks < dim1 ? job->numtasks : dim1;
IntMatrix *array_sizes;

  for (i = 0; i < job->numtasks; i++) begin_j[i] = 0;
  for (i = 0; i < job->numtasks; i++) end_j[i] = 0;
  for (i = 0; i < job->numtasks; i++) begin_ax1[i] = 0;
  for (i = 0; i < job->numtasks; i++) end_ax1[i] = 0;
  mpi_begin_end(begin_j, end_j, dim1, num_proc, job, file);
  mpi_begin_end(begin_ax2, end_ax2, dim1ax, job->numtasks, job, file);
  AllocateIntMatrix(&array_sizes,&job->numtasks,&job->numtasks,job);

  *sendtotal = 0;
  *recvtotal = 0;
  senddispls[0] = 0;
  recvdispls[0] = 0;

  begin_ax1[0] = 0;
  end_ax1[0] = 0;
  for (j1 = begin_j[0]; j1 < end_j[0]; j1++) end_ax1[0] += atoms_ax->bfnnumb_sh[j1];
  for (i = 1; i < job->numtasks; i++) {
    begin_ax1[i] = end_ax1[i - 1];
    end_ax1[i] = end_ax1[i - 1];
    for (j1 = begin_j[i]; j1 < end_j[i]; j1++) end_ax1[i] += atoms_ax->bfnnumb_sh[j1]; 
   }
  for (i = 0; i < job->numtasks; i++) {
    fprintf(file.out,"begin end %3d %3d %3d   %3d %3d\n",i,begin_ax1[i],end_ax1[i],begin_ax2[i],end_ax2[i]);
   }

  i = 0; 
  j = 0;
  pointer = 0;
  ResetIntMatrix(array_sizes);
  do {
      if (end_ax1[i] < end_ax2[j]) {
      array_sizes->a[i][j] += end_ax1[i] - pointer;
      pointer = end_ax1[i];
      i++;
     }
      else if (end_ax1[i] >  end_ax2[j]) {
      array_sizes->a[i][j] += end_ax2[j] - pointer;
      pointer = end_ax2[j];
      j++;
     }
      else if (end_ax1[i] == end_ax2[j]) {
      array_sizes->a[i][j] += end_ax1[i] - pointer;
      pointer = end_ax1[i];
      i++;
      j++;
     }
     fprintf(file.out,"i j %3d %3d    %3d\n",i,j,pointer);
    } while (i < job->numtasks && j < job->numtasks);
  for (i = 0; i < job->numtasks; i++) {
    for (j = 0; j < job->numtasks; j++) {
      fprintf(file.out,"%3d %3d   %3d\n",i,j,array_sizes->a[i][j]);
     }
    }

  for (i = 0; i < job->numtasks; i++) {
    sendcounts[i] = array_sizes->a[job->taskid][i];
    recvcounts[i] = array_sizes->a[i][job->taskid];
   }

      for (i = 0; i < job->numtasks; i++) fprintf(file.out,"%3d %3d %3d\n",i,sendcounts[i],recvcounts[i]);
      DestroyIntMatrix(&array_sizes,job);
  
  //fprintf(file.out,"send recv %3d %3d\n", *sendtotal, *recvtotal);
  //for (i = 0; i < job->numtasks; i++) { sendcounts[i] = 0;
  //for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { sendcounts[i] += atoms_ax->bfnnumb_sh[j1] * nbands * nbands; } }
  for (i = 0; i < job->numtasks - 1; i++) { senddispls[i+1] = sendcounts[i] + *sendtotal; *sendtotal += sendcounts[i]; }
  *sendtotal += sendcounts[job->numtasks - 1];
  //for (i = 0; i < job->numtasks; i++) recvcounts[i] = (end_ax2[i] - begin_ax2[i]) * nbands * nbands;
  for (i = 0; i < job->numtasks - 1; i++) { recvdispls[i+1] = recvcounts[i] + *recvtotal; *recvtotal += recvcounts[i]; }
  *recvtotal += recvcounts[job->numtasks - 1];
  fprintf(file.out,"send recv %3d %3d\n", *sendtotal, *recvtotal);
  for (i = 0; i < job->numtasks; i++) \
  fprintf(file.out,"%3d  %3d %3d  %3d %3d\n",i,sendcounts[i],recvcounts[i],senddispls[i],recvdispls[i]);

}

void ring_size_transfer1(double *buffer1, double *buffer2, int dim1, int dim2, int num_procs, JOB_PARAM *job, FILES file)

{   

int i, j;
double time1, time2;

  time1 = MPI_Wtime();
  if (job->taskid != 0 && job->taskid < num_procs) {
  MPI_Recv(buffer2, dim2, MPI_DOUBLE, job->taskid - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 } 
  if (job->taskid < num_procs) {
  MPI_Send(buffer1, dim1, MPI_DOUBLE, (job->taskid + 1) % num_procs, 0, MPI_COMM_WORLD);
 }
  if (job->taskid == 0) {
  MPI_Recv(buffer2, dim2, MPI_DOUBLE, num_procs - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 }
  time2 = MPI_Wtime() - time1;

}

void ring_rotate_once1(double **send_buffer, int k, int num_procs, int *dim3, int *dims, int *dimr, JOB_PARAM *job, FILES file)

{   

int i, j, r;
int dim4;
double time1, time2;
double *recv_buffer;

  time1 = MPI_Wtime();
  if (job->taskid < num_procs) {
  r = (job->taskid - k) % num_procs;
  dim4 = dimr[r < 0 ? r + num_procs : r];
  AllocateDoubleArray(&recv_buffer,&dim4,job);
  if (job->numtasks == 1)
  for (i = 0; i < dim4; i++) { recv_buffer[i] = (*send_buffer)[i];  }
  else if (job->numtasks > 1)
  ring_size_transfer1(*send_buffer,recv_buffer,*dim3,dim4,num_procs,job,file);
  DestroyDoubleArray(send_buffer,dim3,job);
  *dim3 = dim4;
  AllocateDoubleArray(send_buffer,dim3,job);
  for (i = 0; i < dim4; i++) { (*send_buffer)[i] = recv_buffer[i];  }
  //memcpy(*send_buffer, recv_buffer, dim4 * sizeof(double));
  DestroyDoubleArray(&recv_buffer,&dim4,job);
 }
  time2 += MPI_Wtime() - time1;

}

void initialise_ring_rotate(int num_procs, int *dim3, int *dims, int *dimr, JOB_PARAM *job, FILES file)

{

int i;

  for (i = 1; i < num_procs; i++) dimr[i] = dims[i - 1];
  dimr[0] = dims[num_procs - 1];
 *dim3 = dims[job->taskid];

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

*/
/*

void hamiltonian_in_core_screened_ij_ab(int *j, int *dim_ham, int *begin_j, int *end_j, int *offset_j, int *nbsize, int *mpA, int *nqA, int *nprow_nbsize, int *npcol_nbsize, int* nprow_myrow, int *npcol_mycol, double *Ham_buffer1, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, l, a1, il, jl, j1, I1, I2, m, n, r, s;
int nd6, offset, offset1, offset2;
int *offset_j1;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int dim4, mdim;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
double time1, time2, time3, time4;
double *temp4;
double *m_1, *m_2, *m_buffer1, *m_buffer2;
double *cas_eigenvalues;
FILE *cas_evals;
char zz4[24] = "cas_evalues";

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
  //printf("dim4 %3d\n",dim4);

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("Temp4                                      %10.2f\n",time4);

  AllocateIntArray(&offset_j1,&dim1,job);
  for (i = 0; i < dim1; i++) offset_j1[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j1[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1]; }}

  // ******************************************************************************************
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();
  read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz4, job, file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("Read                                       %10.2f\n",time4);

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with <alpha|occ-occ>                             *
  // ******************************************************************************************
  
  mdim = *mpA * *nqA * job->bse_lim;
  AllocateDoubleArray(&m_buffer1,&mdim,job);
  AllocateDoubleArray(&m_buffer2,&mdim,job);
  ResetDoubleArray(m_buffer1,&mdim);
  ResetDoubleArray(m_buffer2,&mdim);

  time3 = MPI_Wtime();
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l = 0; l < job->bse_lim; l++) {
      for (il = 0; il < *mpA; il++) {
        I1 = *nprow_nbsize * (il / *nbsize) + il % *nbsize + *nprow_myrow;
        m = I1 / nvir;
        n = I1  - m  * nvir + nocc;
        for (jl = 0; jl < *nqA; jl++) {
          I2 = *npcol_nbsize * (jl / *nbsize) + jl % *nbsize + *npcol_mycol;
          r = I2 / nvir;
          s = I2  - r  * nvir + nocc;
          for (a1 = 0; a1 < nd6; a1++) {
            m_buffer1[l * *mpA * *nqA + il * *nqA + jl] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[l * nd6 + offset1 + a1];
           }
          }
         }
        }
       }
//for (i = 0; i < mdim; i++)
//for (il = 0; il < *mpA; il++) {
//for (jl = 0; jl < *nqA; jl++) {
//for (l = 0; l < job->bse_lim; l++) {
//fprintf(file.out,"m_buffer1 %3d %3d %3d %10.6f\n",il, jl, l, m_buffer1[l * *mpA * *nqA + il * *nqA + jl]); }}}
//fprintf(file.out,"\n");

  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l = 0; l < job->bse_lim; l++) {
      for (il = 0; il < *mpA; il++) {
        I1 = *nprow_nbsize * (il / *nbsize) + il % *nbsize + *nprow_myrow;
        m = I1 / nvir;
        n = I1  - m  * nvir + nocc;
        for (jl = 0; jl < *nqA; jl++) {
          I2 = *npcol_nbsize * (jl / *nbsize) + jl % *nbsize + *npcol_mycol;
          r = I2 / nvir;
          s = I2  - r  * nvir + nocc;
          for (a1 = 0; a1 < nd6; a1++) {
            m_buffer2[l * *mpA * *nqA + il * *nqA + jl] += integral_buffer1[offset2 + n * nbands * nd6 + s * nd6 + a1] * \
            temp4[l * nd6 + offset1 + a1];
           }
          }
         }
        }
       }
      time4 = MPI_Wtime() - time3;
      if (job->taskid == 0) printf("Buffers                                    %10.2f\n",time4);

//for (il = 0; il < *mpA; il++) {
//for (jl = 0; jl < *nqA; jl++) {
//for (l = 0; l < job->bse_lim; l++) {
//fprintf(file.out,"m_buffer2 %3d %3d %3d %10.6f\n",il, jl, l, m_buffer2[l * *mpA * *nqA + il * *nqA + jl]); }}}
//fprintf(file.out,"\n");

  AllocateDoubleArray(&m_1,&mdim,job);
  AllocateDoubleArray(&m_2,&mdim,job);
  ResetDoubleArray(m_1,&mdim);
  ResetDoubleArray(m_2,&mdim);

  time3 = MPI_Wtime();
  MPI_Reduce(m_buffer1,m_1,job->bse_lim * dim_ham[*j],MPI_DOUBLE,MPI_SUM,*j,MPI_COMM_WORLD);
  MPI_Reduce(m_buffer2,m_2,job->bse_lim * dim_ham[*j],MPI_DOUBLE,MPI_SUM,*j,MPI_COMM_WORLD);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("MPI_Reduce                                 %10.2f\n",time4);

  DestroyDoubleArray(&m_buffer1,&mdim,job);
  DestroyDoubleArray(&m_buffer2,&mdim,job);

double *Ham_buffer2;
int mpanqa = *mpA * *nqA;
AllocateDoubleArray(&Ham_buffer2,&mpanqa,job);
ResetDoubleArray(Ham_buffer2,&mpanqa);

  if (job->taskid == *j) {
    time3 = MPI_Wtime();
    for (l = 0; l < job->bse_lim; l++) {
      for (il = 0; il < *mpA; il++) {
        for (jl = 0; jl < *nqA; jl++) {
          Ham_buffer1[il + *mpA * jl] += two * m_1[l * *mpA * *nqA + il * *nqA + jl] * m_2[l * *mpA * *nqA + il * *nqA + jl] / \
          cas_eigenvalues[l];
          Ham_buffer2[il + *mpA * jl] += two * m_1[l * *mpA * *nqA + il * *nqA + jl] * m_2[l * *mpA * *nqA + il * *nqA + jl] / \
          cas_eigenvalues[l];
          fprintf(file.out,"%3d %3d %3d %3d %8.4f %8.4f %8.4f\n",*j,il,jl,l,Ham_buffer2[il + *mpA * jl],\
          m_1[l * *mpA * *nqA + il * *nqA + jl], m_2[l * *mpA * *nqA + il * *nqA + jl]);
         }
        }
       }
      time4 = MPI_Wtime() - time3;
      if (job->taskid == 0) printf("Contract                                   %10.2f\n",time4);
      }

printf("%3d\n",mpanqa);
for (i = 0; i < mpanqa; i++)
fprintf(file.out,"Ham2 old %10.6f\n",Ham_buffer2[i]);

  DestroyDoubleArray(&m_1,&mdim,job);
  DestroyDoubleArray(&m_2,&mdim,job);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&temp4,&dim4,job);

}

void hamiltonian_in_core_screened_ij_ab1(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *Ham_buffer1, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, l, a1, il, jl, j1, I1, I2, m, n, r, s;
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
double time1, time2, time3, time4, time5, time6;
double *temp4;
double *M_1, *M_2, *M_buffer1, *M_buffer2;
double *cas_eigenvalues;
FILE *cas_evals;
char zz4[24] = "cas_evalues";

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
  time4 = MPI_Wtime() - time3;
  //if (job->taskid == 0) printf("Temp4                                      %10.2f\n",time4);

  AllocateIntArray(&offset_j1,&dim1,job);
  for (i = 0; i < dim1; i++) offset_j1[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j1[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1]; }}

  // ******************************************************************************************
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();
  read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz4, job, file);
  time4 = MPI_Wtime() - time3;
  //if (job->taskid == 0) printf("Read                                       %10.2f\n",time4);

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with <alpha|occ-occ>                             *
  // ******************************************************************************************
  
  time2 = k_zero;
  time4 = k_zero;
  time6 = k_zero;

  //double *Ham_buffer2;
  //int mpanqa = mpA * nqA;
  //int mpanqa = *mpA * *nqA;
  //AllocateDoubleArray(&Ham_buffer2,&mpanqa,job);
  //ResetDoubleArray(Ham_buffer2,&mpanqa);
  //mpa = 6; nqa = 2; myrow = 0; mycol = job->taskid; nprow = 1; npcol = 3;
  //printf("%3d %3d %3d %3d %3d %3d %3d\n",job->taskid,mpa,nqa,myrow,mycol,nprow,npcol);
  //printf("itr %3d %3d %3d %3d\n",itr,nvir,mpa,nqa);

  for (l1 = 0; l1 < itr; l1++) {

  time1 = MPI_Wtime();
  ////mdim1 = nocc * nocc * nvir * nocc;
  mdim1 = nocc * nocc * nvir;
  AllocateDoubleArray(&M_buffer1,&mdim1,job);
  ResetDoubleArray(M_buffer1,&mdim1);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    ////for (l2 = 0; l2 < nocc * nvir; l2++) {
    for (l2 = 0; l2 < nvir; l2++) {
      for (m = 0; m < nocc; m++) {
        for (r = 0; r < nocc; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
        //M_buffer1[m * nocc * nvir * nocc + r * nvir * nocc + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
            M_buffer1[m * nocc * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
           }
          }
         }
        }
       }
      //for (i = 0; i < mdim1; i++)
      //fprintf(file.out,"M_buffer1 %3d %10.6f\n",l1 * nvir + l2,M_buffer1[i]);
      //fprintf(file.out,"\n");

  AllocateDoubleArray(&M_1,&mdim1,job);
  ResetDoubleArray(M_1,&mdim1);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer1,M_1,mdim1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  DestroyDoubleArray(&M_buffer1,&mdim1,job);

  mdim2 = nvir * nvir * nvir;
  ////mdim2 = nvir * nvir * nvir * nocc;
  AllocateDoubleArray(&M_buffer2,&mdim2,job);
  ResetDoubleArray(M_buffer2,&mdim2);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    ////for (l2 = 0; l2 < nocc * nvir; l2++) {
    for (l2 = 0; l2 < nvir; l2++) {
      for (n = 0; n < nvir; n++) {
        for (s = 0; s < nvir; s++) {
          for (a1 = 0; a1 < nd6; a1++) {
            ////M_buffer2[n * nvir * nvir * nocc + s * nvir * nocc + l2] += \
            integral_buffer1[offset2 + (nocc + n) * nbands * nd6 + (nocc + s) * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1];
            M_buffer2[n * nvir * nvir + s * nvir + l2] += \
            integral_buffer1[offset2 + (nocc + n) * nbands * nd6 + (nocc + s) * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1];
           }
          }
         }
        }
       }
    //for (i = 0; i < mdim2; i++)
    //fprintf(file.out,"M_buffer2 %3d %10.6f\n",l1 * nvir + l2,M_buffer2[i]);
    //fprintf(file.out,"\n");

  AllocateDoubleArray(&M_2,&mdim2,job);
  ResetDoubleArray(M_2,&mdim2);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer2,M_2,mdim2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  DestroyDoubleArray(&M_buffer2,&mdim2,job);
  time2 += MPI_Wtime() - time1;

  time5 = MPI_Wtime();
  for (il = 0; il < mpA; il++) {
  //for (il = 0; il < *mpA; il++) {
    I1 = nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize;
    //I1 = *nprow_nbsize * (il / *nbsize) + il % *nbsize + *nprow_myrow;
    m = I1 / nvir;
    n = I1  - m  * nvir;
    ////for (l2 = 0; l2 < nocc * nvir; l2++) {
    for (jl = 0; jl < nqA; jl++) {
    //for (jl = 0; jl < *nqA; jl++) {
      I2 = npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + mycol) % npcol) * *nbsize;
      //I2 = *npcol_nbsize * (jl / *nbsize) + jl % *nbsize + *npcol_mycol;
      r = I2 / nvir;
      s = I2  - r  * nvir;
      //fprintf(file.out,"I1 I2 %3d %3d  mnrs  %3d %3d %3d %3d  mrns %3d %3d %3d %3d\n",I1,I2,m,n,r,s,m,r,n,s);
      for (l2 = 0; l2 < nvir; l2++) {
      ////Ham_buffer1[il + *mpA * jl] += two * M_1[m * nocc * nvir * nocc + r * nvir * nocc + l2] * \
                                           M_2[n * nvir * nvir * nocc + s * nvir * nocc + l2] / \
                                           cas_eigenvalues[l1 * nvir + l2];
        //Ham_buffer1[il + *mpA * jl] += two * M_1[m * nocc * nvir + r * nvir + l2] * M_2[n * nvir * nvir + s * nvir + l2] / \
        cas_eigenvalues[l1 * nvir + l2];
        Ham_buffer1[il + mpA * jl] += two * M_1[m * nocc * nvir + r * nvir + l2] * M_2[n * nvir * nvir + s * nvir + l2] / \
        cas_eigenvalues[l1 * nvir + l2];
        //Ham_buffer2[il + *mpA * jl] += two * M_1[m * nocc  * nvir + r * nvir + l2] * M_2[n * nvir  * nvir + s * nvir + l2] / \
        cas_eigenvalues[l1 * nvir + l2];
        //Ham_buffer2[il + mpA * jl] += two * M_1[m * nocc  * nvir + r * nvir + l2] * M_2[n * nvir * nvir + s * nvir + l2] / \
        cas_eigenvalues[l1 * nvir + l2];
        //Ham_buffer2[il + *mpA * jl] += two * M_1[m * nocc * nvir * nocc + r * nvir * nocc + l2] * \
                                             M_2[n * nvir * nvir * nocc + s * nvir * nocc + l2] / \
                                             cas_eigenvalues[l1 * nvir + l2];
        //fprintf(file.out,"ab1 %3d %3d %3d %3d   %3d %3d %3d %3d %8.4f %8.4f %8.4f %8.4f\n",j1,il,jl,l2,m,r,nocc+n,nocc+s,\
        Ham_buffer2[il + *mpA * jl],M_1[m * nocc * nvir * nocc + r * nvir * nocc + l2], \
        M_2[n * nvir * nvir * nocc + s * nvir * nocc + l2], cas_eigenvalues[l1 + l2]);
       }
      }
     }
  //fprintf(file.out,"\n");
  time6 += MPI_Wtime() - time5;
  //for (i = 0; i < mdim1; i++)
  //fprintf(file.out,"M_1 %3d %10.6f\n",l1 * nvir + l2,M_1[i]);
  //fprintf(file.out,"\n");
  //for (i = 0; i < mdim2; i++)
  //fprintf(file.out,"M_2 %3d %10.6f\n",l1 * nvir + l2,M_2[i]);
  //fprintf(file.out,"\n");
  DestroyDoubleArray(&M_1,&mdim1,job);
  DestroyDoubleArray(&M_2,&mdim1,job);
 } // close loop on l1

  //for (i = 0; i < mpA * nqA; i++)
  //fprintf(file.out,"Ham2 %10.6f\n",Ham_buffer2[i]);
  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Contract                                   %10.2f\n",time6);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&temp4,&dim4,job);

}
*/

/*
void hamiltonian_in_core_screened_ij_ab3(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *GW_eigenvalues, double *Ham_buffer1, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

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
char zz4[24] = "scf_evalues";
char zz5[24] = "cas_evalues";

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;
  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
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
  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();
  read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);
  read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz5, job, file);
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

////CHP
//int kk  = 0;
////CHP
  time1 = MPI_Wtime();
  mdim = nbands * nbands * nvir;
  //fprintf(file.out,"hamiltonian_in_core_screened_ij_ab3  allocating nbands^2 %7d * nvir %7d = %7d\n", \
  nbands*nbands,nvir, nbands*nbands*nvir); fflush(file.out);
  AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l2 = 0; l2 < nvir; l2++) {
      for (m = 0; m < nbands; m++) {
        for (r = 0; r < nbands; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer[m * nbands * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]; 
////CHP
      //      if (kk < 1000 && (M_buffer[m * nbands * nvir + r * nvir + l2] != M_buffer[m * nbands * nvir + r * nvir + l2] || \
      //          integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] !=  integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] || \
      //          temp4[(l1 * nvir + l2) * nd6 + offset1 + a1] != temp4[(l1 * nvir + l2) * nd6 + offset1 + a1])){
      //          fprintf(file.out,"%3d %3d %3d %3d %3d %f %f %f\n",j1,l2,m,r,a1,M_buffer[m * nbands * nvir + r * nvir + l2],\
      //          integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1],  temp4[(l1 * nvir + l2) * nd6 + offset1 + a1]);
      //          kk++;
      //     }
      //
           }
          }
         }
        }
       }
  time2 += MPI_Wtime() - time1;

  AllocateDoubleArray(&M,&mdim,job);
  //fprintf(file.out,"hamiltonian_in_core_screened_ij_ab3  allocating mdim %7d\n", mdim); fflush(file.out);
  ResetDoubleArray(M,&mdim);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  DestroyDoubleArray(&M_buffer,&mdim,job);

////CHP
//int ii = 0;

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
  //      if (M[i * nbands * nvir + k * nvir + l2] != M[i * nbands * nvir + k * nvir + l2] || \
  //          M[k * nbands * nvir + i * nvir + l2] != M[k * nbands * nvir + i * nvir + l2] || 
  //          cas_eigenvalues[l1 * nvir + l2] != cas_eigenvalues[l1 * nvir + l2]) {
  //          sigma_factor = 0.0;
  //          if (ii < 1000) {
  //          fprintf(file.out,"i %3d k %3d l2 %3d l1 %3d Mik %f Mki %f cas_eval %f\n",\
  //          i,k,l2,l1,M[i * nbands * nvir + k * nvir + l2],M[k * nbands * nvir + i * nvir + l2],cas_eigenvalues[l1 * nvir + l2]);
  //          ii++;
  //       }
  //      }
////CHP
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }

////CHP
//int jj = 0;

  for (i = 0; i < nbands; i++) {
    for (k = nocc; k < nbands; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] - cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;
 //       if (M[i * nbands * nvir + k * nvir + l2] != M[i * nbands * nvir + k * nvir + l2] || \
 //           M[k * nbands * nvir + i * nvir + l2] != M[k * nbands * nvir + i * nvir + l2] || 
 //           cas_eigenvalues[l1 * nvir + l2] != cas_eigenvalues[l1 * nvir + l2]) {
 //           sigma_factor = 0.0;
 //           if (jj < 1000) {
 //           fprintf(file.out,"i %3d k %3d l2 %3d l1 %3d Mik %f Mki %f cas_eval %f\n",\
 //           i,k,l2,l1,M[i * nbands * nvir + k * nvir + l2],M[k * nbands * nvir + i * nvir + l2],cas_eigenvalues[l1 * nvir + l2]);
 //           jj++;
 //        }
 //       }
////CHP
        Sigma[i] += sigma_factor;
        dSigma_dE[i] -= sigma_factor / denom;
       }
      }
     }
  //for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);
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

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);
  if (job->taskid == 0 ) printf("Contract                                   %10.2f\n",time8);
  if (job->taskid == 0) {
  fprintf(file.out,"GW Eigenvalues: %4d RPA states used\n",job->bse_lim);
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
  fflush(file.out);
 }

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
char zz4[24] = "scf_evalues";
char zz5[24] = "cas_evalues";

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
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
  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();
  read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);
  read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz5, job, file);
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
  //mdim = nbands * nbands * nvir;
  //AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  ResetDoubleArray(M,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
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

  //AllocateDoubleArray(&M,&mdim,job);
  //ResetDoubleArray(M,&mdim);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  //DestroyDoubleArray(&M_buffer,&mdim,job);

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;

        //fprintf(file.out,"%3d %3d %3d %10.3e %10.3e %10.3e %10.3e\n",i,k,l2,\
        sigma_factor,M[i * nbands * nvir + k * nvir + l2], M[k * nbands * nvir + i * nvir + l2], denom);

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
  //for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);
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

  //DestroyDoubleArray(&M,&mdim,job);

 } // close loop on l1
  DestroyDoubleArray(&M_buffer,&mdim,job);
  DestroyDoubleArray(&M,&mdim,job);

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);
  if (job->taskid == 0 ) printf("Contract                                   %10.2f\n",time8);
  if (job->taskid == 0) {
  fprintf(file.out,"GW Eigenvalues: %4d RPA states used\n",job->bse_lim);
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
  fflush(file.out);
 }

  DestroyDoubleArray(&scf_eigenvalues,&nbands,job);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&temp4,&dim4,job);

}
*/

/*
void self_energy_diagonal_in_core(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

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
double *cas_eigenvalues, *scf_eigenvalues, *GW_eigenvalues;
double Sigma[nbands], dSigma_dE[nbands], sigma_factor, denom;
char zz4[24] = "scf_evalues";
char zz5[24] = "cas_evalues";

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
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
  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();
  read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);
  read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz5, job, file);
  time4 = MPI_Wtime() - time3;

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with <alpha|occ-occ>                             *
  // ******************************************************************************************
  
  time2 = k_zero;
  time4 = k_zero;
  time6 = k_zero;

  for (i = 0; i < nbands; i++) Sigma[i] = k_zero;
  for (i = 0; i < nbands; i++) dSigma_dE[i] = k_zero;

  for (l1 = 0; l1 < itr; l1++) {

  time1 = MPI_Wtime();
  mdim = nbands * nbands * nvir;
  AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l2 = 0; l2 < nvir; l2++) {
      for (m = 0; m < nbands; m++) {
        for (r = 0; r < nbands; r++) {
          for (a1 = 0; a1 < nd6; a1++) {
            M_buffer[m * nbands * nvir + r * nvir + l2] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1]; // * \
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
  //for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);
  time6 += MPI_Wtime() - time5;

  DestroyDoubleArray(&M,&mdim,job);

 } // close loop on l1

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);
  if (job->taskid == 0) {
  fprintf(file.out,"GW Eigenvalues: %4d RPA states used\n",job->bse_lim);
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
  fflush(file.out);
  char zz4[24] = "GW_evalues";
  FILE *evals = fopen(zz4, "wb");
  fseek(evals,0,SEEK_SET);
  fwrite(GW_eigenvalues,sizeof(double),job->spin_dim * nbands,evals);
  fclose(evals);
 }
  DestroyDoubleArray(&GW_eigenvalues,&nbands,job);
  DestroyDoubleArray(&scf_eigenvalues,&nbands,job);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&temp4,&dim4,job);

}

void self_energy_diagonal_in_core1(int *ictxt, int *nbsize, int *begin_j, int *end_j, int *offset_j, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

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
double *cas_eigenvalues, *scf_eigenvalues, *GW_eigenvalues;
double Sigma[nbands], dSigma_dE[nbands], sigma_factor, denom;
char zz4[24] = "scf_evalues";
char zz5[24] = "cas_evalues";

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
  itr = job->bse_lim / nvir;
  rem = job->bse_lim - itr * nvir;

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
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
  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();
  read_scf_GW_eigenvalues(scf_eigenvalues, fermi->bands[0] - 1, nbands, zz4, job, file);
  read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz5, job, file);
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
  //mdim = nbands * nbands * nvir;
  //AllocateDoubleArray(&M_buffer,&mdim,job);
  ResetDoubleArray(M_buffer,&mdim);
  ResetDoubleArray(M,&mdim);
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
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

  //AllocateDoubleArray(&M,&mdim,job);
  //ResetDoubleArray(M,&mdim);
  time3 = MPI_Wtime();
  MPI_Allreduce(M_buffer,M,mdim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  time4 += MPI_Wtime() - time3;
  //DestroyDoubleArray(&M_buffer,&mdim,job);

  time5 = MPI_Wtime();
  for (i = 0; i < nbands; i++) {
    for (k = 0; k < nocc; k++) {
      for (l2 = 0; l2 < nvir; l2++) {
        denom = scf_eigenvalues[i] - scf_eigenvalues[k] + cas_eigenvalues[l1 * nvir + l2];
        sigma_factor = two * M[i * nbands * nvir + k * nvir + l2] * M[k * nbands * nvir + i * nvir + l2] / denom;

        //fprintf(file.out,"%3d %3d %3d %10.3e %10.3e %10.3e %10.3e\n",i,k,l2,\
        sigma_factor,M[i * nbands * nvir + k * nvir + l2], M[k * nbands * nvir + i * nvir + l2], denom);

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
  //for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);
  time6 += MPI_Wtime() - time5;

 } // close loop on l1
  DestroyDoubleArray(&M_buffer,&mdim,job);
  DestroyDoubleArray(&M,&mdim,job);

  for (i = 0; i < nbands; i++) GW_eigenvalues[i] = scf_eigenvalues[i] + Sigma[i] / (k_one - dSigma_dE[i]);

  if (job->taskid == 0)  printf("Buffers                                    %10.2f\n",time2);
  if (job->taskid == 0)  printf("MPI_Allreduce                              %10.2f\n",time4);
  if (job->taskid == 0 ) printf("Sigma                                      %10.2f\n",time6);
  if (job->taskid == 0) {
  fprintf(file.out,"GW Eigenvalues: %4d RPA states used\n",job->bse_lim);
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
  fflush(file.out);
  char zz4[24] = "GW_evalues";
  FILE *evals = fopen(zz4, "wb");
  fseek(evals,0,SEEK_SET);
  fwrite(GW_eigenvalues,sizeof(double),job->spin_dim * nbands,evals);
  fclose(evals);
 }

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

void hamiltonian_in_core_screened_ij_ab2(int *j, int *dim_ham, int *begin_j, int *end_j, int *offset_j, int *nbsize, int *mpA, int *nqA, int *nprow_nbsize, int *npcol_nbsize, int* nprow_myrow, int *npcol_mycol, double *Ham_buffer1, double *integral_buffer1, double *integral_buffer2, ATOM *atoms_ax, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, l, a1, il, jl, j1, I1, I2, m, n, r, s;
int nd6, offset, offset1, offset2;
int *offset_j1;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int dim4, mdim;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int nocc = fermi->occupied[0];
int nvir = nbands - fermi->occupied[0];
int ntransitions = nocc * nvir;
double time1, time2, time3, time4;
double *temp4;
double *m_1, *m_2, *m_buffer1, *m_buffer2;
double *cas_eigenvalues;
FILE *cas_evals;
char zz4[24] = "cas_evalues";

  // ******************************************************************************************
  // * Generate temp4 array                                                                   *
  // ******************************************************************************************
  
  if (job->bse_lim == 0 || job->bse_lim > ntransitions) job->bse_lim = ntransitions;
  if (job->taskid == 0) printf("BSE vector limit %3d\n",job->bse_lim);
  dim4 = 0;
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) { dim4 += job->bse_lim * atoms_ax->bfnnumb_sh[j1]; }
printf("dim4 %3d\n",dim4);

  if (dim4 == 0) dim4 = 1;
  AllocateDoubleArray(&temp4,&dim4,job);
  ResetDoubleArray(temp4,&dim4);
  time3 = MPI_Wtime();
  generate_temp4_in_core(temp4,integral_buffer2,fermi,atoms_ax,job,file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("Temp4                                      %10.2f\n",time4);

  AllocateIntArray(&offset_j1,&dim1,job);
  for (i = 0; i < dim1; i++) offset_j1[i] = 0;
  for (i = 0; i < job->numtasks; i++) { offset = 0;
  for (j1 = begin_j[i]; j1 < end_j[i]; j1++) { offset_j1[j1] = offset;
  offset += atoms_ax->bfnnumb_sh[j1]; }}

  // ******************************************************************************************
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  AllocateDoubleArray(&cas_eigenvalues,&job->bse_lim,job);
  ResetDoubleArray(cas_eigenvalues,&job->bse_lim);
  time3 = MPI_Wtime();
  read_scf_GW_eigenvalues(cas_eigenvalues, 0, job->bse_lim, zz4, job, file);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("Read                                       %10.2f\n",time4);

  // ******************************************************************************************
  // * Contract V <alpha|beta><beta|occ-vir> with <alpha|occ-occ>                             *
  // ******************************************************************************************
  
  mdim = *mpA * *nqA * job->bse_lim;
  AllocateDoubleArray(&m_buffer1,&mdim,job);
  AllocateDoubleArray(&m_buffer2,&mdim,job);
  ResetDoubleArray(m_buffer1,&mdim);
  ResetDoubleArray(m_buffer2,&mdim);

  time3 = MPI_Wtime();
  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l = 0; l < job->bse_lim; l++) {
      for (il = 0; il < *mpA; il++) {
        I1 = *nprow_nbsize * (il / *nbsize) + il % *nbsize + *nprow_myrow;
        m = I1 / nvir;
        n = I1  - m  * nvir + nocc;
        for (jl = 0; jl < *nqA; jl++) {
          I2 = *npcol_nbsize * (jl / *nbsize) + jl % *nbsize + *npcol_mycol;
          r = I2 / nvir;
          s = I2  - r  * nvir + nocc;
          for (a1 = 0; a1 < nd6; a1++) {
            m_buffer1[l * *mpA * *nqA + il * *nqA + jl] += integral_buffer1[offset2 + m * nbands * nd6 + r * nd6 + a1] * \
            temp4[l * nd6 + offset1 + a1];
           }
          }
         }
        }
       }

  for (j1 = begin_j[job->taskid]; j1 < end_j[job->taskid]; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    offset1 = job->bse_lim * offset_j1[j1];
    offset2 = offset_j[j1];
    for (l = 0; l < job->bse_lim; l++) {
      for (il = 0; il < *mpA; il++) {
        I1 = *nprow_nbsize * (il / *nbsize) + il % *nbsize + *nprow_myrow;
        m = I1 / nvir;
        n = I1  - m  * nvir + nocc;
        for (jl = 0; jl < *nqA; jl++) {
          I2 = *npcol_nbsize * (jl / *nbsize) + jl % *nbsize + *npcol_mycol;
          r = I2 / nvir;
          s = I2  - r  * nvir + nocc;
          for (a1 = 0; a1 < nd6; a1++) {
            m_buffer2[l * *mpA * *nqA + il * *nqA + jl] += integral_buffer1[offset2 + n * nbands * nd6 + s * nd6 + a1] * \
            temp4[l * nd6 + offset1 + a1];
           }
          }
         }
        }
       }
      time4 = MPI_Wtime() - time3;
      if (job->taskid == 0) printf("Buffers                                    %10.2f\n",time4);

  AllocateDoubleArray(&m_1,&mdim,job);
  AllocateDoubleArray(&m_2,&mdim,job);
  ResetDoubleArray(m_1,&mdim);
  ResetDoubleArray(m_2,&mdim);

  time3 = MPI_Wtime();
  MPI_Reduce(m_buffer1,m_1,job->bse_lim * dim_ham[*j],MPI_DOUBLE,MPI_SUM,*j,MPI_COMM_WORLD);
  MPI_Reduce(m_buffer2,m_2,job->bse_lim * dim_ham[*j],MPI_DOUBLE,MPI_SUM,*j,MPI_COMM_WORLD);
  time4 = MPI_Wtime() - time3;
  if (job->taskid == 0) printf("MPI_Reduce                                 %10.2f\n",time4);

  DestroyDoubleArray(&m_buffer1,&mdim,job);
  DestroyDoubleArray(&m_buffer2,&mdim,job);

  if (job->taskid == *j) {
    time3 = MPI_Wtime();
    for (l = 0; l < job->bse_lim; l++) {
      for (il = 0; il < *mpA; il++) {
        for (jl = 0; jl < *nqA; jl++) {
          Ham_buffer1[il + *mpA * jl] += two * m_1[l * *mpA * *nqA + il * *nqA + jl] * m_2[l * *mpA * *nqA + il * *nqA + jl] /\
          cas_eigenvalues[l];
         }
        }
       }
      time4 = MPI_Wtime() - time3;
      if (job->taskid == 0) printf("Contract                                   %10.2f\n",time4);
      }

  DestroyDoubleArray(&m_1,&mdim,job);
  DestroyDoubleArray(&m_2,&mdim,job);
  DestroyDoubleArray(&cas_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&temp4,&dim4,job);

}

void generate_coulomb_matrix_inverse(DoubleMatrix *V_inv, ATOM_TRAN *atom_p, ATOM *atoms, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k;
int Function[8];
int dim, dimg, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int info = 0;
double *eigenvalues;
double time1, time2;
double alpha = k_one, beta = k_zero;
char uplo = 'U';
char jobz = 'V';
char NoTrans = 'N', Trans = 'T';
DoubleMatrix *V, *xtrn1, *eigenvectors;
PAIR_TRAN pair_p;
INT_1E one_ints;

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ; // needed for screened potential 
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 1 ;

  time1 = MPI_Wtime();
  count_density_pairs4(&pair_p, atoms, atom_p, symmetry, R, R_tables, job, file);
  allocate_PAIR_TRAN(&pair_p, atoms, symmetry, R_tables, job, file);
  generate_density_pairs4(&pair_p, atoms, atom_p, symmetry, R, R_tables, job, file);
  print_pairs(&pair_p, atoms, R, job, file);
  array_dimensions(&dim, &dimg, &pair_p, atoms_ax, job, file); // don't change
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  fock_element_1e2(&one_ints, &pair_p, Function, R, G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  AllocateDoubleMatrix(&V,&dim1ax,&dim1ax,job);
  fourier_transform_molecule(&one_ints.Coulomb[0], &V->a[0][0], &pair_p, R, atoms_ax, shells_ax, symmetry, job, file);
  //fprintf(file.out,"V\n");
  //print_double_matrix(V, file);
  free_INT_1E(&one_ints, Function, job, file);
  free_PAIR_TRAN(&pair_p,job);
  AllocateDoubleMatrix(&xtrn1,&dim1ax,&dim1ax,job);
  AllocateDoubleMatrix(&eigenvectors,&dim1ax,&dim1ax,job);
  AllocateDoubleArray(&eigenvalues,&dim1ax,job);
  ResetDoubleMatrix(xtrn1);
  ResetDoubleMatrix(eigenvectors);
  ResetDoubleArray(eigenvalues,&dim1ax);
  DiagonaliseSymmetrical(&V, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);
  DestroyDoubleMatrix(&V,job);
  for (i = 0; i < dim1ax; i++) {
    for (j = 0; j < dim1ax; j++) {
      xtrn1->a[i][j] = eigenvectors->a[i][j] / eigenvalues[i];
     }
    }
  ResetDoubleMatrix(V_inv);
  DoubleGEMM(&Trans, &NoTrans, &alpha, &xtrn1, &eigenvectors, &beta, &V_inv);
  DestroyDoubleMatrix(&xtrn1,job);
  DestroyDoubleMatrix(&eigenvectors,job);
  DestroyDoubleArray(&eigenvalues,&dim1ax,job);
  //fprintf(file.out,"V inverse\n");
  //print_double_matrix(V_inv, file);
  time2 = MPI_Wtime() - time1;
  //if (job->taskid == 0 && job->verbosity > 1) printf("Coulomb Matrix Inverse %10.4f\n",time2);

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
  char buf2[110], xy[14] = "/scf_evec_spk";
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
*/
