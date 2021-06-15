/*
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <xc.h>
#include <unistd.h>
#include <stdlib.h>
#include "mycomplex.h"
#include "mylogical.h"
#include "MATRIX_UTIL.h"
#include "SYMMETRY.h"
#include "TOOLS.h"
#include "CRYSTAL09.h"
#include "PARALLEL.h"
#include "PRINT_UTIL.h"
#include "PLOTTING.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "INTEGRALS.h"
#include "INTEGRALS_TWO_CENTRE.h"
#include "INTEGRALS_TEST.h"
#include "INTEGRALS_THREE_CENTRE.h"
#include "INTEGRALS1.h"
#include "DENSITY_FITTING.h"
#include "ROTATION_OPERATORS.h"
#include "ROTATIONS_MOLECULE.h"
#include "FOURIER_TRANSFORM.h"
#include "SETUP_SYMMETRY.h"
#include "SETUP_RECIPROCAL_LATTICE.h"
#include "ATOM_SCF.h"
#include "DFT.h"
#include "mkl_blas.h"
#include "mkl_pblas.h"
#include "mkl_blacs.h"
#include "mkl_cblas.h"
#include "LIMITS.h"
#include "BUILD_FOCK_MATRIX.h"
*/
#define MKL_Complex16 Complex
#define MKL_INT16 short int
#include <cstring>
#include "mycomplex.h"
#include "myconstants.h"
#include "conversion_factors.h"
#include "mkl_scalapack.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "PRINT_MOLECULE.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "SETUP_SYMMETRY.h"
#include "FOURIER_TRANSFORM.h"
#include "SETUP_RECIPROCAL_LATTICE.h"
#include "PAIRS_QUADS.h"
#include "PARALLEL.h"
#include "KPOINTS.h"
#include "ROTATIONS_MOLECULE.h"
#include "INTEGRALS1.h"
#include "INTEGRALS_2C_MOLECULE.h"
#include "INTEGRALS_2C_CRYSTAL.h"
#include "DENSITY_FITTING_CRYSTAL.h"
#include "ALLOCATE_MEMORY.h"
#include "SCALAPACK.h"
#include "TDHF_CRYSTAL.h"
extern "C" void  Cblacs_pinfo( int*, int*);
extern "C" void  Cblacs_get( int context, int request, int* value);
extern "C" int   Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern "C" int   Cblacs_gridmap( int* context, int *imap, int ldimap , int np_row, int np_col);
extern "C" void  Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern "C" void  Cblacs_gridexit( int context);
extern "C" void  Cblacs_exit( int error_code);
extern "C" void  Cblacs_barrier(int context, char *scope);
extern "C" int   Cblacs_pnum(int, int, int);

using namespace std;

void bse_crystal1(FERMI* fermi, ATOM* atoms, ATOM_TRAN* atom_p, int *numfrag, int *natom, int nat[][2], SHELL* shells, GAUSSIAN* gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Routine calculates BSE and BSE-TDA excitations                                         *
  // ******************************************************************************************
  
  KPOINT_TRAN knet;
  count_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  if (job->kss == 0)      fermi->nkunique = knet.nktot;
  else if (job->kss == 1) fermi->nkunique = knet.unique;
  allocate_k_points(&knet,crystal,job,file);
  generate_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //print_knet(&knet, fermi->is, crystal, job, file);
  allocate_fermi(fermi,atoms,job,file);
  fermi->knet = &knet;
  fermi->nktot = knet.nktot;
  fermi->occupied[0] = fermi->homo[0] - fermi->bands[0] + 1;
  if (job->spin_dim == 2)
  fermi->occupied[fermi->nkunique] = fermi->homo[1] - fermi->bands[2] + 1;

int i, q, dummy;
int band_range_k[24], band_range_kq[24], nband_k[12], nband_kq[12];
int ntransitions;
int count;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int ictxt, nbsize_row, nbsize_col, cblacs_taskid, itemp, nprow, npcol, myrow, mycol, mpA, nqA, izero = 0, ione = 1;
int descA[9];
int info = 0;
int buffer_size1, buffer_size2;
int *dim_send;
Complex *integral_buffer1, *integral_buffer2;
double *bse_eigenvalues;
double time1, time2, time3, time4, time5, time6;
double *GW_eigenvalues, *scf_eigenvalues;
char yy[20] = "./bse_eigenvalues";
FILE *bse_evalues;
Complex *Ham_buffer1, *Ham_buffer2;

  time1 = MPI_Wtime();
  setup_hamiltonian_parameters_zero_q(&ntransitions, band_range_k, band_range_kq, nband_k, nband_kq, fermi, job, file);
  initialise_spk_grid_crystal(&ntransitions, fermi, dummy, &ictxt, &nbsize_row, &nbsize_col, job, file);
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, &nbsize_row, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, &nbsize_col, &mycol, &izero, &npcol);
  itemp = max(1, mpA);
  descinit_(descA, &ntransitions, &ntransitions, &nbsize_row, &nbsize_col, &izero, &izero, &ictxt, &itemp, &info);
  //if (job->taskid == 0) printf("mpA %7d nqA %7d blocksizes %7d %7d ntransitions %7d\n",mpA,nqA,nbsize_row,nbsize_col,ntransitions);

  // ******************************************************************************************
  // * Allocate memory for BSE eigenvectors and eigenvalues                                   *
  // ******************************************************************************************

  buffer_size1 = mpA * nqA;
  if      (job->bse_tda == 0) buffer_size2 = buffer_size1;
  else if (job->bse_tda == 1) buffer_size2 = 1;
  AllocateDoubleArray(&bse_eigenvalues,&ntransitions,job);
  AllocateComplexArray(&Ham_buffer1,&buffer_size1,job);
  AllocateComplexArray(&Ham_buffer2,&buffer_size2,job);
  ResetDoubleArray(bse_eigenvalues,&ntransitions);
  ResetComplexArray(Ham_buffer1,&buffer_size1);
  ResetComplexArray(Ham_buffer2,&buffer_size2);

  // ******************************************************************************************
  // * Generate three centre integrals needed for matrix elements                             *
  // * Generate RPA, BSE or BSE-TDA Hamiltonian                                               *
  // ******************************************************************************************
  
  int nbsize = nbsize_col;

  tdhf_hamiltonian_in_core_crystal_finite_q(&ictxt,&nbsize_row,&nbsize_col,Ham_buffer1,Ham_buffer2,fermi,atom_p,atoms,shells,\
  gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  tdhf_hamiltonian_in_core_crystal_zero_q(&ictxt,&nbsize_row,&nbsize_col,Ham_buffer1,Ham_buffer2,fermi,atom_p,atoms,shells,\
  gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  time2 = MPI_Wtime() - time1;
  MPI_Barrier(MPI_COMM_WORLD);

  time3 = MPI_Wtime();
  switch (job->bse_tda) {

  case 0: // Full-matrix BSE or TDHF

  //diagonalise_bse_hamiltonian_crystal(&ictxt,&nbsize,Ham_buffer1,Ham_buffer2,fermi,job,file);

  break;

  case 1: // Tamm-Dancoff Approximation: TDHF or BSE

   if (job->taskid >= 0 && job->verbosity > 1) {
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"|                                     GW BSE HAMILTONIAN MATRIX (eV)                                      |\n");
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

int j, H = 1;
   for (i = 0; i < mpA; i++) {
     for (j = 0; j < (nqA > 16 ? 16 : nqA); j++) {
double fac = 1.0 * au_to_eV;
       fprintf(file.out,"%6.3f",fac*(Ham_buffer1[i + mpA * j]).real());
       if (((j+1) / H) * H == j+1) fprintf(file.out," ");
      }
       fprintf(file.out,"\n");
       if (((i+1) / H) * H == i+1) fprintf(file.out,"\n");
     }
       fprintf(file.out,"\n");
       fprintf(file.out,"\n");
   for (i = 0; i < mpA; i++) {
     for (j = 0; j < (nqA > 16 ? 16 : nqA); j++) {
double fac = 1.0;
if (i == j) fac = 1.0*au_to_eV;
       fprintf(file.out,"%6.3f",fac*(Ham_buffer1[i + mpA * j]).imag());
       if (((j+1) / H) * H == j+1) fprintf(file.out," ");
      }
       fprintf(file.out,"\n");
       if (((i+1) / H) * H == i+1) fprintf(file.out,"\n");
     }
  } // close if (job->taskid == 0

  diagonalise_bse_tda_hamiltonian_crystal(&ictxt,&nbsize,Ham_buffer1,fermi,job,file);
  //zheevx diagonalise_bse_tda_hamiltonian_crystal2(&ictxt,&nbsize_row,&nbsize_col,Ham_buffer1,fermi,job,file);

  break;

 } // close switch (job->bse_tda)
  time4 = MPI_Wtime() - time3;

  char yy1[20] = "/bse_eigenvalues";
  read_bse_eigenvalues_crystal(bse_eigenvalues, ntransitions, yy1, job, file);
  if (job->taskid == 0) {
  //bse_evalues  = fopen(yy, "rb");
  //fread(bse_eigenvalues, sizeof(double), ntransitions, bse_evalues);
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
  free_k_points(&knet,job);
  DestroyDoubleArray(&bse_eigenvalues,&ntransitions,job);
  DestroyComplexArray(&Ham_buffer1,&buffer_size1,job);
  DestroyComplexArray(&Ham_buffer2,&buffer_size2,job);
  if (job->taskid == 0) printf("bse                                        %10.2f\n",time2+time4);

}

void tdhf_hamiltonian_in_core_crystal_finite_q(int *ictxt, int *nbsize_row, int *nbsize_col, Complex *Ham_buffer1, Complex *Ham_buffer2, FERMI* fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k, l, m, n, p, q, r, s, t, j1, j2, s1;
int nd6, q1;
int kp, lp, pm, op;
int k1, k2, l2;
int dimk, dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int dim1_kunique = dim1 * fermi->nkunique;
int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
int nk[2], nkq[2], nbands[2];
int nocc, ntransitions, ntrans[2];
int band_range_k[24], band_range_kq[24], nband_k[12], nband_kq[12], buffer_size[24];
int a1, a2, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int bra_offset, ket_offset;
int count_bz, count_bz1, dim_bra, dim_ket, Dim_send;
int fbz1, fbz2, k_bz, kq_bz;
int I1, I2;
int mpA, nqA, row, col, nprow, npcol, myrow, mycol, izero = 0;
int MPA[job->numtasks],NQA[job->numtasks];;
int cblacs_taskid, nprow_nbsize, nprow_myrow, npcol_nbsize, npcol_mycol;
int begin_j[job->numtasks], end_j[job->numtasks];
int begin_q1[job->numtasks], end_q1[job->numtasks];
int *dim_ham, *offset_j;
int offset, offset1, offset2;
int band_offset1, band_offset2;
int seek_spin_offset, seek_point_k, seek_point_kq;
int spin_offset_row, spin_offset_col, spin_offset_msg;
int dest_row, dest_col;
int array_count, array_count1;
int p1, p2, offset4;
int block, size, dim_ham_buffer;
double ferm = double (fermi->is[0] * fermi->is[1] * fermi->is[2]);
double exchange_buffer, exchange_energy;
double q0_part;
double factor;
double time1, time2, time3, time4;
double time5, time6, time7, time8;
double time9, time10, time11, time12;
double time13, time14, time15, time16;
double time17, time18, time19, time20;
double time21, time22;
double *GW_eigenvalues, *scf_eigenvalues;
char buf2[110], xy[14] = "/scf_evec";
char zz4[24] = "/evalfile";
char zz5[24] = "GW_evalues";
Complex *Hamiltonian_buffer;
ComplexMatrix *S1, *S2, *S3, *S4, *S_k, *S_kq, *V_q;
ComplexMatrix *scf_0, *scf_1, *scf_k0, *scf_kq0, *scf_k1, *scf_kq1;
MPI_File fh;

  time2  = k_zero;
  time4  = k_zero;
  time6  = k_zero;
  time8  = k_zero;
  time10 = k_zero;
  time12 = k_zero;
  time14 = k_zero;
  time16 = k_zero;
  time18 = k_zero;
  time20 = k_zero;
  time22 = k_zero;

  time1 = MPI_Wtime();

  nbands[1] = 0;
  nbands[0] = fermi->bands[1] - fermi->bands[0] + 1;
  if (job->spin_dim == 2)
  nbands[1] = fermi->bands[3] - fermi->bands[2] + 1;
  nocc = fermi->homo[0] - fermi->bands[0] + 1;
  if (job->spin_dim == 2) nocc += fermi->homo[1] - fermi->bands[2] + 1;
  dimk = (nbands[0] + nbands[1]) * fermi->nkunique;

  // ******************************************************************************************
  // * Open MPI wave function eigenvector file                                                *
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xy);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;

  AllocateDoubleArray(&scf_eigenvalues,&dimk,job);
  AllocateDoubleArray(&GW_eigenvalues,&dimk,job);
  ResetDoubleArray(scf_eigenvalues,&dimk);
  read_scf_GW_eigenvalues_crystal(scf_eigenvalues, fermi, atoms, zz4, job, file);

  // ******************************************************************************************
  // * Allocate q_G array                                                                     *
  // * Allocate dim_ham array                                                                 *
  // * Initialise processors                                                                  *
  // ******************************************************************************************

  MPI_Request request_a, request_b, request_exchange;
  MPI_Status  status;
  MPI_Status  status1;

  ntrans[0] = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1);
  ntrans[1] = 0;
  if (job->spin_dim == 2) 
  ntrans[1] = (fermi->bands[3] - fermi->homo[1]) * (fermi->homo[1] - fermi->bands[2] + 1);

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);

  AllocateIntArray(&dim_ham,&job->numtasks,job);
  setup_hamiltonian_parameters_finite_q(&ntransitions, band_range_k, band_range_kq, nband_k, nband_kq, fermi, job, file);
  setup_procs(begin_j,end_j,begin_q1,end_q1,MPA,NQA,dim_ham,&dim1,ictxt,&ntransitions,nbsize_row,nbsize_col,fermi,file,job);
  mpA = numroc_(&ntransitions, nbsize_row, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize_col, &mycol, &izero, &npcol);

  Q_LATTICE q_G;
  q_G.last_vector = G->last_vector; 
  q_G.max_vector  = G->max_vector;
  allocate_Q_LATTICE(&q_G, job, file);

  // ******************************************************************************************
  // * Generate range selected pairs                                                          *
  // * Generate real space overlap matrix                                                     *
  // ******************************************************************************************

  time19 = MPI_Wtime();
  PAIR_TRAN pair_p;
  pair_p.cutoff = 8.0;
  //fprintf(file.out,"pair_cutoff = %10.4f\n",pair_p.cutoff);
  count_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //print_pairs(&pair_p,atoms,R,job,file);

  INT_1E one_ints, one_ints_buffer, one_ints_buffer1;
  int dim, dimg;
  int Function[8];
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  int vector_size = nbands[0] * nbfn;

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  array_dimensions(&dim, &dimg, &pair_p, atoms, job, file);
  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p,end_p,receive_p,offset_p,&pair_p,atoms,job,file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  fock_element_1e2(&one_ints, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  time20 += MPI_Wtime() - time19;

  // ******************************************************************************************
  // * Send Irecv requests for data for Ham_buffer1 and Ham_buffer2                           *
  // ******************************************************************************************

  time11 = MPI_Wtime();
  issue_MPI_IRecv(&request_a,Ham_buffer1,MPA,NQA,ictxt,nbsize_row,nbsize_col,fermi,crystal,file,job);

  //if (job->bse_tda == 0)
  //issue_MPI_IRecv(&request_b,Ham_buffer2,MPA,NQA,ictxt,nbsize_row,nbsize_col,fermi,crystal,file,job);
  time12 += MPI_Wtime() - time11;

  // ******************************************************************************************
  // * Build Hamiltonian                                                                      *
  // ******************************************************************************************

  AllocateComplexMatrix(&V_q,&dim1ax,&dim1ax,job);
  exchange_buffer = k_zero;

  //time13 = MPI_Wtime();

  int little_q_group_unique, num_q1, nkunique, kq, q1_list[fermi->nkunique],begin_k[fermi->nkunique], end_k[fermi->nkunique];
  KQPOINT_TRAN kq_pair;
  KPOINT_TRAN knet_little_q_group;
  SYMMETRY symmetry_little_q_group;
  select_kq_pairs(&num_q1,q1_list,begin_k,end_k,&knet_little_q_group,fermi,symmetry,crystal,&q_G,G,job,file);
  for (i = 0; i < num_q1; i++) {
  //exchange_buffer = k_zero;
  q1 = q1_list[i];
  q = fermi->knet->ibz[q1];
  nkunique = end_k[i] - begin_k[i];
  //time17 = MPI_Wtime();
  setup_hamiltonian_parameters_finite_q(&ntransitions, band_range_k, band_range_kq, nband_k, nband_kq, fermi, job, file);
  count_little_k_group_operators(q1,&symmetry_little_q_group,symmetry,crystal,fermi->knet,fermi,job,file);
  allocate_SYMMETRY(&symmetry_little_q_group,job,file);
  generate_little_k_group(q1, &symmetry_little_q_group, fermi, fermi->knet, symmetry, crystal, job, file);
  count_k_points(&knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
  allocate_k_points(&knet_little_q_group,crystal,job,file);
  generate_k_points(&knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
  generate_q_lattice(&fermi->knet->oblique[q], &q_G, fermi, G, crystal, job, file);
  allocate_kq_pairs(&kq_pair,fermi,crystal,symmetry,job,file);
  generate_kq_pairs(q1,&kq_pair,crystal,symmetry,fermi,job,file);
  little_q_group_unique = knet_little_q_group.unique;
  //time18 += MPI_Wtime() - time17;

  count_bz = 0; for (k = 0; k < begin_k[i]; k++) count_bz += kq_pair.num[k];

  // ******************************************************************************************
  // * Calculate Overlap matrices                                                             *
  // ******************************************************************************************


  int count = 0;
  for (s1 = 0; s1 < job->spin_dim; s1++) 
  for (k = begin_k[i]; k < end_k[i]; k++)
  for (k1 = 0; k1 < kq_pair.num[k]; k1++) count++;

  ComplexMatrix *SS1[count], *SS2[count], *SS3[count], *SS4[count];
  count = 0;
  array_count = 0;
  array_count1 = 0;
  count_bz1 = count_bz;
  for (block = 0; block < 2 - job->bse_tda; block++) {
    for (s1 = 0; s1 < job->spin_dim; s1++) {
      AllocateComplexMatrix(&scf_0,&nband_k[array_count],&nbfn,job);
      AllocateComplexMatrix(&scf_k0,&nband_k[array_count],&nbfn,job);
      AllocateComplexMatrix(&scf_kq0,&nband_kq[array_count],&nbfn,job);
      AllocateComplexMatrix(&scf_k1,&nband_k[array_count + 1],&nbfn,job);
      AllocateComplexMatrix(&scf_1,&nband_kq[array_count + 1],&nbfn,job);
      AllocateComplexMatrix(&scf_kq1,&nband_kq[array_count + 1],&nbfn,job);

      count_bz = count_bz1;
      bra_offset = 0;
      ket_offset = 0;
      for (k = begin_k[i]; k < end_k[i]; k++) {

      k_bz  = kq_pair.bz1[count_bz];
      kq_bz = kq_pair.bz2[count_bz];
      fbz1 = fermi->knet->fbz[k_bz];
      fbz2 = fermi->knet->fbz[kq_bz];
      seek_spin_offset = 0;
      //if (job->spin_dim == 2 && (((s1 + 1) / 4) * 4 == s1 + 1 || ((s1 + 2) / 4) * 4 == s1 + 2)) 
      seek_spin_offset = s1 * fermi->nkunique * nbfn * nbfn;

      int vector_size_0 = nband_k[array_count] * nbfn;
      int vector_size_1 = nband_k[array_count + 1] * nbfn;
      seek_point_k  = seek_spin_offset + fbz1 * nbfn * nbfn  + band_range_k[array_count1 + 0] * nbfn;
      MPI_File_seek(fh, seek_point_k * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_0->a[0][0], 2 * vector_size_0, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_0->a[0][0],&scf_k0->a[0][0],nband_k[array_count],k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);

      seek_point_kq = seek_spin_offset + fbz2 * nbfn * nbfn  + band_range_kq[array_count1 + 0] * nbfn;
      MPI_File_seek(fh, seek_point_kq * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_0->a[0][0], 2 * vector_size_0, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_0->a[0][0],&scf_kq0->a[0][0],nband_kq[array_count],kq_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);

      seek_point_k  = seek_spin_offset + fbz1 * nbfn * nbfn  + band_range_k[array_count1 + 2] * nbfn;
      MPI_File_seek(fh, seek_point_k * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_1->a[0][0], 2 * vector_size_1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_1->a[0][0],&scf_k1->a[0][0],nband_k[array_count + 1],k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);

      seek_point_kq = seek_spin_offset + fbz2 * nbfn * nbfn  + band_range_kq[array_count1 + 2] * nbfn;
      MPI_File_seek(fh, seek_point_kq * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_1->a[0][0], 2 * vector_size_1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_1->a[0][0],&scf_kq1->a[0][0],nband_kq[array_count + 1],kq_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);

      time13 = MPI_Wtime();
      for (k1 = 0; k1 < kq_pair.num[k]; k1++) {
        AllocateComplexMatrix(&SS1[count],&nband_k[array_count],&nband_k[array_count],job);
        AllocateComplexMatrix(&SS2[count],&nband_k[array_count + 1],&nband_k[array_count + 1],job);
        AllocateComplexMatrix(&SS3[count],&nband_kq[array_count],&nband_kq[array_count],job);
        AllocateComplexMatrix(&SS4[count],&nband_kq[array_count + 1],&nband_kq[array_count + 1],job);
        k_bz  = kq_pair.bz1[count_bz];
        kq_bz = kq_pair.bz2[count_bz];
        band_offset1 = fermi->homo[s1];
        band_offset2 = fermi->bands[2 * s1] - 1;
        overlap_matrices5(SS1[count],SS2[count],scf_k0,scf_k1,&kq_pair,&s1,&band_offset1,&band_offset2,&nband_k[array_count],\
        &nband_k[array_count + 1],&count_bz,&one_ints,pair_p,atoms,atom_p,shells,R,fermi,symmetry,crystal,&fh,file,job);
        overlap_matrices6(SS3[count],SS4[count],scf_kq0,scf_kq1,&kq_pair,&s1,&band_offset1,&band_offset2,&nband_kq[array_count],\
        &nband_kq[array_count + 1],&count_bz,&one_ints,pair_p,atoms,atom_p,shells,R,fermi,symmetry,crystal,&fh,file,job);
        count++;
        count_bz++;
       } // close loop on k1
       time14 += MPI_Wtime() - time13;
      } // close loop on k

        DestroyComplexMatrix(&scf_0,job);
        DestroyComplexMatrix(&scf_1,job);
        DestroyComplexMatrix(&scf_k0,job);
        DestroyComplexMatrix(&scf_kq0,job);
        DestroyComplexMatrix(&scf_k1,job);
        DestroyComplexMatrix(&scf_kq1,job);
        array_count += 2;
        array_count1 += 4;
       }
      }

  // ******************************************************************************************
  // * Calculate Coulomb matrix V_q                                                           *
  // ******************************************************************************************

  time5 = MPI_Wtime();
  ResetComplexMatrix(V_q);
  generate_coulomb_matrix_inverse_complex(V_q,q,fermi,atom_p,atoms,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  //fprintf(file.out,"V_q\n");
  //print_complex_matrix(V_q,file);
  time6 += MPI_Wtime() - time5;

  // ******************************************************************************************
  // * Allocate integral_buffers                                                              *
  // * Calculate <ij|a> matrix elements                                                       *
  // ******************************************************************************************

  time3 = MPI_Wtime();
  Complex *integral_buffers[job->band_dim], *integral_tmp, *bra_integrals, *ket_integrals;
  for (s = 0; s < job->band_dim; s++) {
    buffer_size[s] = 0;
    for (j1 = 0; j1 < dim1; j1++) buffer_size[s] += nband_k[s] * nband_kq[s] * nkunique * atoms_ax->bfnnumb_sh[j1];
    AllocateComplexArray(&integral_buffers[s],&buffer_size[s],job);
    ResetComplexArray(integral_buffers[s],&buffer_size[s]);
    //printf("Buffer_size %3d %3d %3d %3d %3d\n",job->taskid,s,buffer_size[s],nband_k[s],nband_kq[s]);
   } // close s loop

    density_fitting_crystal_contract_integrals(&q1,&begin_k[i],&end_k[i],band_range_k,band_range_kq,&pair_p,&knet_little_q_group,fh,\
    integral_buffers,fermi,atom_p,atoms,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,&symmetry_little_q_group,\
    &q_G,R,R_tables,G,job,file);
    time4 += MPI_Wtime() - time3;
    time7 = MPI_Wtime();
    for (s = 1; s < job->band_dim; s++) {
      AllocateComplexArray(&integral_tmp,&buffer_size[s],job);
      for (t = 0; t < buffer_size[s]; t++) integral_tmp[t] = integral_buffers[s][t];
        ResetComplexArray(integral_buffers[s],&buffer_size[s]);
        for (m = 0; m < nkunique * nband_k[s] * nband_kq[s]; m++) {
          for (a1 = 0; a1 < dim1ax; a1++) {
            for (a2 = 0; a2 < dim1ax; a2++) {
              integral_buffers[s][m * dim1ax + a1] += V_q->a[a1][a2] * conj(integral_tmp[m * dim1ax + a2]);
              //fprintf(file.out,"m %3d a1,2 %3d %3d off2 %3d V_q[a2][a1] %3d %3d  off1 %3d %10.4lf %10.4lf %10.4lf\n",\
              m,a1,a2, m * dim1ax + a1,a1,a2,m*dim1ax + a2,\
              integral_buffers[2 * s + 1][m * dim1ax + a1].real(),V_q->a[a1][a2].real(),integral_tmp[s][m * dim1ax + a2].real());
             }
            }
           }
          DestroyComplexArray(&integral_tmp,&buffer_size[s],job);
          s++; // multiply only the kets by V_q
         }
    time8 += MPI_Wtime() - time7;

  // ******************************************************************************************
  // * Loop over unique k points for this q vector                                            *
  // ******************************************************************************************

time15 = MPI_Wtime();
array_count = 0;
array_count1 = 0;
count = 0;
  for (block = 0; block < 2 - job->bse_tda; block++) {
    for (s1 = 0; s1 < job->spin_dim; s1++) {
      int vector_size_k  = nband_k[array_count]  * nbfn;
      int vector_size_kq = nband_kq[array_count + 1]  * nbfn;
      dim_bra = nband_k[array_count] * nband_kq[array_count] * dim1ax;
      dim_ket = nband_k[array_count + 1] * nband_kq[array_count + 1] * dim1ax;
      dim_ham_buffer = nband_k[array_count] * nband_kq[array_count] * nband_k[array_count + 1] * nband_kq[array_count + 1];
      AllocateComplexArray(&Hamiltonian_buffer,&dim_ham_buffer,job);
      AllocateComplexArray(&bra_integrals,&dim_bra,job);
      AllocateComplexArray(&ket_integrals,&dim_ket,job);
      count_bz = count_bz1;
      bra_offset = 0;
      ket_offset = 0;
      for (k = begin_k[i]; k < end_k[i]; k++) {

      k_bz  = kq_pair.bz1[count_bz];
      kq_bz = kq_pair.bz2[count_bz];
      fbz1 = fermi->knet->fbz[k_bz];
      fbz2 = fermi->knet->fbz[kq_bz];
      //fprintf(file.out,"s1 k num %3d %3d %3d end %3d arr %3d range %3d %3d %3d %3d     %3d %3d %3d %3d %3d %3d \n",\
      s1,k,kq_pair.num[k],end_k[i],array_count,band_range_k[array_count1 + 0],band_range_kq[array_count1 + 0],\
      band_range_k[array_count1 + 2],band_range_kq[array_count1 + 2],fermi->homo[s1],fermi->bands[2 * s1],band_range_kq[4],\
      band_range_kq[5],band_range_kq[6],band_range_kq[7]);

        for (k1 = 0; k1 < kq_pair.num[k]; k1++) {
          time17 = MPI_Wtime();
          ResetComplexArray(bra_integrals,&dim_bra);
          density_fitting_crystal_rotate_integrals(&bra_offset,1,&kq_pair,SS1[count],SS3[count],&integral_buffers[array_count],bra_integrals,\
          &nband_k[array_count],&nband_kq[array_count],&count_bz,atom_p,atoms_ax,shells_ax,symmetry,file,job);
          ResetComplexArray(ket_integrals,&dim_ket);
          density_fitting_crystal_rotate_integrals(&ket_offset,0,&kq_pair,SS2[count],SS4[count],&integral_buffers[array_count + 1],\
          ket_integrals,&nband_k[array_count + 1],&nband_kq[array_count + 1],&count_bz,atom_p,atoms_ax,shells_ax,symmetry,file,job);
          ResetComplexArray(Hamiltonian_buffer,&dim_ham_buffer);
          time18 += MPI_Wtime() - time17;

          if (block == 0 && q == 0) {  //A matrix
          p = 0;
          offset1 = fbz1 * (nbands[0] + nbands[1]) + s1 * nbands[0];
          for (m = 0; m < nband_k[array_count + 1]; m++) {
            for (n = nband_k[array_count + 1]; n < nband_k[array_count + 1] + nband_k[array_count]; n++) {
              //Ham_buffer1[offset + p * mpA + p] += scf_eigenvalues[offset1 + n] - scf_eigenvalues[offset1 + m];
              Hamiltonian_buffer[p + p * nband_k[array_count + 1] * nband_k[array_count]] += \
              scf_eigenvalues[offset1 + n] - scf_eigenvalues[offset1 + m] - job->scissor;
              Hamiltonian_buffer[p + p * nband_k[array_count + 1] * nband_k[array_count]] -= \
              pow(0.75 / pi / crystal->primitive_cell_volume / ferm, 0.3333333) * four * job->scalefac;
              p++;
             }
            }
           } // close if (block

          p1 = 0;
          for (r = 0; r < nband_k[array_count + 1]; r++) {      // v
            for (m = 0; m < nband_k[array_count]; m++) {        // c
              p2 = 0;
              for (s = 0; s < nband_kq[array_count + 1]; s++) { // v'
                for (n = 0; n < nband_kq[array_count]; n++) {   // c'
                  for (a1 = 0; a1 < dim1ax; a1++) {
                    Hamiltonian_buffer[p1 + p2 * nband_k[array_count + 1] * nband_k[array_count]] -= \
                    bra_integrals[m * nband_kq[array_count]     * dim1ax + n * dim1ax + a1] * \
                    ket_integrals[r * nband_kq[array_count + 1] * dim1ax + s * dim1ax + a1] / \
                    (double) (fermi->knet->nktot) * job->scalefac;
                    //fprintf(file.out,"s1 %3d p1,p2 %3d %3d  %3d %14.8f %14.8f    %14.8f\n",s1,p1,p2,a1,\
                    (bra_integrals[m * nband_kq[array_count]     * dim1ax + n * dim1ax + a1]).real(),\
                    (ket_integrals[r * nband_kq[array_count + 1] * dim1ax + s * dim1ax + a1]).real(),\
                    (Hamiltonian_buffer[p1 + p2 * nband_k[array_count + 1] * nband_k[array_count]]).real());
                   } 
                  p2++;
                  } 
                 } 
                p1++;
                } 
               } 
          
          // B matrix
          if (block == 1 && q == 0) { 
           }

           spin_offset_row = s1 * ntrans[0];
           spin_offset_col = s1 * ntrans[0];
           spin_offset_msg = s1 * fermi->nktot * fermi->nktot;
           size = dim_ham_buffer;
           k_bz  = kq_pair.bz1[count_bz];
           kq_bz = kq_pair.bz2[count_bz];
           I1 = spin_offset_row + k_bz  * (ntrans[0] + ntrans[1]);
           I2 = spin_offset_col + kq_bz * (ntrans[0] + ntrans[1]);
           dest_row = (I1 / *nbsize_row) % nprow;
           dest_col = (I2 / *nbsize_col) % npcol;
           cblacs_taskid = Cblacs_pnum(*ictxt,dest_row,dest_col);
           //if (k1==0) printf("%3d sending message %8d to core %3d %3d %14.8f\n", \
           job->taskid, k_bz * fermi->nktot + kq_bz, cblacs_taskid,size,(Hamiltonian_buffer[0]).real());
           //fflush(stdout);
           //fprintf(file.out,"%3d %3d %3d sending message %8d to core %3d %3d %14.8f\n", \
           job->taskid,k,k1, k_bz * fermi->nktot + kq_bz, cblacs_taskid,size,(Hamiltonian_buffer[0]).real());
           //fflush(file.out);
           time9 = MPI_Wtime();
           MPI_Send(Hamiltonian_buffer, 2 * size, MPI_DOUBLE, cblacs_taskid, spin_offset_msg + k_bz * fermi->nktot + kq_bz, MPI_COMM_WORLD);
           time10 += MPI_Wtime() - time9;
           //if(k1==kq_pair.num[k]-1) printf("%3d sent message    %8d to core %3d %9.2e\n", job->taskid, k_bz * fermi->nktot + kq_bz, cblacs_taskid,tim);
           //fprintf(file.out,"%3d %3d %3d sent message    %8d to core %3d %9.2e\n", job->taskid, k,k1,k_bz * fermi->nktot + kq_bz, cblacs_taskid,tim);
           //fflush(stdout);
           //fflush(file.out);
           count_bz++;
           count++;
           } // close loop on k1
           bra_offset += dim_bra;
           ket_offset += dim_ket;
          } // close loop on k
         array_count += 2;
         array_count1 += 4;
         DestroyComplexArray(&bra_integrals,&dim_bra,job);
         DestroyComplexArray(&ket_integrals,&dim_ket,job);
         DestroyComplexArray(&Hamiltonian_buffer,&dim_ham_buffer,job);
        } // close loop on s1
       } // close loop on block

  time16 += MPI_Wtime() - time15;

  if      (job->spin_dim == 1) { array_count = 2; array_count1 = 4; }
  else if (job->spin_dim == 2) { array_count = 4; array_count1 = 8; }
  
  if (job->bse_exc == 1) { 
  for (s1 = 0; s1 < job->spin_dim; s1++) {
    dim_bra = nband_k[array_count] * nband_kq[array_count] * dim1ax;
    dim_ket = nband_k[array_count + 1] * nband_kq[array_count + 1] * dim1ax;
    AllocateComplexArray(&bra_integrals,&dim_bra,job);
    AllocateComplexArray(&ket_integrals,&dim_ket,job);
    AllocateComplexMatrix(&S_k,&nbfn,&nbfn,job);
    AllocateComplexMatrix(&S_kq,&nbfn,&nbfn,job);
    AllocateComplexMatrix(&S1,&nband_k[array_count],&nband_k[array_count],job);
    AllocateComplexMatrix(&S2,&nband_k[array_count + 1],&nband_k[array_count + 1],job);
    AllocateComplexMatrix(&S3,&nband_kq[array_count],&nband_kq[array_count],job);
    AllocateComplexMatrix(&S4,&nband_kq[array_count + 1],&nband_kq[array_count + 1],job);
    ComplexMatrix *scf_0, *scf_1, *scf_k0, *scf_kq0, *scf_k1, *scf_kq1;
    int vector_size_k  = nband_k[array_count]  * nbfn;
    int vector_size_kq = nband_kq[array_count + 1]  * nbfn;
    AllocateComplexMatrix(&scf_0,&nband_k[array_count],&nbfn,job);
    AllocateComplexMatrix(&scf_k0,&nband_k[array_count],&nbfn,job);
    AllocateComplexMatrix(&scf_kq0,&nband_kq[array_count],&nbfn,job);
    AllocateComplexMatrix(&scf_k1,&nband_k[array_count + 1],&nbfn,job);
    AllocateComplexMatrix(&scf_1,&nband_kq[array_count + 1],&nbfn,job);
    AllocateComplexMatrix(&scf_kq1,&nband_kq[array_count + 1],&nbfn,job);

    bra_offset = 0;
    ket_offset = 0;
    count_bz = count_bz1;
    for (k = begin_k[i]; k < end_k[i]; k++) {
 
      k_bz  = kq_pair.bz1[count_bz];
      kq_bz = kq_pair.bz2[count_bz];
      fbz1 = fermi->knet->fbz[k_bz];
      fbz2 = fermi->knet->fbz[kq_bz];
      seek_spin_offset = s1 * fermi->nkunique * nbfn * nbfn;
      int vector_size_0 = nband_k[array_count] * nbfn;
      int vector_size_1 = nband_k[array_count + 1] * nbfn;

      seek_point_k  = seek_spin_offset + fbz1 * nbfn * nbfn  + band_range_k[array_count1] * nbfn;
      MPI_File_seek(fh, seek_point_k * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_0->a[0][0], 2 * vector_size_0, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_0->a[0][0],&scf_k0->a[0][0],nband_k[array_count],k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);

      seek_point_k  = seek_spin_offset + fbz1 * nbfn * nbfn  + band_range_k[array_count1] * nbfn;
      MPI_File_seek(fh, seek_point_k * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_1->a[0][0], 2 * vector_size_1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_1->a[0][0],&scf_k1->a[0][0],nband_k[array_count + 1],k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);

      seek_point_kq = seek_spin_offset + fbz2 * nbfn * nbfn  + band_range_kq[array_count1] * nbfn;
      MPI_File_seek(fh, seek_point_kq * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_0->a[0][0], 2 * vector_size_0, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_0->a[0][0],&scf_kq0->a[0][0],nband_kq[array_count],kq_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);

      seek_point_kq = seek_spin_offset + fbz2 * nbfn * nbfn  + band_range_kq[array_count1] * nbfn;
      MPI_File_seek(fh, seek_point_kq * sizeof(Complex), MPI_SEEK_SET);
      MPI_File_read(fh, &scf_1->a[0][0], 2 * vector_size_1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_1->a[0][0],&scf_kq1->a[0][0],nband_kq[array_count + 1],kq_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);

      for (k1 = 0; k1 < kq_pair.num[k]; k1++) {
        ResetComplexArray(bra_integrals,&dim_bra);
        ResetComplexArray(ket_integrals,&dim_ket);
        k_bz  = kq_pair.bz1[count_bz];
        kq_bz = kq_pair.bz2[count_bz];
        band_offset1 = fermi->bands[2 * s1] - 1;
        band_offset2 = fermi->bands[2 * s1] - 1;
        overlap_matrices5(S1,S2,scf_k0,scf_k1,&kq_pair,&s1,&band_offset1,&band_offset2,&nband_k[array_count],&nband_k[array_count + 1],\
        &count_bz,&one_ints,pair_p,atoms,atom_p,shells,R,fermi,symmetry,crystal,&fh,file,job);
        overlap_matrices6(S3,S4,scf_kq0,scf_kq1,&kq_pair,&s1,&band_offset1,&band_offset2,&nband_kq[array_count],\
        &nband_kq[array_count + 1],&count_bz,&one_ints,pair_p,atoms,atom_p,shells,R,fermi,symmetry,crystal,&fh,file,job);
        density_fitting_crystal_rotate_integrals(&bra_offset,1,&kq_pair,S1,S3,&integral_buffers[array_count],bra_integrals,\
        &nband_k[array_count],&nband_kq[array_count],&count_bz,atom_p,atoms_ax,shells_ax,symmetry,file,job);
        density_fitting_crystal_rotate_integrals(&ket_offset,0,&kq_pair,S2,S4,&integral_buffers[array_count + 1],ket_integrals,\
        &nband_k[array_count + 1],&nband_kq[array_count + 1],&count_bz,atom_p,atoms_ax,shells_ax,symmetry,file,job);
        for (m = 0; m < nband_k[array_count]; m++) {
          for (r = 0; r < nband_kq[array_count]; r++) {
            for (a1 = 0; a1 < dim1ax; a1++) { ;
              exchange_buffer -= k_one / (double) (fermi->knet->nktot * fermi->knet->nktot) * \
             (bra_integrals[m * nband_kq[array_count] * dim1ax + r * dim1ax + a1] * \
              ket_integrals[m * nband_kq[array_count + 1] * dim1ax + r * dim1ax + a1]).real();
              //fprintf(file.out,"s1 %3d %3d %3d %3d   %14.8f %14.8f\n",s1,m,r,a1,\
             (bra_integrals[m * nband_kq[array_count] * dim1ax + r * dim1ax + a1]).real(),\
             (ket_integrals[m * nband_kq[array_count + 1] * dim1ax + r * dim1ax + a1]).real());
              //fprintf(file.out,"arr %3d q,k,k1 %3d %3d %3d integrals %3d %3d %3d %14.8f %14.8f  %14.8f %14.8f\n",array_count,q,k,k1,m,r,a1,\
             (bra_integrals[m * nband_kq[array_count] * dim1ax + r * dim1ax + a1]).real(), \
             (bra_integrals[m * nband_kq[array_count] * dim1ax + r * dim1ax + a1]).imag(), \
             (ket_integrals[m * nband_kq[array_count + 1] * dim1ax + r * dim1ax + a1]).real(),\
             (ket_integrals[m * nband_kq[array_count + 1] * dim1ax + r * dim1ax + a1]).imag());
            }
           }
          } 
         count_bz++;
        }
       bra_offset += dim_bra;
       ket_offset += dim_ket;
      }
      array_count  +=2;
      array_count1 +=4;
      DestroyComplexMatrix(&S1,job);
      DestroyComplexMatrix(&S2,job);
      DestroyComplexMatrix(&S3,job);
      DestroyComplexMatrix(&S4,job);
      DestroyComplexMatrix(&S_k,job);
      DestroyComplexMatrix(&S_kq,job);
      DestroyComplexMatrix(&scf_0,job);
      DestroyComplexMatrix(&scf_1,job);
      DestroyComplexMatrix(&scf_k0,job);
      DestroyComplexMatrix(&scf_kq0,job);
      DestroyComplexMatrix(&scf_k1,job);
      DestroyComplexMatrix(&scf_kq1,job);
      DestroyComplexArray(&bra_integrals,&dim_bra,job);
      DestroyComplexArray(&ket_integrals,&dim_ket,job);
     } // loop on s1
    } //close if (job->bse_exc == 1) { 

  count = 0;
  for (s1 = 0; s1 < job->spin_dim; s1++) {
  for (k = begin_k[i]; k < end_k[i]; k++) {
  for (k1 = 0; k1 < kq_pair.num[k]; k1++) { 
  DestroyComplexMatrix(&SS1[count],job);
  DestroyComplexMatrix(&SS2[count],job);
  DestroyComplexMatrix(&SS3[count],job);
  DestroyComplexMatrix(&SS4[count],job);
  count++;
 }}}

  for (s = 0; s < job->band_dim; s++) 
  DestroyComplexArray(&integral_buffers[s],&buffer_size[s],job);
  free_SYMMETRY(&symmetry_little_q_group,job);
  free_k_points(&knet_little_q_group,job);
  free_kq_pairs(&kq_pair,job);
 } // close i loop
  //time14 += MPI_Wtime() - time13;

  time21 = MPI_Wtime();
  MPI_Reduce(&exchange_buffer,&exchange_energy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  time22 += MPI_Wtime() - time21;
  //printf("%3d %10.4f\n",job->taskid,exchange_buffer);
  q0_part = - nocc / (double) job->spin_dim * pow(0.75 / pi / crystal->primitive_cell_volume / ferm, 0.3333333) * four;
  if (job->taskid == 0 && job->bse_exc == 1) printf("exchange energy                      %16.8f %16.8f %16.8lf\n",\
  exchange_energy / (double) job->spin_dim, q0_part, q0_part + exchange_energy / (double) job->spin_dim); 

  DestroyComplexMatrix(&V_q,job);
  free_INT_1E(&one_ints, Function, job, file);
  free_Q_LATTICE(&q_G,job);
  free_PAIR_TRAN(&pair_p,job);
  DestroyIntArray(&dim_ham,&job->numtasks,job);
  DestroyDoubleArray(&GW_eigenvalues,&dimk,job);
  DestroyDoubleArray(&scf_eigenvalues,&dimk,job);
  MPI_File_close(&fh);

  time2 += MPI_Wtime() - time1;

  if (job->taskid == 1) printf("Total Time                                 %10.2f\n",time2);
  if (job->taskid == 1) printf("Recv                                       %10.2f\n",time12);
  if (job->taskid == 1) printf("Overlap                                    %10.2f\n",time20);
  if (job->taskid == 1) printf("Overlap Matrices                           %10.2f\n",time14);
  if (job->taskid == 1) printf("Rotate integrals                           %10.2f\n",time18);
  //if (job->taskid == 1) printf("Little k group                             %10.2f\n",time18);
  if (job->taskid == 1) printf("Invert Coulomb matrix                      %10.2f\n",time6);
  if (job->taskid == 1) printf("Integrals                                  %10.2f\n",time4);
  if (job->taskid == 1) printf("Contract integrals                         %10.2f\n",time8);
  if (job->taskid == 1) printf("k loop                                     %10.2f\n",time16);
  if (job->taskid == 1) printf("Send                                       %10.2f\n",time10);
  if (job->taskid == 1) printf("Windup                                     %10.2f\n",time22);

}

void tdhf_hamiltonian_in_core_crystal_zero_q(int *ictxt, int *nbsize_row, int *nbsize_col, Complex *Ham_buffer1, Complex *Ham_buffer2, FERMI* fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k, l, m, n, p, q, r, s, t, q1, j1, j2;
int p1, p2, q2, s1, s2;
int kp, lp, pm, op;
int k1, l1, k2, l2;
int dimk, dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int ntransitions;
int ntransition[2], nbands[2];
int band_range_k[24], band_range_kq[24], nband_k[12], nband_kq[12], buffer_size[24];
int a1, a2, nd6, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int count_bz, dim_bra, dim_ket, Dim_send;
int k_bz, kq_bz;
int bra_offset, ket_offset;
int fbz1, fbz2;
int local_row, local_col, op1, op2;
int I1, I2;
int mpA, nqA, row, col, nprow, npcol, myrow, mycol, izero = 0;
int mpa, nqa, MPA[job->numtasks],NQA[job->numtasks];;
int cblacs_taskid, nprow_nbsize, nprow_myrow, npcol_nbsize, npcol_mycol;
int begin_j[job->numtasks], end_j[job->numtasks];
int begin_q1[job->numtasks], end_q1[job->numtasks];
int nkunique, kq, q1_list[fermi->nkunique],begin_k[fermi->nkunique], end_k[fermi->nkunique];
int *dim_ham;
int offset, spin_offset_row, spin_offset_col;
int num_proc;
int dest_row, dest_col;
int array_count;
double time1, time2, time3, time4;
double time5, time6, time7, time8;
double time9, time10, time11, time12;
double time13, time14, time15, time16;
double time17, time18, time19, time20;
double time21, time22;
double *GW_eigenvalues, *scf_eigenvalues;
char buf2[110], xy[14] = "/scf_evec";
char zz4[24] = "/evalfile";
char zz5[24] = "GW_evalues";
Complex *Hamiltonian_buffer;
ComplexMatrix *V_q;
MPI_File fh;

  time2  = k_zero;
  time4  = k_zero;
  time6  = k_zero;
  time8  = k_zero;
  time10 = k_zero;
  time12 = k_zero;
  time14 = k_zero;
  time16 = k_zero;
  time18 = k_zero;
  time20 = k_zero;
  time22 = k_zero;

  time1 = MPI_Wtime();

  nbands[1] = 0;
  nbands[0] = fermi->bands[1] - fermi->bands[0] + 1;
  if (job->spin_dim == 2)
  nbands[1] = fermi->bands[3] - fermi->bands[2] + 1;
  dimk = (nbands[0] + nbands[1]) * fermi->nkunique;
  ntransition[0] = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1);
  ntransition[1] = 0;
  if (job->spin_dim == 2) 
  ntransition[1] = (fermi->bands[3] - fermi->homo[1]) * (fermi->homo[1] - fermi->bands[2] + 1);

  // ******************************************************************************************
  // * Open MPI wave function eigenvector file                                                *
  // * Allocate memory for GW and SCF eigenvalues                                             *
  // * Read eigenvalues from disk                                                             *
  // ******************************************************************************************

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xy);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;

  AllocateDoubleArray(&scf_eigenvalues,&dimk,job);
  AllocateDoubleArray(&GW_eigenvalues,&dimk,job);
  ResetDoubleArray(scf_eigenvalues,&dimk);
  read_scf_GW_eigenvalues_crystal(scf_eigenvalues, fermi, atoms, zz4, job, file);

  // ******************************************************************************************
  // * Allocate q_G array                                                                     *
  // * Allocate dim_ham array                                                                 *
  // * Initialise processors                                                                  *
  // ******************************************************************************************

  AllocateIntArray(&dim_ham,&job->numtasks,job);
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  //job->bse_cou = 1;
  setup_hamiltonian_parameters_zero_q(&ntransitions, band_range_k, band_range_kq, nband_k, nband_kq, fermi, job, file);
  setup_procs(begin_j,end_j,begin_q1,end_q1,MPA,NQA,dim_ham,&dim1,ictxt,&ntransitions,nbsize_row,nbsize_col,fermi,file,job);
  mpA = numroc_(&ntransitions, nbsize_row, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize_col, &mycol, &izero, &npcol);
  //printf("nt %3d %3d %3d %3d %3d %3d %3d %3d\n",job->taskid,ntransitions,*nbsize_row,*nbsize_col,nprow,npcol,mpA,nqA);

  Q_LATTICE q_G;
  q_G.last_vector = G->last_vector; 
  q_G.max_vector  = G->max_vector;
  allocate_Q_LATTICE(&q_G, job, file);

  // ******************************************************************************************
  // * Generate range selected pairs                                                          *
  // * Generate real space overlap matrix                                                     *
  // ******************************************************************************************

  time19 = MPI_Wtime();
  PAIR_TRAN pair_p;
  pair_p.cutoff = 8.0;
  //pair_p.cutoff = 24.0;
  //fprintf(file.out,"pair_cutoff = %10.4f\n",pair_p.cutoff);
  count_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);

  // ******************************************************************************************
  // * Build Hamiltonian                                                                      *
  // ******************************************************************************************

  time13 = MPI_Wtime();

  KQPOINT_TRAN kq_pair;
  KPOINT_TRAN knet_little_q_group;
  SYMMETRY symmetry_little_q_group;

  i = 0; q1_list[0] = 0;
  q1 = q1_list[i];
  q = fermi->knet->ibz[q1];
  time17 = MPI_Wtime();
  setup_hamiltonian_parameters_zero_q(&ntransitions, band_range_k, band_range_kq, nband_k, nband_kq, fermi, job, file);
  count_little_k_group_operators(q1,&symmetry_little_q_group,symmetry,crystal,fermi->knet,fermi,job,file);
  allocate_SYMMETRY(&symmetry_little_q_group,job,file);
  generate_little_k_group(q1, &symmetry_little_q_group, fermi, fermi->knet, symmetry, crystal, job, file);
  count_k_points(&knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
  allocate_k_points(&knet_little_q_group,crystal,job,file);
  generate_k_points(&knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
  generate_q_lattice(&fermi->knet->oblique[q], &q_G, fermi, G, crystal, job, file);
  allocate_kq_pairs(&kq_pair,fermi,crystal,symmetry,job,file);
  generate_kq_pairs(q1,&kq_pair,crystal,symmetry,fermi,job,file);
  time18 += MPI_Wtime() - time17;

  begin_k[0] = 0; end_k[0] = knet_little_q_group.unique;
  nkunique = end_k[i] - begin_k[i];
  count_bz = 0;

  // ******************************************************************************************
  // * Calculate Coulomb matrix V_q                                                           *
  // ******************************************************************************************

  time5 = MPI_Wtime();
  AllocateComplexMatrix(&V_q,&dim1ax,&dim1ax,job);
  ResetComplexMatrix(V_q);
  generate_coulomb_matrix_inverse_complex(V_q,q,fermi,atom_p,atoms,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
  time6 += MPI_Wtime() - time5;

  // ******************************************************************************************
  // * Allocate integral_buffers                                                              *
  // * Calculate <ij|a> matrix elements                                                       *
  // ******************************************************************************************

  time3 = MPI_Wtime();
  Complex *integral_buffers[job->band_dim], *integral_tmp, *bra_integrals, *ket_integrals;
  for (s = 0; s < job->band_dim; s++) {
    buffer_size[s] = 0;
    for (j1 = 0; j1 < dim1; j1++) buffer_size[s] += nband_k[s] * nband_kq[s] * nkunique * atoms_ax->bfnnumb_sh[j1];
    AllocateComplexArray(&integral_buffers[s],&buffer_size[s],job);
    ResetComplexArray(integral_buffers[s],&buffer_size[s]);
    //fprintf(file.out,"Buffer_size %3d %3d %3d %3d %3d %3d\n",job->taskid,s,buffer_size[s],atoms_ax->bfnnumb_sh[0],nband_k[s],nband_kq[s]);
   } // close s loop
    density_fitting_crystal_contract_integrals(&q1,&begin_k[i],&end_k[i],band_range_k,band_range_kq,&pair_p,&knet_little_q_group,fh,\
    integral_buffers,fermi,atom_p,atoms,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,&symmetry_little_q_group,\
    &q_G,R,R_tables,G,job,file);
    time4 += MPI_Wtime() - time3;
    time7 = MPI_Wtime();

    //for (int ii = 0; ii < buffer_size[0]; ii++) fprintf(file.out,"%3d %14.8f %14.8f  %14.8f %14.8f\n",\
    ii,(integral_buffers[0][ii]).real(),(integral_buffers[0][ii]).imag(),(integral_buffers[1][ii]).real(),(integral_buffers[1][ii]).imag());
    //ComplexMatrix *solution;
    //int dim1a = nband_k[4] * nband_kq[4];
    //AllocateComplexMatrix(&solution,&dim1a,&dim1ax,job);
    //wavefunction_product_density_fit_crystal_test1(solution,&q_G,fermi,integral_buffers,&nband_k[4],&nband_kq[4],&pair_p,atom_p,atoms,\
    shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
    //DestroyComplexMatrix(&solution,job);

    for (s = 1; s < job->band_dim; s++) {
      AllocateComplexArray(&integral_tmp,&buffer_size[s],job);
      for (t = 0; t < buffer_size[s]; t++) integral_tmp[t] = integral_buffers[s][t];
        ResetComplexArray(integral_buffers[s],&buffer_size[s]);
        for (m = 0; m < nkunique * nband_k[s] * nband_kq[s]; m++) {
          for (a1 = 0; a1 < dim1ax; a1++) {
            for (a2 = 0; a2 < dim1ax; a2++) {
              integral_buffers[s][m * dim1ax + a1] += V_q->a[a1][a2] * conj(integral_tmp[m * dim1ax + a2]);
              //fprintf(file.out,"s %3d m %3d a1,2 %3d %3d off2 %3d V_q[a2][a1] %3d %3d  off1 %3d %10.4lf %10.4lf %10.4lf\n",\
              s,m,a1,a2, m * dim1ax + a1,a1,a2,m*dim1ax + a2,\
              integral_buffers[s][m * dim1ax + a1].real(),V_q->a[a1][a2].real(),integral_tmp[m * dim1ax + a2].real());
             }
            }
           }
          DestroyComplexArray(&integral_tmp,&buffer_size[s],job);
          s++; // multiply only the kets by V_q
         }
    //for (int ii = 0; ii < buffer_size[0]; ii++) fprintf(file.out,"%3d %14.8f %14.8f  %14.8f %14.8f\n",\
    ii,(integral_buffers[0][ii]).real(),(integral_buffers[0][ii]).imag(),(integral_buffers[1][ii]).real(),(integral_buffers[1][ii]).imag());
    time8 += MPI_Wtime() - time7;
    //for (s = 1; s < job->band_dim; s++) {
      //ResetComplexArray(integral_buffers[s],&buffer_size[s]);
      //for (m = 0; m < nkunique * nband_k[s] * nband_kq[s]; m++) {
        //for (a1 = 0; a1 < dim1ax; a1++) {
          //for (a2 = 0; a2 < dim1ax; a2++) {
            //integral_buffers[s][m * dim1ax + a1] = solution->a[i][j];
            //fprintf(file.out,"m %3d a1,2 %3d %3d off2 %3d V_q[a2][a1] %3d %3d  off1 %3d %10.4lf %10.4lf %10.4lf\n",\
            m,a1,a2, m * dim1ax + a1,a1,a2,m*dim1ax + a2,\
            integral_buffers[2 * s + 1][m * dim1ax + a1].real(),V_q->a[a1][a2].real(),integral_tmp[s][m * dim1ax + a2].real());
           //}
          //}
         //}
        //s++; // multiply only the kets by V_q
       //}

  // ******************************************************************************************
  // * Loop over unique k points for this q vector                                            *
  // ******************************************************************************************

  time15 = MPI_Wtime();

  for (s1 = 0; s1 < job->spin_dim; s1++) {
    for (s2 = 0; s2 < job->spin_dim; s2++) {
//if (s1 == 0 && s2 == 0 || s1 == 1 && s2 == 1) {
      spin_offset_row = s1 * ntransition[0];
      spin_offset_col = s2 * ntransition[0];
      //CHECK
      for (q2 = 0; q2 < fermi->nkunique; q2++) {
        generate_kq_pairs(q2,&kq_pair,crystal,symmetry,fermi,job,file);
        count_bz = 0;
        for (k = 0; k < kq_pair.unique; k++) {
          for (k1 = 0; k1 < kq_pair.num[k]; k1++) {
            k_bz  = kq_pair.bz1[count_bz + k1];
            kq_bz = kq_pair.bz2[count_bz + k1];
            fbz1 = fermi->knet->fbz[k_bz];
            fbz2 = fermi->knet->fbz[kq_bz];
            op1 = symmetry->inverse[fermi->knet->opr[k_bz]];
            op2 = symmetry->inverse[fermi->knet->opr[kq_bz]];
            I1 = spin_offset_row + k_bz  * (ntransition[0] + ntransition[1]);
            I2 = spin_offset_col + kq_bz * (ntransition[0] + ntransition[1]);
            dest_row = (I1 / *nbsize_row) % nprow;
            dest_col = (I2 / *nbsize_col) % npcol;
            cblacs_taskid = Cblacs_pnum(*ictxt,dest_row,dest_col);
            local_row = *nbsize_row * (I1 / (*nbsize_row * nprow)) + I1 % *nbsize_row;
            local_col = *nbsize_col * (I2 / (*nbsize_col * npcol)) + I2 % *nbsize_col;
            offset = local_row  + mpA * local_col;
 
            if (job->taskid == cblacs_taskid) {
            int j1_inv[dim1];
            int j2_inv[dim1];
            int offset_j1_inv[dim1];
            int offset_j2_inv[dim1];
            for (j1 = 0; j1 < dim1;j1++) j1_inv[j1] = -1;
            for (j2 = 0; j2 < dim1;j2++) j2_inv[j2] = -1;
            for (j1 = 0; j1 < dim1;j1++) offset_j1_inv[j1] = 0;
            for (j2 = 0; j2 < dim1;j2++) offset_j2_inv[j2] = 0;
            for (j1 = 0; j1 < dim1; j1++ ) {
              for (j2 = 0; j2 < dim1; j2++ ) {
                if (atom_p->K[j2 * symmetry->number_of_operators + op1] == j1) j1_inv[j1] = j2;
                if (atom_p->K[j2 * symmetry->number_of_operators + op2] == j1) j2_inv[j1] = j2;
               }
              }
            for (j1 = 0; j1 < dim1; j1++ ) {
              for (j2 = 0; j2 < j1_inv[j1]; j2++ ) {
                offset_j1_inv[j1] +=  atoms_ax->bfnnumb_sh[j2];
               }
              }
            for (j1 = 0; j1 < dim1; j1++ ) {
              for (j2 = 0; j2 < j2_inv[j1]; j2++ ) {
                offset_j2_inv[j1] +=  atoms_ax->bfnnumb_sh[j2];
               }
              }
      Complex *integral_buffer1a, *integral_buffer2a;
      bra_offset = fbz1 * nband_k[2 * s1]     * nband_kq[2 * s1]     * dim1ax;
      ket_offset = fbz2 * nband_k[2 * s2 + 1] * nband_kq[2 * s2 + 1] * dim1ax;
      p1 = 0;
      for (m = 0; m < nband_kq[2 * s1]; m++) {               // v
        for (n = 0; n < nband_k[2 * s1]; n++) {              // c
          p2 = 0;
          for (r = 0; r < nband_kq[2 * s2 + 1]; r++) {       // v'
            for (s = 0; s < nband_k[2 * s2 + 1]; s++) {      // c'
              //fprintf(file.out,"mnrs %3d %3d %3d %3d  occ vir   %3d %3d    %3d %3d %3d %3d\n",\
              m,n,r,s,nocc,nvir,nband_k[array_count],nband_k[array_count+1],nband_kq[array_count],nband_kq[array_count+1]);
              for (j1 = 0; j1 < dim1; j1++) {
                nd6 = atoms_ax->bfnnumb_sh[j1];
                AllocateComplexArray(&integral_buffer1a,&nd6,job);
                AllocateComplexArray(&integral_buffer2a,&nd6,job);
                rotate_single(op1,j1,&integral_buffers[2 * s1][bra_offset + offset_j1_inv[j1] + n * nband_kq[2 * s1] * dim1ax + m * dim1ax],\
                integral_buffer1a,atoms_ax,shells_ax,symmetry,job,file);
                rotate_single(op2,j1,&integral_buffers[2 * s2 + 1][ket_offset + offset_j2_inv[j1] + s * \
                nband_kq[2 * s2 + 1] * dim1ax + r * dim1ax],integral_buffer2a,atoms_ax,shells_ax,symmetry,job,file);
                for (a1 = 0; a1 < nd6; a1++) {
                  //Ham_buffer1[offset + p1 + mpA * p2] += two * integral_buffer1a[a1] * integral_buffer2a[a1] / (double) fermi->nktot;
                  Ham_buffer1[offset + p1 + mpA * p2] += (double) job->spin_fac * integral_buffer1a[a1] * integral_buffer2a[a1] / \
                  (double) fermi->nktot;

                  //fprintf(file.out,"k %3d m %3d n %3d r %3d s %3d coul %3d %3d %14.8f %14.8f %14.8f\n",\
                  k,m,n,r,s,bra_offset + offset_j1_inv[j1] + m * nband_kq[2 * s1] * dim1ax + n * dim1ax + a1,\
                  offset + p1 + mpA * p2,\
                  (Ham_buffer1[offset + p1 + mpA * p2]).real() * au_to_eV,\
                  (integral_buffer1a[a1]).imag(), (integral_buffer2a[a1]).imag());

                  //fprintf(file.out,"Q2 k i3 %3d %3d %3d j1 %3d m %3d n %3d r %3d s %3d coul %3d %3d %3d %3d %14.8f %14.8f %14.8f\n",\
                  q2,k,k1,j1,m,n,r,s,bra_offset + offset_j1_inv[j1] + n * nband_kq[2 * s1] * dim1ax + m * dim1ax + a1,p1,p2,\
                  offset + p1 + mpA * p2,(Ham_buffer1[offset + p1 + mpA * p2]).real() * au_to_eV,\
                  (integral_buffer1a[a1]).real(), (integral_buffer2a[a1]).real());
                  //(integral_buffers[0][bra_offset + offset_j1_inv[j1] + m * nband_kq[array_count] * dim1ax + n * dim1ax + a1]).real(),\
                  (integral_buffers[1][ket_offset + offset_j2_inv[j1] + s * nband_kq[array_count + 1] * dim1ax + r * dim1ax + a1]).real());
                 } 
                DestroyComplexArray(&integral_buffer1a,&nd6,job);
                DestroyComplexArray(&integral_buffer2a,&nd6,job);
               } 
              p2++;
             } 
            } 
           p1++;
          } 
         } 
        } // close if (job->taskid == cblacs_taskid)
       } // close loop on k1
      count_bz += kq_pair.num[k];
     } // close loop on k
    } // close loop on q2
//}
   }  // close loop on s1
  }  // close loop on s2
   time8 += MPI_Wtime() - time7;

  // ******************************************************************************************
  // * Compute Coulomb Energy                                                                 *
  // ******************************************************************************************
  time9 = MPI_Wtime();
  double coulomb_buffer = k_zero, coulomb_energy;
  Complex *integral_buffer1a, *integral_buffer2a;

  if      (job->spin_dim == 1) array_count = 2;
  else if (job->spin_dim == 2) array_count = 8;
  
  if (job->bse_cou == 1)
  for (s1 = 0; s1 < job->spin_dim; s1++) {
    for (s2 = 0; s2 < job->spin_dim; s2++) {
      for (q2 = 0; q2 < fermi->nkunique; q2++) {
        generate_kq_pairs(q2,&kq_pair,crystal,symmetry,fermi,job,file);
        count_bz = 0;
        for (k = 0; k < kq_pair.unique; k++) {
          for (k1 = 0; k1 < kq_pair.num[k]; k1++) {
            k_bz  = kq_pair.bz1[count_bz + k1];
            kq_bz = kq_pair.bz2[count_bz + k1];
            fbz1 = fermi->knet->fbz[k_bz];
            fbz2 = fermi->knet->fbz[kq_bz];
            op1 = symmetry->inverse[fermi->knet->opr[k_bz]];
            op2 = symmetry->inverse[fermi->knet->opr[kq_bz]];
            int j1_inv[dim1];
            int j2_inv[dim1];
            int offset_j1_inv[dim1];
            int offset_j2_inv[dim1];
            for (j1 = 0; j1 < dim1;j1++) j1_inv[j1] = -1;
            for (j2 = 0; j2 < dim1;j2++) j2_inv[j2] = -1;
            for (j1 = 0; j1 < dim1;j1++) offset_j1_inv[j1] = 0;
            for (j2 = 0; j2 < dim1;j2++) offset_j2_inv[j2] = 0;
            for (j1 = 0; j1 < atoms->number_of_atoms_in_unit_cell; j1++ ) {
              for (j2 = 0; j2 < atoms->number_of_atoms_in_unit_cell; j2++ ) {
                if (atom_p->K[j2 * symmetry->number_of_operators + op1] == j1) j1_inv[j1] = j2;
                if (atom_p->K[j2 * symmetry->number_of_operators + op2] == j1) j2_inv[j1] = j2;
               }
              }
            for (j1 = 0; j1 < dim1; j1++ ) {
              for (j2 = 0; j2 < j1_inv[j1]; j2++ ) {
                offset_j1_inv[j1] +=  atoms_ax->bfnnumb_sh[j2];
               }
              }
            for (j1 = 0; j1 < dim1; j1++ ) {
              for (j2 = 0; j2 < j2_inv[j1]; j2++ ) {
                offset_j2_inv[j1] +=  atoms_ax->bfnnumb_sh[j2];
               }
              }

  bra_offset = fbz1 * nband_k[array_count + 2 * s1] * nband_kq[array_count + 2 * s1] * dim1ax;
  ket_offset = fbz2 * nband_k[array_count + 2 * s2 + 1] * nband_kq[array_count + 2 * s2 + 1] * dim1ax;
  for (m = 0; m < nband_k[array_count + 2 * s1]; m++) {
    for (r = 0; r < nband_kq[array_count + 2 * s2 + 1]; r++) {
      for (j1 = 0; j1 < dim1; j1++) {
        nd6 = atoms_ax->bfnnumb_sh[j1];
        AllocateComplexArray(&integral_buffer1a,&nd6,job);
        AllocateComplexArray(&integral_buffer2a,&nd6,job);
        rotate_single(op1,j1,&integral_buffers[array_count + 2 * s1][bra_offset + offset_j1_inv[j1] + m * \
        nband_kq[array_count + 2 * s1] * dim1ax + m * dim1ax],integral_buffer1a,atoms_ax,shells_ax,symmetry,job,file);
        rotate_single(op2,j1,&integral_buffers[array_count + 2 * s2 + 1][ket_offset + offset_j2_inv[j1] + r * \
        nband_kq[array_count + 2 * s2 + 1] * dim1ax + r * dim1ax],integral_buffer2a,atoms_ax,shells_ax,symmetry,job,file);
        //rotate_single(op1,j1,\
        &integral_buffers[array_count][bra_offset + offset_j1_inv[j1] + m * nband_kq[array_count] * dim1ax + m * dim1ax],\
        integral_buffer1a,atoms_ax,shells_ax,symmetry,job,file);
        //rotate_single(op2,j1,\
        &integral_buffers[array_count + 1][ket_offset + offset_j2_inv[j1] + r * nband_kq[array_count + 1] * dim1ax + r * dim1ax],\
        integral_buffer2a,atoms_ax,shells_ax,symmetry,job,file);
        //rotate_single(op1, j1, &integral_buffers[2][bra_offset + offset_j1_inv[j1] + m * nband_kq[array_count] * dim1ax + m * dim1ax],\
        integral_buffer1a,atoms_ax,shells_ax,symmetry,job,file);
        //rotate_single(op2, j1, &integral_buffers[3][ket_offset + offset_j2_inv[j1] + r * nband_kq[array_count + 1] * dim1ax + r * dim1ax],\
        integral_buffer2a,atoms_ax,shells_ax,symmetry,job,file);
        for (a1 = 0; a1 < nd6; a1++) {
          coulomb_buffer += job->spin_fac / (double) (fermi->knet->nktot * fermi->knet->nktot) * \
          (integral_buffer1a[a1] * integral_buffer2a[a1]).real();
          //fprintf(file.out,"Q2 k,k1 j1 %3d %3d %3d %3d m %3d r %3d bra,ket %3d %3d coul %14.8f %14.8f  %14.8f  %14.8f\n",\
          q2,k,k1,j1,m,r,bra_offset + offset_j1_inv[j1] + m * nband_kq[array_count + 2 * s1] * dim1ax + m * dim1ax + a1,\
          ket_offset + offset_j2_inv[j1] + r * nband_kq[array_count + 2 * s2 + 1] * dim1ax + r * dim1ax + a1,\
          coulomb_buffer,two / (double) (fermi->knet->nktot * fermi->knet->nktot) * \
          (integral_buffer1a[a1] * integral_buffer2a[a1]).real(),\
          (integral_buffers[array_count + 2 * s1][bra_offset + offset_j1_inv[j1] + m * nband_kq[array_count + 2 * s1] * dim1ax + \
          m * dim1ax + a1]).real(),\
          (integral_buffers[array_count + 2 * s2 + 1][ket_offset + offset_j2_inv[j1] + r * nband_kq[array_count + 2 * s2 + 1] * dim1ax + \
          r * dim1ax + a1]).real());
         } 
          DestroyComplexArray(&integral_buffer1a,&nd6,job);
          DestroyComplexArray(&integral_buffer2a,&nd6,job);
         } // close loop on j1
        } // close loop on r 
       } // close loop on m
      } // close loop on k1
     count_bz += kq_pair.num[k];
    } // close loop on k
   }  // close loop on q2
   //array_count2 += 2;
  }  // close loop on s2
  //array_count1 += 2;
 }  // close loop on s1
  time10 += MPI_Wtime() - time9;

  time16 += MPI_Wtime() - time15;
  for (s = 0; s < job->band_dim; s++) 
  DestroyComplexArray(&integral_buffers[s],&buffer_size[s],job);
  free_SYMMETRY(&symmetry_little_q_group,job);
  free_k_points(&knet_little_q_group,job);
  free_kq_pairs(&kq_pair,job);
  time14 += MPI_Wtime() - time13;

  time21 = MPI_Wtime();
  time22 += MPI_Wtime() - time21;
  if (job->taskid == 0 && job->bse_cou == 1) printf("coulomb energy                       %16.8f\n",coulomb_buffer / (int)job->spin_dim); 

  DestroyComplexMatrix(&V_q,job);
  free_Q_LATTICE(&q_G,job);
  free_PAIR_TRAN(&pair_p,job);
  DestroyIntArray(&dim_ham,&job->numtasks,job);
  DestroyDoubleArray(&GW_eigenvalues,&dimk,job);
  DestroyDoubleArray(&scf_eigenvalues,&dimk,job);
  MPI_File_close(&fh);

  time2 += MPI_Wtime() - time1;

  if (job->taskid == 1) printf("Total Time                                 %10.2f\n",time2);
  if (job->taskid == 1) printf("Recv                                       %10.2f\n",time12);
  if (job->taskid == 1) printf("Overlap                                    %10.2f\n",time20);
  if (job->taskid == 1) printf("q1 loop                                    %10.2f\n",time14);
  if (job->taskid == 1) printf("Little k group                             %10.2f\n",time18);
  if (job->taskid == 1) printf("Invert Coulomb matrix                      %10.2f\n",time6);
  if (job->taskid == 1) printf("Integrals                                  %10.2f\n",time4);
  if (job->taskid == 1) printf("Contract integrals                         %10.2f\n",time8);
  if (job->taskid == 1) printf("k loop                                     %10.2f\n",time16);
  if (job->taskid == 1) printf("Send                                       %10.2f\n",time10);
  if (job->taskid == 1) printf("Windup                                     %10.2f\n",time22);

}

void setup_hamiltonian_parameters_zero_q(int *ntransitions, int *band_range_k, int *band_range_kq, int *nband_k, int *nband_kq, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, j, k, l, s;
//int coulomb = 1, exchange = 1;
int spin_fac, bse_tda_fac, mult_fac, tmp_size = 12;
int band_k_tmp[2 * tmp_size], band_kq_tmp[2 * tmp_size];
int count1, count2;

  if (job->spin_polarisation == 0)      spin_fac = 1;    // unpolarised
  else if (job->spin_polarisation == 1) spin_fac = 4;    // polarised
  if (job->bse_tda == 0)             bse_tda_fac = 2;    // full BSE/TDHF
  else if (job->bse_tda == 1)        bse_tda_fac = 1;    // TDA BSE/TDHF
  if (job->bse_spin == 0)               mult_fac = 1;    // singlet
  else if (job->bse_spin == 1)          mult_fac = 0;    // triplet
  //*band_range_dim = 2 * spin_fac * bse_tda_fac * mult_fac;
  job->band_dim = 2 * spin_fac * bse_tda_fac * mult_fac;
  //if (coulomb == 1) *band_range_dim += spin_fac;
  if (job->bse_cou == 1 && job->spin_dim == 2) job->band_dim += 2 * job->spin_dim;
  if (job->bse_cou == 1 && job->spin_dim == 1) job->band_dim += 2;
  //printf("band %d cou %d\n",job->band_dim,job->bse_cou);

  // A block q == 0 (cv|v'c') = (cv|conj(c'v'|

  band_k_tmp[0] = fermi->homo[0];                 // c up
  band_k_tmp[1] = fermi->bands[1];

  band_k_tmp[2] = fermi->homo[0];                 // c' up
  band_k_tmp[3] = fermi->bands[1];

  band_k_tmp[4] = fermi->homo[1];                 // c down
  band_k_tmp[5] = fermi->bands[3];

  band_k_tmp[6] = fermi->homo[1];                 // c' down
  band_k_tmp[7] = fermi->bands[3];


  band_kq_tmp[0] = fermi->bands[0] - 1;           // v up
  band_kq_tmp[1] = fermi->homo[0];

  band_kq_tmp[2] = fermi->bands[0] - 1;           // v' up
  band_kq_tmp[3] = fermi->homo[0];

  band_kq_tmp[4] = fermi->bands[2] - 1;           // v down
  band_kq_tmp[5] = fermi->homo[1];

  band_kq_tmp[6] = fermi->bands[2] - 1;           // v' down
  band_kq_tmp[7] = fermi->homo[1];

  // B block q == 0 (cv|c'v') = (cv|conj(v'c'|

  band_k_tmp[8] = fermi->homo[0];                 // c up
  band_k_tmp[9] = fermi->bands[1];

  band_k_tmp[10] = fermi->bands[0] - 1;           // v' up
  band_k_tmp[11] = fermi->homo[0];

  band_k_tmp[12] = fermi->homo[1];                // c down
  band_k_tmp[13] = fermi->bands[3];

  band_k_tmp[14] = fermi->bands[2] - 1;           // v' down
  band_k_tmp[15] = fermi->homo[1];


  band_kq_tmp[8] = fermi->bands[0] - 1;           // v up
  band_kq_tmp[9] = fermi->homo[0];

  band_kq_tmp[10] = fermi->homo[0];               // c' up
  band_kq_tmp[11] = fermi->bands[1];

  band_kq_tmp[12] = fermi->bands[2] - 1;          // v down
  band_kq_tmp[13] = fermi->homo[1];

  band_kq_tmp[14] = fermi->homo[1];               // c' down
  band_kq_tmp[15] = fermi->bands[3];
 
  // coulomb q == 0 (vv|v'v') = (vv|conj(v'v'|

  band_k_tmp[16] = fermi->bands[0] - 1;           // v up
  band_k_tmp[17] = fermi->homo[0];               

  band_k_tmp[18] = fermi->bands[0] - 1;           // v' up
  band_k_tmp[19] = fermi->homo[0];               

  band_k_tmp[20] = fermi->bands[2] - 1;           // v down
  band_k_tmp[21] = fermi->homo[1];               

  band_k_tmp[22] = fermi->bands[2] - 1;           // v' down
  band_k_tmp[23] = fermi->homo[1];               


  band_kq_tmp[16] = fermi->bands[0] - 1;          // v up
  band_kq_tmp[17] = fermi->homo[0];               

  band_kq_tmp[18] = fermi->bands[0] - 1;          // v' up
  band_kq_tmp[19] = fermi->homo[0];              

  band_kq_tmp[20] = fermi->bands[2] - 1;          // v down
  band_kq_tmp[21] = fermi->homo[1];               

  band_kq_tmp[22] = fermi->bands[2] - 1;          // v' down
  band_kq_tmp[23] = fermi->homo[1];                 

  count1 = 0;
  count2 = 0;

  for (i = 0; i < 2 - job->bse_tda; i++) {
    if (i == 1) count2 = 8;
    for (j = 0; j < job->spin_dim; j++) {
      for (k = 0; k < job->spin_dim; k++) {
        for (l = 0; l < 4; l++) {
          band_range_k[count1] = band_k_tmp[count2];
          band_range_kq[count1] = band_kq_tmp[count2];
          //printf("Count1 %3d %3d %3d %3d\n",count1,count2,band_range_k[count1],band_range_kq[count1]);
          count1++;
          count2++;
         }
        }
       }
      }

  if (job->bse_cou == 1) {
    count2 = 16;
    for (j = 0; j < job->spin_dim; j++) {
      for (k = 0; k < job->spin_dim; k++) {
        for (l = 0; l < 2 * job->spin_fac; l++) {
          band_range_k[count1] = band_k_tmp[count2];
          band_range_kq[count1] = band_kq_tmp[count2];
          //printf("count1 %3d %3d %3d %3d\n",count1,count2,band_range_k[count1],band_range_kq[count1]);
          count1++;
          count2++;
         }  
        }  
       }   
      }

  for (i = 0; i < job->band_dim; i++) {
    nband_k[i]  = band_range_k[2 * i + 1]  - band_range_k[2 * i];
    nband_kq[i] = band_range_kq[2 * i + 1] - band_range_kq[2 * i];
    //printf("%3d %3d %3d\n",i,nband_k[i],nband_kq[i]);
   }

  *ntransitions  = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1) * fermi->nktot;
  if (job->spin_dim == 2) 
  *ntransitions += (fermi->bands[3] - fermi->homo[1]) * (fermi->homo[1] - fermi->bands[2] + 1) * fermi->nktot;

}

void setup_hamiltonian_parameters_finite_q(int *ntransitions, int *band_range_k, int *band_range_kq, int *nband_k, int *nband_kq, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int i, j, k, l, s;
int spin_fac, bse_tda_fac, mult_fac;
int band_k_tmp[24], band_kq_tmp[24];
int count1, count2;

  //job->bse_exc = 1;
  if (job->spin_polarisation == 0)      spin_fac = 1;    // unpolarised
  else if (job->spin_polarisation == 1) spin_fac = 2;    // polarised
  if (job->bse_tda == 0)             bse_tda_fac = 2;    // full BSE/TDHF
  else if (job->bse_tda == 1)        bse_tda_fac = 1;    // TDA BSE/TDHF
  if (job->bse_spin == 0)               mult_fac = 1;    // singlet
  else if (job->bse_spin == 1)          mult_fac = 1;    // triplet
  job->band_dim = 2 * spin_fac * bse_tda_fac * mult_fac;
  if (job->bse_exc == 1) job->band_dim += 2 * spin_fac;

  // A block q != 0 (cc'|v'v) = (cc'|conj(vv'|

  band_k_tmp[0] = fermi->homo[0];                 // c up
  band_k_tmp[1] = fermi->bands[1];

  band_k_tmp[2] = fermi->bands[0] - 1;            // v up
  band_k_tmp[3] = fermi->homo[0];

  band_k_tmp[4] = fermi->homo[1];                 // c down
  band_k_tmp[5] = fermi->bands[3];

  band_k_tmp[6] = fermi->bands[2] - 1;            // v down
  band_k_tmp[7] = fermi->homo[1];


  band_kq_tmp[0] = fermi->homo[0];                // c' up
  band_kq_tmp[1] = fermi->bands[1];

  band_kq_tmp[2] = fermi->bands[0] - 1;           // v' up
  band_kq_tmp[3] = fermi->homo[0];

  band_kq_tmp[4] = fermi->homo[1];                // c' down
  band_kq_tmp[5] = fermi->bands[3];

  band_kq_tmp[6] = fermi->bands[2] - 1;           // v' down
  band_kq_tmp[7] = fermi->homo[1];

  // B block q != 0 (cv'|c'v) = (cv'|conj(vc'|

  band_k_tmp[8] = fermi->homo[0];                 // c up
  band_k_tmp[9] = fermi->bands[1];

  band_k_tmp[10] = fermi->bands[0] - 1;           // v up
  band_k_tmp[11] = fermi->homo[0];

  band_k_tmp[12] = fermi->homo[1];                // c down
  band_k_tmp[13] = fermi->bands[3];

  band_k_tmp[14] = fermi->bands[2] - 1;           // v down
  band_k_tmp[15] = fermi->occupied[1];


  band_kq_tmp[8] = fermi->bands[0] - 1;           // v' up
  band_kq_tmp[9] = fermi->homo[0];

  band_kq_tmp[10] = fermi->homo[0];               // c' up
  band_kq_tmp[11] = fermi->bands[1];

  band_kq_tmp[12] = fermi->bands[2] - 1;          // v' down
  band_kq_tmp[13] = fermi->homo[1];

  band_kq_tmp[14] = fermi->homo[1];               // c' down
  band_kq_tmp[15] = fermi->bands[3];
 
  // exchange q != 0 (vv'|v'v) = (vv'|conj(vv'|

  band_k_tmp[16] = fermi->bands[0] - 1;           // v up
  band_k_tmp[17] = fermi->homo[0];                 

  band_k_tmp[18] = fermi->bands[0] - 1;           // v up
  band_k_tmp[19] = fermi->homo[0];                 

  band_k_tmp[20] = fermi->bands[2] - 1;           // v down
  band_k_tmp[21] = fermi->homo[1];                 

  band_k_tmp[22] = fermi->bands[2] - 1;           // v down
  band_k_tmp[23] = fermi->homo[1];                 


  band_kq_tmp[16] = fermi->bands[0] - 1;          // v' up
  band_kq_tmp[17] = fermi->homo[0];                 

  band_kq_tmp[18] = fermi->bands[0] - 1;          // v' up
  band_kq_tmp[19] = fermi->homo[0];                 

  band_kq_tmp[20] = fermi->bands[2] - 1;          // v' down
  band_kq_tmp[21] = fermi->homo[1];                 

  band_kq_tmp[22] = fermi->bands[2] - 1;          // v' down
  band_kq_tmp[23] = fermi->homo[1];                 

  count1 = 0;
  count2 = 0;

  for (i = 0; i < 2 - job->bse_tda; i++) {
    if (i == 1) count2 = 8;
    for (j = 0; j < job->spin_dim; j++) {
      for (l = 0; l < 4; l++) {
        band_range_k[count1] = band_k_tmp[count2];
        band_range_kq[count1] = band_kq_tmp[count2];
        //printf("Eount1 %3d %3d %3d %3d\n",count1,count2,band_range_k[count1],band_range_kq[count1]);
        count1++;
        count2++;
       }
      }
     }

  if (job->bse_exc == 1) {
    count2 = 16;
    for (j = 0; j < job->spin_dim; j++) {
      for (l = 0; l < 4; l++) {
        band_range_k[count1] = band_k_tmp[count2];
        band_range_kq[count1] = band_kq_tmp[count2];
        //printf("count1 %3d %3d %3d %3d\n",count1,count2,band_range_k[count1],band_range_kq[count1]);
        count1++;
        count2++;
       }  
      }  
     }   

  for (i = 0; i < job->band_dim; i++) {
    nband_k[i]  = band_range_k[2 * i + 1]  - band_range_k[2 * i];
    nband_kq[i] = band_range_kq[2 * i + 1] - band_range_kq[2 * i];
    //printf("%3d %3d %3d %3d %3d %3d %3d\n",\
    i,nband_k[i],nband_kq[i],fermi->bands[0],fermi->bands[2],fermi->occupied[0],fermi->occupied[fermi->nkunique]);
   }

  *ntransitions  = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1) * fermi->nktot;
  if (job->spin_dim == 2) 
  *ntransitions += (fermi->bands[3] - fermi->homo[1]) * (fermi->homo[1] - fermi->bands[2] + 1) * fermi->nktot;

}

void diagonalise_bse_tda_hamiltonian_crystal(int *ictxt, int *nbsize, Complex *Ham_buffer1, FERMI* fermi, JOB_PARAM *job, FILES file)

{

int ntransitions;
int band_range_k[24], band_range_kq[24], nband_k[12], nband_kq[12];
int buffer_size;
int nprow, npcol, myrow, mycol, mpA, nqA, izero = 0, ione = 1;
int itemp, descA[9];
int info = 0;
int lwork = -1;
double time1, time2, time3, time4;
double *bse_eigenvalues;
double *rwork;
char jobz = 'V';
char uplo = 'U';
char xx[20] = "/bse_eigenvectors";
char yy[20] = "/bse_eigenvalues";
FILE *bse_evalues;
Complex *bse_eigvec, *work;

  setup_hamiltonian_parameters_zero_q(&ntransitions, band_range_k, band_range_kq, nband_k, nband_kq, fermi, job, file);
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(&ntransitions, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(&ntransitions, nbsize, &mycol, &izero, &npcol);
  itemp = max(1, mpA);
  descinit_(descA, &ntransitions, &ntransitions, nbsize, nbsize, &izero, &izero, ictxt, &itemp, &info);

  buffer_size = mpA * nqA;
  AllocateDoubleArray(&bse_eigenvalues,&ntransitions,job);
  AllocateComplexArray(&bse_eigvec,&buffer_size,job);
  ResetDoubleArray(bse_eigenvalues,&ntransitions);
  ResetComplexArray(bse_eigvec,&buffer_size);

  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0)) {
int i, j;
//for (i = 0; i < 8; i++) {
//for (j = 0; j < 8; j++) {
//fprintf(file.out,"%6.3lf ", au_to_eV*(Ham_buffer1[i * ntransitions + j]).real());}fprintf(file.out,"\n");}

 // fprintf(file.out,"\n"); 
 // fprintf(file.out,"\n"); 
 // for (i = 0; i < 6; i++) {
 // for (j = 0; j < 6; j++) {
 // fprintf(file.out,"%6.3lf", (Ham_buffer1[i + mpA * j]).real());
 //}fprintf(file.out,"\n"); }
 // fprintf(file.out,"\n"); 
 // fprintf(file.out,"\n");
 //
 // fprintf(file.out,"\n"); 
 // for (int k = 0; k < 8; k++) {
 // fprintf(file.out,"\n"); 
 // for (i = 0; i < 6; i++) {
 // for (j = 0; j < 6; j++) {
 // fprintf(file.out,"%6.3lf", (Ham_buffer1[k * 6 + i + mpA * j]).real());
 //}fprintf(file.out,"\n"); }
 // fprintf(file.out,"\n"); }
 // fprintf(file.out,"\n");
 //
 //
 //
 // fprintf(file.out,"\n"); 
 // fprintf(file.out,"\n"); 
 // fprintf(file.out,"\n"); 
 // for (i = 0; i < mpA; i++) {
 // for (j = 0; j < 16; j++) {
 // fprintf(file.out,"%8.4lf", (Ham_buffer1[i + mpA * j]).real());
 //}fprintf(file.out,"\n"); }
 // fprintf(file.out,"\n");
 // fprintf(file.out,"\n");
 //
 // fprintf(file.out,"\n"); 
 // fprintf(file.out,"\n"); 
 // for (i = 0; i < mpA; i++) {
 // for (j = 16; j < 32; j++) {
 // fprintf(file.out,"%9.5lf", (Ham_buffer1[i + mpA * j]).real());
 //}fprintf(file.out,"\n"); }
 // fprintf(file.out,"\n");
 // fprintf(file.out,"\n");
 //
 //  fprintf(file.out,"%5.2lf%5.2lf ", (Ham_buffer1[i * ntransitions + j]).real()*au_to_eV,\
 (Ham_buffer1[i * ntransitions + j]).imag()*au_to_eV);
 //}   fprintf(file.out,"\n"); }

  time1 = MPI_Wtime();
  int lwork, lrwork;
  work = (Complex *) calloc(2, sizeof(Complex)) ;
  if (work==NULL) { 
  fprintf(file.out,"ERROR: allocation of work array on core %3d\n",job->taskid); 
  MPI_Finalize();
  exit(1); 
 }
  double *rwork;
  rwork = (double *) calloc(2, sizeof(double)) ;
  if (rwork==NULL) { 
  fprintf(file.out,"ERROR: allocation of rwork array on core %3d\n",job->taskid); 
  MPI_Finalize();
  exit(1); 
 }
  lwork=-1;
  lrwork=-1;
  pzheev_(&jobz,&uplo,&ntransitions,Ham_buffer1,&ione,&ione,descA,bse_eigenvalues,bse_eigvec,&ione,&ione,descA,\
  work,&lwork,rwork,&lrwork,&info);
  lrwork = (int) rwork[0];
  lwork  = (int)  work[0].real();
  free(work);
  free(rwork);
  work = (MKL_Complex16 *) calloc(lwork, sizeof(MKL_Complex16)) ;
  if (work==NULL) { 
  fprintf(file.out,"ERROR: allocation of work array on core %3d\n",job->taskid); 
  MPI_Finalize();
  exit(1); 
 }
  rwork = (double *) calloc(2 * lrwork, sizeof(double)) ;
  if (rwork==NULL) { 
  fprintf(file.out,"ERROR: allocation of rwork array on core %3d\n",job->taskid); 
  MPI_Finalize();
  exit(1); 
 }

  pzheev_(&jobz,&uplo,&ntransitions,Ham_buffer1,&ione,&ione,descA,bse_eigenvalues,bse_eigvec,&ione,&ione,descA,\
  work,&lwork,rwork,&lrwork,&info);
  free(work);
  free(rwork);
  time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("pzheev spk                                 %10.2f\n",time2);
  } // close if (

  if (job->taskid == 0) printf("ntransitions %6d bse_lim %6d\n",ntransitions,job->bse_lim);
  block_cyclic_to_linear_limit_complex(&ntransitions,ictxt,nbsize,job->bse_lim * fermi->nktot,bse_eigvec,xx,job,file);

  char buf1[110];
  strcpy(buf1,file.scf_eigvec);
  strcat(buf1,yy);
  MPI_File hh;
  MPI_File_open(MPI_COMM_WORLD,buf1,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&hh) ;
  MPI_File_seek(hh, 0, MPI_SEEK_SET) ;
  MPI_File_write(hh,bse_eigenvalues,ntransitions,MPI_DOUBLE,MPI_STATUS_IGNORE);
  MPI_File_close(&hh);

  //if (job->taskid == 0) {
  //bse_evalues  = fopen(yy, "wb");
  //fwrite(bse_eigenvalues, sizeof(double), ntransitions, bse_evalues);
  //fflush(bse_evalues);
  //fclose(bse_evalues);
 //}
 
  //DestroyDoubleArray(&bse_eigenvalues,&nt,job);
  //DestroyComplexArray(&bse_eigvec,&buffer_size,job);

}

void overlap_matrices5(ComplexMatrix *S1, ComplexMatrix *S2, ComplexMatrix *scf_eigenvectors_k, ComplexMatrix *scf_eigenvectors_kq, KQPOINT_TRAN *kq_pair, int *s1, int *band_offset1, int *band_offset2, int *nband_k, int *nband_kq, int *count_bz, INT_1E *one_ints, PAIR_TRAN pair_p, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, REAL_LATTICE *R, FERMI *fermi, SYMMETRY *symmetry, CRYSTAL *crystal, MPI_File *fh, FILES file, JOB_PARAM *job)

{

int i, j, k;
int op, k_bz, fbz1;
int vector_size_k, vector_size_kq;
int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
int nk[2];
double time1, time2, time3, time4, time5, time6;
char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
KPOINT_TRAN knet_tmp;
ComplexMatrix *scf_eigenvectors_0, *scf_eigenvectors_1;
ComplexMatrix *scf_eigenvectors_kp, *scf_eigenvectors_kpp;
ComplexMatrix *scf_eigenvectors_kqp, *scf_eigenvectors_kqpp;
ComplexMatrix *tmp_k, *tmp_kq, *S_k;
Complex alpha = Complex(k_one, k_zero);
Complex  beta = Complex(k_zero, k_zero);
//Complex alpha, beta;
//alpha.real() = k_one;
//alpha.imag() = k_zero;
//beta.real() = k_zero;
//beta.imag() = k_zero;

  time1 = MPI_Wtime();

  knet_tmp.nktot = 1;
  knet_tmp.unique = 1;
  vector_size_k  = *nband_k  * nbfn;
  vector_size_kq = *nband_kq * nbfn;

  AllocateComplexMatrix(&S_k,&nbfn,&nbfn,job);
  allocate_k_points(&knet_tmp,crystal,job,file);

  k_bz  = kq_pair->bz1[*count_bz];
  fbz1 = fermi->knet->fbz[k_bz];
  op  = symmetry->inverse[kq_pair->opr[*count_bz]];
  int seek_spin_offset = *s1 * fermi->nkunique * nbfn * nbfn;
//if (*s1 == 1) fprintf(file.out,"seek %3d\n",seek_spin_offset);

  nk[0] = k_bz;
  nk[1] = k_bz;
  fourier_transform(&one_ints->Overlap[0], &S_k->a[0][0], fermi->knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);

  knet_tmp.opr[0] = kq_pair->opr[*count_bz];
  knet_tmp.trs[0] = 0; // fermi->knet->trs[k_bz];
  knet_tmp.cart[0].comp1 = fermi->knet->cart[k_bz].comp1;
  knet_tmp.cart[0].comp2 = fermi->knet->cart[k_bz].comp2;
  knet_tmp.cart[0].comp3 = fermi->knet->cart[k_bz].comp3;

  time5 = MPI_Wtime();
  AllocateComplexMatrix(&tmp_k,&nbfn,nband_k,job);
  AllocateComplexMatrix(&scf_eigenvectors_0,nband_k,&nbfn,job);
  AllocateComplexMatrix(&scf_eigenvectors_kp,nband_k,&nbfn,job);
  AllocateComplexMatrix(&scf_eigenvectors_kpp,nband_k,&nbfn,job);
  ResetComplexMatrix(S1);
  ResetComplexMatrix(tmp_k);
  MPI_File_seek(*fh, (seek_spin_offset + fbz1 * nbfn * nbfn + *band_offset1 * nbfn) * sizeof(Complex), MPI_SEEK_SET);
  //MPI_File_seek(*fh,  (seek_spin_offset + fbz1 * nbfn * nbfn  + (fermi->bands[2 * *s1] - 1) * nbfn) * sizeof(Complex), MPI_SEEK_SET);
  //MPI_File_seek(*fh,  (seek_spin_offset + fbz1 * nbfn * nbfn  + (fermi->homo[0]) * nbfn) * sizeof(Complex), MPI_SEEK_SET);
  MPI_File_read(*fh, &scf_eigenvectors_0->a[0][0], 2 * vector_size_k, MPI_DOUBLE, MPI_STATUS_IGNORE);
  rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_kp->a[0][0],*nband_k,k_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
  rotate_psi(&scf_eigenvectors_k->a[0][0],&scf_eigenvectors_kpp->a[0][0],*nband_k,0,&knet_tmp,atom_p,atoms,R,shells,symmetry,job,file);
  ComplexGEMM3(&NoTrans,&ConjTrans,&nbfn,nband_k,&nbfn,&alpha,&(S_k->a[0]),&nbfn,&(scf_eigenvectors_kpp->a[0]),&nbfn,&beta,&(tmp_k->a[0]),\
  nband_k);
  ComplexGEMM3(&NoTrans,&NoTrans,nband_k,nband_k,&nbfn,&alpha,&(scf_eigenvectors_kp->a[0]),&nbfn,&(tmp_k->a[0]),nband_k,&beta,&(S1->a[0]),\
  nband_k);
//if (*s1 == 2) {
//fprintf(file.out,"%3d %3d %3d\n",fbz1,vector_size_k,seek_spin_offset + fbz1 * nbfn * nbfn + *band_offset1 * nbfn);
//print_complex_matrix(scf_eigenvectors_0,file);
//print_complex_matrix(scf_eigenvectors_k,file);
//print_complex_matrix(scf_eigenvectors_kp,file);
//print_complex_matrix(scf_eigenvectors_kpp,file);
//}

  DestroyComplexMatrix(&tmp_k,job);
  DestroyComplexMatrix(&scf_eigenvectors_0,job);
  DestroyComplexMatrix(&scf_eigenvectors_kp,job);
  DestroyComplexMatrix(&scf_eigenvectors_kpp,job);

  AllocateComplexMatrix(&tmp_kq,&nbfn,nband_kq,job);
  AllocateComplexMatrix(&scf_eigenvectors_1,nband_kq,&nbfn,job);
  AllocateComplexMatrix(&scf_eigenvectors_kqp,nband_kq,&nbfn,job);
  AllocateComplexMatrix(&scf_eigenvectors_kqpp,nband_kq,&nbfn,job);
  ResetComplexMatrix(S2);
  ResetComplexMatrix(tmp_kq);
  rotate_psi(&scf_eigenvectors_kq->a[0][0],&scf_eigenvectors_kqpp->a[0][0],*nband_kq,0,&knet_tmp,atom_p,atoms,R,shells,symmetry,job,file);
  MPI_File_seek(*fh,  (seek_spin_offset + fbz1 * nbfn * nbfn + *band_offset2 * nbfn) * sizeof(Complex), MPI_SEEK_SET) ;
  //MPI_File_seek(*fh,  (seek_spin_offset + fbz1 * nbfn * nbfn  + (fermi->bands[2 * *s1] - 1) * nbfn) * sizeof(Complex), MPI_SEEK_SET) ;
  MPI_File_read(*fh, &scf_eigenvectors_1->a[0][0], 2 * vector_size_kq, MPI_DOUBLE, MPI_STATUS_IGNORE);
  rotate_psi(&scf_eigenvectors_1->a[0][0],&scf_eigenvectors_kqp->a[0][0],*nband_kq,k_bz,fermi->knet,atom_p,atoms,R,shells,\
  symmetry,job,file);
  ComplexGEMM3(&NoTrans,&ConjTrans,&nbfn,nband_kq,&nbfn,&alpha,&(S_k->a[0]),&nbfn,&(scf_eigenvectors_kqpp->a[0]),&nbfn,&beta,\
  &(tmp_kq->a[0]),nband_kq);
  ComplexGEMM3(&NoTrans,&NoTrans,nband_kq,nband_kq,&nbfn,&alpha,&(scf_eigenvectors_kqp->a[0]),&nbfn,&(tmp_kq->a[0]),nband_kq,&beta,\
  &(S2->a[0]),nband_kq);
  DestroyComplexMatrix(&tmp_kq,job);
  DestroyComplexMatrix(&scf_eigenvectors_1,job);
  DestroyComplexMatrix(&scf_eigenvectors_kqp,job);
  DestroyComplexMatrix(&scf_eigenvectors_kqpp,job);
  DestroyComplexMatrix(&S_k,job);
  free_k_points(&knet_tmp,job);

}

void overlap_matrices6(ComplexMatrix *S1, ComplexMatrix *S2, ComplexMatrix *scf_eigenvectors_k, ComplexMatrix *scf_eigenvectors_kq, KQPOINT_TRAN *kq_pair, int *s1, int *band_offset1, int *band_offset2, int *nband_k, int *nband_kq, int *count_bz, INT_1E *one_ints, PAIR_TRAN pair_p, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, REAL_LATTICE *R, FERMI *fermi, SYMMETRY *symmetry, CRYSTAL *crystal, MPI_File *fh, FILES file, JOB_PARAM *job)

{

int i, j, k;
int op, kq_bz, fbz2;
int vector_size_k, vector_size_kq;
int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
int nk[2];
char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
KPOINT_TRAN knet_tmp;
ComplexMatrix *scf_eigenvectors_0, *scf_eigenvectors_1;
ComplexMatrix *scf_eigenvectors_kp, *scf_eigenvectors_kpp;
ComplexMatrix *scf_eigenvectors_kqp, *scf_eigenvectors_kqpp;
ComplexMatrix *tmp_k, *tmp_kq, *S_k;
Complex alpha = Complex(k_one, k_zero);
Complex  beta = Complex(k_zero, k_zero);
//Complex alpha, beta;
//alpha.real() = k_one;
//alpha.imag() = k_zero;
//beta.real() = k_zero;
//beta.imag() = k_zero;

  knet_tmp.nktot = 1;
  knet_tmp.unique = 1;
  //vector_size_k  = nband_k[0]  * nbfn;
  //vector_size_kq = nband_kq[0] * nbfn;
  vector_size_k  = *nband_k  * nbfn;
  vector_size_kq = *nband_kq * nbfn;

  AllocateComplexMatrix(&S_k,&nbfn,&nbfn,job);
  allocate_k_points(&knet_tmp,crystal,job,file);

  kq_bz = kq_pair->bz2[*count_bz];
  fbz2 = fermi->knet->fbz[kq_bz];
  op  = symmetry->inverse[kq_pair->opr[*count_bz]];
  int seek_spin_offset = *s1 * fermi->nkunique * nbfn * nbfn;

  nk[0] = kq_bz;
  nk[1] = kq_bz;
  fourier_transform(&one_ints->Overlap[0], &S_k->a[0][0], fermi->knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
  knet_tmp.opr[0] = kq_pair->opr[*count_bz];
  knet_tmp.trs[0] = 0; // fermi->knet->trs[kq_bz];
  knet_tmp.cart[0].comp1 = fermi->knet->cart[kq_bz].comp1;
  knet_tmp.cart[0].comp2 = fermi->knet->cart[kq_bz].comp2;
  knet_tmp.cart[0].comp3 = fermi->knet->cart[kq_bz].comp3;

  //AllocateComplexMatrix(&tmp_k,&nbfn,&nband_k[0],job);
  //AllocateComplexMatrix(&scf_eigenvectors_0,&nband_k[0],&nbfn,job);
  //AllocateComplexMatrix(&scf_eigenvectors_kp,&nband_k[0],&nbfn,job);
  //AllocateComplexMatrix(&scf_eigenvectors_kpp,&nband_k[0],&nbfn,job);
  AllocateComplexMatrix(&tmp_k,&nbfn,nband_k,job);
  AllocateComplexMatrix(&scf_eigenvectors_0,nband_k,&nbfn,job);
  AllocateComplexMatrix(&scf_eigenvectors_kp,nband_k,&nbfn,job);
  AllocateComplexMatrix(&scf_eigenvectors_kpp,nband_k,&nbfn,job);
  ResetComplexMatrix(S1);
  ResetComplexMatrix(tmp_k);
  MPI_File_seek(*fh,  (seek_spin_offset + fbz2 * nbfn * nbfn + *band_offset1 * nbfn) * sizeof(Complex), MPI_SEEK_SET);
  //MPI_File_seek(*fh,  (seek_spin_offset + fbz2 * nbfn * nbfn  + (fermi->bands[2 * *s1] - 1) * nbfn) * sizeof(Complex), MPI_SEEK_SET);
  //MPI_File_seek(*fh,  (fbz2 * nbfn * nbfn  + (fermi->homo[0]) * nbfn) * sizeof(Complex), MPI_SEEK_SET);
  MPI_File_read(*fh, &scf_eigenvectors_0->a[0][0], 2 * vector_size_k, MPI_DOUBLE, MPI_STATUS_IGNORE);
  rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_kp->a[0][0],*nband_k,kq_bz,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
  rotate_psi(&scf_eigenvectors_k->a[0][0],&scf_eigenvectors_kpp->a[0][0],*nband_k,0,&knet_tmp,atom_p,atoms,R,shells,symmetry,job,file);
  ComplexGEMM3(&NoTrans,&ConjTrans,&nbfn,nband_k,&nbfn,&alpha,&(S_k->a[0]),&nbfn,&(scf_eigenvectors_kp->a[0]),&nbfn,&beta,&(tmp_k->a[0]),\
  nband_k);
  ComplexGEMM3(&NoTrans,&NoTrans,nband_k,nband_k,&nbfn,&alpha,&(scf_eigenvectors_kpp->a[0]),&nbfn,&(tmp_k->a[0]),nband_k,&beta,&(S1->a[0]),\
  nband_k);
  DestroyComplexMatrix(&tmp_k,job);
  DestroyComplexMatrix(&scf_eigenvectors_0,job);
  DestroyComplexMatrix(&scf_eigenvectors_kp,job);

  //AllocateComplexMatrix(&tmp_kq,&nbfn,&nband_kq[0],job);
  //AllocateComplexMatrix(&scf_eigenvectors_1,&nband_kq[0],&nbfn,job);
  //AllocateComplexMatrix(&scf_eigenvectors_kqp,&nband_kq[0],&nbfn,job);
  //AllocateComplexMatrix(&scf_eigenvectors_kqpp,&nband_kq[0],&nbfn,job);
  AllocateComplexMatrix(&tmp_kq,&nbfn,nband_kq,job);
  AllocateComplexMatrix(&scf_eigenvectors_1,nband_kq,&nbfn,job);
  AllocateComplexMatrix(&scf_eigenvectors_kqp,nband_kq,&nbfn,job);
  AllocateComplexMatrix(&scf_eigenvectors_kqpp,nband_kq,&nbfn,job);
  ResetComplexMatrix(S2);
  ResetComplexMatrix(tmp_kq);
  rotate_psi(&scf_eigenvectors_kq->a[0][0],&scf_eigenvectors_kqpp->a[0][0],*nband_kq,0,&knet_tmp,atom_p,atoms,R,shells,symmetry,job,file);
  MPI_File_seek(*fh,  (seek_spin_offset + fbz2 * nbfn * nbfn + *band_offset2 * nbfn) * sizeof(Complex), MPI_SEEK_SET) ;
  //MPI_File_seek(*fh,  (seek_spin_offset + fbz2 * nbfn * nbfn  + (fermi->bands[2 * *s1] - 1) * nbfn) * sizeof(Complex), MPI_SEEK_SET) ;
  MPI_File_read(*fh, &scf_eigenvectors_1->a[0][0], 2 * vector_size_kq, MPI_DOUBLE, MPI_STATUS_IGNORE);
  rotate_psi(&scf_eigenvectors_1->a[0][0],&scf_eigenvectors_kqp->a[0][0],*nband_kq,kq_bz,fermi->knet,atom_p,atoms,R,shells,\
  symmetry,job,file);
  ComplexGEMM3(&NoTrans,&ConjTrans,&nbfn,nband_kq,&nbfn,&alpha,&(S_k->a[0]),&nbfn,&(scf_eigenvectors_kqp->a[0]),&nbfn,&beta,\
  &(tmp_kq->a[0]),nband_kq);
  ComplexGEMM3(&NoTrans,&NoTrans,nband_kq,nband_kq,&nbfn,&alpha,&(scf_eigenvectors_kqpp->a[0]),&nbfn,&(tmp_kq->a[0]),nband_kq,&beta,\
  &(S2->a[0]),nband_kq);
  DestroyComplexMatrix(&tmp_kq,job);
  DestroyComplexMatrix(&scf_eigenvectors_1,job);
  DestroyComplexMatrix(&scf_eigenvectors_kqp,job);

  //printf("nband %3d  %3d %3d\n",*s1,*nband_k,*nband_kq);

  DestroyComplexMatrix(&S_k,job);
  free_k_points(&knet_tmp,job);

}

void generate_kq_pairs(int q1, KQPOINT_TRAN *kq_pairs, CRYSTAL *crystal, SYMMETRY *symmetry, FERMI *fermi, JOB_PARAM *job, FILES file)

{

int nd6;
int k, l, n, q, t, j1, j2, k1, t1, kq, rot_index_k, rot_index_q, rot_index_kq;
int offset1;
int ksize, symmetry_number_of_operators;
int is[3];
int *p_inr;
VECTOR_INT kvec[3], qvec[3];
SYMMETRY symmetry_little_q_group;
KPOINT_TRAN knet_little_q_group;

  is[0] = fermi->is[0];
  is[1] = fermi->is[1];
  is[2] = fermi->is[2];
  q = fermi->knet->ibz[q1];
  qvec[0] = decompose_k_point(is, q, crystal, job, file);
  if (job->kss == 0)      symmetry_number_of_operators = 1;
  else if (job->kss == 1) symmetry_number_of_operators = symmetry->number_of_operators;
  knet_size(&ksize,fermi->is,crystal);
  for (j2 = 0; j2 < ksize * symmetry_number_of_operators; j2++) kq_pairs->ibz[j2] = -1;
  for (j2 = 0; j2 < ksize * symmetry_number_of_operators; j2++) kq_pairs->bz1[j2] = -1;
  for (j2 = 0; j2 < ksize * symmetry_number_of_operators; j2++) kq_pairs->bz2[j2] = -1;
  count_little_k_group_operators(q1,&symmetry_little_q_group,symmetry,crystal,fermi->knet,fermi,job,file);
  allocate_SYMMETRY(&symmetry_little_q_group,job,file);
  generate_little_k_group(q1, &symmetry_little_q_group, fermi, fermi->knet, symmetry, crystal, job, file);
  count_k_points(&knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
  allocate_k_points(&knet_little_q_group,crystal,job,file);
  generate_k_points(&knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
  free_SYMMETRY(&symmetry_little_q_group,job);
  kq_pairs->unique = 0;
  for (k = 0; k < knet_little_q_group.unique; k++) {
    kq_pairs->num[k] = 0;
    k1 = knet_little_q_group.ibz[k];
    kvec[0]  = decompose_k_point(is, k1, crystal, job, file);
    for (t1 = 0; t1 <= job->trs; t1++) {
      for (l = 0; l < symmetry_number_of_operators; l++) {
        p_inr = symmetry->inr + symmetry->inverse[l] * 9;
        rotate_vector_int_latt(p_inr, is, t1, &kvec[0], &kvec[1], crystal, job, file);
        //fprintf(file.out,"kvec    %3d %3d %3d   %3d %3d %3d\n",\
        kvec[0].comp1,kvec[0].comp2,kvec[0].comp3,kvec[1].comp1,kvec[1].comp2,kvec[1].comp3);
        rot_index_k = kvec[1].comp1 * is[1] * is[2] + kvec[1].comp2 * is[2] + kvec[1].comp3;
        //rotate_vector_int(p_inr, &qvec[2], &qvec[1]);
        rotate_vector_int_latt(p_inr, is, t1, &qvec[0], &qvec[1], crystal, job, file);
        rot_index_q  =  qvec[1].comp1 * is[1] * is[2] +  qvec[1].comp2 * is[2] +  qvec[1].comp3;
        rot_index_kq = compose_k_point(fermi->is, kvec[1], qvec[1], crystal, job, file);
        //fprintf(file.out,"qvec    %3d %3d %3d   %3d %3d %3d   %3d %3d  op %3d\n",\
        qvec[0].comp1,qvec[0].comp2,qvec[0].comp3,qvec[1].comp1,qvec[1].comp2,qvec[1].comp3,rot_index_k,rot_index_kq,l);
        for (n = 0; n <= kq_pairs->unique; n++) {
          if (rot_index_k == kq_pairs->bz1[n] && rot_index_kq == kq_pairs->bz2[n]) break;
          if (n == kq_pairs->unique && (kq_pairs->bz1[n] != rot_index_k || kq_pairs->bz2[n] != rot_index_kq)) {
          kq_pairs->opr[n] = l;
          //kq_pairs->ibz[n] = k;
          kq_pairs->ibz[n] = rot_index_q;
          kq_pairs->bz1[n] = rot_index_k;
          kq_pairs->bz2[n] = rot_index_kq;
          kq_pairs->num[k]++;
          (kq_pairs->unique)++;
          //printf("q %3d q1 %3d k %3d rot_k %3d rot_kq %3d num_k %3d total %3d\n",\
          q,q1,k,kq_pairs->bz1[n],kq_pairs->bz2[n],kq_pairs->num[k],kq_pairs->unique);
          //fprintf(file.out,"q %3d q1 %3d k %5d ibz %5d rot_k %5d rot_kq %5d opr %5d num_k %3d total %7d\n", \
          q,q1,k,k1,kq_pairs->bz1[n],kq_pairs->bz2[n],kq_pairs->opr[n],kq_pairs->num[k],kq_pairs->unique); fflush(file.out);
          break;
         }
        } // close loop on n
       } // close loop on l
      } // close loop on t1
       //fprintf(file.out,"q1 %3d k %3d num %3d opr %3d total_pairs %7d\n",\
       q1,k,kq_pairs->num[k],kq_pairs->opr[n],knet_little_q_group.unique);
       //fprintf(file.out,"\n");
       //printf("\n");
     } // close loop on k
       kq_pairs->unique = knet_little_q_group.unique;
       free_k_points(&knet_little_q_group,job);

}

void select_kq_pairs(int *num_q1, int *q1_list, int *begin_k, int *end_k, KPOINT_TRAN *knet_little_q_group, FERMI *fermi, SYMMETRY *symmetry, CRYSTAL *crystal, Q_LATTICE *q_G, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k, q, q1;
int flag, total_kq_pairs;
int begin_kq[job->numtasks], end_kq[job->numtasks];
SYMMETRY symmetry_little_q_group;
KQPOINT_TRAN kq_pair;

  for (i = 0; i < 2; i++) { 
    k = 0;
    *num_q1 = 0;
    total_kq_pairs = 0;
    for (q1 = 0; q1 < fermi->knet->unique; q1++) {
      q = fermi->knet->ibz[q1];
      count_little_k_group_operators(q1,&symmetry_little_q_group,symmetry,crystal,fermi->knet,fermi,job,file);
      allocate_SYMMETRY(&symmetry_little_q_group,job,file);
      generate_little_k_group(q1, &symmetry_little_q_group, fermi, fermi->knet, symmetry, crystal, job, file);
      count_k_points(knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
      allocate_k_points(knet_little_q_group,crystal,job,file);
      generate_k_points(knet_little_q_group,fermi->is,crystal,&symmetry_little_q_group,job,file);
      generate_q_lattice(&fermi->knet->oblique[q], q_G, fermi, G, crystal, job, file);
      allocate_kq_pairs(&kq_pair,fermi,crystal,symmetry,job,file);
      generate_kq_pairs(q1,&kq_pair,crystal,symmetry,fermi,job,file);
      flag = 0;
      for (j = 0; j < knet_little_q_group->unique; j++) {
        if (i == 1 && total_kq_pairs + j >= begin_kq[job->taskid] && total_kq_pairs + j < end_kq[job->taskid]) {
        if (flag == 0) {
        begin_k[k] = j;
        q1_list[k] = q1;
        (*num_q1)++;
       }
        end_k[k] = j + 1;
        flag = 1;
       }
      }
       if (flag == 1) k++;
       total_kq_pairs += knet_little_q_group->unique;
      }
       if (i == 0) mpi_begin_end(begin_kq,end_kq,total_kq_pairs,job->numtasks,job,file);
      }
       MPI_Barrier(MPI_COMM_WORLD);
       //for (i = 0; i < *num_q1; i++) printf("%3d  %3d %3d %3d %3d\n",job->taskid,i,q1_list[i],begin_k[i],end_k[i]);

}

void read_scf_GW_eigenvalues_crystal(double *eigenvalues, FERMI *fermi, ATOM *atoms, char *zz, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Read SCF, GW or CAS eigenvalues from disk                                              *
  // ******************************************************************************************

int i, j, s;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int nbands[2];
int seekpoint;
char buf1[110];
MPI_File eh;

  nbands[0] = fermi->bands[1] - fermi->bands[0] + 1;
  nbands[1] = 0;
  if (job->spin_dim == 2) 
  nbands[1] = fermi->bands[3] - fermi->bands[2] + 1;
  //printf("%s %3d %3d\n",zz,fermi->nkunique,nbands);
  //FILE *evals;
  //int i, s;
  //double *evalues_temp;
  //AllocateDoubleArray(&evalues_temp,&datasize,job);
  //ResetDoubleArray(evalues_temp,&datasize);
  //if (job->taskid == 0) {
  //evals = fopen(zz, "rb");
  //for (i = 0; i < fermi->nkunique; i++) {
  //seekpoint = i * dim1 + fermi->bands[0] - 1;
  //fseek(evals,seekpoint*sizeof(double),SEEK_SET);
  //fread(&eigenvalues[i * nbands],sizeof(double),nbands,evals);
 //}
  //fclose(evals);
 //}
  //MPI_Allreduce(evalues_temp,eigenvalues,job->spin_dim * datasize,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //DestroyDoubleArray(&evalues_temp,&datasize,job);

  strcpy(buf1,file.scf_eigvec);
  strcat(buf1,zz);
  //printf("reading from %s\n",buf1);
  //printf("%3d %3d %3d %3d %3d\n",fermi->nkunique,dim1,fermi->bands[0],nbands[0],nbands[1]);
  MPI_File_open(MPI_COMM_WORLD,buf1,MPI_MODE_RDONLY,MPI_INFO_NULL,&eh) ;
  //for (s = 0; s < job->spin_dim; s++) {
  for (i = 0; i < fermi->nkunique; i++) {
  for (s = 0; s < job->spin_dim; s++) {
  //fprintf(file.out,"Eigenvalues %3d %3d\n",dim1,fermi->bands[0]);
  seekpoint = i * job->spin_dim * dim1 + s * dim1 + fermi->bands[2 * s] - 1;
  //fprintf(file.out,"seek %3d %3d %3d %3d\n",i,s,seekpoint,i * (nbands[0] + nbands[1]) + s * nbands[0]);
  //seekpoint = s * fermi->nkunique * dim1 + i * dim1 + fermi->bands[2 * s] - 1;
  //MPI_File_write(eh, &eigval1[(s * knet.unique + k) * dim1 + fermi->bands[0]-1],value_size,MPI_DOUBLE,MPI_STATUS_IGNORE);
  //fprintf(file.out,"%3d %3d %3d %3d %3d %3d\n",s,i,nbands,fermi->nkunique,seekpoint,fermi->bands[2 * s]);
  //seekpoint = i * dim1 + fermi->bands[2 * s] - 1;
  //seekpoint = i * dim1 + fermi->bands[0] - 1;
  //seekpoint = fermi->knet->fbz[i] * dim1 + fermi->bands[0] - 1;
  MPI_File_seek(eh, seekpoint * sizeof(double), MPI_SEEK_SET) ;
  //MPI_File_read(eh,&eigenvalues[i * nbands],nbands,MPI_DOUBLE,MPI_STATUS_IGNORE);
  //MPI_File_read(eh,&eigenvalues[s * fermi->nkunique * nbands[0] + i * nbands[s]],nbands[s],MPI_DOUBLE,MPI_STATUS_IGNORE);
  MPI_File_read(eh,&eigenvalues[i * (nbands[0] + nbands[1]) + s * nbands[0]], nbands[s], MPI_DOUBLE, MPI_STATUS_IGNORE);
 }
 }
  MPI_File_close(&eh);

  if (job->taskid == 0 && job->verbosity >= 1) {
  fprintf(file.out,"Eigenalues %3d \n",job->spin_dim);
  //for (s = 0; s < job->spin_dim; s++) {
  for (i = 0; i < fermi->nkunique; i++) {
  for (s = 0; s < job->spin_dim; s++) {
  for (j = 0; j < nbands[s]; j++) {
  fprintf(file.out,"%3d %5d %3d %10.4lf\n",i, s, i * (nbands[0] + nbands[1]) + s * nbands[0] + j,eigenvalues[i * (nbands[0] + nbands[1]) + s \
  * nbands[0] + j] * au_to_eV); }
  //fprintf(file.out,"%3d %5d %10.4lf\n",s, i * nbands[s] + j, eigenvalues[s * fermi->nkunique * nbands[0] + i * nbands[s] + j] * au_to_eV); }
  fprintf(file.out,"\n"); } } }

}

void read_bse_eigenvalues_crystal(double *eigenvalues, int size, char *zz, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Read SCF, GW or CAS eigenvalues from disk                                              *
  // ******************************************************************************************

int i;
char buf1[110];
MPI_File eh;

  strcpy(buf1,file.scf_eigvec);
  strcat(buf1,zz);
  MPI_File_open(MPI_COMM_WORLD,buf1,MPI_MODE_RDONLY,MPI_INFO_NULL,&eh) ;
  MPI_File_seek(eh, 0, MPI_SEEK_SET) ;
  MPI_File_read(eh,eigenvalues,size,MPI_DOUBLE,MPI_STATUS_IGNORE);
  MPI_File_close(&eh);

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out,"reading from %s\n",buf1);
  fprintf(file.out,"BSE Eigenvalues\n");
  for (i = 0; i < size; i++)
  fprintf(file.out,"%5d %10.4lf\n",i,eigenvalues[i] * au_to_eV);
 }

}

void issue_MPI_IRecv(MPI_Request *request, Complex* Ham_buffer, int *MPA, int *NQA, int *ictxt, int *nbsize_row, int *nbsize_col, FERMI* fermi, CRYSTAL *crystal, FILES file, JOB_PARAM *job)

{

int i, k, l, kq, q1, q3, q4, s;
int k_bz, q_bz, kq_bz, fbz1;
int k_count, q_count;
int cblacs_taskid, itemp, nprow, npcol, myrow, mycol, mpA, nqA, izero = 0, info = 0, ione = 1;
int local_row, local_col, dest_row, dest_col, I1, I2;
int starts[2], subsize[2], largesize[2];
int spin_offset_row,spin_offset_col,spin_offset_msg,ntransitions[2];
double time1, time2;
VECTOR_INT kvec[2], qvec[2];

  ntransitions[0] = (fermi->bands[1] - fermi->homo[0]) * (fermi->homo[0] - fermi->bands[0] + 1);
  ntransitions[1] = 0;
  if (job->spin_dim == 2) 
  ntransitions[1] = (fermi->bands[3] - fermi->homo[1]) * (fermi->homo[1] - fermi->bands[2] + 1);

  for (s = 0; s < job->spin_dim; s++) {
    spin_offset_row = s * ntransitions[0];
    spin_offset_col = s * ntransitions[0];
    spin_offset_msg = s * fermi->nktot * fermi->nktot;
    time1 = MPI_Wtime();
    Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
    for (i = 0; i < job->numtasks; i++) {
      q_count = 0;
      for (q1 = 0; q1 < fermi->knet->unique; q1++) {
        for (q3 = 0; q3 < fermi->knet->num[q1]; q3++) {
          k_count = 0;
          for (k = 0; k < fermi->knet->unique; k++) {
            for (l = 0; l < fermi->knet->num[k]; l++) {
              kvec[0] = decompose_k_point(fermi->is, k_count, crystal, job, file);
              qvec[0] = decompose_k_point(fermi->is, q_count, crystal, job, file);
              kq  = compose_k_point(fermi->is, kvec[0], qvec[0], crystal, job, file);
              if (job->kss == 0) k_bz = k_count;                  // use for scf_evec_sym
              if (job->kss == 1) k_bz = fermi->knet->bz[k_count]; // use for scf_evec_no_sym
              q_bz = fermi->knet->bz[q_count];
              kvec[0]  = decompose_k_point(fermi->is, k_bz, crystal, job, file);
              qvec[0]  = decompose_k_point(fermi->is, q_bz, crystal, job, file);
              kq_bz = compose_k_point(fermi->is, kvec[0], qvec[0], crystal, job, file);
              I1 = spin_offset_row + k_bz  * (ntransitions[0] + ntransitions[1]);
              I2 = spin_offset_col + kq_bz * (ntransitions[0] + ntransitions[1]);
              dest_row = (I1 / *nbsize_row) % nprow;
              dest_col = (I2 / *nbsize_col) % npcol;
              cblacs_taskid = Cblacs_pnum(*ictxt,dest_row,dest_col);
              //printf("%3d %3d %3d %3d %3d\n",k,l,I1,I2,cblacs_taskid);
              local_row = *nbsize_row * (I1 / (*nbsize_row * nprow)) + I1 % *nbsize_row;
              local_col = *nbsize_col * (I2 / (*nbsize_col * npcol)) + I2 % *nbsize_col;
              //offset = local_row  + mpA * local_col;
              starts[0] = 2 * local_row;
              starts[1] = local_col;
              subsize[0] = 2 * ntransitions[s];
              subsize[1] = 1 * ntransitions[s];
              largesize[0] = 2 * MPA[cblacs_taskid];
              largesize[1] =     NQA[cblacs_taskid];
              //printf("ntrans %3d %3d %3d %3d %3d %3d\n",ntransitions[0],ntransitions[1],subsize[0],subsize[1],largesize[0],largesize[1]);
              MPI_Datatype sub_block;
              MPI_Type_create_subarray(2, largesize, subsize, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE, &sub_block);
              MPI_Type_commit(&sub_block);
              if (job->taskid == cblacs_taskid) { 
              MPI_Irecv(&Ham_buffer[0], 1, sub_block, i, spin_offset_msg + k_bz * fermi->nktot + kq_bz, MPI_COMM_WORLD, &request[0]);
              //printf("task %3d Irecv from %3d tag %3d\n",job->taskid,i,spin_offset_msg + k_bz * fermi->nktot + kq_bz);
             }
              MPI_Type_free(&sub_block);
              k_count++;
              //kq_count++;
             } // close loop on l
            } // close loop on k
           q_count++;
          } // close loop on q3
         } // close loop on q1
        } // close loop on i
       } // close loop on s
        time2 += MPI_Wtime() - time1;

}

/*
void density_fitting_crystal_contract_integrals_constrained(ComplexMatrix *solution1, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int i, j, k, q;
int a1, a2, i1, k1, l1, kp, lp, j0, j1, j2, j3, t2, q2;
int nd4, nd5, nd6b;
int count;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int dim11 = dim1 + 1, dim2 = dim1 * dim1;
int dim1a, dim456, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int bfpos6b, bfposk1, bfposl1;
int index_i, shelposi, gausposi, i4;
int index_j, shelposj, gausposj, j4;
int *ipiv;
int info = 0;
double normfac, expnta;
double mean_absolute_error;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10;
Complex *three_centre_integrals;
Complex orbital_product_norm_unconstrain[dim2], orbital_product_norm_constrained[dim2];
Complex maximum_error;
ComplexMatrix *V_mn_beta, *V_mn_beta1;
ComplexMatrix *solution;
TRIPLE_TRAN triple1;
PAIR_TRAN pair_p;
KPOINT_TRAN knet;
VECTOR_INT qvec;

  //count_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //if (job->kss == 0)      fermi->nkunique = knet.nktot;
  //else if (job->kss == 1) fermi->nkunique = knet.unique;
  //allocate_k_points(&knet,crystal,job,file);
  //generate_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //print_knet(&knet, fermi->is, crystal, job, file);
  //allocate_fermi(fermi,atoms,job,file);
  //fermi->knet = &knet;
  //fermi->nktot = knet.nktot;

  q = 0;
  //q = 7;
  //qvec.comp1 = 0;
  //qvec.comp2 = 0;
  //qvec.comp3 = 1;
  qvec.comp1 = fermi->knet->oblique[q].comp1;
  qvec.comp2 = fermi->knet->oblique[q].comp2;
  qvec.comp3 = fermi->knet->oblique[q].comp3;
  Q_LATTICE q_G;
  q_G.last_vector = G->last_vector;
  q_G.max_vector  = G->max_vector;
  allocate_Q_LATTICE(&q_G, job, file);
  generate_q_lattice(&qvec, &q_G, fermi, G, crystal, job, file);
  //free_Q_LATTICE(&q_G, job);

  // ******************************************************************************************
  // * Generate range selected pairs                                                          *
  // ******************************************************************************************

  pair_p.cutoff = 6.0;
pair_p.cutoff = 6.0;
  fprintf(file.out,"pair_cutoff = %10.4f\n",pair_p.cutoff);
  count_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  //generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  generate_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  print_pairs(&pair_p,atoms,R,job,file);

  // ******************************************************************************************
  // * Generate S_q for q = 0                                                                 *
  // ******************************************************************************************

  INT_1E one_ints, one_ints_buffer;
  int dim, dimg, nk[2];
  int Function[8];
  ComplexMatrix *S_q, *tmp;
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  array_dimensions(&dim, &dimg, &pair_p, atoms, job, file); // don't change
  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p,end_p,receive_p,offset_p,&pair_p,atoms,job,file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  AllocateComplexMatrix(&S_q,&dim1,&dim1,job);

  fock_element_1e1(&one_ints_buffer, dim, &pair_p, pair_p.nump, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  MPI_Allreduce(&one_ints_buffer.Overlap[0],&one_ints.Overlap[0],dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  nk[0] = 0;
  nk[1] = 0;
  printf("%3d %3d %3d \n",fermi->knet->oblique[q].comp1,fermi->knet->oblique[q].comp2, fermi->knet->oblique[q].comp3);

  fourier_transform(&one_ints.Overlap[0], &S_q->a[0][0], fermi->knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);

  // ******************************************************************************************
  // * Generate V_q and v_q for q = 0                                                         *
  // ******************************************************************************************

  ComplexMatrix *V_q, *v_q;
  AllocateComplexMatrix(&V_q,&dim1ax,&dim1ax,job);
  ResetComplexMatrix(V_q);
  generate_coulomb_matrix_inverse_complex(V_q,q,fermi,atom_p,atoms,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,\
  G,job,file);
  fprintf(file.out,"v_matrix inverse\n");
  print_complex_matrix(V_q, file);

  // ******************************************************************************************
  // * Allocate v_matrix for linear equations to solve for c_alpha expansion coefficients     *
  // ******************************************************************************************
 
  int dim11ax = dim1ax + 1;
  ComplexMatrix *v_matrix;

  AllocateComplexMatrix(&v_q,&dim1ax,&dim1ax,job);
  AllocateComplexMatrix(&v_matrix,&dim11ax,&dim11ax,job);
  ResetComplexMatrix(v_q);
  ResetComplexMatrix(v_matrix);

fprintf(file.out,"orbital prod\n");
  two_centre_coulomb1_crystal(v_q, R, &q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  print_complex_matrix(v_q, file);
  for (i = 0; i < dim1ax; i++) {
    for (j = 0; j < dim1ax; j++) {
      v_matrix->a[i][j] = v_q->a[i][j];
     }
    }

  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    shelposi = atoms_ax->shelposn_sh[j1];
    gausposi = atoms_ax->gausposn_sh[j1];
    bfpos6b = atoms_ax->bfnposn_sh[j1];
    for (index_i = shelposi; index_i < shelposi + atoms_ax->nshel_sh[j1]; index_i++) {
      if (shells_ax->type1[index_i] != 1) { //printf("bfpos6b %3d shells %3d\n",bfpos6b,shells_ax->type_sh[index_i]); 
      bfpos6b += shells_ax->type_sh[index_i]; gausposi += shells_ax->ng_sh[index_i]; continue; }
      normfac = k_zero;
      for (i4 = 0; i4 < shells_ax->ng_sh[index_i]; i4++) {
        expnta = gaussians_ax->expo_sh[gausposi + i4];
        normfac += gaussians_ax->sc[gausposi + i4] * pow(pi / expnta, 1.5);
        //normfac = pow(two * pi / expnta, 0.75);
       }
        v_matrix->a[bfpos6b][dim1ax] = normfac;
        v_matrix->a[dim1ax][bfpos6b] = normfac;   // need conjugate here
        gausposi += shells_ax->ng_sh[index_i];
        bfpos6b += shells_ax->type_sh[index_i];
       }
      }

  // ******************************************************************************************
  // * Generate V_q_mn for q = 0 and multiply in V_q inverse                                  *
  // ******************************************************************************************

  time2 = k_zero;
  time4 = k_zero;
  time6 = k_zero;
  time8 = k_zero;
  time10 = k_zero;

  int dimf, dimp;
  Complex *S1;
  ComplexMatrix *V_screen;
  time1 = MPI_Wtime();
  sh_array_dimensions(&dimp,&dimf,&pair_p,atoms,job,file); 
  AllocateComplexArray(&S1,&dimf,job);
  if (dimf != job->dimf) job->dimf = dimf;
  shell_screen_complex(S1,&pair_p,R,G,&q_G,atoms,shells,gaussians,symmetry,crystal,job,file);
  AllocateComplexMatrix(&V_screen,&dim1ax,&dim1ax,job);
  two_centre_coulomb1_crystal(V_screen, R, &q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  time2 = MPI_Wtime() - time1;

  // ******************************************************************************************
  // * Multiply inverse of coulomb matrix into three centre integrals                         *
  // ******************************************************************************************

  AllocateComplexMatrix(&V_mn_beta,&dim2,&dim1ax,job);
  ResetComplexMatrix(V_mn_beta);
  AllocateComplexMatrix(&V_mn_beta1,&dim2,&dim11ax,job);
  ResetComplexMatrix(V_mn_beta1);

    for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
      nd6b = atoms_ax->bfnnumb_sh[j1];
      allocate_TRIPLE_TRAN(&triple1, job, file);
      for (j = 0; j < pair_p.nump; j++) {
        q2 = pair_p.posn[j];
        nd4 = atoms->bfnnumb_sh[pair_p.cell1[q2]];
        nd5 = atoms->bfnnumb_sh[pair_p.cell2[q2]];
        dim456 = nd4 * nd5 * nd6b;
        AllocateComplexArray(&three_centre_integrals,&dim456,job);
        for (j2 = 0; j2 < pair_p.numb[j]; j2++) {
          kp = pair_p.cell1[q2 + j2];
          lp = pair_p.cell2[q2 + j2];
          t2 = pair_p.latt2[q2 + j2];
          //printf("j1  %3d  pair_p %3d %3d  %3d\n",j1,kp,lp,t2);
          j3 = q2 + j2;
          triple1.cell1[0] = kp;
          triple1.cell2[0] = lp;
          triple1.cell3[0] = j1;
          triple1.latt1[0] = 0;
          triple1.latt2[0] = t2;
          triple1.latt3[0] = 0;
          bfpos6b = atoms_ax->bfnposn_sh[j1];
          bfposk1 = atoms->bfnposn_sh[kp];
          bfposl1 = atoms->bfnposn_sh[lp];
 
          time5 = MPI_Wtime();
          int nshells = atoms->nshel_sh[kp] * atoms->nshel_sh[lp] * atoms_ax->nshel_sh[j1];
          int start_index[nshells];
          int skip = 1;
          for (i = 0; i < nshells; i++) start_index[i] = 0;
          shell_screen3(start_index,V_screen,S1,&pair_p,&triple1,atoms,shells,atoms_ax,shells_ax,job,file);
          //for (i = 0; i < nshells; i++) start_index[i] = 1;
          for (i = 0; i < nshells; i++) if (start_index[i] == 1) skip = 0;
          for (i = 0; i < nshells; i++) start_index[i] = 1;
          skip = 0;
          //printf("skip task %3d j1 %3d q1 %3d j, j2 %3d %3d skip %3d  %3d\n",job->taskid,j1,*q1,j,j2,skip,i);
          //fflush(stdout);
          //fprintf(file.out,"skip task %3d j1 %3d q1 %3d j, j2 %3d %3d skip %3d\n",job->taskid,j1,*q1,j,j2,skip);
          //fflush(file.out);
          //if (job->taskid == 0) printf("%3d  %3d out of %3d\n",j2,j+1,pair_p.nump);
          //for (i = 0; i < nshells; i++) fprintf(file.out,"start %3d %3d\n",i,start_index[i]);
          if (skip == 1) continue;
          ////if (triple.latt3[q2+j2] != 0) continue;
          //if (triple.cell1[q2+j2] != triple.cell2[q2+j2] || triple.latt3[q2+j2] != 0) continue;
          ResetComplexArray(three_centre_integrals,&dim456);
          three_centre_coulomb1_reversed2_crystal_test1(0,&triple1,start_index,three_centre_integrals,R,&q_G,G,atoms,shells, \
          gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,job,file);
          //for (i = 0; i < dim456; i++) { fprintf(file.out,"dim456 %3d triple1 %3d %3d %3d  %3d  %14.8f %14.8f\n",\
          dim456,kp,lp,j1,t2,(three_centre_integrals[i]).real(),(three_centre_integrals[i]).imag()); }

          count = 0;
          for (k1 = 0; k1 < nd4; k1++) {
            for (l1 = 0; l1 < nd5; l1++) {
              for (a2 = 0; a2 < nd6b; a2++) {
                V_mn_beta->a[ (bfposk1 + k1) * dim1 + bfposl1 + l1][bfpos6b + a2] += three_centre_integrals[count];
                V_mn_beta1->a[(bfposk1 + k1) * dim1 + bfposl1 + l1][bfpos6b + a2] += three_centre_integrals[count];
                V_mn_beta1->a[(bfposk1 + k1) * dim1 + bfposl1 + l1][dim1ax] = S_q->a[bfposk1 + k1][bfposl1 + l1];
                count++;
               }
              }
             }
            } // close loop on j2
           DestroyComplexArray(&three_centre_integrals,&dim456,job);
          } // close loop on j
         } // close loop on j1

  AllocateComplexMatrix(&solution,&dim2,&dim1ax,job);
  ResetComplexMatrix(solution);

  ResetComplexMatrix(solution1);
  for (i = 0; i < dim1ax; i++) {
    for (j = 0; j < dim1ax; j++) {
      for (k = 0; k < dim2; k++) {
        solution1->a[k][i] += V_q->a[i][j] * V_mn_beta->a[k][j];
       }
      }
     }

  //fprintf(file.out,"v_matrix\n");
  //print_complex_matrix(v_matrix, file);

  //fprintf(file.out,"V_mn_beta\n");
  //print_complex_matrix(V_mn_beta, file);

  //fprintf(file.out,"V_mn_beta1\n");
  //print_complex_matrix(V_mn_beta1, file);

  AllocateIntArray(&ipiv,&dim11ax,job);
  zgesv_(&dim11ax, &dim2, v_matrix->a[0], &dim11ax, ipiv, V_mn_beta1->a[0], &dim11ax, &info);

  //fprintf(file.out,"V_mn_beta solution\n");
  //print_complex_matrix(solution, file);

  fprintf(file.out,"V_mn_beta1 solution\n");
  print_complex_matrix(V_mn_beta1, file);

}

void wavefunction_product_density_fit_crystal_test1(ComplexMatrix *solution1, Q_LATTICE *q_G, FERMI *fermi, Complex **integral_buffer, int *nband_k, int *nband_kq, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k;
int a1, a2, i1, k1, l1, kp, lp, j0, j1, j2, j3, t2, q2;
int nd4, nd5, nd6b;
int count;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int dim11 = dim1 + 1, dim2 = dim1 * dim1;
int dim1a, dim456, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int bfpos6b, bfposk1, bfposl1;
int index_i, shelposi, gausposi, i4;
int index_j, shelposj, gausposj, j4;
int *ipiv;
int info = 0;
double normfac, expnta;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10;
ComplexMatrix *V_mn_beta, *V_mn_beta1;

  time2 = k_zero;
  time4 = k_zero;
  time6 = k_zero;
  time8 = k_zero;
  time10 = k_zero;

  // ******************************************************************************************
  // * Generate S_q for q = 0                                                                 *
  // ******************************************************************************************

  INT_1E one_ints, one_ints_buffer;
  int dim, dimg, nk[2];
  int Function[8];
  ComplexMatrix *S_q, *tmp;
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  array_dimensions(&dim, &dimg, pair_p, atoms, job, file); // don't change
  mpi_begin_end(begin_p,end_p,pair_p->nump,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p,end_p,receive_p,offset_p,pair_p,atoms,job,file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  AllocateComplexMatrix(&S_q,&dim1,&dim1,job);
  fock_element_1e1(&one_ints_buffer, dim, pair_p, pair_p->nump, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  MPI_Allreduce(&one_ints_buffer.Overlap[0],&one_ints.Overlap[0],dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  nk[0] = 0;
  nk[1] = 0;
  fourier_transform(&one_ints.Overlap[0], &S_q->a[0][0], fermi->knet, nk, pair_p, R, atoms, shells, symmetry, job, file);
  //printf("%3d %3d %3d \n",fermi->knet->oblique[q].comp1,fermi->knet->oblique[q].comp2, fermi->knet->oblique[q].comp3);

  // ******************************************************************************************
  // * Allocate v_matrix for linear equations to solve for c_alpha expansion coefficients     *
  // ******************************************************************************************
 
  int dim11ax = dim1ax + 1;
  ComplexMatrix *v_matrix, *v_q;

  AllocateComplexMatrix(&v_q,&dim1ax,&dim1ax,job);
  AllocateComplexMatrix(&v_matrix,&dim11ax,&dim11ax,job);
  ResetComplexMatrix(v_q);
  ResetComplexMatrix(v_matrix);
  two_centre_coulomb1_crystal(v_q, R, q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  for (i = 0; i < dim1ax; i++) {
    for (j = 0; j < dim1ax; j++) {
      v_matrix->a[i][j] = v_q->a[i][j];
     }
    }
  //print_complex_matrix(v_q, file);

  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    shelposi = atoms_ax->shelposn_sh[j1];
    gausposi = atoms_ax->gausposn_sh[j1];
    bfpos6b = atoms_ax->bfnposn_sh[j1];
    for (index_i = shelposi; index_i < shelposi + atoms_ax->nshel_sh[j1]; index_i++) {
      if (shells_ax->type1[index_i] != 1) { //printf("bfpos6b %3d shells %3d\n",bfpos6b,shells_ax->type_sh[index_i]); 
      bfpos6b += shells_ax->type_sh[index_i]; gausposi += shells_ax->ng_sh[index_i]; continue; }
      normfac = k_zero;
      for (i4 = 0; i4 < shells_ax->ng_sh[index_i]; i4++) {
        expnta = gaussians_ax->expo_sh[gausposi + i4];
        normfac += gaussians_ax->sc[gausposi + i4] * pow(pi / expnta, 1.5);
        //normfac = pow(two * pi / expnta, 0.75);
       }
        v_matrix->a[bfpos6b][dim1ax] = normfac;
        v_matrix->a[dim1ax][bfpos6b] = normfac;   // need conjugate here
        gausposi += shells_ax->ng_sh[index_i];
        bfpos6b += shells_ax->type_sh[index_i];
       }
      }

  // ******************************************************************************************
  // * Generate V_q_mn for q = 0 and multiply in V_q inverse                                  *
  // ******************************************************************************************

  AllocateComplexMatrix(&V_mn_beta,&dim2,&dim1ax,job);
  ResetComplexMatrix(V_mn_beta);
  AllocateComplexMatrix(&V_mn_beta1,&dim2,&dim11ax,job);
  ResetComplexMatrix(V_mn_beta1);

  for (k1 = 0; k1 < *nband_k; k1++) {
    for (l1 = 0; l1 < *nband_kq; l1++) {
      for (a1 = 0; a1 < dim1ax; a1++) {
        V_mn_beta1->a[k1 * *nband_kq +l1][a1] = integral_buffer[4][k1 * *nband_kq + l1];
        V_mn_beta1->a[(bfposk1 + k1) * dim1 + bfposl1 + l1][dim1ax] = S_q->a[bfposk1 + k1][bfposl1 + l1];
       }
        //V_mn_beta1->a[k1 * *nband_kq +l1][dim1ax] = S_q[k1 * *nband_kq + l1];
      }
     }

  AllocateIntArray(&ipiv,&dim11ax,job);
  zgesv_(&dim11ax, &dim2, v_matrix->a[0], &dim11ax, ipiv, V_mn_beta1->a[0], &dim11ax, &info);
  ResetComplexMatrix(solution1);
  for (i = 0; i < dim2; i++) {
    for (j = 0; j < dim1ax; j++) {
      solution1->a[i][j] = V_mn_beta1->a[i][j];
     }
    }

  fprintf(file.out,"Variational Fit\n");
  print_complex_matrix(solution1, file);
  fprintf(file.out,"V_mn_beta1\n");
  print_complex_matrix(V_mn_beta1, file);
  DestroyComplexMatrix(&V_mn_beta,job);
  DestroyComplexMatrix(&V_mn_beta1,job);
  DestroyComplexMatrix(&v_q,job);
  DestroyComplexMatrix(&v_matrix,job);
  DestroyComplexMatrix(&S_q,job);
  free_INT_1E(&one_ints, Function, job, file);
  free_INT_1E(&one_ints_buffer, Function, job, file);

}

void orbital_product_density_fit_crystal_test1(ComplexMatrix *solution1, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int i, j, k, q;
int a1, a2, i1, k1, l1, kp, lp, j0, j1, j2, j3, t2, q2;
int nd4, nd5, nd6b;
int count;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int dim11 = dim1 + 1, dim2 = dim1 * dim1;
int dim1a, dim456, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int bfpos6b, bfposk1, bfposl1;
int index_i, shelposi, gausposi, i4;
int index_j, shelposj, gausposj, j4;
int *ipiv;
int info = 0;
double normfac, expnta;
double mean_absolute_error;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10;
Complex *three_centre_integrals;
Complex orbital_product_norm_unconstrain[dim2], orbital_product_norm_constrained[dim2];
Complex maximum_error;
ComplexMatrix *V_mn_beta, *V_mn_beta1;
ComplexMatrix *solution;
TRIPLE_TRAN triple1;
PAIR_TRAN pair_p;
KPOINT_TRAN knet;
VECTOR_INT qvec;

  //count_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //if (job->kss == 0)      fermi->nkunique = knet.nktot;
  //else if (job->kss == 1) fermi->nkunique = knet.unique;
  //allocate_k_points(&knet,crystal,job,file);
  //generate_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //print_knet(&knet, fermi->is, crystal, job, file);
  //allocate_fermi(fermi,atoms,job,file);
  //fermi->knet = &knet;
  //fermi->nktot = knet.nktot;

  q = 0;
  //q = 7;
  //qvec.comp1 = 0;
  //qvec.comp2 = 0;
  //qvec.comp3 = 1;
  qvec.comp1 = fermi->knet->oblique[q].comp1;
  qvec.comp2 = fermi->knet->oblique[q].comp2;
  qvec.comp3 = fermi->knet->oblique[q].comp3;
  Q_LATTICE q_G;
  q_G.last_vector = G->last_vector;
  q_G.max_vector  = G->max_vector;
  allocate_Q_LATTICE(&q_G, job, file);
  generate_q_lattice(&qvec, &q_G, fermi, G, crystal, job, file);
  //free_Q_LATTICE(&q_G, job);

  // ******************************************************************************************
  // * Generate range selected pairs                                                          *
  // ******************************************************************************************

  pair_p.cutoff = 6.0;
pair_p.cutoff = 2.0;
  fprintf(file.out,"pair_cutoff = %10.4f\n",pair_p.cutoff);
  count_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  //generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  generate_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //print_pairs(&pair_p,atoms,R,job,file);

  // ******************************************************************************************
  // * Generate S_q for q = 0                                                                 *
  // ******************************************************************************************

  INT_1E one_ints, one_ints_buffer;
  int dim, dimg, nk[2];
  int Function[8];
  ComplexMatrix *S_q, *tmp;
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  array_dimensions(&dim, &dimg, &pair_p, atoms, job, file); // don't change
  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p,end_p,receive_p,offset_p,&pair_p,atoms,job,file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  AllocateComplexMatrix(&S_q,&dim1,&dim1,job);

  fock_element_1e1(&one_ints_buffer, dim, &pair_p, pair_p.nump, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  MPI_Allreduce(&one_ints_buffer.Overlap[0],&one_ints.Overlap[0],dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  nk[0] = 0;
  nk[1] = 0;
  //printf("%3d %3d %3d \n",fermi->knet->oblique[q].comp1,fermi->knet->oblique[q].comp2, fermi->knet->oblique[q].comp3);

  fourier_transform(&one_ints.Overlap[0], &S_q->a[0][0], fermi->knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);

  // ******************************************************************************************
  // * Generate V_q and v_q for q = 0                                                         *
  // ******************************************************************************************

  ComplexMatrix *V_q, *v_q;
  AllocateComplexMatrix(&V_q,&dim1ax,&dim1ax,job);
  ResetComplexMatrix(V_q);
  generate_coulomb_matrix_inverse_complex(V_q,q,fermi,atom_p,atoms,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,\
  G,job,file);
  //fprintf(file.out,"v_matrix inverse\n");
  //print_complex_matrix(V_q, file);

  // ******************************************************************************************
  // * Allocate v_matrix for linear equations to solve for c_alpha expansion coefficients     *
  // ******************************************************************************************
 
  int dim11ax = dim1ax + 1;
  ComplexMatrix *v_matrix;

  AllocateComplexMatrix(&v_q,&dim1ax,&dim1ax,job);
  AllocateComplexMatrix(&v_matrix,&dim11ax,&dim11ax,job);
  ResetComplexMatrix(v_q);
  ResetComplexMatrix(v_matrix);

//fprintf(file.out,"orbital prod\n");
  two_centre_coulomb1_crystal(v_q, R, &q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  //print_complex_matrix(v_q, file);
  for (i = 0; i < dim1ax; i++) {
    for (j = 0; j < dim1ax; j++) {
      v_matrix->a[i][j] = v_q->a[i][j];
     }
    }

  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    shelposi = atoms_ax->shelposn_sh[j1];
    gausposi = atoms_ax->gausposn_sh[j1];
    bfpos6b = atoms_ax->bfnposn_sh[j1];
    for (index_i = shelposi; index_i < shelposi + atoms_ax->nshel_sh[j1]; index_i++) {
      if (shells_ax->type1[index_i] != 1) { //printf("bfpos6b %3d shells %3d\n",bfpos6b,shells_ax->type_sh[index_i]); 
      bfpos6b += shells_ax->type_sh[index_i]; gausposi += shells_ax->ng_sh[index_i]; continue; }
      normfac = k_zero;
      for (i4 = 0; i4 < shells_ax->ng_sh[index_i]; i4++) {
        expnta = gaussians_ax->expo_sh[gausposi + i4];
        normfac += gaussians_ax->sc[gausposi + i4] * pow(pi / expnta, 1.5);
        //normfac = pow(two * pi / expnta, 0.75);
       }
        v_matrix->a[bfpos6b][dim1ax] = normfac;
        v_matrix->a[dim1ax][bfpos6b] = normfac;   // need conjugate here
        gausposi += shells_ax->ng_sh[index_i];
        bfpos6b += shells_ax->type_sh[index_i];
       }
      }

  // ******************************************************************************************
  // * Generate V_q_mn for q = 0 and multiply in V_q inverse                                  *
  // ******************************************************************************************

  time2 = k_zero;
  time4 = k_zero;
  time6 = k_zero;
  time8 = k_zero;
  time10 = k_zero;

  int dimf, dimp;
  Complex *S1;
  ComplexMatrix *V_screen;
  time1 = MPI_Wtime();
  sh_array_dimensions(&dimp,&dimf,&pair_p,atoms,job,file); 
  AllocateComplexArray(&S1,&dimf,job);
  if (dimf != job->dimf) job->dimf = dimf;
  shell_screen_complex(S1,&pair_p,R,G,&q_G,atoms,shells,gaussians,symmetry,crystal,job,file);
  AllocateComplexMatrix(&V_screen,&dim1ax,&dim1ax,job);
  two_centre_coulomb1_crystal(V_screen, R, &q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  time2 = MPI_Wtime() - time1;

  // ******************************************************************************************
  // * Multiply inverse of coulomb matrix into three centre integrals                         *
  // ******************************************************************************************

  AllocateComplexMatrix(&V_mn_beta,&dim2,&dim1ax,job);
  ResetComplexMatrix(V_mn_beta);
  AllocateComplexMatrix(&V_mn_beta1,&dim2,&dim11ax,job);
  ResetComplexMatrix(V_mn_beta1);

    for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
      nd6b = atoms_ax->bfnnumb_sh[j1];
      allocate_TRIPLE_TRAN(&triple1, job, file);
      for (j = 0; j < pair_p.nump; j++) {
        q2 = pair_p.posn[j];
        nd4 = atoms->bfnnumb_sh[pair_p.cell1[q2]];
        nd5 = atoms->bfnnumb_sh[pair_p.cell2[q2]];
        dim456 = nd4 * nd5 * nd6b;
        AllocateComplexArray(&three_centre_integrals,&dim456,job);
        for (j2 = 0; j2 < pair_p.numb[j]; j2++) {
          kp = pair_p.cell1[q2 + j2];
          lp = pair_p.cell2[q2 + j2];
          t2 = pair_p.latt2[q2 + j2];
          //printf("j1  %3d  pair_p %3d %3d  %3d\n",j1,kp,lp,t2);
          j3 = q2 + j2;
          triple1.cell1[0] = kp;
          triple1.cell2[0] = lp;
          triple1.cell3[0] = j1;
          triple1.latt1[0] = 0;
          triple1.latt2[0] = t2;
          triple1.latt3[0] = 0;
          bfpos6b = atoms_ax->bfnposn_sh[j1];
          bfposk1 = atoms->bfnposn_sh[kp];
          bfposl1 = atoms->bfnposn_sh[lp];
 
          time5 = MPI_Wtime();
          int nshells = atoms->nshel_sh[kp] * atoms->nshel_sh[lp] * atoms_ax->nshel_sh[j1];
          int start_index[nshells];
          int skip = 1;
          for (i = 0; i < nshells; i++) start_index[i] = 0;
          shell_screen3(start_index,V_screen,S1,&pair_p,&triple1,atoms,shells,atoms_ax,shells_ax,job,file);
          //for (i = 0; i < nshells; i++) start_index[i] = 1;
          for (i = 0; i < nshells; i++) if (start_index[i] == 1) skip = 0;
          for (i = 0; i < nshells; i++) start_index[i] = 1;
          skip = 0;
          //printf("skip task %3d j1 %3d q1 %3d j, j2 %3d %3d skip %3d  %3d\n",job->taskid,j1,*q1,j,j2,skip,i);
          //fflush(stdout);
          //fprintf(file.out,"skip task %3d j1 %3d q1 %3d j, j2 %3d %3d skip %3d\n",job->taskid,j1,*q1,j,j2,skip);
          //fflush(file.out);
          //if (job->taskid == 0) printf("%3d  %3d out of %3d\n",j2,j+1,pair_p.nump);
          //for (i = 0; i < nshells; i++) fprintf(file.out,"start %3d %3d\n",i,start_index[i]);
          if (skip == 1) continue;
          ////if (triple.latt3[q2+j2] != 0) continue;
          //if (triple.cell1[q2+j2] != triple.cell2[q2+j2] || triple.latt3[q2+j2] != 0) continue;
          ResetComplexArray(three_centre_integrals,&dim456);
          three_centre_coulomb1_reversed2_crystal_test1(0,&triple1,start_index,three_centre_integrals,R,&q_G,G,atoms,shells, \
          gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,job,file);
          //for (i = 0; i < dim456; i++) { fprintf(file.out,"dim456 %3d triple1 %3d %3d %3d  %3d  %14.8f %14.8f\n",\
          dim456,kp,lp,j1,t2,(three_centre_integrals[i]).real(),(three_centre_integrals[i]).imag()); }

          count = 0;
          for (k1 = 0; k1 < nd4; k1++) {
            for (l1 = 0; l1 < nd5; l1++) {
              for (a2 = 0; a2 < nd6b; a2++) {
                V_mn_beta->a[ (bfposk1 + k1) * dim1 + bfposl1 + l1][bfpos6b + a2] += three_centre_integrals[count];
                V_mn_beta1->a[(bfposk1 + k1) * dim1 + bfposl1 + l1][bfpos6b + a2] += three_centre_integrals[count];
                V_mn_beta1->a[(bfposk1 + k1) * dim1 + bfposl1 + l1][dim1ax] = S_q->a[bfposk1 + k1][bfposl1 + l1];
                count++;
               }
              }
             }
            } // close loop on j2
           DestroyComplexArray(&three_centre_integrals,&dim456,job);
          } // close loop on j
         } // close loop on j1

  AllocateComplexMatrix(&solution,&dim2,&dim1ax,job);
  ResetComplexMatrix(solution);

  ResetComplexMatrix(solution1);
  for (i = 0; i < dim1ax; i++) {
    for (j = 0; j < dim1ax; j++) {
      for (k = 0; k < dim2; k++) {
        solution1->a[k][i] += V_q->a[i][j] * V_mn_beta->a[k][j];
       }
      }
     }

  fprintf(file.out,"Robust Fit\n");
  print_complex_matrix(solution1, file);

  //fprintf(file.out,"v_matrix\n");
  //print_complex_matrix(v_matrix, file);

  //fprintf(file.out,"V_mn_beta\n");
  //print_complex_matrix(V_mn_beta, file);

  //fprintf(file.out,"V_mn_beta1\n");
  //print_complex_matrix(V_mn_beta1, file);

  AllocateIntArray(&ipiv,&dim11ax,job);
  zgesv_(&dim11ax, &dim2, v_matrix->a[0], &dim11ax, ipiv, V_mn_beta1->a[0], &dim11ax, &info);

  ResetComplexMatrix(solution1);
  for (i = 0; i < dim2; i++) {
    for (j = 0; j < dim1ax; j++) {
      solution1->a[i][j] = V_mn_beta1->a[i][j];
     }
    }

  fprintf(file.out,"Variational Fit\n");
  print_complex_matrix(V_mn_beta1, file);

 // int s_basis_functions = 0;
 // mean_absolute_error = k_zero;
 // maximum_error = Complex(k_zero, k_zero);
 // for (k = 0; k < dim2; k++) {
 //   orbital_product_norm_unconstrain[k] = Complex(k_zero, k_zero);
 //   orbital_product_norm_constrained[k] = Complex(k_zero, k_zero);
 //  }
 // count = 0;
 // for (j0 = 0; j0 < atoms_ax->number_of_atoms_in_unit_cell; j0++) {
 // shelposi = atoms_ax->shelposn_sh[j0];
 // gausposi = atoms_ax->gausposn_sh[j0];
 // bfpos6b = atoms_ax->bfnposn_sh[j0];
 // for (index_i = shelposi; index_i < shelposi + atoms_ax->nshel_sh[j0]; index_i++) {
 //   if (shells_ax->type1[index_i] != 1) { //printf("bfpos6b %3d shells %3d\n",bfpos6b,shells_ax->type_sh[index_i]); 
 //   bfpos6b += shells_ax->type_sh[index_i]; gausposi += shells_ax->ng_sh[index_i]; continue; }
 //   normfac = k_zero;
 //   for (i4 = 0; i4 < shells_ax->ng_sh[index_i]; i4++) {
 //     expnta = gaussians_ax->expo_sh[gausposi + i4];
 //     normfac += gaussians_ax->sc[gausposi + i4] * pow(pi / expnta, 1.5);
 //     //printf("%3d %10.4f %10.4f %10.4f %10.4f\n",\
 //     i4,expnta,gaussians_ax->sc[gausposi + i4], gaussians_ax->sc[gausposi + i4] * pow(pi / expnta, 1.5),normfac);
 //     //normfac = pow(two * pi / expnta, 0.75);
 //    }
 //     s_basis_functions++;
 //     //printf("j0 %3d bfpos6b %3d %3d %3d %f %f shells %3d %3d\n",\
 //     j0,bfpos6b,shelposi,gausposi,expnta,normfac,shells_ax->type_sh[index_i],dim2); 
 //     for (k = 0; k < dim2; k++) {
 //       orbital_product_norm_unconstrain[k] += solution->a[k][bfpos6b] * normfac;
 //       orbital_product_norm_constrained[k] += V_mn_beta1->a[k][bfpos6b] * normfac;
 //      }
 //     gausposi += shells_ax->ng_sh[index_i];
 //     bfpos6b += shells_ax->type_sh[index_i];
 //    } // close loop on index_i
 //   } // close loop on j0
 //
 // count = 0;
 // int max_i, max_j;
 // int bin_ij[20][dim1][dim1];
 // int bin_off[20];
 // int rel_error_diag[20];
 // int rel_error_off_diag[20];
 // for (i = 0; i <20 ; i++) { for (j = 0; j < dim1; j++) { for (k = 0; k < dim1; k++) { bin_ij[i][j][k] = 0; }}}
 // for (i = 0; i < 20; i++) { bin_off[i] = 0; rel_error_diag[i] = 0; rel_error_diag[i] = 0; }
 // 
 // for (i = 0; i < dim1; i++) {
 //   for (j = 0; j < dim1; j++) {
 //     if (fabs((S_q->a[i][j]).real()) < 1e-20) (S_q->a[i][j]).real() = k_zero;
 //     fprintf(file.out,"S_q %3d %3d %14.8f %14.8f  %14.8f     %14.8f  %9.2e  %10.4f %3d   %14.8f  %3d\n",\
 //     i,j,(S_q->a[i][j]).real(), (orbital_product_norm_constrained[count]).real(),\
 //     (orbital_product_norm_unconstrain[count]).real(),\
 //     (orbital_product_norm_constrained[count]).real() - (S_q->a[i][j]).real(),\
 //     (orbital_product_norm_unconstrain[count]).real() - (S_q->a[i][j]).real(),\
 //     //(orbital_product_norm_constrained[count]).real() / (S_q->a[i][j]).real(),
 //     floor(log10(fabs((orbital_product_norm_unconstrain[count]).real() - (S_q->a[i][j]).real() - 1e-60))),\
 //     (int)-floor(log10(fabs((orbital_product_norm_unconstrain[count]).real() - (S_q->a[i][j]).real() - 1e-60))),\
 //     (orbital_product_norm_unconstrain[count] - S_q->a[i][j] - 1e-60).real() / (S_q->a[i][j] + 1e-20).real(), \
 //     (int)-floor(log10(fabs((orbital_product_norm_unconstrain[count] - S_q->a[i][j] - 1e-60).real() / (S_q->a[i][j] + 1e-20).real()))));
 //      //fprintf(file.out,"%3d %3d %10.4f\n",i,j,log10(fabs((orbital_product_norm_unconstrain[count]).real() - (S_q->a[i][j]).real())));
 //      if ((int)-floor(log10(fabs((orbital_product_norm_unconstrain[count]).real() - (S_q->a[i][j]).real()))) < 20 && \
 //      (int)-floor(log10(fabs((orbital_product_norm_unconstrain[count]).real() - (S_q->a[i][j]).real()))) >= 0) { \
 //      //if (i == j)
 //      bin_ij[(int)-floor(log10(fabs((orbital_product_norm_unconstrain[count]).real() - (S_q->a[i][j]).real())))][i][j]++; \
 //      if (i != j)
 //      bin_off[(int)-floor(log10(fabs((orbital_product_norm_unconstrain[count]).real() - (S_q->a[i][j]).real())))]++; \
 //     }
 //      if (fabs((orbital_product_norm_unconstrain[count] - S_q->a[i][j]).real()) > fabs((maximum_error).real())) {
 //      maximum_error = orbital_product_norm_unconstrain[count] - S_q->a[i][j];
 //      mean_absolute_error += fabs((orbital_product_norm_unconstrain[count] - S_q->a[i][j]).real());
 //      max_i = i;
 //      max_j = j;
 //     }
 //      count++;
 //    }
 //   }
 //
 // int count_i, count_j, shelmax = 5;
 // int shell_bin[shelmax][shelmax][20];
 // for (i = 0; i < shelmax; i++) { for (j = 0; j < shelmax; j++) { for (k = 0; k < 20; k++) { shell_bin[i][j][k] = 0; }}}
 //
 // shelposi = atoms->shelposn_sh[0];
 // shelposj = atoms->shelposn_sh[0];
 // count_i = 0;
 // for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[0]; index_i++) {
 //   for (i1 = 0; i1 < shells->type_sh[index_i]; i1++) {
 //     count_j = 0;
 //     for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[0]; index_j++) {
 //       if (shells->type_sh[index_i] == 1) i = 0;
 //       if (shells->type_sh[index_i] == 3) i = 1;
 //       if (shells->type_sh[index_i] == 5) i = 2;
 //       if (shells->type_sh[index_i] == 7) i = 3;
 //       if (shells->type_sh[index_i] == 9) i = 4;
 //       if (shells->type_sh[index_j] == 1) j = 0;
 //       if (shells->type_sh[index_j] == 3) j = 1;
 //       if (shells->type_sh[index_j] == 5) j = 2;
 //       if (shells->type_sh[index_j] == 7) j = 3;
 //       if (shells->type_sh[index_j] == 9) j = 4;
 //       for (j1 = 0; j1 < shells->type_sh[index_j]; j1++) {
 //         if (count_i <= count_j)
 //         for (k = 0; k < 20; k++) {
 //           shell_bin[i][j][k] += bin_ij[k][count_i][count_j];
 //          }
 //         count_j++;
 //        }
 //       //printf("%3d %3d  %3d %3d\n",index_i,index_j,shells->type_sh[index_i],shells->type_sh[index_j]);
 //       }
 //      count_i++;
 //     }
 //   }
 //
 // for (k = 0; k < 20; k++) { 
 // fprintf(file.out,\
 // "mag. %3d s|s %3d s|p %3d s|d %3d s|f %3d s|g %3d p|p %3d p|d %3d p|f %3d p|g %3d d|d %3d d|f %3d d|g %3d f|f %3d f|g %3d g|g %3d\n",\
 // k,shell_bin[0][0][k],shell_bin[0][1][k],shell_bin[0][2][k],shell_bin[0][3][k],shell_bin[0][4][k],shell_bin[1][1][k],shell_bin[1][2][k],\
 //   shell_bin[1][3][k],shell_bin[1][4][k],shell_bin[2][2][k],shell_bin[2][3][k],shell_bin[2][4][k],shell_bin[3][3][k],shell_bin[3][4][k],\
 //   shell_bin[4][4][k]);
 //}
 // for (k = 0; k < 20; k++) { fprintf(file.out,"%3d  %3d %3d %3d %3d %3d   %3d %3d %3d %3d   %3d %3d %3d   %3d %3d   %3d\n",\
 // k,shell_bin[0][0][k],shell_bin[0][1][k],shell_bin[0][2][k],shell_bin[0][3][k],shell_bin[0][4][k],shell_bin[1][1][k],shell_bin[1][2][k],\
 //   shell_bin[1][3][k],shell_bin[1][4][k],shell_bin[2][2][k],shell_bin[2][3][k],shell_bin[2][4][k],shell_bin[3][3][k],shell_bin[3][4][k],\
 //   shell_bin[4][4][k]);
 //}
 //   
 // //for (i = 0; i < 20; i++) { fprintf(file.out,"%3d %3d %3d\n",i,bin_ij[i][0][0],bin_off[i]); }
 //
 // fprintf(file.out,"max error %14.8f max_i,j %3d %3d mean absolute error %14.8f %f\n",\
 // (maximum_error).real(), mean_absolute_error, (double) s_basis_functions);
 //
 // free_TRIPLE_TRAN(&triple1,job);
 // DestroyComplexMatrix(&solution,job);
 // DestroyComplexMatrix(&V_mn_beta,job);
 // DestroyComplexMatrix(&V_mn_beta1,job);
 // DestroyComplexArray(&S1,&dimf,job);
 // DestroyComplexMatrix(&V_screen,job);
 // DestroyComplexMatrix(&v_q,job);
 // DestroyComplexMatrix(&v_matrix,job);
 // DestroyComplexMatrix(&V_q,job);
 // DestroyComplexMatrix(&S_q,job);
 // free_INT_1E(&one_ints, Function, job, file);
 // free_INT_1E(&one_ints_buffer, Function, job, file);
 // free_PAIR_TRAN(&pair_p,job);
 // free_Q_LATTICE(&q_G,job);
  //time2 = MPI_Wtime() - time1;
  if (job->taskid == 0) printf("orbital_product_density_fit_crystal_test   %10.2f %10.2f %10.2f %10.2f %10.2f\n",\
  time2,time4,time6,time8,time10);
//MPI_Finalize();
//exit(0);

}
*/
