
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
#include "mycomplex.h"
#include "conversion_factors.h"
#include "myconstants.h"
#include "LIMITS.h"
#include "USER_DATA.h"
#include "SETUP_SYMMETRY.h"
#include "PAIRS_QUADS.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "PARALLEL.h"
#include "DENSITY_MATRIX_MOLECULE.h"
#include "SYMMETRY_ADAPTATION.h"
#include "MATRIX_UTIL.h"
#include "PRINT_MOLECULE.h"
#include "BUILD_FOCK_MATRIX_MOLECULE.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "ROTATIONS_MOLECULE.h"
#include "FOURIER_TRANSFORM.h"
#include "INTEGRALS_2C_MOLECULE.h"
#include "SCF_MOLECULE.h"

using namespace std;

  // ******************************************************************************************
  // * This routine is the SCF loop                                                           *
  // ******************************************************************************************

void scf_molecule(FERMI *fermi, ATOM *atoms, ATOM *atoms_ax, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  int begin_k, end_k;
  int i, j, k, l, n, p, q, s, ii, iii, jjj;
  int dim1, dim2, dim3;
  int dimp_spin, dimf_spin, dimp_read, dimf_read, dimp_read_spin, dimf_read_spin, spin_offset;
  int dim, dimf, dimg, dimp, dim11;
  int count;
  int nk[2], ksize;
  int num_sym = symmetry->number_of_operators;
  int vector_size, value_size;
  int mixing_order = job->mixing_order;
  int dim_salc2[symmetry->number_of_classes];
  int ngrid_points;
  int restart_parameters[6];
  int Function[8];
  int info = 0, izero = 0;
  int *irrep, *p_irrep;
  double total_electrons, total_population, fac;
  double *eigenvalues;
  double time1, time2, time3, time4, time5, time6, time7;
  double homo, lumo, gap;
  double *Fock, *eigval, *eigval1, *p_eigval;
  double *P0, *P1, *F0, *F1;
  double *P2;
  double *delta_P0, *delta_F0, *lse_F0;
  double threshold = 1e-04;
  char NoTrans = 'N', ConjTrans = 'C';
  char uplo = 'U';
  char jobz = 'V';
  char buf1[110], buf2[110], buf4[110];
  char xx[10] = "/evalfile", yy[10] = "/scf_evec", zz[20] = "/new_density_matrix";
  char zz2[24] = "scf_evectors";
  FILE *scf_evectors; 
  Complex alpha = Complex(k_one, k_zero);
  Complex beta = Complex(k_zero, k_zero);
  PAIR_TRAN pair_p;
  KPOINT_TRAN knet;
  INT_1E one_ints, one_ints_buffer, one_ints_buffer1;
  DoubleMatrix *C_i, *C_o, *D;
  ComplexMatrix *Residual[job->spin_dim * job->mixing_order * symmetry->number_of_classes];
  ComplexMatrix *FockHist[job->spin_dim * job->mixing_order * symmetry->number_of_classes];
  ComplexMatrix *eigenvectors, *eigvec;
  MPI_File eh;
  MPI_File fh;
  MPI_File gh;

  strcpy(buf1,file.scf_eigvec);
  strcat(buf1,xx);
  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,yy);
  strcpy(buf4,file.scf_eigvec);
  strcat(buf4,zz);

  // ******************************************************************************************
  // * Count and generate atom pairs for SCF calculation                                      *
  // * Count size of reduced (dimp) and full (dimf) atom pair array                           *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  if (job->taskid == 0 && job->verbosity > 1) 
  print_pairs(&pair_p,atoms,R,job,file);
  sh_array_dimensions(&dimp,&dimf,&pair_p,atoms,job,file); 

  // ******************************************************************************************
  // * Allocate and initialize P and F and matrices for a new SCF calculation                 *
  // ******************************************************************************************

  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;
  AllocateDoubleArray(&Fock,&dimp_spin,job);
  dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  dim2 = dim1 * dim1;
  AllocateDoubleArray(&eigenvalues,&dim1,job);
  ResetDoubleArray(eigenvalues,&dim1);
  AllocateComplexMatrix(&eigenvectors,&dim1,&dim1,job);
  AllocateComplexMatrix(&eigvec,&dim1,&dim1,job);

  if (job->guess_type == 0) {
  AllocateDoubleArray(&P0,&dimp_spin,job);
  AllocateDoubleArray(&F0,&dimf_spin,job);
  AllocateDoubleArray(&delta_P0,&dimp_spin,job);
  AllocateDoubleArray(&delta_F0,&dimf_spin,job);
  ResetDoubleArray(P0,&dimp_spin);
  ResetDoubleArray(F0,&dimf_spin);
  ResetDoubleArray(delta_P0,&dimp_spin);
  ResetDoubleArray(delta_F0,&dimf_spin);
 }
  AllocateDoubleArray(&P1,&dimp_spin,job);
  AllocateDoubleArray(&F1,&dimf_spin,job);
  AllocateDoubleArray(&P2,&dimp_spin,job);

  // ******************************************************************************************
  // * Initial wave function guess                                                            *
  // * Allocate reduced (P) and full (F) density matrices for a restart                       *
  // ******************************************************************************************

  if (job->guess_type == 1) {
    dimp_read = 0;
    dimf_read = 0;
    if (job->taskid == 0) printf("Restart run\n");
    read_density_matrix(fermi,&P0,&dimp_read,&dimf_read,atoms,job,file);
    dimp_read_spin = dimp_read * job->spin_dim;
    dimf_read_spin = dimf_read * job->spin_dim;
    if (dimp != dimp_read || dimf != dimf_read) {
    if (job->taskid == 0 || job->verbosity > 1)
      fprintf(file.out,"ERROR: Reduced/Full density matrices read from disk have sizes %d/%d.\n\
			       Reduced/Full density matrices sizes calculated to be %d/%d\n",dimp_read,dimf_read,dimp,dimf);
      MPI_Finalize();
      exit(1);
    }
    dimp = dimp_read;
    dimf = dimf_read;
    dimp_spin = dimp * job->spin_dim;
    dimf_spin = dimf * job->spin_dim;
    if (job->taskid == 0 || job->verbosity > 1)
    fprintf(file.out,"Reduced density matrix read from disk.\nReduced/Full density matrices have sizes %d/%d.\n\n",dimp,dimf);
    ksize = 1;
    fermi->nktot = ksize;
    allocate_fermi(fermi,atoms,job,file);
    AllocateDoubleArray(&P0,&dimp_read_spin,job);
    AllocateDoubleArray(&F0,&dimf_read_spin,job);
    AllocateDoubleArray(&delta_P0,&dimp_read_spin,job);
    AllocateDoubleArray(&delta_F0,&dimf_read_spin,job);
    ResetDoubleArray(delta_P0,&dimp_read_spin);
    ResetDoubleArray(delta_F0,&dimf_read_spin);
  }

  // ******************************************************************************************
  // * Generate k points for SCF calculation (kpoints == 0) or BAND STRUCTURE (kpoints == 1)  *
  // ******************************************************************************************

    if (job->kpoints == 0) {
    knet.nktot = 1;
    knet.unique = 1;
    fermi->knet = &knet;
    fermi->nkunique = knet.nktot;
    fermi->nktot = knet.nktot;
    allocate_k_points(&knet,crystal,job,file);
    knet.ibz[0] = 0;
    knet.opr[0] = 0;
    knet.num[0] = 1;
    allocate_fermi(fermi,atoms,job,file);
   }

  // ******************************************************************************************
  // * Allocate and initialize P and F and matrices guess_type == 0 new guess_type == 1 old   *
  // ******************************************************************************************

    if (job->guess_type == 0) {
    initial_density_matrix(P0,F0,&pair_p,fermi,atoms,atom_p,shells,gaussians,crystal,symmetry,R,R_tables,G,job,file);
    ResetDoubleArray(F0,&dimf_spin);
    expand_density_matrix(P0,F0,&pair_p,atoms,shells,symmetry,job,file);
   }
    if (job->guess_type == 1) {
    read_density_matrix(fermi,&P0,&dimp_read_spin,&dimf_read_spin,atoms,job,file);
    expand_density_matrix(P0,F0,&pair_p,atoms,shells,symmetry,job,file);
   }

    dim3 = job->spin_dim * knet.unique * atoms->number_of_sh_bfns_in_unit_cell;
    AllocateIntArray(&irrep,&dim3,job);
    AllocateDoubleArray(&eigval,&dim3,job);
    AllocateDoubleArray(&eigval1,&dim3,job);
    int job_parameters[41 + job->spin_dim * fermi->nkunique];
    double job_parameters1[11 + dim3];
    write_JOB_PARAM_array(fermi, job_parameters, job_parameters1, job, file);

    vector_size = (fermi->bands[1] - fermi->bands[0] + 1) * dim1;
    value_size  = (fermi->bands[1] - fermi->bands[0] + 1);
    if (fermi->bands[1] > atoms->number_of_sh_bfns_in_unit_cell || fermi->bands[1] < 1) {
      if (job->taskid == 0)
      fprintf(file.out,"Upper limit in band range %5d must not exceed number of AOs %5d or be less than 1\n", \
      fermi->bands[1],atoms->number_of_sh_bfns_in_unit_cell);
      MPI_Finalize();
      exit(1);
     }
    if (fermi->bands[0] > atoms->number_of_sh_bfns_in_unit_cell || fermi->bands[0] < 1) {
      if (job->taskid == 0)
      fprintf(file.out,"Lower limit in band range %5d must not exceed number of AOs %5d or be less than 1\n", \
      fermi->bands[0],atoms->number_of_sh_bfns_in_unit_cell);
      MPI_Finalize();
      exit(1);
     }

  // ******************************************************************************************
  // * Nuclear repulsion energy                                                               *
  // ******************************************************************************************

  int atm_n0, atm_n1;
  double result2, Rsqrd;
  VECTOR_DOUBLE s_12;
  job->nuc_nuc = k_zero;

  for (atm_n0 = 0; atm_n0 < atoms->number_of_atoms_in_unit_cell; atm_n0++) {
    for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {
         s_12.comp1 = atoms->cell_vector[atm_n0].comp1 - atoms->cell_vector[atm_n1].comp1;
         s_12.comp2 = atoms->cell_vector[atm_n0].comp2 - atoms->cell_vector[atm_n1].comp2;
         s_12.comp3 = atoms->cell_vector[atm_n0].comp3 - atoms->cell_vector[atm_n1].comp3;
         if (atm_n0 != atm_n1) {
         Rsqrd = double_vec_dot(&s_12, &s_12);
         result2 =  k_one / sqrt(Rsqrd);
         job->nuc_nuc += double(atoms->atomic_number[atm_n0] * atoms->atomic_number[atm_n1]) * result2 / two;
        }
       }
      }

  // ******************************************************************************************
  // * Calculate 1e contribution to Fock matrix                                               *
  // ******************************************************************************************

  Function[0] = 1 ;
  Function[1] = 1 ;
  Function[2] = 1 ;
  Function[3] = 0 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  array_dimensions(&dim, &dimg, &pair_p, atoms, job, file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  fock_element_1e2(&one_ints, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);

  // ******************************************************************************************
  // * Check for small eigenvalues of S_k                                                     *
  // ******************************************************************************************

  dim11 = dim1;
  if (job->sym_adapt == 0) {
    n = 0;
    ComplexMatrix *S_k1;
    AllocateComplexMatrix(&S_k1,&dim1,&dim1,job);
      ResetDoubleArray(eigenvalues,&dim1);
      nk[0] = 0;
      nk[1] = 0;
      fourier_transform(&one_ints.Overlap[0], &S_k1->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      DiagonaliseHermitian(&S_k1, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);
      for (ii = 0; ii < atoms->number_of_sh_bfns_in_unit_cell; ii++) {
        if (eigenvalues[ii] < threshold) { 
        n++;
        fprintf(file.out,"small eigenvalue %d %e\n",ii,eigenvalues[ii]); 
       }
      }
    dim11 = dim1 - n;
    DestroyComplexMatrix(&S_k1,job);
    ResetDoubleArray(eigenvalues,&dim1);
   }

  else if (job->sym_adapt == 1) {
    int dim_salc, num_irrep_in_basis[symmetry->number_of_classes];
    int count1;
    double *eigenvalues_S_k[symmetry->number_of_classes];
    ComplexMatrix *S_salc1[symmetry->number_of_classes];
    ComplexMatrix *eigenvectors_S_k[symmetry->number_of_classes];
    count_basis_irrep(num_irrep_in_basis,atoms,atom_p,shells,symmetry,job,file);
    for (iii = 0; iii < symmetry->number_of_classes; iii++) {
      dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
      if (dim_salc == 0) dim_salc = 1; 
      if (num_irrep_in_basis[iii] == 0) continue;
      nk[0] = 0;
      nk[1] = 0;
      AllocateComplexMatrix(&S_salc1[iii],&dim_salc,&dim_salc,job);
      AllocateComplexMatrix(&eigenvectors_S_k[iii],&dim_salc,&dim_salc,job);
      AllocateDoubleArray(&eigenvalues_S_k[iii],&dim_salc,job);
     }

      fourier_transform_salc1(&one_ints.Overlap[0], S_salc1, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job,file);

    for (iii = 0; iii < symmetry->number_of_classes; iii++) {
      if (num_irrep_in_basis[iii] == 0) continue;
      DiagonaliseHermitian(&S_salc1[iii], &eigenvalues_S_k[iii], &eigenvectors_S_k[iii], &jobz, &uplo, &info);
      dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
      if (dim_salc == 0) dim_salc = 1;
      n = 0;
      for (jjj = 0; jjj < dim_salc; jjj++) {
        if (eigenvalues_S_k[iii][jjj] < threshold) { n++; 
        fprintf(file.out,"small eigenvalue irrep %3d %3d %e\n",iii,jjj,eigenvalues_S_k[iii][jjj]); 
       }
       }
        dim_salc2[iii] = dim_salc - n;
        DestroyDoubleArray(&eigenvalues_S_k[iii],&dim_salc,job);
        DestroyComplexMatrix(&eigenvectors_S_k[iii],job);
        DestroyComplexMatrix(&S_salc1[iii],job);
       }
      }

  // ******************************************************************************************
  // * Initialize density matrix mixing                                                       *
  // ******************************************************************************************

  if (job->guess_type == 1) job->mixing_order = 1;
  AllocateDoubleMatrix(&C_i,&job->mixing_order,&dimp_spin,job);
  AllocateDoubleMatrix(&C_o,&job->mixing_order,&dimp_spin,job);
  AllocateDoubleMatrix(&D,&job->mixing_order,&dimp_spin,job);
  int count1, dim_salc, basis_dimension;
  if (job->sym_adapt == 0) basis_dimension = 1;
  else if (job->sym_adapt == 1) basis_dimension = symmetry->number_of_classes;
  int num_irrep_in_basis[basis_dimension];

  if (job->sym_adapt == 0) {
  num_irrep_in_basis[0] = 1;
  count1 = 0;
  for (jjj = 0; jjj < job->spin_dim * job->mixing_order; jjj++) {
    AllocateComplexMatrix(&Residual[count1],&dim11,&dim11,job);
    AllocateComplexMatrix(&FockHist[count1],&dim1,&dim1,job);
    //AllocateComplexMatrix(&FockHist[count1],&dim11,&dim11,job);
    ResetComplexMatrix(Residual[count1]);
    ResetComplexMatrix(FockHist[count1]);
    count1++;
   }
  }

 else if (job->sym_adapt == 1) {
  count_basis_irrep(num_irrep_in_basis,atoms,atom_p,shells,symmetry,job,file);
  count1 = 0;
  for (jjj = 0; jjj < job->spin_dim * job->mixing_order; jjj++) {
    for (iii = 0; iii < symmetry->number_of_classes; iii++) {
      dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
      if (dim_salc == 0) dim_salc = 1;
      AllocateComplexMatrix(&Residual[count1],&dim_salc2[iii],&dim_salc2[iii],job);
      AllocateComplexMatrix(&FockHist[count1],&dim_salc,&dim_salc,job);
      ResetComplexMatrix(Residual[count1]);
      ResetComplexMatrix(FockHist[count1]);
      count1++;
     }
    }
   }

  ResetDoubleMatrix(C_i);
  ResetDoubleMatrix(C_o);
  ResetDoubleMatrix(D);
  for (i = 0; i < dimp_spin; i++)
  delta_P0[i] = P0[i];
  for (i = 0; i < dimf_spin; i++)
  delta_F0[i] = F0[i];

  // ******************************************************************************************
  // * Initialize grid for DFT if needed                                                      *
  // ******************************************************************************************

//  DFT_GRID dft_grid;
//  if (job->xc_num > 0) {
//  time1 = MPI_Wtime();
//  count_dft_savin_grid(&ngrid_points,&pair_p,atoms,crystal,R,symmetry,job,file);
//  allocate_dft_grid(&dft_grid, ngrid_points,atoms,job,file);
//  generate_dft_savin_grid(&dft_grid,&pair_p,atoms,crystal,R,symmetry,job,file);
//  time2 = MPI_Wtime();
//  if (job->taskid == 0)
//  printf("Time to generate %7d dft_grid points %10.4e\n",ngrid_points,(double)(time2 - time1));
// }
//  if (job->xc_num == 0) {
//  ngrid_points = 1;
//  allocate_dft_grid(&dft_grid, ngrid_points, atoms, job, file);
// }

  // ******************************************************************************************
  // * SCF loop begins here                                                                   *
  // ******************************************************************************************

  time3 = MPI_Wtime();

  job->guess_type = 0;

  double *S1;
  AllocateDoubleArray(&S1,&job->dimf,job);
  if (job->scf_coulomb  == 1 || job->scf_exchange == 1)
  fock_matrix_molecule_compute_screening_integrals(S1,&pair_p,R,G,atoms,shells,gaussians,symmetry,crystal,job,file);

  for (job->iter = 1; job->iter <= job->max_cycle; job->iter++) {

    if (job->taskid == 0) 
    printf("ITERATION %d\n",job->iter);

    ResetIntArray(irrep,&dim3);
    ResetDoubleArray(eigval,&dim3);

    time6 = MPI_Wtime();
    time7 = k_zero;

  // ******************************************************************************************
  // * Calculate 2e contribution to Fock matrix                                               *
  // ******************************************************************************************

   time2 = MPI_Wtime();
   build_fock_matrix_molecule(Fock,&one_ints,fermi,&total_electrons,S1,P0,F0,delta_P0,delta_F0,&pair_p,atom_p,atoms, \
   shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
   if (job->taskid == 0) printf("\n");
   time1 = MPI_Wtime() - time2;

   time2 = MPI_Wtime();
   //build_fock_matrix_molecule_no_sym(Fock,&one_ints,fermi,&total_electrons,S1,P0,F0,delta_P0,delta_F0,&pair_p,atom_p,atoms, \
   shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
   if (job->taskid == 0) printf("\n");
   time1 = MPI_Wtime() - time2;

   if (job->taskid == 0 && job->mpp == 0) scf_evectors = fopen(zz2, "wb");
   if (job->values == 1 || job->values == 2)
   MPI_File_open(MPI_COMM_WORLD,buf1,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&eh) ;
   if (job->vectors == 1 || job->vectors == 2)
   MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
   if (job->density == 1 || job->density == 2)
   MPI_File_open(MPI_COMM_WORLD,buf4,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&gh) ;

  // ******************************************************************************************
  // * Analysis of atomic populations and density matrix normalization                        *
  // ******************************************************************************************

   atom_shell_populations2(&one_ints, P0, &pair_p, atoms, shells, symmetry, job, file);
   total_population = print_total_population2(atoms, shells, job, file);
   fac = k_one;
   total_energy(&one_ints,Fock,P0,&pair_p,atoms,job,file);

   if (job->taskid >= 0) {
   fprintf(file.out,"===========================================================================================================\n");
   fprintf(file.out,"| ITERATION             %3d | ENERGY %16.9e | ENERGY CHANGE %9.2e | DENSITY NORM  %9.6lf |\n", \
   job->iter, job->total_energy, job->energy_change, fac);
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   print_atom_populations2(atoms, shells, job, file);
   if (job->verbosity > 0) print_shell_populations2(atoms, shells, job, file);
   if (job->xc_num > 0) {
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"| %12.6e ELECTRONS    |                         |                         |                         |\n", \
   total_electrons);
  }
 }

  // ******************************************************************************************
  // * Loop over k points                                                                     *
  // ******************************************************************************************

  if (job->mpp == 1)
  ResetDoubleArray(P2,&dimp_spin);

  k = 0;
  nk[0] = knet.ibz[k];
  nk[1] = knet.ibz[k];

  if (job->taskid == 0) {

  // ******************************************************************************************
  // * Work in non-symmetry adapted basis                                                     *
  // ******************************************************************************************

  if (job->sym_adapt == 0) {

  // ******************************************************************************************
  // * Find transformation to orthogonalise basis at this k point and store in xtrn1          *
  // ******************************************************************************************

 ComplexMatrix *Cholesky_S_k_inverse;
 int *Cholesky_S_k_pivot;

 job->scf_trans = 0; // set to canonical transformation = 1 is Cholesky

 ComplexMatrix *xtrn1, *F_k0, *F_k1, *P0_k0, *S_k, *xtmp1, *eigenvectors1;
 AllocateComplexMatrix(&S_k,&dim1,&dim1,job);

 fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);

 if (job->scf_trans == 0) { // canonical transformation

 DiagonaliseHermitian(&S_k, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);

 n = 0;
 for (ii = 0; ii < atoms->number_of_sh_bfns_in_unit_cell; ii++) {
 if (job->iter == 1 & eigenvalues[ii] < threshold) { 
 fprintf(file.out,"small eigenvalue %d %d %e\n",k,ii,eigenvalues[ii]); 
   }
    if (eigenvalues[ii] < threshold)  n++; 
   }
    dim11 = dim1 - n;
    if (job->kpoints == 0 && dim11 < fermi->occupied[k]) {
    if (job->taskid == 0)
    fprintf(file.out,"Limited basis set is also linearly dependent at k-point %3d. Stopping.\n",k);
    MPI_Finalize();
    exit(1);
   }

    AllocateComplexMatrix(&xtrn1,&dim11,&dim1,job);

    for (i = n; i < dim1; i++) { 
      for (j = 0; j < dim1; j++) { 
	xtrn1->a[i - n][j] = eigenvectors->a[i][j] / sqrt(*(eigenvalues + i));
       }
      }

  // ******************************************************************************************
  // * Check that overlap matrix is correctly Fourier transformed                             *
  // ******************************************************************************************

  if (job->taskid == 0 && job->verbosity > 1) {
    ComplexMatrix *S_k1, *S_k2;
    AllocateComplexMatrix(&xtmp1,&dim11,&dim1,job);
    AllocateComplexMatrix(&S_k1,&dim1,&dim1,job);
    AllocateComplexMatrix(&S_k2,&dim11,&dim11,job);
    fourier_transform(&one_ints.Overlap[0], &S_k1->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
    ResetComplexMatrix(xtmp1);
    ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1, &S_k1, &beta, &xtmp1);
    ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1, &xtrn1, &beta, &S_k2);
    fprintf(file.out,"S_k1 orthogonalised %d\n",k);
    print_complex_matrix(S_k1, file);
    print_complex_matrix(S_k2, file);
    DestroyComplexMatrix(&S_k1,job);
    DestroyComplexMatrix(&S_k2,job);
    DestroyComplexMatrix(&xtmp1,job);
  }

  } // close if (job->scf_trans == 0)

    else if (job->scf_trans == 1) { // cholesky
      AllocateIntArray(&Cholesky_S_k_pivot,&dim1,job);
      CholeskyHermitian(&S_k,Cholesky_S_k_pivot,&dim11,&threshold);
      AllocateComplexMatrix(&Cholesky_S_k_inverse,&dim11,&dim1,job);
      CholeskyInverseHermitian(&S_k,&Cholesky_S_k_inverse,&dim11);
     }

    AllocateComplexMatrix(&xtmp1,&dim11,&dim1,job);
    AllocateComplexMatrix(&F_k0,&dim1,&dim1,job);
    AllocateComplexMatrix(&F_k1,&dim11,&dim11,job);
    AllocateComplexMatrix(&eigenvectors1,&dim11,&dim11,job);
    ResetComplexMatrix(eigvec);

  for (s = 0; s < job->spin_dim; s++) {

    fourier_transform(&Fock[s * dimp], &F_k0->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);

    if (job->scf_trans == 0) {

    if (job->diis == 1) {
      AllocateComplexMatrix(&P0_k0,&dim1,&dim1,job);
      fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      fourier_transform(&P0[s * dimp], &P0_k0->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      spin_offset = s * job->mixing_order * symmetry->number_of_classes;
      C_DIIS_extrapolation(&F_k0, &P0_k0, &S_k, &xtrn1, &FockHist[spin_offset], &Residual[spin_offset], &dim1, &dim11, job, file);
      DestroyComplexMatrix(&P0_k0,job);
     }
      ResetComplexMatrix(xtmp1);
      ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1, &F_k0, &beta, &xtmp1);
      ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1, &xtrn1, &beta, &F_k1);

      if (job->taskid == 0 && job->verbosity > 1) {
      fprintf(file.out,"Xtrn k %3d spin %3d\n",k,s);
      print_complex_matrix(xtrn1, file);
      fprintf(file.out,"F_k orthogonalised k %3d spin %3d\n",k,s);
      print_complex_matrix(F_k1, file);
     }
    } 

    else if (job->scf_trans == 1) {
   //   dim11 = 0;
   //   for (i = 0; i < dim1; i++) {
   //   if ((S_k->a[i][i]).real() > threshold) dim11++;
   //    }
      ResetComplexMatrix(xtmp1);
      CholeskyPermuteHermitian(&F_k0,Cholesky_S_k_pivot,0,job);
      ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &Cholesky_S_k_inverse, &F_k0, &beta, &xtmp1);
      ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1, &Cholesky_S_k_inverse, &beta, &F_k1);
     }

    time2 = MPI_Wtime();
    ResetComplexMatrix(eigenvectors1);
    p_irrep  =  &irrep[(s * fermi->nkunique + k) * dim1];
    p_eigval = &eigval[(s * fermi->nkunique + k) * dim1];
    DiagonaliseHermitian(&F_k1, &p_eigval, &eigenvectors1, &jobz, &uplo, &info);
    for (i = dim11; i < dim1; i++) p_eigval[i] = 9999.9; // set zero eigenvalues to large number
    for (i = dim11; i < dim1; i++) p_irrep[i] = 0; // only one irrep for no symmetry adaptation

    if (job->scf_trans == 0) {
      ResetComplexMatrix(eigvec);
      ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigenvectors1, &xtrn1, &beta, &eigvec);
     }

    else if (job->scf_trans == 1) {
      ResetComplexMatrix(eigvec);
      ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigenvectors1, &Cholesky_S_k_inverse, &beta, &eigvec);
      CholeskyPermuteHermitian(&eigvec,Cholesky_S_k_pivot,2,job);
     }

    if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"k point %d eigenvalues of F_k\n",k);
    print_complex_eigenvector_matrix3(eigenvectors1, eigenvalues, file);
    fprintf(file.out,"transformed eigenvectors of F_k\n");
    print_complex_matrix(eigvec, file);
   }

    if (job->vectors == 2 && job->mpp == 0) {
    if (crystal->type[0] == 'M' && job->taskid == 0) {
    fwrite(&eigvec->a[fermi->bands[0] - 1][0], sizeof(double), 2 * vector_size, scf_evectors);
    fflush(scf_evectors);
   }
  }

    else if (job->mpp == 1) {
    time2 = MPI_Wtime();
    reduced_density_matrix_molecule_mpp(s,p_irrep,p_eigval,eigvec,fermi,&P2[s * job->dimp],&knet,&fermi->nkunique,R,\
    &pair_p,atom_p,atoms,shells,symmetry,job,file);
    time7 += MPI_Wtime() - time2;
   }

 } // close loop on s

      if (job->scf_trans == 0) {
      DestroyComplexMatrix(&xtrn1,job);
     }
      else if (job->scf_trans == 1) {
      DestroyIntArray(&Cholesky_S_k_pivot,&dim1,job);
      DestroyComplexMatrix(&Cholesky_S_k_inverse,job);
     }
      DestroyComplexMatrix(&xtmp1,job);
      DestroyComplexMatrix(&F_k0,job);
      DestroyComplexMatrix(&F_k1,job);
      DestroyComplexMatrix(&eigenvectors1,job);
      DestroyComplexMatrix(&S_k,job);

    } // close if (job->sym_adapt == 0)

  // ******************************************************************************************
  // * End work in non-symmetry adapted basis                                                 *
  // ******************************************************************************************

  // ******************************************************************************************
  // * Work in symmetry adapted basis                                                         *
  // ******************************************************************************************

    else if (job->sym_adapt == 1) {

    int n1, dim_salc, dim_salc1, num_irrep_in_basis[symmetry->number_of_classes];
    int num_eigval, num_eigval1;
    int *Cholesky_S_k_pivot[symmetry->number_of_classes];
    double *eigenvalues_S_k[symmetry->number_of_classes], *eigenvalues_salc[symmetry->number_of_classes];
    double eigvec_bas;

    job->scf_trans = 0; // use canonical transformation

    count_basis_irrep(num_irrep_in_basis,atoms,atom_p,shells,symmetry,job,file);
    ComplexMatrix *Cholesky_S_k[symmetry->number_of_classes];
    ComplexMatrix *Cholesky_S_k_inverse[symmetry->number_of_classes];
    ComplexMatrix *eigenvectors_S_k[symmetry->number_of_classes];
    ComplexMatrix *F_k0_salc[symmetry->number_of_classes];
    ComplexMatrix *F_k1_salc[symmetry->number_of_classes];
    ComplexMatrix *S_salc1[symmetry->number_of_classes];
    ComplexMatrix *xtmp1_salc[symmetry->number_of_classes];
    ComplexMatrix *eigenvectors_salc[symmetry->number_of_classes];
    ComplexMatrix *eigvec_salc[symmetry->number_of_classes];
    ComplexMatrix *xtrn1_salc[symmetry->number_of_classes];

    ComplexMatrix *P0_salc[symmetry->number_of_classes];
    ComplexMatrix *S0_salc[symmetry->number_of_classes];

  // ******************************************************************************************
  // * Allocate memory                                                                        *
  // ******************************************************************************************

    for (iii = 0; iii < symmetry->number_of_classes; iii++) {
    dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
    if (dim_salc == 0) dim_salc = 1;
    AllocateComplexMatrix(&S_salc1[iii],&dim_salc,&dim_salc,job);
    AllocateComplexMatrix(&Cholesky_S_k[iii],&dim_salc,&dim_salc,job);
    AllocateComplexMatrix(&eigenvectors_S_k[iii],&dim_salc,&dim_salc,job);
    AllocateComplexMatrix(&eigvec_salc[iii],&dim_salc,&dim_salc,job);
    AllocateComplexMatrix(&F_k0_salc[iii],&dim_salc,&dim_salc,job);
    AllocateDoubleArray(&eigenvalues_salc[iii],&dim_salc,job);
    AllocateDoubleArray(&eigenvalues_S_k[iii],&dim_salc,job);
   }
  
  // ******************************************************************************************
  // * scf_trans == 0/1 Canonical/Cholesky transformation of AO basis                         *
  // ******************************************************************************************

    if (job->scf_trans == 0) 
    fourier_transform_salc1(&one_ints.Overlap[0], S_salc1, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job,file);
    if (job->scf_trans == 1) 
    fourier_transform_salc1(&one_ints.Overlap[0], Cholesky_S_k, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job,file);

  // ******************************************************************************************
  // * Loop over symmetry classes                                                             *
  // ******************************************************************************************

    for (iii = 0; iii < symmetry->number_of_classes; iii++) {
      dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
      if (dim_salc == 0) dim_salc = 1;
      if (num_irrep_in_basis[iii] == 0) continue;

  // ******************************************************************************************
  // * Canonical transformation                                                               *
  // ******************************************************************************************

      if (job->scf_trans == 0) { 

	DiagonaliseHermitian(&S_salc1[iii], &eigenvalues_S_k[iii], &eigenvectors_S_k[iii], &jobz, &uplo, &info);
	n1 = 0;
	for (jjj = 0; jjj < dim_salc; jjj++) {
	  //if (eigenvalues_S_k[iii][jjj] < threshold && job->taskid == 0)  \
	  fprintf(file.out,"small eigenvalue %3d irrep %3d eigenvalue %3d %e\n",k,iii,jjj,eigenvalues_S_k[iii][jjj]); 
	  if (eigenvalues_S_k[iii][jjj] < threshold)  n1++; 
	 }
	  dim_salc1 = dim_salc - n1;
	  AllocateComplexMatrix(&xtrn1_salc[iii],&dim_salc1,&dim_salc,job);
	  AllocateComplexMatrix(&xtmp1_salc[iii],&dim_salc1,&dim_salc,job);
	  AllocateComplexMatrix(&F_k1_salc[iii],&dim_salc1,&dim_salc1,job);
	  AllocateComplexMatrix(&eigenvectors_salc[iii],&dim_salc1,&dim_salc1,job);
	  for (i = n1; i < dim_salc; i++) { 
	    for (j = 0; j < dim_salc; j++) { 
	      xtrn1_salc[iii]->a[i - n1][j] = eigenvectors_S_k[iii]->a[i][j] / sqrt(*(eigenvalues_S_k[iii] + i));
	     }
	    }

	   }

  // ******************************************************************************************
  // * Cholesky transformation                                                                *
  // ******************************************************************************************

      else if (job->scf_trans == 1) { 
	AllocateIntArray(&Cholesky_S_k_pivot[iii],&dim_salc,job);
	CholeskyHermitian(&Cholesky_S_k[iii],Cholesky_S_k_pivot[iii],&dim_salc1,&threshold);
	AllocateComplexMatrix(&xtrn1_salc[iii],&dim_salc1,&dim_salc,job);
	AllocateComplexMatrix(&Cholesky_S_k_inverse[iii],&dim_salc1,&dim_salc,job);
	CholeskyInverseHermitian(&Cholesky_S_k[iii],&Cholesky_S_k_inverse[iii],&dim_salc1);
	AllocateComplexMatrix(&xtmp1_salc[iii],&dim_salc1,&dim_salc,job);
	AllocateComplexMatrix(&F_k1_salc[iii],&dim_salc1,&dim_salc1,job);
	AllocateComplexMatrix(&eigenvectors_salc[iii],&dim_salc1,&dim_salc1,job);
       }

      } // close loop on iii

  // ******************************************************************************************
  // * Loop over spin if spin-polarized                                                       *
  // ******************************************************************************************

  for (s = 0; s < job->spin_dim; s++) {

    fourier_transform_salc1(&Fock[s * dimp], F_k0_salc, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job, file);

  // ******************************************************************************************
  // * DIIS Fock matrix mixing in orthogonalised salc basis                                   *
  // ******************************************************************************************

    if (job->diis == 1) {
      for (iii = 0; iii < symmetry->number_of_classes; iii++) {
        dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
        if (dim_salc == 0) dim_salc = 1;
        AllocateComplexMatrix(&P0_salc[iii],&dim_salc,&dim_salc,job);
        AllocateComplexMatrix(&S0_salc[iii],&dim_salc,&dim_salc,job);
       }
        fourier_transform_salc1(&P0[s * dimp], P0_salc, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job, file);
        fourier_transform_salc1(&one_ints.Overlap[0], S0_salc, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job,file);
        spin_offset = s * job->mixing_order * symmetry->number_of_classes;
        C_DIIS_extrapolation_salc(F_k0_salc, P0_salc, S0_salc, xtrn1_salc, &FockHist[spin_offset], &Residual[spin_offset], num_irrep_in_basis, &dim1, dim_salc2,\
        symmetry, job,file);
        //C_DIIS_extrapolation_salc(F_k0_salc, P0_salc, S0_salc, xtrn1_salc, FockHist, Residual, num_irrep_in_basis, &dim1, dim_salc2,\
        symmetry, job,file);
       } 

      num_eigval = 0; num_eigval1 = 0;

      for (iii = 0; iii < symmetry->number_of_classes; iii++) {
	if (num_irrep_in_basis[iii] == 0) continue;
	dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
	ResetComplexMatrix(xtmp1_salc[iii]);
	ResetComplexMatrix(F_k1_salc[iii]);

  // ******************************************************************************************
  // * Canonical transformation                                                               *
  // ******************************************************************************************

	if (job->scf_trans == 0) {
	  n1 = 0;
	  for (jjj = 0; jjj < dim_salc; jjj++) {
	  if (eigenvalues_S_k[iii][jjj] < threshold)  n1++; 
	 }
	  dim_salc1 = dim_salc - n1;
	  ResetComplexMatrix(xtmp1_salc[iii]);
	  ResetComplexMatrix(F_k1_salc[iii]); 
          //print_complex_matrix2(F_k0_salc[iii],0,6,1.0,file);
	  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1_salc[iii], &F_k0_salc[iii], &beta, &xtmp1_salc[iii]);
	  ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1_salc[iii], &xtrn1_salc[iii], &beta, &F_k1_salc[iii]);
	 }

  // ******************************************************************************************
  // * Cholesky transformation                                                                *
  // ******************************************************************************************

	else if (job->scf_trans == 1) {
	  dim_salc1 = 0;
	  for (i = 0; i < dim_salc; i++) {
	  if ((Cholesky_S_k[iii]->a[i][i]).real() > threshold) dim_salc1++;
	 }
	  CholeskyPermuteHermitian(&F_k0_salc[iii],Cholesky_S_k_pivot[iii],0,job);
	  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &Cholesky_S_k_inverse[iii], &F_k0_salc[iii], &beta, &xtmp1_salc[iii]);
	  ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1_salc[iii], &Cholesky_S_k_inverse[iii], &beta, &F_k1_salc[iii]);
	 }

  // ******************************************************************************************
  // * Diagonalise Fock matrix                                                                *
  // ******************************************************************************************

	DiagonaliseHermitian(&F_k1_salc[iii], &eigenvalues_salc[iii], &eigenvectors_salc[iii], &jobz, &uplo, &info);

  // ******************************************************************************************
  // * Canonical back-transformation                                                          *
  // ******************************************************************************************

	if (job->scf_trans == 0) {
	  ResetComplexMatrix(eigvec_salc[iii]);
	  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigenvectors_salc[iii], &xtrn1_salc[iii], &beta, &eigvec_salc[iii]);
	 }

  // ******************************************************************************************
  // * Cholesky back-transformation                                                           *
  // ******************************************************************************************

	else if (job->scf_trans == 1) {

	  ResetComplexMatrix(eigvec_salc[iii]);
	  ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigenvectors_salc[iii], &Cholesky_S_k_inverse[iii], &beta, &eigvec_salc[iii]);
	  CholeskyPermuteHermitian(&eigvec_salc[iii],Cholesky_S_k_pivot[iii],2,job);

	 }

  // ******************************************************************************************
  // * Transform degenerate irreps (n = 2 only) to standard form                              *
  // ******************************************************************************************

	if (symmetry->irp_dim_k[iii] == 2) 
	transform_degenerate_salc(eigvec_salc, eigenvalues_salc, &iii, num_irrep_in_basis, &knet, nk, symmetry, job, file);

  // ******************************************************************************************
  // * Put Fock eigenvalues in eigval array                                                   *
  // ******************************************************************************************

        p_irrep  =  &irrep[(s * fermi->nkunique + k) * dim1];
        p_eigval = &eigval[(s * fermi->nkunique + k) * dim1];
        for (jjj = 0; jjj < dim_salc1; jjj++) p_eigval[num_eigval + jjj] = eigenvalues_salc[iii][jjj];
        for (jjj = 0; jjj < dim_salc; jjj++)  p_irrep[num_eigval + jjj] = iii;
	  for (i = 0; i < dim_salc1; i++) { 
	    eigvec_bas = k_zero;
	    for (j = 0; j < dim_salc1; j++) { 
	      eigvec_bas += fabs((eigvec_salc[iii]->a[i][j]).real());
	     }
	    }
         for (jjj = dim_salc1; jjj < dim_salc; jjj++) p_eigval[num_eigval + jjj] = 9999.9; // set zero eigenvalues to large number
         //for (jjj = 0; jjj < dim_salc; jjj++) fprintf(file.out,"eigenval salc %3d irrep %3d %10.4f\n",\
	 jjj,p_irrep[num_eigval+jjj],eigenvalues_salc[iii][jjj]);
	 num_eigval  += dim_salc;
	 num_eigval1 += dim_salc1;

   } // close loop on iii (symmetry classes)

    if (job->kpoints == 0 && num_eigval1 < fermi->occupied[k]) {
    if (job->taskid == 0)
    fprintf(file.out,"Limited basis set is also linearly dependent at k-point %3d. Stopping.\n",k);
    MPI_Finalize();
    exit(1);
   }

  // ******************************************************************************************
  // * Transform eigenvectors from salc basis to AO basis                                     *
  // ******************************************************************************************

    ResetComplexMatrix(eigvec);
    transform_salc_to_atomic_basis(eigvec,p_eigval,p_irrep,eigvec_salc,&knet,nk,atom_p,&pair_p,R,atoms,shells,symmetry,job,file);
    //print_complex_eigenvector_matrix2(eigvec, p_eigval, 6, num_eigval, 1.0, file);
    //for (iii = 0; iii < dim_salc; iii++) fprintf(file.out,"iii irrep eigval %3d %3d %10.4f\n",iii + 1,p_irrep[iii], au_to_eV * p_eigval[iii]);

  // ******************************************************************************************
  // * Write eigenvectors in AO basis to disk or update density matrix (job->mpp == 0/1)      *
  // ******************************************************************************************

    if (job->vectors == 2 && job->taskid == 0 && job->mpp == 0) {
    fwrite(&eigvec->a[fermi->bands[0] - 1][0], sizeof(double), 2 * vector_size, scf_evectors);
    fflush(scf_evectors);
   }

    else if (job->mpp == 1) {
    time2 = MPI_Wtime();
    reduced_density_matrix_molecule_mpp(s,p_irrep,p_eigval,eigvec,fermi,&P2[s * job->dimp],&knet,\
    &fermi->nkunique,R,&pair_p,atom_p,atoms,shells,symmetry,job,file);
    time7 += MPI_Wtime() - time2;
   }

   if (job->diis == 1) {
     for (iii = 0; iii < symmetry->number_of_classes; iii++) {
       if (num_irrep_in_basis[iii] == 0) continue;
       DestroyComplexMatrix(&P0_salc[iii],job);
       DestroyComplexMatrix(&S0_salc[iii],job);
      }
     }

 } // close loop on s

   for (iii = 0; iii < symmetry->number_of_classes; iii++) {
     if (num_irrep_in_basis[iii] == 0) continue;
      dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
      if (job->scf_trans == 1) {
      DestroyIntArray(&Cholesky_S_k_pivot[iii],&dim_salc,job);
      DestroyComplexMatrix(&Cholesky_S_k_inverse[iii],job);
     }
      DestroyDoubleArray(&eigenvalues_salc[iii],&dim_salc,job);
      DestroyDoubleArray(&eigenvalues_S_k[iii],&dim_salc,job);
      DestroyComplexMatrix(&Cholesky_S_k[iii],job);
      DestroyComplexMatrix(&eigenvectors_S_k[iii],job);
      DestroyComplexMatrix(&xtmp1_salc[iii],job);
      DestroyComplexMatrix(&F_k0_salc[iii],job);
      DestroyComplexMatrix(&F_k1_salc[iii],job);
      DestroyComplexMatrix(&S_salc1[iii],job);
      DestroyComplexMatrix(&eigenvectors_salc[iii],job);
      DestroyComplexMatrix(&eigvec_salc[iii],job);
      DestroyComplexMatrix(&xtrn1_salc[iii],job);
     }

    } // close if (job->sym_adapt == 1)

  // ******************************************************************************************
  // * End work in symmetry adapted basis                                                     *
  // ******************************************************************************************

  } // close if (job->taskid == 0) 

  // k loop closes here in periodic systems

  // ******************************************************************************************
  // * Gather eigenvalues and write to disk if job->values == 2                               *
  // ******************************************************************************************

    ResetDoubleArray(eigval1,&dim3);
    MPI_Allreduce(eigval,eigval1,dim3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    job->values = 2;

    if (job->values == 2) {
      MPI_File_seek(eh, 0, MPI_SEEK_SET) ;
      for (k = 0; k < knet.unique; k++) {
	for (s = 0; s < job->spin_dim; s++) {
	  MPI_File_write(eh, &eigval1[(s * knet.unique + k) * dim1 + fermi->bands[0]-1],value_size,MPI_DOUBLE,MPI_STATUS_IGNORE);
	 }
	}
       }

  // ******************************************************************************************
  // * Calculate fermi->occupied array for spin-polarised molecule                            *
  // ******************************************************************************************

     if (job->spin_polarisation == 1 && crystal->type[0] == 'M' && job->iter > 250) {
       int count1, count2, total_occupied = fermi->occupied[0] + fermi->occupied[1];
       count1 = 0;
       count2 = 0;
       fprintf(file.out,"fermi->occ %3d %3d\n",fermi->occupied[0],fermi->occupied[1]);
       fermi->occupied[0] = 0;
       fermi->occupied[1] = 0;
       do {
         if (eigval1[count1] < eigval1[dim1 + count2]) {
           fermi->occupied[0]++;
           count1++;
          }
         else {
           fermi->occupied[1]++;
           count2++;
          }
           //printf("i %3d s %3d totocc %3d occ[0] %3d occ[1] %3d eigval1 %10.4f eigval2 %10.4f\n",\
           i,s,count1+count2,fermi->occupied[0],fermi->occupied[1],eigval1[count1],eigval1[dim1 + count2]);
         } while (count1 + count2 < total_occupied);
        }

    if (job->taskid == 0 && job->vectors == 2 && crystal->type[0] == 'M') {
    if      (job->spin_dim == 1) printf("Fermi->occupied = %3d\n",fermi->occupied[0]);
    else if (job->spin_dim == 2) printf("Fermi->occupied = %3d %3d\n",fermi->occupied[0], fermi->occupied[1]);
    job->fermi_energy = eigenvalues[fermi->occupied[0] - 1];
   }

    for (i = 0; i < dim3; i++) {
    job_parameters1[11 + i] = fermi->occupation[i];
    //fprintf(file.out,"fermi %d %lf\n",i,fermi->occupation[i]);
   }
    for (i = 0; i < job->spin_dim * fermi->nkunique; i++) {
    job_parameters[41 + i] = fermi->occupied[i];
    //fprintf(file.out,"fermi %d %d\n",i,fermi->occupied[i]);
   }

  // ******************************************************************************************
  // * Calculate new density matrix if needed and write to disk                               *
  // ******************************************************************************************

    if (job->vectors == 2 && (job->density == 1 || job->density == 2)) {

     ResetDoubleArray(P1,&dimp_spin);
     ResetDoubleArray(F1,&dimf_spin);

    if (job->mpp == 0) {
     time2 = MPI_Wtime();
     ResetDoubleArray(P2,&dimp_spin);
     reduced_density_matrix_molecule(p_irrep,p_eigval,eigvec,fermi,P2,&knet,&fermi->nkunique,R,&pair_p,atom_p,atoms,shells,symmetry,job,file);
     MPI_Allreduce(P2,P1,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
     expand_density_matrix(P1,F1,&pair_p,atoms,shells,symmetry,job,file);
     time1 = MPI_Wtime() - time2;
     if (job->taskid == 0)
     printf("Time to generate density matrix %10.4e\n",(double)(time1));
    }

    else if (job->mpp == 1) {
     MPI_Allreduce(P2,P1,dimp_spin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }

  // ******************************************************************************************
  //  * Obtain idempotent density matrix by accepting new density matrix: 1st iter only       *
  // ******************************************************************************************

     double *P0_copy, *F0_copy;
     AllocateDoubleArray(&P0_copy,&dimp_spin,job);
     AllocateDoubleArray(&F0_copy,&dimf_spin,job);
     for (i = 0; i < dimp_spin; i++)
     P0_copy[i] = P0[i];
     for (i = 0; i < dimf_spin; i++)
     F0_copy[i] = F0[i];

     if (job->iter == 1 || (job->diis == 1 && fabs(job->energy_change) < 1e-02)) { 

     for (i = 0; i < dimp_spin; i++)
     P0[i] = P1[i];

     //P0[i] = job->fock_mixing * P0[i] + (k_one - job->fock_mixing) * P1[i];
     for (i = 0; i < dimf_spin; i++) 
     F0[i] = F1[i];
     //F0[i] = job->fock_mixing * F0[i] + (k_one - job->fock_mixing) * F1[i];

    }

  // ******************************************************************************************
  // * Simple mixing                                                                          *
  // ******************************************************************************************

     //if (fabs(job->energy_change) > 5.0e-01) {
     //else if (job->iter > 1) {
     //for (i = 0; i < dimp_spin; i++)
     //P0[i] = job->fock_mixing * P0[i] + (k_one - job->fock_mixing) * P1[i];
     //for (i = 0; i < dimf_spin; i++) 
     //F0[i] = job->fock_mixing * F0[i] + (k_one - job->fock_mixing) * F1[i];
    //}

  // ******************************************************************************************
  // * Pulay mixing                                                                           *
  // ******************************************************************************************

     else if (job->iter > 1 && (job->diis == 0 || fabs(job->energy_change) > 1e-02)) {

     init_mix(C_i, C_o, D, P0, P1);
     mix_density(C_i, C_o, D, P0, job, file);

    }

  // ******************************************************************************************
  // * Calculate full density matrix and density matrix changes                               *
  // ******************************************************************************************

     expand_density_matrix(P0,F0,&pair_p,atoms,shells,symmetry,job,file);

     for (i = 0; i < dimp_spin; i++) {
     delta_P0[i] = P0[i] - P0_copy[i];
    }

     ResetDoubleArray(delta_F0,&dimf_spin);
     expand_density_matrix(delta_P0,delta_F0,&pair_p,atoms,shells,symmetry,job,file);

     for (i = 0; i < dimf_spin; i++) delta_F0[i] = F0[i] - F0_copy[i];
     DestroyDoubleArray(&P0_copy,&dimp_spin,job);
     DestroyDoubleArray(&F0_copy,&dimf_spin,job);
     //print_Fock_matrix_molecule(0,delta_P0, &pair_p, atoms, job, file);
     //print_Fock_matrix_molecule(1,delta_F0, &pair_p, atoms, job, file);
     //print_density_matrix_molecule(0, P0, &pair_p, atoms, job, file);
     //print_density_matrix_molecule(1, F0, &pair_p, atoms, job, file);

     if (job->mpp == 1) MPI_Bcast(P0,dimp_spin,MPI_DOUBLE,0,MPI_COMM_WORLD);

   }

    if (job->density == 2) {
    MPI_File_seek(gh, 0, MPI_SEEK_SET) ;
    MPI_File_write(gh, job_parameters, 41 + job->spin_dim * fermi->nkunique, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write(gh, job_parameters1, 11 + dim3, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_write(gh, &P0[0], dimp_spin, MPI_DOUBLE, MPI_STATUS_IGNORE);
   }

  // ******************************************************************************************
  // * Close files                                                                            *
  // ******************************************************************************************

    if (job->taskid == 0 && job->mpp == 0) fclose(scf_evectors);

    if (job->values  == 1 || job->values  == 2) MPI_File_close(&eh);
    if (job->vectors == 1 || job->vectors == 2) MPI_File_close(&fh);
    if (job->density == 1 || job->density == 2) MPI_File_close(&gh);

    time5 = MPI_Wtime() - time6;

    MPI_Bcast(&job->energy_change,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if ((job->vectors == 2 && job->density == 2) && (fabs(job->energy_change) < job->scf_tol)) break;

 } // close loop on job->iter

  // ******************************************************************************************
  // * Final Fock/Kohn-Sham matrix build in SCF run                                           *
  // ******************************************************************************************

    if (job->vectors == 2 && job->density == 2) { // only scf does this
    if (job->taskid == 0) 
    printf("ITERATION %d\n",job->iter + 1);
    build_fock_matrix_molecule(Fock,&one_ints,fermi,&total_electrons,S1,P0,F0,delta_P0,delta_F0,&pair_p,atom_p,atoms, \
    shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
    time4 = MPI_Wtime();
    if (job->taskid == 0)
    printf("\nTime for SCF calculation %10.4e\n\n",(double)(time4 - time3));
   }

  // ******************************************************************************************
  // * Write SCF eigenvectors to scf_evec MPI file                                            *
  // ******************************************************************************************

  read_write_SCF_eigenvectors(fermi,eigvec,atoms,job,file);

  DestroyDoubleArray(&S1,&job->dimf,job);
   
  //free_dft_grid(&dft_grid,job);
  DestroyDoubleMatrix(&C_i,job);
  DestroyDoubleMatrix(&C_o,job);
  DestroyDoubleMatrix(&D,job);
  if (job->sym_adapt == 1) {
  for (iii = 0; iii < job->spin_dim * job->mixing_order * symmetry->number_of_classes; iii++) {
  DestroyComplexMatrix(&Residual[iii],job);
  DestroyComplexMatrix(&FockHist[iii],job);
 }
 }

  DestroyComplexMatrix(&eigenvectors,job);
  DestroyComplexMatrix(&eigvec,job);
  DestroyDoubleArray(&P2,&dimp_spin,job);
  DestroyDoubleArray(&eigenvalues,&dim1,job);
  DestroyDoubleArray(&eigval,&dim3,job);
  DestroyDoubleArray(&eigval1,&dim3,job);
  DestroyDoubleArray(&Fock,&dimp_spin,job);
  DestroyDoubleArray(&P0,&dimp_spin,job);
  DestroyDoubleArray(&F0,&dimf_spin,job);
  DestroyDoubleArray(&delta_P0,&dimp_spin,job);
  DestroyDoubleArray(&delta_F0,&dimf_spin,job);
  DestroyDoubleArray(&P1,&dimp_spin,job);
  DestroyDoubleArray(&F1,&dimf_spin,job);
  DestroyIntArray(&irrep,&dim3,job);
  free_k_points(&knet,job);
  free_PAIR_TRAN(&pair_p,job);
  free_INT_1E(&one_ints, Function, job, file);

}

void C_DIIS_extrapolation_salc(ComplexMatrix **F_k0_salc, ComplexMatrix **P0_salc, ComplexMatrix **S_k_salc, ComplexMatrix **xtrn_salc, ComplexMatrix **FockHist, ComplexMatrix **Residual, int *num_irrep_in_basis, int *dim1, int *dim_salc2, SYMMETRY *symmetry, JOB_PARAM *job, FILES file) 

{

int i, j, ii, jj, kk, iii;
int info = 0, ione = 1;
int dim_salc;
int nmix = (job->mixing_order < job->iter) ? job->mixing_order : job->iter;
int nmix1 = nmix + 1;
int end, mix_or_not;
int *work;
double *B;
double trace;
char NoTrans = 'N', ConjTrans = 'C';
char uplo = 'U';
char jobz = 'V';
DoubleMatrix *C;
Complex alpha = Complex(k_one, k_zero);
Complex  beta = Complex(k_zero, k_zero);
Complex gamma = Complex(-k_one, k_zero);
ComplexMatrix *Residual_tmp, *xtmp_salc, *xtmp_salc1, *C_element;

  AllocateDoubleMatrix(&C,&nmix1,&nmix1,job);
  AllocateDoubleArray(&B, &nmix1, job);
  work = (int *) malloc(nmix1 * sizeof(int));
  ResetDoubleArray(B,&nmix1);
  ResetDoubleMatrix(C);
  B[nmix] = -k_one;
  for (i = 0; i < nmix; i++) {
    C->a[i][nmix] = -k_one;
    C->a[nmix][i] = -k_one;
    work[i] = 0;
   }
    C->a[nmix][nmix] = k_zero;

  for (iii = 0; iii < symmetry->number_of_classes; iii++) {
    dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
    if (dim_salc == 0) continue; 
    for (kk = job->mixing_order - 1; kk > 0; kk--) {
      for (ii = 0; ii < dim_salc; ii++) {
	for (jj = 0; jj < dim_salc; jj++) {
	  FockHist[iii +  kk * symmetry->number_of_classes]->a[ii][jj] = FockHist[iii + (kk - 1) * symmetry->number_of_classes]->a[ii][jj];
	 }
	}
       }
    for (ii = 0; ii < dim_salc; ii++) {
      for (jj = 0; jj < dim_salc; jj++) {
	FockHist[iii]->a[ii][jj] = F_k0_salc[iii]->a[ii][jj];
       }
      }
    for (kk = job->mixing_order - 1; kk > 0; kk--) {
      for (ii = 0; ii < dim_salc2[iii]; ii++) {
	for (jj = 0; jj < dim_salc2[iii]; jj++) {
	  Residual[iii +  kk * symmetry->number_of_classes]->a[ii][jj] = Residual[iii + (kk - 1) * symmetry->number_of_classes]->a[ii][jj];
	 }
	}
       //print_complex_matrix2(Residual[iii + kk * symmetry->number_of_classes],0,6,1.0,file);
      }

    AllocateComplexMatrix(&C_element,&dim_salc2[iii],&dim_salc2[iii],job);
    AllocateComplexMatrix(&xtmp_salc,&dim_salc2[iii],&dim_salc,job);
    AllocateComplexMatrix(&xtmp_salc1,&dim_salc,&dim_salc,job);
    AllocateComplexMatrix(&Residual_tmp,&dim_salc,&dim_salc,job);
    ResetComplexMatrix(Residual[iii]);
    ResetComplexMatrix(Residual_tmp);
    ResetComplexMatrix(xtmp_salc);
    ResetComplexMatrix(xtmp_salc1);

    ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &F_k0_salc[iii], &P0_salc[iii], &beta, &xtmp_salc1);
    ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp_salc1, &S_k_salc[iii], &beta, &Residual_tmp);
    ResetComplexMatrix(xtmp_salc1);
    ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &S_k_salc[iii], &P0_salc[iii], &beta, &xtmp_salc1);
    ComplexGEMM1(&NoTrans, &ConjTrans, &gamma, &xtmp_salc1, &F_k0_salc[iii], &alpha, &Residual_tmp);
    ResetComplexMatrix(xtmp_salc);
    ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn_salc[iii], &Residual_tmp, &beta, &xtmp_salc);
    ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp_salc, &xtrn_salc[iii], &beta, &Residual[iii]);
    DestroyComplexMatrix(&xtmp_salc,job);
    DestroyComplexMatrix(&xtmp_salc1,job);
    DestroyComplexMatrix(&Residual_tmp,job);

    for (ii = 0; ii < nmix; ii++) {
      for (jj = 0; jj < nmix; jj++) {
	ResetComplexMatrix(C_element);
        //if (job->taskid == 0) printf("ii jj %3d %3d iii %3d %3d %3d %3d %3d %3d\n",ii,jj,iii,dim_salc2[iii],\
        iii+ii * symmetry->number_of_classes,iii + jj * symmetry->number_of_classes,Residual[iii + ii * symmetry->number_of_classes]->iCols,\
        Residual[iii + jj * symmetry->number_of_classes]->iCols);
        ComplexGEMM1(&ConjTrans, &NoTrans, &alpha, &Residual[iii + ii * symmetry->number_of_classes], \
	&Residual[iii + jj * symmetry->number_of_classes], &beta, &C_element);
	//print_complex_matrix2(C_element,0,7,1.0,file);
	trace = k_zero;
        for (kk = 0; kk < dim_salc2[iii]; kk++) trace += (C_element->a[kk][kk]).real();
        C->a[ii][jj] += trace;
       }
      }
     DestroyComplexMatrix(&C_element,job);
    } // close loop on iii

    //print_double_matrix2(C,0,7,1.0,file);
    dgesv_(&nmix1,&ione,C->a[0],&nmix1,work,B,&nmix1,&info);
    mix_or_not = 1;
    for (i = 0; i < nmix; i++) { if (fabs(B[i]) > 20.0 || fabs(job->energy_change) > 0.5) mix_or_not = 0; }
    //fprintf(file.out,"mix %3d %3d %3d %14.8lf %14.8lf",mix_or_not,job->iter,job->mixing_order,B[i],job->energy_change); fprintf(file.out," B\n"); }
    if (mix_or_not == 1) {
      end = 8 - nmix;
      if (end < 0) end = 0;
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"| SIIS MIXING COEFFS");
      for (kk = 0; kk < nmix; kk++) fprintf(file.out," %8.3f",B[kk]);
      for (kk = 0; kk < end; kk++)  fprintf(file.out,"          ");
      fprintf(file.out,"|\n");
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");

    for (iii = 0; iii < symmetry->number_of_classes; iii++) {
      dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
      if (dim_salc == 0) continue; 
      for (ii = 0; ii < dim_salc; ii++) {
        for (jj = 0; jj < dim_salc; jj++) {
          F_k0_salc[iii]->a[ii][jj] *= B[0];
         }
        }
      for (kk = 1; kk < nmix; kk++) {
        for (ii = 0; ii < dim_salc; ii++) {
          for (jj = 0; jj < dim_salc; jj++) {
            F_k0_salc[iii]->a[ii][jj] += B[kk] * FockHist[iii + kk * symmetry->number_of_classes]->a[ii][jj];
           }
          }
         }
        }
       }

    free(work);
    DestroyDoubleArray(&B, &nmix1, job);
    DestroyDoubleMatrix(&C,job);

}

void C_DIIS_extrapolation(ComplexMatrix **F_k0, ComplexMatrix **P0, ComplexMatrix **S_k, ComplexMatrix **xtrn, ComplexMatrix **FockHist, ComplexMatrix **Residual, int *dim1, int* dim11, JOB_PARAM *job, FILES file) 

{

int i, j, ii, jj, kk;
int info = 0, ione = 1;
int nmix = (job->mixing_order < job->iter) ? job->mixing_order : job->iter;
int nmix1 = nmix + 1;
int end, mix_or_not;
int *work;
double *B;
double trace;
char NoTrans = 'N', ConjTrans = 'C';
char uplo = 'U';
char jobz = 'V';
DoubleMatrix *C;
ComplexMatrix *xtmp, *xtmp1, *C_element, *Residual_tmp;
Complex alpha = Complex(k_one, k_zero);
Complex  beta = Complex(k_zero, k_zero);
Complex gamma = Complex(-k_one, k_zero);

  for (kk = job->mixing_order - 1; kk > 0; kk--) {
    for (ii = 0; ii < *dim1; ii++) {
      for (jj = 0; jj < *dim1; jj++) {
        FockHist[kk]->a[ii][jj] = FockHist[kk - 1]->a[ii][jj];
       }
      }
     }

  for (ii = 0; ii < *dim1; ii++) {
    for (jj = 0; jj < *dim1; jj++) {
      FockHist[0]->a[ii][jj] = F_k0[0]->a[ii][jj];
     }
    }

  for (kk = job->mixing_order - 1; kk > 0; kk--) {
    for (ii = 0; ii < *dim11; ii++) {
      for (jj = 0; jj < *dim11; jj++) {
        Residual[kk]->a[ii][jj] = Residual[kk - 1]->a[ii][jj];
       }
      }
     }

    AllocateComplexMatrix(&xtmp,dim11,dim1,job);
    AllocateComplexMatrix(&xtmp1,dim1,dim1,job);
    AllocateComplexMatrix(&Residual_tmp,dim1,dim1,job);
    ResetComplexMatrix(Residual[0]);
    ResetComplexMatrix(Residual_tmp);
    ResetComplexMatrix(xtmp);
    ResetComplexMatrix(xtmp1);
    ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &F_k0[0],&P0[0], &beta, &xtmp1);
    ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1, &S_k[0], &beta, &Residual_tmp);
    ResetComplexMatrix(xtmp1);
    ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &S_k[0], &P0[0], &beta, &xtmp1);
    ComplexGEMM1(&NoTrans, &ConjTrans, &gamma, &xtmp1,&F_k0[0], &alpha, &Residual_tmp);
    ResetComplexMatrix(xtmp);
    ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn[0], &Residual_tmp, &beta, &xtmp);
    ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp, &xtrn[0], &beta, &Residual[0]);
    DestroyComplexMatrix(&xtmp,job);
    DestroyComplexMatrix(&xtmp1,job);

    AllocateComplexMatrix(&C_element,dim11,dim11,job);
    AllocateDoubleMatrix(&C,&nmix1,&nmix1,job);
    AllocateDoubleArray(&B, &nmix1, job);
    work = (int *) malloc(nmix1 * sizeof(int));
    ResetDoubleArray(B,&nmix1);
    B[nmix] = -k_one;
    for (i = 0; i < nmix; i++) {
      C->a[i][nmix] = -k_one;
      C->a[nmix][i] = -k_one;
      work[i] = 0;
     }
      C->a[nmix][nmix] = k_zero;

    for (ii = 0; ii < nmix; ii++) {
      for (jj = 0; jj < nmix; jj++) {
	ResetComplexMatrix(C_element);
        ComplexGEMM1(&ConjTrans, &NoTrans, &alpha, &Residual[ii], &Residual[jj], &beta, &C_element);
	//print_complex_matrix2(C_element,0,6,1.0,file);
	trace = k_zero;
        for (kk = 0; kk < *dim11; kk++) trace += (C_element->a[kk][kk]).real();
	  C->a[ii][jj] = trace;
	 }
	}
       //print_double_matrix2(C,0,7,1.0,file);

    dgesv_(&nmix1,&ione,C->a[0],&nmix1,work,B,&nmix1,&info);
    mix_or_not = 1;
    for (i = 0; i < nmix; i++) { if (fabs(B[i]) > 20.0 || fabs(job->energy_change) > 0.5) mix_or_not = 0; }
    if (mix_or_not == 1) {
      end = 8 - nmix;
      if (end < 0) end = 0;
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"| DIIS MIXING COEFFS");
      for (kk = 0; kk < nmix; kk++) fprintf(file.out," %8.3f",B[kk]);
      for (kk = 0; kk < end; kk++)  fprintf(file.out,"          ");
      fprintf(file.out,"|\n");
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      for (ii = 0; ii < *dim1; ii++) {
        for (jj = 0; jj < *dim1; jj++) {
          F_k0[0]->a[ii][jj] *= B[0];
         }
        }
      for (kk = 1; kk < nmix; kk++) {
        for (ii = 0; ii < *dim1; ii++) {
          for (jj = 0; jj < *dim1; jj++) {
            F_k0[0]->a[ii][jj] += B[kk] * FockHist[kk]->a[ii][jj];
           }
          }
         }
        }

    free(work);
    DestroyDoubleArray(&B, &nmix1, job);
    DestroyDoubleMatrix(&C,job);
    DestroyComplexMatrix(&C_element,job);

}
