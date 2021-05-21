/*
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <mpi.h>
//#include <xc.h>
#include <unistd.h>
#include <stdlib.h>
*/
#include <cstring>
#include "mycomplex.h"
#include "myconstants.h"
#include "conversion_factors.h"
#include "USER_DATA.h"
#include "LIMITS.h"
#include "ALLOCATE_MEMORY.h"
#include "PARALLEL.h"
#include "PAIRS_QUADS.h"
#include "KPOINTS.h"
#include "ROTATIONS_MOLECULE.h"
#include "FOURIER_TRANSFORM.h"
#include "SYMMETRY_ADAPTATION.h"
#include "MATRIX_UTIL.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "PRINT_UTIL.h"
#include "INTEGRALS_2C_MOLECULE.h"
#include "DENSITY_MATRIX_MOLECULE.h"
#include "DENSITY_MATRIX_CRYSTAL.h"
#include "BUILD_FOCK_MATRIX.h"
//#include "SCF_MOLECULE.h"
#include "SCF_CRYSTAL.h"

using namespace std;

  // ******************************************************************************************
  // * This routine is the SCF loop                                                           *
  // ******************************************************************************************

void scf_crystal(FERMI *fermi, ATOM *atoms, ATOM *atoms_ax, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  int begin_k, end_k;
  int i, j, k, l, p, s;
  int dim1, dim2, dim3;
  int dimp_spin, dimf_spin, dimp_read, dimf_read, dimp_read_spin, dimf_read_spin;
  int dim, dimf, dimg, dimp, dim11;
  int count;
  int nk[2], ksize;
  int num_sym = symmetry->number_of_operators;
  int vector_size, value_size;
  int mixing_order = job->mixing_order;
  int ngrid_points;
  int restart_parameters[6];
  int Function[8];
  int info = 0;
  int *irrep, *p_irrep;
  double nuc_nuc;
  double total_electrons, total_population, fac;
  double *eigenvalues;
  double time1, time2, time3, time4, time5, time6, time7;
  double homo, lumo, gap;
  double *Fock, *eigval, *eigval1, *p_eigval;
  double *P0, *P1, *F0, *F1;
  double *P2;
  double *delta_P0, *delta_F0, *lse_F0;
  double threshold = 1e-02;
  //double threshold = 1e-05;
  char NoTrans = 'N', ConjTrans = 'C';
  char uplo = 'U';
  char jobz = 'V';
  char buf1[110], buf2[110], buf4[110];
  char xx[10] = "/evalfile", yy[10] = "/scf_evec", zz[20] = "/new_density_matrix";
  char filename[15] = "/scf_evec";
  Complex alpha = Complex(k_one, k_zero);
  Complex beta = Complex(k_zero, k_zero);
  PAIR_TRAN pair_p;
  KPOINT_TRAN knet;
  INT_1E one_ints, one_ints_buffer, one_ints_buffer1;
  DoubleMatrix *C_i, *C_o, *D;
  ComplexMatrix *Residual[job->mixing_order * symmetry->number_of_classes];
  ComplexMatrix *FockHist[job->mixing_order * symmetry->number_of_classes];
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

  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p,end_p,receive_p,offset_p,&pair_p,atoms,job,file);
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
    knet_size(&ksize,fermi->is,crystal);
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
    count_k_points(&knet,fermi->is,crystal,symmetry,job,file);
    fermi->knet = &knet;
    fermi->nktot = knet.nktot;
    if (job->kss == 0)      fermi->nkunique = knet.nktot;
    else if (job->kss == 1) fermi->nkunique = knet.unique;
    allocate_k_points(&knet,crystal,job,file);
    generate_k_points(&knet,fermi->is,crystal,symmetry,job,file);
    print_knet(&knet, fermi->is, crystal, job, file);
    allocate_fermi(fermi,atoms,job,file);
   }
    if (job->kpoints == 1) {
    knet.unique = fermi->npoints;
    knet.nktot  = fermi->npoints;
    allocate_k_points(&knet,crystal,job,file);
    if (job->taskid == 0 && job->verbosity > 1) fprintf(file.out,"K points read from density matrix\n");
    for (i = 0; i < fermi->npoints; i++) {
      knet.cart[i].comp1 = fermi->knet_list[i].cart.comp1;
      knet.cart[i].comp2 = fermi->knet_list[i].cart.comp2;
      knet.cart[i].comp3 = fermi->knet_list[i].cart.comp3;
      knet.ibz[i] = i;
      if (job->taskid == 0 && job->verbosity > 1) {
        fprintf(file.out,"%3d %3d %10.4lf %10.4lf %10.4lf\n",i,fermi->npoints,knet.cart[i].comp1, knet.cart[i].comp2, knet.cart[i].comp3);
       }
      }
     }

  // ******************************************************************************************
  // * Allocate and initialize P and F and matrices guess_type == 0 new guess_type == 1 old   *
  // ******************************************************************************************

    if (job->guess_type == 0) {
    //initial_density_matrix_crystal3(P0,F0,&pair_p,fermi,atoms,atom_p,shells,gaussians,crystal,symmetry,R,R_tables,G,job,file);
    initial_density_matrix(P0,F0,&pair_p,fermi,atoms,atom_p,shells,gaussians,crystal,symmetry,R,R_tables,G,job,file);
    ResetDoubleArray(F0,&dimf_spin);
    expand_density_matrix(P0,F0,&pair_p,atoms,shells,symmetry,job,file);
   }
    if (job->guess_type == 1) {
    read_density_matrix(fermi,&P0,&dimp_read_spin,&dimf_read_spin,atoms,job,file);
    expand_density_matrix(P0,F0,&pair_p,atoms,shells,symmetry,job,file);
   }

    //if (job->guess_type == 0) 
    //initial_density_matrix_crystal3(P0,F0,&pair_p,fermi,atoms,atom_p,shells,gaussians,crystal,symmetry,R,R_tables,G,job,file);
    //if (job->guess_type == 1) {
    //read_density_matrix(fermi,&P0,&dimp_read_spin,&dimf_read_spin,atoms,job,file);
    //expand_density_matrix(P0,F0,&pair_p,atoms,shells,symmetry,job,file);
   //}

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
  // * Initialize density matrix mixing                                                       *
  // ******************************************************************************************

  if (job->guess_type == 1) job->mixing_order = 1;
  AllocateDoubleMatrix(&C_i,&job->mixing_order,&dimp_spin,job);
  AllocateDoubleMatrix(&C_o,&job->mixing_order,&dimp_spin,job);
  AllocateDoubleMatrix(&D,&job->mixing_order,&dimp_spin,job);
  int iii, jjj, count1, dim_salc;
  if (job->sym_adapt == 1) {
  int num_irrep_in_basis[symmetry->number_of_classes];
  count_basis_irrep(num_irrep_in_basis,atoms,atom_p,shells,symmetry,job,file);
  count1 = 0;
  for (jjj = 0; jjj < job->mixing_order; jjj++) {
    for (iii = 0; iii < symmetry->number_of_classes; iii++) {
      dim_salc = symmetry->irp_dim_k[iii] * num_irrep_in_basis[iii];
      if (dim_salc == 0) dim_salc = 1;
      AllocateComplexMatrix(&Residual[count1],&dim_salc,&dim_salc,job);
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

if (job->spin_dim == 1 && job->guess_type == 0) {
ResetDoubleArray(P0,&dimp_spin);
ResetDoubleArray(F0,&dimf_spin);
ResetDoubleArray(delta_P0,&dimp_spin);
ResetDoubleArray(delta_F0,&dimf_spin);
}

  for (i = 0; i < dimp_spin; i++)
  delta_P0[i] = P0[i];
  for (i = 0; i < dimf_spin; i++)
  delta_F0[i] = F0[i];
  //print_Fock_matrix(F0, &pair_p, atoms, job, file);

  // ******************************************************************************************
  // * Nuclear repulsion energy                                                               *
  // ******************************************************************************************

  nuclear_repulsion_energy(&nuc_nuc, R, G, atoms, symmetry, crystal, job, file);

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

  //2021array_dimensions(&dim, &dimg, &pair_p, atoms, job, file); // don't change
  //2021allocate_INT_1E(&one_ints, dim, Function, job, file);
  //2021allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  //2021allocate_INT_1E(&one_ints_buffer1, dim, Function, job, file);
  allocate_INT_1E(&one_ints, dimp, Function, job, file);
  allocate_INT_1E(&one_ints_buffer, dimp, Function, job, file);
  allocate_INT_1E(&one_ints_buffer1, dimp, Function, job, file);
  time1 = MPI_Wtime();
  fock_element_1e1(&one_ints_buffer, dimp, &pair_p, pair_p.nump, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //2021fock_element_1e1(&one_ints_buffer, dim, &pair_p, pair_p.nump, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  time2 = MPI_Wtime() - time1;
  //MPI_Allreduce(&one_ints_buffer.Fock[0],&one_ints_buffer1.Fock[0],dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //MPI_Allreduce(&one_ints_buffer.Overlap[0],&one_ints_buffer1.Overlap[0],dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //MPI_Allreduce(&one_ints_buffer.Kinetic[0],&one_ints_buffer1.Kinetic[0],dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //MPI_Allreduce(&one_ints_buffer.ElecNuc[0],&one_ints_buffer1.ElecNuc[0],dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&one_ints_buffer.Fock[0],&one_ints_buffer1.Fock[0],dimp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&one_ints_buffer.Overlap[0],&one_ints_buffer1.Overlap[0],dimp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&one_ints_buffer.Kinetic[0],&one_ints_buffer1.Kinetic[0],dimp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&one_ints_buffer.ElecNuc[0],&one_ints_buffer1.ElecNuc[0],dimp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer, Function, job, file);
  MPI_Allgatherv(&one_ints_buffer1.Overlap[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Overlap[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgatherv(&one_ints_buffer1.Fock[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Fock[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgatherv(&one_ints_buffer1.Kinetic[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Kinetic[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgatherv(&one_ints_buffer1.ElecNuc[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.ElecNuc[0],receive_p,\
  offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer1, Function, job, file);

  // ******************************************************************************************
  // * Check for small eigenvalues of S_k                                                     *
  // ******************************************************************************************

  if (job->taskid == 0 && job->verbosity > 1) {
    ComplexMatrix *S_k1;
    AllocateComplexMatrix(&S_k1,&dim1,&dim1,job);
    for (k = 0; k < knet.unique; k++) {
      ResetDoubleArray(eigenvalues,&dim1);
      nk[0] = knet.ibz[k];
      nk[1] = knet.ibz[k];
      fourier_transform(&one_ints.Overlap[0], &S_k1->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      DiagonaliseHermitian(&S_k1, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);
      printf("density matrix %3d %3d %3d %9.2e\n",k,nk[0],nk[1],eigenvalues[0]);
      for (int ii = 0; ii < atoms->number_of_sh_bfns_in_unit_cell; ii++) {
        if (eigenvalues[ii] < threshold) { 
        fprintf(file.out,"small eigenvalue %d %d %e\n",k,ii,eigenvalues[ii]); 
        printf("small eigenvalue %d %d %e\n",k,ii,eigenvalues[ii]); 
       }
      }
     }
    DestroyComplexMatrix(&S_k1,job);
    ResetDoubleArray(eigenvalues,&dim1);
   }

// HEREdouble *Fock_EN;
//AllocateDoubleArray(&Fock_EN,&dimp_spin,job);
//ResetDoubleArray(Fock_EN,&dimp_spin);
//fock_element_elecnuc(Fock_EN,&pair_p,R,G,atoms,shells,gaussians,crystal,job,file);
//DestroyDoubleArray(&Fock_EN,&dimp_spin,job);
//MPI_Finalize();
//exit(0);

/*
  ResetDoubleArray(Fock, &dimp_spin);
  count = 0;
  for (s = 0; s < job->spin_dim; s++) {
    for (i = 0; i < dimp; i++) {
      //Fock[count] += one_ints.Kinetic[i];
      Fock[count] += one_ints.Fock[i];
      count++;
     }
    }

    if (job->taskid == 0 && job->verbosity > 1) {
    fprintf(file.out,"Fock matrix 1e\n");
    count = 0;
    for (s = 0; s < job->spin_dim; s++) {
     for (p = 0; p < pair_p.nump; p++) {
      int q = pair_p.posn[p];
      fprintf(file.out,"spin %3d pair %3d posn %3d\n",s,p,q);
      for (i = 0; i < atoms->bfnnumb_sh[pair_p.cell1[q]]; i++) {
	for (j = 0; j < atoms->bfnnumb_sh[pair_p.cell2[q]]; j++) {
	  fprintf(file.out,"%10.4lf ",Fock[count]);
	  count++;
	 }
	fprintf(file.out,"\n");
       }
      fprintf(file.out,"\n\n");
     }
    fprintf(file.out,"\n\n");
   }
  }
*/

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

    double *S1, *S2;
    AllocateDoubleArray(&S1,&job->dimf,job);
    AllocateDoubleArray(&S2,&job->dimf,job);
    if (crystal->type[0] != 'M') AllocateDoubleArray(&S2,&job->dimf,job);
    if (job->scf_coulomb  == 1 || job->scf_exchange == 1)
    shell_screen_compute_integrals(S1,S2,&pair_p,R,G,atoms,shells,gaussians,symmetry,crystal,job,file);

   job->guess_type = 0;

  for (job->iter = 1; job->iter <= job->max_cycle; job->iter++) {

    if (job->taskid == 0) 
    printf("ITERATION %d\n",job->iter);

//if (job->iter > 15) job->scf_exchange = 1;
//if (job->taskid == 0) printf("job->scf_exchange = 1\n");

    ResetIntArray(irrep,&dim3);
    ResetDoubleArray(eigval,&dim3);

    time6 = MPI_Wtime();
    time7 = k_zero;

  // ******************************************************************************************
  // * Calculate 2e contribution to Fock matrix                                               *
  // ******************************************************************************************

   time2 = MPI_Wtime();
   build_fock_matrix(Fock,&one_ints,fermi,&total_electrons,S1,S2,P0,F0,delta_P0,delta_F0,&pair_p,atom_p,atoms, \
   shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
   //build_fock_matrix(Fock,&dft_grid,&one_ints,fermi,&total_electrons,S1,S2,P0,F0,delta_P0,delta_F0,&pair_p,atom_p,atoms, \
   shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
   //build_fock_matrix_no_sym(Fock,&dft_grid,&one_ints,fermi,&total_electrons,S1,S2,P0,F0,delta_P0,delta_F0,&pair_p,atom_p,atoms, \
   shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
   if (job->taskid == 0) printf("\n");
   time1 = MPI_Wtime() - time2;

   //if (crystal->type[0] != 'M') 
   //scf_evectors = fopen(zz1, "wb");
   //else if (crystal->type[0] == 'M' && job->taskid == 0) {
   //scf_evectors = fopen(zz2, "wb");
   //scf_evalues = fopen(zz4, "wb");
  //}

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
   double total_population = print_total_population2(atoms, shells, job, file);
   //double fac = normalise_density_matrix2(P0, F0, total_population, job, file);
   fac = k_one;
   total_energy(&one_ints,Fock,P0,&pair_p,atoms,job,file);
   if (job->taskid >= 0) {
   //printf("Density matrix normalisation factor %16.10e\n",fac - k_one);
   fprintf(file.out,"===========================================================================================================\n");
   fprintf(file.out,"| ITERATION             %3d | ENERGY %16.9e | ENERGY CHANGE %9.2e | DENSITY NORM  %9.6lf |\n", \
   job->iter,job->total_energy,job->energy_change,fac);
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

  int *counter, winRank, myCounter, increment;
  const int izero = 0;
  MPI_Win win;
  winRank = 0;
  CreateCounter(job->taskid, winRank, &counter, izero, &win, MPI_COMM_WORLD);
  myCounter = 0;

  if ((crystal->type[0] == 'M'  && job->taskid == 0) || (crystal->type[0] != 'M' && (job->taskid > 0 || job->numtasks == 1))) {
  while (myCounter < knet.unique) {

    increment = (knet.unique - myCounter) / job->numtasks / 2;
    if (increment <  1) increment =  1;
    if (increment > 64) increment = 64;
    increment = 1;
    myCounter = GetCounter(winRank, increment , &win);
    begin_k = (myCounter             < knet.unique) ? myCounter : knet.unique;
    end_k   = (myCounter + increment < knet.unique) ? myCounter + increment : knet.unique;
    //printf("proc %d myCounter %d\n",job->taskid,myCounter);
    //printf("proc %d Fock begin %d end %d out of %d\n",job->taskid,begin_k,end_k,fermi->nkunique);

    for (k = begin_k; k < end_k; k++) {

      nk[0] = knet.ibz[k];
      nk[1] = knet.ibz[k];

  // ******************************************************************************************
  // * Work in non-symmetry adapted basis                                                     *
  // ******************************************************************************************

    if (job->sym_adapt == 0) {

    ComplexMatrix *Cholesky_S_k_inverse;
    int *Cholesky_S_k_pivot;
    ComplexMatrix *xtrn1, *F_k0, *F_k1, *S_k, *xtmp1, *eigenvectors1;
    AllocateComplexMatrix(&S_k,&dim1,&dim1,job);

    job->scf_trans = 0; // set to canonical transformation = 1 is Cholesky

    fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], &knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);

  // ******************************************************************************************
  // * Canonical transformation to orthogonalise basis at this k point and store in xtrn1     *
  // ******************************************************************************************

    if (job->scf_trans == 0) { // canonical transformation

    DiagonaliseHermitian(&S_k, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);

    int n = 0;
    for (int ii = 0; ii < atoms->number_of_sh_bfns_in_unit_cell; ii++) {
      if (job->iter == 1 && eigenvalues[ii] < threshold) { 
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

  // ******************************************************************************************
  // * Cholesky transformation to orthogonalise basis at this k point and store in xtrn1      *
  // ******************************************************************************************

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

  // ******************************************************************************************
  // * Canonical transformation of Fock matrix                                                *
  // ******************************************************************************************

    if (job->scf_trans == 0) {
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

  // ******************************************************************************************
  // * Cholesky transformation of Fock matrix                                                 *
  // ******************************************************************************************

    else if (job->scf_trans == 1) {
      //   dim11 = 0;
      //   for (i = 0; i < dim1; i++) {
      //     if ((S_k->a[i][i]).real() > threshold) dim11++;
      //    }
      ResetComplexMatrix(xtmp1);
      CholeskyPermuteHermitian(&F_k0,Cholesky_S_k_pivot,0,job);
      ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &Cholesky_S_k_inverse, &F_k0, &beta, &xtmp1);
      ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1, &Cholesky_S_k_inverse, &beta, &F_k1);
     }

  // ******************************************************************************************
  // * Diagonalise Fock matrix and back-transform by Canonical or Cholesky                    *
  // ******************************************************************************************

    time2 = MPI_Wtime();
    ResetComplexMatrix(eigenvectors1);
    p_eigval = &eigval[(s * knet.unique + k) * dim1];
    DiagonaliseHermitian(&F_k1, &p_eigval, &eigenvectors1, &jobz, &uplo, &info);
    for (i = dim11; i < dim1; i++) p_eigval[i] = 9999.9; // set zero eigenvalues from S_k to large number
    //for (i = 0; i < dim1; i++) fprintf(file.out,"%3d %10.4lf\n",i,p_eigval[i]);

    if (job->scf_trans == 0) {
      ResetComplexMatrix(eigvec);
      ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigenvectors1, &xtrn1, &beta, &eigvec);
     }

    else if (job->scf_trans == 1) {
      ResetComplexMatrix(eigvec);
      ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigenvectors1, &Cholesky_S_k_inverse, &beta, &eigvec);
      CholeskyPermuteHermitian(&eigvec,Cholesky_S_k_pivot,2,job);
     }

    if (job->verbosity > 1) {
      fprintf(file.out,"Fock eigenvectors1 spin  %3d k-point %3d\n",s,k);
      ////print_complex_eigenvector_matrix(eigenvectors1, eigenvalues, file);
      fprintf(file.out,"Back-transformed Fock eigec spin  %3d k-point %3d\n",s,k);
      ////print_complex_eigenvector_matrix(eigvec, eigenvalues, file);
     }

    if (job->vectors == 2 && job->mpp == 0) {
      MPI_File_seek(fh, (s * knet.unique + k) * vector_size * sizeof(Complex), MPI_SEEK_SET) ;
      MPI_File_write(fh, &eigvec->a[fermi->bands[0] - 1][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
      //if (crystal->type[0] == 'M' && job->taskid == 0) {
      //fwrite(&eigvec->a[fermi->bands[0] - 1][0], sizeof(double), 2 * vector_size, scf_evectors);
      //fflush(scf_evectors);
     //}
    }

    else if (job->mpp == 1) {
    time2 = MPI_Wtime();
    reduced_density_matrix_crystal3(fermi,P2,&knet,s,k,eigvec,R,R_tables,&pair_p,atom_p,atoms,shells,symmetry,job,file);
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

//    else if (job->sym_adapt == 2 || crystal->type[0] == 'R') {
//
//    ATOM_TRAN atom_k;
//    SYMMETRY symmetry_little_k_group;
//    REAL_LATTICE_TABLES R_tables_little_k;
//    count_little_k_group_operators(k,&symmetry_little_k_group,symmetry,crystal,&knet,fermi,job,file);
//    //CHP March2019count_little_k_group_operators(k,&symmetry_little_k_group,symmetry,&knet,fermi,job,file);
//    allocate_SYMMETRY(&symmetry_little_k_group,job,file);
//    generate_little_k_group(k, &symmetry_little_k_group, fermi, &knet, symmetry, crystal, job, file);
//    //march 2019 generate_little_k_group(k, &symmetry_little_k_group, fermi, &knet, symmetry, job, file);
//
//    allocate_ATOM_TRAN(&atom_k,atoms,&symmetry_little_k_group,job,file);
//    generate_ATOM_TRAN(crystal,R,&symmetry_little_k_group,atoms,&atom_k,job,file);
//
//    R_tables_little_k.last_vector = R->last_vector;
//    R_tables_little_k.margin_vector = R->margin_vector;
//    R_tables_little_k.max_vector = R->max_vector;
//    allocate_REAL_LATTICE_TABLES(&R_tables_little_k,&symmetry_little_k_group,job,file);
//    generate_real_lattice_tables(crystal,&symmetry_little_k_group,R,&R_tables_little_k,job,file);
//
//    PAIR_TRAN pair_k;
//    fprintf(file.out,"LITTLE GROUP OF K POINT %3d\n",k);
//    printf("LITTLE GROUP OF K POINT %3d %3d\n",job->taskid,k);
//    count_density_pairs4(&pair_k,atoms,&atom_k,&symmetry_little_k_group,R,&R_tables_little_k,job,file);
//    allocate_PAIR_TRAN(&pair_k,atoms,&symmetry_little_k_group,&R_tables_little_k,job,file);
//    generate_density_pairs4(&pair_k,atoms,&atom_k,&symmetry_little_k_group,R,&R_tables_little_k,job,file);
//    if (job->taskid == 0 && job->verbosity > 1) 
 //   print_pairs(&pair_k,atoms,R,job,file);
 //
 //   int iii, jjj, n1, dim_salc, dim_salc1, num_irrep_in_basis[symmetry_little_k_group.number_of_classes];
 //   int num_eigval, num_eigval1;
 //   int  *Cholesky_S_k_pivot[symmetry_little_k_group.number_of_classes];
 //   double  *eigenvalues_S_k[symmetry_little_k_group.number_of_classes];
 //   double *eigenvalues_salc[symmetry_little_k_group.number_of_classes];
 //   double threshold, eigvec_bas;
 //
 //   threshold = 0.0001;
 //   job->scf_trans = 0;
 //
 //   if (crystal->type[0] == 'C')
 //   count_basis_irrep_crystal(num_irrep_in_basis,atoms,&atom_k,&R_tables_little_k,shells,&symmetry_little_k_group,job,file);
 //
//    // FIX LATT IN ABOVE
//    if (crystal->type[0] == 'M')
//    count_basis_irrep(num_irrep_in_basis,atoms,&atom_k,shells,&symmetry_little_k_group,job,file);
//    ComplexMatrix *Cholesky_S_k[symmetry_little_k_group.number_of_classes];
//   ComplexMatrix *Cholesky_S_k_inverse[symmetry_little_k_group.number_of_classes];
//   ComplexMatrix *eigenvectors_S_k[symmetry_little_k_group.number_of_classes];
//   ComplexMatrix *F_k0_salc[symmetry_little_k_group.number_of_classes];
//   ComplexMatrix *F_k1_salc[symmetry_little_k_group.number_of_classes];
 //   ComplexMatrix *S_salc1[symmetry_little_k_group.number_of_classes];
 //   ComplexMatrix *xtmp1_salc[symmetry_little_k_group.number_of_classes];
 //   ComplexMatrix *eigenvectors_salc[symmetry_little_k_group.number_of_classes];
 //   ComplexMatrix *eigvec_salc[symmetry_little_k_group.number_of_classes];
 //   ComplexMatrix *xtrn1_salc[symmetry_little_k_group.number_of_classes];
 //
 // // ******************************************************************************************
 // // * Allocate memory                                                                        *
 // // ******************************************************************************************
 // 
 //    for (iii = 0; iii < symmetry_little_k_group.number_of_classes; iii++) {
 //    dim_salc = symmetry_little_k_group.irp_dim_k[iii] * num_irrep_in_basis[iii];
 //    //printf("im here %3d %3d %3d\n",job->taskid,iii,dim_salc);// }
 //    if (dim_salc == 0) dim_salc = 1;
 //    AllocateComplexMatrix(&S_salc1[iii],&dim_salc,&dim_salc,job);
 //    AllocateComplexMatrix(&Cholesky_S_k[iii],&dim_salc,&dim_salc,job);
 //   AllocateComplexMatrix(&eigenvectors_S_k[iii],&dim_salc,&dim_salc,job);
 //   AllocateComplexMatrix(&eigvec_salc[iii],&dim_salc,&dim_salc,job);
 //   AllocateComplexMatrix(&F_k0_salc[iii],&dim_salc,&dim_salc,job);
 //   AllocateDoubleArray(&eigenvalues_salc[iii],&dim_salc,job);
 //   AllocateDoubleArray(&eigenvalues_S_k[iii],&dim_salc,job);
 //  }
 //
 //   if (job->scf_trans == 0) 
 //   fourier_transform_salc_crystal(&one_ints.Overlap[0], S_salc1, &knet, nk, &atom_k, &pair_k, R, &R_tables_little_k, \
 //   atoms, shells, &symmetry_little_k_group, symmetry, job, file);
 //   //fourier_transform_salc1(&one_ints.Overlap[0], S_salc1, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job,file);
 //   //if (job->scf_trans == 1) 
 //   //fourier_transform_salc1(&one_ints.Overlap[0], Cholesky_S_k, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job,file);
 //   for (iii = 0; iii < symmetry_little_k_group.number_of_classes; iii++) {
 //     dim_salc = symmetry_little_k_group.irp_dim_k[iii] * num_irrep_in_basis[iii];
 //     if (dim_salc == 0) dim_salc = 1;
 //     if (num_irrep_in_basis[iii] == 0) continue;
 //
 //     if (job->scf_trans == 0) { // canonical transformation
 //       DiagonaliseHermitian(&S_salc1[iii], &eigenvalues_S_k[iii], &eigenvectors_S_k[iii], &jobz, &uplo, &info);
 //       n1 = 0;
 //       for (jjj = 0; jjj < dim_salc; jjj++) {
 //         //if (eigenvalues_S_k[iii][jjj] < threshold && job->taskid == 0)  
 //         fprintf(file.out,"small eigenvalue %3d irrep %3d eigenvalue %3d %e\n",k,iii,jjj,eigenvalues_S_k[iii][jjj]); 
 //         if (eigenvalues_S_k[iii][jjj] < threshold)  n1++; 
 //        }
 //         dim_salc1 = dim_salc - n1;
 //         AllocateComplexMatrix(&xtrn1_salc[iii],&dim_salc1,&dim_salc,job);
 //         AllocateComplexMatrix(&xtmp1_salc[iii],&dim_salc1,&dim_salc,job);
 //         AllocateComplexMatrix(&F_k1_salc[iii],&dim_salc1,&dim_salc1,job);
 //         AllocateComplexMatrix(&eigenvectors_salc[iii],&dim_salc1,&dim_salc1,job);
 //         for (i = n1; i < dim_salc; i++) { 
 //   for (j = 0; j < dim_salc; j++) { 
 //             xtrn1_salc[iii]->a[i - n1][j] = eigenvectors_S_k[iii]->a[i][j] / sqrt(*(eigenvalues_S_k[iii] + i));
 //            }
 //           }
 //          }
 //
 //     else if (job->scf_trans == 1) { // cholesky
 //       AllocateIntArray(&Cholesky_S_k_pivot[iii],&dim_salc,job);
 //       CholeskyHermitian(&Cholesky_S_k[iii],Cholesky_S_k_pivot[iii],&dim_salc1,&threshold);
 //       AllocateComplexMatrix(&xtrn1_salc[iii],&dim_salc1,&dim_salc,job);
 //       AllocateComplexMatrix(&Cholesky_S_k_inverse[iii],&dim_salc1,&dim_salc,job);
 //       CholeskyInverseHermitian(&Cholesky_S_k[iii],&Cholesky_S_k_inverse[iii],&dim_salc1);
 //       AllocateComplexMatrix(&xtmp1_salc[iii],&dim_salc1,&dim_salc,job);
 //       AllocateComplexMatrix(&F_k1_salc[iii],&dim_salc1,&dim_salc1,job);
 //       AllocateComplexMatrix(&eigenvectors_salc[iii],&dim_salc1,&dim_salc1,job);
 //      }
 //
 //     } // close loop on iii

  // ******************************************************************************************
  // * Loop over spin if spin-polarized                                                       *
  // ******************************************************************************************

//  for (s = 0; s < job->spin_dim; s++) {

    //fourier_transform_salc1(&Fock[s*dimp], F_k0_salc, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job, file);
    //fourier_transform_salc_crystal(&Fock[s*dimp], F_k0_salc, &knet, nk, &atom_k, &pair_k, R, &R_tables_little_k, \
    atoms, shells, &symmetry_little_k_group, symmetry, job, file);

    //ComplexMatrix *P0_salc[symmetry_little_k_group.number_of_classes];
    //ComplexMatrix *S0_salc[symmetry_little_k_group.number_of_classes];
    //for (iii = 0; iii < symmetry_little_k_group.number_of_classes; iii++) {
      //dim_salc = symmetry_little_k_group.irp_dim_k[iii] * num_irrep_in_basis[iii];
      //if (dim_salc == 0) dim_salc = 1;
      //AllocateComplexMatrix(&P0_salc[iii],&dim_salc,&dim_salc,job);
      //AllocateComplexMatrix(&S0_salc[iii],&dim_salc,&dim_salc,job);
      //ResetComplexMatrix(P0_salc[iii]);
      //ResetComplexMatrix(P0_salc[iii]);
     //}
      //job->diis = 0;
      //if (job->diis == 1) {
      //fourier_transform_salc1(&P0[s*dimp], P0_salc, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job, file);
      //fourier_transform_salc1(&one_ints.Overlap[0], S0_salc, &knet, nk, atom_p, &pair_p, R, atoms, shells, symmetry, job,file);
     //} 
     // NEED TO FREE P0_salc and S0_salc somewhere below

 //     num_eigval = 0, num_eigval1 = 0;
 //
 //     for (iii = 0; iii < symmetry_little_k_group.number_of_classes; iii++) {
 //       if (num_irrep_in_basis[iii] == 0) continue;
 //       dim_salc = symmetry_little_k_group.irp_dim_k[iii] * num_irrep_in_basis[iii];
 //       ResetComplexMatrix(xtmp1_salc[iii]);
 //       ResetComplexMatrix(F_k1_salc[iii]);
 //
 //       if (job->scf_trans == 0) {
 //         n1 = 0;
 //         //for (jjj = 0; jjj < dim_salc; jjj++) {
 //         //if (eigenvalues_S_k[iii][jjj] < threshold)  n1++; 
 //        //}
 //         dim_salc1 = dim_salc - n1;
 //         ResetComplexMatrix(xtmp1_salc[iii]);
 //         ResetComplexMatrix(F_k1_salc[iii]); // CHANGES2017
 //         ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &xtrn1_salc[iii], &F_k0_salc[iii], &beta, &xtmp1_salc[iii]);
 //         ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1_salc[iii], &xtrn1_salc[iii], &beta, &F_k1_salc[iii]);
 //         //if (job->diis == 1)
 //           //C_DIIS_extrapolation(F_k1_salc, S0_salc, P0_salc, FockHist, Residual, num_irrep_in_basis, &iii, symmetry, job, file);
 //        }
 //
 //       else if (job->scf_trans == 1) {
 //         dim_salc1 = 0;
 //         for (i = 0; i < dim_salc; i++) {
 //         if ((Cholesky_S_k[iii]->a[i][i]).real() > threshold) dim_salc1++;
 //        }
 //         CholeskyPermuteHermitian(&F_k0_salc[iii],Cholesky_S_k_pivot[iii],0,job);
 //         ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &Cholesky_S_k_inverse[iii], &F_k0_salc[iii], &beta, &xtmp1_salc[iii]);
 //         ComplexGEMM1(&NoTrans, &ConjTrans, &alpha, &xtmp1_salc[iii], &Cholesky_S_k_inverse[iii], &beta, &F_k1_salc[iii]);
 //        }
 //
 //       //DiagonaliseHermitian(&F_k1_salc[iii], &eigenvalues_salc[iii], &eigenvectors_salc[iii], &jobz, &uplo, &info);
 //
 //       if (job->scf_trans == 0) {
 //         ResetComplexMatrix(eigvec_salc[iii]);
 //         ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigenvectors_salc[iii], &xtrn1_salc[iii], &beta, &eigvec_salc[iii]);
 //        }
 //
 //       else if (job->scf_trans == 1) {
 //         ResetComplexMatrix(eigvec_salc[iii]);
 //         ComplexGEMM1(&NoTrans, &NoTrans, &alpha, &eigenvectors_salc[iii], &Cholesky_S_k_inverse[iii], &beta, &eigvec_salc[iii]);
 //         //print_complex_eigenvector_matrix2(eigvec_salc[iii], eigenvalues_salc[iii],6,12,au_to_eV,file);
 //         CholeskyPermuteHermitian(&eigvec_salc[iii],Cholesky_S_k_pivot[iii],2,job);
 //        }
 //
 //       //if (symmetry_little_k_group.irp_dim_k[iii] == 2) 
 //         //transform_degenerate_salc(eigvec_salc, eigenvalues_salc, &iii, num_irrep_in_basis, &knet, nk, symmetry, job, file);
 //         for (jjj = 0; jjj < dim_salc; jjj++) irrep[num_eigval + jjj] = iii;
 //         for (jjj = 0; jjj < dim_salc1; jjj++) eigval[num_eigval + jjj] = eigenvalues_salc[iii][jjj];
 //         for (i = 0; i < dim_salc1; i++) { 
 //         eigvec_bas = k_zero;
 //           for (j = 0; j < dim_salc1; j++) { 
 //             eigvec_bas += fabs((eigvec_salc[iii]->a[i][j]).real());
 //            }
 //           }
 //        for (jjj = dim_salc1; jjj < dim_salc; jjj++) eigval[num_eigval + jjj] = 9999.0;
 //        num_eigval  += dim_salc;
 //        num_eigval1 += dim_salc1;
 //
 //  } // close loop on iii
 //
 //   if (job->kpoints == 0 && num_eigval1 < fermi->occupied[k]) {
 //   if (job->taskid == 0)
 //   fprintf(file.out,"Limited basis set is also linearly dependent at k-point %3d. Stopping.\n",k);
 //   MPI_Finalize();
 //   exit(1);
 //  }
 //   ResetComplexMatrix(eigvec);

    //transform_salc_to_atomic_basis(eigvec,eigval,irrep,eigvec_salc,&knet,nk,atom_p,&pair_p,R,atoms,shells,symmetry,job,file);
    //print_complex_eigenvector_matrix2(eigvec, eigval, 6, num_eigval, 1.0, file);

    //if (job->vectors == 2 && job->mpp == 0) {
    //MPI_File_seek(fh, (s * knet.unique + k) * vector_size * sizeof(Complex), MPI_SEEK_SET) ;
    //MPI_File_write(fh, &eigvec->a[fermi->bands[0] - 1][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    //if (crystal->type[0] == 'M' && job->taskid == 0) {
    //fwrite(&eigvec->a[fermi->bands[0] - 1][0], sizeof(double), 2 * vector_size, scf_evectors);
    //fflush(scf_evectors);
   //}
  //}
    //else if (job->mpp == 1) {
    //time2 = MPI_Wtime();
    //reduced_density_matrix_crystal3(fermi,P2,&knet,s,k,eigvec,R,R_tables,&pair_p,atom_p,atoms,shells,symmetry,job,file);
    //time7 += MPI_Wtime() - time2;
   //}

// } // close loop on s

 //  for (iii = 0; iii < symmetry_little_k_group.number_of_classes; iii++) {
 //    if (num_irrep_in_basis[iii] == 0) continue;
 //     dim_salc = symmetry_little_k_group.irp_dim_k[iii] * num_irrep_in_basis[iii];
 //     if (job->scf_trans == 1) {
 //     DestroyIntArray(&Cholesky_S_k_pivot[iii],&dim_salc,job);
 //     DestroyComplexMatrix(&Cholesky_S_k_inverse[iii],job);
 //    }
 //     DestroyDoubleArray(&eigenvalues_salc[iii],&dim_salc,job);
 //     DestroyDoubleArray(&eigenvalues_S_k[iii],&dim_salc,job);
 //     DestroyComplexMatrix(&Cholesky_S_k[iii],job);
 //     DestroyComplexMatrix(&eigenvectors_S_k[iii],job);
 //     DestroyComplexMatrix(&xtmp1_salc[iii],job);
 //     DestroyComplexMatrix(&F_k0_salc[iii],job);
 //     DestroyComplexMatrix(&F_k1_salc[iii],job);
 //     DestroyComplexMatrix(&S_salc1[iii],job);
 //     DestroyComplexMatrix(&eigenvectors_salc[iii],job);
 //     DestroyComplexMatrix(&eigvec_salc[iii],job);
 //     DestroyComplexMatrix(&xtrn1_salc[iii],job);
 //    }
 //
 //     free_PAIR_TRAN(&pair_k,job);
 //     free_REAL_LATTICE_TABLES(&R_tables_little_k,job);
 //     free_ATOM_TRAN(&atom_k,job);
 //     free_SYMMETRY(&symmetry_little_k_group,job);

//    } // close if (job->sym_adapt == 2)

  // ******************************************************************************************
  // * End work in symmetry adapted basis                                                     *
  // ******************************************************************************************

    } // close loop over k

   } // close while (
	    
} 
      DestroyCounter(job->taskid, winRank, &win, counter);

  // ******************************************************************************************
  // * Gather eigenvalues and write to disk if job->values == 2                               *
  // ******************************************************************************************

    ResetDoubleArray(eigval1,&dim3);
    MPI_Allreduce(eigval,eigval1,dim3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if (job->values == 2) {
      if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
      MPI_File_seek(eh, 0, MPI_SEEK_SET) ;
      for (k = 0; k < knet.unique; k++) {
	for (s = 0; s < job->spin_dim; s++) {
	  MPI_File_write(eh, &eigval1[(s * knet.unique + k) * dim1 + fermi->bands[0]-1],value_size,MPI_DOUBLE,MPI_STATUS_IGNORE);
	  //if (job->taskid == 0) {
	  if (job->taskid == 0 && k == 0) {
	  fprintf(file.out,"%3d\n",k);
	  for (int eig =  0; eig < 10; eig++)
	  fprintf(file.out,"%10.4lf ", eigval1[(s * knet.unique + k) * dim1 + eig] * au_to_eV);
	  fprintf(file.out,"\n");
	  for (int eig = 10; eig < 20; eig++)
	  fprintf(file.out,"%10.4lf ", eigval1[(s * knet.unique + k) * dim1 + eig] * au_to_eV);
	  fprintf(file.out,"\n");
	  for (int eig = 20; eig < 30; eig++)
	  fprintf(file.out,"%10.4lf ", eigval1[(s * knet.unique + k) * dim1 + eig] * au_to_eV);
	  fprintf(file.out,"\n");
	  for (int eig = 30; eig < 40; eig++)
	  fprintf(file.out,"%10.4lf ", eigval1[(s * knet.unique + k) * dim1 + eig] * au_to_eV);
	  fprintf(file.out,"\n");
	  for (int eig = 40; eig < 50; eig++)
	  fprintf(file.out,"%10.4lf ", eigval1[(s * knet.unique + k) * dim1 + eig] * au_to_eV);
	  fprintf(file.out,"\n");
	  for (int eig = 50; eig < 60; eig++)
	  fprintf(file.out,"%10.4lf ", eigval1[(s * knet.unique + k) * dim1 + eig] * au_to_eV);
	  fprintf(file.out,"\n");
	 }
	 }
	}
       }
      //else if (crystal->type[0] == 'M' && job->taskid == 0) {
	//for (s = 0; s < job->spin_dim; s++) {
	//fwrite(&eigval1[s * dim1 + fermi->bands[0] - 1], sizeof(double), value_size, scf_evalues);
	//printf("%3d %3d %lf %lf %lf %lf\n",s * dim1 + fermi->bands[0] - 1,value_size,eigval1[s * dim1 + 0],eigval1[s * dim1 + 1],\
        eigval1[s * dim1 + 2],eigval1[s * dim1 + 3]);
	//fflush(scf_evalues);
       //}
      //}
     }

  // ******************************************************************************************
  // * Calculate Fermi level if system is periodic                                            *
  // ******************************************************************************************

    if (job->iter > 2 && job->fix_occ == 0 && job->vectors == 2 && (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P'))
    //if (job->fix_occ == 0 && job->vectors == 2 && (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P'))
    //calculate_fermi_level_metal_old(fermi,eigval1,fermi->bands,&knet,fermi->nkunique,fermi->nktot,atoms,job,file);

    if (job->taskid == 0 && job->vectors == 2 && crystal->type[0] == 'M') {
    printf("Fermi->occupied = %3d\n",fermi->occupied[0]);
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
     if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P')
     density_matrix_crystal2(fermi,P1,F1,&knet,filename,&fermi->nkunique,R,R_tables,&pair_p,atom_p,atoms,shells,symmetry,job,file);
     if (crystal->type[0] == 'M')
     density_matrix_molecule2(fermi,P1,F1,&knet,&fermi->nkunique,R,&pair_p,atom_p,atoms,shells,symmetry,job,file);
     time1 = MPI_Wtime() - time2;
     if (job->taskid == 0)
     printf("Time to generate density matrix %10.4e\n",(double)(time1));
    }

     else if (job->mpp == 1) {
     MPI_Reduce(P2,P1,dimp_spin,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     //printf("Time to generate MPP density matrix %10.4e\n",(double)(time7));
    }

     double *P0_copy, *F0_copy;
     AllocateDoubleArray(&P0_copy,&dimp_spin,job);
     AllocateDoubleArray(&F0_copy,&dimf_spin,job);
     for (i = 0; i < dimp_spin; i++)
     P0_copy[i] = P0[i];
     for (i = 0; i < dimf_spin; i++)
     F0_copy[i] = F0[i];

     if (job->iter == 1) { // Obtain idempotent density matrix by accepting 100% of new density matrix at first iteration
     for (i = 0; i < dimp_spin; i++)
     P0[i] = P1[i];
     for (i = 0; i < dimf_spin; i++) 
     F0[i] = F1[i];
    }
 
     //else if (job->iter > 1) {
     //for (i = 0; i < dimp_spin; i++)
     //P0[i] = job->fock_mixing * P0[i] + (k_one - job->fock_mixing) * P1[i];
     //for (i = 0; i < dimf_spin; i++) 
     //F0[i] = job->fock_mixing * F0[i] + (k_one - job->fock_mixing) * F1[i];
    //}
 
     // Pulay mixing
     //else if (fabs(job->energy_change) <= 5.0e-01) {
     else if (job->iter > 1) {
     init_mix(C_i, C_o, D, P0, P1);
     mix_density(C_i, C_o, D, P0, job, file);
    }

  //******************************************************************************************
  //* Calculate full density matrix and density matrix changes                               *
  //******************************************************************************************
 
     if (job->mpp == 1)
     MPI_Bcast(P0,dimp_spin,MPI_DOUBLE,0,MPI_COMM_WORLD);
     expand_density_matrix(P0,F0,&pair_p,atoms,shells,symmetry,job,file);
 
     //expand_density_matrix(P0,F0,&pair_p,atoms,shells,symmetry,job,file);

     for (i = 0; i < dimp_spin; i++) {
     delta_P0[i] = P0[i] - P0_copy[i];
    }

     ResetDoubleArray(delta_F0,&dimf_spin);
     expand_density_matrix(delta_P0,delta_F0,&pair_p,atoms,shells,symmetry,job,file);
  //print_Fock_matrix(delta_F0, &pair_p, atoms, job, file);
     for (i = 0; i < dimf_spin; i++) delta_F0[i] = F0[i] - F0_copy[i];
       DestroyDoubleArray(&P0_copy,&dimp_spin,job);
       DestroyDoubleArray(&F0_copy,&dimf_spin,job);
  //print_Fock_matrix(delta_F0, &pair_p, atoms, job, file);

     //for (i = 0; i < dimp_spin; i++) {
     //delta_P0[i] = (P1[i] - P0[i]) * (k_one - job->fock_mixing);
    //}
     //for (i = 0; i < dimf_spin; i++) {
     //delta_F0[i] = (F1[i] - F0[i]) * (k_one - job->fock_mixing);
    //}
     //density_matrix_largest_shell_elements(lse_F0, delta_F0, &pair_p, &pair_p, atoms, shells, job, file);
   }

    if (job->kpoints == 0 && job->density == 2) {
    MPI_File_seek(gh, 0, MPI_SEEK_SET) ;
    MPI_File_write(gh, job_parameters, 41 + job->spin_dim * fermi->nkunique, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write(gh, job_parameters1, 11 + dim3, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_write(gh, &P0[0], dimp_spin, MPI_DOUBLE, MPI_STATUS_IGNORE);
   }

  // ******************************************************************************************
  // * Close files                                                                            *
  // ******************************************************************************************

    //if (crystal->type[0] == 'M' && job->taskid == 0) {
    //fclose(scf_evectors);
    //fclose(scf_evalues);
   //}

    if (job->values == 1  || job->values == 2)
    MPI_File_close(&eh);
    if (job->vectors == 1 || job->vectors == 2)
    MPI_File_close(&fh);
    if (job->density == 1 || job->density == 2)
    MPI_File_close(&gh);

    time5 = MPI_Wtime() - time6;

    MPI_Bcast(&job->energy_change,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if ((job->vectors == 2 && job->density == 2) && (fabs(job->energy_change) < job->scf_tol)) {
    break;
   }

 } // close loop on job->iter

  // ******************************************************************************************
  // * Final Fock/Kohn-Sham matrix build and diagonalisation in SCF run                       *
  // ******************************************************************************************

    if (job->vectors == 2 && job->density == 2) { // only scf does this
    if (job->taskid == 0) 
    printf("ITERATION %d\n",job->iter + 1);
    build_fock_matrix(Fock,&one_ints,fermi,&total_electrons,S1,S2,P0,F0,delta_P0,delta_F0,&pair_p,atom_p,atoms, \
    shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
    //build_fock_matrix(Fock,&dft_grid,&one_ints,fermi,&total_electrons,S1,S2,P0,F0,delta_P0,delta_F0,&pair_p,atom_p,atoms, \
    shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
    time4 = MPI_Wtime();
    if (job->taskid == 0)
    printf("\nTime for SCF calculation %10.4e\n\n",(double)(time4 - time3));
   }
   
  DestroyDoubleArray(&S1,&job->dimf,job);
  if (crystal->type[0] != 'M') { 
  DestroyDoubleArray(&S2,&job->dimf,job);
 }
  free_INT_1E(&one_ints, Function, job, file);
  ////free_dft_grid(&dft_grid,job);
  if (job->sym_adapt == 1) {
  for (iii = 0; iii < job->mixing_order * symmetry->number_of_classes; iii++) {
  DestroyComplexMatrix(&Residual[iii],job);
  DestroyComplexMatrix(&FockHist[iii],job);
 }
 }
  DestroyIntArray(&irrep,&dim3,job);
  DestroyDoubleArray(&eigval,&dim3,job);
  DestroyDoubleArray(&eigval1,&dim3,job);
  DestroyDoubleMatrix(&C_i,job);
  DestroyDoubleMatrix(&C_o,job);
  DestroyDoubleMatrix(&D,job);
  free_k_points(&knet,job);
  DestroyDoubleArray(&P0,&dimp_spin,job);
  DestroyDoubleArray(&F0,&dimf_spin,job);
  DestroyDoubleArray(&delta_P0,&dimp_spin,job);
  DestroyDoubleArray(&delta_F0,&dimf_spin,job);
  DestroyDoubleArray(&P1,&dimp_spin,job);
  DestroyDoubleArray(&F1,&dimf_spin,job);
  DestroyDoubleArray(&P2,&dimp_spin,job);
  free_PAIR_TRAN(&pair_p,job);
  DestroyComplexMatrix(&eigvec,job);
  DestroyComplexMatrix(&eigenvectors,job);
  DestroyDoubleArray(&eigenvalues,&dim1,job);
  DestroyDoubleArray(&Fock,&dimp_spin,job);

}
