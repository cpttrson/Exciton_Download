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
#include <mkl_scalapack.h>
#include "mycomplex.h"
#include "mylogical.h"
#include "conversion_factors.h"
#include "MATRIX_UTIL.h"
#include "INTEGRALS_2C_MOLECULE.h"
#include "PRINT_MOLECULE.h"
#include "SCALAPACK.h"
*/

#include <cstring>
#include "myconstants.h"
#include "USER_DATA.h"
#include "PARALLEL.h"
#include "PAIRS_QUADS.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "ROTATION_OPERATORS.h"
#include "ROTATIONS_MOLECULE.h"
#include "FOURIER_TRANSFORM.h"
#include "INTEGRALS_2C_MOLECULE.h"
#include "INTEGRALS_3C_MOLECULE.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "DENSITY_FITTING_MOLECULE.h"

using namespace std;

void integrals_molecule_IJ_alpha(double *integral_buffer1, double *integral_buffer2, int *ictxt, int *nbsize, FERMI *fermi, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

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

void contract_integrals_molecule_ij_alpha(int *j1, MPI_File fh, double *integral_buffer, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate <a|ij> intermediates needed for density fitting integrals over centre *j1    *
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
      //rotate_permute_triple_ax_reversed2(&kp,&lp,j1,&op,&pm,reduced_three_centre_integrals, \
      three_centre_integrals,atom_p,atoms,shells,atoms_ax,shells_ax,symmetry,job,file);
      rotate_permute_triple_ij_alpha(&kp,&lp,j1,&op,&pm,reduced_three_centre_integrals, \
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

