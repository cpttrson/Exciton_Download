#include <cstring>
#include "myconstants.h"
#include "USER_DATA.h"
#include "KPOINTS.h"
#include "MATRIX_UTIL.h"
#include "PRINT_MOLECULE.h"
#include "BUILD_FOCK_MATRIX.h"
#include "SETUP_RECIPROCAL_LATTICE.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "ROTATIONS_MOLECULE.h"
#include "INTEGRALS_2C_CRYSTAL.h"
#include "INTEGRALS_3C_CRYSTAL.h"
#include "ALLOCATE_MEMORY.h"
#include "DENSITY_FITTING_CRYSTAL.h"

using namespace std;

void generate_coulomb_matrix_inverse_complex(ComplexMatrix *V_q, int q, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k;
int Function[8];
int dim, dimg, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int info = 0;
int begin_k, end_k;
int nk[2];
double *eigenvalues;
double time1, time2;
char uplo = 'U';
char jobz = 'V';
char NoTrans = 'N', ConjTrans = 'C';
ComplexMatrix *xtrn1, *eigenvectors;
Complex alpha = Complex(k_one, k_zero);
Complex  beta = Complex(k_zero, k_zero);

  time1 = MPI_Wtime();

  Q_LATTICE q_G;
  q_G.last_vector = G->last_vector; 
  q_G.max_vector  = G->max_vector;
  allocate_Q_LATTICE(&q_G, job, file);

  AllocateComplexMatrix(&xtrn1,&dim1ax,&dim1ax,job);
  AllocateComplexMatrix(&eigenvectors,&dim1ax,&dim1ax,job);
  AllocateDoubleArray(&eigenvalues,&dim1ax,job);

  generate_q_lattice(&fermi->knet->oblique[q], &q_G, fermi, G, crystal, job, file);
  integrals_crystal_ij(V_q, R, &q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  //fprintf(file.out,"V_q\n");
  //print_complex_matrix(V_q,file);
  DiagonaliseHermitian(&V_q, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);
  for (i = 0; i < dim1ax; i++) {
  if (eigenvalues[i] < 1e-05) { //printf("Smallest V_q eigenvalue %3d %3d %9.2e\n",q,i,eigenvalues[i]); 
  eigenvalues[i] = 1.0e15; }
  if (eigenvalues[i] < 1e-05) continue;
    for (j = 0; j < dim1ax; j++) {
      xtrn1->a[i][j] = eigenvectors->a[i][j] / eigenvalues[i];
     }
    }
  ResetComplexMatrix(V_q);
  ComplexGEMM1(&ConjTrans, &NoTrans, &alpha, &xtrn1, &eigenvectors, &beta, &V_q);

}

void density_fitting_crystal_contract_integrals(int *q1, int *begin_k, int *end_k, int *band_range_k, int *band_range_kq, PAIR_TRAN *pair_p, KPOINT_TRAN *knet_little_q_group, MPI_File fh, Complex **integral_buffer, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, SYMMETRY *symmetry_little_q_group, Q_LATTICE *q_G, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int b, i, j, k, l, m, n, p, p1, p2, p3, q, s;
int a1, j1, k1, l1, kp, lp, kq, j2, j3, t2;
int nd4, nd5, nd6, nshells;
int k_bz, kq_bz, k_fbz, kq_fbz;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int dim1a, dim456, dimk, dimkq;
int dimf, dimp, dimp_spin, dimf_spin;
int bfposa1, bfposk1, bfposl1;
int op, pm;
int nkunique = *end_k - *begin_k;
int seek_point_k, seek_point_kq, seek_spin_offset;
int vector_size_k, vector_size_kq;
int offset_k, offset_kq;
int count, offset;
int nband_k[job->band_dim], nband_kq[job->band_dim], offset_j1[job->band_dim];
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, time12;
double dot_product;
Complex temp, temp1;
Complex *three_centre_integrals;
Complex *kq_buffer[job->band_dim];
Complex *S1;
ComplexMatrix *V_screen;
ComplexMatrix *scf_eigenvectors_0, *scf_eigenvectors_1;
ComplexMatrix *scf_eigenvectors_k[job->band_dim], *scf_eigenvectors_kq[job->band_dim];
VECTOR_INT kvec[2], qvec[2];
TRIPLE_TRAN triple, triple1;

  time2 = k_zero;
  time4 = k_zero;
  time6 = k_zero;
  time8 = k_zero;
  time10 = k_zero;
  time12 = k_zero;

  // ******************************************************************************************
  // * Calculate screening integrals                                                          *
  // ******************************************************************************************

  time1 = MPI_Wtime();
  sh_array_dimensions(&dimp,&dimf,pair_p,atoms,job,file); 
  dimp_spin = dimp * job->spin_dim;
  dimf_spin = dimf * job->spin_dim;
  AllocateComplexArray(&S1,&dimf_spin,job);
  shell_screen_complex(S1,pair_p,R,G,q_G,atoms,shells,gaussians,symmetry,crystal,job,file);
  AllocateComplexMatrix(&V_screen,&dim1ax,&dim1ax,job);
  integrals_crystal_ij(V_screen, R, q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  //print_pairs(pair_p,atoms,R,job,file);
  time2 = MPI_Wtime() - time1;
 
  // ******************************************************************************************
  // * Read wave functions                                                                    *
  // ******************************************************************************************

  time3 = MPI_Wtime();
  q = fermi->knet->ibz[*q1];
  qvec[0] = decompose_k_point(fermi->is, q, crystal, job, file);
  //fprintf(file.out,"q = %3d\n",q);
  for (s = 0; s < job->band_dim; s++) {
    nband_k[s]  = band_range_k[2 * s + 1]  - band_range_k[2 * s];
    nband_kq[s] = band_range_kq[2 * s + 1] - band_range_kq[2 * s];
    vector_size_k  = nband_k[s]  * dim1;
    vector_size_kq = nband_kq[s] * dim1;
    //fprintf(file.out,"%3d %3d %3d %3d %3d %3d\n",job->band_dim,s,nband_k[s],nband_kq[s],band_range_k[2 * s],band_range_kq[2 * s]);
    dimk  = nkunique * nband_k[s];
    dimkq = nkunique * nband_kq[s];
    AllocateComplexMatrix(&scf_eigenvectors_k[s],&dimk,&dim1,job);
    AllocateComplexMatrix(&scf_eigenvectors_kq[s],&dimkq,&dim1,job);
    AllocateComplexMatrix(&scf_eigenvectors_0,&dimk,&dim1,job);
    AllocateComplexMatrix(&scf_eigenvectors_1,&dimkq,&dim1,job);
    seek_spin_offset = 0;
    if (job->spin_dim == 2 && (((s + 1) / 4) * 4 == s + 1 || ((s + 2) / 4) * 4 == s + 2)) seek_spin_offset = fermi->nkunique * dim1 * dim1;
    //fprintf(file.out,"spin_offset %3d %3d\n",s,seek_spin_offset);
    for (k = *begin_k; k < *end_k; k++) {
      k_bz  = knet_little_q_group->ibz[k]; // little group
      k_fbz = fermi->knet->fbz[k_bz];      // full group
      kvec[0]  = decompose_k_point(fermi->is, k_bz, crystal, job, file);
      kq_bz = compose_k_point(fermi->is, kvec[0], qvec[0], crystal, job, file);
      kq_fbz = fermi->knet->fbz[kq_bz];    // full group
      offset_k  = (k  - *begin_k) * nband_k[s];
      offset_kq = (k  - *begin_k) * nband_kq[s];
      seek_point_k = seek_spin_offset + k_fbz * dim1 * dim1  + band_range_k[2 * s] * dim1;
      MPI_File_seek(fh,  seek_point_k * sizeof(Complex), MPI_SEEK_SET) ;
      MPI_File_read(fh, &scf_eigenvectors_0->a[0][0], 2 * vector_size_k, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_k[s]->a[offset_k][0],nband_k[s],k_bz,fermi->knet,\
      atom_p,atoms,R,shells,symmetry,job,file);

      seek_point_kq = seek_spin_offset + kq_fbz * dim1 * dim1  + band_range_kq[2 * s] * dim1;
      MPI_File_seek(fh, seek_point_kq * sizeof(Complex), MPI_SEEK_SET) ;
      MPI_File_read(fh, &scf_eigenvectors_1->a[0][0], 2 * vector_size_kq, MPI_DOUBLE, MPI_STATUS_IGNORE);
      rotate_psi(&scf_eigenvectors_1->a[0][0],&scf_eigenvectors_kq[s]->a[offset_kq][0],nband_kq[s],kq_bz,fermi->knet,\
      atom_p,atoms,R,shells,symmetry,job,file);

      if (job->taskid == 0 && job->verbosity > 1) {
      fprintf(file.out,"q1 %2d k %2d bandsk,kq %4d %4d k_fbz kq_fbz %2d %2d k_bz kq_bz %2d %2d \n",\
      *q1,k,seek_point_k,seek_point_kq,k_fbz,kq_fbz,k_bz,kq_bz);
      for (i = 0; i < nband_k[s]; i++) {
        for (j = 0; j < dim1; j++) {
          fprintf(file.out,"s %2d k %2d k_fbz kq_fbz %2d %2d k_bz kq_bz %2d %2d %8.2lf %8.2lf   %8.2lf %8.2lf\n",\
          s,k,k_fbz,kq_fbz,k_bz,kq_bz, \
          (scf_eigenvectors_k[s]->a[i][j]).real(),(scf_eigenvectors_k[s]->a[i][j]).imag(), \
          (scf_eigenvectors_kq[s]->a[i][j]).real(),(scf_eigenvectors_kq[s]->a[i][j]).imag());
         }
          fprintf(file.out,"\n");
         }
       }

     } // close loop on k
    DestroyComplexMatrix(&scf_eigenvectors_0,job);
    DestroyComplexMatrix(&scf_eigenvectors_1,job);
  } // close loop on s
    time4 = MPI_Wtime() - time3;

  // ******************************************************************************************
  // * Calculate three-centre integrals and contract with psi-k+q                             *
  // ******************************************************************************************

  triple1.tot = 1;
  allocate_TRIPLE_TRAN(&triple1, job, file);
  for (s = 0; s < job->band_dim; s++) offset_j1[s] = 0;
  for (j1 = 0; j1 < atoms_ax->number_of_atoms_in_unit_cell; j1++) {
    nd6 = atoms_ax->bfnnumb_sh[j1];
    for (s = 0; s < job->band_dim; s++) {
      dim1a = nkunique * nband_kq[s] * dim1 * nd6; 
      AllocateComplexArray(&kq_buffer[s],&dim1a,job);
      ResetComplexArray(kq_buffer[s],&dim1a);
     }
  for (p = 0; p < pair_p->nump; p++) {
    p1 = pair_p->posn[p];
    nd4 = atoms->bfnnumb_sh[pair_p->cell1[p1]];
    nd5 = atoms->bfnnumb_sh[pair_p->cell2[p1]];
    dim456 = nd4 * nd5 * nd6;
    AllocateComplexArray(&three_centre_integrals,&dim456,job);
    for (p2 = 0; p2 < pair_p->numb[p]; p2++) {
      kp = pair_p->cell1[p1 + p2];
      lp = pair_p->cell2[p1 + p2];
      t2 = pair_p->latt2[p1 + p2];
      triple1.cell1[0] = kp;
      triple1.cell2[0] = lp;
      triple1.cell3[0] = j1;
      triple1.latt1[0] = 0;
      triple1.latt2[0] = t2;
      triple1.latt3[0] = 0;
      time5 = MPI_Wtime();
      nshells = atoms->nshel_sh[kp] * atoms->nshel_sh[lp] * atoms_ax->nshel_sh[j1];
      int start_index[nshells];
      for (i = 0; i < nshells; i++) start_index[i] = 0;
      integral_screen_crystal_ija(start_index,V_screen,S1,pair_p,&triple1,atoms,shells,atoms_ax,shells_ax,job,file);
      //shell_screen3(start_index,V_screen,S1,pair_p,&triple1,atoms,shells,atoms_ax,shells_ax,job,file);
      ResetComplexArray(three_centre_integrals,&dim456);
      integrals_crystal_ija(0,&triple1,start_index,three_centre_integrals,R,q_G,G,atoms,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,job,file);
      time6 += MPI_Wtime() - time5;
      nd4 = atoms->bfnnumb_sh[kp];
      nd5 = atoms->bfnnumb_sh[lp];
      bfposk1 = atoms->bfnposn_sh[kp];
      bfposl1 = atoms->bfnposn_sh[lp];
      time7 = MPI_Wtime();
      for (s = 0; s < job->band_dim; s++) {
        for (k = 0; k < nkunique; k++) {
          k_bz = knet_little_q_group->ibz[*begin_k + k];
          kvec[0]  = decompose_k_point(fermi->is, k_bz, crystal, job, file);
          kq_bz = compose_k_point(fermi->is, kvec[0], qvec[0], crystal, job, file);
          dot_product = double_vec_dot(&fermi->knet->cart[kq_bz],&R->vec_ai[t2]); 
          temp = Complex(cos(dot_product), sin(dot_product)); // phase factor for Psi-k+q
          for (m = 0; m < nband_kq[s]; m++) {
            count = 0;
            for (k1 = 0; k1 < nd4; k1++) {
              for (l1 = 0; l1 < nd5; l1++) {
                temp1 = scf_eigenvectors_kq[s]->a[k * nband_kq[s] + m][bfposl1 + l1] * temp; // * Psi-k+q
                for (a1 = 0; a1 < nd6; a1++) {
                  kq_buffer[s][k * dim1 * nband_kq[s] * nd6 + (bfposk1 + k1) * nband_kq[s] * nd6 + m * nd6 + a1] \
                  += three_centre_integrals[count] * temp1; // * Psi-k+q
                  count++;
                 }
                }
               }
              }
             } // close loop on k
            } // close loop on s
          time8 += MPI_Wtime() - time7;
         } // close loop on p2
         DestroyComplexArray(&three_centre_integrals,&dim456,job);
        } //   close loop on p

  // ******************************************************************************************
  // * Contract with psi-k*                                                                   *
  // ******************************************************************************************

  time9 = MPI_Wtime();
  for (s = 0; s < job->band_dim; s++) {
    for (k = 0; k < nkunique; k++) {
      for (l = 0; l < nband_k[s]; l++) {
        for (k1 = 0; k1 < dim1; k1++) {
          for (m = 0; m < nband_kq[s]; m++) {
            for (a1 = 0; a1 < nd6; a1++) {
              integral_buffer[s][k * nband_k[s] * nband_kq[s] * dim1ax + l * nband_kq[s] * dim1ax + m * dim1ax + offset_j1[s] + a1] \
              += kq_buffer[s][k * dim1 * nband_kq[s] * nd6 + k1 * nband_kq[s] * nd6 + m * nd6 + a1] * \
              conj(scf_eigenvectors_k[s]->a[k * nband_k[s] + l][k1]);  // * Psi-k^*
             }
            }
           }
          }
         }
        }
  time10 += MPI_Wtime() - time9;

   for (s = 0; s < job->band_dim; s++) {
     dim1a = nkunique * nband_kq[s] * dim1 * nd6; 
     DestroyComplexArray(&kq_buffer[s],&dim1a,job);
     offset_j1[s] += nd6;
    }
   } // close loop on j1
   free_TRIPLE_TRAN(&triple1,job);
   for (s = 0; s < job->band_dim; s++) {
   DestroyComplexMatrix(&scf_eigenvectors_k[s],job);
   DestroyComplexMatrix(&scf_eigenvectors_kq[s],job);
   } // close loop on s

   DestroyComplexArray(&S1,&dimf_spin,job);
   DestroyComplexMatrix(&V_screen,job);

  // ******************************************************************************************
  // * Print contracted integrals                                                             *
  // ******************************************************************************************

  if (job->taskid == 0 && job->verbosity > 1) 
  for (s = 0; s < job->band_dim; s++) {
    for (k = 0; k < nkunique; k++) {
      k_bz = knet_little_q_group->ibz[*begin_k + k];
      kvec[0]  = decompose_k_point(fermi->is, k_bz, crystal, job, file);
      kq_bz = compose_k_point(fermi->is, kvec[0], qvec[0], crystal, job, file);
      for (l = 0; l < nband_k[s]; l++) {
        for (m = 0; m < nband_kq[s]; m++) {
          for (a1 = 0; a1 < dim1ax; a1++) {
            fprintf(file.out,"VIR2 *q1 %3d s %3d k %3d l,m %3d %3d a1 %3d all  %3d    %14.8lf %14.8lf\n",*q1,s,*begin_k + k,l,m,a1,\
            k * nband_k[s] * nband_kq[s] * dim1ax + l * nband_kq[s] * dim1ax + m * dim1ax + a1,\
            (integral_buffer[s][k * nband_k[s] * nband_kq[s] * dim1ax + l * nband_kq[s] * dim1ax + m * dim1ax + a1]).real(),\
            (integral_buffer[s][k * nband_k[s] * nband_kq[s] * dim1ax + l * nband_kq[s] * dim1ax + m * dim1ax + a1]).imag());
           }
          }
         }
        fprintf(file.out,"\n");
       }
      }

   time12 = MPI_Wtime() - time1;
   if (job->taskid == 0) printf("%3d integrals_occ_occ_vir_vir_crystal          %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",\
   job->taskid,time2,time4,time6,time8,time10,time12);
   fflush(stdout);

}

void density_fitting_crystal_rotate_integrals(int *offset_k, int conjugate, KQPOINT_TRAN *kq_pair, ComplexMatrix *S1, ComplexMatrix *S2, Complex **integral_buffers, Complex *rotated_integrals, int *nband_k, int *nband_kq, int *count_bz, ATOM_TRAN *atom_p, ATOM *atoms_ax, SHELL *shells_ax, SYMMETRY *symmetry, FILES file, JOB_PARAM *job)

{

int a1, j1, k1, k2, l1, l2 ;
int dim1 = atoms_ax->number_of_atoms_in_unit_cell;
int dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int dim_rot;
int k_bz, kq_bz;
int offset, offset_j1;
int nd6, op;
int j1_inv[dim1], j2;
int offset_j1_inv[dim1];
Complex *integral_buffer_rotated;

  k_bz  = kq_pair->bz1[*count_bz];
  kq_bz = kq_pair->bz2[*count_bz];
  dim_rot = *nband_k * *nband_kq * dim1ax;
  AllocateComplexArray(&integral_buffer_rotated,&dim_rot,job);
  op  = symmetry->inverse[kq_pair->opr[*count_bz]];

  for (j1 = 0; j1 < dim1;j1++) offset_j1_inv[j1] = 0;
    for (j1 = 0; j1 < dim1; j1++ ) {
      for (j2 = 0; j2 < dim1; j2++ ) {
        if (atom_p->K[j2 * symmetry->number_of_operators + op] == j1) j1_inv[j1] = j2;
       }
      }
  for (j1 = 0; j1 < dim1; j1++ ) {
    for (j2 = 0; j2 < j1_inv[j1]; j2++ ) {
      offset_j1_inv[j1] += atoms_ax->bfnnumb_sh[j2];
     }
    }

  offset_j1 = 0;
  for (k1 = 0; k1 < *nband_k; k1++) {
    for (l1 = 0; l1 < *nband_kq; l1++) {
      offset = 0;
      for (j1 = 0; j1 < dim1; j1++) { 
        rotate_single(op, j1, &integral_buffers[0][*offset_k + offset_j1_inv[j1] + k1 * *nband_kq * dim1ax + l1 * dim1ax],\
        &integral_buffer_rotated[k1 * *nband_kq * dim1ax + l1 * dim1ax + offset], atoms_ax, shells_ax, symmetry, job, file);
        offset += atoms_ax->bfnnumb_sh[j1];
       } // close loop on j1
      } // close loop on l1
     } // close loop on k1

  Complex overlap;
  ResetComplexArray(rotated_integrals,&dim_rot);
  for (k1 = 0; k1 < *nband_k; k1++) {
    for (l1 = 0; l1 < *nband_kq; l1++) {
      for (k2 = 0; k2 < *nband_k; k2++) {
        for (l2 = 0; l2 < *nband_kq; l2++) {
          if (conjugate == 0)      overlap =      S1->a[k1][k2] * S2->a[l2][l1];
          else if (conjugate == 1) overlap = conj(S1->a[k1][k2] * S2->a[l2][l1]);
          if (overlap.real() * overlap.real() + overlap.imag() * overlap.imag() < 1.0e-09) continue;
          for (a1 = 0; a1 < dim1ax; a1++) {
            rotated_integrals[k1 * *nband_kq * dim1ax + l1 * dim1ax + a1] += \
            integral_buffer_rotated[k2 * *nband_kq * dim1ax + l2 * dim1ax + a1] * overlap;
            //fprintf(file.out,"%3d %3d %3d %3d %3d %3d %3d %f %f\n",k1,l1,k2,l2,a1,*nband_k,*nband_kq,overlap.real(),overlap.imag());
            //fflush(file.out);
            //if (a1 == 22) fprintf(file.out,"%3d %3d %3d %3d %3d %3d  %14.8f %14.8f %14.8f %14.8f\n",conjugate,k1,l1,k2,l2,a1,\
            (rotated_integrals[k1 * *nband_kq * dim1ax + l1 * dim1ax + a1]).real(),\
            (rotated_integrals[k1 * *nband_kq * dim1ax + l1 * dim1ax + a1]).imag(),overlap.real(),overlap.imag());
            //if (a1 == 22) fprintf(file.out,"%3d %3d %3d %3d %3d %3d  %14.8f %14.8f %14.8f %14.8f\n",conjugate,k1,l1,k2,l2,a1,\
            (integral_buffer_rotated[k2 * *nband_kq * dim1ax + l2 * dim1ax + a1]).real(),\
            (integral_buffer_rotated[k2 * *nband_kq * dim1ax + l2 * dim1ax + a1]).imag(),overlap.real(),overlap.imag());
           } // close loop on a1
          } // close loop on l2
         } // close loop on k2
        } // close loop on l1
       } // close loop on k1
    DestroyComplexArray(&integral_buffer_rotated,&dim_rot,job);

}

