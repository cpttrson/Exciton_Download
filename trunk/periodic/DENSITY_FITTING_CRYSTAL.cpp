#include <cstring>
#include "myconstants.h"
#include "USER_DATA.h"
/*
#include "PARALLEL.h"
#include "PAIRS_QUADS.h"
#include "ROTATION_OPERATORS.h"
#include "ROTATIONS_MOLECULE.h"
#include "FOURIER_TRANSFORM.h"
#include "INTEGRALS_2C_MOLECULE.h"
#include "INTEGRALS_3C_MOLECULE.h"
*/
#include "KPOINTS.h"
#include "MATRIX_UTIL.h"
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
int is[3]; //move outside here
double *eigenvalues;
double time1, time2;
char uplo = 'U';
char jobz = 'V';
char NoTrans = 'N', ConjTrans = 'C';
ComplexMatrix *xtrn1, *eigenvectors;
Complex alpha = Complex(k_one, k_zero);
Complex  beta = Complex(k_zero, k_zero);
//Complex alpha, beta;

  Q_LATTICE q_G;
  q_G.last_vector = G->last_vector; 
  q_G.max_vector  = G->max_vector;
  allocate_Q_LATTICE(&q_G, job, file);

  //alpha.real() = k_one;
  //alpha.imag() = k_zero;
  //beta.real() = k_zero;
  //beta.imag() = k_zero;

  AllocateComplexMatrix(&xtrn1,&dim1ax,&dim1ax,job);
  AllocateComplexMatrix(&eigenvectors,&dim1ax,&dim1ax,job);
  AllocateDoubleArray(&eigenvalues,&dim1ax,job);

  time1 = MPI_Wtime();

  generate_q_lattice(&fermi->knet->oblique[q], &q_G, fermi, G, crystal, job, file);
  //two_centre_coulomb1_crystal_test(V_q, R, &q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  ////two_centre_coulomb1_crystal(V_q, R, &q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  DiagonaliseHermitian(&V_q, &eigenvalues, &eigenvectors, &jobz, &uplo, &info);
  for (i = 0; i < dim1ax; i++) {
  /////if (fabs(eigenvalues[i]) < 1e-05) { 
  ////printf("Smallest V_q fabs eigenvalue %3d %3d %9.2e\n",q,i,eigenvalues[i]); eigenvalues[i] = 1.0e10; }
  if (eigenvalues[i] < 1e-05) { /////printf("Smallest V_q eigenvalue %3d %3d %9.2e\n",q,i,eigenvalues[i]); 
  eigenvalues[i] = 1.0e15; }
  if (eigenvalues[i] < 1e-05) continue;
  //if (eigenvalues[i] < 1e-05) { printf("Smallest V_q eigenvalue %3d %3d %9.2e\n",q,i,eigenvalues[i]);  }
    //if (eigenvalues[i] < 1.0e-05) eigenvalues[i] = 100.0;
    for (j = 0; j < dim1ax; j++) {
      xtrn1->a[i][j] = eigenvectors->a[i][j] / eigenvalues[i];
      //xtrn1->a[i][j] = eigenvectors->a[i][j] / sqrt(eigenvalues[i]);
     }
    }
  ResetComplexMatrix(V_q);
  ComplexGEMM1(&ConjTrans, &NoTrans, &alpha, &xtrn1, &eigenvectors, &beta, &V_q);

}

/*
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
*/

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
Complex *three_centre_integrals, *reduced_three_centre_integrals;
Complex *kq_buffer[job->band_dim];
Complex *S1;
ComplexMatrix *V_screen;
ComplexMatrix *scf_eigenvectors_0, *scf_eigenvectors_1;
ComplexMatrix *scf_eigenvectors_k[job->band_dim], *scf_eigenvectors_kq[job->band_dim];
VECTOR_INT kvec[2], qvec[2];
TRIPLE_TRAN triple, triple1;

//INTEGRAL_LIST_COMPLEX integral_list;

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
  two_centre_coulomb1_crystal(V_screen, R, q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
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
  // * Calculate 3-centre integrals and contract with psi-k+q                                 *
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

    //allocate_integral_list_complex(&integral_list, dim456, job, file);
    //AllocateComplexArray(&reduced_three_centre_integrals,&dim456,job);
    AllocateComplexArray(&three_centre_integrals,&dim456,job);
    //ResetComplexArray(reduced_three_centre_integrals,&dim456);
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
      shell_screen3(start_index,V_screen,S1,pair_p,&triple1,atoms,shells,atoms_ax,shells_ax,job,file);
      //for (i = 0; i < nshells; i++) start_index[i] = 1;
      //if (job->taskid == 1) {printf("%2d %2d %2d %2d ",kp,lp,j1,t2);for(i=0;i<nshells;i++) printf("%2d",start_index[i]);printf("\n");}
      ResetComplexArray(three_centre_integrals,&dim456);
      //ResetComplexArray(integral_list.value,&dim456);
      //integral_list.num = 0;
      //three_centre_coulomb1_reversed2_crystal_test(0,&triple1,three_centre_integrals,R,q_G,G,atoms,shells, \
      gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,job,file);
      three_centre_coulomb1_reversed2_crystal_test1(0,&triple1,start_index,three_centre_integrals,R,q_G,G,atoms,shells, \
      gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,job,file);
      ////three_centre_coulomb1_reversed2_crystal_test2(&integral_list,0,&triple1,start_index,three_centre_integrals,R,q_G,G,atoms,shells, \
      gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,job,file);
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
       //   Complex temp2[nband_kq[s]][nd5];
       //
       //   for (m = 0; m < nband_kq[s]; m++) {
       //     for (n = 0; n < integral_list.num; n++) {
       //       l1 = integral_list.j[n];
       //       temp2[m][l1] = scf_eigenvectors_kq[s]->a[k * nband_kq[s] + m][bfposl1 + l1] * temp; // * Psi-k+q
       //      }
       //     }
       //   for (n = 0; n < integral_list.num; n++) {
       //     k1 = integral_list.i[n];
       //     l1 = integral_list.j[n];
       //     a1 = integral_list.k[n];
       //     for (m = 0; m < nband_kq[s]; m++) {
       //       kq_buffer[s][k * dim1 * nband_kq[s] * nd6 + (bfposk1 + k1) * nband_kq[s] * nd6 + m * nd6 + a1] \
       //       += integral_list.value[n] * temp2[m][l1]; // * Psi-k+q
       //      }
       //     }
          for (m = 0; m < nband_kq[s]; m++) {
            //for (n = 0; n < integral_list.num; n++) {
              //k1 = integral_list.i[n];
              //l1 = integral_list.j[n];
              //a1 = integral_list.k[n];
              //temp1 = scf_eigenvectors_kq[s]->a[k * nband_kq[s] + m][bfposl1 + l1] * temp; // * Psi-k+q
              //kq_buffer[s][k * dim1 * nband_kq[s] * nd6 + (bfposk1 + k1) * nband_kq[s] * nd6 + m * nd6 + a1] \
              += integral_list.value[n] * temp1; // * Psi-k+q
             //}
            count = 0;
            for (k1 = 0; k1 < nd4; k1++) {
              for (l1 = 0; l1 < nd5; l1++) {
                temp1 = scf_eigenvectors_kq[s]->a[k * nband_kq[s] + m][bfposl1 + l1] * temp; // * Psi-k+q
                for (a1 = 0; a1 < nd6; a1++) {
                  kq_buffer[s][k * dim1 * nband_kq[s] * nd6 + (bfposk1 + k1) * nband_kq[s] * nd6 + m * nd6 + a1] \
                  += three_centre_integrals[count] * temp1; // * Psi-k+q
                  ////kq_buffer[s][k * dim1 * nband_kq[s] * nd6 + (bfposk1 + k1) * nband_kq[s] * nd6 + m * nd6 + a1] \
                  += three_centre_integrals[count] * scf_eigenvectors_kq[s]->a[k * nband_kq[s] + m][bfposl1 + l1] * temp; // * Psi-k+q
                  count++;
                 }
                }
               }
              }
             } // close loop on k
            } // close loop on s
          time8 += MPI_Wtime() - time7;
         } // close loop on p2
          //DestroyComplexArray(&reduced_three_centre_integrals,&dim456,job);
          DestroyComplexArray(&three_centre_integrals,&dim456,job);

          //free_integral_list_complex(&integral_list,job);

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
   if (job->taskid == 0) printf("%3d integrals_occ_occ_vir_vir2_crystal         %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",\
   job->taskid,time2,time4,time6,time8,time10,time12);
   fflush(stdout);
   if (job->taskid >= 0) fprintf(file.out,"%3d integrals_occ_occ_vir_vir2_crystal         %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",\
   job->taskid,time2,time4,time6,time8,time10,time12);
   fflush(file.out);

}

/*
void integrals_occ_occ_vir_vir2_crystal(int *j1, int *q1, KPOINT_TRAN *knet_little_q_group, MPI_File fh, Complex *integral_buffer2, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, SYMMETRY *symmetry_little_q_group, Q_LATTICE *q_G, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Calculate intermediates needed for density fitting integrals                           *
  // ******************************************************************************************

int i, j, k, l, m, q;
int a1, k1, l1, kp, lp, kq, j2, j3, t2;
int nd4, nd5, nd6;
int k_bz, kq_bz, k_fbz, kq_fbz;
int count, offset, vector_size;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int dim1 = atoms->number_of_sh_bfns_in_unit_cell;
int dimk = knet_little_q_group->unique * job->spin_dim * nbands;
int dim1a, dim456, dim1ax = atoms_ax->number_of_sh_bfns_in_unit_cell;
int op, pm;
int bfposa1, bfposk1, bfposl1;
double time1, time2, time3, time4, time5, time6, time7, time8, time9, time10;
double dot_product;
Complex temp;
Complex *three_centre_integrals, *reduced_three_centre_integrals, *temp1_buffer, *temp2_buffer;
ComplexMatrix *scf_eigenvectors_0, *scf_eigenvectors_1, *scf_eigenvectors_2;
VECTOR_INT kvec[2], qvec[2];
TRIPLE_TRAN triple, triple1;
PAIR_TRAN pair_p;

  pair_p.cutoff = 6.0;
pair_p.cutoff = 8.0;
  fprintf(file.out,"pair_cutoff = %10.4f\n",pair_p.cutoff);
  count_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  //generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  generate_range_selected_pairs(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //print_pairs(&pair_p,atoms,R,job,file);

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
  shell_screen_complex(S1,&pair_p,R,G,q_G,atoms,shells,gaussians,symmetry,crystal,job,file);
  AllocateComplexMatrix(&V_screen,&dim1ax,&dim1ax,job);
  two_centre_coulomb1_crystal(V_screen, R, q_G, atoms_ax, shells_ax, gaussians_ax, crystal, job, file);
  time2 = MPI_Wtime() - time1;

  time3 = MPI_Wtime();
  AllocateComplexMatrix(&scf_eigenvectors_0,&dimk,&dim1,job);
  AllocateComplexMatrix(&scf_eigenvectors_1,&dimk,&dim1,job);
  AllocateComplexMatrix(&scf_eigenvectors_2,&dimk,&dim1,job);

ComplexMatrix *scf_eigenvectors_3;
AllocateComplexMatrix(&scf_eigenvectors_3,&dimk,&dim1,job);

  q = fermi->knet->ibz[*q1];
fprintf(file.out,"q = %3d\n",q);
  qvec[0] = decompose_k_point(fermi->is, q, crystal, job, file);
  for (k = 0; k < knet_little_q_group->unique; k++) {
    k_bz  = knet_little_q_group->ibz[k]; // little group
    k_fbz = fermi->knet->fbz[k_bz];      // full group
//int kbz = fermi->knet->bz[k_bz];
    kvec[0]  = decompose_k_point(fermi->is, k_bz, crystal, job, file);
    kq_bz = compose_k_point(fermi->is, kvec[0], qvec[0], crystal, job, file);
//int kqbz = fermi->knet->bz[kq_bz];
    kq_fbz = fermi->knet->fbz[kq_bz];    // full group
    offset = k * nbands;
    vector_size = nbands * dim1;

    MPI_File_seek(fh,  (k_fbz * dim1 * dim1  + (fermi->bands[0] - 1) * dim1) * sizeof(Complex), MPI_SEEK_SET) ;
    MPI_File_read(fh, &scf_eigenvectors_0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_1->a[offset][0],nbands,k_bz,fermi->knet,\
    atom_p,atoms,R,shells,symmetry,job,file);
    //rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_1->a[offset][0],nbands,kbz,fermi->knet,\
    atom_p,atoms,R,shells,symmetry,job,file);

//rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_3->a[0][0],nbands,0,fermi->knet,\
//atom_p,atoms,R,shells,symmetry,job,file);
//rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_3->a[nbands][0],nbands,2,fermi->knet,\
//atom_p,atoms,R,shells,symmetry,job,file);
//for (i = 0; i < nbands; i++) {
  //for (j = 0; j < dim1; j++) {
    //fprintf(file.out,"k %2d k_fbz kq_fbz %2d %2d k_bz kq_bz %2d %2d %16.10lf %16.10lf   %16.10lf %16.10lf\n",
    //fprintf(file.out,"evec3 k %2d k_fbz kq_fbz %2d %2d k_bz kq_bz %2d %2d %8.2lf %8.2lf   %8.2lf %8.2lf\n",\
    //k,k_fbz,kq_fbz,k_bz,kq_bz, \
    //(scf_eigenvectors_3->a[i][j]).real(),(scf_eigenvectors_3->a[i][j]).imag(), \
    //(scf_eigenvectors_3->a[nbands+i][j]).real(),(scf_eigenvectors_3->a[nbands+i][j]).imag());
   //}
  //fprintf(file.out,"\n");
 //}

    //rotate_psi(&eigenvectors->a[(s * fermi->nkunique + kunique) * dim2][0],&eigenvectors1->a[0][0],nbands,knet->bz[countk],knet,\
    atom_p,atoms,R,shell,symmetry,job,file);

//these rereads are for no symmetry fbz only
if (job->kss == 0) {
MPI_File_seek(fh,  (k_bz * dim1 * dim1  + (fermi->bands[0] - 1) * dim1) * sizeof(Complex), MPI_SEEK_SET) ;
MPI_File_read(fh, &scf_eigenvectors_1->a[offset][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

    MPI_File_seek(fh,  (kq_fbz * dim1 * dim1  + (fermi->bands[0] - 1) * dim1) * sizeof(Complex), MPI_SEEK_SET) ;
    MPI_File_read(fh, &scf_eigenvectors_0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_2->a[offset][0],nbands,kq_bz,fermi->knet,\
    atom_p,atoms,R,shells,symmetry,job,file);
    //rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_2->a[offset][0],nbands,kqbz,fermi->knet,\
    atom_p,atoms,R,shells,symmetry,job,file);

//these rereads are for no symmetry fbz only
if (job->kss == 0) {
MPI_File_seek(fh,  (kq_bz * dim1 * dim1  + (fermi->bands[0] - 1) * dim1) * sizeof(Complex), MPI_SEEK_SET) ;
MPI_File_read(fh, &scf_eigenvectors_2->a[offset][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

    }
  time4 = MPI_Wtime() - time3;

   //for (i=0;i<knet_little_q_group->unique*nbands;i++) { for (int jj=0;jj<dim1;jj++){ fprintf(file.out,"%3d %3d %10.4f %10.4f\n",\
   i,jj,(scf_eigenvectors_2->a[i][jj]).real(), (scf_eigenvectors_2->a[i][jj]).imag());}}

    //rotate_psi(&scf_eigenvectors_0->a[0][0],&scf_eigenvectors_2->a[offset][0],nbands,2,fermi->knet,\
    atom_p,atoms,R,shells,symmetry,job,file);

  nd6 = atoms_ax->bfnnumb_sh[*j1];
  dim1a = knet_little_q_group->unique * nbands * dim1 * nd6; 
  AllocateComplexArray(&temp2_buffer,&dim1a,job);
  ResetComplexArray(temp2_buffer,&dim1a);
  ////count_triples1_reversed(*j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
  ////allocate_TRIPLE_TRAN(&triple, job, file);
allocate_TRIPLE_TRAN(&triple1, job, file);
  ////generate_triples1_reversed(*j1,&triple, atoms, atom_p, symmetry, R, R_tables, job, file);
  //for (j = 0; j < triple.nump; j++) {
for (j = 0; j < pair_p.nump; j++) {
//if (job->taskid == 0) printf("%3d out of %3d\n",j+1,pair_p.nump);
    //q = triple.posn[j];
    //nd4 = atoms->bfnnumb_sh[triple.cell1[q]];
    //nd5 = atoms->bfnnumb_sh[triple.cell2[q]];
int q2;
q2 = pair_p.posn[j];
nd4 = atoms->bfnnumb_sh[pair_p.cell1[q2]];
nd5 = atoms->bfnnumb_sh[pair_p.cell2[q2]];
    dim456 = nd4 * nd5 * nd6;
    AllocateComplexArray(&reduced_three_centre_integrals,&dim456,job);
    AllocateComplexArray(&three_centre_integrals,&dim456,job);
    ResetComplexArray(reduced_three_centre_integrals,&dim456);
    ////three_centre_coulomb1_reversed2(j,&triple,reduced_three_centre_integrals,R,G,atoms,shells, \
    gaussians,atoms_ax, shells_ax,gaussians_ax,crystal,job,file);
    //for (i = 0; i < dim456; i++) fprintf(file.out,"%3d %3d %10.4lf\n",j,i,reduced_three_centre_integrals[i]);
    //for (j2 = 0; j2 < triple.numb[j]; j2++) {
for (j2 = 0; j2 < pair_p.numb[j]; j2++) {
      //kp = triple.cell1[q + j2];
      //lp = triple.cell2[q + j2];
      //t2 = triple.latt2[q + j2];
      kp = pair_p.cell1[q2 + j2];
      lp = pair_p.cell2[q2 + j2];
      t2 = pair_p.latt2[q2 + j2];
      j3 = q2 + j2;
triple1.cell1[0] = kp;
triple1.cell2[0] = lp;
triple1.cell3[0] = *j1;
triple1.latt1[0] = 0;
triple1.latt2[0] = t2;
triple1.latt3[0] = 0;

    time5 = MPI_Wtime();
      int nshells = atoms->nshel_sh[kp] * atoms->nshel_sh[lp] * atoms_ax->nshel_sh[*j1];
      int start_index[nshells];
      int skip = 1;
///      for (i = 0; i < nshells; i++) start_index[i] = 0;
///      shell_screen3(start_index,V_screen,S1,&pair_p,&triple1,atoms,shells,atoms_ax,shells_ax,job,file);
      //for (i = 0; i < nshells; i++) start_index[i] = 1;
      for (i = 0; i < nshells; i++) if (start_index[i] == 1) skip = 0;
for (i = 0; i < nshells; i++) start_index[i] = 1;
skip = 0;
      //printf("skip task %3d j1 %3d q1 %3d j, j2 %3d %3d skip %3d  %3d\n",job->taskid,*j1,*q1,j,j2,skip,i);
      //fflush(stdout);
      //fprintf(file.out,"skip task %3d j1 %3d q1 %3d j, j2 %3d %3d skip %3d\n",job->taskid,*j1,*q1,j,j2,skip);
      //fflush(file.out);
      //if (job->taskid == 0) printf("%3d  %3d out of %3d\n",j2,j+1,pair_p.nump);
      //for (i = 0; i < nshells; i++) fprintf(file.out,"start %3d %3d\n",i,start_index[i]);
      if (skip == 1) continue;
      ////if (triple.latt3[q2+j2] != 0) continue;
      //if (triple.cell1[q2+j2] != triple.cell2[q2+j2] || triple.latt3[q2+j2] != 0) continue;
      //fprintf(file.out,"triple %3d   %3d %3d  %3d  %3d  %3d %3d\n",\
      *j1,triple.cell1[q2+j2],triple.cell2[q2+j2],triple.latt2[q2+j2],triple.cell3[q2+j2],triple.latt1[q2+j2],triple.latt3[q2+j2]);
      ResetComplexArray(three_centre_integrals,&dim456);
//three_centre_coulomb1_reversed2_crystal_test(0,&triple1,three_centre_integrals,R,q_G,G,atoms,shells, \
gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,job,file);
//fprintf(file.out,"T2 %3d %3d %3d   j2 %3d \n",kp,lp,t2,j2);
three_centre_coulomb1_reversed2_crystal_test1(0,&triple1,start_index,three_centre_integrals,R,q_G,G,atoms,shells, \
gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,job,file);
    time6 += MPI_Wtime() - time5;
      ////three_centre_coulomb1_reversed2_crystal(j3,&triple,three_centre_integrals,R,q_G,G,atoms,shells, \
      gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,job,file);
      //for (int i = 0; i < dim456; i++) fprintf(file.out,"%3d %3d %3d   %3d %3d %3d %14.8lf  %14.8lf\n",\
      j,j2,i,kp,lp,t2,(three_centre_integrals[i]).real(),(three_centre_integrals[i]).imag());
      //fprintf(file.out,"j %3d j2 %3d kp %3d lp %3d t2 %3d\n",j,j2,kp,lp,t2);
      ////op = triple.k[q2 + j2];
      ////pm = triple.p[q2 + j2];
      nd4 = atoms->bfnnumb_sh[kp];
      nd5 = atoms->bfnnumb_sh[lp];
      bfposk1 = atoms->bfnposn_sh[kp];
      bfposl1 = atoms->bfnposn_sh[lp];
      //ResetDoubleArray(three_centre_integrals,&dim456);
      //rotate_permute_triple_ax_reversed2(&kp,&lp,j1,&op,&pm,reduced_three_centre_integrals, \
      three_centre_integrals,atom_p,atoms,shells,atoms_ax,shells_ax,symmetry,job,file);
    time7 = MPI_Wtime();
      for (k = 0; k < knet_little_q_group->unique; k++) {
        k_bz = knet_little_q_group->ibz[k];
        kvec[0]  = decompose_k_point(fermi->is, k_bz, crystal, job, file);
        kq_bz = compose_k_point(fermi->is, kvec[0], qvec[0], crystal, job, file);
        dot_product = double_vec_dot(&fermi->knet->cart[kq_bz],&R->vec_ai[t2]); 
//dot_product= 0.0;
        temp = Complex(cos(dot_product), sin(dot_product)); // phase factor for Psi-k+q
        //dot_product = double_vec_dot(&fermi->knet->cart[k_bz],&R->vec_ai[t2]); 
        //temp = Complex(cos(dot_product), -sin(dot_product));   // phase factor for Psi-k^*
        for (m = 0; m < nbands; m++) {
          count = 0;
          for (k1 = 0; k1 < nd4; k1++) {
            for (l1 = 0; l1 < nd5; l1++) {
              for (a1 = 0; a1 < nd6; a1++) {
                //for (m = 0; m < nbands; m++) {
                //temp2_buffer[k * dim1 * nbands * nd6 + (bfposk1 + k1) * nbands * nd6 + m * nd6 + a1] += \
                conj(three_centre_integrals[count]) * scf_eigenvectors_1->a[k * nbands + m][bfposl1 + l1] * temp;
                temp2_buffer[k * dim1 * nbands * nd6 + (bfposk1 + k1) * nbands * nd6 + m * nd6 + a1] += \
                three_centre_integrals[count] * scf_eigenvectors_2->a[k * nbands + m][bfposl1 + l1] * temp; // * Psi-k+q
                count++;
                }
               //count++;
              }
             }
            }
           } // close loop on k
    time8 += MPI_Wtime() - time7;
          //fprintf(file.out,"\n");
           //fprintf(file.out,"\n");
          } // close loop on j2
       DestroyComplexArray(&reduced_three_centre_integrals,&dim456,job);
       DestroyComplexArray(&three_centre_integrals,&dim456,job);
      } // close loop on j

    time9 = MPI_Wtime();
      for (k = 0; k < knet_little_q_group->unique; k++) {
        for (l = 0; l < nbands; l++) {
          for (k1 = 0; k1 < dim1; k1++) {
            for (m = 0; m < nbands; m++) {
              for (a1 = 0; a1 < nd6; a1++) {
                integral_buffer2[k * nbands * nbands * nd6 + l * nbands * nd6 + m * nd6 + a1] += \
                temp2_buffer[k * dim1 * nbands * nd6 + k1 * nbands * nd6 + m * nd6 + a1] * \
                conj(scf_eigenvectors_1->a[k * nbands + l][k1]);  // * Psi-k^*
               }
              }
             }
            }
           }
    time10 += MPI_Wtime() - time9;

      for (k = 0; k < knet_little_q_group->unique; k++) {
        k_bz = knet_little_q_group->ibz[k];
        kvec[0]  = decompose_k_point(fermi->is, k_bz, crystal, job, file);
        kq_bz = compose_k_point(fermi->is, kvec[0], qvec[0], crystal, job, file);
        for (l = 0; l < nbands; l++) {
          for (m = 0; m < nbands; m++) {
            for (a1 = 0; a1 < nd6; a1++) {
              fprintf(file.out,"vir2 j1 %3d *q1 %3d k %3d l,m %3d %3d a1 %3d all  %3d    %14.8lf %14.8lf\n",*j1,*q1,k,l,m,a1,\
              k * nbands * nbands * nd6 + l * nbands * nd6 + m * nd6 + a1,\
              (integral_buffer2[k * nbands * nbands * nd6 + l * nbands * nd6 + m * nd6 + a1]).real(),\
              (integral_buffer2[k * nbands * nbands * nd6 + l * nbands * nd6 + m * nd6 + a1]).imag());
             }
            }
           }
           //fprintf(file.out,"\n");
          }


          ////free_TRIPLE_TRAN(&triple,job);
          free_TRIPLE_TRAN(&triple1,job);
          free_PAIR_TRAN(&pair_p,job);
          DestroyComplexArray(&temp2_buffer,&dim1a,job);
          DestroyComplexMatrix(&scf_eigenvectors_0,job);
          DestroyComplexMatrix(&scf_eigenvectors_1,job);
          DestroyComplexMatrix(&scf_eigenvectors_2,job);
  DestroyComplexArray(&S1,&dimf,job);
  DestroyComplexMatrix(&V_screen,job);
          //time2 = MPI_Wtime() - time1;
          //if (job->taskid == 0) printf("integrals_occ_occ_vir_vir2_crystal         %10.2f\n",time2);
          if (job->taskid == 0) printf("integrals_occ_occ_vir_vir2_crystal         %10.2f %10.2f %10.2f %10.2f %10.2f\n",\
          time2,time4,time6,time8,time10);

}
*/
