#include <cstring>
#include "mycomplex.h"
#include "myconstants.h"
#include "conversion_factors.h"
#include "USER_DATA.h"
#include "PAIRS_QUADS.h"
#include "KPOINTS.h"
#include "TDHF_CRYSTAL.h"
#include "PRINT_MOLECULE.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "FOURIER_TRANSFORM.h"
#include "PARALLEL.h"
#include "INTEGRALS_2C_MOLECULE.h"
#include "INTEGRALS1.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "ROTATIONS_MOLECULE.h"
#include "GW_BSE_MOLECULE.h"
#include "OPTICAL_SPECTRUM.h"

using namespace std;

void optical_spectrum_crystal2(FERMI* fermi, ATOM* atoms, ATOM_TRAN* atom_p, int *numfrag, int *natom, int nat[][2], SHELL* shells, GAUSSIAN* gaussians, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax,CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  //KPOINT_TRAN knet;
  //count_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //if (job->kss == 0)      fermi->nkunique = knet.nktot;
  //else if (job->kss == 1) fermi->nkunique = knet.unique;
  //allocate_k_points(&knet,crystal,job,file);
  //generate_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  //print_knet(&knet, fermi->is, crystal, job, file);
  //allocate_fermi(fermi,atoms,job,file);
  //fermi->knet = &knet;
  //fermi->nktot = knet.nktot;
  //fermi->occupied[0] = fermi->homo[0] - fermi->bands[0] + 1;
  //fermi->occupied[fermi->nkunique] = fermi->homo[1] - fermi->bands[2] + 1;

  KPOINT_TRAN knet;
  count_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  if (job->kss == 0)      fermi->nkunique = knet.nktot;
  else if (job->kss == 1) fermi->nkunique = knet.unique;
  allocate_k_points(&knet,crystal,job,file);
  generate_k_points(&knet,fermi->is,crystal,symmetry,job,file);
  print_knet(&knet, fermi->is, crystal, job, file);
  fermi->knet = &knet;
  fermi->nktot = knet.nktot;

int i, j, k, ii, jj, i1, j1, k1, i2, j2, i3, j3, l, m, n, s, t;
int nk[2];
int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
int ntrans[2], nband[2], ntransitions, nt, nbands;
int vector_size;
int dimk;
int offset, eigval_size, size, seek, count;
double E, E_tot, energy, range = job->energy_range[1] - job->energy_range[0];
double increment = range / (double) (job->npoints + 1);
double *scf_eigenvalues, *GW_eigenvalues, *bse_eigenvalues;
char xx[20] = "/bse_eigenvectors";
char yy[20] = "/bse_eigenvalues";
char zz[20] = "/scf_evec";
char zz4[24] = "/evalfile";
char zz5[24] = "GW_evalues";
char buf1[110], buf2[110];
FILE *bse_evalues, *scf_evalues, *spect;
ComplexMatrix *bse_eigenvectors;
ComplexMatrix *eigvec0, *eigvec1, *eigvec2;
MPI_File fh, gh;

  // ******************************************************************************************
  // * open SCF and BSE eigenvector MPI files and load eigenvalues                            *
  // ******************************************************************************************

  strcpy(buf1,file.scf_eigvec);
  strcat(buf1,zz);
  MPI_File_open(MPI_COMM_WORLD,buf1,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,xx);
  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&gh) ;

  // ******************************************************************************************
  // * Allocate Memory                                                                        *
  // ******************************************************************************************

  int nt_local;
  int dim, dim1, dim3, dim5, dimf;
  int begin_m[job->numtasks], end_m[job->numtasks], receive_m[job->numtasks], offset_m[job->numtasks];
  int begin_p[job->numtasks], end_p[job->numtasks], receive_p[job->numtasks], offset_p[job->numtasks];
  int Function[8] ;
  ComplexMatrix *S_k, *M_k, *tmp, *M_x, *S_x;
  INT_1E one_ints, one_ints_buffer1;
  PAIR_TRAN pair_p;
  nband[0] = fermi->bands[1] - fermi->bands[0] + 1;
  ntrans[0] = (fermi->homo[0] - fermi->bands[0] + 1) * (fermi->bands[1] - fermi->homo[0]);
  ntransitions = ntrans[0];
  nbands = nband[0];
  if (job->spin_dim == 2) {
  nband[1] = fermi->bands[3] - fermi->bands[2] + 1;
  ntrans[1] = (fermi->homo[1] - fermi->bands[2] + 1) * (fermi->bands[3] - fermi->homo[1]);
  ntransitions += ntrans[1];
  nbands += nband[1];
  }
  dimk = fermi->nkunique * job->spin_dim * nbands;
  nt_local = fermi->nktot * job->bse_lim;
  if (job->bse_lim == 0 || fermi->nktot * job->bse_lim > fermi->nktot * ntransitions) nt_local = fermi->nktot * ntransitions;
  //printf("%3d %3d %3d %3d %3d %3d\n",job->bse_lim,ntransitions,fermi->nktot,nt_local,dimk,nbands);
  //if (job->bse_lim == 0) nt_local = fermi->nktot * ntransitions;

  // ******************************************************************************************
  // * Read SCF and BSE eigenvalues from disk                                                 *
  // ******************************************************************************************
 
  AllocateDoubleArray(&scf_eigenvalues,&dimk,job);
  ResetDoubleArray(scf_eigenvalues,&dimk);
  read_scf_GW_eigenvalues_crystal(scf_eigenvalues, fermi, atoms, zz4, job, file);

  eigval_size = ntransitions * fermi->nktot;
  AllocateDoubleArray(&bse_eigenvalues,&eigval_size,job);
  ResetDoubleArray(bse_eigenvalues,&eigval_size);
  read_bse_eigenvalues_crystal(bse_eigenvalues, eigval_size, yy, job, file);

  // ******************************************************************************************
  // * Generate pairs                                                                         *
  // ******************************************************************************************

  count_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  allocate_PAIR_TRAN(&pair_p,atoms,symmetry,R_tables,job,file);
  generate_density_pairs4(&pair_p,atoms,atom_p,symmetry,R,R_tables,job,file);
  //if (job->taskid == 0 && job->print_pairs == 1)
  //print_pairs(&pair_p,atoms,R,job,file);

  // ******************************************************************************************
  // * Generate momentum matrix elements                                                      *
  // ******************************************************************************************

  Function[0] = 0 ;
  Function[1] = 0 ;
  Function[2] = 0 ;
  Function[3] = 1 ;
  Function[4] = 1 ;
  Function[5] = 0 ;
  Function[6] = 0 ;
  Function[7] = 0 ;

  mpi_begin_end(begin_p,end_p,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_pairs(begin_p,end_p,receive_p,offset_p,&pair_p,atoms,job,file);
  mpi_begin_end(begin_m,end_m,pair_p.nump,job->numtasks,job,file);
  mpi_receive_offset_momentum(begin_m, end_m, receive_m, offset_m, &pair_p, atoms, job, file);
  array_dimensions(&dim, &dimf, &pair_p, atoms, job, file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  allocate_INT_1E(&one_ints_buffer1, dim, Function, job, file);
  fock_element_1e(&one_ints_buffer1, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  MPI_Allgatherv(&one_ints_buffer1.Overlap[offset_p[job->taskid]],receive_p[job->taskid],MPI_DOUBLE,&one_ints.Overlap[0],\
  receive_p,offset_p,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgatherv(&one_ints_buffer1.Momentum[offset_m[job->taskid]],receive_m[job->taskid],MPI_DOUBLE,&one_ints.Momentum[0],\
  receive_m,offset_m,MPI_DOUBLE,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer1, Function, job, file);

  // ******************************************************************************************
  // * Allocate memory for optical matrix elements and initialise arrays                      *
  // ******************************************************************************************

  char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
Complex alpha = Complex(k_one, k_zero);
Complex  beta = Complex(k_zero, k_zero);
  //Complex alpha, beta;
  //alpha.real() = k_one;
  //alpha.imag() = k_zero;
  //beta.real() = k_zero;
  //beta.imag() = k_zero;
  double bse_spectrum_buffer[job->spin_dim][2 * job->field_dirs][job->npoints + 1];
  double free_particle_spectrum_buffer[2 * job->spin_dim][job->field_dirs][job->npoints + 1];
  double TRK_sum_rule_buffer[job->field_dirs];
  double TRK_sum_rule1_buffer[job->field_dirs];
  Complex bse_matrix_element_buffer[job->spin_dim][job->field_dirs][nt_local];
  double vol_fac = job->spin_fac * hbar * hbar * hbar * hbar / el / au_to_eV / au_to_eV / au_to_eV / epsilon_0 / m0 / m0 / \
  a0 / a0 / a0 / a0 / a0 / crystal->primitive_cell_volume / fermi->is[0] / fermi->is[1] / fermi->is[2];

  for (s = 0; s < job->spin_dim; s++) {
    for (j3 = 0; j3 < 2 * job->field_dirs; j3++) {
      for (j1 = 0; j1 < job->npoints; j1++) {
        free_particle_spectrum_buffer[s][j3][j1] = k_zero;
       }
      }
     }

  for (s = 0; s < job->spin_dim; s++) {
    for (j3 = 0; j3 < 2 * job->field_dirs; j3++) {
      for (j1 = 0; j1 < job->npoints; j1++) {
        bse_spectrum_buffer[s][j3][j1] = k_zero;
       }
      }
     }

  for (s = 0; s < job->spin_dim; s++) {
    for (j3 = 0; j3 < job->field_dirs; j3++) {
      for (t = 0; t < nt_local; t++) {
        bse_matrix_element_buffer[s][j3][t] = Complex(k_zero, k_zero);
       }
      }
     }

  for (j3 = 0; j3 < job->field_dirs; j3++) {
    TRK_sum_rule_buffer[j3] = k_zero;
   }

  for (j3 = 0; j3 < job->field_dirs; j3++) {
    TRK_sum_rule1_buffer[j3] = k_zero;
   }

  // ******************************************************************************************
  // * Calculate optical spectra and sum rules                                                *
  // ******************************************************************************************

  int begin_k1[job->numtasks], end_k1[job->numtasks];
  mpi_begin_end(begin_k1,end_k1,fermi->knet->nktot,job->numtasks,job,file);

  for (s = 0; s < job->spin_dim; s++) {
    dim = job->spin_dim * fermi->knet->unique * nband[s] ;
    dim3 = atoms->number_of_sh_bfns_in_unit_cell * 3;
    dim5 = nband[s] * 3;
    //nt_local = fermi->nktot * job->bse_lim;
    //if (job->bse_lim == 0) nt_local = fermi->nktot * ntransitions;
    AllocateComplexMatrix(&bse_eigenvectors,&nt_local,&ntrans[s],job);
    double E_trans[ntrans[s]];
    Complex free_particle_matrix_element[job->field_dirs][ntrans[s]];
    AllocateComplexMatrix(&eigvec0,&nband[s],&atoms->number_of_sh_bfns_in_unit_cell,job);
    AllocateComplexMatrix(&eigvec1,&nband[s],&atoms->number_of_sh_bfns_in_unit_cell,job);
    AllocateComplexMatrix(&eigvec2,&nband[s],&atoms->number_of_sh_bfns_in_unit_cell,job);
    AllocateComplexMatrix(&tmp,&atoms->number_of_sh_bfns_in_unit_cell,&nband[s],job);
    AllocateComplexMatrix(&S_k,&atoms->number_of_sh_bfns_in_unit_cell,&atoms->number_of_sh_bfns_in_unit_cell,job);
    AllocateComplexMatrix(&M_k,&dim3,&atoms->number_of_sh_bfns_in_unit_cell,job);
    AllocateComplexMatrix(&S_x,&nband[s],&nband[s],job);
    AllocateComplexMatrix(&M_x,&dim5,&nband[s],job);
    for (j = begin_k1[job->taskid]; j < end_k1[job->taskid]; j++) {
    if (job->taskid == 0) printf("%3d j %3d\n",job->taskid,j);
      k1 = fermi->knet->fbz[j];
      nk[0] = j;
      nk[1] = j;
      fourier_transform_3(&one_ints.Momentum[0], &M_k->a[0][0], fermi->knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      fourier_transform(&one_ints.Overlap[0], &S_k->a[0][0], fermi->knet, nk, &pair_p, R, atoms, shells, symmetry, job, file);
      vector_size = nband[s] * atoms->number_of_sh_bfns_in_unit_cell;
      offset = (fermi->bands[2 * s] - 1) * nbfn * sizeof(Complex);
      MPI_File_seek(fh, (s * fermi->knet->unique + k1) * nbfn * nbfn * sizeof(Complex) + offset, MPI_SEEK_SET) ;
      MPI_File_read(fh, &eigvec0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
      rotate_psi(&eigvec0->a[0][0],&eigvec1->a[0][0],nband[s],j,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
      MPI_File_seek(fh, (s * fermi->knet->unique + k1) * nbfn * nbfn * sizeof(Complex) + offset, MPI_SEEK_SET) ;
      MPI_File_read(fh, &eigvec0->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
      rotate_psi(&eigvec0->a[0][0],&eigvec2->a[0][0],nband[s],j,fermi->knet,atom_p,atoms,R,shells,symmetry,job,file);
      for (i2 = 0; i2 < nband[s]; i2++) {
        for (j2 = 0; j2 < nbfn; j2++) {
          eigvec1->a[i2][j2] *= Complex(-k_one, k_zero);
          //eigvec1->a[i2][j2].imag() *= -k_one;
         }
        }
      ResetComplexMatrix(S_x);
      ResetComplexMatrix(tmp);
      ComplexGEMM3(&NoTrans,&ConjTrans,&nbfn,&nband[s],&nbfn,&alpha,&(S_k->a[0]),&nbfn,&(eigvec2->a[0]),&nbfn,&beta,&(tmp->a[0]),&nband[s]);
      ComplexGEMM3(&NoTrans,&NoTrans,&nband[s],&nband[s],&nbfn,&alpha,&(eigvec2->a[0]),&nbfn,&(tmp->a[0]),&nband[s],&beta,&(S_x->a[0]),\
      &nband[s]);
      ResetComplexMatrix(M_x);
      ResetComplexMatrix(tmp);
      for (j3 = 0; j3 < 3; j3++) {
      dim = j3 * atoms->number_of_sh_bfns_in_unit_cell;
      dim5 = j3 * nband[s];
      ComplexGEMM3(&NoTrans,&Trans,&nbfn,&nband[s],&nbfn,&alpha,&(M_k->a[dim]),&nbfn,&(eigvec2->a[0]),&nbfn,&beta,&(tmp->a[0]),&nband[s]);
      ComplexGEMM3(&NoTrans,&NoTrans,&nband[s],&nband[s],&nbfn,&alpha,&(eigvec1->a[0]),&nbfn,&(tmp->a[0]),&nband[s],&beta,&(M_x->a[dim5]),\
      &nband[s]);
     }
      if (job->verbosity > 1 && j < 6) {
        for (ii = 0; ii < nband[s]; ii++) {
          for (jj = 0; jj < atoms->number_of_sh_bfns_in_unit_cell; jj++) {
            fprintf(file.out,"k %2d %2d %16.10lf %16.10lf   %16.10lf %16.10lf\n",k1,j,(eigvec0->a[ii][jj]).real(),\
            (eigvec0->a[ii][jj]).imag(),(eigvec2->a[ii][jj]).real(),(eigvec2->a[ii][jj]).imag());
           }
          fprintf(file.out,"\n");
         }
      fprintf(file.out,"M_k\n");
      print_complex_matrix(M_k,file);
      fprintf(file.out,"\neigvec1\n");
      print_complex_matrix(eigvec1,file);
      fprintf(file.out,"\neigvec2\n");
      print_complex_matrix(eigvec2,file);
      fprintf(file.out,"tmp\n");
      print_complex_matrix(tmp,file);
      fprintf(file.out,"\nM_x\n");
      print_complex_matrix(M_x,file);
     }

      for (j3 = 0; j3 < job->field_dirs; j3++) {
        for (t = 0; t < ntrans[s]; t++) {
          free_particle_matrix_element[j3][t] = Complex(k_zero, k_zero);
         }
        }

      count = 0;
      for (m = 0; m < fermi->homo[s] - fermi->bands[2 * s] + 1; m++) {
        for (n = fermi->homo[s] - fermi->bands[2 * s] + 1; n < nband[s]; n++) {
          E_trans[count] = scf_eigenvalues[k1 * nbands + s * nband[0] + n] - scf_eigenvalues[k1 * nbands + s * nband[0] + m];
          //fprintf(file.out,"E_t %3d %3d %3d %10.4f %10.4f %10.4f\n",\
          m,n,count,E_trans[count],scf_eigenvalues[k1 * nbands + s * nband[0] + n],scf_eigenvalues[k1 * nbands + s * nband[0] + m]);
          for (j3 = 0; j3 < job->field_dirs; j3++) {
            free_particle_matrix_element[j3][count] = \
           (job->e_field[j3].comp1 * M_x->a[m][n] + \
            job->e_field[j3].comp2 * M_x->a[nband[s] + m][n] + \
            job->e_field[j3].comp3 * M_x->a[2 * nband[s] + m][n]);
           }
          count++; // counter for all transitions
         } // close loop on n
        } // close loop on m

  // ******************************************************************************************
  // * Read BSE eigenvectors from disk                                                        *
  // ******************************************************************************************

  long long lim, memsize;  // still need to fix for large hamiltonians
  memsize = nt_local;
  if (job->bse_tda == 0) {
    //MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
    //for (i = 0; i < job->bse_lim; i++) 
    //MPI_File_read(fh, &bse_eigenvectors->a[job->bse_lim - 1 - i][0], 2 * nt, MPI_DOUBLE, MPI_STATUS_IGNORE);
    //for (t = 0; t < job->bse_lim * fermi->nktot; t++) {
    for (t = 0; t < nt_local; t++) {
      MPI_File_seek(gh, t * ntransitions * sizeof(Complex), MPI_SEEK_SET) ;
      MPI_File_read(gh, &bse_eigenvectors->a[t][0], 2 * memsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
      count = 0;
      for (m = 0; m < ntrans[s]; m++) {
        for (j3 = 0; j3 < job->field_dirs; j3++) {
          bse_matrix_element_buffer[s][j3][t] += free_particle_matrix_element[j3][count] * bse_eigenvectors->a[t][count];
         }
        count++; // counter for all transitions
       } // close loop on m
      for (j1 = 0; j1 < bse_eigenvectors->iRows; j1++) { 
        for (k = 0; k < bse_eigenvectors->iCols; k++) { 
          //fprintf(file.out,"%3d %3d %10.4f\n",j,k,bse_eigenvectors->a[j][k]);
          bse_eigenvectors->a[j1][k] /= two * sqrt(bse_eigenvalues[j1]); 
         }
        } 
       }
      } // close if (job->bse_tda
  else if (job->bse_tda == 1) {
    //long long lim, memsize;  // still need to fix for large hamiltonians
    memsize = ntrans[s];
    for (t = 0; t < nt_local; t++) {
      MPI_File_seek(gh, (t * fermi->nktot * ntransitions + j * ntransitions + s * ntrans[0]) * sizeof(Complex), MPI_SEEK_SET) ;
      MPI_File_read(gh, &bse_eigenvectors->a[t][0], 2 * memsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
      for (j3 = 0; j3 < job->field_dirs; j3++) {
        for (m = 0; m < ntrans[s]; m++) {
          bse_matrix_element_buffer[s][j3][t] += free_particle_matrix_element[j3][m] * bse_eigenvectors->a[t][m];
         } // close loop on m
        } // close loop on j3
       } // close loop on t
       //print_complex_eigenvector_matrix2(bse_eigenvectors, &bse_eigenvalues[j * nt_local], 5, 50, au_to_eV, file);
     } // close else if (job->bse_tda
  for (j1 = 0; j1 < job->npoints; j1++) {
    E = increment * (double) j1 + job->energy_range[0];
    for (j3 = 0; j3 < job->field_dirs; j3++) {
      for (m = 0; m < ntrans[s]; m++) {
        free_particle_spectrum_buffer[s][j3][j1] += vol_fac / E / E * \
        (conj(free_particle_matrix_element[j3][m]) * free_particle_matrix_element[j3][m]).real() * \
        rtpi * exp(-(E_trans[m] - E) * (E_trans[m] - E) / job->linewidth / job->linewidth) / job->linewidth;
       }
      }
     }
    } // close loop over j
  DestroyComplexMatrix(&bse_eigenvectors,job);
  DestroyComplexMatrix(&eigvec0,job);
  DestroyComplexMatrix(&eigvec1,job);
  DestroyComplexMatrix(&eigvec2,job);
  DestroyComplexMatrix(&tmp,job);
  DestroyComplexMatrix(&S_k,job);
  DestroyComplexMatrix(&S_x,job);
  DestroyComplexMatrix(&M_k,job);
  DestroyComplexMatrix(&M_x,job);
 } // close loop on s

  for (s = 0; s < job->spin_dim; s++) {
    for (j3 = 0; j3 < job->field_dirs; j3++) {
      for (j1 = 0; j1 < job->npoints; j1++) {
        //fprintf(file.out,"%3d %3d %3d free %10.4lf \n",s,j3,j1,free_particle_spectrum_buffer[s][j3][j1]);
       }
      }
     }

  double free_particle_spectrum[job->spin_dim][2 * job->field_dirs][job->npoints + 1];
  for (s = 0; s < job->spin_dim; s++) {
    for (j3 = 0; j3 < 2 * job->field_dirs; j3++) { 
      for (j = 0; j < job->npoints + 1; j++) {
        free_particle_spectrum[s][j3][j] = k_zero;
       }
      }
     }
  for (s = 0; s < job->spin_dim; s++) {
    for (j3 = 0; j3 < job->field_dirs; j3++) { // factor of 2 below allows for imaginary, real parts
      MPI_Reduce(&free_particle_spectrum_buffer[s][j3][0],&free_particle_spectrum[s][2 * j3][0],job->npoints + 1, \
      MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     }
    }
  Complex bse_matrix_element[job->spin_dim][job->field_dirs][nt_local];
  for (s = 0; s < job->spin_dim; s++) {
    for (j3 = 0; j3 < job->field_dirs; j3++) {
      for (t = 0; t < nt_local; t++) {
        bse_matrix_element[s][j3][t] = Complex(k_zero, k_zero);
       }
      }
     }
  size = 2 * job->spin_dim * job->field_dirs * nt_local;
  MPI_Reduce(bse_matrix_element_buffer, bse_matrix_element, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (job->taskid == 0) {

  double bse_spectrum[job->spin_dim][2 * job->field_dirs][job->npoints + 1];
  double TRK_sum_rule[job->field_dirs];
  double TRK_sum_rule1[job->field_dirs];
  char field[80];
  for (s = 0; s < job->spin_dim; s++) {
    for (j3 = 0; j3 < job->field_dirs; j3++) {
      for (t = 0; t < job->npoints; t++) {
        //fprintf(file.out,"%3d %3d %3d free %10.4lf \n",s,j3,t,free_particle_spectrum_buffer[s][j3][t]);
      //for (t = 0; t < nt_local; t++) {
        //fprintf(file.out,"%3d %3d bse %10.4lf %10.4lf \n",j3,t,\
        //(bse_matrix_element[s][j3][t]).real(),(bse_matrix_element[s][j3][t]).imag());
       }
      }
     }
  for (s = 0; s < job->spin_dim; s++) {
    for (j3 = 0; j3 < 2 * job->field_dirs; j3++) { 
      for (j = 0; j < job->npoints + 1; j++) {
        bse_spectrum[s][j3][j] = k_zero;
       }
      }
     }

  for (j3 = 0; j3 < job->field_dirs; j3++) TRK_sum_rule[j3]  = k_zero;
  for (j3 = 0; j3 < job->field_dirs; j3++) TRK_sum_rule1[j3] = k_zero;
  sprintf(field," %d     %4.1lf %4.1lf %4.1lf    %4.1lf %4.1lf %4.1lf   %4.1lf %4.1lf %4.1lf",job->field_dirs, \
  job->e_field[0].comp1, job->e_field[0].comp2, job->e_field[0].comp3, \
  job->e_field[1].comp1, job->e_field[1].comp2, job->e_field[1].comp3, \
  job->e_field[2].comp1, job->e_field[2].comp2, job->e_field[2].comp3);
  fprintf(file.out,"|                                      LARGEST OSCILLATOR STRENGTHS                                       |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"| MODE NUMBER   ENERGY (eV) | FIELDS   %66s |\n", field);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  for (s = 0; s < job->spin_dim; s++) {
    for (j1 = 0; j1 < job->npoints; j1++) {
      E = increment * (double) j1 + job->energy_range[0];
      for (j3 = 0; j3 < job->field_dirs; j3++) {
        for (t = 0; t < nt_local; t++) {
          bse_spectrum[s][2 * j3][j1] += vol_fac / E / E * \
         (conj(bse_matrix_element[s][j3][t]) * bse_matrix_element[s][j3][t]).real() * \
          rtpi * exp(-(bse_eigenvalues[t] - E) * (bse_eigenvalues[t] - E) / job->linewidth / job->linewidth) / job->linewidth; // pi * delta
         }
        }
       }
      }
  for (s = 0; s < job->spin_dim; s++) {
    for (t = 0; t < nt_local; t++) {
      for (j3 = 0; j3 < job->field_dirs; j3++) {
        TRK_sum_rule[j3] += job->spin_fac * two / three * \
       (conj(bse_matrix_element[s][j3][t]) * bse_matrix_element[s][j3][t]).real() * bse_eigenvalues[t];
        //fprintf(file.out,"%10.4f ",TRK_sum_rule[j3]);
       }
        //fprintf(file.out,"\n");
      }
        //fprintf(file.out,"\n");
     }

   //for (t = 0; t < nt_local; t++) {
     //if   ((conj(bse_matrix_element[0][t]) * bse_matrix_element[0][t]).real() * bse_eigenvalues[t] > 0.05 || \
           (conj(bse_matrix_element[1][t]) * bse_matrix_element[1][t]).real() * bse_eigenvalues[t] > 0.05 || \
           (conj(bse_matrix_element[2][t]) * bse_matrix_element[2][t]).real() * bse_eigenvalues[t] > 0.05) {
     //fprintf(file.out,"| %4d           %10.4lf ", t + 1, bse_eigenvalues[t] * au_to_eV);
     //for (j3 = 0; j3 < job->field_dirs; j3++) {
       //fprintf(file.out,"|              %10.4lf ", \
       job->spin_fac * two / three * (conj(bse_matrix_element[j3][t]) * bse_matrix_element[j3][t]).real() * \
       bse_eigenvalues[t]);
      //}
     //fprintf(file.out,"|\n");
    //}
   //}

  spect = fopen("Free_particle_optical_spectrum.dat", "w");
  if (spect == NULL) { fprintf(file.out, "cannot open file Free_particle_optical_spectrum.dat in bethe_salpeter\n"); exit(1); }
      for (j = 0; j < job->npoints; j++) {
       energy = job->energy_range[0] + range * (double) j / (double) job->npoints;
        fprintf(spect, "%10.4e   ", energy * au_to_eV);
         for (l = 0; l < job->spin_dim; l++) {
          for (j3 = 0; j3 < 2 * job->field_dirs; j3++) {
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
          for (j3 = 0; j3 < 2 * job->field_dirs; j3++) {
          fprintf(spect, "%12.4e", bse_spectrum[l][j3][j]);
         }
        }
       fprintf(spect,"\n");
      }
     fflush(spect);
     fclose(spect);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"| TRK SUM RULE   %10.4lf |              %10.4lf |              %10.4lf |              %10.4lf |\n", \
  TRK_sum_rule[0] + TRK_sum_rule[1] + TRK_sum_rule[2], TRK_sum_rule[0], TRK_sum_rule[1], TRK_sum_rule[2]);
  fprintf(file.out,"===========================================================================================================\n");
  if (job->taskid == 0 && job->field_dirs < 3)
  fprintf(file.out,"WARNING: Thomas-Reiche-Kun Sum Rule: %3d out of 3 components calculated.\n", job->field_dirs);
  fflush(file.out);

  } // close if (job->taskid == 0

  free_PAIR_TRAN(&pair_p,job);
  free_INT_1E(&one_ints, Function, job, file);
  DestroyDoubleArray(&bse_eigenvalues,&ntransitions,job);
  DestroyDoubleArray(&scf_eigenvalues,&dimk,job);
  MPI_File_close(&fh) ;
  MPI_File_close(&gh) ;
  free_k_points(&knet,job);

}

