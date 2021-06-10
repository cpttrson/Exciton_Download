#include <cstring>
#include "mycomplex.h"
#include "myconstants.h"
#include "conversion_factors.h"
#include "USER_DATA.h"
#include "PAIRS_QUADS.h"
#include "PRINT_MOLECULE.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "FOURIER_TRANSFORM.h"
#include "PARALLEL.h"
#include "INTEGRALS_2C_MOLECULE.h"
#include "INTEGRALS1.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "ROTATIONS_MOLECULE.h"
#include "GW_BSE_MOLECULE.h"
#include "OPTICAL_SPECTRUM_MOLECULE.h"

using namespace std;

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

void dipole_matrix_elements_molecule(DoubleMatrix *M_x, FERMI*fermi, ATOM* atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN* gaussians, CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, i1, j1, i3, j3, s;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int neigenvalues, ntransitions;
int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
int vector_size = nbands * atoms->number_of_sh_bfns_in_unit_cell;
int block_size = vector_size * sizeof(Complex);
int dim, dimf, dimk, dim1 = atoms->number_of_sh_bfns_in_unit_cell, dim3 = nbfn * 3, dim4 = nbands * 3, dim5;
int Function[8];
double alpha = k_one, beta = k_zero;
size_t result;
char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
char zz2[24] = "scf_evectors";
PAIR_TRAN pair_p;
DoubleMatrix *M_k, *tmp;
DoubleMatrix *eigvec;
ComplexMatrix *scf_eigenvectors;
INT_1E one_ints, one_ints_buffer;
FILE *scf_evectors;

  ntransitions = (nbands - fermi->occupied[0]) * fermi->occupied[0];
  if (ntransitions <= 0) {
  if (job->taskid == 0)
  fprintf(file.out,"Error in RPA_molecule. HOMO level not in range of bands. ntransitions = %d\n",ntransitions);
  MPI_Finalize();
  exit(0);
 }

  // ******************************************************************************************
  // * Count and generate pairs for excited state Hamiltonian                                 *
  // ******************************************************************************************

  count_density_pairs4(&pair_p, atoms, atom_p, symmetry, R, R_tables, job, file);
  allocate_PAIR_TRAN(&pair_p, atoms, symmetry, R_tables, job, file);
  generate_density_pairs4(&pair_p, atoms, atom_p, symmetry, R, R_tables, job, file);
  //print_pairs(&pair_p, atoms, R, job, file);
  
  // ******************************************************************************************
  // * Allocate memory for dipole matrix elements in AO basis                                 *
  // ******************************************************************************************

  AllocateDoubleMatrix(&M_k,&dim3,&nbfn,job);
  AllocateDoubleMatrix(&tmp,&nbfn,&nbands,job);
  ResetDoubleMatrix(tmp);

  // ******************************************************************************************
  // * Generate dipole matrix elements                                                        *
  // ******************************************************************************************

  Function[0] = 0;
  Function[1] = 0;
  Function[2] = 0;
  Function[3] = 0;
  Function[4] = 0;
  Function[5] = 0;
  Function[6] = 1;
  Function[7] = 0;

  array_dimensions(&dim, &dimf, &pair_p, atoms, job, file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  fock_element_1e2(&one_ints, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);

  //allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  //fock_element_1e(&one_ints_buffer, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //fock_element_1e1(&one_ints_buffer, dim, &pair_p, pair_p.nump, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //MPI_Allreduce(&one_ints_buffer.Dipole[0],&one_ints.Dipole[0],3*dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //free_INT_1E(&one_ints_buffer, Function, job, file);

  // ******************************************************************************************
  // * Initialise parameters for parallel run                                                 *
  // ******************************************************************************************
  
  //strcpy(buf2,file.directory1);
  //strcat(buf2,yy);
  //MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;

  AllocateDoubleMatrix(&eigvec,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);

  dimk = job->spin_dim * nbands;
  vector_size = dimk * dim1;
  AllocateComplexMatrix(&scf_eigenvectors,&dimk,&dim1,job);

  // ******************************************************************************************
  // * Generate optical matrix elements over Molecular Orbitals                               *
  // ******************************************************************************************
  
  if (job->taskid == 0) {

  //for (s = 0; s < job->spin_dim; s++) {

    rotate_permute_expand_tensor_3(&one_ints.Dipole[0], &M_k->a[0][0], &pair_p, R, atoms, shells, symmetry, job, file);

    if (job->taskid == 0) {
    scf_evectors = fopen(zz2, "rb");
    fseek(scf_evectors, dim1 * (fermi->bands[0] - 1) * sizeof(Complex),SEEK_SET);
    result = fread(&scf_eigenvectors->a[0][0],sizeof(Complex),vector_size,scf_evectors);
    fclose(scf_evectors);
   }

  for (i = 0; i < dimk; i++) {
  for (j = 0; j < dim1; j++) {
  eigvec->a[i][j] = (scf_eigenvectors->a[i][j]).real();
 }
 }
  for (i3 = 0; i3 < 3; i3++) {
  dim  = i3 * atoms->number_of_sh_bfns_in_unit_cell;
  dim5 = i3 * nbands;
  DoubleGEMM3(&NoTrans,&Trans,&nbfn,&nbands,&nbfn,&alpha,&(M_k->a[dim]),&nbfn,&eigvec->a[0],&nbfn,&beta,&(tmp->a[0]),&nbands);
  DoubleGEMM3(&NoTrans,&NoTrans,&nbands,&nbands,&nbfn,&alpha,&eigvec->a[0],&nbfn,&(tmp->a[0]),&nbands,&beta,&(M_x->a[dim5]),&nbands);
 }
 
 //} // close loop on s
 
  //fprintf(file.out,"DIPOLE MATRIX\n");
  //print_double_matrix2(M_x,0,6,1.0,file);

 } // close if (job->taskid

  free_PAIR_TRAN(&pair_p,job);
  free_INT_1E(&one_ints, Function, job, file);

}

void fragment_transition_dipole_moment_molecule(FERMI *fermi, int *numfrag, int *fragment_size, int fragment_atom[][2], DoubleMatrix *bse_eigenvectors, double *bse_eigenvalues, ATOM *atoms, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN* gaussians, CRYSTAL* crystal, SYMMETRY* symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

int i, j, k, i1, j1, i3, j3, l, m, n, s, t;
int nbands = fermi->bands[1] - fermi->bands[0] + 1;
int ntransitions = (nbands - fermi->occupied[0]) * fermi->occupied[0];
int nbfn = atoms->number_of_sh_bfns_in_unit_cell;
int dim, dimf, dimk, dim1 = atoms->number_of_sh_bfns_in_unit_cell, dim3 = nbfn * 3, dim4 = nbands * 3, dim5;
int vector_size;
int Function[8];
int count;
double fragment_DM[*numfrag][3][ntransitions];
double nuclear_dipole[*numfrag][3];
double bse_matrix_element[job->field_dirs][ntransitions];
double free_particle_matrix_element[job->field_dirs][ntransitions];
size_t result;
char zz2[24] = "scf_evectors";
FILE *scf_evectors;
double alpha = k_one, beta = k_zero;
char ConjTrans = 'C', Trans = 'T', NoTrans = 'N';
DoubleMatrix *eigvec, *M_k, *M_x, *tmp;
ComplexMatrix *scf_eigenvectors;
INT_1E one_ints, one_ints_buffer;
PAIR_TRAN pair_p;

  // ******************************************************************************************
  // * Allocate memory for dipole matrix elements in AO basis                                 *
  // ******************************************************************************************

  AllocateDoubleMatrix(&M_k,&dim3,&nbfn,job);

  // ******************************************************************************************
  // * Count and generate pairs for dipole matrix elements                                    *
  // ******************************************************************************************

  count_density_pairs4(&pair_p, atoms, atom_p, symmetry, R, R_tables, job, file);
  allocate_PAIR_TRAN(&pair_p, atoms, symmetry, R_tables, job, file);
  generate_density_pairs4(&pair_p, atoms, atom_p, symmetry, R, R_tables, job, file);
  print_pairs(&pair_p, atoms, R, job, file);
  
  // ******************************************************************************************
  // * Generate dipole matrix elements                                                        *
  // ******************************************************************************************

  Function[0] = 0;
  Function[1] = 0;
  Function[2] = 0;
  Function[3] = 0;
  Function[4] = 0;
  Function[5] = 0;
  Function[6] = 1;
  Function[7] = 0;

  array_dimensions(&dim, &dimf, &pair_p, atoms, job, file);
  allocate_INT_1E(&one_ints, dim, Function, job, file);
  //allocate_INT_1E(&one_ints_buffer, dim, Function, job, file);
  fock_element_1e2(&one_ints, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //fock_element_1e(&one_ints_buffer, dim, &pair_p, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //fock_element_1e1(&one_ints, dim, &pair_p, &pair_p.nump, Function, R, G, atoms, shells, gaussians, crystal, job, file);
  //MPI_Allreduce(&one_ints_buffer.Dipole[0],&one_ints.Dipole[0],3*dim,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //free_INT_1E(&one_ints_buffer, Function, job, file);
  rotate_permute_expand_tensor_3(&one_ints.Dipole[0], &M_k->a[0][0], &pair_p, R, atoms, shells, symmetry, job, file);
  free_PAIR_TRAN(&pair_p,job);
  free_INT_1E(&one_ints, Function, job, file);

  // ******************************************************************************************
  // * Transition dipole moments of constituent molecular fragments                           *
  // ******************************************************************************************

  if (job->taskid == 0) {

  dimk = job->spin_dim * nbands;
  vector_size = dimk * dim1;
  AllocateComplexMatrix(&scf_eigenvectors,&dimk,&dim1,job);
  AllocateDoubleMatrix(&eigvec,&nbands,&atoms->number_of_sh_bfns_in_unit_cell,job);
  AllocateDoubleMatrix(&M_x,&dim4,&nbands,job);
  AllocateDoubleMatrix(&tmp,&nbfn,&nbands,job);
  ResetDoubleMatrix(tmp);

  scf_evectors = fopen(zz2, "rb");
  fseek(scf_evectors, dim1 * (fermi->bands[0] - 1) * sizeof(Complex),SEEK_SET);
  result = fread(&scf_eigenvectors->a[0][0],sizeof(Complex),vector_size,scf_evectors);
  fclose(scf_evectors);

  for (i = 0; i < *numfrag; i++) {
    ResetDoubleMatrix(eigvec);
    ResetDoubleMatrix(M_x);
    nuclear_dipole[i][0] = k_zero;
    nuclear_dipole[i][1] = k_zero;
    nuclear_dipole[i][2] = k_zero;
    for (j = 0; j < fragment_size[i]; j++) {
      for (k = 0; k < dimk; k++) {
        for (l = atoms->bfnposn_sh[fragment_atom[j][i]]; \
          l < atoms->bfnposn_sh[fragment_atom[j][i]] + atoms->bfnnumb_sh[fragment_atom[j][i]]; l++) {
          eigvec->a[k][l] = (scf_eigenvectors->a[k][l]).real();
          //fprintf(file.out,"%3d %3d %3d  %10.4f\n",i,k,l,eigvec->a[k][l]);
         }
        }
       }

 for (i3 = 0; i3 < 3; i3++) {
 dim  = i3 * atoms->number_of_sh_bfns_in_unit_cell;
 dim5 = i3 * nbands;
 DoubleGEMM3(&NoTrans,&Trans,&nbfn,&nbands,&nbfn,&alpha,&(M_k->a[dim]),&nbfn,&eigvec->a[0],&nbfn,&beta,&(tmp->a[0]),&nbands);
 DoubleGEMM3(&NoTrans,&NoTrans,&nbands,&nbands,&nbfn,&alpha,&eigvec->a[0],&nbfn,&(tmp->a[0]),&nbands,&beta,&(M_x->a[dim5]),&nbands);
 }
  //fprintf(file.out,"M_x\n");
  //print_real_matrix2(M_x,0,8,1.0,file);
  for (j3 = 0; j3 < job->field_dirs; j3++) {
    for (t = 0; t < ntransitions; t++) {
      free_particle_matrix_element[j3][t] = k_zero;
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
          for (t = 0; t < ntransitions; t++) {
          bse_matrix_element[j3][t] += free_particle_matrix_element[j3][count] * bse_eigenvectors->a[t][count];
         }
        }
        count++; // counter for all transitions
       } // close loop on n
      } // close loop on m
      for (j3 = 0; j3 < job->field_dirs; j3++) {
        for (t = 0; t < ntransitions; t++) {
          fragment_DM[i][j3][t] = bse_matrix_element[j3][t];
          //fragment_DM[i][j3][t] = free_particle_matrix_element[j3][t];
         }
        }
       } // close loop on i

  DestroyDoubleMatrix(&tmp,job);
  DestroyDoubleMatrix(&M_x,job);
  DestroyDoubleMatrix(&eigvec,job);
  DestroyComplexMatrix(&scf_eigenvectors,job);

 }

  if (job->taskid == 0 && job->verbosity >= 1) {
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|                        GW BSE TRANSITION DIPOLE MOMENTS FOR MOLECULAR FRAGMENTS                         |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"| DIPOLE DIRECTION COSINES  | DIR 1 %5.2f %5.2f %5.2f | DIR 2 %5.2f %5.2f %5.2f | DIR 3 %5.2f %5.2f %5.2f |\n",\
  job->e_field[0].comp1,job->e_field[0].comp2,job->e_field[0].comp3,\
  job->e_field[1].comp1,job->e_field[1].comp2,job->e_field[1].comp3,\
  job->e_field[2].comp1,job->e_field[2].comp2,job->e_field[2].comp3);
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  fprintf(file.out,"|       TRANSITION NO.      |        FRAGMENT 1       |        FRAGMENT 2       |                         |\n");
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  for (t = 0; t < ntransitions; t++) {
  fprintf(file.out,"|         %4d %10.4lf   | %7.4f %7.4f %7.4f | %7.4f %7.4f %7.4f |                         |\n",\
  t+1,bse_eigenvalues[t]*au_to_eV,\
      (fragment_DM[0][0][t] + nuclear_dipole[0][0]) * bohr_to_AA,(fragment_DM[0][1][t] + nuclear_dipole[0][1]) * bohr_to_AA,\
      (fragment_DM[0][2][t] + nuclear_dipole[0][2]) * bohr_to_AA,(fragment_DM[1][0][t] + nuclear_dipole[1][0]) * bohr_to_AA,\
      (fragment_DM[1][1][t] + nuclear_dipole[1][1]) * bohr_to_AA,(fragment_DM[1][2][t] + nuclear_dipole[1][2]) * bohr_to_AA);
 }
  fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 }

  DestroyDoubleMatrix(&M_k,job);

}
