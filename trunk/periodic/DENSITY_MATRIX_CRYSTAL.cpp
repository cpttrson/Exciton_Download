#include <cstring>
#include "mycomplex.h"
#include "conversion_factors.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "ROTATIONS_MOLECULE.h"
#include "PRINT_UTIL.h"
#include "MATRIX_UTIL.h"
#include "DENSITY_MATRIX_MOLECULE.h"
#include "DENSITY_MATRIX_CRYSTAL.h"

using namespace std;

void initial_density_matrix_crystal3(double *P0, double *F0, PAIR_TRAN *pair_p, FERMI *fermi, ATOM *atoms,  ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file){

  int i, j, k, l, m, p, q, r, s, u;
  int dim, dim1, dimk, dimp, dims;
  int bfposi;
  int nd1, nd2;
  int vector_size, offset, nkunique;
  int count;
  int num_sym = symmetry->number_of_operators;
  int taken[atoms->number_of_basis_sets];
  int occupied[job->spin_dim];
  double fac;
  double electron_count_up, electron_count_down, electron_occupation_up, electron_occupation_down;
  double total_electron_occupation_up, total_electron_occupation_down;
  double time1, time2;
  double *eigenvalues;
  double *occupation;
  DoubleMatrix *eigenvectors, *P;

  num_sym = 1;
  nkunique = fermi->nkunique;
  dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  dimk = job->spin_dim * nkunique * atoms->number_of_sh_bfns_in_unit_cell;
  vector_size = dimk * dim1;
  time1 = MPI_Wtime();
  for (i = 0; i < atoms->number_of_basis_sets; i++) taken[i] = -1;
  total_electron_occupation_up = k_zero;
  total_electron_occupation_down = k_zero;
  job->electron_count = k_zero;

  // ******************************************************************************************
  // * Loop over all basis sets and unique atoms                                              *
  // ******************************************************************************************

  for (i = 0; i < atoms->number_of_basis_sets; i++) {
    for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
      if (atoms->basis_set[j] == i) {
        if (taken[i] == -1) {
          electron_occupation_up = k_zero;
          electron_occupation_down = k_zero;
          electron_count_up = k_zero;
          electron_count_down = k_zero;
          dim  = atoms->bfnnumb_sh[j];
          dimp = atoms->magnetic[j] * atoms->bfnnumb_sh[j];
          dims = job->spin_dim * atoms->bfnnumb_sh[j];
          AllocateDoubleMatrix(&P,&dimp,&dim,job);
          AllocateDoubleMatrix(&eigenvectors,&dimp,&dim,job);
          AllocateDoubleArray(&eigenvalues,&dims,job);
          AllocateDoubleArray(&occupation,&dims,job);
          ResetDoubleMatrix(P);
          ResetDoubleMatrix(eigenvectors);
          ResetDoubleArray(eigenvalues,&dims);
          ResetDoubleArray(occupation,&dims);
          atom_scf(atoms,j,eigenvectors,eigenvalues,occupation,occupied,shells,gaussians,job,file);
          if (atoms->magnetic[j] == 1)
            for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
              electron_occupation_up   += occupation[m];
             }
          if (atoms->magnetic[j] == 2 && atoms->spin[j] > 0)
          //if (atoms->magnetic[j] == 2 && atoms->spin[j] == 1)
            for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
              electron_occupation_up   += occupation[m];
              electron_occupation_down += occupation[atoms->bfnnumb_sh[j] + m];
             }
          if (atoms->magnetic[j] == 2 && atoms->spin[j] < 0)
          //if (atoms->magnetic[j] == 2 && atoms->spin[j] == -1)
            for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
              electron_occupation_up   += occupation[atoms->bfnnumb_sh[j] + m];
              electron_occupation_down += occupation[m];
             }
          for (s = 0; s < atoms->magnetic[j]; s++) {
            for (m = 0; m < occupied[s]; m++) {
              if (s == 0) electron_count_up   += k_one;
              if (s == 1) electron_count_down += k_one;
             }
            }
          //for (m = 0; m < dims; m++) fprintf(file.out,"%3d %10.4f\n",m,occupation[m]);
        fprintf(file.out,"%d %d Occup/dn %3d %3d occup/dn %10.4f %10.4f totup/dn %10.4f %10.4f ctup/dn %10.4f %10.4f job->spin_dim %3d\n",\
        atoms->magnetic[j],atoms->spin[j],occupied[0],occupied[1],electron_occupation_up,electron_occupation_down,electron_count_up,\
        electron_count_down,job->spin_dim);

  // ******************************************************************************************
  // * Calculate density matrix for unique atom j                                             *
  // ******************************************************************************************

          for (s = 0; s < atoms->magnetic[j]; s++) {
            offset = s * atoms->bfnnumb_sh[j];
            fac = k_one;
            if (job->spin_dim == 1) fac = electron_occupation_up / electron_count_up;
            //if (s == 0 && occupied[s] > 0) fac = electron_occupation_up   / electron_count_up;
            //if (s == 1 && occupied[s] > 0) fac = electron_occupation_down / electron_count_down;
            //fprintf(file.out,"fac %f %f %f %f\n",electron_occupation_up,electron_count_up,electron_occupation_down,electron_count_down);
            for (k = 0; k < dim; k++) {
              for (l = 0; l < dim; l++) {
                for (m = 0; m < occupied[s]; m++) {
                  P->a[offset + k][l] += fac * eigenvectors->a[m][k] * eigenvectors->a[m][l];
                 }
                }
               }
              }
          if (job->taskid == 0 && job->verbosity > 1) {
            fprintf(file.out,"eigvec2 %3d %3d\n",dims,dim);
            if (job->spin_polarisation == 0)
            print_real_eigenvector_matrix(eigenvectors, eigenvalues, file);
            if (job->spin_polarisation == 1)
            PrintDoubleMatrix(eigenvectors, file);
            fprintf(file.out,"Atom density matrix\n");
            for (k = 0; k < P->iRows; k++) {
            for (l = 0; l < P->iCols; l++) {
            fprintf(file.out,"%6.3f ",P->a[k][l]); }
            fprintf(file.out,"\n"); }
            fflush(file.out);
           }
          taken[i] = 0;
         } // close if (taken[i]

  // ******************************************************************************************
  // * Loop over atoms symmetry equivalent to unique atom j                                   *
  // * Assemble reduced and full density matrices P0 and F0                                   *
  // ******************************************************************************************

  //for (s = 0; s < job->spin_dim; s++) {
  //for (s = 0; s < atoms->magnetic[j]; s++) {
    //if (atoms->spin[j] == -1) u = 1 - s;
    //if (atoms->spin[j] ==  0) u = 0;
    //if (atoms->spin[j] ==  1) u = s;
    //offset = s * atoms->bfnnumb_sh[j];
    //fprintf(file.out,"j %d s %d u %d\n",j,s,u);

    bfposi = 0;
    for (p = 0; p < pair_p->nump; p++) {
      q = pair_p->posn[p];
      nd1 = atoms->bfnnumb_sh[pair_p->cell1[q]];
      nd2 = atoms->bfnnumb_sh[pair_p->cell2[q]];
      if (pair_p->cell1[q] == j && pair_p->cell2[q] == j && pair_p->latt1[q] == 0 && pair_p->latt2[q] == 0) {
        if (job->spin_dim == 1) {
        count = 0;
        for (k = 0; k < dim; k++) {
          for (l = 0; l < dim; l++) {
            //P0[s * job->dimp + bfposi + count] = P->a[offset + k][l];
            P0[bfposi + count] = P->a[k][l];
            F0[bfposi + count] = P->a[k][l];
            count++;
           }
          }
         }
        if (job->spin_dim == 2 && atoms->magnetic[j] == 1) {
        count = 0;
        for (k = 0; k < dim; k++) {
          for (l = 0; l < dim; l++) {
            P0[bfposi + count] = P->a[k][l];
            P0[bfposi + count + job->dimp] = P->a[k][l];
            F0[bfposi + count] = P->a[k][l];
            F0[bfposi + count + job->dimf] = P->a[k][l];
            count++;
           }
          }
         }
        if (job->spin_dim == 2 && atoms->magnetic[j] == 2) {
        count = 0;
        for (k = 0; k < dim; k++) {
          for (l = 0; l < dim; l++) {
            if (atoms->spin[j] > 0) {
            //if (atoms->spin[j] == 1) {
            P0[bfposi + count] = P->a[k][l];
            P0[bfposi + count + job->dimp] = P->a[atoms->bfnnumb_sh[j] + k][l];
            F0[bfposi + count] = P->a[k][l];
            F0[bfposi + count + job->dimf] = P->a[atoms->bfnnumb_sh[j] + k][l];
           }
            else if (atoms->spin[j] < 0) {
            //else if (atoms->spin[j] == -1) {
            P0[bfposi + count] = P->a[atoms->bfnnumb_sh[j] + k][l];
            P0[bfposi + count + job->dimp] = P->a[k][l];
            F0[bfposi + count] = P->a[atoms->bfnnumb_sh[j] + k][l];
            F0[bfposi + count + job->dimf] = P->a[k][l];
           }
            count++;
           }
          }
         }
        }
       bfposi += nd1 * nd2;
      }

    //bfposi = 0;
    //for (p = 0; p < pair_p->tot; p++) {
      //nd1 = atoms->bfnnumb_sh[pair_p->cell1[p]];
      //nd2 = atoms->bfnnumb_sh[pair_p->cell2[p]];
      //if (pair_p->cell1[p] == j && pair_p->cell2[p] == j && pair_p->latt1[p] == 0 && pair_p->latt2[p] == 0) {
        //count = 0;
        //for (k = 0; k < dim; k++) {
          //for (l = 0; l < dim; l++) {
            //F0[s * job->dimf + bfposi + count] = P->a[offset + k][l];
            //count++;
           //}
          //}
         //}
        //bfposi += nd1 * nd2;
       //}
      //} // close loop on s

  // ******************************************************************************************
  // * Count total number of spin up and spin down electrons                                  *
  // ******************************************************************************************

  if (atoms->magnetic[j] == 1)
  for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
    total_electron_occupation_up   += occupation[m] / two;
    total_electron_occupation_down += occupation[m] / two;
    job->electron_count += occupation[m];
   }
  if (atoms->magnetic[j] == 2 && atoms->spin[j] > 0)
  //if (atoms->magnetic[j] == 2 && atoms->spin[j] == 1)
  for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
    total_electron_occupation_up   += occupation[m];
    total_electron_occupation_down += occupation[atoms->bfnnumb_sh[j] + m];
    job->electron_count += occupation[m] + occupation[atoms->bfnnumb_sh[j] + m];
   }
  else if (atoms->magnetic[j] == 2 && atoms->spin[j] < 0)
  //else if (atoms->magnetic[j] == 2 && atoms->spin[j] == -1)
  for (m = 0; m < atoms->bfnnumb_sh[j]; m++) {
    total_electron_occupation_up   += occupation[atoms->bfnnumb_sh[j] + m];
    total_electron_occupation_down += occupation[m];
    job->electron_count += occupation[m] + occupation[atoms->bfnnumb_sh[j] + m];
   }
    //fprintf(file.out,"electron count %e occupationup/dn %10.4f %10.4f\n",\
    job->electron_count,total_electron_occupation_up,total_electron_occupation_down);
   } // close if (atoms
  } // close loop on j
   if (taken[i] == 0) {
     DestroyDoubleMatrix(&eigenvectors,job);
     DestroyDoubleArray(&eigenvalues,&dims,job);
     DestroyDoubleArray(&occupation,&dims,job);
     DestroyDoubleMatrix(&P,job);
     taken[i]++;
    }
   } // close loop on i

  // ******************************************************************************************
  // * Generate fermi->occupied and fermi->occupation arrays                                  *
  // ******************************************************************************************

  for (s = 0; s < job->spin_dim; s++) {
    for (k = 0; k < nkunique; k++) {
      fermi->occupied[s * nkunique + k] = 0;
      for (i = 0; i < fermi->occupied[s]; i++) {
      fermi->occupation[(s * nkunique  + k) * dim1 + i] = k_zero;
     }
    }
   }

  switch (crystal->type[0]) {

    case 'C':
    case 'S':
    case 'P':

    for (s = 0; s < job->spin_dim; s++) {
      if      (s == 0) fermi->occupied[s] = (int)(total_electron_occupation_up   + 0.00001);
      else if (s == 1) fermi->occupied[s] = (int)(total_electron_occupation_down + 0.00001);
      for (k = 0; k < nkunique; k++) {
        for (i = 0; i < fermi->occupied[s]; i++) {
          fermi->occupation[(s * nkunique  + k) * dim1 + i] = k_one;
          //fprintf(file.out,"spin %3d occupied %10.4f\n",s,fermi->occupation[(s * nkunique  + k) * dim1 + i]);
         }
        }
       }

    for (s = 0; s < job->spin_dim; s++) {
      for (k = 0; k < nkunique; k++) {
        //if      (s == 0) fermi->occupied[k]            = 13;
        //else if (s == 1) fermi->occupied[nkunique + k] = 11;
        if      (s == 0) fermi->occupied[k]            = (int)(total_electron_occupation_up   + 0.00001);
        else if (s == 1) fermi->occupied[nkunique + k] = (int)(total_electron_occupation_down + 0.00001);
        //fprintf(file.out,"spin %3d k %3d occupied %3d\n",s,k,fermi->occupied[s * nkunique + k]);
         }
        }

   break;

    case 'M':

    if (job->spin_dim == 1) fermi->occupied[0] = (int)(total_electron_occupation_up + 0.00001);
    else if (job->spin_dim == 2) {
      fermi->occupied[0] = (int)(total_electron_occupation_up   + 0.00001);
      fermi->occupied[1] = (int)(total_electron_occupation_down + 0.00001);
      for (s = 0; s < job->spin_dim; s++) {
        for (i = 0; i < fermi->occupied[s]; i++) {
          fermi->occupation[s * dim1 + i] = k_one;
          //fprintf(file.out,"spin %3d occupied %3d occupation %lf\n",s,fermi->occupied[s],fermi->occupation[s * dim1 + i]);
         }
        }
       }

     break;
  
  } // close switch

  for (i = 0; i < atoms->number_of_basis_sets; i++) {
      if (taken[i] == -1) {
        if (job->taskid == 0 || job->verbosity > 1)
        fprintf(file.out,"Initial guess for atom %d failed in initial_density_matrix_crystal\n",j);
        printf("Initial guess for atom %d failed in initial_density_matrix_crystal\n",j);
        MPI_Finalize();
        exit(0);
      }
     }

  time2 = MPI_Wtime();

}

void density_matrix_crystal2(FERMI *fermi, double *P0, double *F0, KPOINT_TRAN *knet, char *filename, int *nkunique, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

  {

    int i, j, k, p, q, s;
    int dimp_spin, dimf_spin;
    int dim1, dim2, count;

    dimp_spin = job->dimp * job->spin_dim;
    dimf_spin = job->dimf * job->spin_dim;

    if (job->guess_type == 0 && job->density == 2)   // only zero P0 if beginning from atomic wave functions
      for (i=0;i<dimp_spin;i++) { P0[i] = k_zero;}
      for (i=0;i<dimf_spin;i++) { F0[i] = k_zero;}


    if (job->guess_type == 0 && job->density == 2) { 
    reduced_density_matrix_crystal(fermi,P0,knet,filename,R,R_tables,pair_p,atom_p,atoms,shells,symmetry,job,file);
   }

    expand_density_matrix(P0,F0,pair_p,atoms,shells,symmetry,job,file);

    if (job->taskid == 0 && job->verbosity > 1) {
      fprintf(file.out,"density matrix %d\n\n",pair_p->nump);
      count = 0;
      for (s = 0; s < job->spin_dim; s++) {
        for (p = 0; p < pair_p->nump; p++) {
          q = pair_p->posn[p];
            fprintf(file.out,"unique pair %d spin %d\n",p,s);
              for (k = 0; k < pair_p->numb[p]; k++) {
              fprintf(file.out,"%3d   [%d] [%d]   [%d]\n",q+k,pair_p->cell1[q+k],pair_p->cell2[q+k],pair_p->latt2[q+k]);
          for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
            for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
              //fprintf(file.out,"%16.12lf ",F0[count]);
              //fprintf(file.out,"%10.3e ",F0[count]);
              fprintf(file.out,"%6.3f ",F0[count]);
              count++;
              }
            fprintf(file.out,"\n");
           }
            fprintf(file.out,"\n");
           }
          fprintf(file.out,"\n\n");
         }
        fprintf(file.out,"\n\n");
       }
      }

}

void reduced_density_matrix_crystal(FERMI *fermi, double *P, KPOINT_TRAN *knet, char *filename, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shell, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, gj, k, kp, l, m, n, p, q, s;
int dim, dim1, dim2, dimk;
int bfposi, bfposj;
int vector_size;
int count, count1, countk, kunique, dimp_spin = job->dimp * job->spin_dim;
int nbands, bands[4];
double temp;
Complex temp1;
char buf2[110];
ComplexMatrix *eigenvectors, *eigenvectors1;
MPI_File fh;

  strcpy(buf2,file.scf_eigvec);
  strcat(buf2,filename);

  dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  if (job->fix_occ == 0) {
  dim2 = dim1;
 }
  else if (job->fix_occ == 1) {
  s = 0; // need to fix this
  bands[0] = fermi->bands[2 * s] - 1; 
  bands[1] = fermi->bands[2 * s + 1] - 1;
  nbands = bands[1] - bands[0] + 1;
  dim2 = nbands;
 }
  dimk = job->spin_dim * fermi->nkunique * dim2;
  vector_size = dimk * dim1;
  AllocateComplexMatrix(&eigenvectors,&dimk,&dim1,job);
  AllocateComplexMatrix(&eigenvectors1,&dimk,&dim1,job);

  MPI_File_open(MPI_COMM_WORLD,buf2,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) ;
  MPI_File_seek(fh, 0, MPI_SEEK_SET) ;
  MPI_File_read(fh, &eigenvectors->a[0][0], 2 * vector_size, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
  MPI_File_close(&fh);

  count1 = 0;
  for (s = 0; s < job->spin_dim; s++) {
    countk = 0;
    for (k = 0; k < fermi->nkunique; k++) {
      if (job->kss == 0) kunique = countk;
      else if (job->kss == 1) kunique = k;
      for (kp = 0; kp < knet->num[k]; kp++) {
        if (job->fix_occ == 0) {
        bands[0] = 0; 
        bands[1] = fermi->occupied[s * fermi->nkunique + k];
        nbands = bands[1] - bands[0];
        dim2 = dim1;
       }
        else if (job->fix_occ == 1) {
        bands[0] = fermi->bands[2 * s] - 1; 
        bands[1] = fermi->bands[2 * s + 1] - 1;
        nbands = bands[1] - bands[0] + 1;
        dim2 = nbands;
       }
        rotate_psi(&eigenvectors->a[(s * fermi->nkunique + kunique) * dim2][0],&eigenvectors1->a[0][0],nbands,knet->bz[countk],knet,\
        atom_p,atoms,R,shell,symmetry,job,file);
        dim = s * job->dimp;
        for (p = 0; p < pair_p->nump; p++) {
          q = pair_p->posn[p];
          gj = pair_p->latt2[q]; 
          temp = double_vec_dot(&knet->cart[knet->bz[countk]],&R->vec_ai[gj]) ;
          temp1 = (double)job->spin_fac * Complex(cos(temp),sin(temp)) ;
          bfposi = 0;
          bfposj = 0;
          for (l = 0; l < pair_p->cell1[q]; l++)
          bfposi += atoms->bfnnumb_sh[l];
          for (l = 0; l < pair_p->cell2[q]; l++)
          bfposj += atoms->bfnnumb_sh[l];
           for (n = 0; n < nbands; n++) {
            count = 0;
            for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
              for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
                P[dim + count] += (conj(eigenvectors1->a[n][bfposi + i]) * eigenvectors1->a[n][bfposj + j] * \
                fermi->occupation[count1 + n] / (double)fermi->nktot * temp1).real();
                count++;
               }
              }
             }
            dim += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
           } // close loop on p
          countk++;
         } // close loop on kp
        count1 += dim1;
       } // close loop on k
      } // close loop on s

  DestroyComplexMatrix(&eigenvectors,job);
  DestroyComplexMatrix(&eigenvectors1,job);

  if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"reduced density matrix\n");
     count = 0;
     for (s = 0; s < job->spin_dim; s++) {
       for (p = 0; p < pair_p->nump; p++) {
         q = pair_p->posn[p];
         fprintf(file.out,"pair %d spin %d   %3d %3d  %4d\n",p,s,pair_p->cell1[q],pair_p->cell2[q],pair_p->latt2[q]);
         for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
           for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
             fprintf(file.out,"%7.3lf",P[count]);
             count++;
            }
           fprintf(file.out,"\n");
          }
         fprintf(file.out,"\n");
        }
       fprintf(file.out,"\n");
      }
       fflush(file.out);
     }

}

void reduced_density_matrix_crystal3(FERMI *fermi, double *P, KPOINT_TRAN *knet, int s, int k, ComplexMatrix *eigenvectors, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, PAIR_TRAN *pair_p, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shell, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, gj, kp, l, n, p, q;
int dim, dim1;
int bfposi, bfposj;
int count, count1, countk, dimp_spin = job->dimp * job->spin_dim;
int nbands;
double temp;
Complex temp1;
ComplexMatrix *eigenvectors1;

  dim1 = atoms->number_of_sh_bfns_in_unit_cell;
  AllocateComplexMatrix(&eigenvectors1,&dim1,&dim1,job);

  countk = 0;
  for (kp = 0; kp < k; kp++) 
  countk += knet->num[kp];
  count1 = (s * fermi->nkunique + k ) * dim1;
  nbands = fermi->bands[1] - fermi->bands[0] + 1;
  for (kp = 0; kp < knet->num[k]; kp++) {
  //fprintf(file.out,"s %d k %d kp %d countk %d count1 %d nbands %d \n",s,k,kp,countk,count1,nbands);
      rotate_psi(&eigenvectors->a[0][0],&eigenvectors1->a[0][0],nbands,knet->bz[countk],knet,atom_p,atoms,R,shell,symmetry,job,file);
      dim = s * job->dimp;
      for (p = 0; p < pair_p->nump; p++) {
        q = pair_p->posn[p];
        gj = pair_p->latt2[q]; 
        temp = double_vec_dot(&knet->cart[knet->bz[countk]],&R->vec_ai[gj]) ;
        temp1 = (double)job->spin_fac * Complex(cos(temp),sin(temp)) ;
        bfposi = 0;
        bfposj = 0;
        for (l = 0; l < pair_p->cell1[q]; l++)
        bfposi += atoms->bfnnumb_sh[l];
        for (l = 0; l < pair_p->cell2[q]; l++)
        bfposj += atoms->bfnnumb_sh[l];
         for (n = 0; n < nbands; n++) {
          count = 0;
          for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
            for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
              P[dim + count] += (conj(eigenvectors1->a[n][bfposi + i]) * eigenvectors1->a[n][bfposj + j] * fermi->occupation[count1 + n] \
              / (double)fermi->nktot * temp1).real();
              count++;
             }
            }
           }
          dim += atoms->bfnnumb_sh[pair_p->cell1[q]] * atoms->bfnnumb_sh[pair_p->cell2[q]];
         } // close loop on p
       countk++;
       } // close loop on kp

  DestroyComplexMatrix(&eigenvectors1,job);

  if (job->taskid == 0 && job->verbosity > 1) {
     fprintf(file.out,"MPP reduced density matrix\n");
     count = 0;
     for (s = 0; s < job->spin_dim; s++) {
       for (p = 0; p < pair_p->nump; p++) {
         q = pair_p->posn[p];
         fprintf(file.out,"pair %d spin %d   %3d %3d  %4d\n",p,s,pair_p->cell1[q],pair_p->cell2[q],pair_p->latt2[q]);
         for (i = 0; i < atoms->bfnnumb_sh[pair_p->cell1[q]]; i++) {
           for (j = 0; j < atoms->bfnnumb_sh[pair_p->cell2[q]]; j++) {
             fprintf(file.out,"%7.3lf",P[count]);
             count++;
            }
           fprintf(file.out,"\n");
          }
         fprintf(file.out,"\n");
        }
       fprintf(file.out,"\n");
      }
       fflush(file.out);
     }

}

