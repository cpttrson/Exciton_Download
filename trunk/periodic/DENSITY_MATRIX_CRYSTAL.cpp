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

