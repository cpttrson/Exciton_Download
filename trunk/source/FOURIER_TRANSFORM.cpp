#include "mycomplex.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "PRINT_UTIL.h"
#include "MATRIX_UTIL.h"
#include "ROTATIONS_MOLECULE.h"
#include "SYMMETRY_ADAPTATION.h"
#include "ALLOCATE_MEMORY.h"
#include "FOURIER_TRANSFORM.h"

using namespace std;

void fourier_transform(double *matrix_reduced, Complex *matrix_k, KPOINT_TRAN *knet, int nk[2], PAIR_TRAN *pair_p, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int i, j, k, l, m, p, pm, k1, dim2, dim3, dim4 = atoms->number_of_sh_bfns_in_unit_cell * atoms->number_of_sh_bfns_in_unit_cell;
  int atm1, atm2, nd1, nd2, nkpt = nk[1] - nk[0] + 1;
  double *p_matrix_reduced, *p_matrix, *matrix, temp1;
  double time1 = MPI_Wtime();
  Complex *p_matrix_k, temp;

  dim2 = 0;
  for (i = 0; i < pair_p->nump; i++) {
    atm1 = pair_p->cell1[pair_p->posn[i]];
    atm2 = pair_p->cell2[pair_p->posn[i]];
    nd1 = atoms->bfnnumb_sh[atm1];
    nd2 = atoms->bfnnumb_sh[atm2];
    dim2 += nd1 * nd2 * pair_p->numb[i];
  }

  matrix = (double *) malloc(dim2 * sizeof(double));
  if (matrix == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double matrix! \n");
    exit(1);
  }

  p_matrix = matrix;
  for (i = 0; i < dim2; i++) {
    *p_matrix = k_zero;
    p_matrix++;
  }

  p_matrix_k = matrix_k;
  for (i = 0; i < dim4 * nkpt; i++) {
    *p_matrix_k = k_zero;
    p_matrix_k++;
  }

  dim3 = 0;
  for (p = 0; p < pair_p->nump; p++) {
    i = pair_p->posn[p];
    nd1 = atoms->bfnnumb_sh[pair_p->cell1[i]];
    nd2 = atoms->bfnnumb_sh[pair_p->cell2[i]];
    rotate_permute_expand_pair(p, pair_p, &matrix_reduced[dim3], matrix, atoms, shells, symmetry, job, file);
    dim3 += nd1 * nd2;
    k1 = 0;
    for (k = nk[0]; k <= nk[1]; k++) {
      for (j = 0; j < pair_p->numb[p]; j++) {
        nd1 = atoms->bfnnumb_sh[pair_p->cell1[i + j]];
        nd2 = atoms->bfnnumb_sh[pair_p->cell2[i + j]];
        temp1 = -double_vec_dot(&knet->cart[k], &R->vec_ai[pair_p->latt2[i + j]]);
        temp = Complex(cos(temp1), sin(temp1));
        p_matrix   = matrix + nd1 * nd2 * j;
        p_matrix_k = matrix_k + atoms->bfnposn_sh[pair_p->cell1[i + j]] * \
        atoms->number_of_sh_bfns_in_unit_cell + atoms->bfnposn_sh[pair_p->cell2[i + j]] + dim4 * k1;
        for (l = 0; l < nd1; l++) {
          for (m = 0; m < nd2; m++) {
            *p_matrix_k += *p_matrix * temp;
             p_matrix++;
             p_matrix_k++;
            }
             p_matrix_k += atoms->number_of_sh_bfns_in_unit_cell - nd2;
            } // close loop on l
           } // close loop on j
          k1++;
         } // close loop on k
        } // close loop on p

  free(matrix);

  double time2 = MPI_Wtime();
  if (job->verbosity > 1)
  fprintf(file.out, "Time for Fourier Transform %10.4e\n", (double) time2 - time1);

}

void fourier_transform_3(double *matrix_reduced, Complex *matrix_k, KPOINT_TRAN *knet, int nk[2], PAIR_TRAN *pair_p, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{
      int i, j, k, l, m, n, p, dim2, dim3, dim4 = atoms->number_of_sh_bfns_in_unit_cell * atoms->number_of_sh_bfns_in_unit_cell;
      int atm1, atm2, nd1, nd2, k1, pm, nkpt = nk[1] - nk[0] + 1;
      int lim1, lim2;
      double *p_matrix_reduced, *p_matrix, *matrix, temp1 ;
      Complex *p_matrix_k, temp ;
      double time1 = MPI_Wtime();

  // identical to fourier_transform_3 except that it uses rotate_permute_expand_tensor1 instead of new_rotate_tensor1_pair

      dim2 = 0 ;

      for (i = 0; i < pair_p->nump; i++) {
      atm1 = pair_p->cell1[pair_p->posn[i]];
      atm2 = pair_p->cell2[pair_p->posn[i]];
      nd1 = atoms->bfnnumb_sh[atm1];
      nd2 = atoms->bfnnumb_sh[atm2];
      dim2 += nd1 * nd2 * pair_p->numb[i];
     }

      matrix = (double *) malloc(dim2 * 3 * sizeof(double)) ;
      if (matrix==NULL) { fprintf(stderr, "ERROR: not enough memory for double matrix! \n") ; exit(1) ; }

      p_matrix = matrix ;
      for (i = 0; i < dim2 * 3; i++) {
      *p_matrix = k_zero ; p_matrix++ ; }

      p_matrix_k = matrix_k ;
      for ( i = 0; i < dim4 * 3 * nkpt; i++) {
      *p_matrix_k = k_zero ;
       p_matrix_k++ ; }

      dim3 = 0 ;

      for ( p = 0; p < pair_p->nump; p++) {
      i = pair_p->posn[p];
      nd1 = atoms->bfnnumb_sh[pair_p->cell1[i]];
      nd2 = atoms->bfnnumb_sh[pair_p->cell2[i]];
      rotate_permute_expand_tensor1(p, pair_p, &matrix_reduced[3 * dim3], matrix, atoms, shells, symmetry, job, file);

      dim3 += nd1 * nd2 ;

      k1 = 0;
      for (k = nk[0]; k <= nk[1]; k++) {
      for (j = 0; j < pair_p->numb[p]; j++) {
        nd1 = atoms->bfnnumb_sh[pair_p->cell1[i + j]];
        nd2 = atoms->bfnnumb_sh[pair_p->cell2[i + j]];
      //temp1 = -double_vec_dot(&knet->cart[k], &R->vec_ai[pair_p->latt2[i + j]]);
      //temp1 = two;
      temp1 = double_vec_dot(&knet->cart[k], &R->vec_ai[pair_p->latt2[i + j]]);
      // CHECK temp1 = double_vec_dot(&knet->cart[knet->ibz[k]], &R->vec_ai[pair_p->latt2[i + j]]);
      if (pair_p->p[i + j] == 0) 
      temp  =  Complex(sin(temp1),-cos(temp1)) ; // Multiply in -i factor for -i del here
      if (pair_p->p[i + j] == 1) 
      temp  =  Complex(-sin(temp1), cos(temp1)) ; // Multiply in -i factor for -i del here change of sign by hermiticity
      for (n = 0; n < 3; n++) {
      p_matrix      = matrix      + nd1 * nd2 * (3 * j + n) ;
      p_matrix_k = matrix_k + atoms->bfnposn_sh[pair_p->cell1[i + j]] * \
      atoms->number_of_sh_bfns_in_unit_cell + atoms->bfnposn_sh[pair_p->cell2[i + j]] + dim4 * (3 * k1 + n);
      //fprintf(file.out,"%d\n",n);
      //fprintf(file.out,"M_k address %4d %4d %4d %4d\n",atoms->bfnposn_sh[pair_p->cell1[i + j]] * \
      atoms->number_of_sh_bfns_in_unit_cell , atoms->bfnposn_sh[pair_p->cell2[i + j]] , dim4 * (3 * k1 + n),\
      atoms->number_of_sh_bfns_in_unit_cell - nd2);
      for (l = 0; l < nd1; l++) {
      for (m = 0; m < nd2; m++) {
      *p_matrix_k += *p_matrix * temp ;
      //fprintf(file.out,"matrix temp %4d %4d %12.4e %12.4e %12.4e  %12.4e %12.4e\n", \
      l,m,*p_matrix,temp.real(),temp.imag(),(*p_matrix_k).real(),(*p_matrix_k).imag());
      //fprintf(file.out,"%10.4lf",*p_matrix);
       p_matrix++ ;
       p_matrix_k++ ;
      }
      //fprintf(file.out,"\n");
       p_matrix_k += atoms->number_of_sh_bfns_in_unit_cell - nd2;
      } // close loop on l

      } // close loop on n
      } // close loop on j
      k1++;
      } // close loop on k
      } // close loop on p

      free(matrix) ;

}

void fourier_transform_molecule(double *matrix_reduced, double *matrix_k, PAIR_TRAN *pair_p, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  int i, j, l, m, p, k1, dim2, dim3, dim4 = atoms->number_of_sh_bfns_in_unit_cell * atoms->number_of_sh_bfns_in_unit_cell;
  int atm1, atm2, nd1, nd2; //, nkpt = nk[1] - nk[0] + 1;
  double *p_matrix_reduced, *p_matrix, *matrix; //, temp1;
  double *p_matrix_k; //, temp;
  double time1 = MPI_Wtime();

  dim2 = 0;
  for (i = 0; i < pair_p->nump; i++) {
    atm1 = pair_p->cell1[pair_p->posn[i]];
    atm2 = pair_p->cell2[pair_p->posn[i]];
    nd1 = atoms->bfnnumb_sh[atm1];
    nd2 = atoms->bfnnumb_sh[atm2];
    dim2 += nd1 * nd2 * pair_p->numb[i];
  }

  matrix = (double *) malloc(dim2 * sizeof(double));
  if (matrix == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double matrix! \n");
    exit(1);
  }

  for (i = 0; i < dim2; i++) matrix[i] = k_zero;
  for (i = 0; i < dim4; i++) matrix_k[i] = k_zero;
  dim3 = 0;
  for (p = 0; p < pair_p->nump; p++) {
    i = pair_p->posn[p];
    nd1 = atoms->bfnnumb_sh[pair_p->cell1[i]];
    nd2 = atoms->bfnnumb_sh[pair_p->cell2[i]];
    rotate_permute_expand_pair(p, pair_p, &matrix_reduced[dim3], matrix, atoms, shells, symmetry, job, file);
    dim3 += nd1 * nd2;
    k1 = 0;
    for (j = 0; j < pair_p->numb[p]; j++) {
      nd1 = atoms->bfnnumb_sh[pair_p->cell1[i + j]];
      nd2 = atoms->bfnnumb_sh[pair_p->cell2[i + j]];
      p_matrix   = matrix + nd1 * nd2 * j;
      p_matrix_k = matrix_k + atoms->bfnposn_sh[pair_p->cell1[i + j]] * \
      atoms->number_of_sh_bfns_in_unit_cell + atoms->bfnposn_sh[pair_p->cell2[i + j]] + dim4 * k1;
      for (l = 0; l < nd1; l++) {
        for (m = 0; m < nd2; m++) {
          *p_matrix_k += *p_matrix;
           p_matrix++;
           p_matrix_k++;
          }
           p_matrix_k += atoms->number_of_sh_bfns_in_unit_cell - nd2;
          } // close loop on l
         } // close loop on j
        k1++;
       } // close loop on p
  free(matrix);
  double time2 = MPI_Wtime();
  if (job->verbosity > 1)
  fprintf(file.out, "Time for Fourier Transform %10.4e\n", (double) time2 - time1);

}

void rotate_permute_expand_tensor_3(double *matrix_reduced, double *matrix_k, PAIR_TRAN *pair_p, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{
      int i, j, l, m, n, p, dim2, dim3, dim4 = atoms->number_of_sh_bfns_in_unit_cell * atoms->number_of_sh_bfns_in_unit_cell;
      int atm1, atm2, nd1, nd2, pm;
      double *p_matrix_reduced, *p_matrix, *matrix, *p_matrix_k ;
      double time1 = MPI_Wtime();

  // Routine similar to fourier_transform_3 which expands, permutes and rotates pairs for finite systems

      dim2 = 0 ;

      for (i = 0; i < pair_p->nump; i++) {
      atm1 = pair_p->cell1[pair_p->posn[i]];
      atm2 = pair_p->cell2[pair_p->posn[i]];
      nd1 = atoms->bfnnumb_sh[atm1];
      nd2 = atoms->bfnnumb_sh[atm2];
      dim2 += nd1 * nd2 * pair_p->numb[i];
     }

      matrix = (double *) malloc(dim2 * 3 * sizeof(double)) ;
      if (matrix==NULL) { fprintf(stderr, "ERROR: not enough memory for double matrix! \n") ; exit(1) ; }

      p_matrix = matrix ;
      for (i = 0; i < dim2 * 3; i++) {
      *p_matrix = k_zero ; p_matrix++ ; }

      p_matrix_k = matrix_k ;
      for ( i = 0; i < dim4 * 3; i++) {
      *p_matrix_k = k_zero ;
       p_matrix_k++ ; }

      dim3 = 0 ;

      for ( p = 0; p < pair_p->nump; p++) {
      i = pair_p->posn[p];
      nd1 = atoms->bfnnumb_sh[pair_p->cell1[i]];
      nd2 = atoms->bfnnumb_sh[pair_p->cell2[i]];
      rotate_permute_expand_tensor1(p, pair_p, &matrix_reduced[3 * dim3], matrix, atoms, shells, symmetry, job, file);

      dim3 += nd1 * nd2 ;

      //fprintf(file.out,"pair_p %d   [%3d][%3d]  [%3d]\n",i,pair_p->cell1[i],pair_p->cell2[i],pair_p->latt2[i]);
      for (j = 0; j < pair_p->numb[p]; j++) {
      //fprintf(file.out,"pair_p %d   [%3d][%3d]  [%3d]\n",i,pair_p->cell1[i + j],pair_p->cell2[i + j],pair_p->latt2[i + j]);
      nd1 = atoms->bfnnumb_sh[pair_p->cell1[i + j]];
      nd2 = atoms->bfnnumb_sh[pair_p->cell2[i + j]];
      for (n = 0; n < 3; n++) {
      p_matrix      = matrix      + nd1 * nd2 * (3 * j + n) ;
      p_matrix_k = matrix_k + atoms->bfnposn_sh[pair_p->cell1[i + j]] * \
      atoms->number_of_sh_bfns_in_unit_cell + atoms->bfnposn_sh[pair_p->cell2[i + j]] + dim4 * n;
      //fprintf(file.out,"%d\n",n);
      //fprintf(file.out,"M_k address %4d %4d %4d %4d\n",atoms->bfnposn_sh[pair_p->cell1[i + j]] * \
      atoms->number_of_sh_bfns_in_unit_cell , atoms->bfnposn_sh[pair_p->cell2[i + j]] , dim4 * n,\
      atoms->number_of_sh_bfns_in_unit_cell - nd2);
      for (l = 0; l < nd1; l++) {
      for (m = 0; m < nd2; m++) {
      *p_matrix_k += *p_matrix; 
       p_matrix++ ;
       p_matrix_k++ ;
      }
       p_matrix_k += atoms->number_of_sh_bfns_in_unit_cell - nd2;
      } // close loop on l
      } // close loop on n
      } // close loop on j
      } // close loop on p

      free(matrix) ;

}

void fourier_transform_salc1(double *matrix_reduced, ComplexMatrix **S_k, KPOINT_TRAN *knet, int nk[2], ATOM_TRAN *atom_p, PAIR_TRAN *pair_p, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)
{

int i5;
int i, iii, j, k, r, s, s1, t, u, v, p, pp;
int dim2, dim3;
int atm, atm1, atm2, nd1, nd2, nkpt = nk[1] - nk[0] + 1;
int count_i, count_j;
int atm10, atm20, atm11, atm21;
int offset1, offset2;
int num_irrep_in_basis[symmetry->number_of_classes];
double *p_matrix_reduced, *p_matrix, *matrix, temp1;
double time_salc = k_zero, time_mult = k_zero;
double time_salc1, time_mult1;
Complex temp;
SALC salc1, salc2;
DoubleMatrix *salc_rows, *salc_cols;
IntMatrix *num_irrep_in_atom, *salc_count1;

  AllocateIntMatrix(&num_irrep_in_atom, &symmetry->number_of_classes, &atoms->number_of_unique_atoms, job);
  AllocateIntMatrix(&salc_count1, &symmetry->number_of_classes, &atoms->number_of_unique_atoms, job);
  ResetIntMatrix(salc_count1);
  count_atom_irrep(num_irrep_in_atom,atoms,atom_p,shells,symmetry,job,file);
  count_basis_irrep(num_irrep_in_basis,atoms,atom_p,shells,symmetry,job,file);
 
  for(i = 0; i < symmetry->number_of_classes; i++) {
  if (num_irrep_in_basis[i] == 0) continue;
  ResetComplexMatrix(S_k[i]);
 }

  dim2 = 0;
  for (i = 0; i < pair_p->nump; i++) {
    atm1 = pair_p->cell1[pair_p->posn[i]];
    atm2 = pair_p->cell2[pair_p->posn[i]];
    nd1 = atoms->bfnnumb_sh[atm1];
    nd2 = atoms->bfnnumb_sh[atm2];
    dim2 += nd1 * nd2 * pair_p->numb[i];
  }

  matrix = (double *) malloc(dim2 * sizeof(double));
  if (matrix == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double matrix! \n");
    exit(1);
  }

  p_matrix = matrix;
  for (i = 0; i < dim2; i++) {
    *p_matrix = k_zero;
     p_matrix++;
  }

//int salc_count[symmetry->number_of_classes][atoms->number_of_unique_atoms];
  //for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
  //for (s = 0; s < symmetry->number_of_classes; s++) {
  //salc_count[s][atm] = 0;
 //}
 //}
  for (atm = 1; atm < atoms->number_of_unique_atoms; atm++) {
  for (s = 0; s < symmetry->number_of_classes; s++) {
  //salc_count[s][atm] = salc_count[s][atm - 1] + num_irrep_in_atom->a[s][atm - 1] * symmetry->irp_dim_k[s];
  salc_count1->a[s][atm] = salc_count1->a[s][atm - 1] + num_irrep_in_atom->a[s][atm - 1] * symmetry->irp_dim_k[s];
 }
 }

  dim3 = 0;
  for (p = 0; p < pair_p->nump; p++) {
  i5 = pair_p->posn[p];
  nd1 = atoms->bfnnumb_sh[pair_p->cell1[i5]];
  nd2 = atoms->bfnnumb_sh[pair_p->cell2[i5]];
  atm1 = pair_p->cell1[i5];
  atm2 = pair_p->cell2[i5];
  atm10 = atoms->uniq[atm1];
  atm20 = atoms->uniq[atm2];
  atm11 = atom_p->posn[atm10];
  atm21 = atom_p->posn[atm20];
  rotate_permute_expand_pair(p, pair_p, &matrix_reduced[dim3], matrix, atoms, shells, symmetry, job, file);
  for (s = 0; s < symmetry->number_of_classes; s++) {
    offset1 = salc_count1->a[s][atm10];
    offset2 = salc_count1->a[s][atm20];
    count_atom_salc(s,atm1,&salc1,atoms,atom_p,pair_p,shells,symmetry,job,file);
    count_atom_salc(s,atm2,&salc2,atoms,atom_p,pair_p,shells,symmetry,job,file);
    //fprintf(file.out,"salc %3d %3d %3d   %3d %3d %3d\n",\
    s,atm1,salc_count1->a[s][atoms->uniq[atm1]],s,atm2,salc_count1->a[s][atoms->uniq[atm2]]);
    if (salc1.num_salc == 0 || salc2.num_salc == 0) continue;
    allocate_SALC(&salc1,symmetry,job,file);
    allocate_SALC(&salc2,symmetry,job,file);
    time_salc1 = MPI_Wtime();
    generate_atom_salc(s,atm1,&salc1,atoms,atom_p,pair_p,shells,symmetry,job,file);
    generate_atom_salc(s,atm2,&salc2,atoms,atom_p,pair_p,shells,symmetry,job,file);
    time_salc += MPI_Wtime() - time_salc1;
    AllocateDoubleMatrix(&salc_rows,&salc1.num_salc,&nd1,job);
    AllocateDoubleMatrix(&salc_cols,&nd2,&salc2.num_salc,job);
    ResetDoubleMatrix(salc_rows);
    ResetDoubleMatrix(salc_cols);
    for (pp = 0; pp < pair_p->numb[p]; pp++) {
      //temp1 = -double_vec_dot(&knet->cart[k], &R->vec_ai[pair_p->latt2[i5 + pp]]);
      temp1 = k_zero;
      temp = Complex(cos(temp1), sin(temp1));
      //fprintf(file.out,"PERM %3d %3d %3d   %3d %3d\n",p,pp,pair_p->p[i5 + pp],pair_p->cell1[i5 + pp],pair_p->cell2[i5 + pp]);

      if (pair_p->p[i5 + pp] == 0) {

      count_j = 0;
      for (r = 0; r < salc1.num_irp[s]; r++){
        for (t = 0; t < salc1.num_coef[r]; t++){
          salc_rows->a[r][salc1.bfn_posn->a[pair_p->cell1[i5 + pp] - atm11][count_j+t]] = \
          salc1.coeff->a[pair_p->cell1[i5 + pp] - atm11][count_j+t];
         }
        count_j += salc1.num_coef[r];
       }
    count_j = 0;
      for (r = 0; r < salc2.num_irp[s]; r++){
        for (t = 0; t < salc2.num_coef[r]; t++){
          salc_cols->a[salc2.bfn_posn->a[pair_p->cell2[i5 + pp] - atm21][count_j+t]][r] = \
          salc2.coeff->a[pair_p->cell2[i5 + pp]- atm21][count_j+t];
         }
        count_j += salc2.num_coef[r];
       }

      time_mult1 = MPI_Wtime();
      for(r=0;r<salc1.num_irp[s];r++){
        for(u=0;u<salc2.num_irp[s];u++){
          for(s1=0;s1<nd1;s1++){
            for(t=0;t<nd2;t++){
              S_k[s]->a[offset1 + r][offset2 + u] += \
              salc_rows->a[r][s1] * matrix[pp * nd1 * nd2 + s1 * nd2 + t] * salc_cols->a[t][u] * temp;
             }
            }
           }
          }
         time_mult += MPI_Wtime() - time_mult1;

        }

      if (pair_p->p[i5 + pp] == 1) {

      count_j = 0;
      for (r = 0; r < salc1.num_irp[s]; r++){
        for (t = 0; t < salc1.num_coef[r]; t++){
          salc_rows->a[r][salc1.bfn_posn->a[pair_p->cell2[i5 + pp] - atm11][count_j+t]] = \
          salc1.coeff->a[pair_p->cell2[i5 + pp] - atm11][count_j+t];
         }
        count_j += salc1.num_coef[r];
       }
    count_j = 0;
      for (r = 0; r < salc2.num_irp[s]; r++){
        for (t = 0; t < salc2.num_coef[r]; t++){
          salc_cols->a[salc2.bfn_posn->a[pair_p->cell1[i5 + pp] - atm21][count_j+t]][r] = \
          salc2.coeff->a[pair_p->cell1[i5 + pp]- atm21][count_j+t];
         }
        count_j += salc2.num_coef[r];
       }

      time_mult1 = MPI_Wtime();
      for(r=0;r<salc1.num_irp[s];r++){
        for(u=0;u<salc2.num_irp[s];u++){
          for(s1=0;s1<nd1;s1++){
            for(t=0;t<nd2;t++){
              S_k[s]->a[offset2 + u][offset1 + r] += \
              salc_cols->a[t][u] * matrix[pp * nd1 * nd2 + t * nd1 + s1] * salc_rows->a[r][s1] * temp;
             }
            }
           }
          }
         time_mult += MPI_Wtime() - time_mult1;

        }
         } // pp

   DestroyDoubleMatrix(&salc_rows,job);
   DestroyDoubleMatrix(&salc_cols,job);
   free_SALC(&salc2,job);
   free_SALC(&salc1,job);

  } // close loop on s

   dim3 += nd1 * nd2;

  } // close loop on p

  free(matrix);
  DestroyIntMatrix(&salc_count1,job);
  DestroyIntMatrix(&num_irrep_in_atom,job);
  //printf("salc %10.4f mult %10.4f \n",time_salc,time_mult);

}

void fourier_transform_salc_crystal(double *matrix_reduced, ComplexMatrix **S_k, KPOINT_TRAN *knet, int nk[2], ATOM_TRAN *atom_k, PAIR_TRAN *pair_p, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables_little_k, ATOM *atoms, SHELL *shells, SYMMETRY* symmetry_little_k_group, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)
{

int i5;
int i, iii, j, k, r, s, s1, t, u, v, p, pp;
int dim2, dim3;
int atm, atm1, atm2, nd1, nd2, nkpt = nk[1] - nk[0] + 1;
int count_i, count_j;
int atm10, atm20, atm11, atm21;
int offset1, offset2;
int num_irrep_in_basis[symmetry_little_k_group->number_of_classes];
double *p_matrix_reduced, *p_matrix, *matrix, temp1;
Complex temp;
SALC salc1, salc2;
DoubleMatrix *salc_rows, *salc_cols;
IntMatrix *num_irrep_in_atom, *salc_count1;

double time_salc = k_zero, time_mult = k_zero;
double time_salc1, time_mult1;

  AllocateIntMatrix(&num_irrep_in_atom, &symmetry_little_k_group->number_of_classes, &atoms->number_of_unique_atoms, job);
  AllocateIntMatrix(&salc_count1, &symmetry_little_k_group->number_of_classes, &atoms->number_of_unique_atoms, job);
  ResetIntMatrix(salc_count1);
  //count_atom_irrep(num_irrep_in_atom,atoms,atom_p,shells,symmetry,job,file);
  //count_basis_irrep(num_irrep_in_basis,atoms,atom_p,shells,symmetry,job,file);
  count_atom_irrep_crystal(num_irrep_in_atom,atoms,atom_k,R_tables_little_k,shells,symmetry_little_k_group,job,file);
  count_basis_irrep_crystal(num_irrep_in_basis,atoms,atom_k,R_tables_little_k,shells,symmetry_little_k_group,job,file);

  for(i = 0; i < symmetry_little_k_group->number_of_classes; i++) {
  if (num_irrep_in_basis[i] == 0) continue;
  ResetComplexMatrix(S_k[i]);
 }

  dim2 = 0;
  for (i = 0; i < pair_p->nump; i++) {
    atm1 = pair_p->cell1[pair_p->posn[i]];
    atm2 = pair_p->cell2[pair_p->posn[i]];
    nd1 = atoms->bfnnumb_sh[atm1];
    nd2 = atoms->bfnnumb_sh[atm2];
    dim2 += nd1 * nd2 * pair_p->numb[i];
  }

  matrix = (double *) malloc(dim2 * sizeof(double));
  if (matrix == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double matrix! \n");
    exit(1);
  }

  p_matrix = matrix;
  for (i = 0; i < dim2; i++) {
    *p_matrix = k_zero;
     p_matrix++;
  }

  for (atm = 1; atm < atoms->number_of_unique_atoms; atm++) {
  for (s = 0; s < symmetry_little_k_group->number_of_classes; s++) {
  salc_count1->a[s][atm] = salc_count1->a[s][atm - 1] + num_irrep_in_atom->a[s][atm - 1] * symmetry_little_k_group->irp_dim_k[s];
 }
 }

  dim3 = 0;
  for (p = 0; p < pair_p->nump; p++) {
  i5 = pair_p->posn[p];
  nd1 = atoms->bfnnumb_sh[pair_p->cell1[i5]];
  nd2 = atoms->bfnnumb_sh[pair_p->cell2[i5]];
  atm1 = pair_p->cell1[i5];
  atm2 = pair_p->cell2[i5];
  atm10 = atoms->uniq[atm1];
  atm20 = atoms->uniq[atm2];
  atm11 = atom_k->posn[atm10];
  atm21 = atom_k->posn[atm20];
  int ip, latt;
  int array[symmetry_little_k_group->number_of_operators];
  int oper[symmetry_little_k_group->number_of_operators];
  int lat_rot[symmetry_little_k_group->number_of_operators];
  IntMatrix *lattice_tran1, *lattice_tran2;
  ip = pair_p->cell1[i5];
  latt = pair_p->latt1[i5];
  count_atom_cluster(ip,latt,&salc1,array,oper,lat_rot,symmetry_little_k_group,atom_k,R_tables_little_k);
  fprintf(file.out,"count_atom1 %3d\n",salc1.num_atom);
  AllocateIntMatrix(&lattice_tran1, &salc1.num_atom, &symmetry_little_k_group->number_of_operators, job);
  generate_atom_cluster(lattice_tran1,&salc1,array,oper,lat_rot,symmetry_little_k_group);
  ip = pair_p->cell2[i5];
  latt = pair_p->latt2[i5];
  count_atom_cluster(ip,latt,&salc2,array,oper,lat_rot,symmetry_little_k_group,atom_k,R_tables_little_k);
  fprintf(file.out,"count_atom2 %3d\n",salc2.num_atom);
  AllocateIntMatrix(&lattice_tran2, &salc2.num_atom, &symmetry_little_k_group->number_of_operators, job);
  generate_atom_cluster(lattice_tran2,&salc2,array,oper,lat_rot,symmetry_little_k_group);
  rotate_permute_expand_pair(p, pair_p, &matrix_reduced[dim3], matrix, atoms, shells, symmetry, job, file);

/*
  for (pp = 0; pp < pair_p->numb[p]; pp++) {
  fprintf(file.out,"p %3d pp %3d cell2 %3d %10.4lf\n",p,pp,pair_p->latt2[i5 + pp],\
  cos(-double_vec_dot(&knet->cart[nk[0]], &R->vec_ai[pair_p->latt2[i5 + pp]])));
int count1 = 0;
  for (j = 0; j < 2; j++) {
  for (k = 0; k < 2; k++) {
  fprintf(file.out,"%16.10lf ",matrix_reduced[dim3 + count1]);
  count1++;
 }
  fprintf(file.out,"\n");
 }
  fprintf(file.out,"\n");
 }
*/

  for (s = 0; s < symmetry_little_k_group->number_of_classes; s++) {
    offset1 = salc_count1->a[s][atm10];
    offset2 = salc_count1->a[s][atm20];
    count_atom_salc_crystal(s,lattice_tran1,&salc1,atoms,atom_k,pair_p,R_tables_little_k,shells,symmetry_little_k_group,job,file);
    count_atom_salc_crystal(s,lattice_tran2,&salc2,atoms,atom_k,pair_p,R_tables_little_k,shells,symmetry_little_k_group,job,file);
    //fprintf(file.out,"salc %3d %3d %3d   %3d %3d %3d\n",\
    s,atm1,salc_count1->a[s][atoms->uniq[atm1]],s,atm2,salc_count1->a[s][atoms->uniq[atm2]]);
    if (salc1.num_salc == 0 || salc2.num_salc == 0) continue;
    allocate_SALC(&salc1,symmetry_little_k_group,job,file);
    allocate_SALC(&salc2,symmetry_little_k_group,job,file);
    time_salc1 = MPI_Wtime();
    generate_atom_salc_crystal(s,lattice_tran1,&salc1,atoms,atom_k,pair_p,R_tables_little_k,shells,symmetry_little_k_group,job,file);
    generate_atom_salc_crystal(s,lattice_tran2,&salc2,atoms,atom_k,pair_p,R_tables_little_k,shells,symmetry_little_k_group,job,file);

/*
  for (pp = 0; pp < pair_p->numb[p]; pp++) {
  fprintf(file.out,"p %3d pp %3d s %3d num_irp %3d %3d\n",p,pp,s,salc1.num_irp[s],salc2.num_irp[s]);
  count_j = 0;
  for (r = 0; r < salc2.num_irp[s]; r++) {
    fprintf(file.out,"%3d %3d ",s,r);
    for (t = 0; t < salc2.num_coef[r]; t++){
      fprintf(file.out,"%7.3f ",salc2.coeff->a[pair_p->cell2[i5 + pp] - atm21][count_j+t]);
     }
    count_j += salc2.num_coef[r];
    fprintf(file.out,"\n");
   }
   }
*/

    time_salc += MPI_Wtime() - time_salc1;
    AllocateDoubleMatrix(&salc_rows,&salc1.num_salc,&nd1,job);
    AllocateDoubleMatrix(&salc_cols,&nd2,&salc2.num_salc,job);
    ResetDoubleMatrix(salc_rows);
    ResetDoubleMatrix(salc_cols);
    for (pp = 0; pp < pair_p->numb[p]; pp++) {
      temp1 = -double_vec_dot(&knet->cart[nk[0]], &R->vec_ai[pair_p->latt2[i5 + pp]]);
      temp = Complex(cos(temp1), sin(temp1));

      if (pair_p->p[i5 + pp] == 0) {

      count_j = 0;
      for (r = 0; r < salc1.num_irp[s]; r++){
        for (t = 0; t < salc1.num_coef[r]; t++){
          salc_rows->a[r][salc1.bfn_posn->a[pair_p->cell1[i5 + pp] - atm11][count_j+t]] = \
          salc1.coeff->a[pair_p->cell1[i5 + pp] - atm11][count_j+t];
         }
        count_j += salc1.num_coef[r];
       }
    count_j = 0;
      for (r = 0; r < salc2.num_irp[s]; r++){
        for (t = 0; t < salc2.num_coef[r]; t++) {
          salc_cols->a[salc2.bfn_posn->a[pair_p->cell2[i5 + pp] - atm21][count_j + t]][r] = \
          salc2.coeff->a[pair_p->cell2[i5 + pp] - atm21][count_j + t];
         }
        count_j += salc2.num_coef[r];
       }

      time_mult1 = MPI_Wtime();
      for (r = 0; r < salc1.num_irp[s]; r++) {
        for (u = 0; u < salc2.num_irp[s]; u++) {
          for (s1 = 0; s1 < nd1; s1++) {
            for (t = 0; t < nd2; t++) {
              //fprintf(file.out,"s %3d offset1 %3d offset2 %3d r %3d s1 %3d t %d u %3d temp %10.4lf\n",s,offset1,offset2,r,s1,t,u,temp);
              //S_k[s]->a[offset1 + r][offset2 + u] += \
              salc_rows->a[r][s1] * matrix[pp * nd1 * nd2 + s1 * nd2 + t] * salc_cols->a[t][u] * temp;
              S_k[s]->a[offset1 + r][offset2 + u] += \
              salc_rows->a[r][s1] * matrix[pp * nd1 * nd2 + s1 * nd2 + t] * salc_cols->a[t][u] * temp;
             }
            }
           }
          }
         time_mult += MPI_Wtime() - time_mult1;

        }

      if (pair_p->p[i5 + pp] == 1) {

      count_j = 0;
      for (r = 0; r < salc1.num_irp[s]; r++) {
        for (t = 0; t < salc1.num_coef[r]; t++) {
          salc_rows->a[r][salc1.bfn_posn->a[pair_p->cell2[i5 + pp] - atm11][count_j + t]] = \
          salc1.coeff->a[pair_p->cell2[i5 + pp] - atm11][count_j + t];
         }
        count_j += salc1.num_coef[r];
       }
    count_j = 0;
      for (r = 0; r < salc2.num_irp[s]; r++){
        for (t = 0; t < salc2.num_coef[r]; t++){
          salc_cols->a[salc2.bfn_posn->a[pair_p->cell1[i5 + pp] - atm21][count_j + t]][r] = \
          salc2.coeff->a[pair_p->cell1[i5 + pp]- atm21][count_j + t];
         }
        count_j += salc2.num_coef[r];
       }

      time_mult1 = MPI_Wtime();
      for(r=0;r<salc1.num_irp[s];r++){
        for(u=0;u<salc2.num_irp[s];u++){
          for(s1=0;s1<nd1;s1++){
            for(t=0;t<nd2;t++){
              S_k[s]->a[offset2 + u][offset1 + r] += \
              salc_cols->a[t][u] * matrix[pp * nd1 * nd2 + t * nd1 + s1] * salc_rows->a[r][s1] * temp;
             }
            }
           }
          }
         time_mult += MPI_Wtime() - time_mult1;

        }
         } // pp

   DestroyDoubleMatrix(&salc_rows,job);
   DestroyDoubleMatrix(&salc_cols,job);
   free_SALC(&salc2,job);
   free_SALC(&salc1,job);

  } // close loop on s

   dim3 += nd1 * nd2;
   DestroyIntMatrix(&lattice_tran1,job);
   DestroyIntMatrix(&lattice_tran2,job);

  } // close loop on p

  free(matrix);
  DestroyIntMatrix(&salc_count1,job);
  DestroyIntMatrix(&num_irrep_in_atom,job);
  //printf("salc %10.4f mult %10.4f \n",time_salc,time_mult);

}

void count_atom_cluster(int ip, int latt, SALC *salc, int *array, int *oper, int *lat_rot, SYMMETRY *symmetry, ATOM_TRAN *atom_p, REAL_LATTICE_TABLES *R_tables)

{

int i, op;
int atm_rot, O_temp, latt_rot, count_atoms;

  count_atoms = 0;
  for (i = 0; i < symmetry->number_of_operators; i++) array[i] = -1;
  for (op = 0; op < symmetry->number_of_operators; op++) {
    atm_rot  = atom_p->K[ip * symmetry->number_of_operators + op];
    O_temp   = atom_p->O[ip * symmetry->number_of_operators + op];
    latt_rot = R_tables->diffvec[R_tables->lattvec[latt * symmetry->number_of_operators + op] * R_tables->margin_vector + O_temp];
    lat_rot[op] = R_tables->diffvec[R_tables->lattvec[latt * symmetry->number_of_operators + op] * R_tables->margin_vector + O_temp];
    for (i = 0; i <= count_atoms; i++) {
      if (latt_rot == array[i]) break;
      if (latt_rot != array[i] && i == count_atoms) { array[i] = latt_rot; oper[i] = op; count_atoms++; break; }
     }
    }
   salc->num_atom = count_atoms;

}

void generate_atom_cluster(IntMatrix *lattice_tran, SALC *salc, int *array, int *oper, int *lat_rot, SYMMETRY *symmetry)

{

int i, j, op;
int latt_tran[salc->num_atom][symmetry->number_of_operators];
int latt_rotd[salc->num_atom][symmetry->number_of_operators];

  for (i = 0; i < salc->num_atom; i++) {
    for (op = 0; op < symmetry->number_of_operators; op++) {
      latt_tran[i][op] = lat_rot[symmetry->grp_k[oper[i] * symmetry->number_of_operators + op]];
     }
    }
  for (i = 0; i < salc->num_atom; i++) {
    for (j = 0; j < salc->num_atom; j++) {
      for (op = 0; op < symmetry->number_of_operators; op++) {
       if (latt_tran[j][op] == array[i]) lattice_tran->a[j][op] = i; 
      }
     }
    }

}
