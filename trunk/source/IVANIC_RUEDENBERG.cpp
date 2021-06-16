/*
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <fstream>
#include <mpi.h>
#include "mycomplex.h"
#include "mylogical.h"
#include "conversion_factors.h"
#include "LIMITS.h"
#include "CRYSTAL1.h"
*/
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "PRINT_UTIL.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "SYMMETRY_ADAPTATION.h"
#include "IVANIC_RUEDENBERG.h"

using namespace std;

void generate_rotation_operators_ivanic_ruedenberg(SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

int i, j, k, l, m, n, count, memsize;
int *p_ind_i, *p_ind_j, *p_num, *p_ost;
double u, v, w;
double R[2 * job->l_max + 1][2 * job->l_max + 1], Z[5][5];
double *p_irr, *p_rot;

//printf("lmax %3d\n", job->l_max);

  memsize = 9;
  for (i = 2; i <= job->l_max; i++) {
  memsize += (2 * i + 1) * (2 * i + 1);
//printf("memsize %3d %3d\n",i,memsize);
 }
  double S[memsize];

  count = 0;
  p_num = symmetry->num_ij;
  p_ost = symmetry->op_shift; 
  p_ind_i = symmetry->ind_i; 
  p_ind_j = symmetry->ind_j;
  p_rot = symmetry->rot;

  for (k = 0; k < symmetry->number_of_operators; k++) {

   *p_num = 1; // s shell
   *p_ost = count;
   *p_ind_i = 0;
   *p_ind_j = 0;
   *p_rot = k_one;
    //fprintf(file.out,"S %3d %3d\n",*symmetry->op_shift,count);
    //fprintf(file.out,"s rrk %3d %3d %3d %3d %3d\n",k,*p_num,*p_ind_i,*p_ind_j,*p_ost) ;
    p_ind_i++;
    p_ind_j++;
    p_rot++;
    p_num++;
    p_ost++;
    count++;

   *p_num = 0; // s part of sp shell
   *p_ind_i = 0;
   *p_ind_j = 0;
   *p_rot = k_one;
    p_irr = symmetry->irr + k * 9;
    (*p_num)++;
    //fprintf(file.out,"rr %3d %3d %3d\n",*p_ind_i,*p_ind_j,*p_rot) ;
    p_ind_i++;
    p_ind_j++;
    p_rot++;

    for (i = 1; i < 4; i++) { // p part of sp shell
      for (j = 1; j < 4; j++) {
        //fprintf(file.out,"sp p_irr %3d %3d %lf %3d\n",i,j,*p_irr,*p_num) ;
        if (fabs(*(p_irr)) > 0.00001) {
          *p_ind_i = i;
          *p_ind_j = j;
          *p_rot = *p_irr;
          //fprintf(file.out,"rr %3d %3d %3d\n",*p_ind_i,*p_ind_j,*p_rot) ;
          p_ind_i++;
          p_ind_j++;
          p_rot++;
          (*p_num)++;
        } // end if (fabs
        p_irr++;
      }
    }
   *p_ost = count;
    //fprintf(file.out,"SP %3d %3d\n",*symmetry->op_shift,count);
    count += *p_num;
    p_num++;
    p_ost++;

   *p_num = 0;
    p_irr = symmetry->irr + k * 9;
    for (i = 0; i < 3; i++) { // p shell
      for (j = 0; j < 3; j++) {
        //fprintf(file.out,"p p_irr %3d %3d %lf %3d\n",i,j,*p_irr,*p_num) ;
        if (fabs(*(p_irr)) > 0.00001) {
          *p_ind_i = i;
          *p_ind_j = j;
          *p_rot = *p_irr;
          //fprintf(file.out,"rr %3d %3d %10.4lf\n",*p_ind_i,*p_ind_j,*p_rot) ;
          p_ind_i++;
          p_ind_j++;
          p_rot++;
          (*p_num)++;
        } // end if (fabs
        p_irr++;
      }
    }
   *p_ost = count;
    //fprintf(file.out,"P %3d %3d\n",*symmetry->op_shift,count);
    count += *p_num;
    p_num++;
    p_ost++;

    for (i = 0; i < memsize; i++)
    S[i] = k_zero;

    p_irr = symmetry->irr + k * 9;

    S[0] = *(p_irr + 0);  // order is m,m' = 11, 10, 1-1, 01, 00, 0-1, 1-1, -10, -1-1
    S[1] = *(p_irr + 2);
    S[2] = *(p_irr + 1);
    S[3] = *(p_irr + 6);
    S[4] = *(p_irr + 8);
    S[5] = *(p_irr + 7);
    S[6] = *(p_irr + 3);
    S[7] = *(p_irr + 5);
    S[8] = *(p_irr + 4);

    for (l = 2; l <= job->l_max; l++) {

    for (m = -l; m <= l; m++) {
      for (n = -l; n <= l; n++) {
        uvw(l, m, n, &u, &v, &w);
        //printf("lmn %3d %3d %3d U\n",l,m,n);
        u *= U(l, m, n, S);
        //printf("lmn %3d %3d %3d V\n",l,m,n);
        v *= V(l, m, n, S);
        //printf("lmn %3d %3d %3d W\n",l,m,n);
        w *= W(l, m, n, S);
        S[M_offset(l, m, n)] = u + v + w;
        R[l + m][l + n] = u + v + w;
        //fprintf(file.out,"l m n %3d %3d %3d   l+m l+n %3d %3d  offset %3d %10.4lf\n", l,m,n, l+m, l+n, M_offset(l, m, n),u+v+w);
       }
      }

    if (l == 2) { // transform operators for l = 2 to Crystal order z^2, xz, yz, x^2 - y^2, xy
    for (i = 0; i < 5; i++) {
      for (j = 0; j < 5; j++) {
        Z[i][j] = R[i][j];
       }
      }
       crystal_sh_order(Z);
    for (i = 0; i < 5; i++) {
      for (j = 0; j < 5; j++) {
        R[i][j] = Z[i][j];
       }
      }
     }

   *p_num = 0;
    for (j = 0; j < 2 * l + 1; j++) { 
      for (i = 0; i < 2 * l + 1; i++) {
        //fprintf(file.out,"d/f/g p_irr %3d %3d %lf %3d\n",i,j,*p_irr,*p_num) ;
        if (fabs(R[i][j]) > 0.00001) {
          *p_ind_i = i;
          *p_ind_j = j;
          *p_rot = R[i][j];
          //fprintf(file.out,"rr %3d %3d %lf\n",*p_ind_i,*p_ind_j,*p_rot) ;
          p_ind_i++;
          p_ind_j++;
          p_rot++;
          (*p_num)++;
        } // end if (fabs
      }
    }
    //fprintf(file.out,"\n");
   *p_ost = count;
    //fprintf(file.out,"DFG %3d %3d\n",*symmetry->op_shift,count);
    //fprintf(file.out,"l %3d count %3d *p_num %3d\n",l,count,*p_num) ;
    count += *p_num;
    p_num++;
    p_ost++;

    } // close loop on l

  } // close loop on k

  if (job->taskid == 0 && job->verbosity > 1) {
  p_ind_i = symmetry->ind_i;
  p_ind_j = symmetry->ind_j;
  p_rot = symmetry->rot;
  p_ost = symmetry->op_shift;
  p_num = symmetry->num_ij;
  for (k = 0; k < symmetry->number_of_operators; k++) {
    fprintf(file.out, "k = symmetry->%3d\n", k);
    for (i = 0; i < job->l_max + 2; i++) {
      for (j = 0; j < *p_num; j++) {
        fprintf(file.out, "symm %3d ind_i %3d %3d num_ij %3d %3d rot  %10.4lf\n", i, *p_ind_i, *p_ind_j, *p_num,
            *p_ost, *p_rot);
        p_ind_i++;
        p_ind_j++;
        p_rot++;
      }
      fprintf(file.out, "\n");
      p_ost++;
      p_num++;
    }
    fprintf(file.out, "\n");
    fflush(file.out);
  } // close loop on k
 }

}

void crystal_sh_order(double Z[5][5])

{
int i,j, p, q;
double R[5][5], T[5][5];

  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      T[i][j] = k_zero;
      R[i][j] = k_zero;
     }
    }

    T[4][0] = k_one;
    T[2][1] = k_one;
    T[0][2] = k_one;
    T[1][3] = k_one;
    T[3][4] = k_one;

    for (p = 0; p < 5; p++) { 
      for (q = 0; q < 5; q++) { 
        for (i = 0; i < 5; i++) { 
          for (j = 0; j < 5; j++) {
            R[p][q] += T[p][i] * Z[i][j] * T[q][j];
           }
          }
         //printf("%10.4lf ",R[p][q]);
        }
       //printf("\n");
      }
     //printf("\n");

  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      Z[i][j] = R[i][j];
     }
    }

}

void uvw(const int l, const int m, const int n, double *u, double *v, double *w)

{

    double d = delta(m, 0);
    int abs_m = abs(m);
    double denom;

    if (abs(n) == l)
    denom = (2 * l) * (2 * l - 1);

    else
    denom = (l + n) * (l - n);

    *u =  sqrt((l + m) * (l - m) / denom);
    *v =  0.5 * sqrt((1 + d) * (l + abs_m - 1) * (l + abs_m) / denom) * (1 - 2 * d);
    *w = -0.5 * sqrt((l - abs_m - 1) * (l - abs_m) / denom) * (1 - d);

}

double U(const int l, const int m, const int n, double *R)

{

    return (P(0, l, m, n, R));

}

double V(const int l, const int m, const int n, double *R)

{
    if (m == 0) {
        double p0 = P(1, l, 1, n, R);
        double p1 = P(-1, l, -1, n, R);
        return (p0 + p1);
       }

    else if (m > 0) {
        double d = delta(m, 1);
        double p0 = P(1, l, m - 1, n, R);
        double p1 = P(-1, l, -m + 1, n, R);
        return (p0 * sqrt(1 + d) - p1 * (1 - d));
       }

    else  {  // m < 0
        double d = delta(m, -1);
        double p0 = P(1, l, m + 1, n, R);
        double p1 = P(-1, l, -m - 1, n, R);
        return (p0 * (1 - d) + p1 * sqrt(1 + d));
       }

}

double W(const int l, const int m, const int n, double *R)

{

    if (m == 0) {
      return (0.0);
     }

    else if (m > 0) {
      double p0 = P(1, l, m + 1, n, R);
      double p1 = P(-1, l, -m - 1, n, R);
      return (p0 + p1);
     }

    else {  // m < 0
      double p0 = P(1, l, m - 1, n, R);
      double p1 = P(-1, l, -m + 1, n, R);
      return (p0 - p1);
     }

}

double M(const int l, const int m, const int n, double *S)

{

   //  returns m,n element of real spherical harmonic Y^l(m,n)
   //  order of (m,n) indices is l,..,2,1,0,-1,-2,..,-l

int i;
int offset, row_index, col_index;

  offset = 0;
  for (i = 1; i < l; i++)
  offset += (2 * i + 1) * (2 * i + 1);

  row_index = (l - m);
  col_index = (l - n);

  if (abs(m) > l || abs(n) > l) return k_zero;
  else {
  //printf("M offset %3d lmn %3d %3d %3d row col %3d %3d tot %3d %10.4lf\n",offset,l,m,n,row_index,col_index,\
  offset + row_index * (2 * l + 1) + col_index,S[offset + row_index * (2 * l + 1) + col_index]);
  return S[offset + row_index * (2 * l + 1) + col_index];
 }

}

int M_offset(const int l, const int m, const int n)

{

     //returns memory offset of m,n element of real spherical harmonic Y^l(m,n)
     //order of (m,n) indices is l,..,2,1,0,-1,-2,..,-l

int i;
int offset, row_index, col_index;

  offset = 0;
  for (i = 1; i < l; i++)
  offset += (2 * i + 1) * (2 * i + 1);

  row_index = (l - m);
  col_index = (l - n);

  //printf("M offset %3d lmn %3d %3d %3d row col %3d %3d tot %3d\n",offset,l,m,n,row_index,col_index,\
  offset + row_index * (2 * l + 1) + col_index);

  return offset + row_index * (2 * l + 1) + col_index;

}

double P(const int i, const int l, const int a, const int b, double *R)

{

double rip1 = M(1, i,  1, R);
double rim1 = M(1, i, -1, R);
double ri00 = M(1, i,  0, R);

    if (b == -l)
      return (rip1 * M(l - 1, a, -l + 1, R) + rim1 * M(l - 1, a, l - 1, R));

    else if (b == l)
      return (rip1 * M(l - 1, a, l - 1, R) - rim1 * M(l - 1, a, -l + 1, R));

    else
      return (ri00 * M(l - 1, a, b, R));

}

void R_offset(const int l, double *R, double *S)

{

int i, m, n;
int offset, row_index, col_index;

  offset = 0;
  for (i = 1; i < l; i++)
  offset += (2 * i + 1) * (2 * i + 1);

  for (m = -l; m <= l; m++) {
    for (n = -l; n <= l; n++) {
      if (l == 2) {
        if (m ==  0) row_index = 0; // 3z^2 - r^2
        if (m ==  1) row_index = 1; //  xz
        if (m == -1) row_index = 2; //  yz
        if (m ==  2) row_index = 3; //  x^2 - y^2
        if (m == -2) row_index = 4; //  xy
        if (n ==  0) col_index = 0; // 3z^2 - r^2
        if (n ==  1) col_index = 1; //  xz
        if (n == -1) col_index = 2; //  yz
        if (n ==  2) col_index = 3; //  x^2 - y^2
        if (n == -2) col_index = 4; //  xy
       }
        R[offset + row_index * (2 * l + 1) + col_index] = S[M_offset(l, m, n)];
        //printf("%3d %3d %3d %3d %3d %10.4lf %10.4lf\n", \
        l,m,n,offset + row_index * (2 * l + 1) + col_index, M_offset(l, m, n),R[offset + row_index * (2 * l + 1) + col_index], S[M_offset(l, m, n)]);
      }
     }

}

double delta(const int m, const int n)

{
    // Kronecker Delta
    return (m == n ? 1 : 0);
}

void test_rotation_operators(ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // generates a full operator from a reduced operator using permutation and rotation operations

  int i, j, l;
  int op;
  int oppshift1, oppshift2;
  int sheli1, shelj1;
  int index_i, shelposi;
  int *p_i, *p_i1, *p_j, *p_j1;
  double *p_rot1, *p_rot2, *p_F_unit_matrix, *F_rotate, *p_F_rotate, *F_unit_matrix;;

  //for (l = 1; l < 6; l+=2) {

/*
  for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
    shelposi = atoms->shelposn[i];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel[i]; index_i++) {
    fprintf(file.out,"atom %3d shell %3d order %3d\n",i,index_i,shells->ord[index_i]);
    for (op = 0; op < symmetry->number_of_operators; op++) {
    oppshift1 = *(symmetry->op_shift + op * 5 + shells->ord[index_i]);
    sheli1 = *(symmetry->num_ij + op * 5 + shells->ord[index_i]);
    p_i1 = symmetry->ind_i + oppshift1;
    p_i = symmetry->ind_j + oppshift1;
    p_rot1 = symmetry->rot + oppshift1;
    for (j = 0; j < sheli1; j++) {
    fprintf(file.out,"op %3d op_shift %5d i %3d   i1,i  %3d %3d  rot1 %10.4lf\n",op,oppshift1,j,*p_i1,*p_i,*p_rot1);
   }
  }
  fprintf(file.out,"\n");
 }
}
*/

//CHANGES2015 l = 5;
  l = 2 * job->l_max + 1;

  F_rotate = (double *) malloc(l * l * sizeof(double));
  if (F_rotate == NULL) {  fprintf(stderr, "ERROR: not enough memory for double F_rotate! \n"); exit(1); }

  F_unit_matrix = (double *) malloc(l * l * sizeof(double));
  if (F_unit_matrix == NULL) {  fprintf(stderr, "ERROR: not enough memory for double F_unit_matrix! \n"); exit(1); }

    for (j = 0; j < l * l; j++)
    F_unit_matrix[j] = k_zero;
 
    for (j = 0; j < l; j++) 
    F_unit_matrix[j * l + j] = k_one;
    //F_unit_matrix[0] = k_one;
      
    for (op = 0; op < symmetry->number_of_operators; op++) {

    for (j = 0; j < l * l; j++) 
    F_rotate[j] = k_zero;
 
        //oppshift1 = *(symmetry->op_shift + op * 5 + shells->ord[index_i]);
        //oppshift2 = *(symmetry->op_shift + op * 5 + shells->ord[index_j]);
        //sheli1 = *(symmetry->num_ij + op * 5 + shells->ord[index_i]);
        //shelj1 = *(symmetry->num_ij + op * 5 + shells->ord[index_j]);
        //CHANGES2015oppshift1 = *(symmetry->op_shift + op * 5 + 3);
        //oppshift2 = *(symmetry->op_shift + op * 5 + 3);
        //sheli1 = *(symmetry->num_ij + op * 5 + 3);
        //CHANGES2015shelj1 = *(symmetry->num_ij + op * 5 + 3);
        oppshift1 = *(symmetry->op_shift + op * (job->l_max + 2) + job->l_max + 1);
        oppshift2 = *(symmetry->op_shift + op * (job->l_max + 2) + job->l_max + 1);
        sheli1 = *(symmetry->num_ij + op * (job->l_max + 2) + job->l_max + 1);
        shelj1 = *(symmetry->num_ij + op * (job->l_max + 2) + job->l_max + 1);
        //fprintf(file.out,"oppshift %3d %3d %3d %3d %3d %3d\n",l,op,oppshift1,oppshift2,sheli1,shelj1) ;

        p_i1 = symmetry->ind_i + oppshift1;
        p_i = symmetry->ind_j + oppshift1;
        p_rot1 = symmetry->rot + oppshift1;
        for (i = 0; i < sheli1; i++) {
          p_j1 = symmetry->ind_i + oppshift2;
          p_j = symmetry->ind_j + oppshift2;
          p_rot2 = symmetry->rot + oppshift2;
          for (j = 0; j < shelj1; j++) {
            //p_F_unit_matrix = F_unit_matrix + *p_i1 * l + *p_j1;
            //p_F_rotate      = F_rotate      + *p_i  * l + *p_j;
            p_F_unit_matrix = F_unit_matrix + *p_i * l + *p_j;
            p_F_rotate      = F_rotate      + *p_i1  * l + *p_j1;
            //fprintf(file.out,"*p_i,j  %d %d %d %d\n",*p_i1, *p_j1, *p_i, *p_j);
           *p_F_rotate += *p_rot1 * *p_rot2 * *p_F_unit_matrix;
            p_j++;
            p_j1++;
            p_rot2++;
          }
          p_i++;
          p_i1++;
          p_rot1++;
        }

  if (job->taskid ==0 && job->verbosity > 1) {
  fprintf(file.out,"Initial Operator %3d L = %3d\n",op,l);
  for (i = 0; i < l; i++) {
    for (j = 0; j < l; j++) {
      fprintf(file.out,"%16.10lf ",F_unit_matrix[i * l + j]);
     }
      fprintf(file.out,"\n");
    }

  fprintf(file.out,"Operator %3d L = %3d\n",op,l);
  for (i = 0; i < l; i++) {
    for (j = 0; j < l; j++) {
      fprintf(file.out,"%16.10lf ",F_rotate[i * l + j]);
     }
      fprintf(file.out,"\n");
    }
      fprintf(file.out,"\n");
  }
  
   } // close loop on op

  free(F_rotate);
  free(F_unit_matrix);

}

