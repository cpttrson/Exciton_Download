

  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 Alin-Marin Elena                            *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

#include <cstdlib>
#include <cstring>
#include "myconstants.h"
#include "mycomplex.h"
#include "USER_DATA.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "MATRIX_UTIL.h"
#include "LINEAR_ALGEBRA_UTIL.h"

extern "C" {

void zdotc_(Complex *, const int *, Complex *, const int *, Complex *, const int *);

}

using namespace std;

void ComplexDotProd(Complex *c, const int *n, Complex **a, const int *inca, Complex **b, const int *incb)

{

#ifdef ISOF2003

  ComplexDotProdv1(c, n, a, inca, b, incb );

#else

  zdotc_(c,n, *a, inca, *b, incb );

#endif

}

void DoubleDotProd(double *c, const int *n, double **a, const int *inca, double **b, const int *incb)

{

#ifdef ISOF2003

// DoubleDotProdv1(c, n, a, inca, b, incb );

#else

  //ddot_(c, n, *a, inca, *b, incb);

  *c = ddot_(n, *a, inca, *b, incb);

  //printf("%d %d %d %lf %lf %lf\n",*n,*inca,*incb,*(*a),*(*b),*c);

#endif

}

void DoubleGEMM(const char *ta, const char *tb, const double *alpha, DoubleMatrix **ma,  DoubleMatrix **mb, const double *beta, DoubleMatrix **mc)

{
#ifndef ISOF2003

  int ione=1;
  int n,m,k, lda,ldb,ldc;

#endif

#ifdef ISOF2003

  DoubleGEMMv1(&ta,&tb, alpha, ma, mb,beta, mc);

#else

  if (*ta=='N' && *tb =='N'){
   m=(**mb).iCols;
   n= (**ma).iRows;
   k=(**mb).iRows;
   lda=m;
   ldb=k;
   ldc=m;
  }

  if (((*ta =='T')||(*ta=='C')) && (*tb =='N')){
   m=(**mb).iCols;
   n= (**ma).iCols;
   k=(**mb).iRows;
   lda=m;
   ldb=n;
   ldc=m;
  }

  if ((*ta=='N') && ((*tb =='T')||(*tb=='C'))){
   m=(**mb).iRows;
   n=(**ma).iRows;
   k=(**mb).iCols;
   lda=k;
   ldb=k;
   ldc=m;
  }

  dgemm_(tb, ta, &m, &n, &k, alpha, *((**mb).a), &lda, *((**ma).a), &ldb, beta, *((**mc).a), &ldc, &ione, &ione);

#endif

}

void DoubleGEMM3(const char *ta, const char *tb, int *m, int *n, int *k, const double *alpha, double **ma, int *lda,  double **mb, int *ldb, const double *beta, double **mc, int *ldc)

{

#ifndef ISOF2003

  int ione = 1;

#endif

#ifdef ISOF2003

  DoubleGEMM2v1(&ta, &tb, alpha, ma, mb, beta, mc);

#else

   dgemm_(tb, ta, n, m, k, alpha, *mb, ldb, *ma, lda, beta, *mc, ldc, &ione, &ione);

#endif

}

void ComplexGEMM2(const char *ta, const char *tb, int *m, int *n, int *k, const Complex *alpha, ComplexMatrix **ma, int *lda,  ComplexMatrix **mb, int *ldb, const Complex *beta, ComplexMatrix **mc, int *ldc)

{

#ifndef ISOF2003

  int ione = 1;

#endif

#ifdef ISOF2003

  ComplexGEMM2v1(&ta, &tb, alpha, ma, mb, beta, mc);

#else

//  if (*ta=='N' && *tb =='N'){
//   *m = (**mb).iCols;
//   *n = (**ma).iRows;
//   *k = (**mb).iRows;
//   *lda=*m;
//   *ldb=*k;
//   *ldc=*m;
// }
//
// if (*ta=='T' && *tb =='N'){
//  *m=(**mb).iCols;
//  *n= (**ma).iCols;
//  *k=(**mb).iRows;
//  *lda=*m;
//  *ldb=*n;
//  *ldc=*m;
// }
//
//if (*ta=='N' && *tb =='T'){
//  *m=(**mb).iRows;
//  *n=(**ma).iRows;
//  *k=(**mb).iCols;
//  *lda=*k;
//  *ldb=*k;
//  *ldc=*m;
// }

//printf("ta tb %c %c %d %d %d %d %d %d\n",*ta,*tb,*m,*n,*k,*lda,*ldb,*ldc);

   zgemm_(tb, ta, n, m, k, alpha, *((**mb).a), ldb, *((**ma).a), lda, beta, *((**mc).a), ldc, &ione, &ione);

#endif

}

void ComplexGEMM3(const char *ta, const char *tb, int *m, int *n, int *k, const Complex *alpha, Complex **ma, int *lda,  Complex **mb, int *ldb, const Complex *beta, Complex **mc, int *ldc)

{

#ifndef ISOF2003

  int ione = 1;

#endif

#ifdef ISOF2003

  ComplexGEMM2v1(&ta, &tb, alpha, ma, mb, beta, mc);

#else

   zgemm_(tb, ta, n, m, k, alpha, *mb, ldb, *ma, lda, beta, *mc, ldc, &ione, &ione);

#endif

}

void ComplexGEMM1(const char *ta, const char *tb, const Complex *alpha, ComplexMatrix **ma,  ComplexMatrix **mb, const Complex *beta, ComplexMatrix **mc)

{

#ifndef ISOF2003

  int ione=1;
  int n,m,k, lda,ldb,ldc;

#endif

#ifdef ISOF2003

  ComplexGEMM1v1(&ta,&tb, alpha, ma, mb,beta, mc);

#else

  if (*ta=='N' && *tb =='N'){
   m=(**mb).iCols;
   n= (**ma).iRows;
   k=(**mb).iRows;
   lda=m;
   ldb=k;
   ldc=m;
  }

  if (((*ta =='T')||(*ta=='C')) && (*tb =='N')){
   m=(**mb).iCols;
   n= (**ma).iCols;
   k=(**mb).iRows;
   lda=m;
   ldb=n;
   ldc=m;
  }

  if ((*ta=='N') && ((*tb =='T')||(*tb=='C'))){
   m=(**mb).iRows;
   n=(**ma).iRows;
   k=(**mb).iCols;
   lda=k;
   ldb=k;
   ldc=m;
  }

   zgemm_(tb, ta, &m, &n, &k, alpha, *((**mb).a), &lda, *((**ma).a), &ldb, beta, *((**mc).a), &ldc, &ione, &ione);

//printf("ta tb %d %d %d %d %d %d\n",m,n,k,lda,ldb,ldc);

#endif

}

void ComplexGEMM(const Complex *alpha, ComplexMatrix **ma, ComplexMatrix **mb, const Complex *beta, ComplexMatrix **mc)

{

  int ione=1;
  char ta='N',tb='N';

#ifdef ISOF2003

  ComplexGEMMv1(&ta,&tb, alpha, ma, mb,beta, mc);

#else

  zgemm_(&ta, &tb, &(**mb).iCols, &(**ma).iRows, &(**mb).iRows, alpha, *((**mb).a), &(**mb).iCols, *((**ma).a), &(**ma).iCols, beta, *((**mc).a), &(**mc).iCols,&ione,&ione);

#endif

}

void DoubleTRMM3(const char *side, const char *uplo, const char *ta, const char *diag, int *m, int *n, const double *alpha, double **ma, int *lda,  double **mb, int *ldb)

{

   dtrmm_(side, uplo, ta, diag, n, m, alpha, *mb, ldb, *ma, lda);

}

void SVD_DGESDD(DoubleMatrix **a, double **w, DoubleMatrix **eig_left, DoubleMatrix **eig_right, char *job, int *info)

{

 double tmp[1], *work;
 int lwork, liwork, itmp[1], *iwork;
 int ione = 1;

   tmp[0] = 1;
   itmp[0] = 1;
   lwork  = -1;

   dgesdd_(job, &(**a).iRows, &(**a).iCols, *((**a).a), &(**a).iCols, *w, *((**eig_right).a), &(**eig_right).iCols, \
   *((**eig_left).a), &(**eig_left).iCols, tmp, &lwork, itmp, info);

   lwork  = (int)tmp[0];
   liwork = itmp[0];

printf("lwork %3d    %3d %3d %3d %3d\n",lwork,(**a).iRows,(**a).iCols,(**eig_left).iCols,(**eig_right).iCols);

   work = (double *) malloc(lwork * sizeof(double));

  if (work == NULL) {
    *info = -400;
    return;
  }

  iwork = (int *) malloc(8 * (**eig_right).iCols * sizeof(int));
  if (iwork == NULL) {
    *info = -404;
    return;
  }

   dgesdd_(job, &(**a).iRows, &(**a).iCols, *((**a).a), &(**a).iCols, *w, *((**eig_right).a), &(**eig_right).iCols, \
   *((**eig_left).a), &(**eig_left).iCols, work, &lwork, iwork, info);

  free(work);
  free(iwork);

}

void DiagonaliseSymmetrical(DoubleMatrix **a, double **w, DoubleMatrix **eig, char *job, char *uplo, int *info)

{

#ifndef ISOF2003

 double tmp[1], *work,vl,vu,abstol;
 int lwork, liwork,il,iu, itmp[1],m, *iwork, *issupz;
 char rnge='A';
 int ione=1;

#endif

#ifdef ISOF2003

   DiagonaliseSymmetricalv1(a,w,eig,job,uplo,info);

#else

  abstol = k_tol;
  issupz = (int *) malloc(2*(**a).iCols * sizeof(int));
  if (issupz == NULL) {
    *info = -402;
    return;
  }

   lwork=-1;
   liwork=-1;
   dsyevr_(job,&rnge, uplo, &(**a).iCols, *((**a).a), &(**a).iCols,&vl,&vu,&il,&iu,&abstol,&m, *w,\
       *((**eig).a),&(**eig).iCols,issupz,tmp, &lwork,itmp,&liwork, info,&ione,&ione,&ione);

   lwork=(int)tmp[0];
   liwork=itmp[0];
   work = (double *) malloc(lwork * sizeof(double));

  if (work == NULL) {
    *info = -400;
    return;
  }

  iwork = (int *) malloc(liwork * sizeof(int));
  if (iwork == NULL) {
    *info = -404;
    return;
  }

  dsyevr_(job,&rnge, uplo, &(**a).iCols, *((**a).a), &(**a).iCols,&vl,&vu,&il,&iu,&abstol,&m, *w,\
     *((**eig).a),&(**eig).iCols,issupz,work, &lwork,iwork,&liwork, info,&ione,&ione,&ione);

  free(work);
  free(iwork);
  free(issupz);

#endif

}

void DiagonaliseHermitian(ComplexMatrix **a, double **w, ComplexMatrix **eig, char *job, char *uplo, int *info)

{

#ifndef ISOF2003

 double vl,vu,abstol,rtmp[1],*rwork;
 Complex wtmp[1], *work;
 int lwork, liwork,lrwork,il,iu, itmp[1],m, *iwork, *issupz;
 char rnge='A';
 int ione=1;

#endif

#ifdef ISOF2003

    DiagonaliseHermitianv1(a,w,eig,job,uplo,info);

#else

  int i;
  abstol = k_tol;
  issupz = (int *) malloc(2*(**a).iCols * sizeof(int));

  if (issupz == NULL) {
    *info = -402;
    return;
  }

   lwork=-1;
   liwork=-1;
   lrwork=-1;
   zheevr_(job, &rnge, uplo, &(**a).iCols, *((**a).a), &(**a).iCols, &vl, &vu, &il, &iu, &abstol, &m, *w,\
         *((**eig).a), &(**eig).iCols, issupz, wtmp, &lwork, rtmp, &lrwork, itmp, &liwork, info, &ione, &ione, &ione);
   lwork=(int)wtmp[0].real();
   lrwork=(int)rtmp[0];
   liwork=itmp[0];
   work = (Complex *) malloc(lwork * sizeof(Complex));

  if (work == NULL) {
    *info = -400;
    return;
  }

  rwork = (double *) malloc(lrwork * sizeof(double));

  if (rwork == NULL) {
    *info = -406;
    return;
  }

  iwork = (int *) malloc(liwork * sizeof(int));

  if (iwork == NULL) {
    *info = -404;
    return;
  }

  zheevr_(job,&rnge, uplo, &(**a).iCols, *((**a).a), &(**a).iCols,&vl,&vu,&il,&iu,&abstol,&m, *w,\
        *((**eig).a),&(**eig).iCols,issupz,work, &lwork,rwork, &lrwork,iwork,&liwork, info, &ione,&ione,&ione);

  free(work);
  free(iwork);
  free(rwork);
  free(issupz);

#endif

}

void DiagonaliseHermitianX(ComplexMatrix **a, double **w, ComplexMatrix **eig, char *job, char *uplo, int *info)

{

#ifndef ISOF2003

 double vl,vu,abstol,rtmp[1],*rwork;
 Complex wtmp[1], *work;
 int lwork, il, iu, itmp[1],m, *iwork, *ifail;
 //char rnge='V';
//vl = -16.0;
//vu = 10000.0;
 char rnge='A';
 int ione=1;

*info = 0;

#endif

#ifdef ISOF2003

    DiagonaliseHermitianv1(a,w,eig,job,uplo,info);

#else

  int i;
  abstol = k_tol;
  //abstol = ktol;
  ifail = (int *) malloc((**a).iCols * sizeof(int));

  for (i=0;i<(**a).iCols;i++)
  ifail[i] = 0;

  if (ifail == NULL) {
    *info = -402;
    return;
  }

   lwork=-1;
   zheevx_(job, &rnge, uplo, &(**a).iCols, *((**a).a), &(**a).iCols, &vl, &vu, &il, &iu, &abstol, &m, *w,\
         *((**eig).a), &(**eig).iCols, wtmp, &lwork, rtmp, itmp, ifail, info);
   lwork=(int)wtmp[0].real();
   work = (Complex *) malloc(lwork * sizeof(Complex));

  if (work == NULL) {
    *info = -400;
    return;
  }

  rwork = (double *) malloc(7 * (**eig).iCols * sizeof(double));

  if (rwork == NULL) {
    *info = -406;
    return;
  }

  iwork = (int *) malloc(5 * (**eig).iCols * sizeof(int));

  if (iwork == NULL) {
    *info = -404;
    return;
  }

  zheevx_(job,&rnge, uplo, &(**a).iCols, *((**a).a),&(**a).iCols,&vl,&vu,&il,&iu,&abstol,&m,*w,\
        *((**eig).a),&(**eig).iCols,work,&lwork,rwork,iwork,ifail,info);

  free(work);
  free(iwork);
  free(rwork);
  free(ifail);

#endif

}

void DiagonaliseRealGeneral(DoubleMatrix **a, double **w, double **u, DoubleMatrix **eig, int *info)

{

#ifndef ISOF2003

 double tmp[1], *work;
 int lwork, *iwork;
 char jobvl='N', jobvr='V'; // Only right eigenvectors are calculated

#endif

   lwork=-1;
   dgeev_(&jobvl, &jobvr, &(**a).iCols, *((**a).a), &(**a).iCols, *w, *u, *((**eig).a), &(**eig).iCols, *((**eig).a),&(**eig).iCols,\
   tmp, &lwork, info);

   lwork=(int)tmp[0];
   work = (double *) malloc(lwork * sizeof(double));

  if (work == NULL) {
    *info = -400;
    return;
  }

   dgeev_(&jobvl, &jobvr, &(**a).iCols, *((**a).a), &(**a).iCols, *w, *u, *((**eig).a),&(**eig).iCols, *((**eig).a),&(**eig).iCols,\
   work, &lwork, info);

  free(work);

}

void DiagonaliseComplexGeneral(ComplexMatrix **a, Complex **w, ComplexMatrix **eig, int *info)

{

#ifndef ISOF2003

 Complex tmp[1], *work;
 int lwork, *iwork;
 double *rwork;
 char jobvl='N', jobvr='V'; // Only right eigenvectors are calculated

#endif

   lwork=-1;

   rwork = (double *) malloc(2 * (**a).iCols * sizeof(double));

   zgeev_(&jobvl, &jobvr, &(**a).iCols, *((**a).a), &(**a).iCols, *w, *((**eig).a), &(**eig).iCols, *((**eig).a),&(**eig).iCols,\
   tmp, &lwork, rwork, info);

   lwork=(int)tmp[0].real();
   work = (Complex *) malloc(lwork * sizeof(Complex));

  if (work == NULL) {
    *info = -400;
    return;
  }

   zgeev_(&jobvl, &jobvr, &(**a).iCols, *((**a).a), &(**a).iCols, *w, *((**eig).a),&(**eig).iCols, *((**eig).a),&(**eig).iCols,\
   work, &lwork, rwork, info);

  free(work);
  free(rwork);

}

void CholeskySymmetric(DoubleMatrix **S, int *pivot, int *dimension, double *threshold)

{

int dim = (**S).iRows;
int i, j, k, t, max1;
double max;
double a, v[dim], swap[dim];

  for (i = 0; i < dim; i++) pivot[i] = i;
  for (i = 0; i < dim; i++) {
    max = (**S).a[i][i];
    max1 = i;
    for (j = i; j < dim; j++) {
      if ((**S).a[j][j] > max) {
      max = (**S).a[j][j];
      max1 = j;
     }
    }

  if (max1 != i) {
  t = pivot[max1];
  pivot[max1] = pivot[i];
  pivot[i] = t;
  for (k = 0; k < dim; k++) {
    swap[k] = (**S).a[max1][k];
   }
  for (k = 0; k < dim; k++) {
    (**S).a[max1][k] = (**S).a[i][k];
   }
  for (k = 0; k < dim; k++) {
    (**S).a[i][k] = swap[k];
   }
  for (k = 0; k < dim; k++) {
    swap[k] = (**S).a[k][max1];
   }
  for (k = 0; k < dim; k++) {
    (**S).a[k][max1] = (**S).a[k][i];
   }
  for (k = 0; k < dim; k++) {
    (**S).a[k][i] = swap[k];
   }
  } // close if (max1

  if ((**S).a[i][i] <= k_zero) { (**S).a[i][i] = k_zero; continue; }

  a = sqrt((**S).a[i][i]);
  (**S).a[i][i] = a;
  for (j = i + 1; j < dim; j++) {
    (**S).a[j][i] /= a;
    v[j] = (**S).a[j][i];
   }
  for (j = i + 1; j < dim; j++) {
    for (k = i + 1; k < dim; k++) {
      (**S).a[j][k] -= v[j] * v[k];
     }
    }
  } // close loop on i

  *dimension = 0;
  for (i = 0; i < dim; i++) {
    if ((**S).a[i][i] > *threshold) (*dimension)++;
    //printf("%3d %3d %16.10lf %16.10lf\n",i,*dimension,(**S).a[i][i],*threshold);
   }

  for (i = 0; i < dim; i++) {
    for (j = i + 1; j < dim; j++) {
        (**S).a[i][j] = k_zero;
       }
      }

}

void CholeskySymmetricNoPivot(DoubleMatrix **S)

{

int dim = (**S).iRows;
int i, j, k;
double a, v[dim];

  for (i = 0; i < dim; i++) {
  a = sqrt((**S).a[i][i]);

  if ((**S).a[i][i] <= k_zero) { 
  printf("Element %3d negative %f\n",i,(**S).a[i][i]); 
  (**S).a[i][i] = 1e-04; 
  continue; }

  (**S).a[i][i] = a;
  for (j = i + 1; j < dim; j++) {
    (**S).a[j][i] /= a;
    v[j] = (**S).a[j][i];
   }
  for (j = i + 1; j < dim; j++) {
    for (k = i + 1; k < dim; k++) {
      (**S).a[j][k] -= v[j] * v[k];
     }
    }
  } // close loop on i

  for (i = 0; i < dim; i++) {
    for (j = i + 1; j < dim; j++) {
        (**S).a[i][j] = k_zero;
       }
      }

}

void CholeskyInverseSymmetric(DoubleMatrix **S, DoubleMatrix **S_inverse, int *dimension)

{

int dim = (**S_inverse).iRows;
int i, j, k;

  ResetDoubleMatrix(*S_inverse);
  for (i = 0; i < dim; i++) {
    for (j = i; j < dim; j++) {
      if (i == j) (**S_inverse).a[i][j] = k_one;
      for (k = 0; k < j; k++) {
        (**S_inverse).a[j][i] -= (**S).a[j][k] * (**S_inverse).a[k][i];
       }
        (**S_inverse).a[j][i] /= (**S).a[j][j];
      }
     }

}

void CholeskyPermuteSymmetric(DoubleMatrix **S, int *pivot, int transform, JOB_PARAM *job)

{

int dim = (**S).iCols;
int i, j;
double swap[dim];
DoubleMatrix *temp;

 AllocateDoubleMatrix(&temp,&dim,&dim,job);

  // Row and column permutation by permutation matrix in pivot
  // transform == 0 P   . S . P^T
  // transform == 1 P^T . S . P
  // transform == 2 row permutation 
  // transform == 3 column permutation

  switch (transform) {

  case 0 :

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      temp->a[i][j] = (**S).a[pivot[i]][pivot[j]];
     }
    }
  break;

  case 1 :

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      temp->a[pivot[i]][pivot[j]] = (**S).a[i][j];
     }
    }
  break;

  case 2 :

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      temp->a[i][pivot[j]] = (**S).a[i][j]; // row permutation
     }
    }
  break;

  case 3 :

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      temp->a[pivot[i]][j] = (**S).a[i][j]; // column permutation
     }
    }
  break;

  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      (**S).a[i][j] = temp->a[i][j];
     }
    }

  DestroyDoubleMatrix(&temp,job);

}

void CholeskyMultiply(DoubleMatrix **S, DoubleMatrix **S1, char *trans1, DoubleMatrix **S2, char *trans2)

{

  // ******************************************************************************************
  // * Multiply two lower or upper triangular matrices with or without transpose              *
  // ******************************************************************************************

int dimi = (**S1).iRows, dimj = (**S2).iCols;
int i, j, k;

  ResetDoubleMatrix(*S);

  if (*trans1 == 'Y' && *trans2 == 'Y') {
  for (i = 0; i < dimi; i++) {
    for (j = 0; j < dimj; j++) {
      for (k = 0; k < dimj; k++) {
        (**S).a[i][j] += (**S1).a[k][i] * (**S2).a[j][k];
       }
      }
     }
    }

  if (*trans1 == 'Y' && *trans2 == 'N') {
  for (i = 0; i < dimi; i++) {
    for (j = 0; j < dimj; j++) {
      for (k = 0; k < dimj; k++) {
        (**S).a[i][j] += (**S1).a[k][i] * (**S2).a[k][j];
       }
      }
     }
    }

  if (*trans1 == 'N' && *trans2 == 'Y') {
  for (i = 0; i < dimi; i++) {
    for (j = 0; j < dimj; j++) {
      for (k = 0; k < dimj; k++) {
        (**S).a[i][j] += (**S1).a[i][k] * (**S2).a[j][k];
       }
      }
     }
    }

  if (*trans1 == 'N' && *trans2 == 'N') {
  for (i = 0; i < dimi; i++) {
    for (j = 0; j < dimj; j++) {
      for (k = 0; k < dimj; k++) {
        (**S).a[i][j] += (**S1).a[i][k] * (**S2).a[k][j];
       }
      }
     }
    }

}

void CholeskyHermitian(ComplexMatrix **S, int *pivot, int *dimension, double *threshold)

{

int dim = (**S).iRows;
int i, j, k, t, max1;
double max;
Complex a, v[dim], swap[dim];

  for (i = 0; i < dim; i++) pivot[i] = i;
  for (i = 0; i < dim; i++) {
    max = ((**S).a[i][i]).real();
    max1 = i;
    for (j = i; j < dim; j++) {
      if (((**S).a[j][j]).real() > max) {
      max = ((**S).a[j][j]).real();
      max1 = j;
     }
    }
  if (max1 != i) {
  t = pivot[max1];
  pivot[max1] = pivot[i];
  pivot[i] = t;
  for (k = 0; k < dim; k++) {
    swap[k] = (**S).a[max1][k];
   }
  for (k = 0; k < dim; k++) {
    (**S).a[max1][k] = (**S).a[i][k];
   }
  for (k = 0; k < dim; k++) {
    (**S).a[i][k] = swap[k];
   }
  for (k = 0; k < dim; k++) {
    swap[k] = (**S).a[k][max1];
   }
  for (k = 0; k < dim; k++) {
    (**S).a[k][max1] = (**S).a[k][i];
   }
  for (k = 0; k < dim; k++) {
    (**S).a[k][i] = swap[k];
   }
  } // close if (max1
  //if (((**S).a[i][i]).real() <= k_zero) { ((**S).a[i][i]).real() = k_zero; continue; }
  if (real((**S).a[i][i]) <= k_zero) { (**S).a[i][i] = Complex (k_zero, k_zero); continue; }
  a = sqrt((**S).a[i][i]);
  (**S).a[i][i] = a;
  for (j = i + 1; j < dim; j++) {
    (**S).a[j][i] /= a;
    v[j] = (**S).a[j][i];
   }
  for (j = i + 1; j < dim; j++) {
    for (k = i + 1; k < dim; k++) {
      (**S).a[j][k] -= v[j] * v[k];
     }
    }
  } // close loop on i

  *dimension = 0;
  for (i = 0; i < dim; i++) {
    if (((**S).a[i][i]).real() > *threshold) (*dimension)++;
    //printf("%3d %3d %16.10lf %16.10lf\n",i,*dimension,((**S).a[i][i]).real(),*threshold);
   }

  for (i = 0; i < dim; i++) {
    for (j = i + 1; j < dim; j++) {
        //((**S).a[i][j]).real() = k_zero;
        (**S).a[i][j] = Complex (k_zero, k_zero);
       }
      }

 //printf("Full dimension %4d Cholesky dimension %4d\n",dim, *dimension);

}

void CholeskyInverseHermitian(ComplexMatrix **S, ComplexMatrix **S_inverse, int *dimension)

{

int dim = (**S_inverse).iRows;
int i, j, k;

  ResetComplexMatrix(*S_inverse);
  for (i = 0; i < dim; i++) {
    for (j = i; j < dim; j++) {
      //if (i == j) ((**S_inverse).a[i][j]).real() = k_one;
      if (i == j) (**S_inverse).a[i][j] = Complex (k_one, k_zero);
      for (k = 0; k < j; k++) {
        (**S_inverse).a[j][i] -= (**S).a[j][k] * (**S_inverse).a[k][i];
       }
        (**S_inverse).a[j][i] /= (**S).a[j][j];
      }
     }

}

void CholeskyPermuteHermitian(ComplexMatrix **S, int *pivot, int transform, JOB_PARAM *job)

{

int dim = (**S).iCols;
int i, j;
Complex swap[dim];
ComplexMatrix *temp;

 AllocateComplexMatrix(&temp,&dim,&dim,job);

  // Row and column permutation by permutation matrix in pivot
  // transform == 0 P   . S . P^T
  // transform == 1 P^T . S . P
  // transform == 2 row permutation 
  // transform == 3 column permutation

  switch (transform) {

  case 0 :

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      temp->a[i][j] = (**S).a[pivot[i]][pivot[j]];
     }
    }
  break;

  case 1 :

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      temp->a[pivot[i]][pivot[j]] = (**S).a[i][j];
     }
    }
  break;

  case 2 :

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      temp->a[i][pivot[j]] = (**S).a[i][j]; // row permutation
     }
    }
  break;

  case 3 :

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      temp->a[pivot[i]][j] = (**S).a[i][j]; // column permutation
     }
    }
  break;

  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      (**S).a[i][j] = temp->a[i][j];
     }
    }

  DestroyComplexMatrix(&temp,job);

}
