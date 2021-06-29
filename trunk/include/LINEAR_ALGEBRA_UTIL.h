#ifndef LINALGUTILH
#define LINALGUTILH

/*
extern "C" {

  #include "mycomplex.h"
}
*/

extern "C" {

void ComplexDotProd(Complex *, const int *, Complex **, const int *, Complex **, const int *);

void DoubleDotProd(double *, const int *, double **, const int *, double **, const int *);

void DoubleGEMM(const char *, const char *, const double *,DoubleMatrix **,  DoubleMatrix **,const double *, DoubleMatrix **);

void DoubleGEMMv1(const char *, const char *, const double *, DoubleMatrix **,  DoubleMatrix **, const double *, DoubleMatrix **);

void DoubleGEMM3(const char *, const char *, int *, int *, int *, const double *, double **,  int *, double**, int*, const double*, double**, int*);

void DoubleGEMM2v1(const char *, const char *, int *, int *, int *, const double *, double **,  int *, double**, int*, const double*, double**, int*);

void ComplexGEMM(const Complex *, ComplexMatrix **,  ComplexMatrix **, const Complex *, ComplexMatrix **);

void ComplexGEMMv1(const char *, const char *,const Complex *,ComplexMatrix **,  ComplexMatrix **,const Complex *, ComplexMatrix **);

void ComplexGEMM1(const char *, const char *, const Complex *, ComplexMatrix **,  ComplexMatrix **, const Complex *, ComplexMatrix **);

void ComplexGEMM1v1(const char *, const char *, const Complex *, ComplexMatrix **,  ComplexMatrix **, const Complex *, ComplexMatrix **);

void ComplexGEMM2(const char*, const char*, int*, int*, int*, const Complex*, ComplexMatrix**, int*, ComplexMatrix**, int*, const Complex*, ComplexMatrix**, int*);

void ComplexGEMM2v1(const char *, const char *, const Complex *, ComplexMatrix **,  ComplexMatrix **, const Complex *, ComplexMatrix **);

void ComplexGEMM3(const char *, const char *, int *, int *, int *, const Complex *, Complex **,  int *, Complex**, int*, const Complex*, Complex**, int*);

void DoubleTRMM3(const char *, const char *, const char *, const char *, int *, int *, const double *, double **, int *,  double **, int *);

void SVD_DGESDD(DoubleMatrix**, double**, DoubleMatrix**, DoubleMatrix**, char*, int*);

void DiagonaliseSymmetrical(DoubleMatrix **,double **, DoubleMatrix **, char *, char *,int *);

void DiagonaliseHermitian(ComplexMatrix **,double **, ComplexMatrix **, char *, char *,int *);

void DiagonaliseHermitianX(ComplexMatrix **,double **, ComplexMatrix **, char *, char *,int *);

void DiagonaliseRealGeneral(DoubleMatrix **, double **, double **, DoubleMatrix **, int *);

void DiagonaliseComplexGeneral(ComplexMatrix**, Complex**, ComplexMatrix**, int*);

#ifdef ISOF2003

void ComplexDotProdv1(Complex *, const int *, Complex **, const int *, Complex **, const int *);

void DiagonaliseSymmetricalv1(DoubleMatrix **,double **, DoubleMatrix **, const char *,const char *,int *);

void DiagonaliseHermitianv1(ComplexMatrix **,double **, ComplexMatrix **, const char *,const char *,int *);

#else

void zdotc_(Complex *, const int *, Complex *, const int *, Complex *, const int *);

double ddot_(const int *, double *, const int *, double *, const int *);

void dgemm_(const char *, const char *, const int *, const int *,const int *, const double *, double *, const int *, double *, const int *, const double *, double *, const int *, int*, int * );

void cblas_dgemm_(const char *, const char *, const char *, const int *, const int *,const int *, const double *, double *, const int *, double *, const int *, const double *, double *, const int *);

void zgemm_(const char *, const char *, int *, int *, int *, const Complex *, Complex *, int *, Complex *, int *, const Complex *, Complex *, int *, int *, int *);

void dpotrf_(const char*, const int*, double*, const int*, int*);

void dtrmm_(const char *, const char *, const char *, const char *, const int *, const int *, const double *, double *, const int *, double *, const int *);

void dsyevr_(const char *, const char *,const char *,const int *, double *,const int *,const double *, const double *, const int *, const int *, const double *,int *, double *, double *,const int *, int *, double *, int*,int *, const int *, int *, int *, int *, int *);

//void dgesdd_(const char *, const int *, const int *, double *, const int *, const double *, const double *, const int *, const double *, const int *, double *, const int *, int *);
void dgesdd_(const char *, const int *, const int *, double *, const int *, const double *, const double *, const int *, const double *, const int *, double *, const int *, int *, int *);

void zheevr_(const char *, const char *, const char *,const int *, Complex *, const int *, const double *, const double *, const int *, const int *, const double *, int *,double *, Complex *, const int *, const int *, Complex *, const int *, double *, const int *,int *, const int *, int *, int *,int *, int *);

void zheevx_(const char *, const char *, const char *,const int *, Complex *, const int *, const double *, const double *, const int *, const int *, const double *, int *,double *, Complex *, const int *, Complex *, const int *, double *, int *,int *, int *);
//void zheevx_(const char *, const char *, const char *,const int *, Complex *, const int *, const double *, const double *, const int *, const int *, const double *, int *,double *, Complex *, const int *, Complex *, const int *, double *, int *, int *, int *, int *,int *, int *);

void dgesv_(const int*, const int*, double*, const int*, int*, double*, const int*, int*);

void zgesv_(const int*, const int*, Complex*, const int*, int*, Complex*, const int*, int*);

void dgeev_(const char*, const char*, const int*, double*, const int*, double*, double*, double*, const int*, double*, const int*, double *, int *, int *);

void zgeev_(const char*, const char*, const int*, Complex*, const int*, Complex*, Complex*, const int*, Complex*, const int*, Complex *, int *, double *, int *);

void pdgemm_(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const int*, const int*, const double*, const int*, const int*, const int*, const double*, double*, const int*, const int*, const int*);

void pdtrmm_(const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, const int*, const int*, double*, const int*, const int*, const int*);

void CholeskySymmetric(DoubleMatrix**, int*, int*, double*);

void CholeskySymmetricNoPivot(DoubleMatrix**);

void CholeskyInverseSymmetric(DoubleMatrix**, DoubleMatrix**, int*);

void CholeskyPermuteSymmetric(DoubleMatrix**, int*, int, JOB_PARAM*);

void CholeskyMultiply(DoubleMatrix**, DoubleMatrix**, char*, DoubleMatrix**, char*);

void CholeskyHermitian(ComplexMatrix**, int*, int*, double*);

void CholeskyInverseHermitian(ComplexMatrix**, ComplexMatrix**, int*);

void CholeskyPermuteHermitian(ComplexMatrix**, int*, int, JOB_PARAM*);

#endif

}

#endif
