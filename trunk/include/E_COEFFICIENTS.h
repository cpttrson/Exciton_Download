#ifndef ECOEFF
#define ECOEFF

void E_coefficients(int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void E_coefficients1(int,int,int,int,int,int,int,int,int,int,double*,double*,double*,double*,double*,double*,double*,VECTOR_DOUBLE*,REAL_LATTICE*,ATOM*,SHELL*,GAUSSIAN*,JOB_PARAM*,FILES);
 
void E_coefficients_1c(int, int, double*, double*, double*, double*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void E_coefficients_2c(int,int,int,int,int,int,int,int,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,REAL_LATTICE*,ATOM*,SHELL*,GAUSSIAN*,JOB_PARAM*,FILES);

#endif
