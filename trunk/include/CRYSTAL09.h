#ifndef CRYSTAL09H
#define CRYSTAL09H

#if defined(_athlone)
    #define FORT_HEADER 8
    #define FORT_TRAILER 16
#else
    #define FORT_HEADER 4
    #define FORT_TRAILER 8
#endif

void read_SYMMOP_crystal_09(SYMMETRY*, JOB_PARAM*, FILES);

void read_XCBD_crystal_09(KPOINT_TRAN*, CRYSTAL*, JOB_PARAM*, FILES);
//void read_XCBD_crystal_09(int, KPOINT_TRAN*, CRYSTAL*, JOB_PARAM*, FILES);

void write_vectors_write_values_rma3(int *, KPOINT_TRAN *, int *, ATOM *, SYMMETRY*, CRYSTAL*, JOB_PARAM *, FILES);
//void write_vectors_write_values_rma3(int *, int, KPOINT_TRAN *, int *, ATOM *, SYMMETRY*, CRYSTAL*, JOB_PARAM *, FILES);

//void read_write_evectors_crystal_09(KPOINT_TRAN*, int, int [][3], int*, int*, int*, double*, double*, char*, ATOM*, JOB_PARAM*, FILES);
void read_write_evectors_crystal_09(KPOINT_TRAN*, int, int [][3], int*, int*, int*, char*, ATOM*, JOB_PARAM*, FILES file);


void read_evalues_crystal_09(int , int,  int*, int*, int*, int [][3], int [][3], double*, double*, ATOM*, JOB_PARAM*, FILES);

#endif
/*
void read_evectors_crystal_09(KPOINT_TRAN*, int,  int [][3], int*, int*, int*, double*, double*, char*, ATOM*, JOB_PARAM*, FILES);
void get_vectors_values(int *, int *, int *, VECTOR_KNET *, double *, double *, ATOM *, JOB_PARAM *, FILES);
void get_vectors_get_values(int *, int *, int *, VECTOR_KNET *, ComplexMatrix *, double *, ATOM *, JOB_PARAM *, FILES);
void get_values(int *is, int *nk, int *bands, VECTOR_KNET *knet, double *e, ATOM *atoms, JOB_PARAM *, FILES);
void read_evectors(int , int *, int *, double *, double *, char *, ATOM *, JOB_PARAM *, FILES);
void read_vectors(int , int *, int *, double *, ComplexMatrix *, char *, ATOM *, JOB_PARAM *, FILES);
void read_eivalues(int nk, int *bands, int *latt, double *eigenvalues, ATOM *atoms, JOB_PARAM *, FILES);
void read_evalues(int nk, int *bands, int *latt, int kk[][3], double *eig1, double *eig2, ATOM *atoms, JOB_PARAM *, FILES);
void run_properties(char crystal, int kk, int test, JOB_PARAM *, FILES);
void write_vectors_get_values(int *, int, VECTOR_KNET *, int *, double *, ATOM *, JOB_PARAM *, FILES);
void write_vectors_write_values_rma(int *, int, VECTOR_KNET *, int *, double *, ATOM *, JOB_PARAM *, FILES);
void write_vectors_write_values_rma1(int *, int, KPOINT_TRAN *, int *, ATOM *, JOB_PARAM *, FILES);
void write_vectors_write_values_rma2(int *, int, KPOINT_TRAN *, int *, ATOM *, SYMMETRY*, CRYSTAL*, JOB_PARAM *, FILES);
void write_vectors_write_values_rma3(int *, int, KPOINT_TRAN *, int *, ATOM *, SYMMETRY*, CRYSTAL*, JOB_PARAM *, FILES);
void write_vectors_rma1(int *, int, KPOINT_TRAN *, int *, ATOM *, JOB_PARAM *, FILES);
void write_values_rma1(int *, int, KPOINT_TRAN *, int *, ATOM *, JOB_PARAM *, FILES);
*/
