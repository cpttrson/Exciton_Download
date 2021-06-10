#ifndef INTEGRALS3C
#define INTEGRALS3C

void three_centre_overlap1(double*, int, int, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_overlap2(double*, TRIPLE_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_overlap3(double*, TRIPLE_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb(int, double*, TRIPLE_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb1(int, double*, TRIPLE_TRAN*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb1_reversed(int, double*, TRIPLE_TRAN*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb1_reversed1(int, int, PAIR_TRAN*, double*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb1_reversed2(int, TRIPLE_TRAN*, double*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb1_reversed2_crystal(int, TRIPLE_TRAN*, Complex*, REAL_LATTICE*, Q_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb1_reversed2_crystal_test1(int, TRIPLE_TRAN*, int*, Complex*, REAL_LATTICE*, Q_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb1_reversed2_crystal_test2(INTEGRAL_LIST_COMPLEX*, int, TRIPLE_TRAN*, int*, Complex*, REAL_LATTICE*, Q_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void E_coefficients2(int, int, int, int, int, int, int, int, double*, double*, double, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void E_coefficients3a(int, int, int, int, int, int, int, int, int, int, double*, double*, double*, double*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void E_coefficients3b(int, int, int, int, int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void E_coefficients3c(int, int, int, int, int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void E_coefficients3(int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void mcmurchie_davidson_3c(double*, int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, SHELL*, SHELL*, JOB_PARAM*, FILES);

void mcmurchie_davidson_3c_reversed(double*, int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, SHELL*, SHELL*, JOB_PARAM*, FILES);

void mcmurchie_davidson_3c_reversed_complex(Complex*, Complex*, int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, SHELL*, SHELL*, JOB_PARAM*, FILES);

void non_recursive_three_centre_overlap_primitives(double*, int, int, int, int, int, int, int, int, int, int, double*, double*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);;

#endif


