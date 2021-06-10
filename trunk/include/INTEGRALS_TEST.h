#ifndef INTEGRALS_TEST
#define INTEGRALS_TEST

void integrals_molecule_test(FERMI*, ATOM*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void integrals_crystal_coulomb_test(FERMI*, ATOM*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void integrals_crystal_exchange_test(FERMI*, ATOM*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void integrals_crystal_exchange_energy_test(FERMI*, ATOM*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void integrals_crystal_ijkl3_test(INTEGRAL_LIST*, QUAD_TRAN*, int*, int*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_molecule_ijkl1_test(INTEGRAL_LIST*, int, int, int, int, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_exchange_ijkl1_test(INTEGRAL_LIST*, QUAD_TRAN*, int*, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb1_reversed2_test(INTEGRAL_LIST*, double*, TRIPLE_TRAN*, int, int, int, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void three_centre_coulomb1_reversed2_crystal_test(int, TRIPLE_TRAN*, Complex*, REAL_LATTICE*, Q_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void contract_coulomb_integrals_test(int*, int*, int*, int*, DoubleMatrix*, double*, double*, ATOM*, ATOM*, JOB_PARAM*, FILES);

void contract_coulomb_integrals_complex_test(int*, int*, int*, int*, ComplexMatrix*, Complex*, Complex*, ATOM*, ATOM*, JOB_PARAM*, FILES);

void contract_coulomb_integrals_complex_bands_test(int*, int*, int*, int*, ComplexMatrix*, Complex*, Complex*, ATOM*, ATOM*, JOB_PARAM*, FILES);

void two_centre_coulomb1_crystal_test(ComplexMatrix*, REAL_LATTICE*, Q_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void two_centre_exchange_crystal1_test(int, Complex*, PAIR_TRAN*, REAL_LATTICE*, Q_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

#endif
