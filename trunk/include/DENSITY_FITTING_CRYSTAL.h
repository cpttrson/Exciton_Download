
void generate_coulomb_matrix_inverse_complex(ComplexMatrix*, int, FERMI*, ATOM_TRAN*, ATOM*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void density_fitting_crystal_contract_integrals(int*, int*, int*, int*, int*, PAIR_TRAN*, KPOINT_TRAN*, MPI_File, Complex**, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, SYMMETRY*, Q_LATTICE*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void density_fitting_crystal_rotate_integrals(int*, int, KQPOINT_TRAN*, ComplexMatrix*, ComplexMatrix*, Complex**, Complex*, int*, int*, int*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, FILES, JOB_PARAM*);

