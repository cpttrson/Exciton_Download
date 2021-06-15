
void integrals_molecule_ijkl(INTEGRAL_LIST*, int, int, int, int, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, JOB_PARAM*, FILES);
//void integrals_molecule_ijkl(INTEGRAL_LIST*, int, int, int, int, int*, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_molecule_ijij(double*, int, int, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void pack_write_molecule_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void read_unpack_molecule_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

