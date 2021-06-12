
void integrals_crystal_coulomb_ijkl(INTEGRAL_LIST*, QUAD_TRAN*, int*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_crystal_exchange_ijkl(INTEGRAL_LIST*, QUAD_TRAN*, int*, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_crystal_coulomb_screen_ijij(double*, int, int, int, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_crystal_exchange_screen_ijij(double*, int, int, int, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_crystal_screen_complex(Complex*, int, int, int, REAL_LATTICE*, RECIPROCAL_LATTICE*, Q_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void expand_screening_integral_matrix_complex(Complex*, Complex*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void pack_write_integrals_crystal_coulomb_ijkl(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void read_unpack_integrals_crystal_coulomb_ijkl(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void pack_write_integrals_crystal_exchange_ijkl(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void read_unpack_integrals_crystal_exchange_ijkl(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

