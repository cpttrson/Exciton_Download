
void initial_density_matrix(double*, double*, PAIR_TRAN*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void read_density_matrix(FERMI*,double**,int*,int*,ATOM*,JOB_PARAM*,FILES);

void read_JOB_PARAM_array(FERMI*, int*, double*, JOB_PARAM*, FILES);

void write_JOB_PARAM_array(FERMI*, int*, double*, JOB_PARAM*, FILES);

void density_matrix_molecule2(FERMI*, double*, double*, KPOINT_TRAN*, int*, REAL_LATTICE*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void reduced_density_matrix_molecule2(FERMI*, double*, KPOINT_TRAN*, int*, REAL_LATTICE*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void reduced_density_matrix_molecule3(int*, double*, FERMI*, double*, KPOINT_TRAN*, int*, REAL_LATTICE*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void expand_density_matrix(double*, double*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void atom_shell_populations2(INT_1E*, double*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

double print_total_population2(ATOM*, SHELL*, JOB_PARAM*, FILES);

void print_atom_populations2(ATOM*, SHELL*, JOB_PARAM*, FILES);

void print_shell_populations2(ATOM*, SHELL*, JOB_PARAM*, FILES);

void init_mix(DoubleMatrix*, DoubleMatrix*, DoubleMatrix*, double*, double*);

void mix_density(DoubleMatrix*, DoubleMatrix*, DoubleMatrix*, double*, JOB_PARAM*, FILES);

/*
void expand_screening_integral_matrix(double*, double*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void density_matrix_crystal2(FERMI*, double*, double*, KPOINT_TRAN*, char*, int*, REAL_LATTICE*, REAL_LATTICE_TABLES*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void reduced_density_matrix_crystal(FERMI*, double*, KPOINT_TRAN*, char*, REAL_LATTICE*, REAL_LATTICE_TABLES*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*,FILES);

void reduced_density_matrix_crystal3(FERMI*, double*, KPOINT_TRAN*, int, int, ComplexMatrix*, REAL_LATTICE*, REAL_LATTICE_TABLES*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*,FILES);

void reduced_density_matrix_crystal1(FERMI*, double*, KPOINT_TRAN*, char*, REAL_LATTICE*, REAL_LATTICE_TABLES*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*,FILES);

void reduced_density_matrix_crystal2(FERMI*, double*, KPOINT_TRAN*, int*, REAL_LATTICE*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*,FILES);

void expand_density_matrix_complex(Complex*, Complex*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void expand_screening_integral_matrix_complex(Complex*, Complex*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

double normalise_density_matrix2(double*, double*, double, JOB_PARAM*, FILES);
*/

/*
//void C_DIIS_extrapolation(ComplexMatrix**, ComplexMatrix**, ComplexMatrix**, ComplexMatrix**, ComplexMatrix**, int*, int*, SYMMETRY*, JOB_PARAM*, FILES);

void count_density_matrix_shells2(int*, PAIR_TRAN*, ATOM*, JOB_PARAM*, FILES);

void density_matrix_largest_shell_elements2(double*, double*, PAIR_TRAN*, PAIR_TRAN*, ATOM*, SHELL*, JOB_PARAM*, FILES);

void fetch_scf_vectors(ComplexMatrix*, FERMI*, ATOM*, CRYSTAL*, JOB_PARAM*, FILES);

void fetch_scf_eigenvalues(double*, FERMI*, CRYSTAL*, JOB_PARAM*, FILES);
*/

