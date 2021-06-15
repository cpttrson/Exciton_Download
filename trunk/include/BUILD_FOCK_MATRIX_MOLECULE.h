
void build_fock_matrix_molecule(double*, INT_1E*, FERMI*, double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void build_fock_matrix_molecule_no_sym(double*, INT_1E*, FERMI*, double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void total_energy(INT_1E*, double*, double*, PAIR_TRAN*, ATOM*, JOB_PARAM*, FILES);

void total_energy_final(INT_1E*, double*, double*, double*, double*, PAIR_TRAN*, ATOM*, JOB_PARAM*, FILES);

void total_energy_direct(INT_1E*, double*, double*, PAIR_TRAN*, ATOM*, JOB_PARAM*, FILES);

void fock_matrix_molecule_compute_integrals(double*, double*, double*, double*, PAIR_TRAN*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void fock_matrix_molecule_read_integrals(double *Fock_2c, double *Fock_2e, double *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file);

void fock_matrix_molecule_compute_integrals_no_sym(double *Fock_2c, double *Fock_2e, double *S1, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, JOB_PARAM *job, FILES file);

void fock_matrix_molecule_read_integrals_no_sym(double *Fock_2c, double *Fock_2e, double *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file);

void contract_integrals_molecule_ijkl(double *Fock_2c_buffer, double *Fock_2e_buffer, INTEGRAL_LIST *integral_list_molecule, double *F, PAIR_TRAN *pair_p, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, JOB_PARAM *job, FILES file);

void fock_matrix_molecule_compute_screening_integrals(double*, PAIR_TRAN*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void expand_screening_integral_matrix(double*, double*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void integrals_molecule_screen_ijkl(int*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM*, SHELL*, JOB_PARAM*, FILES);

void integrals_molecule_screen_direct_ijkl(int*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void print_Fock_matrix(double *Fock, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file);

