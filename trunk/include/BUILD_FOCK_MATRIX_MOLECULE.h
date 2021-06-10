#ifndef BUILD_FOCK_MATRIXH
#define BUILD_FOCK_MATRIXH

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

void shell_screen_molecule_compute_integrals(double*, PAIR_TRAN*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void expand_screening_integral_matrix(double*, double*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void shell_screen1(int*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM*, SHELL*, JOB_PARAM*, FILES);

void shell_screen_direct(int*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void print_Fock_matrix(double *Fock, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file);

//void coulomb_matrix_crystal_compute_integrals(double *Fock_2c, double *S1, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file);

//void exchange_matrix_crystal_compute_integrals(double *Fock_2e, double *S2, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM_TRAN *atom_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file);

//void coulomb_matrix_crystal_read_integrals(double *Fock_2c, double *F, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file);

//void coulomb_matrix_crystal_compute_integrals_no_sym(double *Fock_2c, double *S1, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file);

//void exchange_matrix_crystal_compute_integrals_no_sym(double *Fock_2e, double *S2, double *F, PAIR_TRAN *pair_p, PAIR_TRAN *pair_q, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file);

//void coulomb_matrix_crystal_read_integrals_no_sym(double *Fock_2c, double *F, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file);

//void exchange_matrix_crystal_read_integrals_no_sym(double *Fock_2e, double *F, PAIR_TRAN *pair_p, ATOM *atoms, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file);

//void contract_coulomb_integrals(double *Fock_2c_buffer, INTEGRAL_LIST *integral_list_coulomb, double *F, PAIR_TRAN *pair_p, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file);

//void contract_exchange_integrals4(double *Fock_2e_buffer, INTEGRAL_LIST *integral_list_exchange, double *F, PAIR_TRAN *pair_p, QUAD_TRAN *quad, ATOM *atoms, SHELL *shells, SYMMETRY *symmetry, REAL_LATTICE_TABLES *R_tables, JOB_PARAM *job, FILES file);

//void shell_screen_complex(Complex*, PAIR_TRAN*, REAL_LATTICE*, RECIPROCAL_LATTICE*, Q_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

//void shell_screen2(int*, double*, PAIR_TRAN*, QUAD_TRAN*, REAL_LATTICE_TABLES*, ATOM*, SHELL*, JOB_PARAM*, FILES);

//void shell_screen3(int*, ComplexMatrix*, Complex*, PAIR_TRAN*, TRIPLE_TRAN*, ATOM*, SHELL*, ATOM*, SHELL*, JOB_PARAM*, FILES);

//void shell_screen_crystal_coulomb_direct(int*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM*, SHELL*, REAL_LATTICE_TABLES*, SYMMETRY*, JOB_PARAM*, FILES);

//void shell_overlap(int*, QUAD_TRAN*, ATOM*, SHELL*, REAL_LATTICE*, JOB_PARAM*, FILES);

#endif
