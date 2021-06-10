#ifndef TDHF
#define TDHF


void bse_crystal1(FERMI*, ATOM*, ATOM_TRAN*_p, int*, int*, int [][2], SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void tdhf_hamiltonian_in_core_crystal_finite_q(int*, int*, int*, Complex*, Complex*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void tdhf_hamiltonian_in_core_crystal_zero_q(int*, int*, int*, Complex*, Complex*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void setup_hamiltonian_parameters_zero_q(int*, int*, int*, int*, int*, FERMI*, JOB_PARAM*, FILES);

void setup_hamiltonian_parameters_finite_q(int*, int*, int*, int*, int*, FERMI*, JOB_PARAM*, FILES);

void diagonalise_bse_tda_hamiltonian_crystal(int*, int*, Complex*, FERMI*, JOB_PARAM*, FILES);

void overlap_matrices5(ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, KQPOINT_TRAN*, int*, int*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, ATOM*, ATOM_TRAN*, SHELL*, REAL_LATTICE*, FERMI*, SYMMETRY*, CRYSTAL*, MPI_File*, FILES, JOB_PARAM*);

void overlap_matrices6(ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, KQPOINT_TRAN*, int*, int*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, ATOM*, ATOM_TRAN*, SHELL*, REAL_LATTICE*, FERMI*, SYMMETRY*, CRYSTAL*, MPI_File*, FILES, JOB_PARAM*);

void generate_kq_pairs(int, KQPOINT_TRAN*, CRYSTAL*, SYMMETRY*, FERMI*, JOB_PARAM*, FILES);

void select_kq_pairs(int*, int*, int*, int*, KPOINT_TRAN*, FERMI*, SYMMETRY*, CRYSTAL*, Q_LATTICE*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void read_scf_GW_eigenvalues_crystal(double*, FERMI*, ATOM*, char*, JOB_PARAM*, FILES);

void read_bse_eigenvalues_crystal(double*, int, char*, JOB_PARAM*, FILES);

void issue_MPI_IRecv(MPI_Request*,Complex*, int*, int*, int*, int*, int*, FERMI*, CRYSTAL*, FILES, JOB_PARAM*);

/*
void optical_spectrum_crystal1(FERMI*, ATOM*, ATOM_TRAN*_p, int*, int*, int [][2], SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void optical_spectrum_crystal2(FERMI*, ATOM*, ATOM_TRAN*_p, int*, int*, int [][2], SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void density_fitting_crystal_contract_integrals(int*, int*, int*, int*, int*, PAIR_TRAN*, KPOINT_TRAN*, MPI_File, Complex**, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, SYMMETRY*, Q_LATTICE*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void density_fitting_crystal_rotate_integrals1(int*, int*, int, KQPOINT_TRAN*, ComplexMatrix*, ComplexMatrix*, Complex**, Complex*, int*, int*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, ATOM*, ATOM_TRAN*, SHELL*, ATOM*, SHELL*, REAL_LATTICE*, FERMI*, CRYSTAL*, SYMMETRY*, MPI_File*, FILES, JOB_PARAM*);
*/

void density_fitting_crystal_rotate_integrals(int*, int, KQPOINT_TRAN*, ComplexMatrix*, ComplexMatrix*, Complex**, Complex*, int*, int*, int*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, FILES, JOB_PARAM*);

/*
void overlap_matrices(int*, ComplexMatrix*, ComplexMatrix*, KQPOINT_TRAN*, int*, int*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, ATOM*, ATOM_TRAN*, SHELL*, ATOM*, REAL_LATTICE*, FERMI*, SYMMETRY*, CRYSTAL*, MPI_File*, FILES, JOB_PARAM*);

void overlap_matrices1(int*, ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, KQPOINT_TRAN*, int*, int*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, ATOM*, ATOM_TRAN*, SHELL*, ATOM*, REAL_LATTICE*, FERMI*, SYMMETRY*, CRYSTAL*, MPI_File*, FILES, JOB_PARAM*);

void overlap_matrices2(ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, KQPOINT_TRAN*, int*, int*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, ATOM*, ATOM_TRAN*, SHELL*, ATOM*, REAL_LATTICE*, FERMI*, SYMMETRY*, CRYSTAL*, MPI_File*, FILES, JOB_PARAM*);

void overlap_matrices3(ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, KQPOINT_TRAN*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, ATOM*, ATOM_TRAN*, SHELL*, REAL_LATTICE*, FERMI*, SYMMETRY*, CRYSTAL*, MPI_File*, FILES, JOB_PARAM*);

void overlap_matrices4(ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, KQPOINT_TRAN*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, ATOM*, ATOM_TRAN*, SHELL*, REAL_LATTICE*, FERMI*, SYMMETRY*, CRYSTAL*, MPI_File*, FILES, JOB_PARAM*);
*/

/*
void overlap_matrices7(ComplexMatrix*, ComplexMatrix*, ComplexMatrix*, KQPOINT_TRAN*, int*, int*, int*, int*, ATOM*, ATOM_TRAN*, SHELL*, REAL_LATTICE*, FERMI*, SYMMETRY*, CRYSTAL*, MPI_File*, FILES, JOB_PARAM*);
*/
/*
void select_kq_zero_pairs(int*, int*, int*, int*, int*, int*, int*, int*, int*, KPOINT_TRAN*, FERMI*, SYMMETRY*, CRYSTAL*, Q_LATTICE*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void generate_coulomb_matrix_inverse_complex(ComplexMatrix*, int, FERMI*, ATOM_TRAN*, ATOM*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

*/
//void block_cyclic_to_linear_limit_complex(int*, int*, int*, int, Complex*, char*, JOB_PARAM*, FILES);

//void initialise_spk_grid_crystal1(int*, FERMI*, int, int*, int*, int*, JOB_PARAM*, FILES);

//void setup_procs(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, FERMI*, FILES, JOB_PARAM*);

/*
void wavefunction_product_density_fit_crystal_test1(ComplexMatrix*, Q_LATTICE*, FERMI*, Complex**, int*, int*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void orbital_product_density_fit_crystal_test1(ComplexMatrix*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void rotate_bra_ket_integrals1(int*, int*, int*, KQPOINT_TRAN*, Complex**, Complex*, Complex*, int*, int*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, ATOM*, ATOM_TRAN*, SHELL*, ATOM*, SHELL*, REAL_LATTICE*, FERMI*, CRYSTAL*, SYMMETRY*, MPI_File*, FILES, JOB_PARAM*);

void q3_loop_spin_polarised(double*, KQPOINT_TRAN*, Complex**, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, INT_1E*, PAIR_TRAN, double*, ATOM*, ATOM_TRAN*, SHELL*, ATOM*, SHELL*, REAL_LATTICE*, FERMI*, CRYSTAL*, SYMMETRY*, MPI_File*, FILES, JOB_PARAM*);

void integrals_occ_occ_vir_vir2_crystal(int*, int*, KPOINT_TRAN*, MPI_File, Complex*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, SYMMETRY*, Q_LATTICE*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void rotate_q_crystal2(int, int*, int*, int*, int*, int*, int*, CRYSTAL*, SYMMETRY*, FERMI*, JOB_PARAM*, FILES);

void rpa_hamiltonian_in_core_crystal1(int*, int*, int*, Complex*, Complex*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void contract_coulomb_integrals_complex(int*, int*, ComplexMatrix*, Complex*, Complex*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void rpa_hamiltonian_in_core_crystal(int*, int*, int*, Complex*, Complex*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);
*/
#endif
