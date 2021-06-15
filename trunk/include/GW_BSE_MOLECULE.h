
void gw_molecule(FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*,CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void bse_molecule(FERMI*, ATOM*, ATOM_TRAN*_p, int*, int*, int [][2], SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void bse_hamiltonian_in_core1(int*, int*, double*, double*, double*, double*, FERMI*, ATOM_TRAN*, ATOM*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void rpa_molecule1(FERMI*, ATOM*, ATOM_TRAN*_p, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

/*
void optical_spectrum_molecule(FERMI*, ATOM*, ATOM_TRAN*_p, int*, int*, int [][2], SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void bse_hamiltonian(int*, int*, double*, double*, FERMI*, ATOM_TRAN*, ATOM*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void bse_hamiltonian_in_core(int*, int*, double*, double*, double*, double*, FERMI*, ATOM_TRAN*, ATOM*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

*/
void diagonalise_bse_hamiltonian(int*, int*, double*, double*, FERMI*, JOB_PARAM*, FILES);

void diagonalise_bse_tda_hamiltonian(int*, int*, double*, FERMI*, JOB_PARAM*, FILES);

void diagonalise_rpa_hamiltonian(double*, int*, int*, FERMI*, ATOM*, ATOM*, JOB_PARAM*, FILES);

/*
void gw_bse_molecule_integrals(int*, int*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*,CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void gw_bse_molecule_integrals_in_core(double*, double*, int*, int*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*,CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

*/
void casida_diagonal_molecule_mpi_spk(int*, int*, int*, int*, FERMI*, double *, double*, JOB_PARAM*, FILES);

/*
void hamiltonian_ia_jb(double*, double, DoubleMatrix*, FERMI*, int*, int*, ATOM*, JOB_PARAM*, FILES);

void hamiltonian_ij_ab(double*, DoubleMatrix*, FERMI*, int*, int*, ATOM*, JOB_PARAM*, FILES);

void hamiltonian_ia_jb1(double*, double, DoubleMatrix*, FERMI*, int*, int*, ATOM*, JOB_PARAM*, FILES);

void hamiltonian_ib_ja(double*, DoubleMatrix*, FERMI*, int*, int*, ATOM*, JOB_PARAM*, FILES);

*/
void hamiltonian_in_core_coulomb_exchange_energy(int*, int*, int*, int*, int*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);

/*
void hamiltonian_screened_ij_ab(double*, DoubleMatrix*, FERMI*, int*, int*, ATOM*, JOB_PARAM*, FILES);

void hamiltonian_screened_ib_ja(double*, DoubleMatrix*, FERMI*, int*, int*, ATOM*, JOB_PARAM*, FILES);

void generate_temp4(double*, DoubleMatrix*, FERMI*, int*, int*, ATOM*, JOB_PARAM*, FILES);

*/
void generate_temp4_in_core(double*, double*, FERMI*, ATOM*, JOB_PARAM*, FILES);

/*
void self_energy(DoubleMatrix*, double*, DoubleMatrix*, FERMI*, int*, int*, ATOM*, JOB_PARAM*, FILES);

void count_ij_ab_array_size(IntMatrix*, int*, int*, int*, int*, int*, int*, int*, int*, FERMI*, JOB_PARAM*, FILES);

void count_ia_jb_array_size(IntMatrix*, int*, int*, int*, int*, int*, int*, int*, int*, FERMI*, JOB_PARAM*, FILES);

void count_ib_ja_array_size(IntMatrix*, int*, int*, int*, int*, int*, int*, int*, int*, FERMI*, JOB_PARAM*, FILES);

void senddata_ij_ab(double*, double*, int*, int*, int*, int*, int*, int*, FERMI*, ATOM*, JOB_PARAM*, FILES);

void senddata_ia_jb(double*, double*, int*, int*, int*, int*, int*, int*, FERMI*, ATOM*, JOB_PARAM*, FILES);

void senddata_ib_ja(double*, double*, int*, int*, int*, int*, int*, int*, FERMI*, ATOM*, JOB_PARAM*, FILES);

void senddata_screened_ij_ab(DoubleMatrix*, DoubleMatrix*, double*, int*, int*, int*, int*, int*, int*, FERMI*, JOB_PARAM*, FILES);

void senddata_screened_ib_ja(DoubleMatrix*, DoubleMatrix*, double*, int*, int*, int*, int*, int*, int*, FERMI*, JOB_PARAM*, FILES);

void recvdata_buffer(double*, IntMatrix*, int*, int*, double*, int*, int*, int*, int*, int*, double*, int*, int*, int*, int*, int*, int*, FERMI*, JOB_PARAM*, FILES);

void integrals_occ_vir_occ_vir(int*, int*, int, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void integrals_occ_occ_vir_vir(int*, int*, int, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void integrals_self_energy(int*, int*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void block_cyclic_to_linear_limit(int*, int*, int*, int, double*, char*, JOB_PARAM*, FILES);

void integrals_occ_occ_vir_vir1(int*, MPI_File, double*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

*/
/*
void rpa_bse_hamiltonian_in_core(double*, int*, int*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void count_arrays(int*, int*, int*, int*, int*, int*, ATOM*, FERMI*, JOB_PARAM*, FILES);

void contract_coulomb_integrals(int*, int*, DoubleMatrix*, double*, double*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void ring_size_transfer1(double*, double*, int, int, int, JOB_PARAM*, FILES);

void ring_rotate_once1(double**, int, int, int*, int*, int*, JOB_PARAM*, FILES);

void initialise_ring_rotate(int, int*, int*, int*, JOB_PARAM*, FILES);

*/
void hamiltonian_in_core_ia_jb(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);

void hamiltonian_in_core_ij_ab(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);

void hamiltonian_in_core_ib_ja(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);

void hamiltonian_in_core_screened_ij_ab(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);

void hamiltonian_in_core_screened_ij_ab1(int*, int*, int*, int*, int*, double*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);

void hamiltonian_in_core_screened_ij_ab3(int*, int*, int*, int*, int*, double*, double*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);

void hamiltonian_in_core_screened_ij_ab_and_ib_ja(int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);

/*
void self_energy_diagonal_in_core(int*, int*, int*, int*, int*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);
*/
void self_energy_diagonal_in_core1(int*, int*, int*, int*, int*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);
/*
void hamiltonian_in_core_screened_ij_ab2(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, ATOM*, FERMI*, JOB_PARAM*, FILES);

void generate_coulomb_matrix_inverse(DoubleMatrix*, ATOM_TRAN*, ATOM*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

*/
void count_integral_buffer_sizes(int*, FERMI*, ATOM*, JOB_PARAM*, FILES);

void read_scf_GW_eigenvalues(double*, int, int, char*, JOB_PARAM*, FILES);

void read_SCF_GW_eigenvalues(double*, int, char*, JOB_PARAM*, FILES);

void read_write_scf_eigenvectors(FERMI*, ATOM*, JOB_PARAM*, FILES);

