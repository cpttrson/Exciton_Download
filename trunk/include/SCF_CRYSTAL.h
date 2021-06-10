#ifndef SCF_CRYSTALH
#define SCF_CRYSTALH

void scf_crystal(FERMI*, ATOM*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

/*
void C_DIIS_extrapolation(ComplexMatrix**, ComplexMatrix**, ComplexMatrix**, ComplexMatrix**, ComplexMatrix**, int*, int*, SYMMETRY*, JOB_PARAM*, FILES);

void exchange_matrix_elements(double*, double*, double*, double*, double*, double*, double*, double*, KPOINT_TRAN*, ComplexMatrix*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void build_kohn_sham_matrix(double*, DFT_GRID*, INT_1E*, FERMI*, double*, double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM*, FILES);

//void build_kohn_sham_matrix(double*, double*, double*, DFT_GRID*, INT_1E*, FERMI*, double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM*, FILES);

void build_final_kohn_sham_matrix(double*, DFT_GRID*, INT_1E*, FERMI*, double*, double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM*, FILES);

//void build_final_kohn_sham_matrix(double*, double*, double*, DFT_GRID*, INT_1E*, FERMI*, double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM*, FILES);

void coulomb_exchange_matrix_crystal(double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void exchange_matrix_crystal1(double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void exchange_matrix_crystal2(double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void exchange_matrix_crystal3(double*, double*, double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void density_fitting_integrals_molecule_ax1(FERMI*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void coulomb_exchange_matrix_crystal_ax2(double*, double*, ComplexMatrix*, FERMI*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void coulomb_exchange_matrix_crystal_ax1(double*, double*, ComplexMatrix*, FERMI*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void coulomb_exchange_matrix_crystal_ax(double*, double*, FERMI*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void coulomb_exchange_matrix_crystal1(double*, double*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void exchange_correlation_matrix(double*, double*, DFT_GRID*, double*, double*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void contract_coulomb_integrals(double*, INTEGRAL_LIST*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void contract_coulomb_integrals1(REAL_LATTICE*, double*, INTEGRAL_LIST*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, REAL_LATTICE_TABLES*,JOB_PARAM*,FILES);

void contract_coulomb_integrals2(double*, INTEGRAL_LIST*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, REAL_LATTICE_TABLES*,JOB_PARAM*,FILES);

void contract_exchange_integrals(double*, INTEGRAL_LIST*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, REAL_LATTICE_TABLES*,JOB_PARAM*,FILES);

void contract_exchange_integrals2(double*, INTEGRAL_LIST*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void contract_exchange_integrals3(double*, INTEGRAL_LIST*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_coulomb_integrals_crystal(PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void generate_coulomb_integrals_molecule(PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void contract_exchange_integrals1(double*, INTEGRAL_LIST*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*,JOB_PARAM*,FILES);

void contract_molecule_integrals(double*, double*, INTEGRAL_LIST*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void contract_molecule_integrals1(double*, double*, INTEGRAL_LIST*, double*, double*, PAIR_TRAN*, QUAD_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void write_JOB_PARAM_array(FERMI*, int*, double*, JOB_PARAM*, FILES);

void read_JOB_PARAM_array(FERMI*, int*, double*, JOB_PARAM*, FILES);
*/

#endif
