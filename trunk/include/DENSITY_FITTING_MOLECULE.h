
void generate_integral_buffers_molecule_ija(double*, double*, int*, int*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*,CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);
//void integrals_molecule_IJ_alpha(double*, double*, int*, int*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*,CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void contract_integrals_molecule_ija(int*, MPI_File, double*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);
//void contract_integrals_molecule_ij_alpha(int*, MPI_File, double*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void generate_coulomb_matrix_inverse(DoubleMatrix*, ATOM_TRAN*, ATOM*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);
//void generate_coulomb_matrix_inverse(DoubleMatrix *V_inv_k1, FERMI *fermi, ATOM_TRAN *atom_p, ATOM *atoms, ATOM *atoms_ax, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file);

void contract_coulomb_integrals(int*, int*, DoubleMatrix*, double*, double*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void initialise_ring_rotate(int, int*, int*, int*, JOB_PARAM*, FILES);

void ring_rotate_once(double**, int, int, int*, int*, int*, JOB_PARAM*, FILES);

void ring_transfer(double*, double*, int, int, int, JOB_PARAM*, FILES);

