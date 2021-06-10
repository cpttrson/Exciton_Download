void count_pairs4(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, JOB_PARAM*, FILES);

void select_pairs4(PAIR_TRAN*, ATOM*, REAL_LATTICE*, JOB_PARAM*, FILES);

void count_density_pairs4(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE *R, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void count_range_selected_pairs(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE *R, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_density_pairs4(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_range_selected_pairs(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void count_triples1(TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_triples1(TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void count_triples1_reversed(int, TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_triples1_reversed(int, TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_molecule_quads(PAIR_TRAN*, QUAD_TRAN* ,int, int, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_c_quads(PAIR_TRAN*, QUAD_TRAN*, int, int, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_e_quads(PAIR_TRAN*, QUAD_TRAN*, int, int, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void print_pairs(PAIR_TRAN*, ATOM*, REAL_LATTICE*, JOB_PARAM*, FILES);

/*
void generate_triples(TRIPLE_TRAN*, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void count_triples2_reversed(int, TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_triples2_reversed(int, TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void count_triples_q(int*, int*, int*, TRIPLE_TRAN*, PAIR_TRAN*, KPOINT_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_triples_q(int*, int*, int*, int*, TRIPLE_TRAN*, PAIR_TRAN*, KPOINT_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);


void generate_c_quads4(PAIR_TRAN*, PAIR_TRAN*,PAIR_TRAN*,int, int, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_e_quads3(PAIR_TRAN*, PAIR_TRAN*,PAIR_TRAN*,int, int, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_e_quads4(PAIR_TRAN*, PAIR_TRAN*,PAIR_TRAN*,int, int, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

int diffvec(int, int, REAL_LATTICE*, JOB_PARAM*, FILES);

void count_density_pairs3(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void count_density_pairs4(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE *R, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_density_pairs3(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_density_pairs4(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

*/
