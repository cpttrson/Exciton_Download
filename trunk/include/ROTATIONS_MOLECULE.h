
void rotate_single(int, int, Complex*, Complex*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_triple_ija(int*, int*, int*, int*, int*, double*, double*, ATOM_TRAN*, ATOM*, SHELL*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);
//void rotate_permute_triple_ij_alpha(int*, int*, int*, int*, int*, double*, double*, ATOM_TRAN*, ATOM*, SHELL*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_expand_pair(int, PAIR_TRAN*, double*, double*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_expand_pair_complex(int, PAIR_TRAN*, Complex*, Complex*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_expand_tensor1(int, PAIR_TRAN*, double*, double*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

//void rotate_permute_expand_tensor_3(double*, double*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_sum_block(double*, double*, int, int, int, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_psi(Complex*, Complex*, int, int, KPOINT_TRAN*, ATOM_TRAN*, ATOM*, REAL_LATTICE*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

