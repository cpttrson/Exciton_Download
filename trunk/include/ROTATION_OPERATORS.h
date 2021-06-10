#ifndef ROTATION_OPERATORSH
#define ROTATION_OPERATORSH
/*
void generate_rotation_operators(SYMMETRY *symmetry, JOB_PARAM *job, FILES file);

void generate_rotation_operators_ivanic_ruedenberg(SYMMETRY *symmetry, JOB_PARAM *job, FILES file);

//void ivanic_ruedenberg(double*, int, JOB_PARAM *job, FILES file);

void crystal_sh_order(double[5][5]);

void uvw(const int, const int, const int, double*, double*, double*);

inline double delta(const int, const int);

void R_offset(const int, double*, double*);

int M_offset(const int, const int, const int);

double M(const int, const int, const int, double*);

double P(const int, const int, const int, const int, double*);

double U(const int, const int, const int, double*);

double V(const int, const int, const int, double*);

double W(const int, const int, const int, double*);

void test_rotation_operators(ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform_salc(double*, ComplexMatrix**, KPOINT_TRAN*, int*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform_salc1(double*, ComplexMatrix**, KPOINT_TRAN*, int*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform_salc_crystal(double*, ComplexMatrix**, KPOINT_TRAN*, int*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE*, REAL_LATTICE_TABLES*, ATOM*, SHELL*, SYMMETRY*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform(double*, Complex*, KPOINT_TRAN*, int*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform_molecule(double*, double*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void fourier_transform_3(double*, Complex*, KPOINT_TRAN*, int*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*,FILES);

void rotate_permute_expand_tensor_3(double*, double*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_psi(Complex*, Complex*, int, int, KPOINT_TRAN*, ATOM_TRAN*, ATOM*, REAL_LATTICE*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);
*/
//CHANGE2014void rotate_psi(Complex*, Complex*, int[2], int, KPOINT_TRAN*, ATOM_TRAN*, ATOM*, REAL_LATTICE*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);
void rotate_psi1(Complex*, Complex*, int[2], int, KPOINT_TRAN*, ATOM_TRAN*, ATOM*, REAL_LATTICE*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_reduced_density2(Complex*, double*, PAIR_TRAN*, int*, KPOINT_TRAN*, int*, int*, REAL_LATTICE*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_expand_pair(int, PAIR_TRAN*, double*, double*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

//void rotate_permute_expand_pair_complex(int, PAIR_TRAN*, Complex*, Complex*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_expand_triple(int, TRIPLE_TRAN*, double*, double*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_expand_triple_ax(int, TRIPLE_TRAN*, double*, double*, ATOM*, SHELL*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_expand_triple_ax_reversed(int, TRIPLE_TRAN*, double*, double*, ATOM*, SHELL*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_triple_ax_reversed(int, int, int, PAIR_TRAN*, double*, double*, ATOM_TRAN*, ATOM*, SHELL*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_triple_ax_reversed2(int*, int*, int*, int*, int*, double*, double*, ATOM_TRAN*, ATOM*, SHELL*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_pair_complex(int*, int*, int*, int*, Complex*, Complex*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_triple_ax_reversed2_complex(int*, int*, int*, int*, int*, Complex*, Complex*, ATOM_TRAN*, ATOM*, SHELL*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

/*
void rotate_single(int, int, Complex*, Complex*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_quad(int , int , int , int , int , double *, double *, ATOM *, SHELL *, SYMMETRY *, JOB_PARAM *, FILES);

void rotate_permute_expand_pair1(int, PAIR_TRAN*, double*, double*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void new_rotate_pair(int, PAIR_TRAN*, double*, double*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_permute_expand_tensor1(int, PAIR_TRAN*, double*, double*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void rotate_sum_block(double*, double*, int, int, int, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

*/
/*
void two_center_cartesian_to_sh_shell(double*, double*, int, int, int, int, int, int, int, int, SHELL*, JOB_PARAM*, FILES);

void two_center_cartesian_to_sh_shell_complex(Complex*, Complex*, int, int, int, int, int, int, int, int, SHELL*, JOB_PARAM*, FILES);

void two_center_cartesian_to_sh_shell1(double*, double*, int, int, int, int, int, int, int, int, SHELL*, JOB_PARAM*, FILES);

void two_center_cartesian_to_sh_shell1_complex(Complex*, Complex*, int, int, int, int, int, int, int, int, SHELL*, JOB_PARAM*,FILES);

void two_center_vector_cartesian_to_sh_shell(double*, double*, int, int, int, int, int, int, int, int, int, int, SHELL*, JOB_PARAM*, FILES);
*/
void three_center_cartesian_to_sh_shell(double*, double*, int, int, int, int, int, int, int, int, int, int, int, int, int, SHELL*, JOB_PARAM*, FILES);

void three_center_cartesian_to_sh_shell_ax(double*, double*, int, int, int, int, int, int, int, int, int, int, int, int, int, SHELL*, SHELL*, JOB_PARAM*, FILES);

void three_center_cartesian_to_sh_shell_ax_reversed(double*, double*, int, int, int, int, int, int, int, int, int, int, int, int, int, SHELL*, SHELL*, JOB_PARAM*, FILES);

void three_center_cartesian_to_sh_shell_ax_reversed_complex(Complex*, Complex*, int, int, int, int, int, int, int, int, int, int, int, int, int, SHELL*, SHELL*, JOB_PARAM*, FILES);

void new_rotate_tensor1_pair(int, PAIR_TRAN*, double*, double*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void count_atom_cluster(int, int, SALC*, int*, int*, int*, SYMMETRY*, ATOM_TRAN*, REAL_LATTICE_TABLES*);

void generate_atom_cluster(IntMatrix*, SALC*, int*, int*, int*, SYMMETRY*);

#endif
