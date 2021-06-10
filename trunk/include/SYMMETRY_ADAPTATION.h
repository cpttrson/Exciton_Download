#ifndef SYMMETRY_ADAPTATIONH
#define SYMMETRY_ADAPTATIONH

void count_atom_irrep(IntMatrix*, ATOM*, ATOM_TRAN*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void count_atom_irrep_crystal(IntMatrix*, ATOM*, ATOM_TRAN*, REAL_LATTICE_TABLES*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void count_basis_irrep(int*, ATOM*, ATOM_TRAN*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void count_basis_irrep_crystal(int*, ATOM*, ATOM_TRAN*, REAL_LATTICE_TABLES*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void count_atom_salc(int, int, SALC*, ATOM*, ATOM_TRAN*, PAIR_TRAN*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void count_atom_salc_crystal(int, IntMatrix*, SALC*, ATOM*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE_TABLES*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);
//void count_atom_salc_crystal(int, int, int, SALC*, ATOM*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE_TABLES*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void generate_atom_salc(int, int, SALC*, ATOM*, ATOM_TRAN*, PAIR_TRAN*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void generate_atom_salc_crystal(int, IntMatrix*, SALC*, ATOM*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE_TABLES*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);
//void generate_atom_salc_crystal(int, int, int, SALC*, ATOM*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE_TABLES*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void count_shell_salc(int, int, int, SALC*, ATOM*, ATOM_TRAN*, PAIR_TRAN*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void generate_shell_salc(int, int, int, SALC*, ATOM*, ATOM_TRAN*, PAIR_TRAN*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void transform_salc_to_atomic_basis(ComplexMatrix*, double*, int*, ComplexMatrix**, KPOINT_TRAN*, int*, ATOM_TRAN*, PAIR_TRAN*, REAL_LATTICE*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void transform_degenerate_salc(ComplexMatrix**, double**, int*, int*, KPOINT_TRAN*, int *nk, SYMMETRY*, JOB_PARAM*, FILES);

#endif
