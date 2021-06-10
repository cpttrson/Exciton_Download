#ifndef DFTH
#define DFTH

void count_dft_savin_grid(int*, PAIR_TRAN*, ATOM*, CRYSTAL*, REAL_LATTICE*, SYMMETRY*, JOB_PARAM*, FILES);

void generate_dft_savin_grid(DFT_GRID*, PAIR_TRAN*, ATOM*, CRYSTAL*, REAL_LATTICE*, SYMMETRY*, JOB_PARAM*, FILES);

void generate_dft_rho_grid(double*, double*, DFT_GRID*, TRIPLE_TRAN*, PAIR_TRAN*, REAL_LATTICE*, REAL_LATTICE_TABLES*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL *crystal, JOB_PARAM*, FILES);

void generate_orbital_grid(DFT_GRID*, double*, double*, TRIPLE_TRAN*, PAIR_TRAN*, REAL_LATTICE*, REAL_LATTICE_TABLES*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, SYMMETRY*, JOB_PARAM*, FILES);

void vxc_grid_dft(double*, DFT_GRID*, double*, JOB_PARAM*, FILES);

void exchange_matrix_dft(double*, double*, DFT_GRID*, double*, double*, TRIPLE_TRAN*, REAL_LATTICE*, REAL_LATTICE_TABLES*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL *crystal, JOB_PARAM*, FILES);

void dft_grid_parameters(int *L_range, double *r_range, JOB_PARAM *job, FILES file);

#endif
