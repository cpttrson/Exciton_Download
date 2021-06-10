#ifndef INTEGRALSH
#define INTEGRALSH
/*
void integrals_crystal_ijkl(INTEGRAL_LIST*, PAIR_TRAN*, int*, PAIR_TRAN*, int*, int*, int*, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_crystal_ijkl3(INTEGRAL_LIST*, QUAD_TRAN*, int*, int*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_molecule_ijkl(INTEGRAL_LIST*, PAIR_TRAN*, int*, PAIR_TRAN*, int*, int*, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_molecule_ijkl1(INTEGRAL_LIST*, int, int, int, int, int*, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

//void integrals_molecule_ijkl1(INTEGRAL_LIST*, int, int, int, int, int*, int*, PAIR_TRAN*, double*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_exchange_ijkl1(INTEGRAL_LIST*, QUAD_TRAN*, int*, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void pack_write_coulomb_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void read_unpack_coulomb_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void pack_write_exchange_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void read_unpack_exchange_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void pack_write_molecule_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void read_unpack_molecule_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);


void integrals_crystal_ijkl1(double*, PAIR_TRAN*, double*, double*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, ATOM_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_crystal_recip_ijkl(INTEGRAL_LIST*, PAIR_TRAN*, int*, PAIR_TRAN*, int*, int*, int*, int*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void A_coefficients(int, int, int, int, double*, double*, double*, double*, SHELL*, GAUSSIAN*, FILES);
*/
/*
void E_coefficients(int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);
*/
/*
void C_coefficients(int*, int*, int*, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, double*, double*, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void C_coefficients1(int*, int*, int*, int, int, int, int, double*, double*, double*, double*, double*, double*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);
void E_coefficients1(int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, VECTOR_DOUBLE*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);
void mcmurchie_davidson(double*, int, int, int, int,  double*, double*, double*, double*, double*, double*, double*, SHELL*, JOB_PARAM*, FILES);

//void mcmurchie_davidson(double*, double*, double*, int, int, int, int,  double*, double*, double*, double*, double*, double*, double*, SHELL*, JOB_PARAM*, FILES);

void mcmurchie_davidson_screen(double*, int, int, int, int, int, double*, double*, double*, double*, SHELL*, JOB_PARAM*, FILES);

void mcmurchie_davidson_screen_complex(Complex*, int, int, int, int, int, double*, double*, double*, Complex*, SHELL*, JOB_PARAM*, FILES);

void four_centre_cartesian_to_sh(double*, double*, int, int, int, int, SHELL*, JOB_PARAM*, FILES);

*/
/*
void integrals_2e_crystal_nosym2(INTEGRAL_LIST *, PAIR_TRAN *, PAIR_TRAN *, int *, int *, int *, REAL_LATTICE *, ATOM *, SHELL *, GAUSSIAN *, SYMMETRY *, CRYSTAL *, JOB_PARAM *, FILES);

void pack_write_molecular_2e_integrals(INTEGRAL_LIST*, FILE*, FILES);

void read_unpack_molecular_2e_integrals(INTEGRAL_LIST*, FILE*, FILES);

void pack_write_molecular_2e_integrals1(INTEGRAL_LIST*, int, int*, FILE*, JOB_PARAM*, FILES);

void read_unpack_molecular_3c_integrals1(INTEGRAL_LIST*, int, FILE*, JOB_PARAM*, FILES);

void pack_write_molecule_3c_integrals(INTEGRAL_LIST*, FILE*, JOB_PARAM*, FILES);

void read_unpack_molecule_3c_integrals(INTEGRAL_LIST*, FILE*, JOB_PARAM*, FILES);
*/
#endif
