#ifndef INTEGRALS1C
#define INTEGRALS1C

void nuclear_repulsion_energy(double*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void fock_element_1e1(INT_1E*, int, PAIR_TRAN*, int, int[], REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void fock_element_elecnuc(double*, PAIR_TRAN*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void Kinetic(double*, int, int, int, int, int, int, int, int, int,  double*, double*, double*, double*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void ElecNuc(double*,int,int,int,int,int,int,int,double*,double*,double*,VECTOR_DOUBLE*,double*,double*,REAL_LATTICE*,RECIPROCAL_LATTICE*,CRYSTAL*,ATOM*,SHELL*,JOB_PARAM*,FILES);

void Momentum(double*, int, int, int, int, int, int, int, int, int, int,  double*, double*, double*, double*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void Overlap(double*, int, int, int, int, int, int, int, double*, double*, double*, double*, SHELL*, JOB_PARAM*, FILES);

void Dipole(double*,int,int,int,int,int,int,int,int,int,int,int,double*,double*,double*,double*,VECTOR_DOUBLE*,ATOM*,SHELL*,GAUSSIAN*,JOB_PARAM*,FILES);

void two_centre_coulomb(double*,int,int,int,int,int,int,int,double*,double*,double*,VECTOR_DOUBLE*,double*,double*,RECIPROCAL_LATTICE*,CRYSTAL*,ATOM*,SHELL*,GAUSSIAN*,JOB_PARAM*,FILES);
//CHANGES2015void Coulomb(double*,int,int,int,int,int,int,int,double*,double*,double*,VECTOR_DOUBLE*,double*,double*,RECIPROCAL_LATTICE*,CRYSTAL*,ATOM*,SHELL*,GAUSSIAN*,JOB_PARAM*,FILES);

//void two_centre_coulomb1(INT_1E*, PAIR_TRAN*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

//void two_centre_coulomb1_crystal(ComplexMatrix*, REAL_LATTICE*, Q_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void two_centre_exchange_crystal(ComplexMatrix*, REAL_LATTICE*, Q_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void two_centre_exchange_crystal1(int, Complex*, PAIR_TRAN*, REAL_LATTICE*, Q_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void two_centre_overlap(INT_1E*, PAIR_TRAN*, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, JOB_PARAM*, FILES);

void E_coefficients1(int,int,int,int,int,int,int,int,int,int,double*,double*,double*,double*,double*,double*,double*,VECTOR_DOUBLE*,REAL_LATTICE*,ATOM*,SHELL*,GAUSSIAN*,JOB_PARAM*,FILES);
 
void E_coefficients_1c(int, int, double*, double*, double*, double*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

void E_coefficients_2c(int,int,int,int,int,int,int,int,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,REAL_LATTICE*,ATOM*,SHELL*,GAUSSIAN*,JOB_PARAM*,FILES);

void mcmurchie_davidson_2c(double*,int,int,int,int,int,double*,double*,double*,double*,double*,double*,double*,SHELL*,JOB_PARAM*,FILES);

void mcmurchie_davidson_2c_complex(Complex*,Complex*,int,int,int,int,int,double*,double*,double*,double*,double*,double*,SHELL*,JOB_PARAM*,FILES);

void integrals_molecule_screen(double*, int, int, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_crystal_exchange_screen(double*, int, int, int, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_crystal_screen(double*, int, int, int, REAL_LATTICE*, RECIPROCAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void integrals_crystal_screen_complex(Complex*, int, int, int, REAL_LATTICE*, RECIPROCAL_LATTICE*, Q_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

#endif
