#ifndef ATOMSCFH
#define  ATOMSCFH

void atom_scf(ATOM *, int, DoubleMatrix *, double *, double *, int *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES);
void init_mix(DoubleMatrix *, DoubleMatrix *, DoubleMatrix *, DoubleMatrix *, DoubleMatrix *);
void mix_density(DoubleMatrix *,DoubleMatrix *,DoubleMatrix *,DoubleMatrix *,double *,double *,double *,int,int*,JOB_PARAM*, FILES);
void initial_density_matrix_atom(DoubleMatrix *, ATOM *, int *, double *, int *, SHELL *, JOB_PARAM *, FILES);
void fock_element_1e_atm(DoubleMatrix*, DoubleMatrix*, int*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);
//void pseudo_1e_atm(DoubleMatrix*, PSEUDO*, int*, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);
//void orbital_grid(double *, double *, double *, double *, double *, int *, int, int, ATOM *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES);
//void orbital_gradient_grid(double*, double*, double*, double*, VECTOR_DOUBLE*, int*, int, int, ATOM*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);
void integrals_2e_atom(INTEGRAL_LIST *, int, ATOM *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES);
void fock_2e_matrix(DoubleMatrix*, INTEGRAL_LIST*, DoubleMatrix*, int, ATOM*, JOB_PARAM*, FILES);
//void kohn_sham_2e_matrix(int *, DoubleMatrix *, INTEGRAL_LIST*, DoubleMatrix *, ATOM *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES);
//void pseudopotential_1e_matrix(DoubleMatrix *, int *, PSEUDO *, ATOM *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES);
//void coulomb_kohn_sham_2e_matrix(int *, DoubleMatrix *, DoubleMatrix*, DoubleMatrix *, ATOM *, SHELL *, GAUSSIAN *, JOB_PARAM *, FILES);
void total_energy_atom(DoubleMatrix*, DoubleMatrix*, DoubleMatrix*, int, double*, double* , ATOM*, JOB_PARAM*, FILES);
void final_total_energy_atom(DoubleMatrix*, DoubleMatrix*, DoubleMatrix*, int, double*, double* , ATOM*, JOB_PARAM*, FILES);
double e_atom(int , int , int , double );
//double e(int , int , int , double , double , double );
//void f000m(double *, double , double , int );
double ftuvn_atom(int , int , int , int , double *);
void g000m(int , double *, double , int );
//void Gauss_Chebyshev_weights_abscissa(int , double *, double *, JOB_PARAM*, FILES);
//void Gauss_Legendre_weights_abscissa(int , double *, double *, JOB_PARAM*, FILES);

#endif
