double e(int, int, int, double, double, double);

double ftuvn(int, int, int, int, double *, VECTOR_DOUBLE);

double cosfactor(int, double);

Complex cosfactor_complex(int, double);

void non_recursive_ftuvn(int, int, double f[][13][13][13], double en[][55], VECTOR_DOUBLE *);

void erf_derivative(double*, int, double, double, RECIPROCAL_LATTICE*);

void erf_exp_derivative(double*, double*, int, double, double, RECIPROCAL_LATTICE*);

void erfc_derivative(double*, int, double, double, RECIPROCAL_LATTICE*);

/*
double e3(int, int, int, int, double, double, double, double);

void non_recursive_e(double[][9][9], int, double, double, double);

void non_recursive_E(double*, double*, int, int, int, double*, double*, double*, SHELL*, FILES);
*/
