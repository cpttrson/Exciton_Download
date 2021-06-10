
void generate_rotation_operators_ivanic_ruedenberg(SYMMETRY*, JOB_PARAM*, FILES);

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

