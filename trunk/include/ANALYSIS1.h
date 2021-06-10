#ifndef ANALYSIS1H
#define  ANALYSIS1H
void fermi_surface(int is[3], int band[2], CRYSTAL *, SYMMETRY *, ATOM *, JOB_PARAM *, FILES);
void state_density(int [3], int [2], int, int, double *, double *, char, ATOM *,  ATOM_TRAN *, ATOM_TRAN *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, RECIPROCAL_LATTICE *G, JOB_PARAM *, FILES);
void state_density1(int [3], int [2], int, int, double *, double *, char, ATOM *,  ATOM_TRAN *, ATOM_TRAN *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, RECIPROCAL_LATTICE *G, JOB_PARAM *, FILES);
void susceptibility(int [3], int [2], int, double [2], VECTOR_DOUBLE*, ATOM *, ATOM_TRAN *, ATOM_TRAN *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, RECIPROCAL_LATTICE *, JOB_PARAM *, FILES);
void susceptibility1(int [3], int [2], int, double [2], VECTOR_DOUBLE*, ATOM *, ATOM_TRAN *, ATOM_TRAN *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, RECIPROCAL_LATTICE *, JOB_PARAM *, FILES);
void susceptibility2(int [3], int [2], int, double [2], VECTOR_DOUBLE*, ATOM *, ATOM_TRAN *, ATOM_TRAN *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, RECIPROCAL_LATTICE *, JOB_PARAM *, FILES);
void susceptibilitya(int [3], int [2], int, double [2], VECTOR_DOUBLE*, ATOM *, ATOM_TRAN *, ATOM_TRAN *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, RECIPROCAL_LATTICE *, JOB_PARAM *, FILES);
void raman_suscep(int [3], int [2], int, double [2], VECTOR_DOUBLE*, ATOM *, ATOM_TRAN *, ATOM_TRAN *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, RECIPROCAL_LATTICE *, JOB_PARAM *, FILES);
void dielectric_function(int [3], int [2], int, double [2], VECTOR_DOUBLE*, ATOM *, ATOM_TRAN *, ATOM_TRAN *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, RECIPROCAL_LATTICE *, JOB_PARAM *, FILES);
//void tetrahedron_spectrum(int, int, int, int is[4], KPT_TRAN *, int b[2], double e[2], double *, ATOM *, CRYSTAL *, JOB_PARAM *, FILES);
void generate_polarisability(Complex *, int[3], int[2], int, PAIR_TRAN *, int, ATOM *, ATOM_TRAN *, ATOM_TRAN *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, RECIPROCAL_LATTICE *, JOB_PARAM *, FILES);
void contract_three_centre_overlap(Complex *, double *, double *, double *, int [3], int [2], int, FERMI *, PAIR_TRAN *, int , PAIR_TRAN *, int , ATOM *, SHELL *, GAUSSIAN *, CRYSTAL *, SYMMETRY *, REAL_LATTICE *, REAL_LATTICE_TABLES *, JOB_PARAM *, FILES);
void polygons_vertices_vol(int [], int *, int *, double *, CRYSTAL *, JOB_PARAM *, FILES file);
#endif
