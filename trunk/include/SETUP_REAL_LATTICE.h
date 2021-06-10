#ifndef SETUP_REAL_LATTICEH
#define SETUP_REAL_LATTICEH

void generate_real_lattice(CRYSTAL*, REAL_LATTICE*, RECIPROCAL_LATTICE*, JOB_PARAM *job, FILES file);

void generate_real_lattice_shells(CRYSTAL*, REAL_LATTICE*, JOB_PARAM*, FILES);

void generate_real_lattice_tables(CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

#endif
